"""Cross-validation style stability checks for GRB model-selection results."""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median, quantiles, stdev
from typing import Any

import numpy as np

from .grb_constants import long_to_short, short_to_long
from .grb_core import GRB
from .grb_stability import ModelComparison
from .safe_good_best import pick_best_single_model

SIMPLE_TO_BB = {
    "PL": "PL_BB",
    "CPL": "CPL_BB",
    "BAND": "BAND_BB",
    "SBPL": "SBPL_BB",
}
DEFAULT_BB_PARAM_NAMES = ["amp_bb", "kt_bb"]
DEFAULT_ELIGIBLE_STATUSES = {"BEST", "SAFE"}
DEFAULT_GRB_NAMES = ["080916C", "140206B", "131014A", "231129C"]
RECORD_FIELDNAMES = [
    "grb",
    "grb_short",
    "interval",
    "interval_kind",
    "simple_model",
    "bb_model",
    "simple_status",
    "bb_status",
    "evaluable",
    "bb_supported",
    "rejection_reason",
    "delta_cstat",
    "delta_aic",
    "delta_bic",
    "lrt_supports_complex",
    "aic_supports_complex",
    "bic_supports_complex",
    "recommended_model",
    "selection_conflict",
    "fold",
    "logo_fold",
]


def load_results(results_path: str | Path = "results.json") -> dict[str, Any]:
    """Load a results.json-style file."""
    with Path(results_path).open("r", encoding="utf-8") as f:
        return json.load(f)


def normalize_grb_name(name: str) -> str:
    """Normalize short or long GRB names to the long results.json key."""
    return short_to_long.get(name, name)


def filter_results_by_grbs(
    data: dict[str, Any],
    grb_names: list[str] | tuple[str, ...] | None,
) -> dict[str, Any]:
    """Filter a results.json-style mapping to the requested GRBs."""
    if grb_names is None:
        return data

    selected = [normalize_grb_name(name) for name in grb_names]
    missing = [name for name in selected if name not in data]
    if missing:
        raise KeyError(f"Requested GRBs not found in results: {', '.join(missing)}")
    return {name: data[name] for name in selected}


def interval_kind(interval_name: str) -> str:
    """Return a coarse interval kind from a results.json interval key."""
    token = interval_name.split()[0]
    if token.startswith("TR"):
        return "TR"
    if token.startswith("BR"):
        return "BR"
    return token


def _status(model_payload: dict[str, Any] | None) -> str | None:
    if not model_payload:
        return None
    status = model_payload.get("_status")
    return str(status).upper() if status is not None else None


def _is_eligible_status(
    model_payload: dict[str, Any] | None,
    eligible_statuses: set[str],
) -> bool:
    return _status(model_payload) in eligible_statuses


def _numeric_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _build_models_for_interval(grb_name: str, interval_name: str, grb_payload: dict[str, Any]):
    grb = GRB.from_dictionary(name=grb_name, grb_data={interval_name: grb_payload[interval_name]})
    return next(iter(grb.intervals)).models


def _base_record(grb_name: str, interval_name: str) -> dict[str, Any]:
    return {
        "grb": grb_name,
        "grb_short": long_to_short.get(grb_name, grb_name),
        "interval": interval_name,
        "interval_kind": interval_kind(interval_name),
        "simple_model": None,
        "bb_model": None,
        "simple_status": None,
        "bb_status": None,
        "evaluable": False,
        "bb_supported": False,
        "rejection_reason": None,
        "delta_cstat": None,
        "delta_aic": None,
        "delta_bic": None,
        "lrt_supports_complex": None,
        "aic_supports_complex": None,
        "bic_supports_complex": None,
        "recommended_model": None,
        "selection_conflict": None,
    }


def collect_bb_selection_records(
    data: dict[str, Any],
    include_marginal: bool = False,
    bb_param_names: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Collect one BB-selection stability record per GRB interval.

    The simple model is the best eligible single model in the interval. Its
    matching +BB model is then evaluated with ModelComparison when available
    and eligible.
    """
    bb_param_names = bb_param_names or DEFAULT_BB_PARAM_NAMES
    eligible_statuses = set(DEFAULT_ELIGIBLE_STATUSES)
    if include_marginal:
        eligible_statuses.add("MARGINAL")

    records: list[dict[str, Any]] = []
    for grb_name in sorted(data):
        grb_payload = data[grb_name]
        for interval_name in sorted(grb_payload):
            interval_payload = grb_payload[interval_name]
            record = _base_record(grb_name, interval_name)

            simple_cstats = {}
            for simple_model in SIMPLE_TO_BB:
                payload = interval_payload.get(simple_model)
                if not _is_eligible_status(payload, eligible_statuses):
                    continue
                cstat = _numeric_or_none(payload.get("c-stat/dof", [None])[0])
                if cstat is not None:
                    simple_cstats[simple_model] = cstat

            if not simple_cstats:
                record["rejection_reason"] = "no eligible simple model"
                records.append(record)
                continue

            simple_model, _ = pick_best_single_model(simple_cstats)
            bb_model = SIMPLE_TO_BB[simple_model]
            simple_payload = interval_payload.get(simple_model)
            bb_payload = interval_payload.get(bb_model)

            record.update(
                {
                    "simple_model": simple_model,
                    "bb_model": bb_model,
                    "simple_status": _status(simple_payload),
                    "bb_status": _status(bb_payload),
                }
            )

            if bb_payload is None:
                record["rejection_reason"] = "matching BB model missing"
                records.append(record)
                continue
            if not _is_eligible_status(bb_payload, eligible_statuses):
                record["rejection_reason"] = f"matching BB model status is {_status(bb_payload)}"
                records.append(record)
                continue

            try:
                model_set = _build_models_for_interval(grb_name, interval_name, grb_payload)
                simple = model_set.get(simple_model)
                complex_ = model_set.get(bb_model)
                comparison = ModelComparison(simple, complex_, bb_param_names=bb_param_names)
                summary = comparison.summary_dict()
                evidence = summary["evidence_summary"]
                info = summary["information_criteria"]
                lrt = summary["lrt"]
            except Exception as exc:
                record["rejection_reason"] = f"comparison failed: {exc}"
                records.append(record)
                continue

            record.update(
                {
                    "evaluable": True,
                    "bb_supported": evidence["recommended_model"] == bb_model,
                    "rejection_reason": None,
                    "delta_cstat": lrt["delta_cstat"],
                    "delta_aic": info["delta_aic"],
                    "delta_bic": info["delta_bic"],
                    "lrt_supports_complex": evidence["lrt_supports_complex"],
                    "aic_supports_complex": evidence["aic_supports_complex"],
                    "bic_supports_complex": evidence["bic_supports_complex"],
                    "recommended_model": evidence["recommended_model"],
                    "selection_conflict": evidence["selection_conflict"],
                }
            )
            records.append(record)

    return records


def make_stratified_interval_folds(
    records: list[dict[str, Any]],
    n_splits: int = 5,
    seed: int = 123,
) -> list[list[int]]:
    """Make deterministic, approximately stratified folds over interval records."""
    if n_splits < 2:
        raise ValueError("n_splits must be at least 2")
    if len(records) < n_splits:
        raise ValueError("n_splits cannot exceed the number of records")

    rng = random.Random(seed)
    strata: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        strata[(record["grb"], record["interval_kind"])].append(index)

    folds: list[list[int]] = [[] for _ in range(n_splits)]
    for key in sorted(strata):
        indices = strata[key]
        rng.shuffle(indices)
        for index in indices:
            target = min(range(n_splits), key=lambda fold_index: len(folds[fold_index]))
            folds[target].append(index)

    for fold in folds:
        fold.sort()
    return folds


def make_leave_one_grb_out_folds(records: list[dict[str, Any]]) -> list[list[int]]:
    """Make one held-out fold per GRB."""
    grb_indices: dict[str, list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        grb_indices[record["grb"]].append(index)
    return [grb_indices[grb] for grb in sorted(grb_indices)]


def _rate(records: list[dict[str, Any]], key: str) -> float | None:
    values = [bool(record[key]) for record in records if record.get(key) is not None]
    if not values:
        return None
    return sum(values) / len(values)


def _mean(records: list[dict[str, Any]], key: str) -> float | None:
    values = [_numeric_or_none(record.get(key)) for record in records]
    values = [value for value in values if value is not None]
    return mean(values) if values else None


def _std(records: list[dict[str, Any]], key: str) -> float | None:
    values = [_numeric_or_none(record.get(key)) for record in records]
    values = [value for value in values if value is not None]
    return stdev(values) if len(values) > 1 else None


def _median(records: list[dict[str, Any]], key: str) -> float | None:
    values = [_numeric_or_none(record.get(key)) for record in records]
    values = [value for value in values if value is not None]
    return median(values) if values else None


def _iqr(records: list[dict[str, Any]], key: str) -> float | None:
    """Q3 - Q1 (linear/inclusive interpolation), robust to skew/outliers that
    distort mean/std on small, heavy-tailed evaluable subsets."""
    values = [_numeric_or_none(record.get(key)) for record in records]
    values = [value for value in values if value is not None]
    if len(values) < 2:
        return None
    q1, _, q3 = quantiles(values, n=4, method="inclusive")
    return q3 - q1


def _summarize_subset(records: list[dict[str, Any]]) -> dict[str, Any]:
    evaluable = [record for record in records if record["evaluable"]]
    rejected = [record for record in records if not record["evaluable"]]
    rejection_counts = Counter(record["rejection_reason"] for record in rejected)
    model_counts = Counter(record["simple_model"] for record in records if record["simple_model"])

    return {
        "n_records": len(records),
        "n_evaluable": len(evaluable),
        "n_rejected": len(rejected),
        "bb_support_rate": _rate(evaluable, "bb_supported"),
        "lrt_support_rate": _rate(evaluable, "lrt_supports_complex"),
        "aic_support_rate": _rate(evaluable, "aic_supports_complex"),
        "bic_support_rate": _rate(evaluable, "bic_supports_complex"),
        "selection_conflict_rate": _rate(evaluable, "selection_conflict"),
        "mean_delta_cstat": _mean(evaluable, "delta_cstat"),
        "std_delta_cstat": _std(evaluable, "delta_cstat"),
        "median_delta_cstat": _median(evaluable, "delta_cstat"),
        "iqr_delta_cstat": _iqr(evaluable, "delta_cstat"),
        "mean_delta_aic": _mean(evaluable, "delta_aic"),
        "std_delta_aic": _std(evaluable, "delta_aic"),
        "mean_delta_bic": _mean(evaluable, "delta_bic"),
        "std_delta_bic": _std(evaluable, "delta_bic"),
        "rejection_counts": dict(rejection_counts),
        "simple_model_counts": dict(model_counts),
    }


def summarize_cross_validation(
    records: list[dict[str, Any]],
    folds: list[list[int]],
    seed: int | None = None,
    fold_groups: list[str] | None = None,
) -> dict[str, Any]:
    """Summarize fold-level and aggregate BB-selection stability metrics."""
    fold_summaries = []
    for fold_number, indices in enumerate(folds, start=1):
        held_out = [records[index] for index in indices]
        train_count = len(records) - len(held_out)
        fold_summary = _summarize_subset(held_out)
        fold_summary.update(
            {
                "fold": fold_number,
                "train_count": train_count,
                "held_out_indices": indices,
            }
        )
        if fold_groups is not None:
            fold_summary["held_out_group"] = fold_groups[fold_number - 1]
        fold_summaries.append(fold_summary)

    evaluable_fold_rates = [
        fold["bb_support_rate"]
        for fold in fold_summaries
        if fold["bb_support_rate"] is not None
    ]

    min_evaluable = min((fold["n_evaluable"] for fold in fold_summaries), default=0)
    summary = {
        "n_splits": len(folds),
        "overall": _summarize_subset(records),
        "folds": fold_summaries,
        "fold_bb_support_rate_mean": mean(evaluable_fold_rates) if evaluable_fold_rates else None,
        "fold_bb_support_rate_std": stdev(evaluable_fold_rates) if len(evaluable_fold_rates) > 1 else None,
    }
    if seed is not None:
        summary["seed"] = seed
    if min_evaluable < 5:
        summary["variance_note"] = (
            "At least one held-out fold has fewer than 5 evaluable records; "
            "fold-rate variance is sensitive to small denominators."
        )
    return summary


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_safe(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_json_safe(item) for item in value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_cross_validation_outputs(
    output_dir: str | Path,
    summary: dict[str, Any],
    records: list[dict[str, Any]],
) -> tuple[Path, Path]:
    """Write JSON summary and per-interval CSV records."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    summary_path = output_path / "bb_selection_cv_summary.json"
    records_path = output_path / "bb_selection_cv_records.csv"
    summary["output_files"] = {
        "summary": str(summary_path),
        "records": str(records_path),
    }

    disk_summary = dict(summary)
    disk_summary.pop("records", None)
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(_json_safe(disk_summary), f, indent=4)

    with records_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=RECORD_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    return summary_path, records_path


def run_bb_selection_cross_validation(
    results_path: str | Path = "results.json",
    n_splits: int = 5,
    seed: int = 123,
    output_dir: str | Path | None = None,
    include_marginal: bool = False,
    grb_names: list[str] | tuple[str, ...] | None = DEFAULT_GRB_NAMES,
) -> dict[str, Any]:
    """Run interval k-fold plus leave-one-GRB-out BB-selection stability checks."""
    data = filter_results_by_grbs(load_results(results_path), grb_names)
    records = collect_bb_selection_records(data, include_marginal=include_marginal)
    folds = make_stratified_interval_folds(records, n_splits=n_splits, seed=seed)
    logo_folds = make_leave_one_grb_out_folds(records)
    logo_groups = []

    for fold_number, indices in enumerate(folds, start=1):
        for index in indices:
            records[index]["fold"] = fold_number

    for fold_number, indices in enumerate(logo_folds, start=1):
        held_out_grbs = sorted({records[index]["grb"] for index in indices})
        logo_groups.append(held_out_grbs[0])
        for index in indices:
            records[index]["logo_fold"] = fold_number

    summary = summarize_cross_validation(records, folds, seed=seed)
    summary["analysis"] = "interval_kfold"
    summary["grb_names"] = [long_to_short.get(name, name) for name in data]
    summary["grb_keys"] = list(data)
    summary["leave_one_grb_out"] = summarize_cross_validation(
        records,
        logo_folds,
        fold_groups=logo_groups,
    )
    summary["records"] = records

    if output_dir is not None:
        write_cross_validation_outputs(output_dir, summary, records)

    return _json_safe(summary)


def _format_rate(value: float | None) -> str:
    return "n/a" if value is None else f"{100.0 * value:.1f}%"


def _format_stat(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.2f}"


def print_summary(summary: dict[str, Any]) -> None:
    """Print a compact CLI summary."""
    overall = summary["overall"]
    logo = summary.get("leave_one_grb_out", {})
    logo_overall = logo.get("overall", {})
    print(f"BB-selection {summary['n_splits']}-fold interval cross-validation")
    print(
        f"records: {overall['n_records']} total, {overall['n_evaluable']} evaluable, {overall['n_rejected']} rejected")
    print(f"overall BB support rate: {_format_rate(overall['bb_support_rate'])}")
    print(
        f"delta_cstat: mean={_format_stat(overall['mean_delta_cstat'])} "
        f"std={_format_stat(overall['std_delta_cstat'])} "
        f"median={_format_stat(overall['median_delta_cstat'])} "
        f"iqr={_format_stat(overall['iqr_delta_cstat'])}"
    )
    print(f"mean fold BB support rate: {_format_rate(summary['fold_bb_support_rate_mean'])}")
    print(f"fold BB support std: {_format_rate(summary['fold_bb_support_rate_std'])}")
    if summary.get("variance_note"):
        print(f"note: {summary['variance_note']}")
    if logo:
        print("leave-one-GRB-out")
        print(f"  held-out GRBs: {logo['n_splits']}")
        print(f"  overall BB support rate: {_format_rate(logo_overall.get('bb_support_rate'))}")
        print(f"  mean held-out BB support rate: {_format_rate(logo.get('fold_bb_support_rate_mean'))}")
        print(f"  held-out BB support std: {_format_rate(logo.get('fold_bb_support_rate_std'))}")
        if logo.get("variance_note"):
            print(f"  note: {logo['variance_note']}")
    if overall["rejection_counts"]:
        print("rejections:")
        for reason, count in sorted(overall["rejection_counts"].items()):
            print(f"  {reason}: {count}")
    if summary.get("output_files"):
        print("outputs:")
        for label, path in summary["output_files"].items():
            print(f"  {label}: {path}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run BB-selection stability checks.")
    parser.add_argument("--results", default="results.json", help="Path to results.json")
    parser.add_argument("--folds", type=int, default=5, help="Number of interval-level k-fold folds")
    parser.add_argument("--seed", type=int, default=123, help="Deterministic fold seed")
    parser.add_argument("--output-dir", default=None, help="Directory for JSON/CSV outputs")
    parser.add_argument(
        "--grbs",
        nargs="+",
        default=DEFAULT_GRB_NAMES,
        help="Short or long GRB names to evaluate",
    )
    parser.add_argument(
        "--include-marginal",
        action="store_true",
        help="Treat MARGINAL models as eligible in addition to BEST/SAFE",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = run_bb_selection_cross_validation(
        results_path=args.results,
        n_splits=args.folds,
        seed=args.seed,
        output_dir=args.output_dir,
        include_marginal=args.include_marginal,
        grb_names=args.grbs,
    )
    print_summary(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
