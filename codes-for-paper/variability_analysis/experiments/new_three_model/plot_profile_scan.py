"""Reporting for the section-5 chi^2 profile scan (profile_scan.py): saves the per-r0 grid and
the decision-rule summary as CSV, and the chi2-vs-r0 profile as PDF+PNG -- this project's
generated-output convention (CSV via df.to_csv, figures via grb_research.update_style() +
grb_research.grb_utils.save_fig, both written into the calling script's own folder; see e.g.
experiments/window_sensitivity_GRB140206275/window_sensitivity_common.py).

Split out from profile_scan.py (which stays self-contained, no matplotlib/grb_research import)
so the core algorithm doesn't carry plotting dependencies.
"""
import sys
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[3]  # new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import update_style  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

from profile_scan import FLAT_ZONE_DELTA_CHI2, R0SelectionResult  # noqa: E402

update_style()


def profile_scan_dataframe(result: R0SelectionResult) -> pd.DataFrame:
    """One row per r0 grid point: the raw scan the spec's plot/decision rule is built from."""
    scan = result.scan
    zone_lo, zone_hi = result.classification.flat_zone
    return pd.DataFrame(
        {
            "r0": scan.r0_grid,
            "success": [r.success for r in scan.results],
            "chi2": scan.chi2,
            "amplitude": [r.amplitude for r in scan.results],
            "t_peak_s": scan.t_peak,
            "t_v_s": scan.t_v,
            "in_flat_zone": (scan.r0_grid >= zone_lo) & (scan.r0_grid <= zone_hi),
        }
    )


def profile_scan_summary_dataframe(result: R0SelectionResult, label: str) -> pd.DataFrame:
    """One row: the section-5 decision-rule outcome for this pulse (shape, chosen r0, Delta-chi2)."""
    c = result.classification
    d = result.delta_chi2
    row = {
        "label": label,
        "shape": c.shape,
        "r0_at_min": c.r0_at_min,
        "chi2_min_3param": c.chi2_min,
        "flat_zone_lo": c.flat_zone[0],
        "flat_zone_hi": c.flat_zone[1],
        "flat_zone_span": c.flat_zone_span,
        "r0_chosen": result.r0_chosen,
        "r0_choice_reason": result.r0_choice_reason,
        "chi2_4param_best": d["chi2_4param"],
        "delta_chi2": d["delta_chi2"],
        "delta_chi2_verdict": d["verdict"],
    }
    if result.fit4 is not None:
        row.update(
            {
                "fit4_amplitude": result.fit4["amplitude"],
                "fit4_t_peak_s": result.fit4["t_peak"],
                "fit4_t_v_s": result.fit4["t_v"],
                "fit4_r": result.fit4["r"],
            }
        )
    return pd.DataFrame([row])


def save_profile_scan_outputs(result: R0SelectionResult, out_dir=HERE, label: str = "pulse") -> dict:
    """Write r0_profile_scan_<label>.csv (grid), r0_profile_scan_summary_<label>.csv (decision),
    and r0_profile_scan_<label>.pdf/.png (chi2 vs r0, log x-axis, per spec section 5 step 3).
    Returns a dict of the written paths.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    grid_csv = out_dir / f"r0_profile_scan_{label}.csv"
    profile_scan_dataframe(result).to_csv(grid_csv, index=False)

    summary_csv = out_dir / f"r0_profile_scan_summary_{label}.csv"
    profile_scan_summary_dataframe(result, label).to_csv(summary_csv, index=False)

    fig_stem = out_dir / f"r0_profile_scan_{label}"
    _plot(result, label, fig_stem)

    return {
        "grid_csv": grid_csv,
        "summary_csv": summary_csv,
        "fig_pdf": Path(f"{fig_stem}.pdf"),
        "fig_png": Path(f"{fig_stem}.png"),
    }


def _plot(result: R0SelectionResult, label: str, fig_stem: Path) -> None:
    scan = result.scan
    c = result.classification

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(scan.r0_grid, scan.chi2, marker="o", ms=3, color="tab:blue", label=r"$\chi^2_{3\mathrm{param}}(r_0)$")
    ax.axhline(
        c.chi2_min + FLAT_ZONE_DELTA_CHI2,
        color="grey",
        ls="--",
        lw=0.8,
        label=rf"$\chi^2_\mathrm{{min}} + {FLAT_ZONE_DELTA_CHI2:.0f}$",
    )
    ax.axvspan(
        c.flat_zone[0],
        c.flat_zone[1],
        color="tab:blue",
        alpha=0.12,
        label="flat zone" if c.shape == "flat" else "near-min zone",
    )
    ax.axvline(result.r0_chosen, color="tab:red", ls=":", lw=1.2, label=rf"$r_0$ chosen = {result.r0_chosen:.3g}")
    if result.fit4 is not None:
        ax.axhline(
            result.fit4["chi2"],
            color="tab:green",
            ls="-.",
            lw=0.8,
            label=rf"best 4-param $\chi^2$ = {result.fit4['chi2']:.3g}",
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"$r_0 = \tau_1/\tau_2$")
    ax.set_ylabel(r"$\chi^2$")
    ax.set_title(f"{label}: {c.shape}-bottom profile (span x{c.flat_zone_span:.1f})")
    ax.legend(fontsize=7, loc="best")

    save_fig(fig, fig_stem)
