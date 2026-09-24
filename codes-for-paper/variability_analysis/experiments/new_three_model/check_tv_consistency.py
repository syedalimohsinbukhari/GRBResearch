"""Section-6 consistency check (norris_3param_spec.md): validate that the pipeline's existing
`tv_value()` (../../norris_fit.py) matches the width function the 3-param spec's f(r) is built to
invert, BEFORE any 3-param wrapper code is written. If this script fails, f(r) in the spec must be
rebuilt as the exact inverse of whatever `tv_value()` actually computes, or every 3-param result will
be subtly biased in (t_peak, t_v)-space (spec section 6).

Checks:
  1. Single worked example from the spec: tau1=85.6, tau2=0.453 -> t_v=1.40710035884524 (also spec T3).
  2. `tv_value()` (imported, the actual pipeline function) reproduces the spec's formula
         t_v = (tau2/2) * sqrt( (ln2 + 2*sqrt(tau1/tau2))^2 - 4*tau1/tau2 )
     to machine precision, evaluated independently (not by reading tv_value's source and copying it).
  3. `tv_value()` evaluated at the archived point-estimate (tau1, tau2) is consistent with the archived
     t_v_s column in every norris_fit_results*.csv under codes-for-paper/variability_analysis/.
     NOTE: archived t_v_s is NOT tv_value(tau1, tau2) -- it is the MC-propagated *median* over the
     fit covariance (see tv_mc_summary() in ../../norris_fit.py), which differs from the point estimate
     by O(1 MC-sigma) for skewed/degenerate pulses (the archived t_v_err_lower_s/upper_s columns are
     that MC spread). So this check passes if the point estimate falls within N sigma of the archived
     median, using the archived asymmetric error bars, not by exact equality.

This script does not modify anything outside this folder. It imports tv_value from norris_fit.py
read-only (sys.path insert, matching the convention in ../normalized_vs_unnormalized_fit/fit_comparison.py).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VARIABILITY_DIR = HERE.parent.parent  # codes-for-paper/variability_analysis/
PROJECT_ROOT = VARIABILITY_DIR.parent.parent  # GRBResearchWork/
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from norris_fit import tv_value  # noqa: E402
from grb_research import get_rng, seed_from_name  # noqa: E402

LN2 = np.log(2)
SEED = seed_from_name(__file__)  # project RNG convention -- src/grb_research/SEEDING.md


def spec_tv_value(tau1, tau2):
    """t_v formula transcribed independently from norris_3param_spec.md section 1/6."""
    tau1 = np.asarray(tau1, dtype=float)
    tau2 = np.asarray(tau2, dtype=float)
    r = tau1 / tau2
    return (tau2 / 2) * np.sqrt((LN2 + 2 * np.sqrt(r)) ** 2 - 4 * r)


def check_worked_example():
    tau1, tau2 = 85.6, 0.453
    expected = 1.40710035884524
    got_pipeline = float(tv_value(tau1, tau2))
    got_spec = float(spec_tv_value(tau1, tau2))
    print(f"[1] worked example tau1={tau1}, tau2={tau2}")
    print(f"    expected      = {expected!r}")
    print(f"    tv_value()    = {got_pipeline!r}  (diff {got_pipeline - expected:.3e})")
    print(f"    spec formula  = {got_spec!r}  (diff {got_spec - expected:.3e})")
    ok = np.isclose(got_pipeline, expected, atol=1e-6) and np.isclose(got_spec, expected, atol=1e-6)
    print(f"    -> {'PASS' if ok else 'FAIL'}")
    return ok


def check_pipeline_matches_spec_formula():
    rng = get_rng(seed=SEED)
    tau1 = rng.uniform(1e-3, 200.0, size=2000)
    tau2 = rng.uniform(1e-3, 200.0, size=2000)
    got_pipeline = tv_value(tau1, tau2)
    got_spec = spec_tv_value(tau1, tau2)
    max_abs_diff = np.max(np.abs(got_pipeline - got_spec))
    max_rel_diff = np.max(np.abs(got_pipeline - got_spec) / np.abs(got_spec))
    print(f"[2] tv_value() vs spec formula, 2000 random (tau1, tau2) in [1e-3, 200]")
    print(f"    max abs diff = {max_abs_diff:.3e}, max rel diff = {max_rel_diff:.3e}")
    ok = max_abs_diff < 1e-9
    print(f"    -> {'PASS' if ok else 'FAIL'}")
    return ok


def check_archived_fits(sigma_tol: float = 5.0):
    csv_paths = sorted(VARIABILITY_DIR.rglob("norris_fit_results*.csv"))
    print(f"[3] archived fit-result CSVs found: {len(csv_paths)}  (pass tol: within {sigma_tol} MC-sigma)")
    all_ok = True
    total_rows = 0
    worst_sigma = 0.0
    required = {"tau1", "tau2", "t_v_s", "t_v_err_lower_s", "t_v_err_upper_s"}
    for path in csv_paths:
        df = pd.read_csv(path)
        if not required.issubset(df.columns):
            print(f"    {path.relative_to(PROJECT_ROOT)}: missing required columns, skipped")
            continue
        point = tv_value(df["tau1"].to_numpy(), df["tau2"].to_numpy())
        archived = df["t_v_s"].to_numpy()
        diff = point - archived
        err = np.where(diff > 0, df["t_v_err_upper_s"].to_numpy(), df["t_v_err_lower_s"].to_numpy())
        n_sigma = np.abs(diff) / err
        n_bad = int(np.sum(n_sigma > sigma_tol))
        total_rows += len(df)
        worst_sigma = max(worst_sigma, float(np.max(n_sigma)))
        status = "PASS" if n_bad == 0 else "FAIL"
        print(f"    {path.relative_to(PROJECT_ROOT)}: {len(df)} rows, max |diff|/sigma = "
              f"{n_sigma.max():.3f} -> {status}")
        all_ok &= n_bad == 0
    print(f"    total rows checked: {total_rows}, worst-case deviation: {worst_sigma:.3f} sigma")
    print("    (deviation is expected: t_v_s is the MC-median of tv_value() over the fit covariance,")
    print("     not tv_value() at the point estimate -- see tv_mc_summary() in ../../norris_fit.py;")
    print("     the largest deviations should coincide with the lowest mc_kept_fraction rows.)")
    return all_ok


if __name__ == "__main__":
    results = [check_worked_example(), check_pipeline_matches_spec_formula(), check_archived_fits()]
    print()
    if all(results):
        print("ALL CHECKS PASSED: tv_value() matches the spec's width function. "
              "f(r) in section 1/2 of norris_3param_spec.md may be used as-is.")
    else:
        print("CHECKS FAILED: tv_value() diverges from the spec formula somewhere above. "
              "f(r) must be rebuilt as the exact inverse of the pipeline's actual tv_value() "
              "before writing the 3-param wrapper (spec section 6).")
        sys.exit(1)
