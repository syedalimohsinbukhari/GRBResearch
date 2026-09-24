"""Test for a candidate hidden peak near t=2s in GRB231129C, reported by the user: found only
once or twice across a ~2-hour multi-start search on a separate machine, with no known initial
parameters. That's exactly what this folder's single-pulse pipeline doesn't need -- bounds_seeding's
neutral seed is argmax-based, not P0-based, so it can probe a candidate pulse with zero prior
guess at its own parameters.

Method: take the 6-pulse joint fit's own residual (data minus the current best 6-pulse model) as
the "y" to search, in a window around the excess actually seen there (checked first, directly
against the light curve: two adjacent bins at t=1.856s and t=1.984s sit at +2.15 and +3.43 sigma
above the 6-pulse model, immediately followed by a compensating dip as pulse 2's rise gets pulled
to partly absorb it -- adjacent-bin correlated excess, not an isolated one-bin outlier). Then run
this folder's full single-pulse pipeline (bounds/seed -> profile scan -> finalize -> audit) on
that residual, exactly as if it were pulse 7, and report honestly whether a real feature emerges.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
GRB231129C_DIR = NEW_THREE_MODEL_DIR / "GRB231129C"
PROJECT_ROOT = NEW_THREE_MODEL_DIR.parents[3]
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(GRB231129C_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from bounds_seeding import median_dt, pulse3_bounds_and_seed  # noqa: E402
from joint_pulse3 import flatten_seed, make_joint_pulse3  # noqa: E402
from load_data import load_light_curve  # noqa: E402
from plot_diagnostics import plot_fit_overlay  # noqa: E402
from profile_scan import chi_square, select_r0  # noqa: E402
from report_pulse import finalize_pulse  # noqa: E402
from audit_checklist import run_audit  # noqa: E402
from pulse3 import make_pulse3  # noqa: E402
from grb_research import get_rng, seed_from_name  # noqa: E402

JOINT_WINDOW = (-2.0, 20.0)
CANDIDATE_WINDOW = (0.5, 2.3)  # brackets the excess (t=1.86-1.98s) without reaching into pulse 2's
# own steep rise/peak region (t_peak=2.47s), where residual structure is about pulse 2's own shape
# mismatch rather than a separate feature.
SEED = seed_from_name(__file__)


def main():
    print("=== Candidate hidden peak near t=2s: single-pulse test on the 6-pulse residual ===")
    t_full, y_full, sigma_full = load_light_curve()
    mask = (t_full >= JOINT_WINDOW[0]) & (t_full <= JOINT_WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]

    df = pd.read_csv(HERE / "GRB231129C_joint_results.csv")
    r0_list = df["r0"].tolist()
    p0_flat = flatten_seed(list(zip(df["amplitude"], df["t_peak_s"], df["t_v_s"])))
    model_6pulse = make_joint_pulse3(r0_list)
    y_fit_6pulse = model_6pulse(t_w, *p0_flat)
    residual_full = y_w - y_fit_6pulse

    cmask = (t_w >= CANDIDATE_WINDOW[0]) & (t_w <= CANDIDATE_WINDOW[1])
    t_c, y_c, sigma_c = t_w[cmask], residual_full[cmask], sigma_w[cmask]
    dt = median_dt(t_c)
    print(f"candidate window: {CANDIDATE_WINDOW}, {len(t_c)} bins, dt={dt:.4f}")

    lb, ub, seed = pulse3_bounds_and_seed(t_c, y_c, dt=dt)
    print(f"neutral seed (argmax-based, NO manual P0): A0={seed[0]:.2f}, t_peak0={seed[1]:.4f}, t_v0={seed[2]:.4f}")
    print(f"bounds: A in [{lb[0]:.2f},{ub[0]:.2f}], t_peak in [{lb[1]:.4f},{ub[1]:.4f}], t_v in [{lb[2]:.4f},{ub[2]:.4f}]")

    print("\n--- chi2 with vs. without a candidate pulse (null test) ---")
    chi2_null = chi_square(y_c, np.zeros_like(y_c), sigma=sigma_c)
    print(f"chi2_null (no pulse, residual should be pure noise) = {chi2_null:.3f} over {len(y_c)} bins "
          f"(expect ~{len(y_c)} if residual really is just noise)")

    print("\n--- Section 5: chi2 profile scan over r0 ---")
    selection = select_r0(t_c, y_c, sigma=sigma_c, dt=dt)
    c = selection.classification
    d = selection.delta_chi2
    print(f"shape={c.shape}  r0_at_min={c.r0_at_min:.4g}  flat_zone={c.flat_zone}  span=x{c.flat_zone_span:.2f}")
    print(f"chi2_min_3param={c.chi2_min:.3f}  delta_chi2 vs 4param={d['delta_chi2']}  verdict={d['verdict']}")
    print(f"r0_chosen={selection.r0_chosen:.4g} ({selection.r0_choice_reason})")

    chi2_with_pulse = c.chi2_min
    print(f"\nchi2_null (no pulse) = {chi2_null:.3f}")
    print(f"chi2_with_candidate_pulse (best r0) = {chi2_with_pulse:.3f}")
    print(f"delta_chi2 (null - with_pulse) = {chi2_null - chi2_with_pulse:.3f} "
          f"(a real pulse costs 3 parameters; >~7-10 would usually justify adding them)")

    print("\n--- Section 4: finalize ---")
    reported = finalize_pulse(t_c, y_c, selection, sigma=sigma_c, dt=dt)
    print(f"A={reported.amplitude:.2f}+/-{reported.amplitude_err:.2f}")
    print(f"t_peak={reported.t_peak:.4f}+/-{reported.t_peak_err:.4f} (stat={reported.t_peak_err_stat:.4f}, sys={reported.t_peak_err_sys:.4f})")
    print(f"t_v={reported.t_v:.4f}+/-{reported.t_v_err:.4f} (stat={reported.t_v_err_stat:.4f}, sys={reported.t_v_err_sys:.4f})")
    print(f"r0_flat_zone={reported.r0_flat_zone}  r0_refined_zone={reported.r0_refined_zone}")

    print("\n--- Section 8: audit ---")
    rng = get_rng(seed=SEED)
    report = run_audit(t_c, y_c, selection, reported, rng, sigma=sigma_c, dt=dt)
    print(str(report))
    print(f"\naudit overall: {'ALL PASS' if report.all_passed else 'SOME CHECKS FAILED'}")

    model3 = make_pulse3(reported.r0)
    y_pulse_fit = model3(t_c, reported.amplitude, reported.t_peak, reported.t_v)
    paths = plot_fit_overlay(
        t_c, y_c, y_pulse_fit, out_dir=str(HERE), label="candidate_hidden_peak", sigma=sigma_c,
        extra_title=f"shape={c.shape}, r0={reported.r0:.3g}, t_peak={reported.t_peak:.3f}+/-{reported.t_peak_err:.3f}, "
                    f"A={reported.amplitude:.1f}+/-{reported.amplitude_err:.1f}",
    )
    print(f"\nwrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print("\n=== Verdict ===")
    a_sigma = reported.amplitude / reported.amplitude_err if reported.amplitude_err > 0 else float("nan")
    print(f"amplitude significance (A / A_err): {a_sigma:.2f} sigma")
    if a_sigma > 3 and (chi2_null - chi2_with_pulse) > 7:
        print("-> looks like a real, if weak, feature: positive amplitude at >3 sigma, meaningful chi2 improvement.")
    elif a_sigma > 2:
        print("-> marginal: positive and roughly consistent with a real feature, but not strong enough to call")
        print("   a confident detection on this data alone. Matches 'only converged once or twice' -- weak signal,")
        print("   easy for a generic multi-start search to miss or land on a spurious alternative instead.")
    else:
        print("-> not supported by this data: consistent with noise once a proper bounded fit is done.")


if __name__ == "__main__":
    main()
