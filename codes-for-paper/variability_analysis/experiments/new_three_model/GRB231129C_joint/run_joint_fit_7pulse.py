"""GRB231129C, joint 7-pulse test: adds the candidate hidden peak near t=2s (test_hidden_peak.py's
finding -- Delta-chi2=7.6, t_v pinned at the 2*dt resolution floor when fit in isolation against a
frozen residual) as a genuine 7th pulse, fit SIMULTANEOUSLY with the other 6 against the raw data.

Why this is a different, better test than test_hidden_peak.py's isolated fit: there, the candidate
was fit against a residual built by subtracting the FIXED 6-pulse model -- it never got a chance to
trade off against pulse 2 (t_peak=2.47s, right next door). Here, all 7 pulses' parameters are free
together, so if the candidate is actually just pulse 2's own rise being slightly the wrong shape
(a modeling artifact, not a separate feature), the joint fit can reveal that by pulse 2 absorbing
it instead of leaving a stable, independent 7th component. Resolution ceiling is real and already
confirmed (64ms bins = 0.064s = the finest available for this burst, per user) -- if the candidate
still wants to be sub-2-bin-wide here, that's the light curve's floor, not a fitting failure.

Reuses run_joint_fit.py's already-tested joint machinery directly (fit_joint, the three audit
checks, plot_joint_overlay) rather than reimplementing it.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_joint_fit import (  # noqa: E402
    EDGE_REL_TOL,
    WIDEN_FACTORS,
    WINDOW,
    check_edge_pinning,
    check_multistart_reproducibility,
    check_window_widening,
    fit_joint,
    load_prior_single_pulse_results,
)
from bounds_seeding import median_dt  # noqa: E402
from joint_pulse3 import flatten_seed, unflatten_params  # noqa: E402
from load_data import load_archived_pulses, load_light_curve  # noqa: E402
from plot_joint_fit import plot_joint_overlay  # noqa: E402
from profile_scan import chi_square  # noqa: E402
from grb_research import get_rng, seed_from_name  # noqa: E402

SEED = seed_from_name(__file__)
# Candidate pulse 7's seed, taken from test_hidden_peak.py's own isolated single-pulse result
# (r0_chosen, amplitude, t_peak, t_v) -- reused as a starting point, not re-derived here.
CANDIDATE_R0 = 31.622776601683793
CANDIDATE_SEED = (1049.29, 1.9050, 0.1500)  # A, t_peak, t_v -- nudged just above the isolated
# test's own t_v=0.1280 (exactly the 2*dt floor there), since starting a bounded optimizer
# exactly on its own boundary is numerically fragile (floating-point can put the seed a hair
# outside [lb, ub]); the joint fit is free to push it back down to the floor again on its own.


def main():
    print("=== GRB231129C: JOINT 7-pulse fit (6 archived pulses + candidate hidden peak near t=2s) ===")
    t_full, y_full, sigma_full = load_light_curve()
    archived = load_archived_pulses()
    prior = load_prior_single_pulse_results()

    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_window, y_window, sigma_window = t_full[mask], y_full[mask], sigma_full[mask]
    dt = median_dt(t_window)

    r0_list_6 = prior["r0_chosen"].tolist()
    p0_6 = list(zip(prior["amplitude"], prior["t_peak_s"], prior["t_v_s"]))
    r0_list_7 = r0_list_6 + [CANDIDATE_R0]
    p0_7 = p0_6 + [CANDIDATE_SEED]
    p0_flat_7 = flatten_seed(p0_7)
    pulse_index_list = prior["pulse_index"].tolist() + [7]
    n_pulses = 7

    print(f"pulses: {pulse_index_list} (pulse 7 = candidate, seed A={CANDIDATE_SEED[0]:.1f}, "
          f"t_peak={CANDIDATE_SEED[1]:.4f}, t_v={CANDIDATE_SEED[2]:.4f}, r0={CANDIDATE_R0:.3g})")

    print("\n--- Joint 7-pulse fit ---")
    model7, popt7, pcov7, chi2_7 = fit_joint(t_window, y_window, sigma_window, r0_list_7, p0_flat_7, dt)
    params7 = unflatten_params(popt7, n_pulses)
    print(f"chi2_7pulse (21 free params) = {chi2_7:.4f}")

    print("\n--- Compare against the already-run 6-pulse joint fit ---")
    prior_results = pd.read_csv(HERE / "GRB231129C_joint_results.csv")
    r0_list_6_check = prior_results["r0"].tolist()
    p0_flat_6 = flatten_seed(list(zip(prior_results["amplitude"], prior_results["t_peak_s"], prior_results["t_v_s"])))
    model6, popt6, pcov6, chi2_6 = fit_joint(t_window, y_window, sigma_window, r0_list_6_check, p0_flat_6, dt)
    print(f"chi2_6pulse (18 free params) = {chi2_6:.4f}")
    print(f"delta_chi2 (6pulse - 7pulse) = {chi2_6 - chi2_7:.4f} for 3 added params "
          f"(compare to the isolated single-pulse test's delta_chi2=7.592)")

    print("\n--- Pulse 7 (candidate) in the joint context ---")
    a7, tp7, tv7 = params7[-1]
    block7 = pcov7[18:21, 18:21]
    a7_err, tp7_err, tv7_err = np.sqrt(np.diag(block7))
    print(f"A={a7:.2f}+/-{a7_err:.2f}  t_peak={tp7:.4f}+/-{tp7_err:.4f}  t_v={tv7:.4f}+/-{tv7_err:.4f}")

    print("\n--- Pulse 2 and pulse 3, 6-pulse vs. 7-pulse (does the candidate come from their shape?) ---")
    for idx, pulse_num in [(1, 2), (2, 3)]:  # position 1 = pulse 2, position 2 = pulse 3 in the sorted list
        a6, tp6, tv6 = unflatten_params(popt6, 6)[idx]
        a7_, tp7_, tv7_ = params7[idx]
        print(f"  pulse {pulse_num}: 6-pulse A={a6:.1f} t_peak={tp6:.4f} t_v={tv6:.4f}  |  "
              f"7-pulse A={a7_:.1f} t_peak={tp7_:.4f} t_v={tv7_:.4f}  |  "
              f"delta A={a7_-a6:+.1f} delta t_peak={tp7_-tp6:+.5f} delta t_v={tv7_-tv6:+.5f}")

    print("\n--- Audit (joint analogs of section 8), 7-pulse fit ---")
    from joint_pulse3 import joint_bounds

    lb7, ub7 = joint_bounds(t_window, dt, n_pulses)
    pinned = check_edge_pinning(popt7, lb7, ub7, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")

    rng = get_rng(seed=SEED)
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(
        t_window, y_window, sigma_window, r0_list_7, p0_flat_7, dt, rng
    )
    print(f"multi-start reproducibility: {'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")
    print(f"  per-pulse spread (A, t_peak, t_v): pulse 7 (candidate) = {repro_spread[18:21]}")

    widen_ok, widen_results = check_window_widening(t_full, y_full, sigma_full, r0_list_7, popt7, dt)
    print(f"window-widening invariance (x1.5, x2.0): {'PASS' if widen_ok else 'FAIL'}")
    for factor, params in widen_results.items():
        a_w, tp_w, tv_w = params[-1]
        print(f"  x{factor}: pulse 7 -> A={a_w:.1f} t_peak={tp_w:.4f} t_v={tv_w:.4f} "
              f"(base: A={a7:.1f} t_peak={tp7:.4f} t_v={tv7:.4f})")

    audit_all_passed = (not pinned) and repro_ok and widen_ok
    print(f"\naudit overall: {'ALL PASS' if audit_all_passed else 'SOME CHECKS FAILED'}")

    print("\n--- Plot ---")
    paths = plot_joint_overlay(
        t_window, y_window, sigma_window, model7, popt7, out_dir=str(HERE), label="GRB231129C_joint_7pulse",
        pulse_indices=pulse_index_list,
    )
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print("\n=== Verdict ===")
    a7_sigma = a7 / a7_err if a7_err > 0 else float("nan")
    print(f"pulse 7 amplitude significance: {a7_sigma:.2f} sigma")
    tv_at_floor = np.isclose(tv7, 2 * dt, rtol=1e-3) or np.isclose(popt7[20], lb7[20], rtol=EDGE_REL_TOL)
    print(f"pulse 7 t_v at the 2*dt resolution floor ({2*dt:.4f}): {tv_at_floor}")


if __name__ == "__main__":
    main()
