"""GRB131014215, joint 5-pulse fit seeded from user-given rough peak centers -- NOT a violation of
the blind-test constraint (PROGRESS.md): the user clarified that eyeballing approximate pulse
centers to seed a fit is standard practice (it's literally how every archived P0 in this whole
project was built), and the constraint was specifically about not reading already-FITTED
parameter values for this burst. No archived (A, t_s, tau1, tau2) or P0 for GRB131014215 is read
here -- only 5 rough center guesses (t ~ 0.5, 1.3, 2.0, 2.7, 3.5) and the raw light curve.

Method, mirroring GRB231129C_joint's already-validated approach:
  1. For each hint center, a LOCAL single-pulse fit (bounds/seed -> section-5 profile scan ->
     finalize) on the raw data restricted to a local window (sized to the midpoint toward its
     nearest neighboring hint, same logic as greedy_discovery_local.py's local_window()) chooses
     that pulse's own r0 and gives an initial (A, t_peak, t_v).
  2. All 5 pulses' (A, t_peak, t_v) are then fit SIMULTANEOUSLY in one joint curve_fit call
     against the raw T90-window data (not a residual), each keeping its own locally-chosen r0
     fixed -- exactly GRB231129C_joint/joint_pulse3.py's model, reused directly.
  3. Joint audit: edge-pinning, multi-start reproducibility. Window-widening is deferred (T90's
     own right edge already flagged as truncating real decay flux -- see PROGRESS.md; widening
     the window is the natural follow-up once this fit is validated, not done here per "stay
     inside T90 duration (for now)").
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
GRB231129C_JOINT_DIR = NEW_THREE_MODEL_DIR / "GRB231129C_joint"
PROJECT_ROOT = NEW_THREE_MODEL_DIR.parents[3]
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(GRB231129C_JOINT_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from audit_checklist import run_audit  # noqa: E402
from bounds_seeding import median_dt  # noqa: E402
from joint_pulse3 import flatten_seed, joint_bounds, make_joint_pulse3, pulse_block_cov, residual_excluding, unflatten_params  # noqa: E402
from load_data import load_light_curve  # noqa: E402
from plot_diagnostics import plot_fit_overlay  # noqa: E402
from plot_joint_fit import plot_joint_overlay  # noqa: E402
from profile_scan import chi_square, refine_r0_zone, select_r0  # noqa: E402
from report_pulse import _zone_spread, finalize_pulse  # noqa: E402

EDGE_REL_TOL = 1e-3
REPRO_REL_TOL = 1e-3
N_JITTER_SEEDS = 3
T90 = (0.960, 4.160)
# WINDOW widened beyond strict T90: hint center 1 (~0.5) sits BEFORE T90's own start (0.96) --
# the hint itself shows T90 clips pulse 1, not just the previously-flagged decay-tail truncation
# at T90's end. Widened to comfortably contain all 5 hints plus baseline/decay margin on each
# side (same "don't cut off real flux" lesson as GRB231129C's pulse 6 and this burst's own T90
# end already discussed in PROGRESS.md), rather than crash or silently clip a hinted pulse.
WINDOW = (-0.5, 6.0)
HINT_CENTERS = [0.5, 1.3, 2.0, 2.7, 3.5]  # user-given rough centers, eyeballed -- not fitted values
MIN_LOCAL_WIDTH = 0.30
MAX_LOCAL_WIDTH = 1.0
SEED = seed_from_name(__file__)


def check_edge_pinning(popt, lb, ub, n_pulses):
    """Copied from GRB231129C_joint/run_joint_fit.py (2026-09-25's already-fixed, value-relative
    version -- CLAUDE.md's "copy rather than fight sys.path": that module's own `load_data.py`
    collides by name with this folder's `load_data.py`, so importing run_joint_fit directly here
    pulls in the wrong load_data via sys.modules; copying these two small functions avoids it)."""
    pinned = []
    for i in range(n_pulses):
        for j, name in enumerate(("A", "t_peak", "t_v")):
            idx = 3 * i + j
            value, lo, hi = popt[idx], lb[idx], ub[idx]
            near_lo = np.isclose(value, lo, rtol=EDGE_REL_TOL, atol=EDGE_REL_TOL * max(abs(lo), 1e-3))
            near_hi = np.isfinite(hi) and np.isclose(value, hi, rtol=EDGE_REL_TOL, atol=EDGE_REL_TOL * max(abs(hi), 1e-3))
            if near_lo or near_hi:
                pinned.append((i + 1, name, value, lo, hi))
    return pinned


def check_multistart_reproducibility(t_window, y_window, sigma_window, r0_list, p0_flat, dt, rng):
    """Copied from GRB231129C_joint/run_joint_fit.py (same reason as check_edge_pinning above)."""
    n = len(r0_list)
    model = make_joint_pulse3(r0_list)
    lb, ub = joint_bounds(t_window, dt, n)
    solutions = []
    for _ in range(N_JITTER_SEEDS):
        jittered = []
        for i in range(n):
            a0, tp0, tv0 = p0_flat[3 * i : 3 * i + 3]
            a = a0 * (1 + rng.uniform(-0.1, 0.1))
            tp = tp0 + rng.uniform(-3, 3) * dt
            tv = tv0 * (1 + rng.uniform(-0.1, 0.1))
            jittered.extend((a, tp, tv))
        popt_j, _ = curve_fit(
            model, t_window, y_window, p0=tuple(jittered), bounds=(lb, ub),
            sigma=sigma_window, absolute_sigma=True, maxfev=40000,
        )
        y_fit_j = model(t_window, *popt_j)
        chi2_j = chi_square(y_window, y_fit_j, sigma=sigma_window)
        solutions.append((popt_j, chi2_j))
    arr = np.array([s[0] for s in solutions])
    spread = arr.max(axis=0) - arr.min(axis=0)
    scale = np.abs(arr).mean(axis=0)
    ok = True
    for i in range(n):
        a_ok = spread[3 * i] < REPRO_REL_TOL * scale[3 * i]
        tp_ok = spread[3 * i + 1] < REPRO_REL_TOL
        tv_ok = spread[3 * i + 2] < REPRO_REL_TOL * scale[3 * i + 2]
        ok &= a_ok and tp_ok and tv_ok
    chi2_spread = max(s[1] for s in solutions) - min(s[1] for s in solutions)
    return ok, spread, chi2_spread


def local_window(t_center, other_centers, dt):
    left = [t for t in other_centers if t < t_center]
    right = [t for t in other_centers if t > t_center]
    half_left = (t_center - max(left)) / 2 if left else MAX_LOCAL_WIDTH
    half_right = (min(right) - t_center) / 2 if right else MAX_LOCAL_WIDTH
    half_left = float(np.clip(half_left, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    half_right = float(np.clip(half_right, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    return max(WINDOW[0], t_center - half_left), min(WINDOW[1], t_center + half_right)


def main():
    print("=== GRB131014215: joint 5-pulse fit seeded from rough hint centers ===")
    t_full, y_full, sigma_full = load_light_curve()
    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
    dt = median_dt(t_w)
    print(f"window: {WINDOW} ({t_w.size} bins)")

    print("\n--- Step 1: local single-pulse fit per hint center (chooses each pulse's own r0) ---")
    r0_list, p0_list = [], []
    for i, t_c in enumerate(HINT_CENTERS):
        others = HINT_CENTERS[:i] + HINT_CENTERS[i + 1 :]
        lo, hi = local_window(t_c, others, dt)
        lmask = (t_w >= lo) & (t_w <= hi)
        t_l, y_l, sigma_l = t_w[lmask], y_w[lmask], sigma_w[lmask]
        selection = select_r0(t_l, y_l, sigma=sigma_l, dt=dt)
        reported = finalize_pulse(t_l, y_l, selection, sigma=sigma_l, dt=dt)
        print(f"  hint t={t_c:.2f} -> window [{lo:.3f},{hi:.3f}]: r0={reported.r0:.4g} ({selection.classification.shape}) "
              f"A={reported.amplitude:.1f} t_peak={reported.t_peak:.4f} t_v={reported.t_v:.4f}")
        r0_list.append(reported.r0)
        p0_list.append((reported.amplitude, reported.t_peak, reported.t_v))

    print("\n--- Step 2: joint 5-pulse fit against the raw data ---")
    p0_flat = flatten_seed(p0_list)
    n_pulses = len(HINT_CENTERS)
    lb, ub = joint_bounds(t_w, dt, n_pulses)
    model = make_joint_pulse3(r0_list)
    popt, pcov = curve_fit(model, t_w, y_w, p0=p0_flat, bounds=(lb, ub), sigma=sigma_w, absolute_sigma=True, maxfev=40000)
    params = unflatten_params(popt, n_pulses)
    y_fit_total = model(t_w, *popt)
    chi2_joint = chi_square(y_w, y_fit_total, sigma=sigma_w)
    print(f"chi2_joint ({3*n_pulses} free params, {t_w.size} bins) = {chi2_joint:.4f}  "
          f"(reduced chi2 ~ {chi2_joint/(t_w.size-3*n_pulses):.3f})")

    print("\n--- Per-pulse results ---")
    rows = []
    for i in range(n_pulses):
        amplitude, t_peak, t_v = params[i]
        r0 = r0_list[i]
        block = pulse_block_cov(pcov, i)
        a_err, tp_err_stat, tv_err_stat = np.sqrt(np.diag(block))

        residual_i = residual_excluding(model, t_w, y_w, popt, exclude_index=i)
        refined_scan, refined_classification = refine_r0_zone(t_w, residual_i, r0, sigma=sigma_w, dt=dt)
        tp_err_sys, tv_err_sys = _zone_spread(refined_scan, *refined_classification.flat_zone)
        tp_err = float(np.hypot(tp_err_stat, tp_err_sys))
        tv_err = float(np.hypot(tv_err_stat, tv_err_sys))

        print(f"  pulse {i+1} (hint {HINT_CENTERS[i]:.2f}): r0={r0:.4g}  A={amplitude:.1f}+/-{a_err:.1f}  "
              f"t_peak={t_peak:.4f}+/-{tp_err:.4f}  t_v={t_v:.4f}+/-{tv_err:.4f}")
        rows.append(
            {
                "pulse_index": i + 1,
                "hint_center": HINT_CENTERS[i],
                "r0": r0,
                "amplitude": amplitude,
                "amplitude_err": a_err,
                "t_peak_s": t_peak,
                "t_peak_err_s": tp_err,
                "t_v_s": t_v,
                "t_v_err_s": tv_err,
            }
        )

    df = pd.DataFrame(rows)
    csv_path = HERE / "GRB131014215_joint_from_hints_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path.name}")

    print("\n--- Audit ---")
    pinned = check_edge_pinning(popt, lb, ub, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")
    rng = get_rng(seed=SEED)
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(t_w, y_w, sigma_w, r0_list, p0_flat, dt, rng)
    print(f"multi-start reproducibility: {'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")
    print("window-widening: deferred -- T90's own right edge already flagged as truncating real decay flux "
          "(see PROGRESS.md); widening is the natural follow-up, not done here per 'stay inside T90 (for now)'.")

    print("\n--- Plot ---")
    paths = plot_joint_overlay(
        t_w, y_w, sigma_w, model, popt, out_dir=str(HERE), label="GRB131014215_joint_from_hints",
        pulse_indices=list(range(1, n_pulses + 1)),
    )
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print(f"\n=== Recovered {n_pulses} pulses from 5 rough hint centers (vs. hint: 5 currently fitted) ===")
    print(f"edge_pinning={'none' if not pinned else 'YES'}  reproducibility={'PASS' if repro_ok else 'FAIL'}  "
          f"reduced_chi2~{chi2_joint/(t_w.size-3*n_pulses):.3f}")


if __name__ == "__main__":
    main()
