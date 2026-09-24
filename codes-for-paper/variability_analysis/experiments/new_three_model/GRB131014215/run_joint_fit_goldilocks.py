"""GRB131014215, joint 5-pulse fit over a "goldilocks" window: wide enough to contain the real
decay tail, not so wide it hands a tail-heavy Norris pulse room to wander into pure noise
(PROGRESS.md's "Full-span test" finding -- the full 478s span made reproducibility WORSE, not
better, because this burst's baseline beyond ~t=10-11s isn't flat zero-mean noise the way
GRB231129C's was). WINDOW = (-1, 11)s: comfortably contains T90 (0.96-4.16), the 5 hint centers,
and the directly-measured real decay tail (mean significance ~1-3 sigma per second-wide bin out to
t=10-11s, checked earlier in this session), stopping before the region that's plausibly pure noise.
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

from bounds_seeding import median_dt  # noqa: E402
from joint_pulse3 import flatten_seed, joint_bounds, make_joint_pulse3, pulse_block_cov, residual_excluding, unflatten_params  # noqa: E402
from load_data import load_light_curve  # noqa: E402
from plot_joint_fit import plot_joint_overlay  # noqa: E402
from profile_scan import chi_square, refine_r0_zone, select_r0  # noqa: E402
from report_pulse import _zone_spread, finalize_pulse  # noqa: E402

WINDOW = (-1.0, 11.0)
T90 = (0.960, 4.160)
PLOT_PAD_S = 2.0
HINT_CENTERS = [0.5, 1.3, 2.0, 2.7, 3.5]
MIN_LOCAL_WIDTH = 0.30
MAX_LOCAL_WIDTH = 1.0
LOCAL_HINT_WINDOW = (-0.5, 6.0)  # unchanged from joint_fit_from_hints.py -- per-pulse r0 choice
# doesn't need the tail; only the JOINT fit's own window changes here.
EDGE_REL_TOL = 1e-3
REPRO_REL_TOL = 1e-3
N_JITTER_SEEDS = 3
SEED = seed_from_name(__file__)


def check_edge_pinning(popt, lb, ub, n_pulses):
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
        ok &= spread[3 * i] < REPRO_REL_TOL * scale[3 * i]
        ok &= spread[3 * i + 1] < REPRO_REL_TOL
        ok &= spread[3 * i + 2] < REPRO_REL_TOL * scale[3 * i + 2]
    chi2_spread = max(s[1] for s in solutions) - min(s[1] for s in solutions)
    return ok, spread, chi2_spread


def local_window(t_center, other_centers, dt):
    left = [t for t in other_centers if t < t_center]
    right = [t for t in other_centers if t > t_center]
    half_left = (t_center - max(left)) / 2 if left else MAX_LOCAL_WIDTH
    half_right = (min(right) - t_center) / 2 if right else MAX_LOCAL_WIDTH
    half_left = float(np.clip(half_left, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    half_right = float(np.clip(half_right, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    return max(LOCAL_HINT_WINDOW[0], t_center - half_left), min(LOCAL_HINT_WINDOW[1], t_center + half_right)


def main():
    print(f"=== GRB131014215: JOINT 5-pulse fit, goldilocks window {WINDOW} ===")
    t_full, y_full, sigma_full = load_light_curve()
    dt = median_dt(t_full)
    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
    print(f"window: {WINDOW} ({t_w.size} bins)")

    print("\n--- Step 1: local single-pulse fit per hint center (unchanged from joint_fit_from_hints.py) ---")
    mask_local = (t_full >= LOCAL_HINT_WINDOW[0]) & (t_full <= LOCAL_HINT_WINDOW[1])
    t_lw, y_lw, sigma_lw = t_full[mask_local], y_full[mask_local], sigma_full[mask_local]
    r0_list, p0_list = [], []
    for i, t_c in enumerate(HINT_CENTERS):
        others = HINT_CENTERS[:i] + HINT_CENTERS[i + 1 :]
        lo, hi = local_window(t_c, others, dt)
        lmask = (t_lw >= lo) & (t_lw <= hi)
        t_l, y_l, sigma_l = t_lw[lmask], y_lw[lmask], sigma_lw[lmask]
        selection = select_r0(t_l, y_l, sigma=sigma_l, dt=dt)
        reported = finalize_pulse(t_l, y_l, selection, sigma=sigma_l, dt=dt)
        r0_list.append(reported.r0)
        p0_list.append((reported.amplitude, reported.t_peak, reported.t_v))
    n_pulses = len(HINT_CENTERS)
    p0_flat = flatten_seed(p0_list)

    print("\n--- Step 2: joint fit over the goldilocks window ---")
    lb, ub = joint_bounds(t_w, dt, n_pulses)
    model = make_joint_pulse3(r0_list)
    popt, pcov = curve_fit(model, t_w, y_w, p0=p0_flat, bounds=(lb, ub), sigma=sigma_w, absolute_sigma=True, maxfev=60000)
    params = unflatten_params(popt, n_pulses)
    y_fit_total = model(t_w, *popt)
    chi2_w = chi_square(y_w, y_fit_total, sigma=sigma_w)
    print(f"chi2 = {chi2_w:.4f}  reduced ~ {chi2_w/(t_w.size-3*n_pulses):.3f}")

    print("\n--- Compare across all three windows tried ---")
    windowed = pd.read_csv(HERE / "GRB131014215_joint_from_hints_results.csv")
    fullspan = pd.read_csv(HERE / "GRB131014215_joint_fullspan_results.csv")
    print(f"{'pulse':>6} {'t_peak(-0.5,6)':>15} {'t_peak(full)':>13} {'t_peak(goldi)':>14}   "
          f"{'t_v(-0.5,6)':>12} {'t_v(full)':>10} {'t_v(goldi)':>11}")
    for i in range(n_pulses):
        rw, rf = windowed.iloc[i], fullspan.iloc[i]
        a, tp, tv = params[i]
        print(f"{i+1:>6} {rw['t_peak_s']:>15.4f} {rf['t_peak_s']:>13.4f} {tp:>14.4f}   "
              f"{rw['t_v_s']:>12.4f} {rf['t_v_s']:>10.4f} {tv:>11.4f}")

    print("\n--- Per-pulse results (goldilocks) ---")
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
        print(f"  pulse {i+1}: A={amplitude:.1f}+/-{a_err:.1f}  t_peak={t_peak:.4f}+/-{tp_err:.4f}  t_v={t_v:.4f}+/-{tv_err:.4f}  r0={r0:.4g}")
        rows.append({"pulse_index": i + 1, "amplitude": amplitude, "amplitude_err": a_err, "t_peak_s": t_peak,
                      "t_peak_err_s": tp_err, "t_v_s": t_v, "t_v_err_s": tv_err, "r0": r0})

    df = pd.DataFrame(rows)
    csv_path = HERE / "GRB131014215_joint_goldilocks_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path.name}")

    print("\n--- Audit ---")
    pinned = check_edge_pinning(popt, lb, ub, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")
    rng = get_rng(seed=SEED)
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(t_w, y_w, sigma_w, r0_list, p0_flat, dt, rng)
    print(f"multi-start reproducibility: {'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")

    print("\n--- Plot (fit: goldilocks window; display boxed to T90 +/- 2s) ---")
    plot_lo, plot_hi = T90[0] - PLOT_PAD_S, T90[1] + PLOT_PAD_S
    box_mask = (t_full >= plot_lo) & (t_full <= plot_hi)
    paths = plot_joint_overlay(
        t_full[box_mask], y_full[box_mask], sigma_full[box_mask], model, popt, out_dir=str(HERE),
        label="GRB131014215_joint_goldilocks", pulse_indices=list(range(1, n_pulses + 1)),
    )
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print(f"\n=== Verdict: reproducibility {'PASS' if repro_ok else 'STILL FAILING'} with the goldilocks window ===")


if __name__ == "__main__":
    main()
