"""GRB131014215, joint 5-pulse fit with ITERATED r0: fixes the real lever identified in
PROGRESS.md's "Goldilocks window: also failed, differently" entry -- three different window sizes
all failed multi-start reproducibility, monotonically WORSE as the window grew (1733 -> 3781 ->
4060), which means window size was never the actual problem. All three reused the same r0 per
pulse, chosen once from narrow LOCAL windows on the raw data -- and this burst is packed tightly
enough (5 pulses in ~3s) that those local windows likely still contain real neighbor contamination
in several of the r0 choices. A wrong, fixed r0 forces the joint (A, t_peak, t_v) fit to contort
itself to compensate; a wider window just gives it more room to find different bad contortions.

Fix: alternate between (a) a joint fit at the current r0's, and (b) re-deriving each pulse's r0
via the full section-5 `select_r0` on the JOINT-FIT-CONSISTENT residual (other 4 pulses' current
best fit subtracted from the data), not the original isolated local window -- the same
"iterate once if the joint fit moves things" idea already validated on GRB231129C
(`GRB231129C/run_pipeline.py`'s pass1->pass2), generalized here to multiple alternating rounds
since this burst's initial r0's are apparently badly contaminated, not just slightly off.

Window: the goldilocks (-1, 11)s (physically the right choice -- contains the real decay tail
without reaching into pure noise), combined with the r0 fix this time rather than instead of it.
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
# Extended per user (2026-09-25): pulse 6 ~5.6s (low amplitude, broad -- the tail-anchoring pulse
# pulse 1 was wrongly standing in for) and a sharp pulse, originally hinted ~3.0s -- pulled back to
# ~2.2s after the 2.7/3.0 pair (0.3s apart) turned out too degenerate to resolve stably (identity
# swap, failed reproducibility). 7 total, matching the original "potentially 7" hint. Sorted so
# local_window()'s neighbor logic (nearest still-unprocessed candidate) works regardless of order.
HINT_CENTERS = sorted([0.5, 1.3, 2.0, 2.2, 2.7, 3.5, 5.6])
MIN_LOCAL_WIDTH = 0.30
MAX_LOCAL_WIDTH = 1.0
LOCAL_HINT_WINDOW = (-0.5, 7.0)  # widened from 6.0 to give pulse 6 (~5.6s) room for its own local fit
N_R0_ITERATIONS = 4
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
            a = float(np.clip(a0 * (1 + rng.uniform(-0.1, 0.1)), lb[3 * i], ub[3 * i]))
            tp = float(np.clip(tp0 + rng.uniform(-3, 3) * dt, lb[3 * i + 1], ub[3 * i + 1]))
            tv = float(np.clip(tv0 * (1 + rng.uniform(-0.1, 0.1)), lb[3 * i + 2], ub[3 * i + 2]))
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


def fit_joint(t_w, y_w, sigma_w, r0_list, p0_flat, dt):
    n = len(r0_list)
    model = make_joint_pulse3(r0_list)
    lb, ub = joint_bounds(t_w, dt, n)
    popt, pcov = curve_fit(model, t_w, y_w, p0=p0_flat, bounds=(lb, ub), sigma=sigma_w, absolute_sigma=True, maxfev=60000)
    return model, popt, pcov


def main():
    print(f"=== GRB131014215: JOINT 5-pulse fit with ITERATED r0, window {WINDOW} ===")
    t_full, y_full, sigma_full = load_light_curve()
    dt = median_dt(t_full)
    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
    n_pulses = len(HINT_CENTERS)
    print(f"window: {WINDOW} ({t_w.size} bins)")

    print("\n--- Pass 0: local single-pulse fit per hint center (initial r0, from raw data) ---")
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
    print(f"initial r0's: {[f'{r:.4g}' for r in r0_list]}")

    p0_flat = flatten_seed(p0_list)
    model, popt, pcov = fit_joint(t_w, y_w, sigma_w, r0_list, p0_flat, dt)
    params = unflatten_params(popt, n_pulses)
    chi2 = chi_square(y_w, model(t_w, *popt), sigma=sigma_w)
    print(f"pass 0 joint fit: chi2={chi2:.2f}")

    history = [{"pass": 0, "r0": list(r0_list), "params": list(params), "chi2": chi2}]

    for pass_num in range(1, N_R0_ITERATIONS + 1):
        print(f"\n--- Pass {pass_num}: re-derive each pulse's r0 from the JOINT-FIT-CONSISTENT residual ---")
        new_r0_list = []
        for i in range(n_pulses):
            residual_i = residual_excluding(model, t_w, y_w, popt, exclude_index=i)
            selection_i = select_r0(t_w, residual_i, sigma=sigma_w, dt=dt)
            new_r0 = selection_i.r0_chosen
            print(f"  pulse {i+1}: r0 {r0_list[i]:.4g} -> {new_r0:.4g} ({selection_i.classification.shape}, "
                  f"delta_chi2={selection_i.delta_chi2['delta_chi2']})")
            new_r0_list.append(new_r0)

        r0_change = max(abs(np.log10(new_r0_list[i]) - np.log10(r0_list[i])) for i in range(n_pulses))
        r0_list = new_r0_list
        p0_flat = tuple(popt)  # seed the next joint fit from where we already are

        model, popt, pcov = fit_joint(t_w, y_w, sigma_w, r0_list, p0_flat, dt)
        params = unflatten_params(popt, n_pulses)
        chi2 = chi_square(y_w, model(t_w, *popt), sigma=sigma_w)
        print(f"pass {pass_num} joint fit: chi2={chi2:.2f}  (max log10(r0) change: {r0_change:.3f})")
        history.append({"pass": pass_num, "r0": list(r0_list), "params": list(params), "chi2": chi2})

        if r0_change < 0.05:  # r0's stopped moving meaningfully (< ~12% change)
            print(f"r0's converged (max change {r0_change:.3f} < 0.05 in log10) -- stopping early")
            break

    print(f"\n=== Final (pass {history[-1]['pass']}) per-pulse results ===")
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
    csv_path = HERE / "GRB131014215_joint_iterated_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path.name}")

    hist_df = pd.DataFrame([{"pass": h["pass"], "chi2": h["chi2"], **{f"r0_pulse{i+1}": h["r0"][i] for i in range(n_pulses)}} for h in history])
    hist_path = HERE / "GRB131014215_joint_iterated_history.csv"
    hist_df.to_csv(hist_path, index=False)
    print(f"wrote {hist_path.name}")
    print(hist_df.to_string(index=False))

    print("\n--- Audit (final pass) ---")
    lb, ub = joint_bounds(t_w, dt, n_pulses)
    pinned = check_edge_pinning(popt, lb, ub, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")
    rng = get_rng(seed=SEED)
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(t_w, y_w, sigma_w, r0_list, tuple(popt), dt, rng)
    print(f"multi-start reproducibility: {'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")

    print("\n--- Plot (display boxed to T90 +/- 2s) ---")
    plot_lo, plot_hi = T90[0] - PLOT_PAD_S, T90[1] + PLOT_PAD_S
    box_mask = (t_full >= plot_lo) & (t_full <= plot_hi)
    paths = plot_joint_overlay(
        t_full[box_mask], y_full[box_mask], sigma_full[box_mask], model, popt, out_dir=str(HERE),
        label="GRB131014215_joint_iterated", pulse_indices=list(range(1, n_pulses + 1)),
    )
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print(f"\n=== Verdict: reproducibility {'PASS -- r0 iteration fixed it' if repro_ok else 'STILL FAILING even with iterated r0'} ===")


if __name__ == "__main__":
    main()
