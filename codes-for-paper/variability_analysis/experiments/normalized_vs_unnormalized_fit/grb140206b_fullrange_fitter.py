"""Separate fitter for GRB140206B at the FULL light-curve range (x.min()/x.max()), following up on
grb140206b_window_scan.py: that scan found the full-range window resolves pulse 7's t_s pinning
(the narrow/wide windows both pinned it against their own edge) but destabilizes pulse 6 instead --
the unnormalized (explicit-x_scale) fit's pulse 6 collapses to tau2=0.00045, the same
near-delta-function pathology that got GRB080916C's pulse 5 dropped from production.

This script re-seeds ONLY pulse 6, changing (A, t_s, tau1, tau2) from COMPLEX_P0's original
(0.23, 23, 83, 1.1). tau1/tau2 reset to a neutral (1, 1) seed, the same "uninformative seed"
reproducibility check already used elsewhere in this project (GRB131014A pulses 1/2, GRB080916C's
norris2-1/norris2-2, both in variability_analysis.md) to test whether a pulse's converged solution is
real (reproduces from a different starting point) or an artifact of where the seed happened to start.
t_s is left at the original 23. A=0.23 (unchanged) was tried first and rejected: it converges, but the
pulse drifts to t_s~5-17s (out of its intended ~23-28s region) and two other pulses pick up new
degeneracies as a side effect. **A=0.2 is the value that actually fixes it** (see comparison.md's
"Recommended configuration" section) -- current value below. Every other pulse's seed is unchanged
from COMPLEX_P0.

Both fitting methods (normalized NorrisFitter, unnormalized explicit-x_scale least_squares) are run,
reusing fit_comparison.py's functions directly (import, not copy).

Does not modify fit_comparison.py, grb140206b_window_scan.py, or anything outside this folder.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fit_comparison import (
    BURSTS, LC_DIR, HERE,
    summed_nai_curve, fit_normalized, fit_unnormalized_xscale, total_model,
)
from norris_fit import t_peak, tv_value
from grb_research import update_style, LINE_WIDTH
from grb_research.grb_utils import save_fig

GRB_DIR_NAME = "GRB140206275"
PAPER_NAME = "GRB140206B"
MAX_NFEV = BURSTS[GRB_DIR_NAME]["max_nfev"]

# Full light-curve range, read directly from the data (not hardcoded) -- same convention as
# grb140206b_window_scan.py / window_sensitivity_GRB140206275/window_sensitivity.py's
# "widest_full_range".
_t_full, _, _ = summed_nai_curve(LC_DIR / GRB_DIR_NAME, (-np.inf, np.inf))
T_MIN, T_MAX = float(_t_full.min()), float(_t_full.max())
WINDOW = (T_MIN, T_MAX)

# COMPLEX_P0 from fitter_GRB140206275.py, pulse 6 re-seeded: A=0.2 (was 0.23 -- see module docstring),
# t_s unchanged (23), (tau1, tau2) reset from the original (83, 1.1) to a neutral (1, 1) -- the
# uninformative-seed reproducibility check.
P0 = [
    (0.13, -0.3, 0.12, 1.68),
    (0.35, 4.0, 6, 14),
    (0.5, 11, 5, 1.4),
    (0.5, 28, 2, 1),
    (0.23, 24, 0.3, 1.434),
    (0.23, 23, 1, 1),  # pulse 6, re-seeded -- was (0.23, 23, 83, 1.1)
    (0.05, -0.9, 3400, 4.4),
]
N_PULSES = len(P0)


def main():
    t, y_raw, dat_nai = summed_nai_curve(LC_DIR / GRB_DIR_NAME, WINDOW)
    y_max = float(np.max(y_raw))
    print(f"window = {WINDOW}  (t.min()={T_MIN:.3f}, t.max()={T_MAX:.3f})  y_max_cts_per_s={y_max:.4f}")

    params_norm = fit_normalized(t, y_raw, y_max, P0, MAX_NFEV)
    model_norm_raw = total_model(t, params_norm, N_PULSES) * y_max
    sse_norm = float(np.sum((model_norm_raw - y_raw) ** 2))

    res_raw = fit_unnormalized_xscale(t, y_raw, y_max, P0, N_PULSES, MAX_NFEV)
    params_raw = res_raw.x
    model_raw = total_model(t, params_raw, N_PULSES)
    sse_raw = float(np.sum((model_raw - y_raw) ** 2))

    print(f"\nconverged: normalized=True  unnormalized_xscale={res_raw.success}")
    print(f"SSE ratio (unnorm/norm): {sse_raw / sse_norm:.4f}")

    rows = []
    print(f"\n{'pulse':>5}  {'A_norm':>8}  {'A_cts/s (n)':>12}  {'t_s (n)':>9}  {'tau1 (n)':>10}  {'tau2 (n)':>10}  "
          f"{'A_cts/s (u)':>12}  {'t_s (u)':>9}  {'tau1 (u)':>10}  {'tau2 (u)':>10}  {'t_peak (n)':>10}  {'t_peak (u)':>10}  {'t_v (n)':>8}  {'t_v (u)':>8}")
    for i in range(N_PULSES):
        A_n, ts_n, tau1_n, tau2_n = params_norm[i * 4:(i + 1) * 4]
        A_r, ts_r, tau1_r, tau2_r = params_raw[i * 4:(i + 1) * 4]
        tp_n, tp_r = t_peak(ts_n, tau1_n, tau2_n), t_peak(ts_r, tau1_r, tau2_r)
        tv_n, tv_r = tv_value(tau1_n, tau2_n), tv_value(tau1_r, tau2_r)
        print(f"{i + 1:5d}  {A_n:8.4f}  {A_n * y_max:12.4f}  {ts_n:9.4f}  {tau1_n:10.4f}  {tau2_n:10.6f}  "
              f"{A_r:12.4f}  {ts_r:9.4f}  {tau1_r:10.4f}  {tau2_r:10.6f}  {tp_n:10.4f}  {tp_r:10.4f}  {tv_n:8.4f}  {tv_r:8.4f}")
        rows.append({
            "grb_name": PAPER_NAME, "window_min_s": WINDOW[0], "window_max_s": WINDOW[1],
            "pulse_index": i + 1, "reseeded": (i + 1 == 6), "y_max_cts_per_s": y_max,
            "normalized_A_norm": A_n, "normalized_A_cts_per_s": A_n * y_max,
            "normalized_t_s": ts_n, "normalized_tau1": tau1_n, "normalized_tau2": tau2_n,
            "normalized_t_peak_s": tp_n, "normalized_t_v_s": tv_n,
            "unnorm_xscale_A_cts_per_s": A_r,
            "unnorm_xscale_t_s": ts_r, "unnorm_xscale_tau1": tau1_r, "unnorm_xscale_tau2": tau2_r,
            "unnorm_xscale_t_peak_s": tp_r, "unnorm_xscale_t_v_s": tv_r,
            "sse_normalized_rescaled": sse_norm, "sse_unnorm_xscale": sse_raw,
            "sse_ratio_unnorm_over_norm": sse_raw / sse_norm,
        })

    df = pd.DataFrame(rows)
    csv_path = HERE / "grb140206b_fullrange_pulse6_reseed_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path}")

    update_style()
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(t, y_raw, color="0.6", lw=LINE_WIDTH * 0.6,
            label=f"10-400 keV NaI ({'+'.join(dat_nai)}, summed)\nBackground subtracted, full range")
    ax.plot(t, model_norm_raw, color="tab:blue", lw=LINE_WIDTH, label="Normalized fit total")
    ax.plot(t, model_raw, color="tab:red", ls="--", lw=LINE_WIDTH, label="Unnormalized fit total (explicit x_scale)")
    ax.set_xlabel("Time since trigger [s]")
    ax.set_ylabel("Count rate [counts/s]")
    ax.set_title(f"{PAPER_NAME}: full-range fit, pulse 6 re-seeded to (A, t_s, 1, 1)")
    ax.legend(fontsize="small")
    fig_path = HERE / "grb140206b_fullrange_pulse6_reseed"
    save_fig(fig, fig_path)
    print(f"wrote {fig_path}.png/.pdf")

    return df


if __name__ == "__main__":
    main()
