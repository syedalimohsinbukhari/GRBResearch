"""Follow-up to fit_comparison.py, GRB140206B only: does pulse 7's t_s pinning (see comparison.md)
resolve if the fit window is widened all the way to the light curve's own x.min()/x.max(), rather than
just the "wide_-20_300" window already used as this burst's robustness check in
../window_sensitivity_GRB140206275/window_sensitivity.py?

Three windows, both fitting methods (normalized NorrisFitter, unnormalized explicit-x_scale
least_squares), same P0/bounds/max_nfev as fit_comparison.py -- reuses that script's functions
directly (import, not copy) so the two experiments can't silently drift apart on the fitting logic
itself, only on the window.

Does not modify fit_comparison.py or anything outside this folder.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fit_comparison import (
    BURSTS, ENERGY_LOW, ENERGY_HIGH, LC_DIR, HERE,
    summed_nai_curve, fit_normalized, fit_unnormalized_xscale, total_model,
)
from norris_fit import t_peak, tv_value
from grb_research import update_style, LINE_WIDTH
from grb_research.grb_utils import save_fig

GRB_DIR_NAME = "GRB140206275"
CFG = BURSTS[GRB_DIR_NAME]

# Full light-curve range, read directly from one detector's own .dat file (not hardcoded) -- same
# "widest_full_range" convention as ../window_sensitivity_GRB140206275/window_sensitivity.py.
_t_full, _, _ = summed_nai_curve(LC_DIR / GRB_DIR_NAME, (-np.inf, np.inf))
T_MIN, T_MAX = float(_t_full.min()), float(_t_full.max())

WINDOWS = {
    "narrow_-1_160": (-1, 160),
    "wide_-20_300": (-20, 300),
    "full_x.min()_x.max()": (T_MIN, T_MAX),
}


def run_window(window_label: str, window: tuple):
    t, y_raw, dat_nai = summed_nai_curve(LC_DIR / GRB_DIR_NAME, window)
    y_max = float(np.max(y_raw))
    n_pulses = len(CFG["p0"])

    params_norm = fit_normalized(t, y_raw, y_max, CFG["p0"], CFG["max_nfev"])
    model_norm_raw = total_model(t, params_norm, n_pulses) * y_max
    sse_norm = float(np.sum((model_norm_raw - y_raw) ** 2))

    res_raw = fit_unnormalized_xscale(t, y_raw, y_max, CFG["p0"], n_pulses, CFG["max_nfev"])
    params_raw = res_raw.x
    model_raw = total_model(t, params_raw, n_pulses)
    sse_raw = float(np.sum((model_raw - y_raw) ** 2))

    rows = []
    for i in range(n_pulses):
        A_n, ts_n, tau1_n, tau2_n = params_norm[i * 4:(i + 1) * 4]
        A_r, ts_r, tau1_r, tau2_r = params_raw[i * 4:(i + 1) * 4]
        rows.append({
            "window": window_label, "window_min_s": window[0], "window_max_s": window[1],
            "pulse_index": i + 1, "y_max_cts_per_s": y_max,
            "normalized_A_norm": A_n, "normalized_A_cts_per_s": A_n * y_max,
            "normalized_t_s": ts_n, "normalized_tau1": tau1_n, "normalized_tau2": tau2_n,
            "normalized_t_peak_s": t_peak(ts_n, tau1_n, tau2_n), "normalized_t_v_s": tv_value(tau1_n, tau2_n),
            "unnorm_xscale_A_cts_per_s": A_r,
            "unnorm_xscale_t_s": ts_r, "unnorm_xscale_tau1": tau1_r, "unnorm_xscale_tau2": tau2_r,
            "unnorm_xscale_t_peak_s": t_peak(ts_r, tau1_r, tau2_r), "unnorm_xscale_t_v_s": tv_value(tau1_r, tau2_r),
            "sse_normalized_rescaled": sse_norm, "sse_unnorm_xscale": sse_raw,
            "sse_ratio_unnorm_over_norm": sse_raw / sse_norm,
        })
    return pd.DataFrame(rows), (t, y_raw, model_norm_raw, model_raw)


def main():
    all_rows = []
    curves = {}
    for label, window in WINDOWS.items():
        print(f"=== {label} = {window} ===")
        df, curve = run_window(label, window)
        all_rows.append(df)
        curves[label] = curve
        print(f"  SSE ratio (unnorm_xscale / normalized): {df['sse_ratio_unnorm_over_norm'].iloc[0]:.4f}")
        p6 = df[df.pulse_index == 6].iloc[0]
        print(f"  pulse 6: t_s_norm={p6.normalized_t_s:.3f}  t_s_unnorm={p6.unnorm_xscale_t_s:.3f}"
              f"  t_peak_norm={p6.normalized_t_peak_s:.3f}  t_peak_unnorm={p6.unnorm_xscale_t_peak_s:.3f}"
              f"  t_v_norm={p6.normalized_t_v_s:.3f}  t_v_unnorm={p6.unnorm_xscale_t_v_s:.3f}")
        p7 = df[df.pulse_index == 7].iloc[0]
        print(f"  pulse 7: t_s_norm={p7.normalized_t_s:.3f}  t_s_unnorm={p7.unnorm_xscale_t_s:.3f}"
              f"  t_peak_norm={p7.normalized_t_peak_s:.3f}  t_peak_unnorm={p7.unnorm_xscale_t_peak_s:.3f}"
              f"  t_v_norm={p7.normalized_t_v_s:.3f}  t_v_unnorm={p7.unnorm_xscale_t_v_s:.3f}")

    combined = pd.concat(all_rows, ignore_index=True)
    csv_path = HERE / "grb140206b_window_scan_results.csv"
    combined.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path} ({len(combined)} rows)")

    # --- Plot: full light curve (widest range) with all three windows' normalized-fit total models
    # overlaid, restricted to each window's own domain -- shows how much of the burst each window
    # actually constrains the fit against.
    update_style()
    fig, ax = plt.subplots(figsize=(12, 6))
    t_full, y_full, _ = summed_nai_curve(LC_DIR / GRB_DIR_NAME, (T_MIN - 1, T_MAX + 1))
    ax.plot(t_full, y_full, color="0.75", lw=LINE_WIDTH * 0.5,
            label=f"10-400 keV NaI ({'+'.join(sorted(f.stem for f in (LC_DIR / GRB_DIR_NAME).glob('*.dat') if 'n' in f.stem))}, summed)\nBackground subtracted")
    colors = {"narrow_-1_160": "tab:blue", "wide_-20_300": "tab:orange", "full_x.min()_x.max()": "tab:green"}
    for label in WINDOWS:
        t_w, y_w, model_norm_w, model_raw_w = curves[label]
        ax.plot(t_w, model_norm_w, color=colors[label], lw=LINE_WIDTH, label=f"{label} (normalized fit)")
    ax.set_xlabel("Time since trigger [s]")
    ax.set_ylabel("Count rate [counts/s]")
    ax.set_title("GRB140206B: normalized-fit total model across three fit windows")
    ax.legend(fontsize="small")
    fig_path = HERE / "grb140206b_window_scan"
    save_fig(fig, fig_path)
    print(f"wrote {fig_path}.png/.pdf")

    return combined


if __name__ == "__main__":
    main()
