r"""
Experiment: same window-widening pulse-count diagnostic as
experiments/window_sensitivity_GRB231129779/window_sensitivity.py, applied here to GRB131014A's
5-pulse fit (fitter_CLAUDE_GRB131014215.py) -- flagged as not-yet-tested in variability_analysis.md's
"Known limitations" section for pulses 1/2.

Mechanism under test:

NorrisFitter.fit_boundaries() (codes-for-paper/variability_analysis/norris_fit.py:49-50) sets t_s's
lower bound to x_values.min() -- i.e., the fit window's own left edge is t_s's box constraint. Widening
the window directly widens how far the optimizer can push t_s along the t_s<->tau1 degeneracy.

Controlled comparison: SAME p0 (the narrow-window run's own converged parameters, read from
norris_fit_results_GRB131014215.csv) fit to the SAME light curve, with ONLY the window bound changed
(narrow [-1,10] -- the fitter_CLAUDE_GRB131014215.py default -- vs wide [-10,20], matching the
GRB231129C experiment's window choice for direct comparability).
"""

import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from light_curves import lightcurve_data
from norris_fit import NorrisFitter, t_peak, tv_mc_summary

from grb_research import update_style
from grb_research.grb_utils import save_fig

update_style()

PROJECT_ROOT = Path(__file__).resolve().parents[4]
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB131014215"
ENERGY_LOW, ENERGY_HIGH = 10, 400

dat_NaI = [f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f]
t_full, r1, b1 = lightcurve_data(f"{LC_DIR}/{dat_NaI[0]}.dat", ENERGY_LOW, ENERGY_HIGH)
y_full = r1 - b1

# Third, most extreme broadening (added 2026-09-16, per GRB140206B's own widest_full_range check):
# push both edges all the way to the light curve's own np.min(t)/np.max(t) -- the widest box
# NorrisFitter.fit_boundaries() can ever present for this data.
T_MIN, T_MAX = float(np.min(t_full)), float(np.max(t_full))

WINDOWS = {
    "narrow_-1_10": (-1, 10),
    "wide_-10_20": (-10, 20),
    "widest_full_range": (T_MIN, T_MAX),
}

# Converged parameters from the narrow-window run (norris_fit_results_GRB131014215.csv), used as p0
# for *all three* windows -- same starting point, only the window bound differs.
P0 = [
    (0.10764203288945787, -0.6801299026283388, 4.736463927778451, 0.3111610774862506),
    (0.2745269047316325, 0.046601440392526616, 25.344683548125637, 0.05606570780848216),
    (0.9341178707469664, 1.1832335152138767, 1.012741437487432, 0.3358238707570263),
    (0.45875295632710456, 2.4705042509819695, 0.04738191980389777, 1.099443921612677),
    (0.2904758094203136, 2.658479250500718, 6.545740339696772, 0.13247748600998696),
]
N_PULSES = len(P0)

rows = []
for window_name, (start, stop) in WINDOWS.items():
    mask = np.logical_and(t_full > start, t_full < stop)
    t_w, y_w = t_full[mask], y_full[mask]
    y_max = np.max(y_w)
    y_norm = y_w / y_max

    nf = NorrisFitter(t_w, y_norm)
    nf.fit(p0=P0)

    nf.plot_fit(
        show_individuals=True,
        x_label="Time since trigger [s]",
        y_label="Normalized count rate",
        data_label="10-400 keV NaI\nBackground Subtracted",
        title=f"GRB131014A window-sensitivity check: {window_name} [{start}, {stop}] s",
    )
    fig_path = Path(__file__).parent / f"window_sensitivity_{window_name}"
    save_fig(plt.gcf(), fig_path)

    print(f"\n=== window {window_name} ({start}, {stop}), y_max={y_max:.2f} cts/s ===")
    for i in range(N_PULSES):
        A, ts, tau1, tau2 = nf.params[i * 4 : (i + 1) * 4]
        tp = t_peak(ts, tau1, tau2)
        mc = tv_mc_summary(nf, pulse_index=i + 1)
        print(
            f"  pulse {i + 1}: A={A:.4f} ts={ts:9.4f} tau1={tau1:10.4f} tau2={tau2:8.4f}  "
            f"t_peak={tp:8.4f}  t_v={mc['t_v_s']:.4f} +{mc['t_v_err_upper_s']:.4f} "
            f"-{mc['t_v_err_lower_s']:.4f}  kept={mc['kept_fraction']:.3f}"
        )
        rows.append(
            {
                "window": window_name,
                "window_start_s": start,
                "window_stop_s": stop,
                "pulse_index": i + 1,
                "A_norm": A,
                "t_s": ts,
                "tau1": tau1,
                "tau2": tau2,
                "t_peak_s": tp,
                "t_v_s": mc["t_v_s"],
                "t_v_err_lower_s": mc["t_v_err_lower_s"],
                "t_v_err_upper_s": mc["t_v_err_upper_s"],
                "mc_kept_fraction": mc["kept_fraction"],
                "n_samples": mc["n_samples"],
                "seed": mc["seed"],
            }
        )

df = pd.DataFrame(rows)
csv_path = Path(__file__).parent / "window_sensitivity_results.csv"
df.to_csv(csv_path, index=False)
print(f"\nwrote {csv_path} ({len(df)} rows)")

# --- Side-by-side comparison per pulse, across all three windows
WINDOW_ORDER = ["narrow_-1_10", "wide_-10_20", "widest_full_range"]
by_window = {w: df[df["window"] == w].set_index("pulse_index") for w in WINDOW_ORDER}
narrow = by_window["narrow_-1_10"]
print("\n=== narrow vs wide vs widest, per pulse ===")
for i in range(1, N_PULSES + 1):
    n = narrow.loc[i]
    parts = [f"pulse {i}: t_s {n.t_s:9.3f}"]
    for w in WINDOW_ORDER[1:]:
        row = by_window[w].loc[i]
        dt_v_pct = 100 * abs(row.t_v_s - n.t_v_s) / n.t_v_s
        parts.append(
            f" -> {row.t_s:9.3f} (t_v={row.t_v_s:.4f}, kept={row.mc_kept_fraction:.3f}, "
            f"dt_v={dt_v_pct:.1f}% vs narrow)"
        )
    print("".join(parts))
