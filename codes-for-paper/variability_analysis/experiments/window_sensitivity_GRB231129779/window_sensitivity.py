"""Experiment: does the Norris-fit *window bound* alone (holding p0, data, and pulse count
fixed) move the fitted parameters and t_v the way the user observed comparing
norris_fit_GRB231129779.png ([-1,10]) against norris_fit_GRB231129779__bkp.png ([-10,20])?

Mechanism under test: NorrisFitter.fit_boundaries() (variability_timescale/norris_fit.py:49-53)
sets t_s's lower bound to x_values.min() -- i.e. the fit window's own left edge is t_s's box
constraint. Combined with the already-established near-total t_s<->tau1 degeneracy (this
session's correlation matrices, TR4, the (tau,xi) reparam test), widening the window directly
widens how far the optimizer can push t_s along that degenerate ridge.

Controlled comparison: SAME p0 (the narrow-window run's own converged parameters, read off
norris_fit_GRB231129779.png's legend) fit to the SAME light curve, with ONLY the window bound
changed. Isolates the window-bound effect from a seed-choice effect (already characterized
separately in variability_analysis.md's guess-sensitivity test).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from variability_timescale.light_curves import lightcurve_data
from variability_timescale.norris_fit import NorrisFitter, t_peak, tv_mc_summary

PROJECT_ROOT = Path(__file__).resolve().parents[4]
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB231129779"
ENERGY_LOW, ENERGY_HIGH = 10, 400

WINDOWS = {
    "narrow_-1_10": (-1, 10),
    "wide_-10_20": (-10, 20),
}

# Converged parameters from the narrow-window run (norris_fit_GRB231129779.png's own legend),
# used as p0 for *both* windows -- same starting point, only the window bound differs.
P0 = [
    (0.570, -0.449, 1.854, 0.637),
    (0.493, 0.127, 3.912, 0.381),
    (0.708, -0.960, 75.728, 0.154),
    (0.472, 2.759, 1.367, 1.470),
    (0.125, 4.609, 0.429, 2.149)
]
N_PULSES = len(P0)

import os

dat_NaI = [f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f]
nai_data = [lightcurve_data(f"{LC_DIR}/{d}.dat", ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]
t_full, r1, b1 = nai_data[0]
y_full = r1 - b1

rows = []
fitters = {}
for window_name, (start, stop) in WINDOWS.items():
    mask = np.logical_and(t_full > start, t_full < stop)
    t_w, y_w = t_full[mask], y_full[mask]
    y_max = np.max(y_w)
    y_norm = y_w / y_max

    nf = NorrisFitter(t_w, y_norm)
    nf.fit(p0=P0)
    nf.plot_fit(show_individuals=True)
    plt.show()
    fitters[window_name] = nf

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

# --- Side-by-side comparison per pulse
print("\n=== narrow vs wide, per pulse ===")
narrow = df[df["window"] == "narrow_-1_10"].set_index("pulse_index")
wide = df[df["window"] == "wide_-10_20"].set_index("pulse_index")
for i in range(1, N_PULSES + 1):
    n, w = narrow.loc[i], wide.loc[i]
    print(
        f"pulse {i}: t_s {n.t_s:9.3f} -> {w.t_s:9.3f}   tau1 {n.tau1:9.3f} -> {w.tau1:9.3f}   "
        f"t_peak {n.t_peak_s:.3f} -> {w.t_peak_s:.3f}   "
        f"t_v {n.t_v_s:.4f}+{n.t_v_err_upper_s:.4f}-{n.t_v_err_lower_s:.4f} (kept={n.mc_kept_fraction:.3f})"
        f"  ->  {w.t_v_s:.4f}+{w.t_v_err_upper_s:.4f}-{w.t_v_err_lower_s:.4f} (kept={w.mc_kept_fraction:.3f})"
    )
