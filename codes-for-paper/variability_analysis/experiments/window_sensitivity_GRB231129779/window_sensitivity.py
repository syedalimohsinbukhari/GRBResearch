r"""
Experiment:

Does the Norris-fit *window bound* alone (holding p0, data, and pulse count fixed) move the fitted parameters and
:math:`t_v` the way the user observed comparing
`norris_fit_GRB231129779.png ([-1,10])` against norris_fit_GRB231129779__bkp.png ([-10,20])?

Mechanism under test:

NorrisFitter.fit_boundaries() (variability_timescale/norris_fit.py:49-53) sets :math:`t_s`'s lower bound to
x_values.min() -- i.e., the fit window's own left edge is :math:`t_s`'s box constraint.
Combined with the already-established near-total :math:`t_s \leftrightarrow \tau_1` degeneracy (this session's
correlation matrices, TR4, the (:math:`\tau`, :math:`\xi`) reparam test), widening the window directly
widens how far the optimizer can push t_s along that degenerate ridge.

Controlled comparison: SAME p0 (the narrow-window run's own converged parameters, read off
norris_fit_GRB231129779.png's legend) fit to the SAME light curve, with ONLY the window bound changed.
Isolates the window-bound effect from a seed-choice effect (already characterized separately in
variability_analysis.md's guess-sensitivity test).
"""

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
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB231129779"
ENERGY_LOW, ENERGY_HIGH = 10, 400

import os

dat_NaI = sorted(f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f)  # order no longer load-bearing -- see BUG-23 fix below; kept sorted for deterministic logging only
nai_data = [lightcurve_data(f"{LC_DIR}/{d}.dat", ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]

# BUG-23 fix (2026-09-23, user decision): sum all of the burst's NaI detectors' background-subtracted
# count rates raw, with no per-detector normalization -- matching fitter_GRB231129779.py's own fix and the
# existing "n3+n4 summed" ROOT cross-check precedent (variability_analysis.md). This structurally closes
# BUG-23: there's no longer a single detector to mis-pick via os.listdir()[0]. Detector time grids are
# confirmed identical before summing, not assumed. (GRB231129C's own os.listdir()[0] pick happened to be
# correct already, per BUGS.md, but this removes the coincidence, not just the bug.)
t_full = nai_data[0][0]
for _det, (_t, _, _) in zip(dat_NaI[1:], nai_data[1:]):
    assert np.array_equal(t_full, _t), f"{_det}'s time grid differs from {dat_NaI[0]}'s -- cannot sum"
r1 = np.sum([r for _, r, _ in nai_data], axis=0)
b1 = np.sum([b for _, _, b in nai_data], axis=0)
y_full = r1 - b1

# Third, most extreme broadening (added 2026-09-16, per GRB140206B's own widest_full_range check):
# push both edges all the way to the light curve's own np.min(t)/np.max(t) -- the widest box
# NorrisFitter.fit_boundaries() can ever present for this data.
T_MIN, T_MAX = float(np.min(t_full)), float(np.max(t_full))

WINDOWS = {
    # "narrow_-1_10": (-1, 10),
    # "wide_-10_20": (-10, 20),
    "widest_full_range": (T_MIN, T_MAX),
}

# Converged parameters from the narrow-window run (norris_fit_GRB231129779.png's own legend),
# used as p0 for *all three* windows -- same starting point, only the window bound differs.
P0 = [
    (0.570, -0.449, 1.854, 0.637),
    (0.493, 0.127, 3.912, 0.381),
    (0.708, -0.960, 75.728, 0.154),
    (0.472, 2.759, 1.367, 1.470),
    (0.125, 4.609, 0.429, 2.149)
]
N_PULSES = len(P0)

rows = []
fitters = {}
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
        title=f"GRB231129C window-sensitivity check: {window_name} [{start}, {stop}] s",
    )
    fig_path = Path(__file__).parent / f"window_sensitivity_{window_name}"
    save_fig(plt.gcf(), fig_path)
    fitters[window_name] = nf

    print(f"\n=== window {window_name} ({start}, {stop}), y_max={y_max:.2f} cts/s ===")
    for i in range(N_PULSES):
        A, ts, tau1, tau2 = nf.params[i * 4: (i + 1) * 4]
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
