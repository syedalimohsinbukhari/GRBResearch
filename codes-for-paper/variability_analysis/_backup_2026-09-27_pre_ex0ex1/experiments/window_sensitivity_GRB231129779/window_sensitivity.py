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

from grb_research import update_style, seed_from_name
from grb_research.grb_utils import save_fig

update_style()

# Deterministic, per-script seed -- project convention (SEEDING.md); added 2026-09-26 (BUG-25),
# every tv_mc_summary() call below used to omit seed= entirely, silently falling back to
# norris_fit.py's own hardcoded SEED=12345 default, which that module no longer has.
SEED = seed_from_name(__file__)

PROJECT_ROOT = Path(__file__).resolve().parents[4]
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB231129779"
ENERGY_LOW, ENERGY_HIGH = 10, 400

import os

dat_NaI = sorted(
    f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f
)  # order no longer load-bearing -- see BUG-23 fix below; kept sorted for deterministic logging only
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
    "narrow_-1_10": (-1, 10),
    "wide_-10_20": (-10, 20),
    "widest_full_range": (T_MIN, T_MAX)
}

# Re-seeded 2026-09-27 from this file's own pre-BUG-23 widest_full_range converged result
# (window_sensitivity_results.csv, 2026-09-16 run) instead of the narrow-window converged
# parameters -- the narrow seed made the optimizer travel from a tiny [-1,10]s box all the way to
# the full [-138,476]s one (NorrisFitter.fit_boundaries() ties t_s's lower bound to the window's
# own left edge), which took >20 min of CPU time without converging. Shape/timing parameters
# transfer across the BUG-23 detector-summing fix fine since the fit always runs on the
# peak-normalized curve (y_norm = y / y.max()), not raw counts -- only the summing changes S/N,
# not pulse shape. Only WINDOWS' one active entry (widest_full_range) uses this P0 currently.
P0_THIN = [
    (0.570, -0.449, 1.854, 0.637),
    (0.493, 0.127, 3.912, 0.381),
    (0.708, -0.960, 75.728, 0.154),
    (0.472, 2.759, 1.367, 1.470),
    (0.125, 4.609, 0.429, 2.149)
]

P0_WIDE = [
    (0.7, -1.14, 90, 0.14),
    (0.570, -0.42, 1.6, 0.7),
    (0.45, 0.16, 3.6, 0.38),
    (0.47, 2.5, 2.5, 0.8),
    (0.1, 4.6, 0.26, 2.5)
]

P0_WIDEST = [
    (0.5, -0.35, 1.2, 0.7),
    (0.7, -5.8, 1065, 0.05),
    (0.5, 0.218, 2.66, 0.4),
    (0.3, 3.5, 1.3, 0.8),
    (0.4, 2.7, 1.7, 1.0),
    (0.3, 4.6, 0.2, 2.57),
]

get_P0 = {
    "narrow_-1_10": P0_THIN,
    "wide_-10_20": P0_WIDE,
    "widest_full_range": P0_WIDEST
}

rows = []
fitters = {}
for window_name, (start, stop) in WINDOWS.items():
    mask = np.logical_and(t_full > start, t_full < stop)
    t_w, y_w = t_full[mask], y_full[mask]
    y_max = np.max(y_w)
    y_norm = y_w / y_max

    P0 = get_P0[window_name]
    N_PULSES = len(P0)

    # max_iterations=20000 (was the NorrisFitter default of 5000) -- same BUG-23 convention as every
    # sibling window_sensitivity_GRB*.py: the summed-detector curve needs more iterations than the
    # single-detector-tuned default.
    nf = NorrisFitter(t_w, y_norm, max_iterations=20000)
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

    print(f"=== window {window_name} ({start}, {stop}), y_max={y_max:.2f} cts/s ===")
    for i in range(N_PULSES):
        A, ts, tau1, tau2 = nf.params[i * 4 : (i + 1) * 4]
        tp = t_peak(ts, tau1, tau2)
        mc = tv_mc_summary(nf, pulse_index=i + 1, seed=SEED)
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

# --- Side-by-side comparison per pulse, across whichever windows actually ran above.
# Was hardcoded to all three (narrow/wide/widest) regardless of WINDOWS -- crashed here once
# narrow/wide were commented out above (2026-09-27), since by_window["narrow_-1_10"] was an empty
# DataFrame with nothing to .loc[i] into. Now derived from WINDOWS itself, and skipped outright
# with fewer than two windows (nothing to compare).
WINDOW_ORDER = [w for w in ["narrow_-1_10", "wide_-10_20", "widest_full_range"] if w in WINDOWS]
if len(WINDOW_ORDER) < 2:
    print(f"\nOnly one window ({WINDOW_ORDER[0]}) was run -- nothing to compare, skipping.")
else:
    by_window = {w: df[df["window"] == w].set_index("pulse_index") for w in WINDOW_ORDER}
    narrow = by_window[WINDOW_ORDER[0]]
    print(f"\n=== {' vs '.join(WINDOW_ORDER)}, per pulse ===")
    for i in range(1, N_PULSES + 1):
        n = narrow.loc[i]
        parts = [f"pulse {i}: t_s {n.t_s:9.3f}"]
        for w in WINDOW_ORDER[1:]:
            row = by_window[w].loc[i]
            dt_v_pct = 100 * abs(row.t_v_s - n.t_v_s) / n.t_v_s
            parts.append(
                f" -> {row.t_s:9.3f} (t_v={row.t_v_s:.4f}, kept={row.mc_kept_fraction:.3f}, "
                f"dt_v={dt_v_pct:.1f}% vs {WINDOW_ORDER[0]})"
            )
        print("".join(parts))
