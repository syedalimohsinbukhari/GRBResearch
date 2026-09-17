r"""
Experiment: same window-widening pulse-count diagnostic as
experiments/window_sensitivity_GRB{131014215,140206275,231129779}/, applied here to GRB080916C's
6-vs-7-pulse question flagged in fitter.py (the user's own working file): a commented-out 7-pulse p0
adds one extra pulse (A=0.1, t_s=20, tau1=9, tau2=1) between the live 6-pulse model's 4th pulse
(t_s=1.3, a broad tau1=43/tau2=14 component) and 5th pulse (t_s=52), sitting inside TR3's window
(15.040-55.296s) -- the same region already flagged in PHASE5_TV_PLAN.md as "near-featureless" around
TR3's own 27.4 GeV defining photon at t=40.5s.

Mechanism under test (same as the other three bursts):

NorrisFitter.fit_boundaries() (codes-for-paper/variability_analysis/norris_fit.py:49-50) sets t_s's
lower bound to x_values.min() -- i.e., the fit window's own left edge is t_s's box constraint. Widening
the window directly widens how far the optimizer can push t_s along the t_s<->tau1 degeneracy. Per the
methodological conclusion already established (variability_analysis.md): if a pulse's mc_kept_fraction
collapses under widening in one pulse-count model but not another, that's evidence the collapsing
model is under-specified.

Stage 1: fit both models (their own raw p0 from fitter.py, not yet converged) over the narrow window
[-1, 70]s (fitter.py's own default). Stage 2: refit both at wide/widest windows using each model's own
narrow-window converged parameters as p0 (isolating the window-bound effect from a seed-choice effect,
same as every other burst's check).
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
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB080916009"
ENERGY_LOW, ENERGY_HIGH = 10, 400

dat_NaI = [f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f]
t_full, r1, b1 = lightcurve_data(f"{LC_DIR}/{dat_NaI[0]}.dat", ENERGY_LOW, ENERGY_HIGH)
y_full = r1 - b1

T_MIN, T_MAX = float(np.min(t_full)), float(np.max(t_full))

WINDOWS = {
    "narrow_-1_70": (-1, 70),
    "wide_-20_150": (-20, 150),
    "widest_full_range": (T_MIN, T_MAX),
}

# Raw p0 from fitter.py -- SIX pulses (the file's currently-live, uncommented fit).
P0_6 = [
    (0.4, -0.7, 2.34, 0.471),
    (0.7, 0.3, 0.88, 6),
    (0.2, 5.3, 0.6, 6),
    (0.3, 1.3, 43, 14),
    (0.3, 52, 18, 0.7),
    (0.3, 61, 0.3, 3),
]

# Raw p0 from fitter.py -- SEVEN pulses (the file's commented-out alternative): identical to P0_6 except
# for one extra pulse (A=0.1, t_s=20, tau1=9, tau2=1) inserted, sitting inside TR3's window.
P0_7 = [
    (0.4, -0.7, 2.34, 0.471),
    (0.7, 0.3, 0.88, 6),
    (0.2, 5.3, 0.6, 6),
    (0.3, 1.3, 43, 14),
    (0.1, 20, 9, 1),
    (0.3, 52, 18, 0.7),
    (0.3, 61, 0.3, 3),
]

MODELS_RAW = {"SIX": P0_6, "SEVEN": P0_7}

def fit_window(p0, start, stop):
    mask = np.logical_and(t_full > start, t_full < stop)
    t_w, y_w = t_full[mask], y_full[mask]
    y_max = np.max(y_w)
    y_norm = y_w / y_max
    nf = NorrisFitter(t_w, y_norm)
    nf.fit(p0=p0)
    return nf, y_max

def summarize(nf, n_pulses, label):
    rows = []
    print(f"\n=== {label} ({n_pulses} pulses) ===")
    for i in range(n_pulses):
        A, ts, tau1, tau2 = nf.params[i * 4 : (i + 1) * 4]
        tp = t_peak(ts, tau1, tau2)
        mc = tv_mc_summary(nf, pulse_index=i + 1)
        print(
            f"  pulse {i + 1}: A={A:.4f} ts={ts:9.4f} tau1={tau1:10.4f} tau2={tau2:8.4f}  "
            f"t_peak={tp:9.4f}  t_v={mc['t_v_s']:.4f} +{mc['t_v_err_upper_s']:.4f} "
            f"-{mc['t_v_err_lower_s']:.4f}  kept={mc['kept_fraction']:.3f}"
        )
        rows.append(
            {"pulse_index": i + 1, "A_norm": A, "t_s": ts, "tau1": tau1, "tau2": tau2,
             "t_peak_s": tp, "t_v_s": mc["t_v_s"], "t_v_err_lower_s": mc["t_v_err_lower_s"],
             "t_v_err_upper_s": mc["t_v_err_upper_s"], "mc_kept_fraction": mc["kept_fraction"],
             "n_samples": mc["n_samples"], "seed": mc["seed"]}
        )
    return rows

# --- Stage 1: narrow-window fit from raw p0, both models.
start, stop = WINDOWS["narrow_-1_70"]
converged_p0 = {}
all_rows = []
for model_name, p0 in MODELS_RAW.items():
    nf, y_max = fit_window(p0, start, stop)
    rows = summarize(nf, len(p0), f"{model_name}-PULSE, narrow_-1_70 (raw p0)")
    for r in rows:
        r.update({"model": model_name, "window": "narrow_-1_70", "window_start_s": start, "window_stop_s": stop})
    all_rows.extend(rows)
    converged_p0[model_name] = [tuple(nf.params[i * 4 : (i + 1) * 4]) for i in range(len(p0))]

    nf.plot_fit(
        show_individuals=True,
        x_label="Time since trigger [s]",
        y_label="Normalized count rate",
        data_label="10-400 keV NaI\nBackground Subtracted",
        title=f"GRB080916C {model_name}-pulse window-sensitivity check: narrow_-1_70 [{start}, {stop}] s",
    )
    fig_path = Path(__file__).parent / f"window_sensitivity_{model_name.lower()}_narrow_-1_70"
    save_fig(plt.gcf(), fig_path)

# --- Stage 2: wide/widest, using each model's own narrow-window converged parameters as p0.
for model_name, p0 in converged_p0.items():
    n_pulses = len(p0)
    for window_name in ["wide_-20_150", "widest_full_range"]:
        start, stop = WINDOWS[window_name]
        nf, y_max = fit_window(p0, start, stop)
        rows = summarize(nf, n_pulses, f"{model_name}-PULSE, {window_name}")
        for r in rows:
            r.update({"model": model_name, "window": window_name, "window_start_s": start, "window_stop_s": stop})
        all_rows.extend(rows)

        nf.plot_fit(
            show_individuals=True,
            x_label="Time since trigger [s]",
            y_label="Normalized count rate",
            data_label="10-400 keV NaI\nBackground Subtracted",
            title=f"GRB080916C {model_name}-pulse window-sensitivity check: {window_name} [{start}, {stop}] s",
        )
        fig_path = Path(__file__).parent / f"window_sensitivity_{model_name.lower()}_{window_name}"
        save_fig(plt.gcf(), fig_path)

df = pd.DataFrame(all_rows)
csv_path = Path(__file__).parent / "window_sensitivity_results.csv"
df.to_csv(csv_path, index=False)
print(f"\nwrote {csv_path} ({len(df)} rows)")

# --- Side-by-side comparison per pulse, per model, across all three windows
WINDOW_ORDER = ["narrow_-1_70", "wide_-20_150", "widest_full_range"]
for model_name, p0 in converged_p0.items():
    n_pulses = len(p0)
    print(f"\n=== {model_name}-PULSE: narrow vs wide vs widest, per pulse ===")
    sub = df[df["model"] == model_name]
    by_window = {w: sub[sub["window"] == w].set_index("pulse_index") for w in WINDOW_ORDER}
    narrow = by_window["narrow_-1_70"]
    for i in range(1, n_pulses + 1):
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
