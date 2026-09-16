r"""
Experiment: same window-widening pulse-count diagnostic as
experiments/window_sensitivity_GRB231129779/ and experiments/window_sensitivity_GRB131014215/,
applied here to GRB140206B -- the one burst in this project's manual-Norris-fit track with two
live candidate decompositions (fitter_GRB140206275_simple.py's 5-pulse SIMPLE model and
fitter_GRB140206275.py's 7-pulse COMPLEX model), so both are checked.

Mechanism under test:

NorrisFitter.fit_boundaries() (codes-for-paper/variability_analysis/norris_fit.py:49-50) sets t_s's
lower bound to x_values.min() -- i.e., the fit window's own left edge is t_s's box constraint. Widening
the window directly widens how far the optimizer can push t_s along the t_s<->tau1 degeneracy.

Controlled comparison: SAME p0 per model (that model's own narrow-window converged parameters, read
from norris_fit_results_GRB140206275.csv / norris_fit_results_GRB140206275_simple.csv) fit to the SAME
light curve, with ONLY the window bound changed, across THREE window widths. Narrow window matches
fitter_GRB140206275.py's own default (START1=-1, END1=160, covering all seven episodes through TR6's
154.240s end); wide window extends both edges well beyond that; widest pushes both edges all the way to
the light curve's own np.min(t)/np.max(t) -- the most extreme box NorrisFitter.fit_boundaries() can ever
present, added specifically to track where SIMPLE pulse 5's t_s (already seen pinned to the wide window's
left edge) ends up when the box constraint is removed as far as the data allows.

The COMPLEX model's own docstring in fitter_GRB140206275.py already flags that splitting the t~23-28s
region into three pulses (COMPLEX pulses 4/5/6) costs reliability relative to SIMPLE's single broad
pulse there (kept ~1.0 for SIMPLE's pulse 4 vs 0.54-0.67 for COMPLEX's 4/5/6) -- a deliberate,
photon-anchored-resolution-vs-stability tradeoff, not a mistake. This diagnostic checks whether that
already-poor COMPLEX reliability collapses further under window-widening (real trouble) or holds
(a legitimate tradeoff, not evidence of a wrong pulse count).
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
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB140206275"
ENERGY_LOW, ENERGY_HIGH = 10, 400

dat_NaI = [f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f]
t_full, r1, b1 = lightcurve_data(f"{LC_DIR}/{dat_NaI[0]}.dat", ENERGY_LOW, ENERGY_HIGH)
y_full = r1 - b1

# Third, most extreme broadening: the fit window's left/right edges pushed all the way out to the
# light curve's own np.min(t)/np.max(t) -- the absolute widest box NorrisFitter.fit_boundaries() can
# ever see for this data, so t_s's lower bound (x_values.min()) is as unconstrained as it can get.
# Added specifically to track pulse 5's t_s (SIMPLE: the broad ~3440s-tau1 pedestal, whose t_s was
# already seen pinned to the wide window's own left edge, -0.96 -> -19.97, in the two-window run).
T_MIN, T_MAX = float(np.min(t_full)), float(np.max(t_full))

WINDOWS = {
    "narrow_-1_160": (-1, 160),
    "wide_-20_300": (-20, 300),
    "widest_full_range": (T_MIN, T_MAX),
}

# Converged parameters from each model's narrow-window run (norris_fit_results_GRB140206275*.csv),
# used as p0 for *both* windows -- same starting point, only the window bound differs.
MODELS = {
    "SIMPLE": [
        (0.13821359231007244, -0.31068299950172557, 0.0786997853378577, 1.7179861641051446),
        (0.35672897543605114, 4.137546935144213, 5.976466850654499, 13.932789445290558),
        (0.5506659874868615, 11.076896556433432, 5.533197632966581, 1.4756390342157466),
        (0.2323517995016442, 26.73173459922519, 2.535166747797389, 3.786359998026647),
        (0.050309612777835425, -0.9518746258166284, 3439.7640412860924, 4.462751768530366),
    ],
    "COMPLEX": [
        (0.13715836098096673, -0.3082778877479421, 0.07280291527795804, 1.7779333409054536),
        (0.3485522723217019, 4.358824486438564, 5.423853787049376, 13.077114856867508),
        (0.5572846305208764, 11.202577021209974, 4.484415169647049, 1.6621123907230426),
        (0.2771673583293726, 23.06415137350121, 46.192192621651564, 1.0927353385142597),
        (0.0635676507074934, 23.081058799329575, 1.0587226502022993, 0.9800958462581607),
        (0.07553834749655068, 26.10103412529695, 68.9065407772954, 2.016129355318782),
        (0.05008650875532962, -0.6470589696141678, 3280.0479468729336, 4.642381720487216),
    ],
}

rows = []
for model_name, p0 in MODELS.items():
    n_pulses = len(p0)
    for window_name, (start, stop) in WINDOWS.items():
        mask = np.logical_and(t_full > start, t_full < stop)
        t_w, y_w = t_full[mask], y_full[mask]
        y_max = np.max(y_w)
        y_norm = y_w / y_max

        nf = NorrisFitter(t_w, y_norm)
        nf.fit(p0=p0)

        nf.plot_fit(
            show_individuals=True,
            x_label="Time since trigger [s]",
            y_label="Normalized count rate",
            data_label="10-400 keV NaI\nBackground Subtracted",
            title=f"GRB140206B {model_name} window-sensitivity check: {window_name} [{start}, {stop}] s",
        )
        fig_path = Path(__file__).parent / f"window_sensitivity_{model_name.lower()}_{window_name}"
        save_fig(plt.gcf(), fig_path)

        print(f"\n=== {model_name}, window {window_name} ({start}, {stop}), y_max={y_max:.2f} cts/s ===")
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
                {
                    "model": model_name,
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

# --- Side-by-side comparison per pulse, per model, across all three windows
WINDOW_ORDER = ["narrow_-1_160", "wide_-20_300", "widest_full_range"]
for model_name, p0 in MODELS.items():
    n_pulses = len(p0)
    print(f"\n=== {model_name}: narrow vs wide vs widest, per pulse ===")
    sub = df[df["model"] == model_name]
    by_window = {w: sub[sub["window"] == w].set_index("pulse_index") for w in WINDOW_ORDER}
    narrow = by_window["narrow_-1_160"]
    for i in range(1, n_pulses + 1):
        n = narrow.loc[i]
        parts = [f"pulse {i}: t_s {n.t_s:9.3f}"]
        t_v_ref = n.t_v_s
        for w in WINDOW_ORDER[1:]:
            row = by_window[w].loc[i]
            dt_v_pct = 100 * abs(row.t_v_s - t_v_ref) / t_v_ref
            parts.append(
                f" -> {row.t_s:9.3f} (t_v={row.t_v_s:.4f}, kept={row.mc_kept_fraction:.3f}, "
                f"dt_v={dt_v_pct:.1f}% vs narrow)"
            )
        print("".join(parts))

# --- Pulse 5 specifically (per the user's request): where does its t_s edge land across all three
# windows, for both models? Pulse 5 is SIMPLE's broad ~3440s-tau1 pedestal (already seen pinned to
# the wide window's left edge); in COMPLEX, pulse 5 is one of the three narrow 23-28s-region pulses
# (t_s~23.08, tau1~1.1) -- a different physical component, reported for completeness since the request
# was "both simple and complex". COMPLEX's own broad pedestal is pulse 7, reported alongside for context.
print("\n=== Pulse 5 t_s edge-tracking across all three windows ===")
for model_name in MODELS:
    sub = df[df["model"] == model_name]
    print(f"\n{model_name} pulse 5 (t_s per window):")
    for w in WINDOW_ORDER:
        start, stop = WINDOWS[w]
        row = sub[(sub["window"] == w) & (sub["pulse_index"] == 5)].iloc[0]
        print(f"  {w:22s} [{start:9.3f}, {stop:9.3f}]  t_s={row.t_s:10.3f}  tau1={row.tau1:10.3f}  kept={row.mc_kept_fraction:.3f}")
    if model_name == "COMPLEX":
        print(f"{model_name} pulse 7 (the actual broad pedestal, for context):")
        for w in WINDOW_ORDER:
            start, stop = WINDOWS[w]
            row = sub[(sub["window"] == w) & (sub["pulse_index"] == 7)].iloc[0]
            print(f"  {w:22s} [{start:9.3f}, {stop:9.3f}]  t_s={row.t_s:10.3f}  tau1={row.tau1:10.3f}  kept={row.mc_kept_fraction:.3f}")
