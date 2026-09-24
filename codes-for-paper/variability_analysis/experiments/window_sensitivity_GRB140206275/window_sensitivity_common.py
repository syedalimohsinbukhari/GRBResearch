r"""
Shared light-curve loading and per-model fitting/plotting logic for GRB140206B's window-sensitivity check, split out of
the combined window_sensitivity.py (2026-09-24) so the SIMPLE and COMPLEX decompositions can be run as two separate,
independent scripts (window_sensitivity_simple.py, window_sensitivity_complex.py) instead of one script always running
both.

Same window-widening pulse-count diagnostic as experiments/window_sensitivity_GRB231129779/ and
experiments/window_sensitivity_GRB131014215/ -- see this folder's original window_sensitivity.py docstring (kept below,
unedited) for the full rationale; only the driving control flow moved.

---

Experiment: same window-widening pulse-count diagnostic as experiments/window_sensitivity_GRB231129779/ and
experiments/window_sensitivity_GRB131014215/, applied here to GRB140206B -- the one burst in this project's
manual-Norris-fit track with two live candidate decompositions (fitter_GRB140206275_simple.py's 5-pulse SIMPLE model and
fitter_GRB140206275.py's 7-pulse COMPLEX model), so both are checked.

Mechanism under test:

NorrisFitter.fit_boundaries() (codes-for-paper/variability_analysis/norris_fit.py:49-50) sets t_s's lower bound to
x_values.min() -- i.e., the fit window's own left edge is t_s's box constraint.
Widening the window directly widens how far the optimizer can push t_s along the t_s<->tau1 degeneracy.

Controlled comparison: SAME p0 per model (that model's own narrow-window converged parameters, read
from norris_fit_results_GRB140206275.csv / norris_fit_results_GRB140206275_simple.csv) fit to the SAME
light curve, with ONLY the window bound changed, across THREE window widths. Narrow window matches
fitter_GRB140206275.py's own default (START1=-1, END1=160, covering all seven episodes through TR6's
154.240s end); wide window extends both edges well beyond that; widest pushes both edges all the way to
the light curve's own np.min(t)/np.max(t) -- the most extreme box NorrisFitter.fit_boundaries() can ever
present, added specifically to track where SIMPLE pulse 5's t_s (already seen pinned to the wide window's
left edge) ends up when the box constraint is removed as far as the data allows.

The COMPLEX model's own docstring in fitter_GRB140206275.py already flags that splitting the t~23-28s region into three
pulses (COMPLEX pulses 4/5/6) costs reliability relative to SIMPLE's single broad pulse there (kept ~1.0 for SIMPLE's
pulse 4 vs 0.54-0.67 for COMPLEX's 4/5/6) -- a deliberate, photon-anchored-resolution-vs-stability tradeoff, not a
mistake.
This diagnostic checks whether that already-poor COMPLEX reliability collapses further under window-widening (real
trouble) or holds (a legitimate tradeoff, not evidence of a wrong pulse count).
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

dat_NaI = sorted(f.split(".")[0] for f in os.listdir(LC_DIR) if f.endswith(".dat") and "n" in f)  # order no longer load-bearing -- see BUG-23 fix below; kept sorted for deterministic logging only

# BUG-23 fix (2026-09-23, user decision): sum all of the burst's NaI detectors' background-subtracted
# count rates raw, with no per-detector normalization -- matching fitter_GRB140206275.py's own fix and the
# existing "n3+n4 summed" ROOT cross-check precedent (variability_analysis.md). This structurally closes
# BUG-23: there's no longer a single detector to mis-pick via os.listdir()[0]. Detector time grids are
# confirmed identical before summing, not assumed.
nai_data = [lightcurve_data(f"{LC_DIR}/{d}.dat", ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]
t_full = nai_data[0][0]
for _det, (_t, _, _) in zip(dat_NaI[1:], nai_data[1:]):
    assert np.array_equal(t_full, _t), f"{_det}'s time grid differs from {dat_NaI[0]}'s -- cannot sum"
r1 = np.sum([r for _, r, _ in nai_data], axis=0)
b1 = np.sum([b for _, _, b in nai_data], axis=0)
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
WINDOW_ORDER = ["narrow_-1_160", "wide_-20_300", "widest_full_range"]


def run_model(model_name: str, p0: list, also_track: list[int] | None = None) -> pd.DataFrame:
    """Fit exactly one model (SIMPLE or COMPLEX) across all three windows, write that model's own
    results CSV + one plot per window, print the narrow/wide/widest per-pulse comparison, and the
    pulse-5 t_s edge-tracking (plus any extra pulses in `also_track`, e.g. COMPLEX's pulse 7 pedestal).
    """
    n_pulses = len(p0)
    rows = []
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
            A, ts, tau1, tau2 = nf.params[i * 4: (i + 1) * 4]
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
    csv_path = Path(__file__).parent / f"window_sensitivity_results_{model_name.lower()}.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path} ({len(df)} rows)")

    # --- Side-by-side comparison per pulse, across all three windows
    print(f"\n=== {model_name}: narrow vs wide vs widest, per pulse ===")
    by_window = {w: df[df["window"] == w].set_index("pulse_index") for w in WINDOW_ORDER}
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

    # --- Pulse 5 t_s edge-tracking across all three windows, plus any extra pulses in also_track
    # (e.g. COMPLEX's pulse 7, the actual broad pedestal there -- see this module's docstring/the
    # original window_sensitivity.py for why pulse 5 means different physical objects in each model).
    print(f"\n=== Pulse 5 t_s edge-tracking across all three windows ({model_name}) ===")
    print(f"\n{model_name} pulse 5 (t_s per window):")
    for w in WINDOW_ORDER:
        start, stop = WINDOWS[w]
        row = df[(df["window"] == w) & (df["pulse_index"] == 5)].iloc[0]
        print(f"  {w:22s} [{start:9.3f}, {stop:9.3f}]  t_s={row.t_s:10.3f}  tau1={row.tau1:10.3f}  kept={row.mc_kept_fraction:.3f}")
    for pulse_idx in also_track or []:
        print(f"{model_name} pulse {pulse_idx} (for context):")
        for w in WINDOW_ORDER:
            start, stop = WINDOWS[w]
            row = df[(df["window"] == w) & (df["pulse_index"] == pulse_idx)].iloc[0]
            print(f"  {w:22s} [{start:9.3f}, {stop:9.3f}]  t_s={row.t_s:10.3f}  tau1={row.tau1:10.3f}  kept={row.mc_kept_fraction:.3f}")

    return df
