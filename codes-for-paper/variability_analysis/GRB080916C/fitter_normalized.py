"""GRB080916C, NORMALIZED method, fit over the light curve's own full x.min()/x.max() range.

Same NorrisFitter/pymultifit path as ../fitter.py (peak-normalize y, fit, rescale A back afterward), same
6-pulse P0, same LAT-photon assignment and episode-matching logic -- only the fit window changed, from
../fitter.py's production (-1, 70) to this burst's full range (see _common.py's FIT_WINDOW). Promotes the
full-range finding already validated in
experiments/normalized_vs_unnormalized_fit/fit_comparison_results.csv (t_peak/t_v agree to <0.3% between
methods for 5 of 6 pulses) into a full per-burst output: CSV + PNG + PDF, MC-propagated t_v, LAT photon
overlay -- everything ../fitter.py produces, not just the point-estimate comparison fit_comparison.py
made.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from _common import (  # noqa: E402 -- sets up sys.path for norris_fit below, must import first
    T_FULL, Y_RAW_FULL, Y_MAX_CTS_PER_S, FIT_WINDOW, DAT_NAI, GRB_080916C, GRB_PAPER_NAME,
    N_PULSES, P0, MAX_NFEV, T05, T95, PLOT_PAD_S, assign_pulses, photon_summary_for, build_results_df,
    add_lat_photon_overlay, save_fig,
)
from norris_fit import NorrisFitter  # noqa: E402

y_norm = Y_RAW_FULL / Y_MAX_CTS_PER_S

nf = NorrisFitter(T_FULL, y_norm, max_iterations=MAX_NFEV)
nf.fit(p0=P0)

photon_pulse = assign_pulses(nf.params, N_PULSES)
photon_summary = photon_summary_for(photon_pulse, N_PULSES)

results_df = build_results_df(
    method_label="normalized", fit_result=nf, n_pulses=N_PULSES, y_max_cts_per_s=Y_MAX_CTS_PER_S,
    amplitude_is_physical=False, photon_summary=photon_summary,
)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_080916C.name}_normalized.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")
print(f"fit window: {FIT_WINDOW} (full x.min()/x.max())")

# --- Plot
fig, ax = plt.subplots(figsize=(13, 6.5))
data_label = f"10-400 keV NaI ({'+'.join(DAT_NAI)}, summed)\nBackground Subtracted"
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Normalized count rate",
    data_label=data_label,
    title=f"{GRB_PAPER_NAME}: normalized fit, full range {FIT_WINDOW[0]:.1f}-{FIT_WINDOW[1]:.1f}s",
    axis=ax
)
# fig = plt.gcf()
# ax = plt.gca()
add_lat_photon_overlay(ax, y_top_data=y_norm.max())

# Plot-axis bounding only (T05-25s .. T95+25s) -- does not change what was fitted, the fit itself still
# ran on the full FIT_WINDOW range above.
ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_080916C.name}_normalized"
save_fig(fig, fig_path)
print(f"wrote {fig_path}.png / .pdf")
