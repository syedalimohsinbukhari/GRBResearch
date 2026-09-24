"""GRB140206B, NORMALIZED method, fit over the light curve's own full x.min()/x.max() range.

Same NorrisFitter/pymultifit path as ../fitter_GRB140206275.py (peak-normalize y, fit, rescale A back
afterward), same 7-pulse COMPLEX P0, same dominant-flux LAT-photon assignment and episode-matching logic
-- only the fit window changed, from that file's production (-1, 160) to this burst's true full range
(see _common.py's FIT_WINDOW and its module docstring's "Known risk" note about pulse 7). Same treatment
as GRB080916C/fitter_normalized.py -- see that file's docstring for the full rationale.

Not run yet, per instruction -- write only.
"""
from pathlib import Path

import matplotlib.pyplot as plt

from _common import (  # noqa: E402 -- sets up sys.path for norris_fit below, must import first
    T_FULL, Y_RAW_FULL, Y_MAX_CTS_PER_S, FIT_WINDOW, DAT_NAI, GRB_140206, GRB_PAPER_NAME,
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
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_140206.name}_normalized.csv"
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
    title=f"{GRB_PAPER_NAME}: normalized fit (COMPLEX), full range {FIT_WINDOW[0]:.1f}-{FIT_WINDOW[1]:.1f}s",
    axis=ax,
)
add_lat_photon_overlay(ax, y_top_data=y_norm.max())

# Plot-axis bounding only (T05-{PLOT_PAD_S}s .. T95+{PLOT_PAD_S}s) -- does not change what was fitted.
ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_140206.name}_normalized"
save_fig(fig, fig_path)
print(f"wrote {fig_path}.png / .pdf")
