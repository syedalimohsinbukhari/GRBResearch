"""GRB140206B, NORMALIZED method, fit over the light curve's own full x.min()/x.max() range --
HISTORICAL RECORD ONLY, the P0 decomposition (7 pulses: near-t=0 pulse out, t=15s pulse in)
superseded by fitter_normalized.py's Q0 (also 7 pulses, but the other way around -- user decision,
2026-09-26; see _common.py's P0/Q0 comments for the numeric comparison). Kept on disk under
*_p0-suffixed filenames rather than deleted, per the user's request that both configurations remain
available for future reference. Suffix is "_p0", not "_7pulse", because both configurations have 7
pulses -- see fitter_normalized.py for the canonical (Q0) counterpart.

Otherwise identical to fitter_normalized.py -- same NorrisFitter/pymultifit path, same dominant-flux
LAT-photon assignment and episode-matching logic, same full-range fit window. Only the seed (P0 vs
Q0) and the output filenames differ.
"""
from pathlib import Path

import matplotlib.pyplot as plt

from _common import (  # noqa: E402 -- sets up sys.path for norris_fit below, must import first
    T_FULL, Y_RAW_FULL, Y_MAX_CTS_PER_S, FIT_WINDOW, DAT_NAI, GRB_140206, GRB_PAPER_NAME,
    N_PULSES_P0, P0, MAX_NFEV, T05, T95, PLOT_PAD_S, assign_pulses, photon_summary_for, build_results_df,
    add_lat_photon_overlay, save_fig,
)
from norris_fit import NorrisFitter  # noqa: E402

N_PULSES = N_PULSES_P0
y_norm = Y_RAW_FULL / Y_MAX_CTS_PER_S

nf = NorrisFitter(T_FULL, y_norm, max_iterations=MAX_NFEV)
nf.fit(p0=P0)

photon_pulse = assign_pulses(nf.params, N_PULSES)
photon_summary = photon_summary_for(photon_pulse, N_PULSES)

results_df = build_results_df(
    method_label="normalized", fit_result=nf, n_pulses=N_PULSES, y_max_cts_per_s=Y_MAX_CTS_PER_S,
    amplitude_is_physical=False, photon_summary=photon_summary,
)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_140206.name}_normalized_p0.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")
print(f"fit window: {FIT_WINDOW} (full x.min()/x.max())")

# --- Plot
fig, ax = plt.subplots(figsize=(13, 6.5))
data_label = f"10-400 keV NaI\n({'+'.join(DAT_NAI)})"
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Normalized count rate",
    data_label=data_label,
    title=" ",
    axis=ax,
)
add_lat_photon_overlay(ax, y_top_data=y_norm.max())

# Plot-axis bounding only (T05-{PLOT_PAD_S}s .. T95+{PLOT_PAD_S}s) -- does not change what was fitted.
ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_140206.name}_normalized_p0"
save_fig(fig, fig_path)
print(f"wrote {fig_path}.png / .pdf")
