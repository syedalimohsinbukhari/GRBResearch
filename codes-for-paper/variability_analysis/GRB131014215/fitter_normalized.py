"""GRB131014A, NORMALIZED method, fit over the light curve's own full x.min()/x.max() range.

Same NorrisFitter/pymultifit path as ../fitter_CLAUDE_GRB131014215.py (peak-normalize y, fit, rescale A
back afterward), same 5-pulse P0, same nearest-preceding-onset LAT-photon assignment and episode-matching
logic -- only the fit window changed, from that file's production (-1, 10) to this burst's full range
(see _common.py's FIT_WINDOW). Same treatment as GRB080916C/fitter_normalized.py -- see that file's
docstring for the full rationale.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pymultifit.fitters.backend import BaseFitter
from scipy.signal import find_peaks

from _common import (  # noqa: E402 -- sets up sys.path for norris_fit below, must import first
    T_FULL, Y_RAW_FULL, Y_MAX_CTS_PER_S, FIT_WINDOW, DAT_NAI, GRB_131014, GRB_PAPER_NAME,
    N_PULSES, P0, MAX_NFEV, T05, T95, PLOT_PAD_S, assign_pulses, photon_summary_for, build_results_df,
    add_lat_photon_overlay, save_fig,
)
from norris_fit import NorrisFitter  # noqa: E402

y_norm = Y_RAW_FULL / Y_MAX_CTS_PER_S

nf: BaseFitter = NorrisFitter(T_FULL, y_norm, max_iterations=MAX_NFEV)

P0 = [
    (0.1, -1.2, 35, 0.07),
    (0.3, -2.57, 155, 0.1),
    (0.6, 1.2, 0.8, 0.3),
    (0.5, 2.4, 0.04, 1),
    (0.27, 2.66, 7, 0.1)
]

# nf.dry_run()
# plt.show()
nf.fit(p0=P0)


def seed_next_pulse_from_residual(t, residual, tau1_guess=1.0, tau2_guess=1.0, t_range=None):
    """Locate the most prominent unmodeled peak in `residual` and return a Norris
    (amplitude, t_s, tau1, tau2) seed whose t_peak = t_s + sqrt(tau1*tau2) lands on it.
    tau1/tau2 are left at 1 (this project's convention for a fresh, unrefined seed --
    see _common.py's P0) and t_s is solved backwards from the detected peak location.
    """
    mask = np.ones_like(t, dtype=bool) if t_range is None else (t >= t_range[0]) & (t <= t_range[1])
    t_masked, res_masked = t[mask], residual[mask]
    peaks, props = find_peaks(res_masked, prominence=0.05 * res_masked.max())
    if len(peaks) == 0:
        best = np.argmax(res_masked)
    else:
        best = peaks[np.argmax(props["prominences"])]
    t_peak = t_masked[best]
    amplitude = res_masked[best]
    t_s = t_peak - np.sqrt(tau1_guess * tau2_guess)
    return (amplitude, t_s, tau1_guess, tau2_guess), t_peak


# residual = y_norm - nf.get_fitted_curve()
# seed5, t_peak5 = seed_next_pulse_from_residual(T_FULL, residual, t_range=(T05 - 5, T95 + 10))
# print(f"residual-seeded 5th pulse: A={seed5[0]:.3f}, t_s={seed5[1]:.3f}, "
#       f"tau1={seed5[2]}, tau2={seed5[3]} -> t_peak={t_peak5:.3f}")

# P0 = P0 + [seed5]
nf.fit(p0=P0)
nf.plot_fit(show_individuals=True)
plt.xlim(-10, 10)
plt.show()

# photon_pulse = assign_pulses(nf.params, N_PULSES)
# photon_summary = photon_summary_for(photon_pulse, N_PULSES)
#
# results_df = build_results_df(
#     method_label="normalized", fit_result=nf, n_pulses=N_PULSES, y_max_cts_per_s=Y_MAX_CTS_PER_S,
#     amplitude_is_physical=False, photon_summary=photon_summary,
# )
# csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_131014.name}_normalized.csv"
# results_df.to_csv(csv_path, index=False)
# print(f"wrote {csv_path} ({len(results_df)} rows)")
# print(f"fit window: {FIT_WINDOW} (full x.min()/x.max())")
#
# # --- Plot
# fig, ax = plt.subplots(figsize=(13, 6.5))
# data_label = f"10-400 keV NaI ({'+'.join(DAT_NAI)}, summed)\nBackground Subtracted"
# nf.plot_fit(
#     show_individuals=True,
#     x_label="Time since trigger [s]",
#     y_label="Normalized count rate",
#     data_label=data_label,
#     title=f"{GRB_PAPER_NAME}: normalized fit, full range {FIT_WINDOW[0]:.1f}-{FIT_WINDOW[1]:.1f}s",
#     axis=ax,
# )
# add_lat_photon_overlay(ax, y_top_data=y_norm.max())
#
# # Plot-axis bounding only (T05-{PLOT_PAD_S}s .. T95+{PLOT_PAD_S}s) -- does not change what was fitted.
# ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)
#
# fig_path = Path(__file__).parent / f"norris_fitted_{GRB_131014.name}_normalized"
# save_fig(fig, fig_path)
# print(f"wrote {fig_path}.png / .pdf")
