"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from variability_timescale.light_curves import lightcurve_data
from variability_timescale.norris_fit import NorrisFitter

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_080916C = LC_DIR / "GRB140206275"
GRB_PAPER_NAME = "GRB140206B"  # for GRBPlotStyle colour + output naming; GRB_131014 above is a path, left as-is
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064
PAD_FRACTION = 0.20
PAD_MIN_S = 0.5
NULL_REL_ERR_THRESHOLD = 1.0  # >100% relative error on tau1 or tau2 => no_reliable_fit
T_ARR_MATCH_TOLERANCE_S = 1.0

START1 = -1
END1 = 160

# Episode boundaries for GRB131014215, read directly from results.json (not eyeballed):
# "EX0 -0.192_2.432", "TR1 0.960_2.432", "TR2 2.432_4.160", "EX1 2.432_6.976".
# EPISODE_BOUNDS = {
#     "EX0": (-0.192, 2.432),
#     "TR1": (0.960, 2.432),
#     "TR2": (2.432, 4.160),
#     "EX1": (2.432, 6.976),
# }
# BOUNDARY_COLORS = {"EX0": "tab:red", "TR1": "tab:green", "TR2": "tab:purple", "EX1": "tab:orange"}

# update_style()

dat = [f for f in os.listdir(f"{GRB_080916C}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = [i.split(".")[0] for i in dat if "n" in i]

nai_data = [lightcurve_data(f"{GRB_080916C}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

t1, r1, b1 = nai_data[0]
# t2, r2, b2 = nai_data[1]

mask_ = np.logical_and(t1 > START1, t1 < END1)

y = (r1 - b1)[mask_]
Y_MAX_CTS_PER_S = np.max(y)  # kept for rescaling A back to physical units after the fit -- see below
y /= Y_MAX_CTS_PER_S

nf = NorrisFitter(t1[mask_], y)

f, ax = plt.subplots(figsize=(12, 8))

nf.fit(
    p0=[
        (0.13, -0.3, 0.12, 1.68),
        (0.35, 4., 6, 14),
        (0.5, 11, 5, 1.4),
        (0.5, 28, 2, 1),
        (0.23, 24, 0.3, 1.434),
        (0.23, 23, 83, 1.1),
        (0.05, -0.9, 3400, 4.4),
    ]
)

# nf.fit(
#     p0=[
#         (0.13, -0.3, 0.12, 1.68),
#         (0.35, 4., 6, 14),
#         (0.5, 11, 5, 1.4),
#         (0.23, 26, 2.5, 3.7),
#         (0.05, -0.96, 3400, 4.4),
#     ]
# )


nf.plot_fit(show_individuals=True, axis=ax)
# plt.savefig(f"norris_fit_{GRB_080916C.name}.png", dpi=SAVE_DPI)
# plt.close()
plt.show()

# parameters = nf.params
# covariance = nf.covariance
# n_pulses = 5
#
# # --- Results CSV: one row per (pulse, episode) match.
# # A pulse can legitimately belong to more than one episode -- EX0 fully contains TR1's window and EX1 fully contains
# # TR2's, per CLAUDE.md's EX0/EX1 definition -- so matching is "does this pulse's t_peak fall inside this episode's
# # boundary", not a one-to-one assignment.
# # t_v is MC-propagated through the fit's own covariance via tv_mc_summary(), the same convention as every other Norris fit in this project.
# rows = []
# for i in range(n_pulses):
#     A, ts, tau1, tau2 = parameters[i * 4: (i + 1) * 4]
#     tp = t_peak(ts, tau1, tau2)
#     mc = tv_mc_summary(nf, pulse_index=i + 1)
#     matched = [name for name, (start, end) in EPISODE_BOUNDS.items() if start <= tp <= end]
#     for episode in matched or [None]:
#         rows.append(
#             {
#                 "grb_name": GRB_PAPER_NAME,
#                 "episode": episode,
#                 "pulse_index": i + 1,
#                 "n_pulses_total": n_pulses,
#                 "fit_window_start_s": START1,
#                 "fit_window_end_s": END1,
#                 "energy_low_keV": ENERGY_LOW,
#                 "energy_high_keV": ENERGY_HIGH,
#                 "A_norm": A,
#                 "A_cts_per_s": A * Y_MAX_CTS_PER_S,
#                 "y_max_cts_per_s": Y_MAX_CTS_PER_S,
#                 "t_s": ts,
#                 "tau1": tau1,
#                 "tau2": tau2,
#                 "t_peak_s": tp,
#                 "t_v_s": mc["t_v_s"],
#                 "t_v_err_lower_s": mc["t_v_err_lower_s"],
#                 "t_v_err_upper_s": mc["t_v_err_upper_s"],
#                 "t_v_definition": (
#                     "t_v=(tau2/2)*sqrt((ln2+2*sqrt(tau1/tau2))^2-4*tau1/tau2), "
#                     "Bukhari et al. 2022 Adv.Space Res. eq.10, attributed to Norris et al. 2005; = FWHM/2"
#                 ),
#                 "mc_kept_fraction": mc["kept_fraction"],
#                 "n_samples": mc["n_samples"],
#                 "seed": mc["seed"],
#             }
#         )
#
# results_df = pd.DataFrame(rows)
# csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_080916C.name}.csv"
# results_df.to_csv(csv_path, index=False)
# print(f"wrote {csv_path} ({len(results_df)} rows)")
#
# # --- Plot
# data_label = f"10-400 keV NaI\nBackground Subtracted"
# nf.plot_fit(
#     show_individuals=True,
#     x_label="Time since trigger [s]",
#     y_label="Count rate [counts/s]",
#     data_label=data_label,
# )
# fig = plt.gcf()
# ax = plt.gca()
#
# fig_path = Path(__file__).parent / f".norris_fitted_{GRB_080916C.name}"
# save_fig(fig, fig_path)
