"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import update_style
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import NorrisFitter, t_peak, tv_mc_summary

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_080916C = LC_DIR / "GRB080916009"
GRB_PAPER_NAME = "GRB080916C"  # for GRBPlotStyle colour + output naming; GRB_131014 above is a path, left as-is
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064
PAD_FRACTION = 0.20
PAD_MIN_S = 0.5
NULL_REL_ERR_THRESHOLD = 1.0  # >100% relative error on tau1 or tau2 => no_reliable_fit
T_ARR_MATCH_TOLERANCE_S = 1.0

START1 = -1
END1 = 70

EPISODE_BOUNDS = {
    "TR1": (1.280, 4.864),
    "TR2": (4.864, 15.040),
    "TR3": (15.040, 55.296),
    "TR4": (55.296, 59.52),
    "TR5": (59.52, 64.256),
}
BOUNDARY_COLORS = {"TR1": "tab:red", "TR2": "tab:green", "TR3": "tab:purple", "TX4": "tab:orange",
                   "TR5": "tab:brown"}

update_style()

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

# 8-pulse decomposition, recovered 2026-09-16 (the earlier hand-typed p0 above lost norris3/norris4 and
# duplicated norris2 while editing). norris2-1/norris2-2 split off norris2's overshadowing broad tail
# (tau2~6.7) into two sharper pulses that the data's own residual (6-pulse fit minus data) shows real,
# coherent structure at: +0.065..0.079 near t=1.41-1.54s and +0.072..0.086 near t=2.43-2.62s (separated by
# a -0.08..-0.10 dip near t=1.92-2.05s), confirmed not noise. These two were NOT reachable by fitting all
# 8 pulses jointly from raw hand-guesses -- norris2's tail absorbs the flux first and the optimizer
# collapses the new pulses to degenerate near-zero-width spikes (tau2~1e-3, kept_fraction=0). The fix:
# seed norris2-1/norris2-2 at the RESIDUAL's own amplitude (~0.08), not the raw light curve's amplitude
# (~0.9), and seed every other pulse from the already-converged 6-pulse fit (norris1/norris2/norris3/
# norris4/norris5/norris6 below), not fresh hand-guesses -- this converges cleanly.
# Neutral-seed reproducibility check (tau1/tau2 = 0.3/0.3 for both new pulses instead of the values below):
# norris2-1 reproduces to t_peak=1.480s, t_v=0.049s (vs. 1.466s/0.051s here) -- reads as real, same
# criterion used for GRB131014A's pulses 3/4/5. norris2-2 reproduces to t_peak=2.462s, t_v=0.177s (vs
# 2.309s/0.148s here) -- a looser match, "probably real, not yet as tight as norris2-1" -- one more seed
# trial before trusting its exact numbers, per that same GRB131014A precedent.
P0 = [
    (0.4304, -0.6943, 2.1212, 0.4932),   # norris1: precursor
    (0.7757, 0.3243, 0.8285, 6.7078),    # norris2: broad envelope/tail (the one overshadowing 3/4/5)
    (0.08, 1.35, 0.15, 0.15),            # norris2-1: NEW, residual-seeded, t_peak lands ~1.47s
    (0.08, 2.15, 0.20, 0.30),            # norris2-2: NEW, residual-seeded, t_peak lands ~2.31s
    (0.2565, 5.3055, 0.6752, 0.4701),    # norris3: the ~5.7-5.9s peak
    (0.3556, 4.7780, 30.6423, 14.3745),  # norris4: broad TR3 pedestal
    (0.1511, 54.1091, 13.0910, 0.8433),  # norris5
    (0.1539, 61.5605, 0.3919, 2.9220),   # norris6
]
nf.fit(p0=P0)

# nf.fit(
#     p0=[
#         (0.4, -0.7, 2.34, 0.471),
#         (0.7, 0.3, 0.88, 6),
#         (0.2, 5.3, 0.6, 6),
#         (0.3, 1.3, 43, 14),
#         (0.1, 20, 9, 1),
#         (0.3, 52, 18, 0.7),
#         (0.3, 61, 0.3, 3),
#     ]
# )

nf.plot_fit(show_individuals=True, axis=ax)
plt.show()

parameters = nf.params
covariance = nf.covariance
n_pulses = len(P0)

# --- Results CSV: one row per (pulse, episode) match.
# A pulse can legitimately belong to more than one episode -- EX0 fully contains TR1's window and EX1 fully contains
# TR2's, per CLAUDE.md's EX0/EX1 definition -- so matching is "does this pulse's t_peak fall inside this episode's
# boundary", not a one-to-one assignment.
# t_v is MC-propagated through the fit's own covariance via tv_mc_summary(), the same convention as every other Norris fit in this project.
rows = []
for i in range(n_pulses):
    A, ts, tau1, tau2 = parameters[i * 4: (i + 1) * 4]
    tp = t_peak(ts, tau1, tau2)
    mc = tv_mc_summary(nf, pulse_index=i + 1)
    matched = [name for name, (start, end) in EPISODE_BOUNDS.items() if start <= tp <= end]
    for episode in matched or [None]:
        rows.append(
            {
                "grb_name": GRB_PAPER_NAME,
                "episode": episode,
                "pulse_index": i + 1,
                "n_pulses_total": n_pulses,
                "fit_window_start_s": START1,
                "fit_window_end_s": END1,
                "energy_low_keV": ENERGY_LOW,
                "energy_high_keV": ENERGY_HIGH,
                "A_norm": A,
                "A_cts_per_s": A * Y_MAX_CTS_PER_S,
                "y_max_cts_per_s": Y_MAX_CTS_PER_S,
                "t_s": ts,
                "tau1": tau1,
                "tau2": tau2,
                "t_peak_s": tp,
                "t_v_s": mc["t_v_s"],
                "t_v_err_lower_s": mc["t_v_err_lower_s"],
                "t_v_err_upper_s": mc["t_v_err_upper_s"],
                "t_v_definition": (
                    "t_v=(tau2/2)*sqrt((ln2+2*sqrt(tau1/tau2))^2-4*tau1/tau2), "
                    "Bukhari et al. 2022 Adv.Space Res. eq.10, attributed to Norris et al. 2005; = FWHM/2"
                ),
                "mc_kept_fraction": mc["kept_fraction"],
                "n_samples": mc["n_samples"],
                "seed": mc["seed"],
            }
        )

results_df = pd.DataFrame(rows)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_080916C.name}.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")

# --- Plot
data_label = f"10-400 keV NaI\nBackground Subtracted"
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Count rate [counts/s]",
    data_label=data_label,
)
fig = plt.gcf()
ax = plt.gca()

fig_path = Path(__file__).parent / f".norris_fitted_{GRB_080916C.name}"
save_fig(fig, fig_path)
