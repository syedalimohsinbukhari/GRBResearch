"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from grb_research import update_style, seed_from_name, MARKER_SIZE
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import NorrisFitter, t_peak, tv_mc_summary

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_131014 = LC_DIR / "GRB131014215"
GRB_PAPER_NAME = "GRB131014A"  # for GRBPlotStyle color + output naming; GRB_131014 above is a path, left as-is
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064
PAD_FRACTION = 0.20
PAD_MIN_S = 0.5
NULL_REL_ERR_THRESHOLD = 1.0  # >100% relative error on tau1 or tau2 => no_reliable_fit
T_ARR_MATCH_TOLERANCE_S = 1.0

START1 = -1
END1 = 10

# Episode boundaries for GRB131014215, read directly from results.json (not eyeballed):
# "EX0 -0.192_2.432", "TR1 0.960_2.432", "TR2 2.432_4.160", "EX1 2.432_6.976".
EPISODE_BOUNDS = {
    "EX0": (-0.192, 2.432),
    "TR1": (0.960, 2.432),
    "TR2": (2.432, 4.160),
    "EX1": (2.432, 6.976),
}
BOUNDARY_COLORS = {"EX0": "tab:red", "TR1": "tab:green", "TR2": "tab:purple", "EX1": "tab:orange"}

update_style()

dat = [f for f in os.listdir(f"{GRB_131014}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = [i.split(".")[0] for i in dat if "n" in i]

nai_data = [lightcurve_data(f"{GRB_131014}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

t1, r1, b1 = nai_data[0]
# t2, r2, b2 = nai_data[1]

mask_ = np.logical_and(t1 > START1, t1 < END1)

y = (r1 - b1)[mask_]
Y_MAX_CTS_PER_S = np.max(y)  # kept for rescaling A back to physical units after the fit -- see below
y /= Y_MAX_CTS_PER_S

nf = NorrisFitter(t1[mask_], y)
# nf = NorrisFitter(t1, (r1 - b1), max_iterations=10000)
# nf.dry_run()
# plt.show()
nf.fit(
    p0=[
        (0.115, -0.9, 1, 1),
        (0.2, -0.8, 1, 1),
        (0.9, 1.19, 1, 1),
        (0.4, 2.4, 1, 1),
        (0.3, 2.71, 1, 1),
    ]
)

# actual scale fit
# nf.fit(
#     p0=[
#         (4960, -0.86, 3.2, 1.15),
#         (9644, -0.9, 179, 0.02),
#         (37.3e3, 1.233, 0.756, 0.357),
#         (18400, 2.411, 0.161, 1.066),
#         (6200, 2.6, 6, 0.1),
#     ]
# )

parameters = nf.params
covariance = nf.covariance
n_pulses = 5

# Deterministic, per-script seed -- project convention (SEEDING.md): derived from this file's own
# name + MASTER_SEED, never a bare literal, so it's reproducible without hardcoding.
SEED = seed_from_name(__file__)

# --- LAT photon list:
# every individual photon (energy, arrival time), not just the one "defining" (max-energy, TS-significant) photon per
# episode already used elsewhere (lorentz_results.csv, EPISODE_PHOTONS in variability_timescale.py).
# FITS-reading convention copied from light_curves/make_lightcurve.py (not re-derived) -- lat[1].data with ENERGY/TIME
# columns and a per-source probability column named after the source (here "GRB131014A", matching the FITS file's own
# TTYPE51).
# T0_MET_S is read directly from LAT_analysis/014__GRB131014215/Ep1__0.960_2.432/GRB131014A_fit_results_0.96_2.432.txt's
# own "T_0, 403420143.2" line -- not re-derived -- and cross-checked: the FITS file's own GTI (403420144.16-403420147.36)
# equals T0_MET_S+0.96 to T0_MET_S+4.16, i.e., exactly the T90 window (0.960-4.160s per results.json),
# confirming this is the right T_0 for this file.
LAT_FITS = Path(__file__).parent / "GRB131014215_lat.fits"  # copied from light_curves/GRB131014215/lat.fits
T0_MET_S = 403420143.2

lat = fits.open(LAT_FITS)[1].data
photon_energy_MeV = np.asarray(lat["ENERGY"])
photon_t_arr_s = np.asarray(lat["TIME"]) - T0_MET_S
photon_prob = np.asarray(lat[GRB_PAPER_NAME])
_order = np.argsort(photon_t_arr_s)
photon_energy_MeV, photon_t_arr_s, photon_prob = photon_energy_MeV[_order], photon_t_arr_s[_order], photon_prob[_order]

# Assign each photon to the pulse whose onset (t_s) it most recently follows -- the largest t_s that is still <= the
# photon's arrival time. Physically motivated (a photon belongs to whichever pulse has "turned on" and not yet been
# superseded by the next pulse's onset), not arbitrary, and it reproduces the manual identification already verified
# against this data: the 445.28 MeV photon (t=2.533s) assigns to pulse 4 (t_s~2.47 is the largest t_s <= 2.533, since
# pulse 5's t_s~2.66-2.76 hasn't happened yet), and the 552.79 MeV photon (t=4.111s) assigns to pulse 5.
# All seven TR1-window photons (t=1.85-2.01s) assign to pulse 3 this way.
pulse_ts_values = [parameters[i * 4 + 1] for i in range(n_pulses)]


def assign_pulse(t_arr: float):
    candidates = [i for i, ts in enumerate(pulse_ts_values) if ts <= t_arr]
    return (max(candidates, key=lambda i: pulse_ts_values[i]) + 1) if candidates else None  # 1-indexed


photon_pulse = [assign_pulse(t) for t in photon_t_arr_s]

print(f"\n{'t_arr_s':>10}  {'E_MeV':>10}  {'prob':>7}  {'pulse':>6}")
for t, e, p, pulse_i in zip(photon_t_arr_s, photon_energy_MeV, photon_prob, photon_pulse):
    print(f"{t:10.5f}  {e:10.2f}  {p:7.4f}  {pulse_i!s:>6}")

# Per-pulse photon summary: count + the highest-energy assigned photon (mirrors the existing
# per-episode "defining photon" convention -- max-E is what feeds Gamma_min).
photon_summary = {}
for i in range(1, n_pulses + 1):
    idx = [k for k, pi in enumerate(photon_pulse) if pi == i]
    if not idx:
        photon_summary[i] = {"n_photons": 0, "photon_e_max_MeV": None, "photon_t_arr_at_e_max_s": None}
        continue
    best = max(idx, key=lambda k: photon_energy_MeV[k])
    photon_summary[i] = {
        "n_photons": len(idx),
        "photon_e_max_MeV": photon_energy_MeV[best],
        "photon_t_arr_at_e_max_s": photon_t_arr_s[best],
    }

# --- Results CSV: one row per (pulse, episode) match.
# A pulse can legitimately belong to more than one episode -- EX0 fully contains TR1's window and EX1 fully contains
# TR2's, per CLAUDE.md's EX0/EX1 definition -- so matching is "does this pulse's t_peak fall inside this episode's
# boundary", not a one-to-one assignment.
# t_v is MC-propagated through the fit's own covariance via tv_mc_summary(), the same convention as every other Norris fit in this project.
rows = []
for i in range(n_pulses):
    A, ts, tau1, tau2 = parameters[i * 4: (i + 1) * 4]
    tp = t_peak(ts, tau1, tau2)
    mc = tv_mc_summary(nf, pulse_index=i + 1, seed=SEED)
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
                "n_lat_photons_assigned": photon_summary[i + 1]["n_photons"],
                "photon_e_max_MeV": photon_summary[i + 1]["photon_e_max_MeV"],
                "photon_t_arr_at_e_max_s": photon_summary[i + 1]["photon_t_arr_at_e_max_s"],
            }
        )

results_df = pd.DataFrame(rows)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_131014.name}.csv"
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

# LAT photons on a twin y-axis: arrival time (shared x-axis) vs energy, red hollow circles --
# matches make_lightcurve.py's own red-hollow-circle convention for LAT photons.
ax_photon = ax.twinx()
ax_photon.scatter(
    photon_t_arr_s, photon_energy_MeV,
    marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE ** 2,
    label="LAT photons",
)
ax_photon.set_ylabel("Photon energy [MeV]")

# Align both y-axes' zero and top *and* every gridline in between: (0, data_max) alone isn't
# enough, since matplotlib's default locator still picks a different number of "nice" ticks
# per axis independently (e.g. 6 on the left, 7 on the right), so gridlines wouldn't actually
# coincide. Forcing identical fractional tick positions (N evenly-spaced ticks from 0 to each
# axis's own max) guarantees every gridline lines up horizontally, regardless of whether the
# resulting right-axis numbers are "round".
N_YTICKS = 6
Y_BUFFER_FRAC = 0.05  # of each axis's own data max, both top and bottom -- keeps points off the plot edges
y_top = y.max()
photon_top = photon_energy_MeV.max()

# Ticks stay at the *unpadded* 0..top fractional positions -- only the axis limits (the
# viewing window) grow by Y_BUFFER_FRAC of each axis's own top. Since the padding fraction is
# identical on both axes, every tick still lands at the same fractional height on both,
# regardless of their different absolute scales (counts/s vs MeV).
ax.set_ylim(-Y_BUFFER_FRAC * y_top, y_top * (1 + Y_BUFFER_FRAC))
ax_photon.set_ylim(-Y_BUFFER_FRAC * photon_top, photon_top * (1 + Y_BUFFER_FRAC))
ax.set_yticks(np.linspace(0, y_top, N_YTICKS))
ax_photon.set_yticks(np.linspace(0, photon_top, N_YTICKS))

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_131014.name}"
save_fig(fig, fig_path)
