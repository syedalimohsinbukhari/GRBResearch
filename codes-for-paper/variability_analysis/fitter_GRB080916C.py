"""Created on Sep 24 2026 -- GRB080916C, full x.min()/x.max() range, normalized method.

Standalone, self-contained copy of GRB080916C/fitter_normalized.py's current behavior, promoted to this
folder's standard per-burst fitter_*.py convention (own script, own light-curve loading, no import from
GRB080916C/_common.py -- CLAUDE.md's "copy rather than fight sys.path": that folder's shared module was
built for two sibling drivers living together, not for reaching into from outside it).

Same NorrisFitter/pymultifit path, same 6-pulse P0, same LAT-photon assignment and episode-matching logic
as fitter.py in this folder -- the only change is the fit window, now the light curve's own full
x.min()/x.max() instead of fitter.py's production (-1, 70), and the plot axis bounded to
T05-{PLOT_PAD_S}s..T95+{PLOT_PAD_S}s so the (much wider) fit window doesn't make the actual burst
illegible. Already validated: experiments/normalized_vs_unnormalized_fit/fit_comparison_results.csv and
GRB080916C/norris_fit_results_GRB080916009_normalized.csv show this converges cleanly and agrees with the
unnormalized method to <0.3% on t_peak/t_v for 5 of 6 pulses (pulse 5/TR4 is the known exception -- see
that folder's discussion, a pre-existing tau1<->tau2 degeneracy, not new here).

Outputs are suffixed `_fullrange` so they never collide with fitter.py's own production
(-1, 70)-window CSV/plot in this same folder.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.ticker import AutoMinorLocator

from grb_research import update_style, MARKER_SIZE
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import NorrisFitter, norris_pulse, t_peak, tv_mc_summary

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_080916C = LC_DIR / "GRB080916009"
GRB_PAPER_NAME = "GRB080916C"
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064

EPISODE_BOUNDS = {
    "TR1": (1.280, 4.864),
    "TR2": (4.864, 15.040),
    "TR3": (15.040, 55.296),
    "TR4": (55.296, 59.52),
    "TR5": (59.52, 64.256),
}

# T05/T95 (T90 start/end) -- plot-axis bounding only, does not change what is fitted (that always runs on
# the light curve's own full x.min()/x.max(), computed below from the data itself). Same values as
# GRB080916C/_common.py.
T05, T95 = 1.280, 64.256
PLOT_PAD_S = 23.0

# LAT photon overlay -- only photons with E > 1 GeV, same convention and same T0_MET_S as fitter.py
# (copied, not re-derived; see that file's own comment for the GTI cross-check that validated this T_0).
LAT_FITS = Path(__file__).parent / "GRB080916009_lat.fits"
T0_MET_S = 243216766.62
PHOTON_E_MIN_MEV = 1000  # > 1 GeV

update_style()

dat = [f for f in os.listdir(f"{GRB_080916C}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = sorted(i for i in dat if "n" in i)  # order no longer load-bearing, kept for deterministic logging only

nai_data = [lightcurve_data(f"{GRB_080916C}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

# BUG-23 fix convention (see fitter.py): sum all the burst's NaI detectors' background-subtracted count
# rates raw, no per-detector normalization. Detector time grids confirmed identical before summing.
t1 = nai_data[0][0]
for _det, (_t, _, _) in zip(dat_NaI[1:], nai_data[1:]):
    assert np.array_equal(t1, _t), f"{_det}'s time grid differs from {dat_NaI[0]}'s -- cannot sum"
r1 = np.sum([r for _, r, _ in nai_data], axis=0)
b1 = np.sum([b for _, _, b in nai_data], axis=0)

# Full x.min()/x.max() range -- was START1=-1, END1=70 in fitter.py. No mask needed: every point in the
# summed light curve is used.
START1, END1 = float(t1.min()), float(t1.max())

y = r1 - b1
Y_MAX_CTS_PER_S = np.max(y)  # kept for rescaling A back to physical units after the fit
y = y / Y_MAX_CTS_PER_S

# max_iterations=20000 (was the NorrisFitter default of 5000) -- same reason as fitter.py: the
# summed-detector curve needs more iterations than the default budget under the same seed.
nf = NorrisFitter(t1, y, max_iterations=20000)

# Same 6-pulse decomposition as fitter.py (2026-09-23 revert) -- already independently confirmed to
# converge cleanly at the full range (experiments/normalized_vs_unnormalized_fit/fit_comparison_results.csv,
# GRB080916C/norris_fit_results_GRB080916009_normalized.csv), so no new seed search was needed here.
P0 = [
    (0.4, -0.8, 2.7, 0.46),
    (0.8, 0.3, 0.7, 7),
    (0.2, 4.8, 3.12, 0.3),
    (0.3, 5.3, 28, 14),
    (0.1, 52, 85, 0.4),
    (0.1, 61, 0.6, 2.5),
]
nf.fit(p0=P0)

parameters = nf.params
covariance = nf.covariance
n_pulses = len(P0)

# --- LAT photon list: only E > 1 GeV -- same convention as fitter.py.
lat = fits.open(LAT_FITS)[1].data
photon_energy_MeV = np.asarray(lat["ENERGY"])
photon_t_arr_s = np.asarray(lat["TIME"]) - T0_MET_S
_e_mask = photon_energy_MeV > PHOTON_E_MIN_MEV
photon_energy_MeV, photon_t_arr_s = photon_energy_MeV[_e_mask], photon_t_arr_s[_e_mask]
_order = np.argsort(photon_t_arr_s)
photon_energy_MeV, photon_t_arr_s = photon_energy_MeV[_order], photon_t_arr_s[_order]

# Nearest-preceding-active-onset (temporal-proximity) photon assignment -- same rule, same ACTIVE_THRESHOLD_FRAC,
# as fitter.py. See that file's own comment for the full derivation history.
ACTIVE_THRESHOLD_FRAC = 0.01
pulse_params = [tuple(parameters[i * 4: (i + 1) * 4]) for i in range(n_pulses)]
pulse_peak_values = [
    norris_pulse(np.array([t_peak(ts, tau1, tau2)]), p)[0] for p, (A, ts, tau1, tau2) in zip(pulse_params, pulse_params)
]


def assign_pulse(t_arr: float):
    candidates = []
    for i, p in enumerate(pulse_params):
        ts = p[1]
        if ts > t_arr:
            continue
        val = norris_pulse(np.array([t_arr]), p)[0]
        if pulse_peak_values[i] > 0 and (val / pulse_peak_values[i]) >= ACTIVE_THRESHOLD_FRAC:
            candidates.append(i)
    return (max(candidates, key=lambda i: pulse_params[i][1]) + 1) if candidates else None  # 1-indexed


photon_pulse = [assign_pulse(t) for t in photon_t_arr_s]

# print(f"\n{'t_arr_s':>10}  {'E_MeV':>10}  {'pulse':>6}")
# for t, e, pulse_i in zip(photon_t_arr_s, photon_energy_MeV, photon_pulse):
#     print(f"{t:10.5f}  {e:10.2f}  {pulse_i!s:>6}")

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

# --- Results CSV: one row per (pulse, episode) match. Suffixed _fullrange -- see module docstring.
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
                "n_lat_photons_assigned": photon_summary[i + 1]["n_photons"],
                "photon_e_max_MeV": photon_summary[i + 1]["photon_e_max_MeV"],
                "photon_t_arr_at_e_max_s": photon_summary[i + 1]["photon_t_arr_at_e_max_s"],
            }
        )

results_df = pd.DataFrame(rows)
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_080916C.name}_fullrange.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")

# --- Plot
data_label = f"10-400 keV NaI\n({'+'.join(dat_NaI)})"
fig, ax = plt.subplots(figsize=(12, 7.5))
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Count rate [counts/s]",
    data_label=data_label,
    title=" ",
    axis=ax,
)

# LAT photons on a twin y-axis: arrival time (shared x-axis) vs energy, red hollow circles -- same
# convention as fitter.py.
ax_photon = ax.twinx()
ax_photon.scatter(
    photon_t_arr_s, photon_energy_MeV,
    marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE ** 2,
    label="LAT photons (E > 1 GeV)",
)
ax_photon.set_ylabel("Photon energy [MeV]")

# Align both y-axes' zero and top and every gridline in between -- same approach as fitter.py.
N_YTICKS = 5
Y_BUFFER_FRAC = 0.25
y_top = y.max()
photon_top = np.round(photon_energy_MeV.max(), -4)

ax.set_ylim(-Y_BUFFER_FRAC * y_top, y_top * (1 + Y_BUFFER_FRAC))
ax_photon.set_ylim(-Y_BUFFER_FRAC * photon_top, photon_top * (1 + Y_BUFFER_FRAC))
ax.set_yticks(np.linspace(0, y_top, N_YTICKS))
ax_photon.set_yticks(np.linspace(0, photon_top, N_YTICKS))

# Explicit, matching minor-tick subdivision for both axes -- AutoMinorLocator's default heuristic picks a
# subdivision count per axis from its own major-tick step size (0.25 here vs 7500 on the photon axis), so
# the two twinned axes' minor ticks/gridlines drift out of alignment even though the 5 major ticks line
# up. Forcing the same explicit n on both keeps them visually matched.
MINOR_TICKS_PER_MAJOR = 5
ax.yaxis.set_minor_locator(AutoMinorLocator(MINOR_TICKS_PER_MAJOR))
ax_photon.yaxis.set_minor_locator(AutoMinorLocator(MINOR_TICKS_PER_MAJOR))

# Plot-axis bounding only (T05-{PLOT_PAD_S}s .. T95+{PLOT_PAD_S}s) -- does not change what was fitted,
# the fit itself ran on the full (START1, END1) range above.
ax.set_xlim(T05 - PLOT_PAD_S, T95 + PLOT_PAD_S)

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_080916C.name}_fullrange"
save_fig(fig, fig_path)
print(f"wrote {fig_path}.png / .pdf")
