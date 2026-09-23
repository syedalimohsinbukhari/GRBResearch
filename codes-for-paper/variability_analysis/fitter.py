"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from grb_research import update_style, MARKER_SIZE
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import NorrisFitter, norris_pulse, t_peak, tv_mc_summary

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

# LAT photon overlay -- only photons with E > 1 GeV, per user request; same FITS-reading convention as
# fitter_GRB231129779.py (lat[1].data with ENERGY/TIME columns, per-source probability column named after
# the source). T0_MET_S is read directly from
# LAT_analysis/018__GRB080916009/Ep1A__m0.128_4.864/GRB080916C_fit_results_-0.128_4.864.txt's own "T_0"
# line -- not re-derived -- and cross-checked: the FITS file's own GTI (243216767.9-243216830.877) equals
# T0_MET_S+1.280 to T0_MET_S+64.257, i.e. exactly the T90 window (1.280-64.256s per results.json/EPISODE_BOUNDS
# above), confirming this is the right T_0 for this file.
LAT_FITS = Path(__file__).parent / "GRB080916009_lat.fits"  # copied from light_curves/GRB080916009/lat.fits, per CLAUDE.md's "copy rather than fight sys.path" convention
T0_MET_S = 243216766.62
PHOTON_E_MIN_MEV = 1000  # > 1 GeV

update_style()

dat = [f for f in os.listdir(f"{GRB_080916C}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = sorted(i for i in dat if "n" in i)  # order no longer load-bearing -- see BUG-23 fix below; kept sorted for deterministic logging only

nai_data = [lightcurve_data(f"{GRB_080916C}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

# BUG-23 fix (2026-09-23, user decision): sum all of the burst's NaI detectors' background-subtracted
# count rates raw, with no per-detector normalization -- each detector's own effective area/viewing angle
# is trusted to weight its own contribution, matching the one existing precedent in this project (the
# GRB080916C ROOT cross-check's "n3+n4 summed" light curve, variability_analysis.md). This also
# structurally closes BUG-23: os.listdir()[0]/sorted()[0] previously had to guess *which* single detector
# was "the" documented one (a guess that silently broke after a resync for three of the four bursts);
# summing every NaI detector in dat_NaI removes that pick entirely, so ordering can no longer matter.
# Detector time grids are confirmed identical before summing, not assumed.
t1 = nai_data[0][0]
for _det, (_t, _, _) in zip(dat_NaI[1:], nai_data[1:]):
    assert np.array_equal(t1, _t), f"{_det}'s time grid differs from {dat_NaI[0]}'s -- cannot sum"
r1 = np.sum([r for _, r, _ in nai_data], axis=0)
b1 = np.sum([b for _, _, b in nai_data], axis=0)

mask_ = np.logical_and(t1 > START1, t1 < END1)

y = (r1 - b1)[mask_]
Y_MAX_CTS_PER_S = np.max(y)  # kept for rescaling A back to physical units after the fit -- see below
y /= Y_MAX_CTS_PER_S


# max_iterations=20000 (was the NorrisFitter default of 5000) -- same reason as the GRB131014A fix: the
# summed 2-detector curve (peak ~3242 cts/s vs the single-detector ~1880) needs more iterations than the
# default budget under the same seed.
nf = NorrisFitter(t1[mask_], y, max_iterations=20000)

f, ax = plt.subplots(figsize=(12, 8))

# Reverted to a 6-pulse decomposition, 2026-09-23 (user call), dropping the 7th pulse (was A=0.1, t_s=20,
# tau1=9, tau2=1, sitting inside TR3's window) that the BUG-23 detector-summing fix above exposed as
# degenerate on the summed curve: it converged to tau2=0.0087 (a near-delta-function spike, not the broad
# feature it used to fit), with mc_kept_fraction collapsing from 0.52 to 0.05. Confirmed this wasn't a seed
# or iteration-budget artifact (both this file's original seed and one reconverged from the old
# single-detector fit landed on the same degenerate optimum, at any iteration budget tried up to 50000),
# and confirmed against the raw summed light curve that t~20.5s has no real standout feature -- the region
# sits at ~1650-1700 cts/s against a ~1000-1800 cts/s noisy floor throughout, i.e. the pulse was fitting
# two noisy bins, not resolving a genuine sub-feature. Dropping it was verified not to disturb the other
# six pulses: every other pulse's t_peak/t_v/A matches the 7-pulse fit's corresponding pulse to within
# ordinary t_s<->tau1-degeneracy noise (largest t_peak shift 0.03s, largest t_v shift 0.02s), and total SSE
# rises only ~0.4% (7.202e7 vs 7.173e7 physical-units counts/s^2) -- the signature of removing a pulse that
# was fitting noise, not a real feature. This is the same P0_6 already independently used as the baseline
# in experiments/window_sensitivity_GRB080916009/window_sensitivity.py (see "6-vs-7-pulse question" in
# variability_analysis.md for the original, pre-BUG-23, single-detector-curve history of that question --
# this reopens it on different grounds, not a re-litigation of the same finding).
P0 = [
    (0.4, -0.7, 2.34, 0.471),
    (0.7, 0.3, 0.88, 6),
    (0.2, 5.3, 0.6, 6),
    (0.3, 1.3, 43, 14),
    (0.3, 52, 18, 0.7),
    (0.3, 61, 0.3, 3),
]
nf.fit(p0=P0)

nf.plot_fit(show_individuals=True, axis=ax)
plt.show()

parameters = nf.params
covariance = nf.covariance
n_pulses = len(P0)

# --- LAT photon list: only E > 1 GeV, per user request -- not every photon (unlike
# fitter_GRB231129779.py's overlay, which plots the full list). Same FITS-reading convention as that file.
lat = fits.open(LAT_FITS)[1].data
photon_energy_MeV = np.asarray(lat["ENERGY"])
photon_t_arr_s = np.asarray(lat["TIME"]) - T0_MET_S
_e_mask = photon_energy_MeV > PHOTON_E_MIN_MEV
photon_energy_MeV, photon_t_arr_s = photon_energy_MeV[_e_mask], photon_t_arr_s[_e_mask]
_order = np.argsort(photon_t_arr_s)
photon_energy_MeV, photon_t_arr_s = photon_energy_MeV[_order], photon_t_arr_s[_order]

# Assign each photon to the pulse whose onset (t_s) it most recently follows, among pulses that are still
# "active" at the photon's arrival time -- nearest-preceding-onset (temporal-proximity), same rule as
# fitter_CLAUDE_GRB131014215.py, chosen over the dominant-flux rule (fitter_GRB231129779.py/
# fitter_GRB140206275.py) per user decision 2026-09-22: the two nearest E>1GeV photons here (6.072s,
# 6.857s) sit only 0.20s/0.99s after norris3's peak (5.869s) but 3.43s/4.21s after norris2's peak (2.646s),
# and fall inside TR2 (norris3's own episode) not TR1 (norris2's) -- dominant-flux picks norris2 anyway,
# since norris2's broad tail (tau2~6.17) still out-predicts norris3's small, narrow pulse at that time, the
# same "overshadowing" behavior already documented for norris2 against the (now-reverted) 8-pulse split.
# Temporal proximity favors norris3 for these two photons, and that's what was asked for here.
#
# ACTIVE_THRESHOLD_FRAC refinement, user decision 2026-09-22: the unrestricted rule (any t_s <= t_arr,
# ignoring whether that pulse still has anything to say) assigned every photon from 6.072s up to pulse 5's
# onset (19.729s) to pulse 3 -- including three (10.215s, 16.538s, 16.798s) that had already decayed to
# ~0% of norris3's own peak while sitting deep inside pulse 4's broad TR3 pedestal, because pulse 4's
# earlier onset (t_s=1.010s) can never "supersede" pulse 3's later one under pure onset-recency. A pulse
# is now only a candidate if its own predicted value is >= ACTIVE_THRESHOLD_FRAC of its own peak value at
# that instant -- among remaining (still-active) candidates, the nearest-preceding-onset tie-break is
# unchanged. 1% was picked for a wide margin, not fine-tuned: norris3's own value falls from 11.4% of its
# own peak at 7.445s to 0.076% at 10.215s, a >2-orders-of-magnitude drop, so the exact threshold value
# doesn't matter anywhere between roughly 0.1% and 10% -- the three far photons move to pulse 4 either way.
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

print(f"\n{'t_arr_s':>10}  {'E_MeV':>10}  {'pulse':>6}")
for t, e, pulse_i in zip(photon_t_arr_s, photon_energy_MeV, photon_pulse):
    print(f"{t:10.5f}  {e:10.2f}  {pulse_i!s:>6}")

# Per-pulse photon summary: count + the highest-energy assigned photon (mirrors the existing per-episode
# "defining photon" convention -- max-E is what feeds Gamma_min).
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
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_080916C.name}.csv"
results_df.to_csv(csv_path, index=False)
print(f"wrote {csv_path} ({len(results_df)} rows)")

# --- Plot
data_label = f"10-400 keV NaI ({'+'.join(dat_NaI)}, summed)\nBackground Subtracted"
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Count rate [counts/s]",
    data_label=data_label,
)
fig = plt.gcf()
ax = plt.gca()

# LAT photons on a twin y-axis: arrival time (shared x-axis) vs energy, red hollow circles --
# matches make_lightcurve.py's own red-hollow-circle convention for LAT photons, same as
# fitter_GRB231129779.py's overlay.
ax_photon = ax.twinx()
ax_photon.scatter(
    photon_t_arr_s, photon_energy_MeV,
    marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE**2,
    label="LAT photons (E > 1 GeV)",
)
ax_photon.set_ylabel("Photon energy [MeV]")

# Align both y-axes' zero and top and every gridline in between -- same approach as
# fitter_GRB231129779.py: force identical fractional tick positions so gridlines coincide
# regardless of each axis's own (different) scale.
N_YTICKS = 6
Y_BUFFER_FRAC = 0.05  # of each axis's own data max, both top and bottom -- keeps points off the plot edges
y_top = y.max()
photon_top = photon_energy_MeV.max()

ax.set_ylim(-Y_BUFFER_FRAC * y_top, y_top * (1 + Y_BUFFER_FRAC))
ax_photon.set_ylim(-Y_BUFFER_FRAC * photon_top, photon_top * (1 + Y_BUFFER_FRAC))
ax.set_yticks(np.linspace(0, y_top, N_YTICKS))
ax_photon.set_yticks(np.linspace(0, photon_top, N_YTICKS))

fig_path = Path(__file__).parent / f".norris_fitted_{GRB_080916C.name}"
save_fig(fig, fig_path)
