"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from grb_research import update_style, seed_from_name, MARKER_SIZE
from grb_research.grb_utils import save_fig
from variability_timescale.light_curves import lightcurve_data
from variability_timescale.norris_fit import NorrisFitter, norris_pulse, t_peak, tv_mc_summary

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_140206 = LC_DIR / "GRB140206275"
GRB_PAPER_NAME = "GRB140206B"  # for GRBPlotStyle color + output naming; GRB_140206 above is a path, left as-is
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064

START1 = -1
END1 = 160

# Episode boundaries for GRB140206275, read directly from results.json (not eyeballed):
# "EX0 4.288_11.072", "TR1 7.488_11.072", "TR2 11.072_20.032", "TR3 20.032_26.752",
# "TR4 26.752_59.776", "TR5 59.776_100.032", "TR6 100.032_154.240". No trailing excess (EX1)
# for this burst -- matches CLAUDE.md's note that GRB140206B lacks one.
EPISODE_BOUNDS = {
    "EX0": (4.288, 11.072),
    "TR1": (7.488, 11.072),
    "TR2": (11.072, 20.032),
    "TR3": (20.032, 26.752),
    "TR4": (26.752, 59.776),
    "TR5": (59.776, 100.032),
    "TR6": (100.032, 154.240),
}

update_style()

dat = [f for f in os.listdir(f"{GRB_140206}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = [i.split(".")[0] for i in dat if "n" in i]

nai_data = [lightcurve_data(f"{GRB_140206}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

t1, r1, b1 = nai_data[0]

mask_ = np.logical_and(t1 > START1, t1 < END1)

y_raw = (r1 - b1)[mask_]
Y_MAX_CTS_PER_S = np.max(y_raw)  # kept for rescaling A back to physical units after the fit
y = y_raw / Y_MAX_CTS_PER_S

SEED = seed_from_name(__file__)

# --- Two candidate decompositions, both fitted by the user in fitter.py before this file
# existed: SIMPLE_P0 (5 pulses) treats the t~23-28s region as one broad pulse; COMPLEX_P0
# (7 pulses) splits that same region into three (t_s~23, 24, 28).
#
# User's actual reason for COMPLEX, not overfitting for its own sake: give the 753.11 MeV
# photon at t=23.998s (TR3's own defining photon, already in lorentz_results.csv) its own
# dedicated peak (pulse 5 here) instead of leaving it buried in a single broad pulse's tail.
# The diagnostic below shows this costs real reliability in that region (kept drops from
# SIMPLE's ~1.0 to COMPLEX's 0.54-0.67) -- a deliberate, known tradeoff (photon-anchored
# resolution vs. numerical stability), not a mistake. Both models are kept as separate files
# (this one, and fitter_GRB140206275_simple.py) rather than picking one as "the" answer.
SIMPLE_P0 = [
    (0.13, -0.3, 0.12, 1.68),
    (0.35, 4.0, 6, 14),
    (0.5, 11, 5, 1.4),
    (0.23, 26, 2.5, 3.7),
    (0.05, -0.96, 3400, 4.4),
]
COMPLEX_P0 = [
    (0.13, -0.3, 0.12, 1.68),
    (0.35, 4.0, 6, 14),
    (0.5, 11, 5, 1.4),
    (0.5, 28, 2, 1),
    (0.23, 24, 0.3, 1.434),
    (0.23, 23, 83, 1.1),
    (0.05, -0.9, 3400, 4.4),
]

nf_simple = NorrisFitter(t1[mask_], y)
nf_simple.fit(p0=SIMPLE_P0)

nf_complex = NorrisFitter(t1[mask_], y)
nf_complex.fit(p0=COMPLEX_P0)


def summarize(nf, n_pulses, label):
    print(f"\n=== {label} ({n_pulses} pulses) ===")
    rows = []
    for i in range(n_pulses):
        A, ts, tau1, tau2 = nf.params[i * 4 : (i + 1) * 4]
        tp = t_peak(ts, tau1, tau2)
        mc = tv_mc_summary(nf, pulse_index=i + 1, seed=SEED)
        print(
            f"  pulse {i + 1}: A={A:.4f} ts={ts:9.3f} tau1={tau1:9.3f} tau2={tau2:8.3f}  "
            f"t_peak={tp:9.3f}  t_v={mc['t_v_s']:.4f} +{mc['t_v_err_upper_s']:.4f} "
            f"-{mc['t_v_err_lower_s']:.4f}  kept={mc['kept_fraction']:.3f}"
        )
        rows.append(
            {"pulse_index": i + 1, "A": A, "t_s": ts, "tau1": tau1, "tau2": tau2, "t_peak_s": tp, **mc}
        )
    return rows


rows_simple = summarize(nf_simple, len(SIMPLE_P0), "SIMPLE")
rows_complex = summarize(nf_complex, len(COMPLEX_P0), "COMPLEX")

# --- Diagnostic: match each SIMPLE pulse to its nearest-t_peak COMPLEX counterpart(s) and
# compare t_v/kept_fraction. A simple pulse whose t_peak sits between several close complex
# pulses is exactly the "does the split earn its place" question -- reported, not resolved
# here, since that's a judgement call same as the GRB231129C case.
diag_rows = []
for rs in rows_simple:
    nearest = min(rows_complex, key=lambda rc: abs(rc["t_peak_s"] - rs["t_peak_s"]))
    diag_rows.append(
        {
            "simple_pulse": rs["pulse_index"],
            "simple_t_peak_s": rs["t_peak_s"],
            "simple_t_v_s": rs["t_v_s"],
            "simple_kept": rs["kept_fraction"],
            "nearest_complex_pulse": nearest["pulse_index"],
            "complex_t_peak_s": nearest["t_peak_s"],
            "complex_t_v_s": nearest["t_v_s"],
            "complex_kept": nearest["kept_fraction"],
            "t_peak_offset_s": abs(nearest["t_peak_s"] - rs["t_peak_s"]),
        }
    )
diag_df = pd.DataFrame(diag_rows)
diag_path = Path(__file__).parent / f"norris_fit_diagnostics_{GRB_140206.name}.csv"
diag_df.to_csv(diag_path, index=False)
print(f"\n=== SIMPLE vs COMPLEX diagnostic (nearest-t_peak match) ===")
print(diag_df.to_string(index=False))
print(f"wrote {diag_path}")

# --- Use COMPLEX (the finer-grained, user-preferred decomposition) for the main per-episode
# results, same schema as the other two fitter_*.py files in this folder.
nf = nf_complex
parameters = nf.params
n_pulses = len(COMPLEX_P0)

# --- LAT photon list: dominant-flux assignment (see fitter_GRB231129779.py for why this
# replaced the earlier nearest-preceding-onset rule -- that rule misattributed photons to a
# fast pulse turning on just before a slower, still-dominant earlier pulse peaks).
LAT_FITS = Path(__file__).parent / "GRB140206275_lat.fits"  # copied from light_curves/GRB140206275/GRB140206Bfiltered_gti_gtsrcprob_7.488_154.176.fits
T0_MET_S = 413361375.84  # LAT_analysis/007__GRB140206275/Ep1__7.488_11.072/*_fit_results_*.txt's own "T_0" line

lat = fits.open(LAT_FITS)[1].data
photon_energy_MeV = np.asarray(lat["ENERGY"])
photon_t_arr_s = np.asarray(lat["TIME"]) - T0_MET_S
photon_prob = np.asarray(lat[GRB_PAPER_NAME])
_order = np.argsort(photon_t_arr_s)
photon_energy_MeV, photon_t_arr_s, photon_prob = photon_energy_MeV[_order], photon_t_arr_s[_order], photon_prob[_order]


def assign_pulse(t_arr: float):
    fluxes = [norris_pulse(np.array([t_arr]), parameters[i * 4 : (i + 1) * 4])[0] for i in range(n_pulses)]
    return (int(np.argmax(fluxes)) + 1) if max(fluxes) > 0 else None


photon_pulse = [assign_pulse(t) for t in photon_t_arr_s]

print(f"\n{'t_arr_s':>10}  {'E_MeV':>10}  {'prob':>7}  {'pulse':>6}")
for t, e, p, pulse_i in zip(photon_t_arr_s, photon_energy_MeV, photon_prob, photon_pulse):
    print(f"{t:10.5f}  {e:10.2f}  {p:7.4f}  {pulse_i!s:>6}")

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

# --- Results CSV: one row per (pulse, episode) match, COMPLEX model only.
rows = []
for i in range(n_pulses):
    A, ts, tau1, tau2 = parameters[i * 4 : (i + 1) * 4]
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
csv_path = Path(__file__).parent / f"norris_fit_results_{GRB_140206.name}.csv"
results_df.to_csv(csv_path, index=False)
print(f"\nwrote {csv_path} ({len(results_df)} rows)")

# --- Plot (COMPLEX model)
data_label = f"10-400 keV NaI\nBackground Subtracted"
nf.plot_fit(
    show_individuals=True,
    x_label="Time since trigger [s]",
    y_label="Count rate [counts/s]",
    data_label=data_label,
)
fig = plt.gcf()
ax = plt.gca()

ax_photon = ax.twinx()
ax_photon.scatter(
    photon_t_arr_s, photon_energy_MeV,
    marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE**2,
    label="LAT photons",
)
ax_photon.set_ylabel("Photon energy [MeV]")

N_YTICKS = 6
Y_BUFFER_FRAC = 0.05
y_top = y.max()
photon_top = photon_energy_MeV.max()

ax.set_ylim(-Y_BUFFER_FRAC * y_top, y_top * (1 + Y_BUFFER_FRAC))
ax_photon.set_ylim(-Y_BUFFER_FRAC * photon_top, photon_top * (1 + Y_BUFFER_FRAC))
ax.set_yticks(np.linspace(0, y_top, N_YTICKS))
ax_photon.set_yticks(np.linspace(0, photon_top, N_YTICKS))

fig_path = Path(__file__).parent / f"norris_fitted_{GRB_140206.name}"
save_fig(fig, fig_path)
