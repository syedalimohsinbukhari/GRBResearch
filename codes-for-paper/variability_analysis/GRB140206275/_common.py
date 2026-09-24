"""Shared data loading and result-building for GRB140206B's normalized-vs-unnormalized full-range Norris
fit, split into fitter_normalized.py and fitter_unnormalized.py. Same treatment, same structure, as
GRB080916C/_common.py -- see that file's docstring for the full rationale; fourth burst through the same
pipeline, adapted for GRB140206B's own episode bounds, P0, and (per that burst's own
fitter_GRB140206275.py) dominant-flux photon assignment -- same rule as GRB231129C/_common.py, not
GRB080916C's nearest-preceding-active-onset or GRB131014A's plain nearest-preceding-onset.

Not a package import across unrelated folders (CLAUDE.md's "copy rather than fight sys.path" is about
that) -- GRB140206275/ is a child of variability_analysis/, same precedent as the other three sibling
folders and experiments/normalized_vs_unnormalized_fit/fit_comparison.py before them.

T05/T95 verified directly against results.json (not copied from a comment without checking) --
GRB140206275's own "T90 7.488_154.240" entry, matching TR1's start and TR6's (the last TR) end exactly,
per CLAUDE.md's own T90 = first-TR-start-to-last-TR-end definition.

Model choice: COMPLEX (7-pulse), not SIMPLE (5-pulse) -- ../fitter_GRB140206275.py itself keeps both live
as separate candidate decompositions but uses COMPLEX for its own main per-episode CSV ("Use COMPLEX (the
finer-grained, user-preferred decomposition) for the main per-episode results" -- that file's own
comment), so this promotes the same choice, not a new one. P0 (COMPLEX_P0) and max_iterations=20000 (that
file itself omits max_iterations, defaulting to 5000 -- 20000 is carried over from every other burst's
full-range treatment and from experiments/normalized_vs_unnormalized_fit/fit_comparison.py's own
per-burst config, which already needed it at a widened window for this exact burst).

Known risk, not yet resolved here: fit_comparison.py's GRB140206275 entry deliberately used a *partial*
widening, (-20, 300), not the true full x.min()/x.max(), specifically because pulse 7 (the broad pedestal)
was landing with t_s pinned against the (-1, 160) window's own left edge -- NorrisFitter.fit_boundaries()
ties t_s's lower bound to the window's left edge, so a pinned t_s means the box constraint, not the data,
decided that parameter. This folder uses the TRUE full range instead, for consistency with the other
three bursts' treatment -- expect pulses 4/5/6 (the fragile 23-28s cluster, already flagged in
window_sensitivity_GRB140206275/ and full_range_all_bursts_check.py) to be the ones most likely affected,
not a new finding if so.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent  # .../variability_analysis/GRB140206275/
VARIABILITY_DIR = HERE.parent  # .../variability_analysis/
CODES_FOR_PAPER_DIR = VARIABILITY_DIR.parent  # .../codes-for-paper/
PROJECT_ROOT = CODES_FOR_PAPER_DIR.parent  # .../GRBResearchWork/
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from light_curves import lightcurve_data  # noqa: E402
from norris_fit import norris_pulse, t_peak, tv_mc_summary  # noqa: E402
from grb_research import update_style, MARKER_SIZE, LINE_WIDTH  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

update_style()

LC_DIR = PROJECT_ROOT / "light_curves"
GRB_140206 = LC_DIR / "GRB140206275"
GRB_PAPER_NAME = "GRB140206B"
ENERGY_LOW, ENERGY_HIGH = 10, 400
MAX_NFEV = 20000

# x_scale for the unnormalized method's least_squares call -- same values and same reasoning as the other
# three bursts' _common.py / experiments/normalized_vs_unnormalized_fit/fit_comparison.py.
X_SCALE_TS_TAU1 = 5.0
X_SCALE_TAU2 = 1.0

# Episode boundaries for GRB140206275, read directly from results.json (not eyeballed) -- copied from
# ../fitter_GRB140206275.py: EX0/TR1..TR6. No trailing excess (EX1) for this burst.
EPISODE_BOUNDS = {
    "EX0": (4.288, 11.072),
    "TR1": (7.488, 11.072),
    "TR2": (11.072, 20.032),
    "TR3": (20.032, 26.752),
    "TR4": (26.752, 59.776),
    "TR5": (59.776, 100.032),
    "TR6": (100.032, 154.240),
}

# T05/T95 (T90 start/end) -- plot-axis bounding only, does not change what was fitted. Verified directly
# against results.json's own "T90 7.488_154.240" entry (see module docstring).
T05, T95 = 7.488, 154.240
PLOT_PAD_S = 25.0

# COMPLEX 7-pulse seed, copied verbatim from ../fitter_GRB140206275.py's COMPLEX_P0 -- see module
# docstring for why COMPLEX, not SIMPLE.
P0 = [
    (0.13, -0.3, 0.12, 1.68),
    (0.35, 4.0, 6, 14),
    (0.5, 11, 5, 1.4),
    (0.5, 28, 2, 1),
    # (0.3, 19, 1, 1),
    (0.23, 24, 0.3, 1.434),
    (0.23, 30, 83, 1.1),
    (0.05, -0.9, 2900, 4.4),
]
N_PULSES = len(P0)

# LAT photon overlay -- every individual photon (energy, arrival time), no energy floor -- same
# convention as ../fitter_GRB140206275.py. T0_MET_S copied, not re-derived; see that file's own comment
# for the T_0 source line.
LAT_FITS = VARIABILITY_DIR / "GRB140206275_lat.fits"
T0_MET_S = 413361375.84


def _load_summed_light_curve():
    dat_NaI = sorted(f.stem for f in GRB_140206.glob("*.dat") if "n" in f.stem)
    nai_data = [lightcurve_data(str(GRB_140206 / f"{d}.dat"), ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]
    t = nai_data[0][0]
    for det, (t_i, _, _) in zip(dat_NaI[1:], nai_data[1:]):
        assert np.array_equal(t, t_i), f"{det}'s time grid differs from {dat_NaI[0]}'s -- cannot sum"
    r = np.sum([r for _, r, _ in nai_data], axis=0)
    b = np.sum([b for _, _, b in nai_data], axis=0)
    return t, r - b, dat_NaI


def _load_lat_photons():
    lat = fits.open(LAT_FITS)[1].data
    energy = np.asarray(lat["ENERGY"])
    t_arr = np.asarray(lat["TIME"]) - T0_MET_S
    order = np.argsort(t_arr)
    return energy[order], t_arr[order]


# Loaded once at import time, shared by both driver scripts -- same light curve, same photon list,
# regardless of fitting method. Full range: no window mask applied (was (-1, 160) in
# ../fitter_GRB140206275.py).
T_FULL, Y_RAW_FULL, DAT_NAI = _load_summed_light_curve()
Y_MAX_CTS_PER_S = float(np.max(Y_RAW_FULL))
FIT_WINDOW = (float(T_FULL.min()), float(T_FULL.max()))
PHOTON_ENERGY_MEV, PHOTON_T_ARR_S = _load_lat_photons()


class FitResult:
    """Minimal duck-typed stand-in for a NorrisFitter -- see GRB080916C/_common.py's own docstring for
    the same class; identical here, kept per-folder per this project's "copy rather than fight sys.path"
    convention for small, self-contained pieces."""

    def __init__(self, params, covariance, n_par=4):
        self.params = np.asarray(params)
        self.covariance = np.asarray(covariance)
        self.n_par = n_par


def assign_pulses(parameters, n_pulses):
    """Dominant-flux photon assignment: each photon goes to whichever pulse has the largest
    model-predicted flux at its arrival time -- same rule as ../fitter_GRB140206275.py and
    GRB231129C/_common.py (see either file's own comment for why the nearest-preceding-onset rule broke
    down on bursts with a fast pulse turning on just before a slower, still-dominant one peaks)."""

    def assign_one(t_arr: float):
        fluxes = [norris_pulse(np.array([t_arr]), parameters[i * 4:(i + 1) * 4])[0] for i in range(n_pulses)]
        return (int(np.argmax(fluxes)) + 1) if max(fluxes) > 0 else None  # 1-indexed

    return [assign_one(t) for t in PHOTON_T_ARR_S]


def photon_summary_for(photon_pulse, n_pulses):
    summary = {}
    for i in range(1, n_pulses + 1):
        idx = [k for k, pi in enumerate(photon_pulse) if pi == i]
        if not idx:
            summary[i] = {"n_photons": 0, "photon_e_max_MeV": None, "photon_t_arr_at_e_max_s": None}
            continue
        best = max(idx, key=lambda k: PHOTON_ENERGY_MEV[k])
        summary[i] = {
            "n_photons": len(idx),
            "photon_e_max_MeV": PHOTON_ENERGY_MEV[best],
            "photon_t_arr_at_e_max_s": PHOTON_T_ARR_S[best],
        }
    return summary


def build_results_df(method_label: str, fit_result: FitResult, n_pulses: int, y_max_cts_per_s: float,
                      amplitude_is_physical: bool, photon_summary: dict) -> pd.DataFrame:
    """One row per (pulse, episode) match, same schema as ../fitter_GRB140206275.py's CSV plus a
    `fit_method` column and an `A_norm` column derived either way (A/y_max) -- see
    GRB080916C/_common.py's build_results_df for the identical convention."""
    rows = []
    for i in range(n_pulses):
        A, ts, tau1, tau2 = fit_result.params[i * 4: (i + 1) * 4]
        A_cts_per_s = A if amplitude_is_physical else A * y_max_cts_per_s
        tp = t_peak(ts, tau1, tau2)
        mc = tv_mc_summary(fit_result, pulse_index=i + 1)
        matched = [name for name, (start, end) in EPISODE_BOUNDS.items() if start <= tp <= end]
        for episode in matched or [None]:
            rows.append({
                "grb_name": GRB_PAPER_NAME,
                "fit_method": method_label,
                "episode": episode,
                "pulse_index": i + 1,
                "n_pulses_total": n_pulses,
                "fit_window_start_s": FIT_WINDOW[0],
                "fit_window_end_s": FIT_WINDOW[1],
                "energy_low_keV": ENERGY_LOW,
                "energy_high_keV": ENERGY_HIGH,
                "detectors": "+".join(DAT_NAI),
                "A_norm": A_cts_per_s / y_max_cts_per_s,
                "A_cts_per_s": A_cts_per_s,
                "y_max_cts_per_s": y_max_cts_per_s,
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
            })
    return pd.DataFrame(rows)


def add_lat_photon_overlay(ax, y_top_data):
    """Twin-axis LAT photon overlay (red hollow circles, gridlines + matching minor ticks aligned), same
    convention as ../fitter_GRB140206275.py and GRB080916C/_common.py's own version (including the
    AutoMinorLocator fix applied there)."""
    from matplotlib.ticker import AutoMinorLocator

    ax_photon = ax.twinx()
    ax_photon.scatter(
        PHOTON_T_ARR_S, PHOTON_ENERGY_MEV,
        marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE ** 2,
        label="LAT photons",
    )
    ax_photon.set_ylabel("Photon energy [MeV]")

    n_yticks = 6
    y_buffer_frac = 0.05
    photon_top = np.round(PHOTON_ENERGY_MEV.max(), -3)

    ax.set_ylim(-y_buffer_frac * y_top_data, y_top_data * (1 + y_buffer_frac))
    ax_photon.set_ylim(-y_buffer_frac * photon_top, photon_top * (1 + y_buffer_frac))
    ax.set_yticks(np.linspace(0, y_top_data, n_yticks))
    ax_photon.set_yticks(np.linspace(0, photon_top, n_yticks))

    minor_ticks_per_major = 5
    ax.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    ax_photon.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    return ax_photon
