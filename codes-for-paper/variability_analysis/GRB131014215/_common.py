"""Shared data loading and result-building for GRB131014A's normalized-vs-unnormalized full-range Norris
fit, split into fitter_normalized.py and fitter_unnormalized.py (2026-09-24). Same treatment, same
structure, as GRB080916C/_common.py -- see that file's docstring for the full rationale; third burst
through the same pipeline, adapted for GRB131014A's own episode bounds, P0, and (per that burst's own
fitter_CLAUDE_GRB131014215.py) the simpler nearest-preceding-onset photon assignment -- no
ACTIVE_THRESHOLD_FRAC "still active" refinement, that was introduced later specifically for GRB080916C
and was never backported to this burst's own production script, so this doesn't add it either.

Not a package import across unrelated folders (CLAUDE.md's "copy rather than fight sys.path" is about
that) -- GRB131014215/ is a child of variability_analysis/, same precedent as GRB080916C/_common.py,
GRB231129C/_common.py, and experiments/normalized_vs_unnormalized_fit/fit_comparison.py before it.

T05/T95 verified directly against results.json (not copied from a comment without checking) --
GRB131014215's own "T90 0.960_4.160" entry, matching TR1's start and TR2's end exactly, per CLAUDE.md's
own T90 = first-TR-start-to-last-TR-end definition.

P0 (the 5-pulse seed) and max_iterations=20000 are copied verbatim from ../fitter_CLAUDE_GRB131014215.py
-- see that file's own comment for why 20000 is needed (the summed 3-detector curve, peak ~130,030 cts/s,
needs more iterations than the single-detector-tuned default of 5000).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent  # .../variability_analysis/GRB131014215/
VARIABILITY_DIR = HERE.parent  # .../variability_analysis/
CODES_FOR_PAPER_DIR = VARIABILITY_DIR.parent  # .../codes-for-paper/
PROJECT_ROOT = CODES_FOR_PAPER_DIR.parent  # .../GRBResearchWork/
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from light_curves import lightcurve_data  # noqa: E402
from norris_fit import t_peak, tv_mc_summary  # noqa: E402
from grb_research import update_style, MARKER_SIZE, LINE_WIDTH  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

update_style()

LC_DIR = PROJECT_ROOT / "light_curves"
GRB_131014 = LC_DIR / "GRB131014215"
GRB_PAPER_NAME = "GRB131014A"
ENERGY_LOW, ENERGY_HIGH = 10, 400
MAX_NFEV = 20000

# x_scale for the unnormalized method's least_squares call -- same values and same reasoning as
# GRB080916C/_common.py / experiments/normalized_vs_unnormalized_fit/fit_comparison.py.
X_SCALE_TS_TAU1 = 5.0
X_SCALE_TAU2 = 1.0

# Episode boundaries for GRB131014215, read directly from results.json (not eyeballed) -- copied from
# ../fitter_CLAUDE_GRB131014215.py: "EX0 -0.192_2.432", "TR1 0.960_2.432", "TR2 2.432_4.160", "EX1 2.432_6.976".
EPISODE_BOUNDS = {
    "EX0": (-0.192, 2.432),
    "TR1": (0.960, 2.432),
    "TR2": (2.432, 4.160),
    "EX1": (2.432, 6.976),
}

# T05/T95 (T90 start/end) -- plot-axis bounding only, does not change what was fitted. Verified directly
# against results.json's own "T90 0.960_4.160" entry (not copied from another file's comment without
# checking -- see module docstring).
T05, T95 = 0.960, 4.160
PLOT_PAD_S = 25.0

# 5-pulse seed, copied verbatim from ../fitter_CLAUDE_GRB131014215.py.
P0 = [
    (0.115, -0.9, 1, 1),
    (0.2, -0.8, 1, 1),
    (0.9, 1.19, 1, 1),
    (0.4, 2.4, 1, 1),
    (0.3, 2.71, 1, 1),
]
N_PULSES = len(P0)

# LAT photon overlay -- every individual photon (energy, arrival time), no energy floor -- same
# convention as ../fitter_CLAUDE_GRB131014215.py. T0_MET_S copied, not re-derived; see that file's own
# comment for the GTI cross-check that validated it.
LAT_FITS = VARIABILITY_DIR / "GRB131014215_lat.fits"
T0_MET_S = 403420143.2


def _load_summed_light_curve():
    dat_NaI = sorted(f.stem for f in GRB_131014.glob("*.dat") if "n" in f.stem)
    nai_data = [lightcurve_data(str(GRB_131014 / f"{d}.dat"), ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]
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
# regardless of fitting method. Full range: no window mask applied (was (-1, 10) in
# ../fitter_CLAUDE_GRB131014215.py).
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
    """Nearest-preceding-onset photon assignment -- each photon goes to whichever pulse has the largest
    t_s that is still <= the photon's arrival time, no "still active" filter. Same rule as
    ../fitter_CLAUDE_GRB131014215.py (copied verbatim, not GRB080916C's later ACTIVE_THRESHOLD_FRAC
    refinement -- that was introduced specifically for GRB080916C and never backported here). Unit
    agnostic: only t_s (always seconds) is used, amplitude never enters."""
    pulse_ts_values = [parameters[i * 4 + 1] for i in range(n_pulses)]

    def assign_one(t_arr: float):
        candidates = [i for i, ts in enumerate(pulse_ts_values) if ts <= t_arr]
        return (max(candidates, key=lambda i: pulse_ts_values[i]) + 1) if candidates else None  # 1-indexed

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
    """One row per (pulse, episode) match, same schema as ../fitter_CLAUDE_GRB131014215.py's CSV plus a
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
    convention as ../fitter_CLAUDE_GRB131014215.py and GRB080916C/_common.py's own version (including the
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
    photon_top = PHOTON_ENERGY_MEV.max()

    ax.set_ylim(-y_buffer_frac * y_top_data, y_top_data * (1 + y_buffer_frac))
    ax_photon.set_ylim(-y_buffer_frac * photon_top, photon_top * (1 + y_buffer_frac))
    ax.set_yticks(np.linspace(0, y_top_data, n_yticks))
    ax_photon.set_yticks(np.linspace(0, photon_top, n_yticks))

    minor_ticks_per_major = 5
    ax.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    ax_photon.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    return ax_photon
