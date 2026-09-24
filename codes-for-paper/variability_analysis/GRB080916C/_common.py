"""Shared data loading and result-building for GRB080916C's normalized-vs-unnormalized full-range Norris
fit, split into fitter_normalized.py and fitter_unnormalized.py (2026-09-24).

Promotes experiments/normalized_vs_unnormalized_fit/fit_comparison.py's GRB080916C full-range finding
(SSE agrees to <0.002% between the two methods, t_peak/t_v agree to <0.3% for 5 of 6 pulses -- see that
folder's fit_comparison_results.csv) out of the experiments/ sandbox and into a proper per-burst output
with the full production machinery: episode matching, LAT-photon-to-pulse assignment, MC-propagated t_v,
CSV + PNG + PDF -- the same outputs ../fitter.py produces, just fit over the light curve's own full
x.min()/x.max() range instead of the production (-1, 70) window, and for both methods rather than one.

Not a package import across unrelated folders (CLAUDE.md's "copy rather than fight sys.path" is about
that) -- GRB080916C/ is a child of variability_analysis/, and experiments/normalized_vs_unnormalized_fit/
fit_comparison.py already set the precedent of reaching into variability_analysis/ for its light_curves.py
and norris_fit.py via sys.path.insert rather than re-copying them a third time; this follows the same
pattern one level closer.

P0 (the 6-pulse seed) is copied verbatim from ../fitter.py, same reasoning as that file's own comment:
the 7th pulse (t_s~20, inside TR3) was dropped 2026-09-23 after it converged to a degenerate
near-delta-function spike on the summed-detector curve -- see BUGS.md BUG-23 and variability_analysis.md.
Already independently confirmed to converge cleanly at the full range in both methods
(fit_comparison_results.csv), so no new seed search was needed here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.ticker import AutoMinorLocator

HERE = Path(__file__).resolve().parent  # .../variability_analysis/GRB080916C/
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
GRB_080916C = LC_DIR / "GRB080916009"
GRB_PAPER_NAME = "GRB080916C"
ENERGY_LOW, ENERGY_HIGH = 10, 400
MAX_NFEV = 20000

# x_scale for the unnormalized method's least_squares call -- see fit_comparison.py's own docstring for
# why an explicit array is needed (x_scale='jac' was tried and rejected there) and the split between
# t_s/tau1 (~5s natural scale) and tau2 (~1s), copied verbatim.
X_SCALE_TS_TAU1 = 5.0
X_SCALE_TAU2 = 1.0

EPISODE_BOUNDS = {
    "TR1": (1.280, 4.864),
    "TR2": (4.864, 15.040),
    "TR3": (15.040, 55.296),
    "TR4": (55.296, 59.52),
    "TR5": (59.52, 64.256),
}

# T05/T95 (T90 start/end), same values as experiments/normalized_vs_unnormalized_fit/
# full_range_single_burst.py's T90 dict -- plot-axis bounding only, does not change what was fitted. The
# fit itself still runs on the full FIT_WINDOW range regardless; only the drivers' ax.set_xlim() uses this.
T05, T95 = 1.280, 64.256
PLOT_PAD_S = 23.0

# Same 6-pulse decomposition as ../fitter.py (2026-09-23 revert) -- see module docstring.
P0 = [
    (0.4, -0.7, 2.34, 0.471),
    (0.7, 0.3, 0.88, 6),
    (0.2, 5.3, 0.6, 6),
    (0.3, 1.3, 43, 14),
    (0.3, 52, 18, 0.7),
    (0.3, 61, 0.3, 3),
]
N_PULSES = len(P0)

# LAT photon overlay -- only E > 1 GeV, same convention and same T0_MET_S as ../fitter.py (copied, not
# re-derived; see that file's own comment for the GTI cross-check that validated this T_0). The FITS file
# itself is read from its existing location in the parent folder rather than copied a second time -- a
# data file, not code, so the "copy rather than fight sys.path" convention (about import gymnastics)
# doesn't apply the same way LC_DIR itself is referenced from PROJECT_ROOT without being copied.
LAT_FITS = VARIABILITY_DIR / "GRB080916009_lat.fits"
T0_MET_S = 243216766.62
PHOTON_E_MIN_MEV = 1000  # > 1 GeV

ACTIVE_THRESHOLD_FRAC = 0.01  # see ../fitter.py's own comment for the derivation of this value


class FitResult:
    """Minimal duck-typed stand-in for a NorrisFitter, so norris_fit.tv_mc_summary()
    (which only reads .params/.covariance/.n_par) works identically for the unnormalized method's bare
    scipy.optimize.least_squares result as it does for the normalized method's real NorrisFitter instance
    -- one MC-propagation implementation for both, not two."""

    def __init__(self, params, covariance, n_par=4):
        self.params = np.asarray(params)
        self.covariance = np.asarray(covariance)
        self.n_par = n_par


def _load_summed_light_curve():
    dat_NaI = sorted(f.stem for f in GRB_080916C.glob("*.dat") if "n" in f.stem)
    nai_data = [lightcurve_data(str(GRB_080916C / f"{d}.dat"), ENERGY_LOW, ENERGY_HIGH) for d in dat_NaI]
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
    mask = energy > PHOTON_E_MIN_MEV
    energy, t_arr = energy[mask], t_arr[mask]
    order = np.argsort(t_arr)
    return energy[order], t_arr[order]


# Loaded once at import time, shared by both driver scripts -- same light curve, same photon list,
# regardless of fitting method. Full range: no window mask applied (was (-1, 70) in ../fitter.py).
T_FULL, Y_RAW_FULL, DAT_NAI = _load_summed_light_curve()
Y_MAX_CTS_PER_S = float(np.max(Y_RAW_FULL))
FIT_WINDOW = (float(T_FULL.min()), float(T_FULL.max()))
PHOTON_ENERGY_MEV, PHOTON_T_ARR_S = _load_lat_photons()


def assign_pulses(parameters, n_pulses):
    """Nearest-preceding-active-onset photon-to-pulse assignment, same rule and same
    ACTIVE_THRESHOLD_FRAC as ../fitter.py (copied verbatim -- see that file's own comment for the full
    reasoning and the 2026-09-22 refinement history). Unit-agnostic: the amplitude A only appears inside
    a per-pulse ratio (val / pulse_peak_values[i]), so this works identically whether `parameters` are
    peak-normalized (fractional) or physical (counts/s) -- t_s/tau1/tau2 are always in seconds either way.
    """
    pulse_params = [tuple(parameters[i * 4:(i + 1) * 4]) for i in range(n_pulses)]
    pulse_peak_values = [
        norris_pulse(np.array([t_peak(ts, tau1, tau2)]), p)[0]
        for p, (A, ts, tau1, tau2) in zip(pulse_params, pulse_params)
    ]

    def assign_one(t_arr: float):
        candidates = []
        for i, p in enumerate(pulse_params):
            ts = p[1]
            if ts > t_arr:
                continue
            val = norris_pulse(np.array([t_arr]), p)[0]
            if pulse_peak_values[i] > 0 and (val / pulse_peak_values[i]) >= ACTIVE_THRESHOLD_FRAC:
                candidates.append(i)
        return (max(candidates, key=lambda i: pulse_params[i][1]) + 1) if candidates else None  # 1-indexed

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
    """One row per (pulse, episode) match, same schema as ../fitter.py's CSV plus a `fit_method` column
    and an `A_norm` column derived either way (A/y_max) so the two methods' CSVs are directly comparable
    column-for-column, not just individually self-consistent."""
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
    """Twin-axis LAT photon overlay (red hollow circles, gridlines aligned), same convention as
    ../fitter.py -- factored out since both driver scripts' plots need it identically."""
    ax_photon = ax.twinx()
    ax_photon.scatter(
        PHOTON_T_ARR_S, PHOTON_ENERGY_MEV,
        marker="o", facecolors="none", edgecolors="red", linewidths=1.2, s=MARKER_SIZE ** 2,
        label="LAT photons (E > 1 GeV)",
    )
    ax_photon.set_ylabel("Photon energy [MeV]")

    n_yticks = 6
    y_buffer_frac = 0.25
    photon_top = PHOTON_ENERGY_MEV.max()

    ax.set_ylim(-y_buffer_frac * y_top_data, y_top_data * (1 + y_buffer_frac))
    ax_photon.set_ylim(-y_buffer_frac * photon_top, photon_top * (1 + y_buffer_frac))
    ax.set_yticks(np.linspace(0, y_top_data, n_yticks))
    ax_photon.set_yticks(np.linspace(0, photon_top, n_yticks))

    # Explicit, matching minor-tick subdivision for both axes -- see fitter_GRB080916C.py's own comment
    # (same fix, applied there first): AutoMinorLocator's default heuristic picks a subdivision count per
    # axis from its own major-tick step size, so the two twinned axes' minor ticks/gridlines drift out of
    # alignment even when the major ticks line up. Forcing the same explicit n on both keeps them matched.
    minor_ticks_per_major = 5
    ax.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    ax_photon.yaxis.set_minor_locator(AutoMinorLocator(minor_ticks_per_major))
    return ax_photon
