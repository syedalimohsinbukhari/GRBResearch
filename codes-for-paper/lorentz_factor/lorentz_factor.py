"""Minimum Lorentz factor Calculator
===================================
Computes Gamma_min using the gamma-gamma opacity method.

Reference: Lithwick & Sari (2001), ApJ, 555, 540
           Abdo et al. (2009), Science, 323, 1688

Formula (Lithwick & Sari 2001, Limit A, Table 1 with redshift corrections):

    tau_hat = 2.1e11 * [(d_L/7Gpc)^2 * (0.511)^(-alpha+1) * f_1] / [(delta_T/0.1s) * (alpha-1)]

    Gamma_min = tau_hat^(1/(2a+2))
                * (E_max/0.511)^((a-1)/(2a+2))
                * (1+z)^((a-1)/(a+1))

    where:
        alpha    = Lithwick & Sari photon index (positive) = -beta (high-energy index)
        f_1      = photon flux at 1 MeV [ph/cm^2/s/MeV]
        E_max    = highest LAT photon energy [MeV]; 0.511 MeV = m_e c^2
        delta_T  = variability timescale [s]
        z        = redshift
        d_L      = luminosity distance [cm]

Spectral parameters are read from ``results.json`` through the ``grb_research`` class API
(``prepare_grbs`` -> ``GRB`` -> ``Model``) rather than being restated here, so this script cannot drift out of sync
with the fitted-model database.
Only quantities that do not live in ``results.json`` — the LAT photon properties and the redshifts — are tabulated below.

delta_T (t_v) precedence, per episode (see ``variability_timescale()``):
    1. A literature override, if one is entered in ``VARIABILITY_TIMESCALE``;
    2. A measured value from the Norris-pulse fits in ``codes-for-paper/variability_analysis/`` (Phase 5), if one exists
       for that episode and passes the ``MC_KEPT_FRACTION_MIN`` quality gate -- see ``load_norris_tv()`` and
       ``lorentz_factor.md``.  The pulse is the one the fitters themselves assign the episode's LAT photon to
       (nearest preceding onset among still-active pulses; ``lorentz_factor.md`` section 15);
    3. Otherwise the episode's own duration, an upper bound on the true t_v (conservative by construction, the original
       convention).

Outputs:
    - lorentz_results.csv   — all computed values
    - lorentz_table.tex     — LaTeX table for paper, rendered by generate_lorentz_table.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from grb_research import draw_model_samples, find_project_root, get_rng, prepare_grbs, seed_from_name
from grb_research.grb_constants import model_n_pars, N_SAMPLES
from grb_research.grb_enums import GRBModelsCombinations as gmC
from grb_research.grb_sed import MODEL_MAP
from grb_research.grb_seds import band_function, cutoff_powerlaw, powerlaw, smoothly_broken_power_law

# Photon-spectrum function for each single component, keyed as in SpectralModels.
SED_FUNCTIONS = {
    gmC.PL: powerlaw,
    gmC.CPL: cutoff_powerlaw,
    gmC.BAND: band_function,
    gmC.SBPL: smoothly_broken_power_law,
}

# ─── Cosmology ────────────────────────────────────────────────────────────────
# Fana Dirirsa et al. (2019) — the single cosmology used across this project.
# Must stay reconciled with tex_files/section-5-data-analysis.tex and
# grb_calculations.py::mc_e_iso_sampler; do not introduce a third value.
H0 = 69.6
OM0 = 0.286
cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)

# ─── Sample ───────────────────────────────────────────────────────────────────
# The four GRBs of this paper. Everything else is keyed off this list.
GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

TEX_NAMES = {
    "080916C": r"\grbzeroeightzeroninesixteenC",
    "131014A": r"\grbthirteentenfourteenA",
    "140206B": r"\grbfourteenzerotwozerosixB",
    "231129C": r"\grbtwentythreeeleventwentynineC",
}

E_1MEV_KEV = 1000.0

# Parameter holding the high-energy photon index, per continuum model.
# The Lithwick & Sari index is its negative.
HIGH_ENERGY_INDEX_PARAM = {gmC.BAND: "index2_band", gmC.SBPL: "index2_sbpl", gmC.PL: "index1_pl"}

# Per-episode LAT products, read from LAT_analysis/lat_photons.csv, which is built from the gtburst/gtlike output files by LAT_analysis/txt_to_csv.py.
# These values used to be transcribed by hand from the paper's appendix table, so an edit to the table left this script silently stale -- BUGS.md OBS-08.
#
# Photon energies are MeV (the source field is "P > 0.9 Max (E) MeV"; the LAT selection floor here is 100 MeV).
# See BUGS.md BUG-18.
LAT_PHOTONS_CSV = find_project_root() / "LAT_analysis" / "lat_photons.csv"

# LAT detection below this TS is not a secure association, so the derived Gamma_min is flagged rather than tabulated alongside the secure ones.
# This is the threshold quoted in the table caption and in the appendix footnotes.
TS_SECURE_DETECTION = 25.0

# Below this, the highest-energy photon gives too weak a pair-production constraint to report as a
# Limit A/B bound. GRB080916C only -- see lorentz_factor.md sec 12.
PHOTON_E_MIN_MEV = 1000.0

# Kept despite a sub-GeV photon: both are 2 of the 3 BB-inclusive episodes (with T90) the existing
# thermal-vs-opacity comparison (Figure 9, sec 8.6) depends on; dropping them would silently
# invalidate an already-published comparison. TR4/TR5/EX1 have no such stake and are cut cleanly.
PHOTON_ENERGY_CUT_EXEMPT_EPISODES = {"EX0", "TR1"}

# Measured variability timescale from the manual Norris-pulse fits in codes-for-paper/variability_analysis/
# (Phase 5, PHASE5_TV_PLAN.md / NORRIS_TV_LORENTZ_INTEGRATION_PLAN.md). Read as plain CSVs rather than an
# import across folders, per CLAUDE.md's "copy rather than fight sys.path" convention.
VARIABILITY_ANALYSIS_DIR = find_project_root() / "codes-for-paper" / "variability_analysis"

# A Norris-measured t_v is only trusted if its MC propagation kept at least this fraction of draws;
# variability_analysis.md flags several pulses below this as "unresolved" (e.g. GRB131014A pulses 1-2,
# 0.10-0.21). Below the threshold the episode falls back to duration, exactly as before -- see load_norris_tv().
MC_KEPT_FRACTION_MIN = 0.5


def load_lat_photons(csv_path=LAT_PHOTONS_CSV):
    """Read the LAT table, returning per-episode photon data and the weak detections.

    Returns
    -------
    photons          : {short_name: {episode: (E_max [MeV], arrival time [s])}}
    low_significance : {(short_name, episode)} for every episode with TS < 25
    """
    table = pd.read_csv(csv_path)
    photons: dict[str, dict[str, tuple[float, float]]] = {}
    low_significance: set[tuple[str, str]] = set()

    for row in table.itertuples(index=False):
        short_name = row.grb_name.removeprefix("GRB")
        photons.setdefault(short_name, {})[row.episode] = (row.e_max_MeV, row.t_arr_s)
        if row.ts < TS_SECURE_DETECTION:
            low_significance.add((short_name, row.episode))

    return photons, low_significance


LAT_PHOTONS, LOW_SIGNIFICANCE = load_lat_photons()


# Photon-to-pulse rule of the Norris fitters (variability_analysis/shared_utilities.py::assign_pulse), reimplemented here
# rather than imported across folders (this project's "copy rather than fight sys.path" convention).  A photon at time t
# is assigned to the pulse with the latest onset t_s <= t among the pulses still "active" at t, i.e. whose own value at t
# is at least ACTIVE_THRESHOLD_FRAC of its own peak.  ACTIVE_THRESHOLD_FRAC must match shared_utilities.py's constant;
# verified 2026-10-02 against all 15 photon assignments the fitters recorded in norris_fit_results_GRB*.csv.
ACTIVE_THRESHOLD_FRAC = 0.01


def _norris_value(t, amplitude, t_s, tau1, tau2):
    """Norris pulse value at ``t`` (zero before onset); same form as ``norris_fit.norris_pulse``."""
    if t <= t_s:
        return 0.0
    return amplitude * np.exp(2.0 * np.sqrt(tau1 / tau2)) * np.exp(-tau1 / (t - t_s) - (t - t_s) / tau2)


def assigned_pulse(pulses, t_arr):
    """Pulse index the fitters' rule assigns a photon at ``t_arr`` to, or ``None`` if no pulse is active then.

    ``pulses`` is a DataFrame with one row per distinct pulse (index = ``pulse_index``) and the columns
    ``A_cts_per_s``, ``t_s``, ``tau1``, ``tau2``, ``t_peak_s``.
    """
    active = []
    for index, pulse in pulses.iterrows():
        if pulse.t_s > t_arr:
            continue
        peak = _norris_value(pulse.t_peak_s, pulse.A_cts_per_s, pulse.t_s, pulse.tau1, pulse.tau2)
        value = _norris_value(t_arr, pulse.A_cts_per_s, pulse.t_s, pulse.tau1, pulse.tau2)
        if peak > 0 and value / peak >= ACTIVE_THRESHOLD_FRAC:
            active.append(index)
    return max(active, key=lambda j: pulses.loc[j].t_s) if active else None


def load_norris_tv(photons):
    """Per-episode measured t_v from the manual Norris-pulse fits, one row selected per episode.

    Several episodes have more than one candidate pulse (their window overlaps a neighbouring
    episode's -- e.g. EX0 is a strict superset of TR1's window, so a pulse belonging to TR1 also
    falls inside EX0).  The pulse that defines the episode's t_v is the one the Norris fitters themselves assign that
    episode's Gamma_min-defining LAT photon (``t_arr_s``, from `photons`) to: the pulse with the latest onset among
    those still active at ``t_arr_s`` (see ``assigned_pulse()``).  This replaced a nearest-``t_peak`` rule on
    2026-10-02 (user decision: the fitter is the primary selector); it changes five GRB231129C/GRB140206B episodes
    and none of GRB080916C's tabulated ones -- ``lorentz_factor.md`` section 15.  If the assigned pulse is not one of
    the episode's own candidate rows (or no pulse is active), the episode falls back to its duration with the reason
    recorded, never silently.

    A selected candidate is only used if its `mc_kept_fraction` is at least `MC_KEPT_FRACTION_MIN`;
    otherwise the rejection is recorded (`rejected_reason`) and the caller falls back to duration --
    per `CLAUDE.md`'s "flag rather than silently use" convention, this is never a silent drop.

    Parameters
    ----------
    photons :
        `LAT_PHOTONS`-shaped: `{short_name: {episode: (E_max_MeV, t_arr_s)}}`.

    Returns
    -------
    dict :
        `{(short_name, episode): {"t_v_s", "t_v_err_lower_s", "t_v_err_upper_s", "pulse_index",
        "mc_kept_fraction", "rejected_reason"}}`. `t_v_s` is `None` when the only candidate(s)
        failed the quality gate.
    """
    result = {}
    for short_name, episode_photons in photons.items():
        csv_path = VARIABILITY_ANALYSIS_DIR / f"norris_fit_results_GRB{short_name}.csv"
        if not csv_path.exists():
            continue

        table = pd.read_csv(csv_path)
        table = table[table["episode"].notna()]
        # One row per distinct pulse of the burst (a pulse repeats across every episode window it falls in); the fitters'
        # assignment rule looks at all of a burst's pulses, not just those of one episode.
        pulses = table.drop_duplicates("pulse_index").set_index("pulse_index").sort_index()

        for episode, t_arr in ((ep, t) for ep, (_, t) in episode_photons.items()):
            candidates = table[table["episode"] == episode]
            if candidates.empty:
                continue

            assigned = assigned_pulse(pulses, t_arr)
            row_match = candidates[candidates["pulse_index"] == assigned] if assigned is not None else candidates.iloc[:0]
            if row_match.empty:
                result[(short_name, episode)] = {
                    "t_v_s": None,
                    "t_v_err_lower_s": None,
                    "t_v_err_upper_s": None,
                    "pulse_index": None if assigned is None else int(assigned),
                    "mc_kept_fraction": None,
                    "rejected_reason": (
                        "no pulse is active at the photon's arrival" if assigned is None
                        else f"assigned pulse {int(assigned)} is not one of this episode's pulse rows"
                    ),
                }
                continue

            row = row_match.iloc[0]
            kept_fraction = float(row["mc_kept_fraction"])

            if kept_fraction < MC_KEPT_FRACTION_MIN:
                result[(short_name, episode)] = {
                    "t_v_s": None,
                    "t_v_err_lower_s": None,
                    "t_v_err_upper_s": None,
                    "pulse_index": int(row["pulse_index"]),
                    "mc_kept_fraction": kept_fraction,
                    "rejected_reason": f"mc_kept_fraction={kept_fraction:.4f} < {MC_KEPT_FRACTION_MIN}",
                }
            else:
                result[(short_name, episode)] = {
                    "t_v_s": float(row["t_v_s"]),
                    "t_v_err_lower_s": float(row["t_v_err_lower_s"]),
                    "t_v_err_upper_s": float(row["t_v_err_upper_s"]),
                    "pulse_index": int(row["pulse_index"]),
                    "mc_kept_fraction": kept_fraction,
                    "rejected_reason": None,
                }

    return result


NORRIS_TV = load_norris_tv(LAT_PHOTONS)

# Spectroscopic redshifts; None means Gamma_min cannot be computed.
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Variability timescale precedence: a literature override (this map) beats a measured Norris t_v
# (NORRIS_TV, Phase 5), which beats the episode's own duration -- the original, still-conservative
# fallback for any episode Phase 5 doesn't cover (T90, or a low-confidence Norris candidate).
# Populate this map to override an episode with a published value; mixing sources within one table
# is fine now -- build_latex_table() marks the source per cell (see NORRIS_TV_LORENTZ_INTEGRATION_PLAN.md).
VARIABILITY_TIMESCALE: dict = {}

# Monte-Carlo settings, matching the rest of the project.
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)
PERCENTILES = (16.0, 50.0, 84.0)


# ─── Spectral quantities from the fitted model ───────────────────────────────


def continuum_photon_flux(model_name, values, energy_kev):
    """Non-thermal continuum photon flux [ph/cm^2/s/keV] for one or many draws.

    Any blackbody component is dropped: the gamma-gamma opacity is set by the
    power-law photons that pair-produce with the LAT photon, not by the thermal
    component. Parameters are sliced in declaration order exactly as
    ``SpectralModels._evaluate_components`` does.

    ``values`` may be ``(n_pars,)`` or ``(n, n_pars)``; the result is ``(n,)``
    for a scalar energy.
    """
    values = np.atleast_2d(np.asarray(values, dtype=float))
    energy = np.atleast_1d(np.asarray(energy_kev, dtype=float))[None, :]

    key = gmC(model_name.lower())
    components = MODEL_MAP.get(key, (key,))

    total = np.zeros((values.shape[0], energy.shape[1]))
    idx = 0
    for component in components:
        n_pars = model_n_pars[component]
        pars = [values[:, idx + offset, None] for offset in range(n_pars)]
        idx += n_pars
        if component is gmC.BB:
            continue
        total = total + SED_FUNCTIONS[component](energy, *pars)

    return total[:, 0]


def f1_from_values(model_name, values):
    """Continuum photon flux at 1 MeV [ph/cm^2/s/MeV] for one or many draws."""
    return continuum_photon_flux(model_name, values, E_1MEV_KEV) * 1000.0


def high_energy_index_param_name(model):
    """Name of the parameter holding the continuum's high-energy index."""
    key = gmC(model.name.lower())
    continuum = key
    if key in MODEL_MAP:
        curved = [c for c in MODEL_MAP[key] if c in (gmC.BAND, gmC.SBPL)]
        continuum = curved[0] if curved else next(c for c in MODEL_MAP[key] if c is not gmC.BB)
    return HIGH_ENERGY_INDEX_PARAM.get(continuum)


def high_energy_index(model):
    """High-energy photon index beta of the continuum, read from the fitted model."""
    param = high_energy_index_param_name(model)
    return model.get_parameter_value(param) if param else None


# ─── Gamma_min formula ───────────────────────────────────────────────────────


def compute_tau_hat(alpha_LS, f_1, delta_T_s, z):
    """
    Dimensionless optical-depth quantity tau-hat (Lithwick & Sari 2001, eq. 4 / eq. 9).

    Shared by both Limit A (this module) and Limit B (``lorentz_factor_limit_b.py``) --
    the two limits differ only in how tau-hat is combined with E_max and z afterwards,
    not in tau-hat itself, so it is computed once here to avoid the formula drifting
    between the two files.

    Parameters
    ----------
    alpha_LS :
        L&S photon index (positive) = -beta (high-energy index)
    f_1 :
        photon flux at 1 MeV [ph/cm^2/s/MeV]
    delta_T_s :
        variability timescale [s]
    z :
        redshift

    Returns
    -------
    tau_hat : float
    """
    d_L_cm = cosmo.luminosity_distance(z).cgs.value
    d_7Gpc = d_L_cm / (7.0 * 3.0857e27)

    return 2.1e11 * d_7Gpc**2 * 0.511 ** (-alpha_LS + 1) * f_1 / ((delta_T_s / 0.1) * (alpha_LS - 1))


def compute_gamma_min(alpha_LS, f_1, E_max_MeV, delta_T_s, z):
    """
    Compute a minimum Lorentz factor (Lithwick & Sari 2001, Limit A).

    Parameters
    ----------
    alpha_LS :
        L&S photon index (positive) = -beta (high-energy index)
    f_1 :
        photon flux at 1 MeV [ph/cm^2/s/MeV]
    E_max_MeV :
        highest LAT photon energy [MeV]
    delta_T_s :
        variability timescale [s]
    z :
        redshift

    Returns
    -------
    gamma_min : float
    tau_hat   : float
    """
    tau_hat = compute_tau_hat(alpha_LS, f_1, delta_T_s, z)

    e1 = 1.0 / (2 * alpha_LS + 2)
    e2 = (alpha_LS - 1) / (2 * alpha_LS + 2)
    e3 = (alpha_LS - 1) / (alpha_LS + 1)

    # 0.511 MeV = m_e c^2, matching the MeV convention used for f_1 and tau_hat above.
    gamma_min = tau_hat**e1 * (E_max_MeV / 0.511) ** e2 * (1 + z) ** e3
    return gamma_min, tau_hat


# ─── COMPUTE ─────────────────────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def variability_timescale(short_name, episode, interval):
    """Variability timescale for an episode, its source, and its errors if measured.

    Precedence: a literature override (`VARIABILITY_TIMESCALE`) beats a measured Norris-fit t_v
    (`NORRIS_TV`, Phase 5), which beats the episode's own duration. The duration is an upper bound
    on the true variability timescale, so falling back to it keeps Gamma_min a conservative lower
    limit rather than an optimistic one; Gamma_min depends only weakly on t_v regardless of source,
    as delta_T^(-1/(2*alpha+2)).

    Returns
    -------
    t_v_s : float
    source : {"literature", "norris", "duration"}
    err_lower_s, err_upper_s : float or None
        Only set when source is "norris" -- a literature value and a duration are both treated as
        exact (no uncertainty to propagate); see `lorentz_factor.md` for the split-normal
        approximation used to turn these into per-draw resamples in `main()`.
    """
    key = (short_name, episode)
    if key in VARIABILITY_TIMESCALE:
        return VARIABILITY_TIMESCALE[key], "literature", None, None

    norris = NORRIS_TV.get(key)
    if norris is not None and norris["t_v_s"] is not None:
        return norris["t_v_s"], "norris", norris["t_v_err_lower_s"], norris["t_v_err_upper_s"]

    return interval.end - interval.start, "duration", None, None


def sample_split_normal(median, err_lower, err_upper, size, rng):
    """Draw `size` samples from a two-piece (split) normal built from an asymmetric 1-sigma interval.

    An approximation, not a reproduction of the Norris fit's own MC draws (which aren't persisted --
    only their 16/50/84 percentiles are). Used to propagate a measured t_v's own uncertainty into
    Gamma_min's MC loop; see `lorentz_factor.md` for why this approximation was chosen over
    re-running the Norris fit here. Draws are clipped to stay positive, since t_v enters
    `compute_gamma_min` as a divisor.
    """
    u = rng.standard_normal(size)
    sigma = np.where(u < 0.0, err_lower, err_upper)
    return np.clip(median + u * sigma, 1e-6, None)


def main():
    """Compute Gamma_min for every episode with LAT coverage, then tabulate."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    results = []
    header = f"{'GRB':<12}{'Ep.':<6}{'model':<10}{'alpha':>7}{'f_1':>12}{'dT [s]':>9}{'Gamma_min':>21}"
    print(f"\n{header}\n" + "-" * len(header))

    for short_name, grb in zip(GRB_LIST, grb_objects):
        redshift = REDSHIFTS[short_name]
        photons = LAT_PHOTONS[short_name]

        for model in grb.get_all_best_models():
            episode = episode_label(model.interval)
            if episode not in photons:
                continue

            e_max_mev, t_arr = photons[episode]
            if (
                short_name == "080916C"
                and e_max_mev < PHOTON_E_MIN_MEV
                and episode not in PHOTON_ENERGY_CUT_EXEMPT_EPISODES
            ):
                print(f"GRB{short_name:<9}{episode:<6}  skipped: E_max={e_max_mev:.1f} MeV < {PHOTON_E_MIN_MEV:.0f} MeV cut")
                continue

            delta_t, delta_t_source, delta_t_err_lo, delta_t_err_hi = variability_timescale(
                short_name, episode, model.interval
            )
            beta = high_energy_index(model)

            if redshift is None or beta is None or beta >= -1.0:
                # No redshift or no usable high-energy index: alpha_LS <= 1 makes the Lithwick & Sari expression singular.
                f_1 = alpha_ls = tau_hat = gamma = gamma_lo = gamma_hi = None
                print(f"GRB{short_name:<9}{episode:<6}{model.name:<10}{'—':>7}{'—':>12}{delta_t:>9.3f}{'—':>11}")
            else:
                f_1 = f1_from_values(model.name, [p.value for p in model.parameters])[0]
                alpha_ls = -beta
                gamma, tau_hat = compute_gamma_min(alpha_ls, f_1, e_max_mev, delta_t, redshift)

                # Statistical uncertainty: propagate the fit covariance through both f_1 and alpha, which are correlated because they come from the same spectral parameters.
                samples = draw_model_samples(model, n_samples=N_SAMPLES, rng=rng)
                index_position = [q.name for q in model.parameters].index(high_energy_index_param_name(model))
                alpha_draws = -samples[:, index_position]
                f1_draws = f1_from_values(model.name, samples)

                # When t_v is a measured Norris value, resample it too (split-normal approximation from
                # its own asymmetric error) so Gamma_min's error bars reflect its uncertainty as well,
                # not just f_1/alpha's -- a duration or literature t_v is treated as exact, as before.
                if delta_t_source == "norris":
                    delta_t_for_draws = sample_split_normal(delta_t, delta_t_err_lo, delta_t_err_hi, N_SAMPLES, rng)
                else:
                    delta_t_for_draws = delta_t

                usable = np.isfinite(f1_draws) & (f1_draws > 0) & (alpha_draws > 1.0)
                delta_t_usable = delta_t_for_draws[usable] if delta_t_source == "norris" else delta_t
                gamma_draws, _ = compute_gamma_min(alpha_draws[usable], f1_draws[usable], e_max_mev, delta_t_usable, redshift)
                gamma_draws = gamma_draws[np.isfinite(gamma_draws)]
                lo, med, hi = np.percentile(gamma_draws, PERCENTILES)
                gamma_lo, gamma_hi = med - lo, hi - med

                flag = "*" if (short_name, episode) in LOW_SIGNIFICANCE else ""
                print(
                    f"GRB{short_name:<9}{episode:<6}{model.name:<10}{alpha_ls:>7.3f}"
                    f"{f_1:>12.4e}{delta_t:>9.3f}{gamma:>8.0f} -{gamma_lo:<5.1f}+{gamma_hi:<5.1f}{flag}"
                )

            results.append(
                {
                    "GRB": f"GRB{short_name}",
                    "tex_name": TEX_NAMES[short_name],
                    "episode": episode,
                    "H0": H0,
                    "Om0": OM0,
                    "z": redshift,
                    "E_max_MeV": e_max_mev,
                    "t_arr_s": t_arr,
                    "t_v_s": delta_t,
                    "t_v_source": delta_t_source,
                    "t_v_err_lower_s": delta_t_err_lo,
                    "t_v_err_upper_s": delta_t_err_hi,
                    "model": model.name,
                    "beta": beta,
                    "alpha_LS": alpha_ls,
                    "f_1": f_1,
                    "tau_hat": tau_hat,
                    "Gamma_min": gamma,
                    "Gamma_min_err_lower": gamma_lo,
                    "Gamma_min_err_upper": gamma_hi,
                    "n_samples": N_SAMPLES,
                    "seed": SEED,
                    "low_significance": (short_name, episode) in LOW_SIGNIFICANCE,
                }
            )

    df = pd.DataFrame(results)
    df.drop(columns=["tex_name"]).to_csv("lorentz_results.csv", index=False)
    print("\nSaved: lorentz_results.csv")

    subprocess.run([sys.executable, str(Path(__file__).parent / "generate_lorentz_table.py")], check=True)


if __name__ == "__main__":
    main()
