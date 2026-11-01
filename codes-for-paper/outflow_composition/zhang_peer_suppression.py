"""
Minimum Poynting-flux fraction that hides the photosphere — Zhang & Pe'er (2009), for the episodes without a blackbody
====================================================================================================================
For an episode whose best-fit model has no thermal component, a pure baryonic fireball with wind luminosity
L_w >= L_gamma would have released a bright photospheric blackbody that the data do not show. Zhang & Pe'er turn that
non-detection into a lower limit on the magnetisation sigma = L_P / L_b at the photosphere.

Reference
---------
Zhang & Pe'er (2009), ApJ 700, L65 (arXiv:0904.2943), section 3:

    L_th = L_w                                      eta >  eta_*  (R_ph < R_c)                           (3)
    kT_ob = (L_w / 4 pi R_0^2 c a)^(1/4) / (1 + z)  = kT_max                                              (text)
    eta_* = (L_w sigma_T / 8 pi m_p c^3 R_0)^(1/4),   R_0 = c * dt                                       (text)

With a Poynting flux, "the above photosphere derivations can be modified by replacing L_w by L_w / (1 + sigma)":

    L_th = L_w / (1 + sigma),     kT_ob = kT_max * (1 + sigma)^(-1/4)

and sigma_min is the smallest sigma for which this blackbody lies below the observed spectrum. This script uses the
eta >= eta_* branch (the paper's thick curve), L_w = L_gamma (100% efficiency, as the paper does), and
R_0 = c * t_v / (1 + z) with t_v the observer-frame variability timescale.

Decisions (user, 2026-10-24: "go with those defaults")
-----------------------------------------------------
    * eta >= eta_* branch; the paper's eta < eta_* curves are not computed (they need an extra parameter).
    * R_0 = c t_v / (1 + z). ``validate_r0_convention()`` shows this is what the paper's own numbers imply: it gives
      eta_* = 812 for the quoted kT_max = 50 keV (quoted 850), while R_0 = c t_v gives 1235.
    * L_gamma is integrated over 10 keV - 10 GeV in the observer frame, the band of Abdo et al. (2009)'s fluence that
      Zhang & Pe'er's L_gamma comes from. ``validate_epoch_b()`` shows this is what reproduces their numbers: with
      Abdo's own epoch-b spectrum it gives kT_max = 53 keV (quoted 50), eta_* = 866 (850) and sigma = 22 (about 20);
      the project's rest-frame 1 - 10^4 keV band (``rest_1_1e4``, the Amati E_iso band) gives L_gamma 2.5x lower and
      sigma = 7.7. Both are tabulated, the Zhang & Pe'er band as the primary.
    * "Below the observed spectrum" is tested over a band: the primary band is 10 keV - 10 MeV (the GBM fit band); the
      sensitivity band extends to 1 keV (the Band extrapolation drawn in the paper's Fig. 2).

Why this is a root-find and not a closed form: the blackbody's peak flux falls as (1+sigma)^-1 and its peak energy as
(1+sigma)^(-1/4), so the largest ratio blackbody / observed spectrum over the band has no closed-form inverse.
The ratio is strictly decreasing in sigma at every energy (``check_monotonic()``), so the root is unique and a
bisection is exact. A draw whose ratio is already <= 1 at sigma = 0 needs no suppression: sigma_min = 0, and the
fraction of draws needing suppression is reported. A draw that cannot be bracketed raises -- it is never defaulted.

Outputs:
    zhang_peer_suppression.csv -- one row per non-blackbody episode per redshift per comparison band
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from grb_research import (
    component_energy_fluxes,
    draw_model_samples,
    find_project_root,
    get_rng,
    prepare_grbs,
    seed_from_name,
)
from grb_research.grb_constants import model_n_pars, N_SAMPLES, kev_to_erg
from grb_research.grb_enums import GRBModelsCombinations as gmC
from grb_research.grb_sed import MODEL_MAP
from grb_research.grb_seds import band_function, cutoff_powerlaw, powerlaw, smoothly_broken_power_law

# ─── Physical constants (CGS) ────────────────────────────────────────────────

SIGMA_T = 6.6524587321e-25  # Thomson cross section [cm^2]
M_P = 1.67262192369e-24  # proton mass [g]
C_LIGHT = 2.99792458e10  # speed of light [cm/s]
A_RAD = 7.565733250e-15  # radiation constant [erg cm^-3 K^-4]
K_B_KEV = 8.617333262e-8  # Boltzmann constant [keV K^-1]

# ─── Cosmology ───────────────────────────────────────────────────────────────
# Fana Dirirsa et al. (2019) — the single cosmology used across this project.
H0 = 69.6
OM0 = 0.286
COSMO = FlatLambdaCDM(H0=H0, Om0=OM0)

# ─── Configuration ───────────────────────────────────────────────────────────

GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

# Spectroscopic redshifts; None means unmeasured and therefore tabulated at Z_VALUES.
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Assumed-redshift grid for the three bursts without a measured z (CLAUDE.md "Output conventions").
Z_VALUES = (1.0, 2.0, 3.0, 5.0, 7.0)

# Photon-spectrum function for each single component, keyed as in SpectralModels.
SED_FUNCTIONS = {
    gmC.PL: powerlaw,
    gmC.CPL: cutoff_powerlaw,
    gmC.BAND: band_function,
    gmC.SBPL: smoothly_broken_power_law,
}

# Band over which L_gamma is integrated, as (frame, low, high) in keV. "zp09" is Zhang & Pe'er's (Abdo et al. 2009's
# fluence band, observer frame); "rest_1_1e4" is the project's E_iso band in the Amati analysis (observer-frame limits
# divided by 1 + z).
# The comparison of the blackbody with the observed spectrum is always in the observer frame.
L_BANDS = {"zp09": ("observer", 10.0, 1.0e7), "rest_1_1e4": ("rest", 1.0, 1.0e4)}
N_L_GRID = 1000

# Observer-frame comparison bands [keV]: "fit" is the GBM fit band, "extrapolated" adds the Band extrapolation to 1 keV.
COMPARISON_BANDS = {"fit": (10.0, 1.0e4), "extrapolated": (1.0, 1.0e4)}
N_COMPARE = 200  # log-spaced energies per comparison band

# Bisection on u = ln(1 + sigma): bracket [0, U_MAX] (1 + sigma up to e^U_MAX), then N_BISECT halvings (final width
# U_MAX / 2^N_BISECT ~ 3e-7 in ln(1 + sigma)).
U_MAX = float(np.log(1.0e8))
N_BISECT = 26

# Wien limit beyond which x^4 / (e^x - 1) is zero in double precision.
X_MAX = 700.0

PERCENTILES = (16.0, 50.0, 84.0)
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)

# Reproduction target: Zhang & Pe'er's epoch b of GRB080916C (dt_ob = 0.5 s, kT_max = 50 keV, eta_* = 850, sigma ~ 20)
# with Abdo et al. (2009) Table 1's Band parameters for 3.58 - 7.68 s (A in ph cm^-2 s^-1 keV^-1 at 100 keV, E_peak in
# keV) and
# its 50 - 300 keV photon flux, which checks the transcription of the table.
EPOCH_B_DT_OB_S = 0.5
EPOCH_B_KT_MAX_KEV = 50.0
EPOCH_B_ETA_STAR = 850.0
EPOCH_B_SIGMA = 20.0
EPOCH_B_Z = 4.35
ABDO_EPOCH_B = {"amp": 35e-3, "e_peak": 1170.0, "alpha": -1.02, "beta": -2.21, "photon_flux_50_300": 5.63}

# ─── Zhang & Pe'er (2009) relations ──────────────────────────────────────────


def kt_max_kev(l_w, r0_cm, z):
    """Observed temperature of the eta > eta_* fireball photosphere [keV]."""
    return K_B_KEV * (l_w / (4 * np.pi * r0_cm**2 * C_LIGHT * A_RAD)) ** 0.25 / (1 + z)


def eta_star(l_w, r0_cm):
    """Critical dimensionless entropy eta_*."""
    return (l_w * SIGMA_T / (8 * np.pi * M_P * C_LIGHT**3 * r0_cm)) ** 0.25


def base_radius_cm(t_v_s, z):
    """R_0 = c t_v / (1 + z), t_v in the observer frame."""
    return C_LIGHT * np.asarray(t_v_s) / (1 + z)


def blackbody_nufnu(energy_kev, kt_kev, flux):
    """nu F_nu [same unit as ``flux``] of a blackbody of temperature ``kt_kev`` and bolometric flux ``flux``.

    nu F_nu = F * (15 / pi^4) * x^4 / (e^x - 1),  x = E / kT. Arrays broadcast; the Wien tail beyond ``X_MAX`` is zero.
    """
    x = np.asarray(energy_kev) / np.asarray(kt_kev)
    x = np.minimum(x, X_MAX)
    return np.asarray(flux) * (15.0 / np.pi**4) * x**4 / np.expm1(x)


def max_ratio(u, kt0_kev, flux0, obs_nufnu, energy_kev):
    """max over energy of (suppressed blackbody) / (observed nu F_nu), for ln(1+sigma) = ``u``; shapes (n,)."""
    shrink = np.exp(-u)
    kt = kt0_kev * shrink**0.25
    ratio = blackbody_nufnu(energy_kev[None, :], kt[:, None], (flux0 * shrink)[:, None]) / obs_nufnu
    return ratio.max(axis=1), ratio.argmax(axis=1)


def solve_one_plus_sigma(kt0_kev, flux0, obs_nufnu, energy_kev):
    """Smallest 1 + sigma for which the blackbody lies at or below the observed spectrum, per draw.

    Returns (one_plus_sigma, binding_index). ``one_plus_sigma`` is exactly 1 where no suppression is needed. Raises if
    any draw's ratio is still above 1 at 1 + sigma = e^U_MAX.
    """
    n = kt0_kev.size
    start, _ = max_ratio(np.zeros(n), kt0_kev, flux0, obs_nufnu, energy_kev)
    needs = start > 1.0

    top, _ = max_ratio(np.full(n, U_MAX), kt0_kev, flux0, obs_nufnu, energy_kev)
    if np.any(needs & (top > 1.0)):
        raise RuntimeError(f"{int(np.sum(needs & (top > 1.0)))} draws not bracketed by 1 + sigma = {np.exp(U_MAX):.1e}")

    low, high = np.zeros(n), np.full(n, U_MAX)
    for _ in range(N_BISECT):
        mid = 0.5 * (low + high)
        above, _ = max_ratio(mid, kt0_kev, flux0, obs_nufnu, energy_kev)
        low = np.where(above > 1.0, mid, low)
        high = np.where(above > 1.0, high, mid)

    u = np.where(needs, high, 0.0)
    _, binding = max_ratio(u, kt0_kev, flux0, obs_nufnu, energy_kev)
    return np.exp(u), binding


# ─── Validation ──────────────────────────────────────────────────────────────


def validate_blackbody_normalisation(verbose=True):
    """The blackbody nu F_nu integrates to the bolometric flux it was given (|error| < 1e-3)."""
    energy = np.logspace(-3, 4, 20_000)
    worst = 0.0
    for kt in (1.0, 10.0, 50.0):
        integral = np.trapezoid(blackbody_nufnu(energy, kt, 1.0), np.log(energy))
        worst = max(worst, abs(integral - 1.0))
    if verbose:
        print(f"Blackbody normalisation: max |integral - 1| = {worst:.2e} over kT = 1, 10, 50 keV")
    assert worst < 1e-3, f"blackbody nu F_nu does not integrate to its flux ({worst:.2e})"
    return worst


def check_monotonic(verbose=True):
    """The max ratio is strictly decreasing in ln(1+sigma), over a grid of (kT, F, spectrum) -- so the root is unique.

    The analytic derivative of ln[(1+sigma)^-1 x^4 / (e^x - 1)] with x ~ (1+sigma)^(1/4) is
    -(x/4) e^x / (e^x - 1) / (1 + sigma) < 0 at every energy; this confirms it numerically for several spectra.
    """
    energy = np.logspace(0, 4, N_COMPARE)
    u = np.linspace(0.0, 12.0, 400)
    worst = -np.inf
    for kt0 in (3.0, 25.0, 200.0):
        for alpha, beta, e_peak in ((-0.5, -2.5, 300.0), (-1.0, -3.0, 80.0), (-1.5, -2.2, 1000.0)):
            obs = energy**2 * band_function(energy, 1.0, e_peak, alpha, beta)
            ratios = np.array(
                [max_ratio(np.array([uu]), np.array([kt0]), np.array([1.0e-6]), obs[None, :], energy)[0][0] for uu in u]
            )
            worst = max(worst, np.max(np.diff(ratios) / ratios[:-1]))
    if verbose:
        print(f"Monotonicity: largest relative step of the max ratio along ln(1+sigma) = {worst:.2e} (must be < 0)")
    assert worst < 0, "max ratio is not strictly decreasing in sigma"
    return worst


def check_solver(verbose=True):
    """At the returned root the max ratio is 1; with a spectrum far above the blackbody the answer is exactly 1."""
    energy = np.logspace(0, 4, N_COMPARE)
    obs = energy**2 * band_function(energy, 1.0, 300.0, -0.7, -2.4)
    kt0 = np.array([5.0, 30.0, 120.0])
    flux0 = np.array([2.0, 5.0, 40.0]) * obs.max() * 10.0
    root, _ = solve_one_plus_sigma(kt0, flux0, np.tile(obs, (3, 1)), energy)
    at_root, _ = max_ratio(np.log(root), kt0, flux0, np.tile(obs, (3, 1)), energy)
    quiet, _ = solve_one_plus_sigma(kt0, flux0 * 1e-9, np.tile(obs, (3, 1)), energy)
    if verbose:
        print(f"Solver: ratio at the root = {np.array2string(at_root, precision=6)} (1 expected); "
              f"1+sigma = {np.array2string(root, precision=1)}; suppression-free case gives {quiet}")
    assert np.allclose(at_root, 1.0, atol=1e-6), "solver root does not sit at ratio = 1"
    assert np.all(quiet == 1.0), "a spectrum above the blackbody must give 1 + sigma = 1 exactly"


def validate_r0_convention(verbose=True):
    """Which R_0 convention reproduces the paper's pair (kT_max = 50 keV, eta_* = 850) for delta t_ob = 0.5 s, z = 4.35.

    Inverts kT_max for L_w under each convention and reports the eta_* that follows. The (1+z) convention is within
    10% of the quoted value; dropping it is off by ~45%.
    """
    z = 4.35
    out = {}
    for label, r0 in (("c dt/(1+z)", base_radius_cm(EPOCH_B_DT_OB_S, z)), ("c dt", C_LIGHT * EPOCH_B_DT_OB_S)):
        l_w = 4 * np.pi * r0**2 * C_LIGHT * A_RAD * (EPOCH_B_KT_MAX_KEV / K_B_KEV * (1 + z)) ** 4
        out[label] = (l_w, eta_star(l_w, r0))
    if verbose:
        print(f"R_0 convention (quoted: kT_max = {EPOCH_B_KT_MAX_KEV:.0f} keV, eta_* = {EPOCH_B_ETA_STAR:.0f}):")
        for label, (l_w, eta) in out.items():
            print(f"  R_0 = {label:<11} -> L_w = {l_w:.2e} erg/s, eta_* = {eta:.0f}")
    assert abs(out["c dt/(1+z)"][1] / EPOCH_B_ETA_STAR - 1) < 0.10, "(1+z) convention fails to reproduce eta_*"
    assert abs(out["c dt"][1] / EPOCH_B_ETA_STAR - 1) > 0.25, "convention test no longer discriminates"
    return out


# ─── Per-episode computation ─────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def observed_nufnu(model_name, values, energy_kev):
    """Observed E^2 N(E) [erg cm^-2 s^-1] of the non-thermal components, shape (n, n_energy).

    Any blackbody component is dropped (these episodes have none; kept so the function is safe for any model).
    Parameters are sliced in declaration order, as ``SpectralModels._evaluate_components`` does.
    """
    values = np.atleast_2d(np.asarray(values, dtype=float))
    energy = np.asarray(energy_kev, dtype=float)[None, :]
    key = gmC(model_name.lower())
    total = np.zeros((values.shape[0], energy.shape[1]))
    idx = 0
    for component in MODEL_MAP.get(key, (key,)):
        n_pars = model_n_pars[component]
        pars = [values[:, idx + offset, None] for offset in range(n_pars)]
        idx += n_pars
        if component is not gmC.BB:
            total += SED_FUNCTIONS[component](energy, *pars)
    return total * energy**2 * kev_to_erg


def l_band_energies(l_band, z):
    """Observer-frame energy grid [keV] that realises ``L_BANDS[l_band]`` at redshift ``z``."""
    frame, low, high = L_BANDS[l_band]
    scale = 1.0 / (1 + z) if frame == "rest" else 1.0
    return np.logspace(np.log10(low * scale), np.log10(high * scale), N_L_GRID)


def isotropic_luminosity(model_name, samples, z, l_band):
    """L_gamma [erg/s] over ``L_BANDS[l_band]``: 4 pi d_L^2 * F[band], one value per draw."""
    d_l = COSMO.luminosity_distance(z).cgs.value
    _, total = component_energy_fluxes(model_name, samples, l_band_energies(l_band, z))
    return 4 * np.pi * d_l**2 * total * kev_to_erg


def bolometric_flux(l_w, z):
    """Observed bolometric flux of a luminosity L: F = L / (4 pi d_L^2)."""
    return l_w / (4 * np.pi * COSMO.luminosity_distance(z).cgs.value ** 2)


def sample_split_normal(median, err_lower, err_upper, size, rng_):
    """Two-piece normal draws from an asymmetric 1-sigma interval, clipped positive (copied from lorentz_factor.py)."""
    u = rng_.standard_normal(size)
    sigma = np.where(u < 0.0, err_lower, err_upper)
    return np.clip(median + u * sigma, 1e-6, None)


def load_variability_timescales():
    """{(grb_name, episode): (t_v [s], source, err_lower, err_upper)} from lorentz_results.csv, the one t_v record."""
    frame = pd.read_csv(find_project_root() / "codes-for-paper" / "lorentz_factor" / "lorentz_results.csv")
    return {
        (r.GRB, r.episode): (r.t_v_s, r.t_v_source, r.t_v_err_lower_s, r.t_v_err_upper_s) for r in frame.itertuples()
    }


def summarise(values):
    """(median, minus, plus) from the 16/50/84 percentiles."""
    lo, med, hi = np.percentile(values, PERCENTILES)
    return med, med - lo, hi - med


def suppression_bound(model, samples, z, t_v_draws, band, l_band):
    """Draw-wise (1 + sigma_min, L_gamma, kT_max, eta_*, binding energy) of one episode at one z and band."""
    l_gamma = isotropic_luminosity(model.name, samples, z, l_band)
    r0 = base_radius_cm(t_v_draws, z)
    kt0 = kt_max_kev(l_gamma, r0, z)
    flux0 = bolometric_flux(l_gamma, z)

    energy = np.logspace(*np.log10(band), N_COMPARE)
    obs = observed_nufnu(model.name, samples, energy)

    usable = np.isfinite(l_gamma) & (l_gamma > 0) & np.all(np.isfinite(obs) & (obs > 0), axis=1)
    one_plus, binding = solve_one_plus_sigma(kt0[usable], flux0[usable], obs[usable], energy)
    return one_plus, l_gamma[usable], kt0[usable], eta_star(l_gamma[usable], r0[usable]), energy[binding]


def collect_results():
    """Suppression bound for every BEST model without a blackbody, every assumed z, both comparison bands."""
    root = find_project_root()
    timescales = load_variability_timescales()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    rows = []
    for short_name, grb in zip(GRB_LIST, grb_objects):
        grb_name = f"GRB{short_name}"
        measured_z = REDSHIFTS[short_name]
        z_grid, z_source = ((measured_z,), "spectroscopic") if measured_z is not None else (Z_VALUES, "assumed")

        for model in grb.get_all_best_models():
            if "BB" in model.name:
                continue
            episode = episode_label(model.interval)
            key = (grb_name, episode)
            if key in timescales:
                t_v, t_v_source, err_lo, err_hi = timescales[key]
            else:
                # No LAT photon, hence no Norris t_v: the duration, an upper bound on the true timescale.
                t_v, t_v_source, err_lo, err_hi = model.interval.end - model.interval.start, "duration", np.nan, np.nan

            samples = draw_model_samples(model, n_samples=N_SAMPLES, rng=rng)
            if t_v_source == "norris":
                t_v_draws = sample_split_normal(t_v, err_lo, err_hi, samples.shape[0], rng)
            else:
                t_v_draws = np.full(samples.shape[0], t_v)

            for z in z_grid:
                for l_band, band_name in itertools.product(L_BANDS, COMPARISON_BANDS):
                    band = COMPARISON_BANDS[band_name]
                    one_plus, l_gamma, kt0, eta_s, binding = suppression_bound(
                        model, samples, z, t_v_draws, band, l_band
                    )
                    needs = one_plus > 1.0
                    rows.append(
                        {
                            "grb_name": grb_name,
                            "episode": episode,
                            "model_name": model.name,
                            "z": float(z),
                            "z_source": z_source,
                            "H0": H0,
                            "Om0": OM0,
                            "branch": "eta>=eta_star",
                            "r0_convention": "c*t_v/(1+z)",
                            "l_w_over_l_gamma": 1.0,
                            "l_band": l_band,
                            "l_band_frame": L_BANDS[l_band][0],
                            "l_band_lo_keV": L_BANDS[l_band][1],
                            "l_band_hi_keV": L_BANDS[l_band][2],
                            "band": band_name,
                            "band_lo_keV": band[0],
                            "band_hi_keV": band[1],
                            "t_v_s": t_v,
                            "t_v_source": t_v_source,
                            "t_v_is_duration_fallback": t_v_source == "duration",
                            "needs_suppression_fraction": float(np.mean(needs)),
                            "binding_energy_keV": float(np.median(binding[needs])) if needs.any() else np.nan,
                            "binding_at_lower_edge_fraction": (
                                float(np.mean(np.isclose(binding[needs], band[0]))) if needs.any() else np.nan
                            ),
                            "summary": {
                                "one_plus_sigma_min": summarise(one_plus),
                                "l_gamma_erg_s": summarise(l_gamma),
                                "kt_max_keV": summarise(kt0),
                                "eta_star": summarise(eta_s),
                            },
                        }
                    )
            print(
                f"  {grb_name:<11}{episode:<5}{model.name:<10}t_v={t_v:8.3f} s ({t_v_source:<8}) "
                f"1+sigma_min(zp09, fit band, last z) = {rows[-4]['summary']['one_plus_sigma_min'][0]:9.2f}"
            )
    return rows


def validate_epoch_b(verbose=True):
    """Reproduce Zhang & Pe'er's epoch b of GRB080916C from Abdo et al. (2009) Table 1's own spectrum.

    Checks, in order: the table's Band parameters reproduce its quoted 50 - 300 keV photon flux (transcription); the
    resulting L_gamma with R_0 = c dt/(1+z) gives kT_max and eta_* near the paper's 50 keV and 850; and the solver's
    sigma lands near their "sigma ~ 20" (their text also says 15 - 20 for the other two cases).
    """
    amp, e_peak, alpha, beta = (ABDO_EPOCH_B[k] for k in ("amp", "e_peak", "alpha", "beta"))
    z = EPOCH_B_Z
    values = np.array([[amp, e_peak, alpha, beta]])

    energy = np.logspace(np.log10(50.0), np.log10(300.0), 4000)
    photon_flux = np.trapezoid(band_function(energy, amp, e_peak, alpha, beta), energy)
    flux_ratio = photon_flux / ABDO_EPOCH_B["photon_flux_50_300"]

    r0 = base_radius_cm(EPOCH_B_DT_OB_S, z)
    out = {}
    for l_band in L_BANDS:
        l_gamma = float(isotropic_luminosity("band", values, z, l_band)[0])
        kt0, eta_s = kt_max_kev(l_gamma, r0, z), eta_star(l_gamma, r0)
        sigma = {}
        for band_name, band in COMPARISON_BANDS.items():
            grid = np.logspace(*np.log10(band), N_COMPARE)
            obs = observed_nufnu("band", values, grid)
            one_plus, _ = solve_one_plus_sigma(np.array([kt0]), np.array([bolometric_flux(l_gamma, z)]), obs, grid)
            sigma[band_name] = float(one_plus[0]) - 1.0
        out[l_band] = (l_gamma, kt0, eta_s, sigma)

    if verbose:
        print(f"Epoch b (Abdo et al. 2009 Table 1: 3.58 - 7.68 s, z = {z}, dt_ob = {EPOCH_B_DT_OB_S} s):")
        print(f"  Table 1 Band parameters give a 50 - 300 keV photon flux of {photon_flux:.3f} ph/cm^2/s "
              f"(quoted {ABDO_EPOCH_B['photon_flux_50_300']}), ratio {flux_ratio:.3f}")
        print(f"  paper: kT_max = {EPOCH_B_KT_MAX_KEV:.0f} keV, eta_* = {EPOCH_B_ETA_STAR:.0f}, "
              f"sigma ~ {EPOCH_B_SIGMA:.0f}")
        for l_band, (l_gamma, kt0, eta_s, sigma) in out.items():
            print(f"  L_gamma band {l_band:<11} L = {l_gamma:.2e} erg/s, kT_max = {kt0:5.1f} keV, "
                  f"eta_* = {eta_s:4.0f}, "
                  f"sigma_min = {sigma['fit']:5.1f} (fit band), {sigma['extrapolated']:5.1f} (extrapolated)")

    l_gamma, kt0, eta_s, sigma = out["zp09"]
    assert abs(flux_ratio - 1) < 0.03, "Abdo Table 1 parameters do not reproduce the quoted photon flux"
    assert abs(kt0 / EPOCH_B_KT_MAX_KEV - 1) < 0.10, "kT_max does not reproduce the paper's 50 keV"
    assert abs(eta_s / EPOCH_B_ETA_STAR - 1) < 0.10, "eta_* does not reproduce the paper's 850"
    assert 0.75 * EPOCH_B_SIGMA < sigma["fit"] < 1.25 * EPOCH_B_SIGMA, "sigma does not reproduce the paper's ~20"
    return out


# ─── Outputs ─────────────────────────────────────────────────────────────────

QUANTITIES = (
    ("one_plus_sigma_min", "one_plus_sigma_min", ""),
    ("l_gamma_erg_s", "L_gamma_erg_s", ""),
    ("kt_max_keV", "kT_max_keV", ""),
    ("eta_star", "eta_star", ""),
)


def write_csv(rows, path="zhang_peer_suppression.csv"):
    """One self-describing row per episode per z per band."""
    records = []
    for row in rows:
        record = {k: v for k, v in row.items() if k != "summary"}
        for key, column, unit in QUANTITIES:
            median, lower, upper = row["summary"][key]
            record[f"{column}{unit}"] = median
            record[f"{column}_err_lower{unit}"] = lower
            record[f"{column}_err_upper{unit}"] = upper
        record["n_samples"] = N_SAMPLES
        record["seed"] = SEED
        records.append(record)
    frame = pd.DataFrame(records)
    frame.to_csv(path, index=False)
    print(f"Saved: {path}  ({len(frame)} rows)")
    return frame


def main():
    """Validate, then bound sigma for every non-blackbody episode."""
    validate_blackbody_normalisation()
    check_monotonic()
    check_solver()
    validate_r0_convention()
    validate_epoch_b()
    print(f"\nZhang & Pe'er (2009) suppression bound ({N_SAMPLES} MC samples)\n")
    rows = collect_results()
    frame = write_csv(rows)
    fallback = frame[frame.t_v_is_duration_fallback].drop_duplicates(["grb_name", "episode"])
    print(f"\n{len(frame)} rows; t_v duration fallback: {', '.join(fallback.grb_name + ' ' + fallback.episode)}")


if __name__ == "__main__":
    main()
