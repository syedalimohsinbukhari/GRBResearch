"""
Hybrid-outflow inversion for the BB-detected episodes — Gao & Zhang (2015), "top-down" approach
===============================================================================================
Turns an episode's observed blackbody temperature ``kT_ob`` and flux ``F_BB`` (plus the total flux ``F_ob``) into the
central-engine pair (eta, 1 + sigma_0), the photospheric radius and the photospheric Lorentz factor, for an *assumed*
jet-base radius r_0.

Reference
---------
Gao & Zhang (2015), ApJ 801, 103 (arXiv:1409.3584), section 4:

    * Table 2  -- observational criteria that select the photospheric regime (II, III, V or VI),
    * eqs. (38), (39), (40) -- the inversion in regimes II, III/VI and V.

Regimes I and IV (r_ph < r_ra) are degenerate and cannot be inverted (their own statement, as in Pe'er+07 for
r_ph <= r_c). Their top-down formulae exist only for the slow-acceleration index delta = 1/3 (their section 6).

Inputs, as in the paper:  kT_ob [keV], F_BB [erg s^-1 cm^-2], f_th = F_BB / F_ob, z, d_L and two constants taken
as known, f_gamma = L_gamma / L_w (0.5, the paper's and Li+24's fiducial) and r_0 (assumed, swept: see
``R0_GZ2015_CM`` and ``R0_ANCHOR_MULTIPLIERS``).

Every row of ``CRITERIA`` and ``INVERSIONS`` below is a literal copy of one printed row of Table 2 or of eqs. (38)-(40),
written as ``coefficient * (1+z)^a * (kT/kT0)^b * (F_BB/F0)^c * r_0,9^d * (f_th/0.1)^e * f_gamma^f * d_L,28^g``, so the
transcription can be audited line against the PDF.

Validation
----------
``validate_round_trip()`` generates (kT_ob, F_BB) from the paper's *bottom-up* equations (24)-(30) with Table 1's
regime criteria, then recovers (eta, 1 + sigma_0, r_ph, Gamma_ph) with the top-down code. The two sets of equations
are printed separately in the paper, so agreement tests the transcription of both, and Table 2 against Table 1.
Run it before trusting any output here.

Outputs:
    gao_zhang_inversion.csv         -- one row per BEST-BB episode per redshift per r_0 trial (the main Phase 6 table)
    gao_zhang_inversion_partner.csv -- the same for the BB *partner* of episodes whose BEST model has no blackbody
                                       (the BEST continuum X plus a blackbody, X_BB; SAFE or MARGINAL only), reported
                                       separately because the data do not require the blackbody (user, 2026-10-24)
"""

from __future__ import annotations

from dataclasses import dataclass

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
from grb_research.grb_constants import N_GRID, N_SAMPLES, kev_to_erg
from grb_research.grb_enums import GoodnessOfFit, GRBModelsCombinations as gmC

# ─── Cosmology ───────────────────────────────────────────────────────────────
# Fana Dirirsa et al. (2019) — the single cosmology used across this project.
H0 = 69.6
OM0 = 0.286
COSMO = FlatLambdaCDM(H0=H0, Om0=OM0)

# ─── Configuration ───────────────────────────────────────────────────────────

GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

# Spectroscopic redshifts; None means unmeasured and therefore tabulated at Z_VALUES.
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Assumed-redshift grid for the three bursts without a measured z, unified across the project on 2026-10-02
# (CLAUDE.md "Output conventions"; z = 2 is the fiducial).
Z_VALUES = (1.0, 2.0, 3.0, 5.0, 7.0)
Z_FIDUCIAL = 2.0

# Constants taken as known (Gao & Zhang 2015, section 4; Li et al. 2024 fix f_gamma = 0.5 as well).
F_GAMMA = 0.5
DELTA = 1.0 / 3.0  # the only value for which the top-down formulae exist

# Jet-base radius trials [cm]. Two sets, user decision 2026-10-24:
#   "GZ2015"          -- Gao & Zhang's own case-study values for GRB110721A (and Li+24's), 1e8, 1e9, 1e10 cm;
#   "project_anchor"  -- the pure-fireball r_0 threshold times 0.1, 1, 10. The threshold is Pe'er+07's r_0
#                        (Y = 1, median over the burst's BB episodes at the same z, evaluated in this script's flux
#                        band by ``compute_pure_fireball_anchors``) rescaled to
#                        the radiative efficiency used here (r_0 ∝ Y^(-3/2), Y = 1/f_gamma), i.e. times f_gamma^(3/2).
#                        By ``check_pure_fireball_closure()`` it is the r_0 at which this inversion returns sigma_0 = 0,
#                        so this set is a scan around the pure-fireball solution, not an independent input.
R0_GZ2015_CM = (1.0e8, 1.0e9, 1.0e10)
R0_ANCHOR_MULTIPLIERS = (0.1, 1.0, 10.0)

# Observer-frame band F_BB and F_ob are integrated over: 10 keV - 10 GeV (user decision, 2026-10-25), the same band as
# the stage 2 primary L_gamma. This is wider than pe_er_photosphere.py's 1 keV - 10 MeV (Phase 1/2, unchanged), so F_BB,
# F_ob and f_th here are *not* the quantities of the published Pe'er+07 table; the r_0 anchor is therefore recomputed
# in this band (``compute_pure_fireball_anchors``) instead of read from pe_er_photosphere.csv.
E_MIN_KEV, E_MAX_KEV = 10.0, 1.0e7
PERCENTILES = (16.0, 50.0, 84.0)

# Geometric prefactors from Pe'er et al. (2007), eqs. (1) and (5); copied from pe_er_photosphere.py (cross-folder
# import convention in CLAUDE.md), needed to evaluate their r_0 in this script's band.
PREFACTOR_R = 1.06
PREFACTOR_T = 1.48
SIGMA_SB = 5.670374419e-5  # Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]
KEV_TO_KELVIN = 1.160451812e7  # K per keV

# 1 + sigma_0 below 1 is unphysical (sigma_0 < 0). The closure test against Pe'er+07 recovers 1 + sigma_0 = 0.98 at
# sigma_0 = 0, from the rounding of the paper's coefficients, so values within this tolerance of 1 count as physical.
SIGMA_FLOOR_TOLERANCE = 0.03

SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)
# Independent child stream for the partner models, so adding them leaves the main table's draws untouched.
partner_rng = get_rng(seed=SEED).spawn(1)[0]

# The paper's blackbody-selection threshold (Delta C-stat, Delta k = 2); every partner model sits below it, which is
# why it is a partner and not BEST. Recorded per row so the separate table can say so.
BB_THRESHOLD_DELTA_CSTAT = 28.74

REGIMES = ("II", "III", "V", "VI")


# ─── Table-driven power laws, one per printed row ───────────────────────────


@dataclass(frozen=True)
class PowerLaw:
    """coefficient * (1+z)^pz * (kT/kt0)^pk * (F_BB/f0)^pf * r_0,9^pr * (f_th/0.1)^pt * f_gamma^pg * d_L,28^pd."""

    coefficient: float
    kt0: float
    f0: float
    pz: float
    pk: float
    pf: float
    pr: float
    pt: float
    pg: float
    pd: float

    def __call__(self, z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28):
        return (
            self.coefficient
            * (1.0 + z) ** self.pz
            * (kt / self.kt0) ** self.pk
            * (flux_bb / self.f0) ** self.pf
            * r0_9**self.pr
            * (f_th / 0.1) ** self.pt
            * f_gamma**self.pg
            * d_l28**self.pd
        )


# Table 2 (non-dissipation rows). Each regime lists its criteria as (power law, required sign), ">" meaning the printed
# expression must exceed 1. Regime III and VI share their first row and differ only by its sign.
_III_VI_A = PowerLaw(8.28, 30, 1e-7, -3 / 2, -3 / 2, 5 / 8, -1, 5 / 4, 5 / 4, 5 / 4)
CRITERIA = {
    "II": (
        (PowerLaw(14.8, 50, 1e-8, 1 / 4, 1 / 4, 3 / 16, -1 / 8, 1 / 2, 1 / 2, 3 / 8), ">"),
        (PowerLaw(0.24, 50, 1e-8, -3, -3, 3 / 4, -3 / 2, 0, 0, 3 / 2), ">"),
        (PowerLaw(1.43e-5, 50, 1e-8, -7, -7, 7 / 4, -7 / 2, 3, 3, 7 / 2), "<"),
    ),
    "III": (
        (_III_VI_A, ">"),
        (PowerLaw(9.42e-2, 30, 1e-7, -14 / 3, -14 / 3, 7 / 6, -7 / 3, 2, 2, 7 / 3), ">"),
    ),
    "V": (
        (PowerLaw(41.4, 10, 1e-9, 1 / 2, 1 / 2, 3 / 8, -1 / 4, 1, 1, 3 / 4), "<"),
        (PowerLaw(5.28, 10, 1e-9, -3, -3, 3 / 4, -3 / 2, 0, 0, 3 / 2), ">"),
        (PowerLaw(1.16e-5, 10, 1e-9, -8, -8, 1, -3, 1, 1, 2), "<"),
    ),
    "VI": (
        (_III_VI_A, "<"),
        (PowerLaw(5.63e-3, 30, 1e-7, -8 / 3, -8 / 3, 1 / 3, -1, 1 / 3, 1 / 3, 2 / 3), ">"),
    ),
}

# Eqs. (38) (regime II), (39) (III and VI) and (40) (V): quantity -> power law.
# eta, r_ph [cm] and Gamma_ph are as printed; "one_plus_sigma0" is 1 + sigma_0.
INVERSIONS = {
    "II": {
        "one_plus_sigma0": PowerLaw(25.5, 50, 1e-8, 4 / 3, 4 / 3, -1 / 3, 2 / 3, -1, -1, -2 / 3),
        "eta": PowerLaw(74.8, 50, 1e-8, 11 / 12, 11 / 12, 1 / 48, 5 / 24, 0, 0, 1 / 24),
        "r_ph": PowerLaw(1.78e10, 50, 1e-8, -25 / 12, -25 / 12, 37 / 48, -7 / 24, 0, 0, 37 / 24),
        "gamma_ph": PowerLaw(46.4, 50, 1e-8, -1 / 12, -1 / 12, 13 / 48, -7 / 24, 0, 0, 13 / 24),
    },
    "III": {
        "one_plus_sigma0": PowerLaw(5.99, 30, 1e-7, 4 / 3, 4 / 3, -1 / 3, 2 / 3, -1, -1, -2 / 3),
        "eta": PowerLaw(20.3, 30, 1e-7, -5 / 6, -5 / 6, 11 / 24, -2 / 3, 3 / 4, 3 / 4, 11 / 12),
        "r_ph": PowerLaw(4.09e11, 30, 1e-7, -3 / 2, -3 / 2, 5 / 8, 0, -1 / 4, -1 / 4, 5 / 4),
        "gamma_ph": PowerLaw(121.3, 30, 1e-7, 1 / 2, 1 / 2, 1 / 8, 0, -1 / 4, -1 / 4, 1 / 4),
    },
    "V": {
        "one_plus_sigma0": PowerLaw(6.43, 10, 1e-9, 4 / 3, 4 / 3, -1 / 3, 2 / 3, -1, -1, -2 / 3),
        "eta": PowerLaw(105.0, 10, 1e-9, 7 / 6, 7 / 6, 5 / 24, 1 / 12, 1 / 2, 1 / 2, 5 / 12),
        "r_ph": PowerLaw(4.62e10, 10, 1e-9, -13 / 6, -13 / 6, 17 / 24, -1 / 4, -1 / 6, -1 / 6, 17 / 12),
        "gamma_ph": PowerLaw(15.3, 10, 1e-9, -1 / 6, -1 / 6, 5 / 24, -1 / 4, -1 / 6, -1 / 6, 5 / 12),
    },
}
INVERSIONS["VI"] = INVERSIONS["III"]  # "regime VI has the identical scalings as regime III" (paper, after eq. 40)


# ─── Regime selection and inversion ──────────────────────────────────────────


def regime_masks(z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28):
    """Boolean array per regime: True where every Table 2 criterion of that regime holds (broadcasts over draws)."""
    masks = {}
    for regime, rows in CRITERIA.items():
        ok = np.ones(np.broadcast(z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28).shape, dtype=bool)
        for law, sign in rows:
            value = law(z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28)
            ok &= (value > 1.0) if sign == ">" else (value < 1.0)
        masks[regime] = ok
    return masks


def classify(masks):
    """Regime label per element: one of REGIMES, "none" (no regime passes) or "ambiguous" (more than one passes).

    Table 2's criteria carry prefactors rounded to 3 figures, so two adjacent regimes can both pass within ~1% of a
    boundary; such points are reported as ambiguous rather than assigned to either.
    """
    stack = np.stack([masks[r] for r in REGIMES])
    count = stack.sum(axis=0)
    labels = np.full(count.shape, "none", dtype=object)
    labels[count > 1] = "ambiguous"
    for index, regime in enumerate(REGIMES):
        labels[(count == 1) & stack[index]] = regime
    return labels


def invert(regime, z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28):
    """(1 + sigma_0, eta, r_ph, Gamma_ph) in one regime, eqs. (38)-(40)."""
    return {
        name: law(z, kt, flux_bb, r0_9, f_th, f_gamma, d_l28) for name, law in INVERSIONS[regime].items()
    }


# ─── Validation: bottom-up eqs. (24)-(30) against the top-down code ──────────

# Table 1 threshold, 7.22e5 L_w,52^(3/5) r_0,9^(-3/5), and the bottom-up kT_ob, F_BB, r_ph, Gamma_ph per regime.
# Written with eta_1 = eta/10, eta_2 = eta/100, (1+sigma_0)_1 = (1+sigma_0)/10, (1+sigma_0)_2 = (1+sigma_0)/100 (paper).


def bottom_up_regime(lw52, r0_9, eta, one_plus_sigma0):
    """Table 1 (non-dissipation): regime from the central-engine parameters, "I" and "IV" being degenerate."""
    threshold = 7.22e5 * lw52**0.6 * r0_9**-0.6
    s = one_plus_sigma0
    if eta > s**0.5:
        if eta**2.4 * s**0.6 > threshold:
            return "I"
        return "III" if eta**2.4 * s**3.6 < threshold else "II"
    if eta**1.2 * s**1.2 > threshold:
        return "IV"
    return "VI" if eta**3.2 * s**3.2 < threshold else "V"


def bottom_up_observables(regime, lw52, r0_9, eta, one_plus_sigma0, z, d_l28):
    """kT_ob [keV], F_BB [erg s^-1 cm^-2], r_ph [cm], Gamma_ph from eqs. (25), (26)/(30) and (29)."""
    e1, e2 = eta / 10.0, eta / 100.0
    s1, s2 = one_plus_sigma0 / 10.0, one_plus_sigma0 / 100.0
    if regime == "II":
        kt = 64.8 / (1 + z) * lw52 ** (-1 / 60) * r0_9 ** (-7 / 30) * e2 ** (16 / 15) * s2 ** (1 / 60)
        flux = 1.24e-7 * lw52 ** (11 / 15) * r0_9 ** (4 / 15) * e2 ** (16 / 15) * s2 ** (-11 / 15) / d_l28**2
        r_ph = 7.22e10 * lw52 ** (3 / 5) * r0_9 ** (2 / 5) * e2 ** (-7 / 5) * s2 ** (-3 / 5)
        gamma = 89.7 * lw52 ** (1 / 5) * r0_9 ** (-1 / 5) * e2 ** (1 / 5) * s2 ** (-1 / 5)
    elif regime in ("III", "VI"):
        kt = 6.65 / (1 + z) * lw52 ** (-5 / 12) * r0_9 ** (1 / 6) * e1 ** (8 / 3) * s1 ** (29 / 12)
        flux = 7.15e-8 * lw52 ** (1 / 3) * r0_9 ** (2 / 3) * e1 ** (8 / 3) * s1 ** (5 / 3) / d_l28**2
        r_ph = 5.81e12 * lw52 * e1**-3 * s1**-3
        gamma = 100.0 * e1 * s1
    elif regime == "V":
        kt = 19.0 / (1 + z) * lw52 ** (-1 / 60) * r0_9 ** (-7 / 30) * e2 ** (8 / 15) * s2 ** (17 / 60)
        flux = 3.63e-8 * lw52 ** (11 / 15) * r0_9 ** (4 / 15) * e2 ** (8 / 15) * s2 ** (-7 / 15) / d_l28**2
        r_ph = 2.46e11 * lw52 ** (3 / 5) * r0_9 ** (2 / 5) * e2 ** (-13 / 15) * s2 ** (-13 / 15)
        gamma = 48.5 * lw52 ** (1 / 5) * r0_9 ** (-1 / 5) * e2 ** (-1 / 15) * s2 ** (-1 / 15)
    else:
        raise ValueError(f"regime {regime} is degenerate and has no top-down inversion")
    return kt, flux, r_ph, gamma


def validate_round_trip(n_points=40_000, verbose=True, seed=20241024):
    """Recover (eta, 1+sigma_0, r_ph, Gamma_ph) from bottom-up observables, over random central-engine parameters.

    Returns {regime: (n_points, max relative error per quantity)}. Rounded 3-figure prefactors in the paper set the
    floor at a few tenths of a percent; classification disagreements within ~1% of a regime boundary are expected.
    """
    generator = np.random.default_rng(seed)
    worst = {}
    counts = {}
    mismatches = 0
    degenerate = 0
    for _ in range(n_points):
        lw52 = 10 ** generator.uniform(-1, 3.5)
        r0_9 = 10 ** generator.uniform(-1.5, 1.5)
        eta = 10 ** generator.uniform(0.3, 3.3)
        s = 10 ** generator.uniform(0, 2.7)
        z = generator.uniform(0.5, 7.0)
        d_l28 = generator.uniform(0.3, 6.0)
        regime = bottom_up_regime(lw52, r0_9, eta, s)
        if regime in ("I", "IV"):
            degenerate += 1
            continue
        kt, flux_bb, r_ph, gamma = bottom_up_observables(regime, lw52, r0_9, eta, s, z, d_l28)
        f_th = flux_bb / (F_GAMMA * lw52 * 1e52 / (4 * np.pi * (d_l28 * 1e28) ** 2))  # F_ob = f_gamma L_w / 4 pi d_L^2
        label = str(classify(regime_masks(z, kt, flux_bb, r0_9, f_th, F_GAMMA, d_l28)))
        same_family = {label, regime} == {"III", "VI"}
        if label != regime and not same_family:
            mismatches += 1
            continue
        out = invert(label, z, kt, flux_bb, r0_9, f_th, F_GAMMA, d_l28)
        errors = np.array(
            [
                abs(out["one_plus_sigma0"] / s - 1),
                abs(out["eta"] / eta - 1),
                abs(out["r_ph"] / r_ph - 1),
                abs(out["gamma_ph"] / gamma - 1),
            ]
        )
        worst[regime] = np.maximum(worst.get(regime, 0.0), errors)
        counts[regime] = counts.get(regime, 0) + 1

    if verbose:
        print("Round-trip validation (bottom-up eqs. 24-30 -> top-down eqs. 38-40, f_gamma = 0.5):")
        header = f"  {'regime':<8}{'n':>7}{'max |d(1+s0)|':>16}{'max |d eta|':>14}{'max |d r_ph|':>15}"
        print(header + f"{'max |d Gamma|':>15}")
        for regime in sorted(worst):
            e = worst[regime]
            print(f"  {regime:<8}{counts[regime]:>7}{e[0]:>16.2e}{e[1]:>14.2e}{e[2]:>15.2e}{e[3]:>15.2e}")
        print(
            f"  degenerate (Table 1 regimes I/IV, not invertible): {degenerate};"
            f" Table 1 vs Table 2 classification disagreements: {mismatches}\n"
        )
    return {r: (counts[r], worst[r]) for r in worst}, mismatches


def check_boundary_continuity(verbose=True):
    """Regimes II and III must give the same bottom-up observables on their shared boundary r_ph = r_c.

    Both are evaluated at eta^(12/5) (1+sigma_0)^(18/5) = 7.22e5 L_w,52^(3/5) r_0,9^(-3/5) (Table 1); the two
    closed forms are independent expressions, so their agreement is a real test of the transcription, not of the
    interpolation. Likewise V/VI at eta^(16/5) (1+sigma_0)^(16/5) = the same threshold.
    """
    results = {}
    # eta ranges per pair: V/VI needs eta < (1+sigma_0)^(1/2) on a boundary with eta (1+sigma_0) fixed, i.e. small eta.
    pairs = ((("II", "III"), 12 / 5, 18 / 5, (1.0, 3.0)), (("V", "VI"), 16 / 5, 16 / 5, (0.2, 1.2)))
    for pair, exponent_eta, exponent_s, (log_eta_lo, log_eta_hi) in pairs:
        worst = 0.0
        used = 0
        generator = np.random.default_rng(7)
        for _ in range(2000):
            lw52 = 10 ** generator.uniform(-1, 3)
            r0_9 = 10 ** generator.uniform(-1, 1)
            eta = 10 ** generator.uniform(log_eta_lo, log_eta_hi)
            threshold = 7.22e5 * lw52**0.6 * r0_9**-0.6
            s = (threshold / eta**exponent_eta) ** (1 / exponent_s)
            if s < 1.0 or (pair == ("II", "III") and eta <= s**0.5) or (pair == ("V", "VI") and eta >= s**0.5):
                continue
            z, d_l28 = 2.0, 1.5
            first = bottom_up_observables(pair[0], lw52, r0_9, eta, s, z, d_l28)
            second = bottom_up_observables(pair[1], lw52, r0_9, eta, s, z, d_l28)
            worst = max(worst, max(abs(a / b - 1) for a, b in zip(first, second)))
            used += 1
        results[pair] = (worst, used)
    if verbose:
        print("Regime-boundary continuity of the bottom-up observables (kT_ob, F_BB, r_ph, Gamma_ph):")
        for pair, (worst, used) in results.items():
            print(f"  {pair[0]}/{pair[1]}: max relative difference {worst:.2e} over {used} boundary points")
        print()
    return results


def check_li_r0_scaling(verbose=True):
    """1 + sigma_0 ∝ r_0^(2/3) in every regime, against Li et al. (2024) Table 1, GRB110721A first slice.

    That slice (kT = 34.419 keV, f_th = 0.0854) lists 1 + sigma_0 = 4.784 at r_0 = 1e8 cm and 22.206 at 1e9 cm.
    F_ob and z are not tabulated there, so absolute values cannot be reproduced, but the ratio is independent of them.
    """
    ratio = 22.206 / 4.784
    expected = 10 ** (2 / 3)
    if verbose:
        label = "Li+24 Table 1 (GRB110721A) r_0 scaling"
        print(f"{label}: (1+s0)[1e9]/(1+s0)[1e8] = {ratio:.4f}, 10^(2/3) = {expected:.4f}")
        print()
    return ratio, expected


def check_pure_fireball_closure(verbose=True, z_checked=(2.0, 4.35)):
    """At sigma_0 = 0 the inversion must reproduce Pe'er+07 (the project's own, independently coded calculation).

    With f_gamma = 1 (Y = 1) and r_0 set to each episode's Pe'er r_0 from pe_er_photosphere.csv, the top-down code is
    fed that episode's (kT_ob, F_BB = R^2 sigma_SB T^4, f_BB) and must return 1 + sigma_0 = 1, eta = Gamma and the same
    r_ph. Residuals of ~2% are expected: the paper's coefficients are rounded to three figures and Pe'er's use 1.06 and
    1.48 exactly. Returns the rows as a DataFrame.
    """
    sigma_sb, kev_to_kelvin = 5.670374419e-5, 1.160451812e7
    frame = load_pure_fireball_frame()
    records = []
    for _, r in frame[frame.z.isin(z_checked)].iterrows():
        kt = r.kt_bb_keV
        flux_bb = r.script_R**2 * sigma_sb * (kt * kev_to_kelvin) ** 4
        d_l28 = COSMO.luminosity_distance(r.z).cgs.value / 1e28
        args = (r.z, kt, flux_bb, r.r0_cm / 1e9, r.f_bb, 1.0, d_l28)
        label = str(classify(regime_masks(*args)))
        out = invert(label, *args) if label in REGIMES else None
        records.append(
            {
                "grb_name": r.grb_name,
                "episode": r.episode,
                "z": r.z,
                "regime": label,
                "one_plus_sigma0": float(out["one_plus_sigma0"]) if out else np.nan,
                "eta_over_Gamma_peer": float(out["eta"]) / r.Gamma if out else np.nan,
                "r_ph_over_peer": float(out["r_ph"]) / r.r_ph_cm if out else np.nan,
            }
        )
    result = pd.DataFrame(records)
    if verbose:
        print("Pure-fireball closure (f_gamma = 1, r_0 = Pe'er+07 r_0 at Y = 1): expect 1+sigma_0 = 1, ratios = 1")
        print(f"  episodes {len(result)}; regimes {result.regime.value_counts().to_dict()}")
        for column in ("one_plus_sigma0", "eta_over_Gamma_peer", "r_ph_over_peer"):
            print(f"  {column:<22} min {result[column].min():.4f}  max {result[column].max():.4f}")
        print()
    return result


# ─── Per-episode computation ─────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def interval_samples(model, n_samples=N_SAMPLES, rng=None):
    """MC draws of (F_BB, F_total, kT) for one interval, in CGS.

    Same recipe as pe_er_photosphere.py's ``interval_samples`` (copied, per the cross-folder-import convention in
    CLAUDE.md): one set of parameter draws, so the three quantities stay correlated.
    """
    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)
    samples = draw_model_samples(model, n_samples=n_samples, rng=rng)
    fluxes, total = component_energy_fluxes(model.name, samples, energy)
    flux_bb = fluxes[gmC.BB] * kev_to_erg
    flux_total = total * kev_to_erg
    kt_index = [p.name for p in model.parameters].index("kt_bb")
    kt_kev = samples[:, kt_index]
    valid = np.isfinite(flux_bb) & np.isfinite(flux_total) & (flux_bb > 0) & (kt_kev > 0) & (flux_bb < flux_total)
    return flux_bb[valid], flux_total[valid], kt_kev[valid]


def summarise(values):
    """(median, minus, plus) from the 16/50/84 percentiles."""
    lo, med, hi = np.percentile(values, PERCENTILES)
    return med, med - lo, hi - med


def load_pure_fireball_frame():
    """pe_er_photosphere.csv, asserted to be the Y = 1 calculation the anchor rescaling below assumes."""
    path = find_project_root() / "codes-for-paper" / "photospheric_radius" / "pe_er_photosphere.csv"
    frame = pd.read_csv(path)
    assert (frame.Y_ratio == 1.0).all(), "r_0 anchor rescaling assumes pe_er_photosphere.csv was computed at Y = 1"
    return frame


def compute_pure_fireball_anchors(grb_objects):
    """Median Pe'er+07 r_0 [cm] at Y = 1 per (burst, z) over the burst's BB episodes, in this script's flux band.

    Pe'er+07 eq. (5), r_0 = 4^1.5 / (1.48^6 1.06^4) * d_L / (1+z)^2 * f_BB^1.5 * R, evaluated on fresh MC draws (own
    stream, so the main table's draws are untouched) as the median over draws per episode and then the median over the
    burst's BB episodes, the same recipe pe_er_photosphere.py used for the 1 keV - 10 MeV values.
    """
    anchor_rng = get_rng(seed=SEED).spawn(2)[1]
    prefactor = 4**1.5 / (PREFACTOR_T**6 * PREFACTOR_R**4)
    per_episode = {}
    for short_name, grb in zip(GRB_LIST, grb_objects):
        grb_name = f"GRB{short_name}"
        z_grid, _ = redshift_grid(short_name)
        for model in grb.get_all_best_models():
            if "BB" not in model.name:
                continue
            flux_bb, flux_total, kt = interval_samples(model, rng=anchor_rng)
            r_value = np.sqrt(flux_bb / (SIGMA_SB * (kt * KEV_TO_KELVIN) ** 4))
            f_bb = flux_bb / flux_total
            for z in z_grid:
                d_l = COSMO.luminosity_distance(z).cgs.value
                r_zero = prefactor * d_l / (1 + z) ** 2 * f_bb**1.5 * r_value
                per_episode.setdefault((grb_name, round(float(z), 6)), []).append(float(np.median(r_zero)))
    return {key: float(np.median(values)) for key, values in per_episode.items()}


def r0_trials(grb_name, z, anchors):
    """[(set name, r_0 [cm], multiplier or NaN)] for one burst at one redshift."""
    trials = [("GZ2015", r0, np.nan) for r0 in R0_GZ2015_CM]
    key = (grb_name, round(float(z), 6))
    if key not in anchors:
        raise KeyError(f"no pure-fireball r_0 computed for {key}")
    threshold = anchors[key] * F_GAMMA**1.5  # Pe'er+07 r_0 at Y = 1 -> at Y = 1/f_gamma, see R0 comment above
    trials += [("project_anchor", threshold * m, m) for m in R0_ANCHOR_MULTIPLIERS]
    return trials


def evaluate_trial(flux_bb, flux_total, kt, z, r0_cm):
    """Regime and inversion for one (z, r_0) over the MC draws.

    The regime of the row is the one selected by the median inputs; the summaries use the draws that fall in that
    same regime, and ``regime_fraction`` reports how many that is, so a row whose draws straddle a boundary is visible
    rather than silently averaged across two formulae.
    """
    d_l28 = COSMO.luminosity_distance(z).cgs.value / 1e28
    r0_9 = r0_cm / 1e9
    f_th = flux_bb / flux_total

    median_label = str(
        classify(
            regime_masks(z, np.median(kt), np.median(flux_bb), r0_9, np.median(f_th), F_GAMMA, d_l28)
        )
    )
    draw_labels = classify(regime_masks(z, kt, flux_bb, r0_9, f_th, F_GAMMA, d_l28))
    if median_label not in REGIMES:
        return median_label, float(np.mean(draw_labels == median_label)), None

    in_regime = draw_labels == median_label
    out = invert(median_label, z, kt[in_regime], flux_bb[in_regime], r0_9, f_th[in_regime], F_GAMMA, d_l28)
    return median_label, float(np.mean(in_regime)), {name: summarise(values) for name, values in out.items()}


def trial_rows(model, grb_name, z_grid, z_source, anchors, rng_, extra=None):
    """Rows for one blackbody-bearing model: every assumed z and every r_0 trial. ``extra`` columns are appended."""
    flux_bb, flux_total, kt = interval_samples(model, rng=rng_)
    rows = []
    for z in z_grid:
        for set_name, r0_cm, multiplier in r0_trials(grb_name, z, anchors):
            regime, fraction, summary = evaluate_trial(flux_bb, flux_total, kt, z, r0_cm)
            rows.append(
                {
                    "grb_name": grb_name,
                    "episode": episode_label(model.interval),
                    "model_name": model.name,
                    "z": float(z),
                    "z_source": z_source,
                    "H0": H0,
                    "Om0": OM0,
                    "e_min_keV": E_MIN_KEV,
                    "e_max_keV": E_MAX_KEV,
                    "f_gamma": F_GAMMA,
                    "delta": DELTA,
                    "r0_set": set_name,
                    "r0_cm": r0_cm,
                    "r0_multiplier": multiplier,
                    "regime": regime,
                    "regime_fraction": fraction,
                    "physical": bool(
                        summary is not None and summary["one_plus_sigma0"][0] >= 1.0 - SIGMA_FLOOR_TOLERANCE
                    ),
                    "kt_bb_keV": float(np.median(kt)),
                    "f_th": float(np.median(flux_bb / flux_total)),
                    "summary": summary,
                    **(extra or {}),
                }
            )
    return rows


def redshift_grid(short_name):
    """(assumed or measured z values, their source) for a burst."""
    measured_z = REDSHIFTS[short_name]
    return ((measured_z,), "spectroscopic") if measured_z is not None else (Z_VALUES, "assumed")


def collect_results():
    """Gao & Zhang inversion for every BB-inclusive BEST model, every assumed z and every r_0 trial."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)
    anchors = compute_pure_fireball_anchors(grb_objects)

    rows = []
    for short_name, grb in zip(GRB_LIST, grb_objects):
        z_grid, z_source = redshift_grid(short_name)
        for model in grb.get_all_best_models():
            if "BB" in model.name:
                rows += trial_rows(model, f"GRB{short_name}", z_grid, z_source, anchors, rng)
    return rows


def partner_pairs(grb):
    """(BEST model, its X_BB partner or None) for every episode whose BEST model has no blackbody.

    The partner is the BEST continuum X plus a blackbody -- ``X_BB`` -- and nothing else (SBPL's partner is SBPL_BB, not
    BAND_BB): the comparison is within the BEST model's own continuum, not across continua (user, 2026-10-24).
    """
    for interval in grb.intervals:
        best = next((m for m in interval.models if m.status is GoodnessOfFit.BEST), None)
        if best is None or "BB" in best.name:
            continue
        yield best, interval.models.get(f"{best.name}_BB")


def collect_partner_results():
    """Inversion for the SAFE or MARGINAL ``X_BB`` partner of every non-BB BEST episode; UNSAFE or missing: skipped."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)
    anchors = compute_pure_fireball_anchors(grb_objects)

    rows = []
    for short_name, grb in zip(GRB_LIST, grb_objects):
        grb_name = f"GRB{short_name}"
        z_grid, z_source = redshift_grid(short_name)
        for best, partner in partner_pairs(grb):
            label = f"{grb_name} {episode_label(best.interval)}"
            if partner is None:
                print(f"  {label:<18} {best.name}_BB: not in results.json, skipped")
                continue
            if partner.status is GoodnessOfFit.UNSAFE:
                print(f"  {label:<18} {partner.name}: UNSAFE, excluded")
                continue
            delta = best.cstat - partner.cstat
            extra = {
                "best_model": best.name,
                "bb_status": partner.status.value,
                "delta_cstat_vs_best": float(delta),
                "delta_dof_vs_best": int(best.dof - partner.dof),
                "bb_threshold_delta_cstat": BB_THRESHOLD_DELTA_CSTAT,
                "bb_required": bool(delta >= BB_THRESHOLD_DELTA_CSTAT),
            }
            print(f"  {label:<18} {partner.name:<9} {partner.status.value:<9} dC-stat vs {best.name} = {delta:5.1f}")
            rows += trial_rows(partner, grb_name, z_grid, z_source, anchors, partner_rng, extra)
    return rows


# ─── Outputs ─────────────────────────────────────────────────────────────────

QUANTITIES = (
    ("one_plus_sigma0", "one_plus_sigma0", ""),
    ("eta", "eta", ""),
    ("r_ph", "r_ph", "_cm"),
    ("gamma_ph", "Gamma_ph", ""),
)


def write_csv(rows, path="gao_zhang_inversion.csv"):
    """One self-describing row per episode per z per r_0 trial; unsolvable trials keep their row with NaN values."""
    records = []
    for row in rows:
        record = {k: v for k, v in row.items() if k != "summary"}
        for key, column, unit in QUANTITIES:
            median, lower, upper = row["summary"][key] if row["summary"] else (np.nan, np.nan, np.nan)
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
    """Validate the transcription, then invert every BEST-BB episode and, separately, the BB partner models."""
    validate_round_trip()
    check_boundary_continuity()
    check_li_r0_scaling()
    check_pure_fireball_closure()
    print(f"Gao & Zhang (2015) top-down inversion (f_gamma = {F_GAMMA}, delta = {DELTA:.3f}, {N_SAMPLES} MC samples)\n")
    rows = collect_results()
    frame = write_csv(rows)
    solved = frame[frame.regime.isin(REGIMES)]
    print(f"\n{len(solved)} of {len(frame)} trials have a regime; regimes: {frame.regime.value_counts().to_dict()}")

    print("\nBB partner models of the episodes whose BEST model has no blackbody (reported separately):")
    partner_frame = write_csv(collect_partner_results(), path="gao_zhang_inversion_partner.csv")
    solved = partner_frame[partner_frame.regime.isin(REGIMES)]
    print(f"{len(solved)} of {len(partner_frame)} partner trials have a regime; "
          f"regimes: {partner_frame.regime.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
