"""
Photospheric Radius and Bulk Lorentz Factor — Pe'er et al. (2007) method
========================================================================
Derives the size at the base of the flow r_0, the outflow bulk Lorentz factor
Gamma, the photospheric radius r_ph and the saturation radius r_s from the
observed temperature and thermal flux of the blackbody component.

Reference
---------
Pe'er, Ryde, Wijers, Meszaros & Rees (2007), ApJ 664, L1,
"A new method of determining the initial size and Lorentz factor of gamma-ray
burst fireballs using a thermal emission component".

Equations, in their numbering:

    R = (F_BB^ob / (sigma T_ob^4))^(1/2) = 1.06 (1+z)^2 r_ph / (d_L Gamma)   (1)

    r_0 (r_ph < r_s) = [1 / 1.06] [d_L / (1+z)^2] R                          (2)

    eta = [1.06 (1+z)^2 d_L Y F^ob sigma_T / (2 m_p c^3 R)]^(1/4)            (4)

    r_0 (r_ph > r_s) = 4^(3/2) / (1.48^6 1.06^4) [d_L / (1+z)^2]
                       (F_BB^ob / (Y F^ob))^(3/2) R                          (5)

with Gamma = eta in the coasting phase, r_s = eta r_0, and
Y = epsilon / epsilon_gamma >= 1 the ratio of total fireball energy to the
energy emitted in gamma-rays.

Y is not measurable from the prompt spectrum alone, so — exactly as Pe'er et al.
do — results are quoted at Y = 1 with the scaling made explicit:

    Gamma ∝ Y^(1/4)     r_0 ∝ Y^(-3/2)     r_ph ∝ Y^(1/4)     r_s ∝ Y^(-5/4)

Validation
----------
``validate_against_peer2007()`` reproduces their worked example for GRB 970828
from the published inputs. Run it before trusting any output here.

Redshift
--------
Only GRB080916C has a spectroscopic redshift (z = 4.35). For the other three
bursts z is swept over a log grid spanning the bulk of the long-GRB
distribution, with a fiducial value highlighted, so the redshift sensitivity is
shown rather than an unmeasured value being invented.

Outputs:
    pe_er_photosphere.csv   -- one row per interval per redshift
    pe_er_photosphere.png   -- r_0(z) and Gamma(z), with the fixed-z burst marked
    pe_er_photosphere.pdf
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from grb_research import (
    EpisodeMarkerResolver,
    GRBPlotStyle,
    component_energy_fluxes,
    draw_model_samples,
    find_project_root,
    get_rng,
    prepare_grbs,
    seed_from_name,
    update_style,
)
from grb_research.grb_constants import LEGEND_FONT_SIZE, LINE_WIDTH, SAVE_DPI, kev_to_erg, N_SAMPLES, N_GRID
from grb_research.grb_enums import GRBModelsCombinations as gmC

# ─── Physical constants (CGS) ────────────────────────────────────────────────

SIGMA_SB = 5.670374419e-5  # Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]
SIGMA_T = 6.6524587321e-25  # Thomson cross section [cm^2]
M_P = 1.67262192369e-24  # proton mass [g]
C_LIGHT = 2.99792458e10  # speed of light [cm/s]
KEV_TO_KELVIN = 1.160451812e7  # K per keV

# Geometric prefactors from Pe'er et al. (2007); see their eqs. (1) and (5).
PREFACTOR_R = 1.06
PREFACTOR_T = 1.48

# ─── Cosmology ───────────────────────────────────────────────────────────────
# Fana Dirirsa et al. (2019) — the single cosmology used across this project.
H0 = 69.6
OM0 = 0.286
COSMO = FlatLambdaCDM(H0=H0, Om0=OM0)

# ─── Configuration ───────────────────────────────────────────────────────────

GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

# Spectroscopic redshifts; None means unmeasured and therefore swept.
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Redshift sweep for the bursts without a measured z:
# the bulk of the long-GRB distribution, with z = 2 (close to the long-GRB median) as the fiducial.
Z_MIN, Z_MAX, Z_POINTS = 0.5, 5.0, 25
Z_FIDUCIAL = 2.0

# Ratio of total fireball energy to energy radiated in gamma-rays. Not
# measurable here, so fixed at unity and reported as an explicit scaling.
Y_RATIO = 1.0

# Per-GRB T90 marker, matching amati_relationship.py and bb_flux_fraction.py.
T90_MARKERS = ["o", "s", "X", "D"]

# Line style per episode, so overlapping curves of one burst stay separable.
EPISODE_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

# Edge width for hollow (BB-augmented) markers.
# Not covered by grb_constants, which sets marker size but not edge weight.
MARKER_EDGE_WIDTH = 1.4

# Observer-frame integration band and MC settings, matching Phase 1 exactly so
# f_BB here is the same quantity as in bb_flux_fraction.py.
E_MIN_KEV, E_MAX_KEV = 1.0, 1.0e4
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)
PERCENTILES = (16.0, 50.0, 84.0)


@dataclass
class PhotosphereResult:
    """Derived fireball parameters for one interval at one redshift."""

    grb_name: str
    episode: str
    model_name: str
    interval_kind: object  # the TimeInterval, for per-episode marker resolution
    z: float
    z_source: str
    kt_bb: tuple
    f_bb: tuple
    script_r: float
    gamma: tuple
    r_zero: tuple
    r_ph: tuple
    r_sat: tuple
    saturated_fraction: float


# ─── Pe'er et al. (2007) relations ───────────────────────────────────────────


def script_r(flux_bb, kt_kev):
    """R = (F_BB / sigma T^4)^(1/2), Pe'er+07 eq. (1). Dimensionless."""
    temperature_k = kt_kev * KEV_TO_KELVIN
    return np.sqrt(flux_bb / (SIGMA_SB * temperature_k ** 4))


def lorentz_factor(flux_total, r_value, z, d_l, y_ratio=Y_RATIO):
    """Gamma = eta, Pe'er+07 eq. (4). Coasting phase, so Gamma equals eta."""
    numerator = PREFACTOR_R * (1 + z) ** 2 * d_l * y_ratio * flux_total * SIGMA_T
    return (numerator / (2 * M_P * C_LIGHT ** 3 * r_value)) ** 0.25


def base_radius_unsaturated(f_bb, r_value, z, d_l, y_ratio=Y_RATIO):
    """r_0 for r_ph > r_s, Pe'er+07 eq. (5) [cm]."""
    prefactor = 4 ** 1.5 / (PREFACTOR_T ** 6 * PREFACTOR_R ** 4)
    return prefactor * d_l / (1 + z) ** 2 * (f_bb / y_ratio) ** 1.5 * r_value


def base_radius_saturated(r_value, z, d_l):
    """r_0 for r_ph < r_s, Pe'er+07 eq. (2) [cm]."""
    return d_l * r_value / (PREFACTOR_R * (1 + z) ** 2)


def photospheric_radius(r_value, gamma, z, d_l):
    """r_ph, inverted from Pe'er+07 eq. (1) [cm]."""
    return r_value * d_l * gamma / (PREFACTOR_R * (1 + z) ** 2)


def validate_against_peer2007(verbose=True):
    """Reproduce Pe'er+07's GRB 970828 worked example from their published inputs.

    Their §3 quotes T_ob = 78.5 keV, R = 1.88e-19, F_BB/F_ob = 0.64,
    d_L = 1.94e28 cm at z = 0.9578, and reports
    Gamma = 305 +/- 28, r_0 = (2.9 +/- 1.8)e8 cm, r_ph = 2.7e11 cm,
    r_s = 9.0e10 cm (all at Y = 1).

    Returns the computed (gamma, r_0, r_ph, r_s).
    """
    z, kt, r_value, f_bb, d_l = 0.9578, 78.5, 1.88e-19, 0.64, 1.94e28
    flux_bb = r_value ** 2 * SIGMA_SB * (kt * KEV_TO_KELVIN) ** 4
    flux_total = flux_bb / f_bb

    gamma = lorentz_factor(flux_total, r_value, z, d_l)
    r_zero = base_radius_unsaturated(f_bb, r_value, z, d_l)
    r_ph = photospheric_radius(r_value, gamma, z, d_l)
    r_sat = gamma * r_zero

    if verbose:
        print("Validation against Pe'er et al. (2007), GRB 970828:")
        print(f"  {'quantity':<12}{'computed':>13}{'published':>16}")
        print(f"  {'Gamma':<12}{gamma:>13.0f}{'305 +/- 28':>16}")
        print(f"  {'r_0 [cm]':<12}{r_zero:>13.3e}{'(2.9+/-1.8)e8':>16}")
        print(f"  {'r_ph [cm]':<12}{r_ph:>13.3e}{'2.7e11':>16}")
        print(f"  {'r_s [cm]':<12}{r_sat:>13.3e}{'9.0e10':>16}")
        print()

    return gamma, r_zero, r_ph, r_sat


# ─── Per-interval computation ────────────────────────────────────────────────


def interval_samples(model, n_samples=N_SAMPLES, rng=None):
    """MC draws of (F_BB, F_total, kT) for one interval, in CGS.

    All three come from the *same* parameter draws, so their errors stay
    correlated when combined into R, Gamma and r_0.
    """
    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)
    samples = draw_model_samples(model, n_samples=n_samples, rng=rng)

    fluxes, total = component_energy_fluxes(model.name, samples, energy)
    flux_bb = fluxes[gmC.BB] * kev_to_erg
    flux_total = total * kev_to_erg

    kt_index = [p.name for p in model.parameters].index("kt_bb")
    kt_kev = samples[:, kt_index]

    valid = np.isfinite(flux_bb) & np.isfinite(flux_total) & (flux_bb > 0) & (kt_kev > 0)
    valid &= flux_bb < flux_total
    return flux_bb[valid], flux_total[valid], kt_kev[valid]


def summarise(values):
    """(median, minus, plus) from the 16/50/84 percentiles."""
    lo, med, hi = np.percentile(values, PERCENTILES)
    return med, med - lo, hi - med


def evaluate_at_redshift(flux_bb, flux_total, kt_kev, z):
    """Derived fireball quantities for a redshift, over the MC draws."""
    d_l = COSMO.luminosity_distance(z).cgs.value
    r_value = script_r(flux_bb, kt_kev)
    f_bb = flux_bb / flux_total

    gamma = lorentz_factor(flux_total, r_value, z, d_l)
    r_zero = base_radius_unsaturated(f_bb, r_value, z, d_l)
    r_ph = photospheric_radius(r_value, gamma, z, d_l)
    r_sat = gamma * r_zero

    # Fraction of draws in the saturated regime, where eq. (5) does not apply
    # and only an upper bound on r_0 is available via eq. (2).
    saturated = float(np.mean(r_ph < r_sat))

    return r_value, f_bb, gamma, r_zero, r_ph, r_sat, saturated


def collect_results():
    """Compute the Pe'er quantities for every BB-inclusive BEST model."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    rows = []
    for short_name, grb in zip(GRB_LIST, grb_objects):
        measured_z = REDSHIFTS[short_name]
        if measured_z is not None:
            z_grid, z_source = np.array([measured_z]), "spectroscopic"
        else:
            z_grid = np.logspace(np.log10(Z_MIN), np.log10(Z_MAX), Z_POINTS)
            z_grid = np.unique(np.append(z_grid, Z_FIDUCIAL))
            z_source = "swept"

        for model in grb.get_all_best_models():
            if "BB" not in model.name:
                continue

            flux_bb, flux_total, kt_kev = interval_samples(model, rng=rng)

            for z in z_grid:
                r_value, f_bb, gamma, r_zero, r_ph, r_sat, saturated = evaluate_at_redshift(
                    flux_bb, flux_total, kt_kev, z
                )
                rows.append(
                    PhotosphereResult(
                        grb_name=f"GRB{short_name}",
                        episode=episode_label(model.interval),
                        model_name=model.name,
                        interval_kind=model.interval,
                        z=float(z),
                        z_source=z_source,
                        kt_bb=summarise(kt_kev),
                        f_bb=summarise(f_bb),
                        script_r=float(np.median(r_value)),
                        gamma=summarise(gamma),
                        r_zero=summarise(r_zero),
                        r_ph=summarise(r_ph),
                        r_sat=summarise(r_sat),
                        saturated_fraction=saturated,
                    )
                )

            label = episode_label(model.interval)
            reference = [r for r in rows if r.grb_name.endswith(short_name) and r.episode == label]
            shown = min(reference, key=lambda r: abs(r.z - (measured_z or Z_FIDUCIAL)))
            print(
                f"  GRB{short_name:<9} {shown.episode:<5} {model.name:<12} z={shown.z:<5.2f} "
                f"Gamma={shown.gamma[0]:6.1f}  r_0={shown.r_zero[0]:.3e} cm  "
                f"r_ph={shown.r_ph[0]:.3e} cm  saturated={shown.saturated_fraction * 100:.0f}%"
            )

    return rows


def episode_order(label):
    """Sort key putting episodes in temporal order: T90, EX0, TR1..TRn, EX1."""
    fixed = {"T90": 0, "EX0": 1, "EX1": 90}
    if label in fixed:
        return fixed[label]
    if label.startswith("TR"):
        return 10 + int(label[2:])
    return 99


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


# ─── Outputs ─────────────────────────────────────────────────────────────────


def write_csv(rows, path="pe_er_photosphere.csv"):
    """Write one self-describing row per interval per redshift."""
    frame = pd.DataFrame(
        [
            {
                "grb_name": r.grb_name,
                "episode": r.episode,
                "model_name": r.model_name,
                "z": r.z,
                "z_source": r.z_source,
                "H0": H0,
                "Om0": OM0,
                "Y_ratio": Y_RATIO,
                "e_min_keV": E_MIN_KEV,
                "e_max_keV": E_MAX_KEV,
                "kt_bb_keV": r.kt_bb[0],
                "kt_bb_err_lower_keV": r.kt_bb[1],
                "kt_bb_err_upper_keV": r.kt_bb[2],
                "f_bb": r.f_bb[0],
                "f_bb_err_lower": r.f_bb[1],
                "f_bb_err_upper": r.f_bb[2],
                "script_R": r.script_r,
                "Gamma": r.gamma[0],
                "Gamma_err_lower": r.gamma[1],
                "Gamma_err_upper": r.gamma[2],
                "r0_cm": r.r_zero[0],
                "r0_err_lower_cm": r.r_zero[1],
                "r0_err_upper_cm": r.r_zero[2],
                "r_ph_cm": r.r_ph[0],
                "r_ph_err_lower_cm": r.r_ph[1],
                "r_ph_err_upper_cm": r.r_ph[2],
                "r_s_cm": r.r_sat[0],
                "r_s_err_lower_cm": r.r_sat[1],
                "r_s_err_upper_cm": r.r_sat[2],
                "saturated_fraction": r.saturated_fraction,
                "n_samples": N_SAMPLES,
                "seed": SEED,
            }
            for r in rows
        ]
    )
    frame.to_csv(path, index=False)
    print(f"\nSaved: {path}  ({len(frame)} rows)")
    return frame


def make_plot(rows, path_stem="pe_er_photosphere"):
    """r_0(z) and Gamma(z): swept bursts as curves, the fixed-z burst as points."""
    update_style()

    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(12.5, 5.0))

    for axis, key, axis_label in ((axes[0], "r_zero", r"$r_0$ [cm]"), (axes[1], "gamma", r"$\Gamma$")):
        for short_name in GRB_LIST:
            grb = f"GRB{short_name}"
            colour = GRBPlotStyle.GRB_COLORS[grb]
            subset = [r for r in rows if r.grb_name == grb]
            if not subset:
                continue

            resolver = EpisodeMarkerResolver(t90_marker=T90_MARKERS[GRB_LIST.index(short_name)])
            episodes = sorted({r.episode for r in subset}, key=episode_order)
            for position, episode in enumerate(episodes):
                track = sorted((r for r in subset if r.episode == episode), key=lambda r: r.z)
                z_values = np.array([r.z for r in track])
                medians = np.array([getattr(r, key)[0] for r in track])
                lower = np.array([getattr(r, key)[1] for r in track])
                upper = np.array([getattr(r, key)[2] for r in track])

                marker = resolver.resolve(track[0].interval_kind)
                m_name = r"$_\text{" + f'{track[0].model_name.replace("_", "+")}' + r"}$"
                series_label = f"{grb} {episode}{m_name}"

                if track[0].z_source == "spectroscopic":
                    # Single measured redshift: one point per episode.
                    axis.errorbar(
                        z_values,
                        medians,
                        yerr=[lower, upper],
                        marker=marker,
                        markerfacecolor="none",
                        markeredgewidth=MARKER_EDGE_WIDTH,
                        color=colour,
                        linestyle="none",
                        capsize=3,
                        label=series_label,
                    )
                else:
                    # Swept redshift: a curve per episode, distinguished by line style, with a marker at the fiducial z
                    # so the episode is identifiable without tracing the line.
                    style = EPISODE_LINESTYLES[position % len(EPISODE_LINESTYLES)]
                    fiducial = int(np.argmin(np.abs(z_values - Z_FIDUCIAL)))
                    # The marker must ride on the *labeled* artist, otherwise the legend handle is a bare line,
                    # and the per-episode marker -- which is what identifies TR1 vs. TR2 vs. EX0 -- never appears.
                    axis.plot(
                        z_values,
                        medians,
                        color=colour,
                        linewidth=LINE_WIDTH,
                        linestyle=style,
                        marker=marker,
                        markevery=[fiducial],
                        markerfacecolor="none",
                        markeredgewidth=MARKER_EDGE_WIDTH,
                        label=series_label,
                    )
                    axis.fill_between(z_values, medians - lower, medians + upper, color=colour, alpha=0.12)

        axis.axvline(Z_FIDUCIAL, color="0.4", linestyle=":", linewidth=LINE_WIDTH, zorder=-10)
        axis.set_xscale("log")
        axis.set_yscale("log")
        # A log axis over 0.5-5 otherwise labels only the single decade tick.
        z_ticks = [0.5, 1.0, 2.0, 3.0, 5.0]
        axis.set_xticks(z_ticks, minor=False)
        axis.set_xticklabels([f"{t:g}" for t in z_ticks])
        axis.set_xticks([], minor=True)
        axis.set_xlabel(r"Redshift $z$")
        axis.set_ylabel(axis_label)

    axes[0].set_title(r"Base radius $r_0$ (at $Y=1$; $r_0 \propto Y^{-3/2}$)")
    axes[1].set_title(r"Lorentz factor $\Gamma$ (at $Y=1$; $\Gamma \propto Y^{1/4}$)")

    # One shared legend in a single column outside the panels, on the left:
    # with 13 series an in-axes legend unavoidably covers data (it hid the
    # GRB140206B track entirely), and a single column keeps the burst/episode
    # grouping readable top-to-bottom.
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="center left", ncols=1, fontsize=LEGEND_FONT_SIZE, frameon=True, bbox_to_anchor=(0.05, 0.5)
    )

    plt.tight_layout(rect=(0.235, 0, 1, 1))

    for extension in ("png", "pdf"):
        plt.savefig(f"./{path_stem}.{extension}", dpi=SAVE_DPI)
    plt.close()
    print(f"Saved: {path_stem}.png / .pdf")


def main():
    """Validate against Pe'er+2007, then compute, tabulate and plot the sample."""
    validate_against_peer2007()
    print(f"Computing Pe'er+2007 fireball parameters (Y = {Y_RATIO:g}, {N_SAMPLES} MC samples)\n")
    rows = collect_results()
    write_csv(rows)
    make_plot(rows)


if __name__ == "__main__":
    main()
