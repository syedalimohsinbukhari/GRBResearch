r"""
Blackbody Fractional Flux Contribution
======================================
Computes f_BB = F_BB / F_total, the fraction of the observed energy flux carried by the thermal (blackbody) component,
for every interval whose BEST-fitting model includes a BB component.

Definition
----------
Following Pe'er et al. (2007, ApJ 664, L1), the quantity that enters the photospheric-radius / Lorentz-factor method is
the ratio of the *observed* thermal flux to the *observed* total (thermal and non-thermal) gamma-ray flux,

    .. math::

        f_BB = F^ob_BB / F^ob → (Pe'er+07)

With :math:`F^ob_BB` "integrated over all frequencies".
Both fluxes are therefore observer-frame and bolometric.
Pe'er et al. quote :math:`f_BB = 0.64 +/- 0.20` for GRB970828 at its break time, which sets the scale to expect.

Because :math:`f_BB` is a ratio of two fluxes evaluated over the same band at the same epoch, it is independent of
redshift and of the assumed cosmology *within one episode* -- neither enters the observer-frame calculation below.
This does *not* make :math:`f_BB` comparable *across* bursts at different redshifts: a fixed observer-frame band is a
different rest-frame band for each burst, and since the BB peaks near 3.92 kT, this changes how much of the BB's
rest-frame shape falls inside vs outside a shared comparison window.
So a second, rest-frame :math:`f_BB` is also computed here (1 keV - 10 MeV rest-frame, matching the rest-frame band
already used for :math:`S_bol/E_iso` elsewhere in this work), using each burst's actual or assumed redshift --
z = 4.35 (measured) for GRB080916C, swept over [0.5, 5.0] with z = 2 marked as the fiducial for the other three,
mirroring the precedent already established for this exact problem in photospheric_radius/`pe_er_photosphere.py`.
See bb_fraction.md Sec 2.2b.

The observer-frame band is 1 keV - 10 MeV, matching the bolometric band used for fluence elsewhere in this work.
The blackbody is fully contained within it: the script checks each interval's band-integrated BB flux (in both bands)
against the analytic all-frequency result

    .. math::
        F_BB(0, \infty) = A_bb \frac{(kT\pi)^4}{15}

And reports the captured fraction, so the approximation to Pe'er's "all frequencies" is quantified rather than assumed.

Uncertainty propagation via Monte-Carlo method:
Parameters are drawn from the fit covariance and pushed through the model-specific physical constraints via
ModelResampler (the same machinery mc_spectra_sampler uses), then f_BB is recomputed per draw and summarized by
percentiles.
The same draws are reused for the observer-frame and every rest-frame/redshift evaluation, so the two are directly
comparable rather than adding independent MC noise.

Outputs:
    bb_flux_fraction.csv -- one row per BB-inclusive interval per redshift
    bb_flux_fraction.png / .pdf -- f_BB (observer-frame) vs. time, one panel per GRB
    bb_flux_fraction_rest_vs_z.png / .pdf -- f_BB (rest-frame) vs. redshift
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotez import plot_errorband, ErrorBandConfig

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
    TITLE_FONT_SIZE,
)
from grb_research.grb_constants import LEGEND_FONT_SIZE, LINE_WIDTH, kev_to_erg, N_SAMPLES, N_GRID
from grb_research.grb_enums import GRBModelsCombinations as gmC
from grb_research.grb_utils import save_fig

# --- Configuration -----------------------------------------------------------

GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

# Per-GRB T90 marker, matching amati_relationship.py so the two figures agree.
T90_MARKERS = ["o", "s", "X", "D"]

# Line style per episode in the rest-vs-z plot, so overlapping curves of one
# burst stay separable (matches pe_er_photosphere.py's convention).
EPISODE_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

# Observer-frame integration band, fixed, log10(keV): 1 keV - 10 MeV.
E_MIN_KEV, E_MAX_KEV = 1.0, 1.0e4

# Rest-frame integration band, log10(keV): 1 keV - 10 MeV, matching the rest-frame band already used for
# S_bol / E_iso elsewhere in this work.
E_MIN_KEV_REST, E_MAX_KEV_REST = 1.0, 1.0e4

# Spectroscopic redshifts; None means unmeasured and therefore swept.
# Matches photospheric_radius/pe_er_photosphere.py's REDSHIFTS exactly, since this is the same problem (three of four
# bursts lack a spectroscopic z).
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Redshift sweep for the bursts without a measured z:
# The bulk of the long-GRB distribution, with z = 2 (close to the long-GRB median) as the fiducial -- same grid as
# pe_er_photosphere.py, so the two are comparable.
Z_MIN, Z_MAX, Z_POINTS = 0.5, 5.0, 25
Z_FIDUCIAL = 2.0

# Cosmology -- Fana Dirirsa et al. (2019), the single cosmology used across this project.
# Inert for the observer-frame columns (f_BB is redshift/cosmology independent per-episode there); reported per
# CLAUDE.md now that z is a genuine input for the rest-frame columns.
H0 = 69.6
OM0 = 0.286

# Monte-Carlo draws evaluated per vectorized block.
# Bounds peak memory at roughly chunk * n_grid * 8 bytes per component array.
CHUNK_SIZE = 2_000

# MC settings — the project-wide convention (amati_relationship.py:39-41).
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)

# Percentiles for the median and the 1-sigma asymmetric interval.
PERCENTILES = (16.0, 50.0, 84.0)

# Edge width for hollow (BB-augmented) markers.
# Not covered by grb_constants, which sets marker size but not edge weight.
MARKER_EDGE_WIDTH = 1.4


@dataclass
class FractionResult:
    """f_BB (observer- and rest-frame, at one redshift) for a single interval."""

    grb_name: str
    episode: str
    marker: str
    model_name: str
    t_start: float
    t_stop: float
    kt_bb: float
    kt_bb_err: float
    f_bb: float
    f_bb_lo: float
    f_bb_hi: float
    flux_bb: tuple
    flux_total: tuple
    bb_captured: float
    z: float
    z_source: str
    f_bb_rest: float
    f_bb_rest_lo: float
    f_bb_rest_hi: float
    bb_captured_rest: float


# --- Spectral integration ----------------------------------------------------


def split_energy_flux(model_name, values, energy):
    """BB and total energy flux for one or many parameter draws [keV/cm^2/s].

    Thin wrapper over ``component_energy_fluxes``; the decomposition itself lives in the package, so Phase 2 shares
    exactly this code path.
    """
    fluxes, total = component_energy_fluxes(model_name, values, energy, chunk=CHUNK_SIZE)
    return fluxes[gmC.BB], total


def analytic_bb_bolometric_flux(amp_bb, kt_bb):
    r"""All-frequency blackbody energy flux [keV/cm^2/s].

    For :math:`N(E) = A * E^2 / (\exp(E/kT) - 1)`, we get

    .. math::

        \int_0^\infty E N(E) dE = A (kT)^4\Gamma(4)\zeta(4) = A \frac{(\pi kT)^4}{15}.
    """
    return amp_bb * kt_bb ** 4 * np.pi ** 4 / 15.0


# --- Monte Carlo -------------------------------------------------------------


def compute_fraction(model, grb_name, n_samples=N_SAMPLES, *, rng):
    """Compute f_BB (observer- and rest-frame, at every z) for one BB-inclusive model.

    Returns a list of ``FractionResult``, one per redshift in this burst's z-grid (a single spectroscopic z, or the
    swept grid).
    The observer-frame values are identical across the list -- they don't depend on z -- only the rest-frame values and
    the z/z_source columns vary.
    """
    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)

    best_values = np.array([p.value for p in model.parameters])
    best_bb, best_total = split_energy_flux(model.name, best_values, energy)
    flux_bb_best = float(best_bb[0])

    samples = draw_model_samples(model, n_samples=n_samples, rng=rng)
    sample_bb, sample_total = split_energy_flux(model.name, samples, energy)
    fractions = sample_bb / sample_total

    # Physically f_BB lies in [0, 1]; draws outside it come from the Gaussian tails of amp_bb and are discarded rather
    # than clipped, so the percentiles are not artificially piled up at the boundaries.
    usable = np.isfinite(fractions) & (fractions > 0) & (fractions < 1)
    lo, med, hi = np.percentile(fractions[usable], PERCENTILES)

    def summarise(values):
        """(median, minus, plus) from the 16/50/84 percentiles."""
        low, mid, high = np.percentile(values, PERCENTILES)
        return mid, mid - low, high - mid

    # Fluxes are summarized over the same accepted draws as f_BB, so all three columns of the table describe one
    # consistent sample.
    flux_bb = tuple(v * kev_to_erg for v in summarise(sample_bb[usable]))
    flux_total = tuple(v * kev_to_erg for v in summarise(sample_total[usable]))

    kt_bb = model.get_parameter_value("kt_bb")
    amp_bb = model.get_parameter_value("amp_bb")
    kt_err = next((p.error for p in model.parameters if p.name == "kt_bb"), np.nan)
    bb_bolometric = analytic_bb_bolometric_flux(amp_bb, kt_bb)
    captured = flux_bb_best / bb_bolometric

    common = dict(
        grb_name=f"GRB{grb_name}",
        episode=episode_label(model.interval),
        marker="",  # assigned by the caller, which knows the GRB's T90 marker
        model_name=model.name,
        t_start=model.interval.start,
        t_stop=model.interval.end,
        kt_bb=kt_bb,
        kt_bb_err=kt_err,
        f_bb=med,
        f_bb_lo=med - lo,
        f_bb_hi=hi - med,
        flux_bb=flux_bb,
        flux_total=flux_total,
        bb_captured=captured,
    )

    measured_z = REDSHIFTS[grb_name]
    if measured_z is not None:
        z_grid, z_source = np.array([measured_z]), "spectroscopic"
    else:
        z_grid = np.logspace(np.log10(Z_MIN), np.log10(Z_MAX), Z_POINTS)
        z_grid = np.unique(np.append(z_grid, Z_FIDUCIAL))
        z_source = "swept"

    # The rest-frame band [E_MIN_KEV_REST, E_MAX_KEV_REST] is fixed; the fitted (observer-frame) spectrum must be
    # evaluated at the observed energies that correspond to it, E_rest / (1 + z).
    # Since logspace(a, b) / (1 + z) == logspace(a - s, b - s) with s = log10(1 + z), the shift is applied to the
    # exponents directly, mirroring mc_e_iso_sampler (grb_calculations.py:394-396) -- pairing a rest-frame grid
    # with an observed-frame variable underestimated E_iso by ~3.7x at z=4.35 once already (BUGS.md, BUG-11);
    # the same care applies here.
    #
    # Every z's shifted grid is stacked into one (n_z, N_GRID) array and integrated in a single batched call rather
    # than one call per z: component_energy_fluxes/_component_energy_flux_block accept a 2D energy grid for exactly
    # this, verified bit-identical to looping per z (see that function's docstring). This is the dominant cost of
    # this script -- profiling showed ~96% of compute_fraction's time was re-running the full N_SAMPLES x N_GRID
    # integration once per swept z; batching removes the per-call Python/chunking overhead of doing that N_z times.
    redshift_shifts = np.log10(1 + z_grid)
    e_rest_as_observed_2d = np.stack(
        [
            np.logspace(np.log10(E_MIN_KEV_REST) - shift, np.log10(E_MAX_KEV_REST) - shift, N_GRID)
            for shift in redshift_shifts
        ]
    )

    rest_best_bb_2d, _ = split_energy_flux(model.name, best_values, e_rest_as_observed_2d)
    captured_rest_all = rest_best_bb_2d[0] / bb_bolometric  # shape (n_z,)

    rest_sample_bb_2d, rest_sample_total_2d = split_energy_flux(model.name, samples, e_rest_as_observed_2d)
    rest_fractions_2d = rest_sample_bb_2d / rest_sample_total_2d  # shape (n_samples, n_z)

    results = []
    for i, z in enumerate(z_grid):
        rest_fractions = rest_fractions_2d[:, i]
        usable_rest = np.isfinite(rest_fractions) & (rest_fractions > 0) & (rest_fractions < 1)
        rest_lo, rest_med, rest_hi = np.percentile(rest_fractions[usable_rest], PERCENTILES)

        results.append(
            FractionResult(
                z=float(z),
                z_source=z_source,
                f_bb_rest=rest_med,
                f_bb_rest_lo=rest_med - rest_lo,
                f_bb_rest_hi=rest_hi - rest_med,
                bb_captured_rest=float(captured_rest_all[i]),
                **common,
            )
        )

    return results


# --- Driver ------------------------------------------------------------------


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def collect_results():
    """Compute f_BB for every BB-inclusive BEST model in the sample, at every redshift."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    rows = []
    for index, (short_name, grb) in enumerate(zip(GRB_LIST, grb_objects)):
        resolver = EpisodeMarkerResolver(t90_marker=T90_MARKERS[index])
        for model in grb.get_all_best_models():
            if "BB" not in model.name:
                continue
            results = compute_fraction(model, short_name, rng=rng)
            marker = resolver.resolve(model.interval)
            for result in results:
                result.marker = marker
            rows.extend(results)

            shown = min(results, key=lambda r: abs(r.z - (REDSHIFTS[short_name] or Z_FIDUCIAL)))
            print(
                f"  GRB{short_name:<9} {shown.episode:<5} {model.name:<12} z={shown.z:<5.2f} ({shown.z_source})   "
                f"f_BB(obs) = {shown.f_bb:6.4f} -{shown.f_bb_lo:6.4f} +{shown.f_bb_hi:6.4f}   "
                f"f_BB(rest) = {shown.f_bb_rest:6.4f} -{shown.f_bb_rest_lo:6.4f} +{shown.f_bb_rest_hi:6.4f}"
            )
    return rows


def _reference_rows(rows):
    """One row per (grb, episode, model): the spectroscopic z, or the row nearest the fiducial z.

    Observer-frame values don't depend on z, so any single row per group carries the correct f_bb/flux_bb/bb_captured
    for the time-series plot; this only avoids re-plotting the same point once per swept z.
    Mirrors pe_er_photosphere.py's table-collapse logic.
    """
    groups: dict = {}
    for r in rows:
        groups.setdefault((r.grb_name, r.episode, r.model_name), []).append(r)
    return [min(group, key=lambda r: abs(r.z - Z_FIDUCIAL)) for group in groups.values()]


def write_csv(rows, path="bb_flux_fraction.csv"):
    """Write the results table, one self-describing row per interval per redshift."""
    frame = pd.DataFrame(
        [
            {
                "grb_name": r.grb_name,
                "episode": r.episode,
                "model_name": r.model_name,
                "t_start_s": r.t_start,
                "t_stop_s": r.t_stop,
                "e_min_keV": E_MIN_KEV,
                "e_max_keV": E_MAX_KEV,
                "frame": "observer",
                "kt_bb_keV": r.kt_bb,
                "kt_bb_err_keV": r.kt_bb_err,
                "f_bb": r.f_bb,
                "f_bb_err_lower": r.f_bb_lo,
                "f_bb_err_upper": r.f_bb_hi,
                "flux_bb_erg_cm2_s": r.flux_bb[0],
                "flux_bb_err_lower_erg_cm2_s": r.flux_bb[1],
                "flux_bb_err_upper_erg_cm2_s": r.flux_bb[2],
                "flux_total_erg_cm2_s": r.flux_total[0],
                "flux_total_err_lower_erg_cm2_s": r.flux_total[1],
                "flux_total_err_upper_erg_cm2_s": r.flux_total[2],
                "bb_fraction_captured_in_band": r.bb_captured,
                "z": r.z,
                "z_source": r.z_source,
                "H0": H0,
                "Om0": OM0,
                "e_min_keV_rest": E_MIN_KEV_REST,
                "e_max_keV_rest": E_MAX_KEV_REST,
                "f_bb_rest": r.f_bb_rest,
                "f_bb_rest_err_lower": r.f_bb_rest_lo,
                "f_bb_rest_err_upper": r.f_bb_rest_hi,
                "bb_fraction_captured_in_band_rest": r.bb_captured_rest,
                "n_samples": N_SAMPLES,
                "seed": SEED,
            }
            for r in rows
        ]
    )
    frame.to_csv(path, index=False)
    print(f"\nSaved: {path}  ({len(frame)} rows)")
    return frame


def make_plot(rows, path_stem="bb_flux_fraction"):
    """f_BB (observer-frame) against interval mid-time, one panel per GRB."""
    update_style()

    reference = _reference_rows(rows)
    grbs = [f"GRB{g}" for g in GRB_LIST]
    _, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9))
    axes = axes.flatten()

    # A single y-range across all panels, so the burst-to-burst contrast in f_BB is readable at a glance rather than
    # rescaled away per panel.
    highest = max(r.f_bb + r.f_bb_hi for r in reference)
    y_top = np.ceil(highest * 20.0) / 20.0 + 0.05

    for i, (grb, axis) in enumerate(zip(grbs, axes)):
        colour = GRBPlotStyle.GRB_COLORS[grb]
        subset = [r for r in reference if r.grb_name == grb]

        for r in subset:
            r: FractionResult
            mid = 0.5 * (r.t_start + r.t_stop)
            half = 0.5 * (r.t_stop - r.t_start)
            m_name = r"$_\text{" + f'{r.model_name.replace("_", "+")}' + r"}$"
            if r.episode != "T90":
                axis.errorbar(
                    mid,
                    r.f_bb,
                    xerr=half,
                    yerr=[[r.f_bb_lo], [r.f_bb_hi]],
                    marker=r.marker,
                    color=colour,
                    # Hollow markers denote BB-augmented models, per the project convention;
                    # construction BB-augments every point here.
                    markerfacecolor="none",
                    markeredgewidth=MARKER_EDGE_WIDTH,
                    linestyle="none",
                    capsize=3,
                    label=f"{r.episode}{m_name}",
                )
            else:
                plot_errorband(
                    [r.t_start, r.t_stop],
                    [r.f_bb, r.f_bb],
                    r.f_bb - r.f_bb_hi,
                    r.f_bb + r.f_bb_lo,
                    line=True,
                    line_config={"c": "k", "ls": "--", "lw": 2},
                    band_config=ErrorBandConfig(color="k", alpha=0.15),
                    data_label=f"{r.episode}{m_name}",
                    axis=axis,
                    plot_title="",
                    x_label="",
                    y_label="",
                )

        # axis.set_title(grb)
        axis.set_ylim(0, y_top)
        if subset:
            axis.legend(
                loc="upper right", title=grb, ncols=1, fontsize=LEGEND_FONT_SIZE, title_fontsize=TITLE_FONT_SIZE
            )

    for axis in axes[2:]:
        axis.set_xlabel(r"Time since $T_0$ [s]")
    for axis in axes[::2]:
        axis.set_ylabel(r"$f_\mathrm{BB} = F_\mathrm{BB}/F_\mathrm{total}$")

    save_fig(_, path_stem)


def make_rest_vs_z_plot(rows, path_stem="bb_flux_fraction_rest_vs_z"):
    r""":math:`f_BB` (rest-frame) vs. redshift: swept bursts as curves, GRB080916C as fixed points.

    Mirrors pe_er_photosphere.py's r_0(z)/Gamma(z) sweep plot -- same grid, same fiducial-marker convention --
    to show the full redshift dependence for the three bursts without a spectroscopic z, rather than committing to one
    assumed value.
    """
    update_style()
    figure, axis = plt.subplots(figsize=(10.0, 6.5))

    for index, short_name in enumerate(GRB_LIST):
        grb = f"GRB{short_name}"
        colour = GRBPlotStyle.GRB_COLORS[grb]
        subset = [r for r in rows if r.grb_name == grb]
        if not subset:
            continue

        episodes = sorted({r.episode for r in subset})
        for position, episode in enumerate(episodes):
            track = sorted((r for r in subset if r.episode == episode), key=lambda r: r.z)
            z_values = np.array([r.z for r in track])
            medians = np.array([r.f_bb_rest for r in track])
            lower = np.array([r.f_bb_rest_lo for r in track])
            upper = np.array([r.f_bb_rest_hi for r in track])

            first = track[0]
            m_name = r"$_\text{" + f'{first.model_name.replace("_", "+")}' + r"}$"
            series_label = f"{grb} {episode}{m_name}"

            if first.z_source == "spectroscopic":
                # Single measured redshift: one point per episode.
                axis.errorbar(
                    z_values,
                    medians,
                    yerr=[lower, upper],
                    marker=first.marker,
                    markerfacecolor="none",
                    markeredgewidth=MARKER_EDGE_WIDTH,
                    color=colour,
                    linestyle="none",
                    capsize=3,
                    label=series_label,
                )
            else:
                # Swept redshift: a curve per episode, distinguished by line style, with a marker at the fiducial z
                # so the episode is identifiable without tracing the line. The marker must ride on the *labeled*
                # artist, otherwise the legend handle is a bare line and the per-episode marker never appears.
                style = EPISODE_LINESTYLES[position % len(EPISODE_LINESTYLES)]
                fiducial = int(np.argmin(np.abs(z_values - Z_FIDUCIAL)))
                axis.plot(
                    z_values,
                    medians,
                    color=colour,
                    linewidth=LINE_WIDTH,
                    linestyle=style,
                    marker=first.marker,
                    markevery=[fiducial],
                    markerfacecolor="none",
                    markeredgewidth=MARKER_EDGE_WIDTH,
                    label=series_label,
                )
                axis.fill_between(z_values, medians - lower, medians + upper, color=colour, alpha=0.12)

    axis.axvline(Z_FIDUCIAL, color="0.4", linestyle=":", linewidth=LINE_WIDTH, zorder=-10)
    axis.set_xscale("log")
    z_ticks = [0.5, 1.0, 2.0, 3.0, 5.0]
    axis.set_xticks(z_ticks, minor=False)
    axis.set_xticklabels([f"{t:g}" for t in z_ticks])
    axis.set_xticks([], minor=True)
    axis.set_xlabel(r"Redshift $[z]$")
    axis.set_ylabel(r"$f_\mathrm{BB}^\mathrm{rest} = F_\mathrm{BB}/F_\mathrm{total}$" + "\nRest-frame band")

    axis.legend(loc="upper right", fontsize=LEGEND_FONT_SIZE, bbox_to_anchor=(1.5, 0.85))

    save_fig(figure, path_stem)
    print(f"Saved: {path_stem}.png / .pdf")


def main():
    print(
        "Computing blackbody flux fractions "
        f"(observer {E_MIN_KEV:g}-{E_MAX_KEV:g} keV; rest-frame {E_MIN_KEV_REST:g}-{E_MAX_KEV_REST:g} keV, "
        f"z swept {Z_MIN:g}-{Z_MAX:g} for unmeasured-z bursts, fiducial z={Z_FIDUCIAL:g}; {N_SAMPLES} MC samples)\n"
    )
    rows = collect_results()
    write_csv(rows)
    make_plot(rows)
    make_rest_vs_z_plot(rows)


if __name__ == "__main__":
    main()
