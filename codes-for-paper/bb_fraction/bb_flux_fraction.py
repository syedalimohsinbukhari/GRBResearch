"""
Blackbody Fractional Flux Contribution
======================================
Computes f_BB = F_BB / F_total, the fraction of the observed energy flux carried
by the thermal (blackbody) component, for every interval whose BEST-fitting
model includes a BB component.

Definition
----------
Following Pe'er et al. (2007, ApJ 664, L1), the quantity that enters the
photospheric-radius / Lorentz-factor method is the ratio of the *observed*
thermal flux to the *observed* total (thermal + non-thermal) gamma-ray flux,

    f_BB = F^ob_BB / F^ob                                              (Pe'er+07)

with F^ob_BB "integrated over all frequencies". Both fluxes are therefore
observer-frame and bolometric. Pe'er et al. quote f_BB = 0.64 +/- 0.20 for
GRB 970828 at its break time, which sets the scale to expect.

Because f_BB is a ratio of two fluxes evaluated over the same band at the same
epoch, it is independent of redshift and of the assumed cosmology -- neither
enters this calculation, so neither is reported as a per-row assumption. What
*is* an assumption is the integration band, which is reported explicitly.

The band is the observer-frame 1 keV - 10 MeV used for the bolometric fluence
elsewhere in this work. The blackbody is fully contained within it: the script
checks each interval's band-integrated BB flux against the analytic all-
frequency result

    F_BB(0, inf) = amp_bb * (kT)^4 * pi^4 / 15

and reports the captured fraction, so the approximation to Pe'er's
"all frequencies" is quantified rather than assumed.

Uncertainties are propagated by Monte Carlo: parameters are drawn from the fit
covariance and pushed through the model-specific physical constraints via
ModelResampler (the same machinery mc_spectra_sampler uses), then f_BB is
recomputed per draw and summarised by percentiles.

Outputs:
    bb_flux_fraction.csv   -- one row per BB-inclusive interval
    bb_flux_fraction.png   -- f_BB vs time, one panel per GRB
    bb_flux_fraction.pdf
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import (
    EpisodeMarkerResolver,
    GRBPlotStyle,
    component_energy_fluxes,
    draw_model_samples,
    find_project_root,
    prepare_grbs,
    update_style,
)
from grb_research.grb_constants import LEGEND_FONT_SIZE, SAVE_DPI, kev_to_erg
from grb_research.grb_enums import GRBModelsCombinations as gmC

# ─── Configuration ───────────────────────────────────────────────────────────

GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

# Per-GRB T90 marker, matching amati_relationship.py so the two figures agree.
T90_MARKERS = ["o", "s", "X", "D"]

# Observer-frame integration band, log10(keV): 1 keV - 10 MeV.
E_MIN_KEV, E_MAX_KEV = 1.0, 1.0e4
N_GRID = 1_000

# Monte-Carlo draws evaluated per vectorised block. Bounds peak memory at
# roughly chunk * n_grid * 8 bytes per component array.
CHUNK_SIZE = 2_000

# MC settings — the project-wide convention (amati_relationship.py:39-41).
N_SAMPLES = 10_000
SEED = 12345

# Percentiles for the median and the 1-sigma asymmetric interval.
PERCENTILES = (16.0, 50.0, 84.0)

# Edge width for hollow (BB-augmented) markers. Not covered by grb_constants,
# which sets marker size but not edge weight.
MARKER_EDGE_WIDTH = 1.4


@dataclass
class FractionResult:
    """f_BB and supporting quantities for a single interval."""

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


# ─── Spectral integration ────────────────────────────────────────────────────


def split_energy_flux(model_name, values, energy):
    """BB and total energy flux for one or many parameter draws [keV/cm^2/s].

    Thin wrapper over ``component_energy_fluxes``; the decomposition itself
    lives in the package so Phase 2 shares exactly this code path.
    """
    fluxes, total = component_energy_fluxes(model_name, values, energy, chunk=CHUNK_SIZE)
    return fluxes[gmC.BB], total


def analytic_bb_bolometric_flux(amp_bb, kt_bb):
    """All-frequency blackbody energy flux [keV/cm^2/s].

    For N(E) = amp * E^2 / (exp(E/kT) - 1),

        int_0^inf E N(E) dE = amp * kT^4 * Gamma(4) * zeta(4) = amp * kT^4 * pi^4 / 15.
    """
    return amp_bb * kt_bb**4 * np.pi**4 / 15.0


# ─── Monte Carlo ─────────────────────────────────────────────────────────────


def compute_fraction(model, grb_name, n_samples=N_SAMPLES, seed=SEED, rng=None):
    """Compute f_BB and its MC uncertainty for one BB-inclusive model."""
    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)
    rng_instance = rng if rng is not None else np.random.default_rng(seed)

    best_values = np.array([p.value for p in model.parameters])
    best_bb, best_total = split_energy_flux(model.name, best_values, energy)
    flux_bb_best = float(best_bb[0])

    samples = draw_model_samples(model, n_samples=n_samples, rng=rng_instance)
    sample_bb, sample_total = split_energy_flux(model.name, samples, energy)
    fractions = sample_bb / sample_total

    # Physically f_BB lies in [0, 1]; draws outside it come from the Gaussian
    # tails of amp_bb and are discarded rather than clipped, so the percentiles
    # are not artificially piled up at the boundaries.
    usable = np.isfinite(fractions) & (fractions > 0) & (fractions < 1)
    fractions = fractions[usable]
    lo, med, hi = np.percentile(fractions, PERCENTILES)

    def summarise(values):
        """(median, minus, plus) from the 16/50/84 percentiles."""
        low, mid, high = np.percentile(values, PERCENTILES)
        return mid, mid - low, high - mid

    # Fluxes are summarised over the same accepted draws as f_BB, so all three
    # columns of the table describe one consistent sample.
    flux_bb = tuple(v * kev_to_erg for v in summarise(sample_bb[usable]))
    flux_total = tuple(v * kev_to_erg for v in summarise(sample_total[usable]))

    kt_bb = model.get_parameter_value("kt_bb")
    amp_bb = model.get_parameter_value("amp_bb")
    kt_err = next((p.error for p in model.parameters if p.name == "kt_bb"), np.nan)
    captured = flux_bb_best / analytic_bb_bolometric_flux(amp_bb, kt_bb)

    return FractionResult(
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


# ─── Driver ──────────────────────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def collect_results():
    """Compute f_BB for every BB-inclusive BEST model in the sample."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    rows = []
    for index, (short_name, grb) in enumerate(zip(GRB_LIST, grb_objects)):
        resolver = EpisodeMarkerResolver(t90_marker=T90_MARKERS[index])
        for model in grb.get_all_best_models():
            if "BB" not in model.name:
                continue
            result = compute_fraction(model, short_name)
            result.marker = resolver.resolve(model.interval)
            rows.append(result)
            print(
                f"  GRB{short_name:<9} {result.episode:<5} {model.name:<12} "
                f"f_BB = {result.f_bb:6.4f} -{result.f_bb_lo:6.4f} +{result.f_bb_hi:6.4f}   "
                f"(BB captured in band: {result.bb_captured * 100:5.2f}%)"
            )
    return rows


def write_csv(rows, path="bb_flux_fraction.csv"):
    """Write the results table, one self-describing row per interval."""
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
                "n_samples": N_SAMPLES,
                "seed": SEED,
            }
            for r in rows
        ]
    )
    frame.to_csv(path, index=False)
    print(f"\nSaved: {path}")
    return frame


def make_plot(rows, path_stem="bb_flux_fraction"):
    """f_BB against interval mid-time, one panel per GRB."""
    update_style()

    grbs = [f"GRB{g}" for g in GRB_LIST]
    _, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9))
    axes = axes.flatten()

    # A single y-range across all panels, so the burst-to-burst contrast in
    # f_BB is readable at a glance rather than rescaled away per panel.
    highest = max(r.f_bb + r.f_bb_hi for r in rows)
    y_top = np.ceil(highest * 20.0) / 20.0 + 0.05

    for i, (grb, axis) in enumerate(zip(grbs, axes)):
        colour = GRBPlotStyle.GRB_COLORS[grb]
        subset = [r for r in rows if r.grb_name == grb]

        for r in subset:
            mid = 0.5 * (r.t_start + r.t_stop)
            half = 0.5 * (r.t_stop - r.t_start)
            axis.errorbar(
                mid,
                r.f_bb,
                xerr=half,
                yerr=[[r.f_bb_lo], [r.f_bb_hi]],
                marker=r.marker,
                color=colour,
                # Hollow markers denote BB-augmented models, per the project
                # convention; every point here is BB-augmented by construction.
                markerfacecolor="none",
                markeredgewidth=MARKER_EDGE_WIDTH,
                linestyle="none",
                capsize=3,
                label=f"{r.episode} ({r.model_name})",
            )

        axis.set_title(grb)
        axis.set_ylim(0, y_top)
        if subset:
            axis.legend(loc="upper right", ncols=1, fontsize=LEGEND_FONT_SIZE)

    for axis in axes[2:]:
        axis.set_xlabel(r"Time since $T_0$ [s]")
    for axis in axes[::2]:
        axis.set_ylabel(r"$f_\mathrm{BB} = F_\mathrm{BB}/F_\mathrm{total}$")

    plt.tight_layout()
    for extension in ("png", "pdf"):
        plt.savefig(f"./{path_stem}.{extension}", dpi=SAVE_DPI)
    plt.close()
    print(f"Saved: {path_stem}.png / .pdf")


def main():
    print("Computing blackbody flux fractions "
          f"({E_MIN_KEV:g}-{E_MAX_KEV:g} keV, observer frame, {N_SAMPLES} MC samples)\n")
    rows = collect_results()
    write_csv(rows)
    make_plot(rows)


if __name__ == "__main__":
    main()
