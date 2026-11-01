"""
Diagnostic figures for the Zhang & Pe'er (2009) suppression bound
==================================================================
Two figures that show what ``zhang_peer_suppression.py`` computes, drawn from its CSV and the stored fits:

    zhang_peer_epoch_b             -- the reproduction check: Abdo et al. (2009) Table 1's epoch-b Band spectrum of
                                      GRB 080916C, the fireball blackbody at sigma = 0 and at sigma_min, for two choices
                                      of the L_gamma band.
    zhang_peer_suppression_spectra -- one panel per episode whose BEST model has no blackbody: the best-fit nu F_nu, the
                                      fireball blackbody at sigma = 0 and at the median sigma_min of the primary
                                      configuration.

Both are validation figures (not yet in the paper). The best-fit curves use the stored best-fit parameters; the
blackbody curves use the median kT_max, L_gamma and 1 + sigma_min of the CSV, so a curve can differ from the per-draw
solutions by ordinary Monte Carlo spread. Nothing here is random, so there is no seed.

Run ``zhang_peer_suppression.py`` first (reads its CSV).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import GRBPlotStyle, episode_order, find_project_root, prepare_grbs, update_style
from grb_research.grb_constants import LEGEND_FONT_SIZE, LINE_WIDTH, kev_to_erg
from grb_research.grb_seds import band_function
from grb_research.grb_utils import save_fig

import zhang_peer_suppression as zp

# ─── Configuration ───────────────────────────────────────────────────────────

CSV_PATH = "zhang_peer_suppression.csv"
PRIMARY_L_BAND = "zp09"
PRIMARY_TEST_BAND = "fit"
Z_FIDUCIAL = 2.0  # panel redshift for the bursts without a measured z

# Plot window [keV]: the test bands (1 or 10 keV to 10 MeV) with a margin either side.
E_PLOT_MIN_KEV, E_PLOT_MAX_KEV = 0.5, 2.0e4
N_PLOT = 600

# Okabe-Ito colours for the two L_gamma bands of the epoch-b figure (distinct from the burst colours).
L_BAND_COLOURS = {"zp09": "#D55E00", "rest_1_1e4": "#56B4E9"}
L_BAND_LABELS = {"zp09": r"$L_\gamma$: 10 keV$-$10 GeV", "rest_1_1e4": r"$L_\gamma$: rest-frame 1$-$10$^4$ keV"}
# Drawn thick and thin so both suppressed curves stay visible where they nearly coincide.
L_BAND_WIDTHS = {"zp09": LINE_WIDTH, "rest_1_1e4": 2.5 * LINE_WIDTH}
BB_ZERO_COLOUR = "0.45"
MARGIN_DEX = 4.0  # panel y-range: this many decades below the best-fit peak
EDGE_LIMITED_FRACTION = 0.5

# ─── Helpers ─────────────────────────────────────────────────────────────────


def energy_grid():
    """Observer-frame plotting grid [keV]."""
    return np.logspace(np.log10(E_PLOT_MIN_KEV), np.log10(E_PLOT_MAX_KEV), N_PLOT)


def blackbody_curve(energy, kt_max, flux0, one_plus_sigma):
    """nu F_nu [erg cm^-2 s^-1] of the fireball blackbody suppressed by 1 + sigma."""
    shrink = 1.0 / one_plus_sigma
    return zp.blackbody_nufnu(energy, kt_max * shrink**0.25, flux0 * shrink)


def episode_axis(axis, energy, nufnu_obs, curves, test_lo):
    """Common axis formatting for the spectrum panels."""
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.axvline(test_lo, color="0.6", linestyle=":", linewidth=LINE_WIDTH, zorder=-10)
    top = max(nufnu_obs.max(), *(c.max() for c in curves))
    axis.set_ylim(nufnu_obs.max() * 10.0**-MARGIN_DEX, 3.0 * top)
    axis.set_xlim(energy[0], energy[-1])


# ─── Figure 1: epoch-b reproduction ──────────────────────────────────────────


def make_epoch_b_figure(path_stem="zhang_peer_epoch_b"):
    """Abdo et al. epoch b with the fireball blackbody at sigma = 0 and sigma_min, for both L_gamma bands."""
    update_style()
    z = zp.EPOCH_B_Z
    ab = zp.ABDO_EPOCH_B
    values = np.array([[ab["amp"], ab["e_peak"], ab["alpha"], ab["beta"]]])
    energy = energy_grid()
    observed = zp.observed_nufnu("band", values, energy)[0]

    figure, axis = plt.subplots(figsize=(7, 6.6))
    colour = GRBPlotStyle.GRB_COLORS["GRB080916C"]
    axis.plot(
        energy, observed, color=colour, linewidth=2 * LINE_WIDTH,
        label=r"GRB080916C epoch b (Abdo et al. 2009, Band)",
    )

    r0 = zp.base_radius_cm(zp.EPOCH_B_DT_OB_S, z)
    test_grid = np.logspace(1, 4, zp.N_COMPARE)
    test_obs = zp.observed_nufnu("band", values, test_grid)
    curves = [observed]
    for l_band, tint in L_BAND_COLOURS.items():
        l_gamma = float(zp.isotropic_luminosity("band", values, z, l_band)[0])
        kt_max = zp.kt_max_kev(l_gamma, r0, z)
        flux0 = zp.bolometric_flux(l_gamma, z)
        one_plus, binding = zp.solve_one_plus_sigma(np.array([kt_max]), np.array([flux0]), test_obs, test_grid)
        sigma = float(one_plus[0]) - 1.0

        unsuppressed = blackbody_curve(energy, kt_max, flux0, 1.0)
        suppressed = blackbody_curve(energy, kt_max, flux0, float(one_plus[0]))
        curves += [unsuppressed, suppressed]
        axis.plot(energy, unsuppressed, color=tint, linestyle="--", linewidth=LINE_WIDTH)
        axis.plot(
            energy, suppressed, color=tint, linestyle="-", linewidth=L_BAND_WIDTHS[l_band],
            zorder=3 if l_band == "zp09" else 2,
            label=(
                rf"{L_BAND_LABELS[l_band]}: $kT_\mathrm{{max}}$ = {kt_max:.0f} keV, "
                rf"$\sigma_\mathrm{{min}}$ = {sigma:.1f}"
            ),
        )
        axis.plot(
            test_grid[binding[0]], suppressed[np.argmin(np.abs(energy - test_grid[binding[0]]))],
            marker="o", color=tint, markerfacecolor="none", zorder=4,
        )

    axis.plot([], [], color="0.3", linestyle="--", linewidth=LINE_WIDTH, label=r"fireball blackbody, $\sigma=0$")
    axis.plot(
        [], [], color="0.3", linestyle="-", linewidth=LINE_WIDTH, label=r"suppressed, $\sigma=\sigma_\mathrm{min}$"
    )
    episode_axis(axis, energy, observed, curves, test_lo=zp.COMPARISON_BANDS[PRIMARY_TEST_BAND][0])
    axis.set_xlabel("Energy [keV]")
    axis.set_ylabel(r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    axis.set_title(
        rf"$z={z}$, $\delta t_\mathrm{{ob}}={zp.EPOCH_B_DT_OB_S}$ s; Zhang & Pe'er: "
        rf"$kT_\mathrm{{max}}\approx{zp.EPOCH_B_KT_MAX_KEV:.0f}$ keV, $\sigma\approx{zp.EPOCH_B_SIGMA:.0f}$"
    )
    axis.legend(fontsize=LEGEND_FONT_SIZE, loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=True)
    save_fig(figure, path_stem)


# ─── Figure 2: one panel per non-BB episode ──────────────────────────────────


def best_fit_models():
    """{(grb_name, episode): model} for every BEST model with no blackbody."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=zp.GRB_LIST, result_file=root / "results.json", get_best=True)
    models = {}
    for short_name, grb in zip(zp.GRB_LIST, grb_objects):
        for model in grb.get_all_best_models():
            if "BB" not in model.name:
                models[(f"GRB{short_name}", zp.episode_label(model.interval))] = model
    return models


def make_spectra_figure(path_stem="zhang_peer_suppression_spectra"):
    """Best-fit spectrum with the fireball blackbody at sigma = 0 and sigma_min, per non-BB episode."""
    update_style()
    frame = pd.read_csv(CSV_PATH)
    frame = frame[(frame.l_band == PRIMARY_L_BAND) & (frame.band == PRIMARY_TEST_BAND)]
    models = best_fit_models()

    keys = sorted(models, key=lambda k: (k[0], episode_order(k[1])))
    n_columns = 4
    n_rows = int(np.ceil((len(keys) + 1) / n_columns))
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(3.4 * n_columns, 2.9 * n_rows))
    axes = axes.flatten()
    energy = energy_grid()
    test_lo = zp.COMPARISON_BANDS[PRIMARY_TEST_BAND][0]

    for axis, (grb, episode) in zip(axes, keys):
        model = models[(grb, episode)]
        subset = frame[(frame.grb_name == grb) & (frame.episode == episode)]
        z_panel = float(subset.z.iloc[0]) if subset.z_source.iloc[0] == "spectroscopic" else Z_FIDUCIAL
        row = subset.iloc[(subset.z - z_panel).abs().argmin()]

        values = np.array([[p.value for p in model.parameters]])
        observed = zp.observed_nufnu(model.name, values, energy)[0]
        flux0 = zp.bolometric_flux(row.L_gamma_erg_s, row.z)
        unsuppressed = blackbody_curve(energy, row.kT_max_keV, flux0, 1.0)
        suppressed = blackbody_curve(energy, row.kT_max_keV, flux0, row.one_plus_sigma_min)

        axis.plot(energy, observed, color=GRBPlotStyle.GRB_COLORS[grb], linewidth=2 * LINE_WIDTH)
        axis.plot(energy, unsuppressed, color=BB_ZERO_COLOUR, linestyle="--", linewidth=LINE_WIDTH)
        axis.plot(energy, suppressed, color="k", linestyle="-", linewidth=LINE_WIDTH)
        episode_axis(axis, energy, observed, [unsuppressed, suppressed], test_lo)

        edge = row.binding_at_lower_edge_fraction > EDGE_LIMITED_FRACTION
        axis.set_title(f"{grb} {episode} ({model.name}), $z$ = {row.z:g}", fontsize=LEGEND_FONT_SIZE)
        axis.text(
            0.97, 0.95,
            rf"$1+\sigma_\mathrm{{min}}$ = {row.one_plus_sigma_min:.1f}" + ("\n(edge-limited)" if edge else ""),
            transform=axis.transAxes, fontsize=LEGEND_FONT_SIZE - 1, va="top", ha="right",
        )

    for axis in axes[len(keys):]:
        axis.axis("off")
    legend_axis = axes[len(keys)]
    legend_axis.plot([], [], color="0.2", linewidth=2 * LINE_WIDTH, label="best-fit spectrum (colour: burst)")
    legend_axis.plot(
        [], [], color=BB_ZERO_COLOUR, linestyle="--", linewidth=LINE_WIDTH, label=r"fireball blackbody, $\sigma=0$"
    )
    legend_axis.plot([], [], color="k", linewidth=LINE_WIDTH, label=r"suppressed, $\sigma=\sigma_\mathrm{min}$")
    legend_axis.plot(
        [], [], color="0.6", linestyle=":", linewidth=LINE_WIDTH, label=f"test-band edge ({test_lo:g} keV)"
    )
    legend_axis.legend(loc="center", fontsize=LEGEND_FONT_SIZE, frameon=True)

    for axis in axes[: len(keys)]:
        axis.set_xlabel("Energy [keV]")
    for axis in axes[: len(keys) : n_columns]:
        axis.set_ylabel(r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    save_fig(figure, path_stem)


def main():
    """Draw both validation figures."""
    make_epoch_b_figure()
    make_spectra_figure()


if __name__ == "__main__":
    main()
