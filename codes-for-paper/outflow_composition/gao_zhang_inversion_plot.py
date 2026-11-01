"""
Figure 3 — outflow magnetisation 1 + sigma_0 against the assumed jet-base radius r_0 (Gao & Zhang 2015 inversion)
================================================================================================================
Two figures drawn from the CSVs of ``gao_zhang_inversion.py``:

    gao_zhang_inversion         -- the main table: the episodes whose BEST model contains a blackbody, one panel per
                                   burst.
    gao_zhang_inversion_partner -- the SAFE/MARGINAL ``X_BB`` partner of episodes whose BEST model has no blackbody,
                                   reported separately because the data do not require the blackbody
                                   (outflow_composition.md sec. 3a).

Per episode one curve through its six r_0 trials (the two trial sets, ``GZ2015`` and ``project_anchor``, share the
axis), median and 16-84% band, at the measured redshift (GRB080916C) or the fiducial z = 2 (the other bursts: the z
dependence at fixed r_0 is modest, outflow_composition.md sec. 5). Horizontal line: 1 + sigma_0 = 1 (sigma_0 = 0, a
pure fireball); vertical dotted line: the burst's pure-fireball r_0 threshold (median over its BEST-BB episodes of
Pe'er+07 r_0 x f_gamma^1.5, the anchor of the second trial set; the same line is drawn on the partner figure). A grey
cross on the lower axis edge marks a trial r_0 at which at least one episode has no valid regime (the inversion does
not exist there).

Run ``gao_zhang_inversion.py`` first. Nothing here is random, so there is no seed.
"""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import EpisodeMarkerResolver, GRBPlotStyle, episode_order, update_style
from grb_research.grb_constants import LEGEND_FONT_SIZE, LINE_WIDTH
from grb_research.grb_utils import save_fig

# ─── Configuration ───────────────────────────────────────────────────────────

MAIN_CSV = "gao_zhang_inversion.csv"
PARTNER_CSV = "gao_zhang_inversion_partner.csv"
GRB_LIST = ["GRB080916C", "GRB131014A", "GRB140206B", "GRB231129C"]
Z_FIDUCIAL = 2.0

# Per-GRB T90 marker, matching pe_er_photosphere.py, amati_relationship.py and bb_flux_fraction.py.
T90_MARKERS = {"GRB080916C": "o", "GRB131014A": "s", "GRB140206B": "X", "GRB231129C": "D"}
EPISODE_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

# Edge width for hollow (BB-bearing) markers; not covered by grb_constants.
MARKER_EDGE_WIDTH = 1.4
# Y-axis padding: factor below the smallest and above the largest plotted value.
Y_PAD = 2.0
NONE_COLOUR = "0.45"
NONE_MARKER_LIFT = 1.35  # the cross sits this factor above the lower axis edge

# ─── Helpers ─────────────────────────────────────────────────────────────────


def interval_stub(label):
    """Minimal stand-in for a TimeInterval (``kind.name`` and ``index``), so the project's marker map can be reused."""
    if label in ("T90", "EX0", "EX1"):
        return SimpleNamespace(kind=SimpleNamespace(name=label), index=0)
    name = label.rstrip("0123456789")  # "TR" or "SP"
    return SimpleNamespace(kind=SimpleNamespace(name=name), index=int(label[len(name):]))


def panel_rows(frame, grb):
    """Rows of one burst at its measured redshift or the fiducial one."""
    subset = frame[frame.grb_name == grb]
    if subset.empty:
        return subset
    z = float(subset.z.iloc[0]) if subset.z_source.iloc[0] == "spectroscopic" else Z_FIDUCIAL
    return subset[np.isclose(subset.z, z)]


def draw_panel(axis, subset, grb, partner):
    """One burst: a curve per episode, the sigma_0 = 0 line, the r_0 threshold and the no-valid-regime crosses."""
    colour = GRBPlotStyle.GRB_COLORS[grb]
    resolver = EpisodeMarkerResolver(t90_marker=T90_MARKERS[grb])
    episodes = sorted(subset.episode.unique(), key=episode_order)

    plotted = []
    for position, episode in enumerate(episodes):
        track = subset[subset.episode == episode].sort_values("r0_cm")
        solved = track[track.one_plus_sigma0.notna()]
        model = track.model_name.iloc[0]
        label = f"{episode} ({model}" + (f", {track.bb_status.iloc[0]})" if partner else ")")
        if solved.empty:
            axis.plot([], [], color=colour, label=label + ", no valid regime")
            continue
        median = solved.one_plus_sigma0.to_numpy()
        lower = solved.one_plus_sigma0_err_lower.to_numpy()
        upper = solved.one_plus_sigma0_err_upper.to_numpy()
        axis.plot(
            solved.r0_cm, median, color=colour, linewidth=LINE_WIDTH,
            linestyle=EPISODE_LINESTYLES[position % len(EPISODE_LINESTYLES)],
            marker=resolver.resolve(interval_stub(episode)), markerfacecolor="none",
            markeredgewidth=MARKER_EDGE_WIDTH, label=label,
        )
        axis.fill_between(solved.r0_cm, median - lower, median + upper, color=colour, alpha=0.12)
        plotted += [(median - lower).min(), (median + upper).max()]

    axis.axhline(1.0, color="0.3", linestyle="-", linewidth=LINE_WIDTH, zorder=-10)
    threshold = subset[(subset.r0_set == "project_anchor") & np.isclose(subset.r0_multiplier, 1.0)].r0_cm.median()
    axis.axvline(threshold, color="0.4", linestyle=":", linewidth=LINE_WIDTH, zorder=-10)

    low = min(plotted + [1.0]) / Y_PAD
    high = max(plotted + [1.0]) * Y_PAD
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_ylim(low, high)

    no_regime = subset[subset.one_plus_sigma0.isna()].r0_cm.unique()
    if len(no_regime):
        axis.plot(no_regime, np.full(len(no_regime), low * NONE_MARKER_LIFT), linestyle="none", marker="x",
                  color=NONE_COLOUR, markeredgewidth=MARKER_EDGE_WIDTH)
    return float(subset.z.iloc[0])


def make_figure(csv_path, path_stem, partner):
    """One panel per burst present in ``csv_path``."""
    update_style()
    frame = pd.read_csv(csv_path)
    bursts = [g for g in GRB_LIST if g in set(frame.grb_name)]
    # One spare cell holds the key to the shared symbols.
    n_columns = 3 if len(bursts) > 3 else 2
    n_rows = int(np.ceil((len(bursts) + 1) / n_columns))
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(5.4 * n_columns, 4.4 * n_rows), squeeze=False)
    axes = axes.flatten()

    for axis, grb in zip(axes, bursts):
        subset = panel_rows(frame, grb)
        z = draw_panel(axis, subset, grb, partner)
        z_text = f"$z={z:g}$" + ("" if subset.z_source.iloc[0] == "spectroscopic" else ", assumed")
        axis.set_title(f"{grb} ({z_text})")
        axis.set_xlabel(r"Assumed $r_0$ [cm]")
        axis.set_ylabel(r"$1+\sigma_0$")
        axis.legend(fontsize=LEGEND_FONT_SIZE - 2, loc="upper left", frameon=True)

    for axis in axes[len(bursts):]:
        axis.axis("off")
    # The shared symbols are explained once, in the first spare cell.
    key = axes[len(bursts)]
    key.plot([], [], color="0.3", linewidth=LINE_WIDTH, label=r"$1+\sigma_0=1$ ($\sigma_0=0$)")
    key.plot(
        [], [], color="0.4", linestyle=":", linewidth=LINE_WIDTH,
        label=r"pure-fireball $r_0$ (median, BEST-BB episodes)",
    )
    key.plot([], [], color=NONE_COLOUR, linestyle="none", marker="x", markeredgewidth=MARKER_EDGE_WIDTH,
             label=r"no valid regime for $\geq1$ episode")
    key.legend(loc="center", fontsize=LEGEND_FONT_SIZE, frameon=True)
    save_fig(figure, path_stem)


def main():
    """Draw the main figure and the separately-reported partner figure."""
    make_figure(MAIN_CSV, "gao_zhang_inversion", partner=False)
    make_figure(PARTNER_CSV, "gao_zhang_inversion_partner", partner=True)


if __name__ == "__main__":
    main()
