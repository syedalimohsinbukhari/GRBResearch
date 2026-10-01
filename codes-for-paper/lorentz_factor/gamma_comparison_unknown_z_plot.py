"""
Gamma_min (Limit A, Limit B) vs. thermal Gamma over redshift -- the three bursts without a measured redshift
==============================================================================================================
Redshift-sweep counterpart of ``gamma_comparison_plot.py`` (GRB080916C, one episode per x position).  For
GRB131014A, GRB140206B and GRB231129C no redshift is measured, so each panel puts redshift on the x-axis, and all three
quantities are curves over z = 0.5-7 with their 1-sigma bands:

    - the thermal Gamma (Pe'er 2007, Y = 1) from ``pe_er_photosphere.csv``;
    - Limit A and Limit B from ``lorentz_curves_unknown_z.csv`` / ``lorentz_curves_limit_b_unknown_z.csv`` -- the same
      calculation as the tables, on a dense z grid from each episode's existing draw set.  The markers sit on the curves at
      the four tabulated redshifts z = 1, 3, 5, 7 (currently disabled -- lines only) and are exactly the values in ``tab:lorentz_unknown_z`` /
      ``tab:lorentz_limit_b_unknown_z``.

Only episodes with a blackbody component (the only ones with a thermal Gamma) are shown, so every episode drawn has all
three quantities; the non-thermal-only episodes are in the tables.

This script computes nothing new: it reads CSVs produced by ``lorentz_factor_unknown_z.py``,
``lorentz_factor_limit_b_unknown_z.py`` and ``../photospheric_radius/pe_er_photosphere.py``.  Episode markers come from
the live ``TimeInterval`` objects via ``EpisodeMarkerResolver``, as in ``gamma_comparison_plot.py``.  No Monte Carlo of
its own, so no seed.

Outputs:
    gamma_comparison_unknown_z.png
    gamma_comparison_unknown_z.pdf
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import EpisodeMarkerResolver, episode_order, find_project_root, prepare_grbs, update_style
from grb_research.grb_constants import (
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LINE_WIDTH,
    MARKER_SIZE,
    TITLE_FONT_SIZE,
)
from grb_research.grb_utils import save_fig

from lorentz_factor import episode_label

GRB_LIST = ["131014A", "140206B", "231129C"]
Z_VALUES = (1, 3, 5, 7)

ROOT = find_project_root()
LORENTZ_DIR = ROOT / "codes-for-paper" / "lorentz_factor"
PHOTOSPHERIC_CSV = ROOT / "codes-for-paper" / "photospheric_radius" / "pe_er_photosphere.csv"

# Same method colours as gamma_comparison_plot.py (Okabe-Ito), so the two figures read as one family.
METHOD_COLORS = {"limit_a": "#0072B2", "limit_b": "#D55E00", "thermal": "#009E73"}
METHOD_LABELS = {
    "limit_a": "$\\Gamma_{\\min}$ (Limit A)",
    "limit_b": "$\\Gamma_{\\min,B}$ (Limit B)",
    "thermal": "$\\Gamma$ (thermal), Pe'er 2007, $Y{=}1$",
}
# Line style is keyed by the episode label, not by position within a burst, so the same episode has the same style in every
# panel (a position-based list shifts GRB140206B's EX0/TR1 by one style, as it has no T90 blackbody episode).
EPISODE_LINESTYLES = {"T90": "-", "EX0": "--", "TR1": "-.", "TR2": ":", "EX1": (0, (5, 1, 1, 1, 1, 1))}
MARKER_EDGE_WIDTH = 1.4
T90_MARKER = "o"


def load_data():
    """Read the three CSVs (limit curves on a dense z grid, thermal curves), restricted to thermally-detected episodes."""
    limit_a = pd.read_csv(LORENTZ_DIR / "lorentz_curves_unknown_z.csv")
    limit_b = pd.read_csv(LORENTZ_DIR / "lorentz_curves_limit_b_unknown_z.csv")
    thermal = pd.read_csv(PHOTOSPHERIC_CSV)

    names = [f"GRB{s}" for s in GRB_LIST]
    thermal = thermal[thermal.grb_name.isin(names) & (thermal.Y_ratio == 1.0)]
    limit_a = limit_a[limit_a.GRB.isin(names)]
    limit_b = limit_b[limit_b.GRB.isin(names)]

    # Keep only (burst, episode) pairs that have a thermal Gamma.
    keys = set(zip(thermal.grb_name, thermal.episode))
    limit_a = limit_a[[(g, e) in keys for g, e in zip(limit_a.GRB, limit_a.episode)]]
    limit_b = limit_b[[(g, e) in keys for g, e in zip(limit_b.GRB, limit_b.episode)]]
    return limit_a, limit_b, thermal


def collect_intervals():
    """Real TimeInterval objects per (burst, episode), for marker identity only."""
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=ROOT / "results.json", get_best=True)
    return {
        (f"GRB{short}", episode_label(m.interval)): m.interval
        for short, grb in zip(GRB_LIST, grb_objects)
        for m in grb.get_all_best_models()
    }


def make_plot(limit_a, limit_b, thermal, path_stem="gamma_comparison_unknown_z"):
    """One panel per burst (2x2 grid, the fourth axes holding the method legend): thermal Gamma curves over z, Limit A/B markers at the assumed z."""
    update_style()

    figure, grid = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
    axes = grid.ravel()
    legend_axis = axes[3]  # no fourth burst: this axes carries the method legend
    intervals = collect_intervals()
    resolver = EpisodeMarkerResolver(t90_marker=T90_MARKER)
    models = {(r.GRB, r.episode): r.model for r in limit_a.drop_duplicates(["GRB", "episode"]).itertuples()}

    for axis, short in zip(axes[:3], GRB_LIST):
        grb = f"GRB{short}"
        episodes = sorted(thermal[thermal.grb_name == grb].episode.unique(), key=episode_order)
        handles = []

        for episode in episodes:
            marker = resolver.resolve(intervals[(grb, episode)])  # only used by the disabled markers below
            style = EPISODE_LINESTYLES[episode]

            hollow = "BB" in models[(grb, episode)].upper()
            handles.append(
                mlines.Line2D(
                    [], [], color="#444444", linestyle=style, linewidth=LINE_WIDTH,
                    # marker=marker,
                    # markersize=MARKER_SIZE, markerfacecolor="none" if hollow else "#444444",
                    # markeredgecolor="#444444", markeredgewidth=MARKER_EDGE_WIDTH,
                    label=f"{episode}" + r"$_\text{" + models[(grb, episode)].replace("_", "+") + r"}$",
                )
            )

            curve = thermal[(thermal.grb_name == grb) & (thermal.episode == episode)].sort_values("z")
            axis.plot(curve.z, curve.Gamma, color=METHOD_COLORS["thermal"], linewidth=LINE_WIDTH, linestyle=style)
            axis.fill_between(
                curve.z,
                curve.Gamma - curve.Gamma_err_lower,
                curve.Gamma + curve.Gamma_err_upper,
                color=METHOD_COLORS["thermal"],
                alpha=0.12,
            )

            for key, frame, column in (("limit_a", limit_a, "Gamma_min_A"), ("limit_b", limit_b, "Gamma_min_B")):
                rows = frame[(frame.GRB == grb) & (frame.episode == episode)].sort_values("z")
                if rows.empty:
                    continue
                color = METHOD_COLORS[key]
                axis.plot(rows.z, rows[column], color=color, linewidth=LINE_WIDTH, linestyle=style)
                axis.fill_between(
                    rows.z, rows[column] - rows[f"{column}_err_lower"], rows[column] + rows[f"{column}_err_upper"],
                    color=color, alpha=0.12,
                )
                # Lines only for now: the markers at the four tabulated redshifts (exact curve points, identifying the episode)
                # are disabled, not removed -- uncomment this block and the marker kwargs in the legend handle below to restore.
                # at_table = rows[np.isclose(rows.z.to_numpy()[:, None], Z_VALUES).any(axis=1)]
                # axis.plot(
                #     at_table.z, at_table[column], linestyle="none", marker=marker, markersize=MARKER_SIZE,
                #     markerfacecolor="none" if hollow else color, markeredgecolor=color,
                #     markeredgewidth=MARKER_EDGE_WIDTH, zorder=3,
                # )

        # Legend per panel (episode sets differ by burst), inside its own frame, in the empty low-Gamma corner.
        axis.legend(
            handles=handles, loc="best", fontsize=LEGEND_FONT_SIZE, title="BEST model",
            title_fontsize=LEGEND_FONT_SIZE, frameon=True,
        )
        axis.set_yscale("log")
        # axis.set_ylim(4, 3000)
        # axis.set_xlim(0.3, 7.7)
        axis.set_xticks(range(1, 8, 2))
        axis.xaxis.minorticks_off()
        axis.set_title(f"GRB{short}", fontsize=TITLE_FONT_SIZE)

    # The top-right panel has no axes below it (the fourth one is legend-only), so it carries its own x tick labels/label.
    axes[1].tick_params(labelbottom=True)
    for axis in (axes[1], axes[2]):
        axis.set_xlabel("assumed redshift $z$", fontsize=LABEL_FONT_SIZE)
    for axis in (axes[0], axes[2]):
        axis.set_ylabel(r"$\Gamma$", fontsize=LABEL_FONT_SIZE)

    method_handles = [
        mlines.Line2D(
            [], [], color=METHOD_COLORS[k], linestyle="none" if k != "thermal" else "-",
            marker="o" if k != "thermal" else None, markerfacecolor=METHOD_COLORS[k], markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH, label=METHOD_LABELS[k],
        )
        for k in ("limit_a", "limit_b", "thermal")
    ]
    legend_axis.axis("off")
    legend_axis.legend(handles=method_handles, loc="center", fontsize=LEGEND_FONT_SIZE, frameon=True)
    save_fig(figure, path_stem)


def main():
    limit_a, limit_b, thermal = load_data()
    make_plot(limit_a, limit_b, thermal)


if __name__ == "__main__":
    main()
