"""
Gamma_min (Limit A, Limit B) vs. thermal Gamma -- comparison plot
=====================================================================
Per-episode comparison of the two gamma-gamma-opacity lower limits on the bulk Lorentz factor (Lithwick & Sari 2001, Limits A and B) against the model-dependent thermal Gamma from the Pe'er (2007) photospheric method, for GRB080916C -- the only burst
with a confirmed redshift, and so the only one for which any of the three quantities can be computed.

This script does not compute anything new: it reads the three CSVs already produced by ``lorentz_factor.py``, ``lorentz_factor_limit_b.py`` and ``../photospheric_radius/pe_er_photosphere.py``, each already a self-contained, defensible source in its own right (CLAUDE.md's CSV standard).

Episode markers and ordering are taken from the *live* ``TimeInterval`` objects (via ``prepare_grbs``), not reparsed from the CSV's string labels.
The marker shapes are always the ones ``EpisodeMarkerResolver`` itself assigns -- no duplicated dispatch logic to drift out of sync with it.

Outputs:
    gamma_comparison.png
    gamma_comparison.pdf
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd

from grb_research import EpisodeMarkerResolver, find_project_root, prepare_grbs, update_style
from grb_research.grb_constants import LABEL_FONT_SIZE, LEGEND_FONT_SIZE, LINE_WIDTH, MARKER_SIZE, SAVE_DPI, \
    TITLE_FONT_SIZE, CAP_SIZE

from lorentz_factor import episode_label

GRB_SHORT = "080916C"
GRB_NAME = f"GRB{GRB_SHORT}"

ROOT = find_project_root()
LORENTZ_DIR = ROOT / "codes-for-paper" / "lorentz_factor"
PHOTOSPHERIC_CSV = ROOT / "codes-for-paper" / "photospheric_radius" / "pe_er_photosphere.csv"

# Method identity -> colour (Okabe-Ito, colourblind-safe). Deliberately not
# GRBPlotStyle.GRB_COLORS: this figure shows a single burst, so colour here
# encodes *method*, not GRB.
METHOD_COLORS = {
    "limit_a": "#0072B2",  # blue
    "limit_b": "#D55E00",  # vermillion
    "thermal": "#009E73",  # bluish green
}
METHOD_LABELS = {
    "limit_a": "$\\Gamma_{\\min}$ (Limit A)\n$\\gamma\\gamma$ annihilation",
    "limit_b": "$\\Gamma_{\\min,B}$ (Limit B)\nCompton on $e^{\\pm}$",
    "thermal": "$\\Gamma$ (thermal)\nPe'er 2007, $Y{=}1$",
}

# Not covered by grb_constants -- matches photospheric_radius/pe_er_photosphere.py's
# own MARKER_EDGE_WIDTH exactly, for visual consistency across the paper's figures.
MARKER_EDGE_WIDTH = 1.4
T90_MARKER = "o"  # matches pe_er_photosphere.py's GRB080916C marker (T90_MARKERS[0])


def episode_order(label):
    """Sort key putting episodes in temporal order: T90, EX0, TR1..TRn, EX1."""
    fixed = {"T90": 0, "EX0": 1, "EX1": 90}
    if label in fixed:
        return fixed[label]
    if label.startswith("TR"):
        return 10 + int(label[2:])
    return 99


def is_bb_model(model_name):
    """True if the named model includes a blackbody component."""
    return "BB" in model_name.upper()


def collect_intervals():
    """Real TimeInterval objects per episode of GRB080916C, for marker/order only."""
    _, _, grb_objects, _ = prepare_grbs(grb_list=[GRB_SHORT], result_file=ROOT / "results.json", get_best=True)
    grb = grb_objects[0]
    return {episode_label(m.interval): m.interval for m in grb.get_all_best_models()}


def load_data():
    """Read the three already-computed, self-contained CSVs, restricted to GRB080916C."""
    limit_a = pd.read_csv(LORENTZ_DIR / "lorentz_results.csv")
    limit_b = pd.read_csv(LORENTZ_DIR / "lorentz_results_limit_b.csv")
    thermal = pd.read_csv(PHOTOSPHERIC_CSV)

    limit_a = limit_a[limit_a.GRB == GRB_NAME]
    limit_b = limit_b[limit_b.GRB == GRB_NAME]
    thermal = thermal[(thermal.grb_name == GRB_NAME) & (thermal.Y_ratio == 1.0)]

    return limit_a, limit_b, thermal


def make_plot(limit_a, limit_b, thermal, path_stem="gamma_comparison"):
    """Episode on the x-axis, Gamma (log-scaled) on the y-axis, one series per method."""
    update_style()

    figure, axis = plt.subplots(figsize=(10, 5.5))

    intervals = collect_intervals()
    resolver = EpisodeMarkerResolver(t90_marker=T90_MARKER)
    episodes = sorted(limit_a["episode"].unique(), key=episode_order)
    x_positions = {episode: i for i, episode in enumerate(episodes)}

    for key, frame, model_col, y_col in (
        ("limit_a", limit_a, "model", "Gamma_min"),
        ("limit_b", limit_b, "model", "Gamma_min_B"),
        ("thermal", thermal, "model_name", "Gamma"),
    ):
        color = METHOD_COLORS[key]
        for _, row in frame.iterrows():
            episode = row["episode"]
            hollow = is_bb_model(row[model_col])
            median, lower, upper = row[y_col], row[f"{y_col}_err_lower"], row[f"{y_col}_err_upper"]

            axis.errorbar(
                x_positions[episode], median, yerr=[[lower], [upper]],
                marker=resolver.resolve(intervals[episode]), markersize=MARKER_SIZE * 1.4,
                markerfacecolor="none" if hollow else color, markeredgecolor=color,
                markeredgewidth=MARKER_EDGE_WIDTH, color=color, linestyle="none",
                capsize=CAP_SIZE, elinewidth=LINE_WIDTH, zorder=3,
            )

    axis.set_yscale("log")
    axis.margins(y=0.15)  # headroom so no marker (esp. the high thermal-Gamma points) touches the frame
    axis.set_xticks(range(len(episodes)))
    axis.xaxis.minorticks_off()
    axis.set_xticklabels(episodes, fontsize=LABEL_FONT_SIZE)
    axis.set_xlabel("Episode", fontsize=LABEL_FONT_SIZE)
    axis.set_ylabel(r"$\Gamma$", fontsize=LABEL_FONT_SIZE)
    axis.set_title(rf"{GRB_NAME}: $\gamma\gamma$-opacity limits vs. Thermal $\Gamma$", fontsize=TITLE_FONT_SIZE)

    # Shrink the axes (not the whole canvas) to leave room on the right for the two out-of-axes legends --
    # the Method legend's labels are multi-line/long, so the legends need real width, but the figure itself doesn't
    # need to keep growing.
    figure.subplots_adjust(right=0.62)

    method_handles = [
        mlines.Line2D(
            [], [], linestyle="none", marker="o", color=METHOD_COLORS[k],
            markerfacecolor=METHOD_COLORS[k], markersize=MARKER_SIZE * 1.4, label=METHOD_LABELS[k],
        )
        for k in ("limit_a", "limit_b", "thermal")
    ]

    # Marker legend: episode + best-fit model, keyed off Limit A's rows since every
    # episode with LAT coverage appears there. Legend entries identify episode *and*
    # model, per CLAUDE.md's convention, never a bare episode label.
    episode_models = limit_a.set_index("episode")["model"].to_dict()
    # print(f'{episode_models=}')
    episode_handles = [
        mlines.Line2D(
            [], [], linestyle="none", marker=resolver.resolve(intervals[episode]), color="#444444",
            markerfacecolor="none" if is_bb_model(episode_models[episode]) else "#444444",
            markeredgecolor="#444444", markeredgewidth=MARKER_EDGE_WIDTH, markersize=MARKER_SIZE * 1.4,
            label=f"{episode}" + r"$_\text{" + f'{episode_models[episode].replace("_", "+")}' + r"}$",
        )
        for episode in episodes
    ]

    # Both legends sit outside the axes, not "upper left"/"upper right" inside it --
    # an in-axes legend previously hid real data in this project (BUGS.md BUG-15,
    # pe_er_photosphere.py), and the same collision showed up here in review: the
    # thermal-Gamma points are the highest values on this log axis, so an in-axes
    # "upper left" box sat directly on top of them.
    legend1 = axis.legend(
        handles=method_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        fontsize=LEGEND_FONT_SIZE, frameon=True,
    )
    axis.add_artist(legend1)
    legend2 = axis.legend(
        handles=episode_handles, bbox_to_anchor=(1.06, 0.55),
        fontsize=LEGEND_FONT_SIZE, title=f"BEST model", title_fontsize=LEGEND_FONT_SIZE,
        frameon=True, ncols=1,
    )

    # No tight_layout() here: it conflicts with the tight savefig bbox (both fight
    # over how much room the out-of-axes legends get).
    #
    # update_style() already sets rcParams["savefig.bbox"] = "tight", but that alone
    # still clipped legend1 during review: a legend added via axis.add_artist() (as
    # opposed to the axes' own current legend) is not reliably picked up by
    # matplotlib's automatic tight-bbox artist search. Passing bbox_extra_artists
    # explicitly is the documented fix -- verified below by re-inspecting the actual
    # saved pixels, not just re-running the script.
    for extension in ("png", "pdf"):
        plt.savefig(
            f"./{path_stem}.{extension}", dpi=SAVE_DPI, bbox_inches="tight",
            bbox_extra_artists=(legend1, legend2),
        )
    plt.close()
    print(f"Saved: {path_stem}.png / .pdf")


def main():
    limit_a, limit_b, thermal = load_data()
    make_plot(limit_a, limit_b, thermal)


if __name__ == "__main__":
    main()
