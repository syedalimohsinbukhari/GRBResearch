"""
Three-Way Refit Comparison for GRB131014A (Joint / GBM-only / GBM-only, NaI >= 40 keV)
========================================================================================

Extends the two-way GBM-only robustness check (`gbm_only_refit.py`, weakness #6 in
`review-resolution.md`) with a third fit variant that answers weakness #7
("No instrumental-artifact / low-energy threshold robustness test"): the user re-ran the
GBM-only RMFIT fit for all five GRB131014A episodes with the NaI lower-energy bound raised
from the standard threshold to 40 keV (NaI band 40-900 keV), specifically to check whether
NaI response systematics at this burst's extreme geometry (ROI = 12 deg, zenith = 100 deg)
could be contaminating the low-energy end of the spectrum and manufacturing a spurious
blackbody. Added to results.json as "GRB131014215GBM40keV", alongside the existing
"GRB131014215" (joint GBM+LAT) and "GRB131014215GBM" (GBM-only, standard NaI range) keys.

This script compares all three fits per episode, deliberately as a *separate* output from
`gbm_only_refit.py` -- it does not overwrite that script's CSV or figures, since those are
already embedded in the paper (section-6-discussion.tex, weakness #6 write-up) describing a
two-way comparison; changing their content in place would silently invalidate text that
already quotes specific two-way numbers. Folding a third-fit narrative into the paper is a
separate, later decision (see gbm_only_refit.md Sec 6 for the reasoning and current status).

Shared helper functions (episode labeling, the f_BB Monte Carlo pipeline, the
Delta-C-stat threshold check, the paired-ratio fractional-difference calculation) are
imported from `gbm_only_refit.py` rather than re-derived, per this project's own precedent
of factoring out identical math shared between adjacent analyses instead of duplicating it
(compute_tau_hat, shared between lorentz_factor.py's Limit A and lorentz_factor_limit_b.py).

Method
------
Identical to gbm_only_refit.py's, applied pairwise across three fits instead of two:
- Model selection status ("_status": BEST) read directly from results.json.
- Independent BASE -> BASE+BB Delta-C-stat >= 28.74 check for the GBM-only and GBM-only
  (NaI >= 40 keV) fits (the joint fit's own status is established elsewhere in the paper's
  pipeline and not re-litigated here).
- f_BB (observer-frame, bolometric, 1 keV - 10 MeV) recomputed for all three fits with the
  same Monte Carlo machinery as bb_fraction/bb_flux_fraction.py.
- Fractional differences computed for two *adjacent* comparisons rather than all three
  pairwise combinations (see gbm_only_refit.md Sec 6 for why): GBM-only vs Joint (isolates
  the effect of dropping LAT -- already answered by gbm_only_refit.py, repeated here so this
  CSV is self-contained) and GBM-only(NaI>=40keV) vs GBM-only (isolates the effect of raising
  the NaI threshold, holding LAT-exclusion fixed -- the new weakness-#7 question). The third
  edge (GBM40 vs Joint, the combined effect) is not separately computed; it is recoverable by
  composing the two adjacent ratios.

Outputs
-------
    gbm_only_refit_3way_comparison.csv              -- one row per shared episode
    gbm_only_refit_3way_kt_comparison.png / .pdf     -- kT_BB, all three fits (primary)
    gbm_only_refit_3way_delta_cstat.png / .pdf       -- BB significance, all three fits (primary)
    gbm_only_refit_3way_fractional_diff.png / .pdf   -- % difference, two adjacent comparisons (secondary)
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gbm_only_refit as two_way
from grb_research import (
    EpisodeMarkerResolver,
    find_project_root,
    get_rng,
    seed_from_name,
    update_style,
    LEGEND_TITLE_FONT_SIZE,
)
from grb_research.grb_constants import (
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LINE_WIDTH,
    MARKER_SIZE,
    TICK_FONT_SIZE,
    CAP_SIZE,
    N_SAMPLES,
)
from grb_research.grb_core import GRB
from grb_research.grb_utils import save_fig

# ─── Configuration ───────────────────────────────────────────────────────────

JOINT_KEY = two_way.JOINT_KEY
GBM_KEY = two_way.GBM_KEY
GBM40_KEY = "GRB131014215GBM40keV"
GRB_DISPLAY_NAME = two_way.GRB_DISPLAY_NAME

T90_MARKER = two_way.T90_MARKER
EPISODE_ORDER = two_way.EPISODE_ORDER
DELTA_CSTAT_THRESHOLD = two_way.DELTA_CSTAT_THRESHOLD

SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)
PERCENTILES = (16.0, 50.0, 84.0)

# Fit-type colors, extending gbm_only_refit.py's red/black with a third, distinguishable hue.
# Blue reads clearly against the white/grid background update_style() sets and is not used
# elsewhere in this folder's figures, so it is unambiguous as "the NaI>=40keV fit" throughout.
JOINT_COLOR = two_way.JOINT_COLOR
GBM_COLOR = two_way.GBM_COLOR
GBM40_COLOR = "blue"

FIT_LABELS = {
    "joint": "Joint (GBM+LAT)",
    "gbm": "GBM-only",
    "gbm40": "GBM-only (NaI 40–900 keV)",
}

MARKER_EDGE_WIDTH = two_way.MARKER_EDGE_WIDTH

# Shorthand for per-point model identification on x-tick labels, e.g. BAND_BB -> B+BB.
MODEL_SHORTHAND = {"BAND_BB": "B+BB", "CPL_BB": "C+BB", "SBPL_BB": "S+BB", "PL_BB": "P+BB"}


def _short(model_name: str) -> str:
    return MODEL_SHORTHAND.get(model_name, model_name.replace("_", "+"))


# ─── Row assembly ────────────────────────────────────────────────────────────


@dataclass
class ComparisonRow3:
    grb_name: str
    episode: str
    t_start_s: float
    t_stop_s: float
    model_name_joint: str
    model_name_gbm: str
    model_name_gbm40: str
    models_agree_joint_gbm: bool
    models_agree_joint_gbm40: bool
    models_agree_gbm_gbm40: bool
    models_agree_all3: bool
    kt_bb_keV_joint: float
    kt_bb_err_keV_joint: float
    kt_bb_keV_gbm: float
    kt_bb_err_keV_gbm: float
    kt_bb_keV_gbm40: float
    kt_bb_err_keV_gbm40: float
    kt_consistent_1sigma_joint_gbm: object
    kt_consistent_1sigma_joint_gbm40: object
    kt_consistent_1sigma_gbm_gbm40: object
    f_bb_joint: float
    f_bb_err_lower_joint: float
    f_bb_err_upper_joint: float
    f_bb_gbm: float
    f_bb_err_lower_gbm: float
    f_bb_err_upper_gbm: float
    f_bb_gbm40: float
    f_bb_err_lower_gbm40: float
    f_bb_err_upper_gbm40: float
    cstat_joint: float
    dof_joint: int
    cstat_gbm: float
    dof_gbm: int
    cstat_gbm40: float
    dof_gbm40: int
    delta_cstat_joint: object
    delta_cstat_gbm: object
    delta_cstat_gbm40: object
    delta_cstat_threshold: float
    bb_required_joint: object
    bb_required_gbm: object
    bb_required_gbm40: object
    any_bb_step_clears_threshold_gbm: object
    recorded_winner_is_bb_gbm: object
    status_consistent_with_threshold_gbm: object
    any_bb_step_clears_threshold_gbm40: object
    recorded_winner_is_bb_gbm40: object
    status_consistent_with_threshold_gbm40: object
    kt_bb_diff_pct_gbm_vs_joint: float
    kt_bb_diff_pct_err_lower_gbm_vs_joint: float
    kt_bb_diff_pct_err_upper_gbm_vs_joint: float
    f_bb_diff_pct_gbm_vs_joint: float
    f_bb_diff_pct_err_lower_gbm_vs_joint: float
    f_bb_diff_pct_err_upper_gbm_vs_joint: float
    kt_bb_diff_pct_gbm40_vs_gbm: float
    kt_bb_diff_pct_err_lower_gbm40_vs_gbm: float
    kt_bb_diff_pct_err_upper_gbm40_vs_gbm: float
    f_bb_diff_pct_gbm40_vs_gbm: float
    f_bb_diff_pct_err_lower_gbm40_vs_gbm: float
    f_bb_diff_pct_err_upper_gbm40_vs_gbm: float
    n_samples: int
    seed: int


def _bb_fields(model, rng):
    if "BB" not in model.name:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    kt = model.get_parameter_value("kt_bb")
    kt_err = next(p.error for p in model.parameters if p.name == "kt_bb")
    f_med, f_lo, f_hi = two_way.compute_f_bb(model, rng=rng)
    return (kt, kt_err, f_med, f_lo, f_hi)


def build_row(ep, joint_interval, gbm_interval, gbm40_interval, rng):
    joint_best = joint_interval.models.best
    gbm_best = gbm_interval.models.best
    gbm40_best = gbm40_interval.models.best

    models_agree_joint_gbm = joint_best.name == gbm_best.name
    models_agree_joint_gbm40 = joint_best.name == gbm40_best.name
    models_agree_gbm_gbm40 = gbm_best.name == gbm40_best.name
    models_agree_all3 = models_agree_joint_gbm and models_agree_joint_gbm40

    kt_j, kt_j_err, f_j, f_j_lo, f_j_hi = _bb_fields(joint_best, rng)
    kt_g, kt_g_err, f_g, f_g_lo, f_g_hi = _bb_fields(gbm_best, rng)
    kt_g40, kt_g40_err, f_g40, f_g40_lo, f_g40_hi = _bb_fields(gbm40_best, rng)

    def _consistent(a, a_err, b, b_err):
        if not (np.isfinite(a) and np.isfinite(b)):
            return np.nan
        return bool(abs(a - b) <= np.hypot(a_err, b_err))

    kt_consistent_joint_gbm = _consistent(kt_j, kt_j_err, kt_g, kt_g_err)
    kt_consistent_joint_gbm40 = _consistent(kt_j, kt_j_err, kt_g40, kt_g40_err)
    kt_consistent_gbm_gbm40 = _consistent(kt_g, kt_g_err, kt_g40, kt_g40_err)

    delta_joint = two_way.base_bb_delta_cstat(joint_interval, joint_best.name) if "BB" in joint_best.name else None
    delta_gbm = two_way.base_bb_delta_cstat(gbm_interval, gbm_best.name) if "BB" in gbm_best.name else None
    delta_gbm40 = (
        two_way.base_bb_delta_cstat(gbm40_interval, gbm40_best.name) if "BB" in gbm40_best.name else None
    )

    bb_required_joint = None if delta_joint is None else bool(delta_joint >= DELTA_CSTAT_THRESHOLD)
    bb_required_gbm = None if delta_gbm is None else bool(delta_gbm >= DELTA_CSTAT_THRESHOLD)
    bb_required_gbm40 = None if delta_gbm40 is None else bool(delta_gbm40 >= DELTA_CSTAT_THRESHOLD)

    any_bb_clears_gbm = two_way.any_bb_step_clears_threshold(gbm_interval)
    recorded_winner_is_bb_gbm = "BB" in gbm_best.name
    status_consistent_gbm = any_bb_clears_gbm == recorded_winner_is_bb_gbm

    any_bb_clears_gbm40 = two_way.any_bb_step_clears_threshold(gbm40_interval)
    recorded_winner_is_bb_gbm40 = "BB" in gbm40_best.name
    status_consistent_gbm40 = any_bb_clears_gbm40 == recorded_winner_is_bb_gbm40

    # Adjacent-comparison fractional differences (see module docstring for why these two edges
    # and not all three pairwise combinations): Joint -> GBM-only isolates dropping LAT;
    # GBM-only -> GBM-only(NaI>=40keV) isolates raising the NaI threshold.
    (
        kt_diff_gj, kt_diff_gj_lo, kt_diff_gj_hi,
        f_diff_gj, f_diff_gj_lo, f_diff_gj_hi,
    ) = two_way.paired_ratio_pct(joint_best, gbm_best, rng=rng)
    (
        kt_diff_g40g, kt_diff_g40g_lo, kt_diff_g40g_hi,
        f_diff_g40g, f_diff_g40g_lo, f_diff_g40g_hi,
    ) = two_way.paired_ratio_pct(gbm_best, gbm40_best, rng=rng)

    return ComparisonRow3(
        grb_name=GRB_DISPLAY_NAME,
        episode=ep,
        t_start_s=joint_interval.start,
        t_stop_s=joint_interval.end,
        model_name_joint=joint_best.name,
        model_name_gbm=gbm_best.name,
        model_name_gbm40=gbm40_best.name,
        models_agree_joint_gbm=models_agree_joint_gbm,
        models_agree_joint_gbm40=models_agree_joint_gbm40,
        models_agree_gbm_gbm40=models_agree_gbm_gbm40,
        models_agree_all3=models_agree_all3,
        kt_bb_keV_joint=kt_j,
        kt_bb_err_keV_joint=kt_j_err,
        kt_bb_keV_gbm=kt_g,
        kt_bb_err_keV_gbm=kt_g_err,
        kt_bb_keV_gbm40=kt_g40,
        kt_bb_err_keV_gbm40=kt_g40_err,
        kt_consistent_1sigma_joint_gbm=kt_consistent_joint_gbm,
        kt_consistent_1sigma_joint_gbm40=kt_consistent_joint_gbm40,
        kt_consistent_1sigma_gbm_gbm40=kt_consistent_gbm_gbm40,
        f_bb_joint=f_j,
        f_bb_err_lower_joint=f_j_lo,
        f_bb_err_upper_joint=f_j_hi,
        f_bb_gbm=f_g,
        f_bb_err_lower_gbm=f_g_lo,
        f_bb_err_upper_gbm=f_g_hi,
        f_bb_gbm40=f_g40,
        f_bb_err_lower_gbm40=f_g40_lo,
        f_bb_err_upper_gbm40=f_g40_hi,
        cstat_joint=joint_best.cstat,
        dof_joint=joint_best.dof,
        cstat_gbm=gbm_best.cstat,
        dof_gbm=gbm_best.dof,
        cstat_gbm40=gbm40_best.cstat,
        dof_gbm40=gbm40_best.dof,
        delta_cstat_joint=delta_joint,
        delta_cstat_gbm=delta_gbm,
        delta_cstat_gbm40=delta_gbm40,
        delta_cstat_threshold=DELTA_CSTAT_THRESHOLD,
        bb_required_joint=bb_required_joint,
        bb_required_gbm=bb_required_gbm,
        bb_required_gbm40=bb_required_gbm40,
        any_bb_step_clears_threshold_gbm=any_bb_clears_gbm,
        recorded_winner_is_bb_gbm=recorded_winner_is_bb_gbm,
        status_consistent_with_threshold_gbm=status_consistent_gbm,
        any_bb_step_clears_threshold_gbm40=any_bb_clears_gbm40,
        recorded_winner_is_bb_gbm40=recorded_winner_is_bb_gbm40,
        status_consistent_with_threshold_gbm40=status_consistent_gbm40,
        kt_bb_diff_pct_gbm_vs_joint=kt_diff_gj,
        kt_bb_diff_pct_err_lower_gbm_vs_joint=kt_diff_gj_lo,
        kt_bb_diff_pct_err_upper_gbm_vs_joint=kt_diff_gj_hi,
        f_bb_diff_pct_gbm_vs_joint=f_diff_gj,
        f_bb_diff_pct_err_lower_gbm_vs_joint=f_diff_gj_lo,
        f_bb_diff_pct_err_upper_gbm_vs_joint=f_diff_gj_hi,
        kt_bb_diff_pct_gbm40_vs_gbm=kt_diff_g40g,
        kt_bb_diff_pct_err_lower_gbm40_vs_gbm=kt_diff_g40g_lo,
        kt_bb_diff_pct_err_upper_gbm40_vs_gbm=kt_diff_g40g_hi,
        f_bb_diff_pct_gbm40_vs_gbm=f_diff_g40g,
        f_bb_diff_pct_err_lower_gbm40_vs_gbm=f_diff_g40g_lo,
        f_bb_diff_pct_err_upper_gbm40_vs_gbm=f_diff_g40g_hi,
        n_samples=N_SAMPLES,
        seed=SEED,
    )


# ─── Plots ───────────────────────────────────────────────────────────────────

WIDTH = 0.22


def _tick_label(row):
    """Episode tick label: single model line if all three fits agree, else one line per fit.

    Extends gbm_only_refit.py's convention of carrying the model name on the x-tick instead of
    in the legend (still true here: color/legend encodes fit type, not model) -- but that
    convention assumed models_agree was True throughout, which no longer holds with a third fit
    in the mix (GRB131014A EX1: BAND_BB in Joint/GBM-only, SBPL_BB in the NaI>=40keV refit).
    When the three fits disagree, all three model names are shown rather than picking one, so
    the figure never silently hides a model-selection disagreement it would otherwise show only
    in the CSV.
    """
    if row["models_agree_all3"]:
        return f"{row['episode']}\n({_short(row['model_name_joint'])})"
    return (
        f"{row['episode']}\n"
        f"J:{_short(row['model_name_joint'])} G:{_short(row['model_name_gbm'])} "
        f"G40:{_short(row['model_name_gbm40'])}"
    )


def make_kt_plot(df, episode_markers, path_stub="gbm_only_refit_3way_kt_comparison"):
    update_style()
    plotted = df.dropna(subset=["kt_bb_keV_joint"])
    fig, ax = plt.subplots()
    x = np.arange(len(plotted))

    series = [
        ("kt_bb_keV_joint", "kt_bb_err_keV_joint", -WIDTH, JOINT_COLOR),
        ("kt_bb_keV_gbm", "kt_bb_err_keV_gbm", 0.0, GBM_COLOR),
        ("kt_bb_keV_gbm40", "kt_bb_err_keV_gbm40", WIDTH, GBM40_COLOR),
    ]
    for i, (_, row) in enumerate(plotted.iterrows()):
        marker = episode_markers[row["episode"]]
        for val_col, err_col, offset, color in series:
            if not np.isfinite(row[val_col]):
                continue
            ax.errorbar(
                x[i] + offset,
                row[val_col],
                yerr=row[err_col],
                fmt=marker,
                mfc="none",
                mec=color,
                ecolor=color,
                capsize=CAP_SIZE,
                markersize=MARKER_SIZE,
                markeredgewidth=MARKER_EDGE_WIDTH,
                linewidth=LINE_WIDTH,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([_tick_label(row) for _, row in plotted.iterrows()], fontsize=TICK_FONT_SIZE)
    ax.set_ylabel(r"$kT_\mathrm{BB}$ [keV]", fontsize=LABEL_FONT_SIZE)
    ax.tick_params(labelsize=TICK_FONT_SIZE)

    fit_handles = [
        mlines.Line2D(
            [], [], marker="o", mfc="none", mec=color, color=color, linestyle="None",
            markeredgewidth=MARKER_EDGE_WIDTH, markersize=MARKER_SIZE, label=FIT_LABELS[key],
        )
        for key, color in (("joint", JOINT_COLOR), ("gbm", GBM_COLOR), ("gbm40", GBM40_COLOR))
    ]
    fit_legend = ax.legend(
        handles=fit_handles, title="Fit type", fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_FONT_SIZE, loc="best",
    )
    ax.add_artist(fit_legend)
    save_fig(fig, path_stub)


def make_delta_cstat_plot(df, episode_markers, path_stub="gbm_only_refit_3way_delta_cstat"):
    update_style()
    plotted = df.dropna(subset=["delta_cstat_joint", "delta_cstat_gbm", "delta_cstat_gbm40"])
    fig, ax = plt.subplots()
    x = np.arange(len(plotted))

    ax.axhline(DELTA_CSTAT_THRESHOLD, color="black", linewidth=LINE_WIDTH, linestyle="--")

    series = [
        ("delta_cstat_joint", -WIDTH, JOINT_COLOR),
        ("delta_cstat_gbm", 0.0, GBM_COLOR),
        ("delta_cstat_gbm40", WIDTH, GBM40_COLOR),
    ]
    for i, (_, row) in enumerate(plotted.iterrows()):
        marker = episode_markers[row["episode"]]
        for val_col, offset, color in series:
            ax.plot(
                x[i] + offset, row[val_col], marker=marker, mfc="none", mec=color, color=color,
                mew=MARKER_EDGE_WIDTH, markersize=MARKER_SIZE, linestyle="None",
            )

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([_tick_label(row) for _, row in plotted.iterrows()], fontsize=TICK_FONT_SIZE)
    ax.set_ylabel(r"$\Delta C\mathrm{-stat}$ (BASE $\to$ BASE+BB)", fontsize=LABEL_FONT_SIZE)
    ax.tick_params(labelsize=TICK_FONT_SIZE)

    fit_handles = [
        mlines.Line2D(
            [], [], marker="o", mfc="none", mec=color, mew=MARKER_EDGE_WIDTH, color=color,
            linestyle="None", markersize=MARKER_SIZE, label=FIT_LABELS[key],
        )
        for key, color in (("joint", JOINT_COLOR), ("gbm", GBM_COLOR), ("gbm40", GBM40_COLOR))
    ]
    fit_handles.append(
        mlines.Line2D([], [], color="black", linewidth=LINE_WIDTH, linestyle="--",
                       label=f"threshold = {DELTA_CSTAT_THRESHOLD}")
    )
    fit_legend = ax.legend(
        handles=fit_handles, title="Fit type", fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_FONT_SIZE, loc="best",
    )
    ax.add_artist(fit_legend)
    save_fig(fig, path_stub)


def make_fractional_diff_plot(df, episode_markers, path_stub="gbm_only_refit_3way_fractional_diff"):
    """Two adjacent-comparison fractional differences, stacked (kT_BB, f_BB) x 2 columns.

    Left column: GBM-only vs Joint (dropping LAT) -- reproduces gbm_only_refit.py's own
    fractional-diff result for cross-reference. Right column: GBM-only(NaI>=40keV) vs GBM-only
    (raising the NaI threshold) -- the new weakness-#7 question this script exists to answer.
    Same two-panel-per-comparison layout as gbm_only_refit.py, doubled into a 2x2 grid rather
    than overlaying both comparisons on one set of axes, since the two ratios have different
    reference fits and overlaying them on a shared x-offset would invite misreading one
    comparison's error bar as commensurate with the other's.
    """
    update_style()
    plotted = df.dropna(subset=["kt_bb_keV_joint", "kt_bb_keV_gbm", "kt_bb_keV_gbm40"])
    fig, axes = plt.subplots(2, 2, sharex=True, figsize=(9, 6))
    x = np.arange(len(plotted))

    columns = [
        (
            0, GBM_COLOR, "GBM-only vs Joint",
            "kt_bb_diff_pct_gbm_vs_joint", "kt_bb_diff_pct_err_lower_gbm_vs_joint",
            "kt_bb_diff_pct_err_upper_gbm_vs_joint",
            "f_bb_diff_pct_gbm_vs_joint", "f_bb_diff_pct_err_lower_gbm_vs_joint",
            "f_bb_diff_pct_err_upper_gbm_vs_joint",
        ),
        (
            1, GBM40_COLOR, "GBM-only (NaI 40–900 keV) vs GBM-only",
            "kt_bb_diff_pct_gbm40_vs_gbm", "kt_bb_diff_pct_err_lower_gbm40_vs_gbm",
            "kt_bb_diff_pct_err_upper_gbm40_vs_gbm",
            "f_bb_diff_pct_gbm40_vs_gbm", "f_bb_diff_pct_err_lower_gbm40_vs_gbm",
            "f_bb_diff_pct_err_upper_gbm40_vs_gbm",
        ),
    ]

    for col_idx, color, title, kt_col, kt_lo, kt_hi, f_col, f_lo, f_hi in columns:
        for row_idx, (val_col, lo_col, hi_col, ylabel) in enumerate(
            [(kt_col, kt_lo, kt_hi, r"$\Delta kT_\mathrm{BB}$ [%]"), (f_col, f_lo, f_hi, r"$\Delta f_\mathrm{BB}$ [%]")]
        ):
            ax = axes[row_idx, col_idx]
            ax.axhline(0.0, color="black", linewidth=LINE_WIDTH, linestyle="--")
            for i, (_, row) in enumerate(plotted.iterrows()):
                marker = episode_markers[row["episode"]]
                ax.errorbar(
                    x[i], row[val_col], yerr=[[row[lo_col]], [row[hi_col]]], fmt=marker,
                    mfc="none", mec=color, ecolor=color, capsize=CAP_SIZE, markersize=MARKER_SIZE,
                    markeredgewidth=MARKER_EDGE_WIDTH, linewidth=LINE_WIDTH,
                )
            if row_idx == 0:
                ax.set_title(title, fontsize=LABEL_FONT_SIZE)
            if col_idx == 0:
                ax.set_ylabel(ylabel, fontsize=LABEL_FONT_SIZE)
            ax.tick_params(labelsize=TICK_FONT_SIZE)

    for ax in axes[1, :]:
        ax.set_xticks(x)
        ax.set_xticklabels(plotted["episode"], fontsize=TICK_FONT_SIZE)

    fig.tight_layout()
    save_fig(fig, path_stub)


# ─── Main ────────────────────────────────────────────────────────────────────


def main():
    root = find_project_root()
    with open(root / "results.json") as f:
        data = json.load(f)

    joint = GRB.from_dictionary(name=JOINT_KEY, grb_data=data[JOINT_KEY])
    gbm = GRB.from_dictionary(name=GBM_KEY, grb_data=data[GBM_KEY])
    gbm40 = GRB.from_dictionary(name=GBM40_KEY, grb_data=data[GBM40_KEY])

    joint_eps = {two_way.episode_label(i): i for i in joint.intervals}
    gbm_eps = {two_way.episode_label(i): i for i in gbm.intervals}
    gbm40_eps = {two_way.episode_label(i): i for i in gbm40.intervals}

    rows = []
    for ep in EPISODE_ORDER:
        if ep not in joint_eps or ep not in gbm_eps or ep not in gbm40_eps:
            continue
        rows.append(build_row(ep, joint_eps[ep], gbm_eps[ep], gbm40_eps[ep], rng))

    df = pd.DataFrame([r.__dict__ for r in rows])
    df.to_csv("gbm_only_refit_3way_comparison.csv", index=False)

    resolver = EpisodeMarkerResolver(t90_marker=T90_MARKER)
    episode_markers = {ep: resolver.resolve(joint_eps[ep]) for ep in EPISODE_ORDER if ep in joint_eps}

    make_kt_plot(df, episode_markers)
    make_delta_cstat_plot(df, episode_markers)
    make_fractional_diff_plot(df, episode_markers)

    pd.set_option("display.width", 220)
    print(df[[
        "episode", "model_name_joint", "model_name_gbm", "model_name_gbm40", "models_agree_all3",
        "kt_bb_keV_joint", "kt_bb_keV_gbm", "kt_bb_keV_gbm40",
        "kt_consistent_1sigma_joint_gbm", "kt_consistent_1sigma_joint_gbm40", "kt_consistent_1sigma_gbm_gbm40",
        "status_consistent_with_threshold_gbm", "status_consistent_with_threshold_gbm40",
    ]])
    print(df[[
        "episode", "kt_bb_diff_pct_gbm_vs_joint", "kt_bb_diff_pct_gbm40_vs_gbm",
        "f_bb_diff_pct_gbm_vs_joint", "f_bb_diff_pct_gbm40_vs_gbm",
    ]])


if __name__ == "__main__":
    main()
