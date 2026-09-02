"""
GBM-only Refit Robustness Check for GRB131014A
===============================================

Weakness #6 (review-resolution.md) asks whether the blackbody (BB) detections
in GRB131014A's five episodes could be an artefact of including Fermi-LAT data
in the joint GBM+LAT fit, given the burst's large (~70 deg) LAT off-axis angle
for its entire T90. A joint fit cannot rule this out on its own: if LAT
systematics biased the continuum parameters through the shared likelihood,
the same joint fit would not reveal it.

The direct test is to drop LAT entirely and refit GBM (NaI+BGO) data alone.
The user ran this refit in RMFIT and added the result to results.json under
the key "GRB131014215GBM", alongside the existing joint-fit key
"GRB131014215". This script compares the two: for each of the five shared
episodes, does the same model win, does BB survive on GBM data alone, and are
the fitted kT_BB and f_BB consistent between the two fits.

Method
------
- Model selection status ("_status": BEST/SAFE/UNSAFE) is read directly from
  results.json, as recorded by the RMFIT-based pipeline used throughout this
  project -- this script does not re-derive it.
- As an independent check, every BASE -> BASE+BB step available in an
  episode's fit output is re-verified from the raw C-stat values, using the
  same Delta-C-stat >= 28.74 threshold (Delta k = 2) used throughout the
  paper (section-1-introduction.tex): if any such step clears the threshold,
  a BB-augmented model should be the recorded winner. This catches any
  disagreement between the recorded "_status" and that rule -- see
  gbm_only_refit.md Sec 2 for why this matters for the T90 episode.
- f_BB (observer-frame, bolometric) is recomputed for the GBM-only fit with
  the same Monte Carlo machinery as codes-for-paper/bb_fraction, so the two
  columns are on identical footing. The energy band, MC settings and
  computation (`component_energy_fluxes` + `draw_model_samples`, then
  percentiles of the accepted draws) are copied from
  codes-for-paper/bb_fraction/bb_flux_fraction.py, since f_BB is not
  reconstructable from the fit files without them (cross-folder-import
  convention, CLAUDE.md).

Outputs
-------
    gbm_only_refit_comparison.csv              -- one row per shared episode
    gbm_only_refit_kt_comparison.png / .pdf    -- kT_BB, joint vs GBM-only (primary)
    gbm_only_refit_delta_cstat.png / .pdf      -- BB significance, joint vs GBM-only (primary)
    gbm_only_refit_fractional_diff.png / .pdf  -- % difference (GBM-only vs Joint) in kT_BB and f_BB (secondary)
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import (
    EpisodeMarkerResolver,
    component_energy_fluxes,
    draw_model_samples,
    find_project_root,
    update_style, LEGEND_TITLE_FONT_SIZE,
)
from grb_research.grb_constants import (
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LINE_WIDTH,
    MARKER_SIZE,
    TICK_FONT_SIZE, CAP_SIZE,
)
from grb_research.grb_core import GRB
from grb_research.grb_enums import GRBModelsCombinations as gmC
from grb_research.grb_utils import save_fig

# ─── Configuration ───────────────────────────────────────────────────────────

JOINT_KEY = "GRB131014215"
GBM_KEY = "GRB131014215GBM"
GRB_DISPLAY_NAME = "GRB131014A"

# T90 marker for GRB131014A, matching bb_flux_fraction.py's T90_MARKERS[1]
# (GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]), so this burst's
# T90 point looks the same here as in the rest of the paper's figures.
T90_MARKER = "s"

# Observer-frame bolometric band, matching bb_fraction/bb_flux_fraction.py --
# f_BB is redshift/cosmology-independent per-episode in this band (established
# there), so no z or cosmology column is carried in the output CSV.
E_MIN_KEV, E_MAX_KEV = 1.0, 1.0e4
N_GRID = 1_000
CHUNK_SIZE = 2_000

# MC settings -- the project-wide convention (amati_relationship.py:39-41).
N_SAMPLES = 10_000
SEED = 12345
PERCENTILES = (16.0, 50.0, 84.0)

# Delta-C-stat threshold for a BASE -> BASE+BB step (Delta k = 2), the same
# value used throughout the paper (section-1-introduction.tex).
DELTA_CSTAT_THRESHOLD = 28.74

# Each BB-augmented model's non-BB base model, for the BASE -> BASE+BB check.
BASE_OF = {"BAND_BB": "BAND", "CPL_BB": "CPL", "SBPL_BB": "SBPL", "PL_BB": "PL"}

# Episodes in time order: leading excess, the two time-resolved intervals, the trailing excess, with T90 first as the
# time-integrated reference.
EPISODE_ORDER = ["T90", "EX0", "TR1", "TR2", "EX1"]

# Fit-type color (per user request, 2026-09-01): yellow/black rather than GRBPlotStyle's per-GRB color, since both
# series are the *same* burst here -- the color dimension in this folder encodes fit type, not GRB identity.
# "gold" rather than pure #FFFF00: legible against the white/grid background update_style() sets, while still reading
# as yellow.
#
# changed to red for better contrast (Sep 1, 2026 - 23:33:20)
JOINT_COLOR = "red"
GBM_COLOR = "black"

MARKER_EDGE_WIDTH = 1.4


# ─── Helpers ─────────────────────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1. Copied from bb_flux_fraction.py."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def analytic_bb_bolometric_flux(amp_bb, kt_bb):
    """All-frequency blackbody energy flux [keV/cm^2/s]. Copied from bb_flux_fraction.py."""
    return amp_bb * kt_bb ** 4 * np.pi ** 4 / 15.0


def split_energy_flux(model_name, values, energy):
    """BB and total energy flux for one or many parameter draws [keV/cm^2/s]."""
    fluxes, total = component_energy_fluxes(model_name, values, energy, chunk=CHUNK_SIZE)
    return fluxes[gmC.BB], total


def compute_f_bb(model, n_samples=N_SAMPLES, seed=SEED):
    """Median and 1-sigma asymmetric interval of f_BB for a BB-augmented model."""
    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)
    samples = draw_model_samples(model, n_samples=n_samples, seed=seed)
    sample_bb, sample_total = split_energy_flux(model.name, samples, energy)
    fractions = sample_bb / sample_total
    usable = np.isfinite(fractions) & (fractions > 0) & (fractions < 1)
    lo, med, hi = np.percentile(fractions[usable], PERCENTILES)
    return med, med - lo, hi - med


def paired_ratio_pct(joint_model, gbm_model, n_samples=N_SAMPLES, seed=SEED):
    """Paired-MC fractional difference, (GBM-only - Joint) / Joint, in kT_BB and f_BB.

    Draws n_samples parameter vectors from *each* fit's own covariance matrix
    (draw_model_samples, which already applies that model's physical-validity
    resampling) and forms the ratio per draw, rather than linearizing two
    marginal errors together (the previous version's fractional_diff_pct).
    This picks up whatever asymmetry/non-Gaussianity each fit's own posterior
    has, at the cost of needing the same functional model on both sides (only
    called when joint_model.name == gbm_model.name, true for all 5 episodes
    here per the #6 result).

    The two draws use *different* seeds (seed and seed+1), not the same one.
    A same-seed (common-random-numbers) version was tried first and rejected
    after validation (see gbm_only_refit.md Sec 3): because both fits use the
    same model form with similar covariance structure, sharing a seed made
    the two draws' kt_bb nearly perfectly correlated (r = 0.9996, measured),
    collapsing the ratio's spread to an artificially tight ~0.15%-wide
    interval -- about 20x tighter than either an independent-seed MC run or
    the original linearized (delta-method) estimate, both of which agree with
    each other at the few-percent level. That tightness was a seed-sharing
    artifact, not a real reduction in uncertainty, so independent seeds are
    used instead; the two fits are separate RMFIT optimizations with no
    actual joint covariance available.

    Returns (kt_diff_pct, kt_diff_lo, kt_diff_hi, f_diff_pct, f_diff_lo, f_diff_hi):
    the paired ratio distribution's median and its (median - p16, p84 - median)
    asymmetric interval, for kT_BB and f_BB respectively.
    """
    if "BB" not in joint_model.name or "BB" not in gbm_model.name:
        return (np.nan,) * 6

    energy = np.logspace(np.log10(E_MIN_KEV), np.log10(E_MAX_KEV), N_GRID)

    joint_samples = draw_model_samples(joint_model, n_samples=n_samples, seed=seed)
    gbm_samples = draw_model_samples(gbm_model, n_samples=n_samples, seed=seed + 1)

    kt_idx_joint = [p.name for p in joint_model.parameters].index("kt_bb")
    kt_idx_gbm = [p.name for p in gbm_model.parameters].index("kt_bb")
    kt_diff = (gbm_samples[:, kt_idx_gbm] / joint_samples[:, kt_idx_joint] - 1.0) * 100.0
    kt_lo, kt_med, kt_hi = np.percentile(kt_diff, PERCENTILES)

    f_joint_bb, f_joint_total = split_energy_flux(joint_model.name, joint_samples, energy)
    f_gbm_bb, f_gbm_total = split_energy_flux(gbm_model.name, gbm_samples, energy)
    f_joint = f_joint_bb / f_joint_total
    f_gbm = f_gbm_bb / f_gbm_total

    usable = np.isfinite(f_joint) & (f_joint > 0) & (f_joint < 1) & np.isfinite(f_gbm) & (f_gbm > 0) & (f_gbm < 1)
    f_lo, f_med, f_hi = np.percentile((f_gbm[usable] / f_joint[usable] - 1.0) * 100.0, PERCENTILES)

    return kt_med, kt_med - kt_lo, kt_hi - kt_med, f_med, f_med - f_lo, f_hi - f_med


def base_bb_delta_cstat(interval, bb_model_name):
    """C-stat(base) - C-stat(base+BB) for one interval, or None if either is missing."""
    base_name = BASE_OF.get(bb_model_name)
    base_model = interval.models.get(base_name) if base_name else None
    bb_model = interval.models.get(bb_model_name)
    if base_model is None or bb_model is None:
        return None
    return base_model.cstat - bb_model.cstat


def any_bb_step_clears_threshold(interval):
    """Whether any fitted BASE -> BASE+BB pair in this interval clears DELTA_CSTAT_THRESHOLD.

    Independent of which model is recorded as "_status": BEST -- this checks
    every base/BB pair actually present in the fit output, so it flags a
    recorded non-BB winner even when the specific BB extension of *that* base
    model failed to fit (see gbm_only_refit.md Sec 2, the T90 case).
    """
    for bb_name, base_name in BASE_OF.items():
        delta = base_bb_delta_cstat(interval, bb_name)
        if delta is not None and delta >= DELTA_CSTAT_THRESHOLD:
            return True
    return False


@dataclass
class ComparisonRow:
    grb_name: str
    episode: str
    t_start_s: float
    t_stop_s: float
    model_name_joint: str
    model_name_gbm: str
    models_agree: bool
    kt_bb_keV_joint: float
    kt_bb_err_keV_joint: float
    kt_bb_keV_gbm: float
    kt_bb_err_keV_gbm: float
    kt_consistent_1sigma: object
    f_bb_joint: float
    f_bb_err_lower_joint: float
    f_bb_err_upper_joint: float
    f_bb_gbm: float
    f_bb_err_lower_gbm: float
    f_bb_err_upper_gbm: float
    cstat_joint: float
    dof_joint: int
    cstat_gbm: float
    dof_gbm: int
    delta_cstat_joint: object
    delta_cstat_gbm: object
    delta_cstat_threshold: float
    bb_required_joint: object
    bb_required_gbm: object
    any_bb_step_clears_threshold_gbm: object
    recorded_winner_is_bb_gbm: object
    status_consistent_with_threshold_gbm: object
    kt_bb_diff_pct: float
    kt_bb_diff_pct_err_lower: float
    kt_bb_diff_pct_err_upper: float
    f_bb_diff_pct: float
    f_bb_diff_pct_err_lower: float
    f_bb_diff_pct_err_upper: float
    n_samples: int
    seed: int


def build_row(ep, joint_interval, gbm_interval):
    joint_best = joint_interval.models.best
    gbm_best = gbm_interval.models.best

    models_agree = joint_best.name == gbm_best.name

    def bb_fields(model):
        if "BB" not in model.name:
            return (np.nan, np.nan, np.nan, np.nan, np.nan)
        kt = model.get_parameter_value("kt_bb")
        kt_err = next(p.error for p in model.parameters if p.name == "kt_bb")
        f_med, f_lo, f_hi = compute_f_bb(model)
        return (kt, kt_err, f_med, f_lo, f_hi)

    kt_j, kt_j_err, f_j, f_j_lo, f_j_hi = bb_fields(joint_best)
    kt_g, kt_g_err, f_g, f_g_lo, f_g_hi = bb_fields(gbm_best)

    kt_consistent = np.nan
    if np.isfinite(kt_j) and np.isfinite(kt_g):
        combined_err = np.hypot(kt_j_err, kt_g_err)
        kt_consistent = bool(abs(kt_j - kt_g) <= combined_err)

    delta_joint = base_bb_delta_cstat(joint_interval, joint_best.name) if "BB" in joint_best.name else None
    delta_gbm = base_bb_delta_cstat(gbm_interval, gbm_best.name) if "BB" in gbm_best.name else None

    bb_required_joint = None if delta_joint is None else bool(delta_joint >= DELTA_CSTAT_THRESHOLD)
    bb_required_gbm = None if delta_gbm is None else bool(delta_gbm >= DELTA_CSTAT_THRESHOLD)

    # Independent sanity check on the GBM-only fit:
    # If any fitted base/BB pair in this interval clears the threshold, the recorded winner should be a BB model
    # (some BB model -- not necessarily the specific pair that cleared it, since a different BB model may fit better
    # still).
    any_bb_clears_gbm = any_bb_step_clears_threshold(gbm_interval)
    recorded_winner_is_bb_gbm = "BB" in gbm_best.name
    status_consistent_gbm = any_bb_clears_gbm == recorded_winner_is_bb_gbm

    # Fractional difference, GBM-only relative to Joint
    # How much does dropping LAT actually move each parameter, not just whether it's consistent within error.
    #
    # Paired-MC (see paired_ratio_pct docstring):
    # Draws each fit's own parameter samples and takes percentiles of the ratio distribution directly, rather than
    # linearizing two marginal errors -- properly asymmetric, and reflects each fit's own posterior shape instead of a
    # Gaussian approximation.
    kt_diff_pct, kt_diff_lo, kt_diff_hi, f_diff_pct, f_diff_lo, f_diff_hi = paired_ratio_pct(joint_best, gbm_best)

    return ComparisonRow(
        grb_name=GRB_DISPLAY_NAME,
        episode=ep,
        t_start_s=joint_interval.start,
        t_stop_s=joint_interval.end,
        model_name_joint=joint_best.name,
        model_name_gbm=gbm_best.name,
        models_agree=models_agree,
        kt_bb_keV_joint=kt_j,
        kt_bb_err_keV_joint=kt_j_err,
        kt_bb_keV_gbm=kt_g,
        kt_bb_err_keV_gbm=kt_g_err,
        kt_consistent_1sigma=kt_consistent,
        f_bb_joint=f_j,
        f_bb_err_lower_joint=f_j_lo,
        f_bb_err_upper_joint=f_j_hi,
        f_bb_gbm=f_g,
        f_bb_err_lower_gbm=f_g_lo,
        f_bb_err_upper_gbm=f_g_hi,
        cstat_joint=joint_best.cstat,
        dof_joint=joint_best.dof,
        cstat_gbm=gbm_best.cstat,
        dof_gbm=gbm_best.dof,
        delta_cstat_joint=delta_joint,
        delta_cstat_gbm=delta_gbm,
        delta_cstat_threshold=DELTA_CSTAT_THRESHOLD,
        bb_required_joint=bb_required_joint,
        bb_required_gbm=bb_required_gbm,
        any_bb_step_clears_threshold_gbm=any_bb_clears_gbm,
        recorded_winner_is_bb_gbm=recorded_winner_is_bb_gbm,
        status_consistent_with_threshold_gbm=status_consistent_gbm,
        kt_bb_diff_pct=kt_diff_pct,
        kt_bb_diff_pct_err_lower=kt_diff_lo,
        kt_bb_diff_pct_err_upper=kt_diff_hi,
        f_bb_diff_pct=f_diff_pct,
        f_bb_diff_pct_err_lower=f_diff_lo,
        f_bb_diff_pct_err_upper=f_diff_hi,
        n_samples=N_SAMPLES,
        seed=SEED,
    )


def make_plot(df, episode_markers, path_stub="gbm_only_refit_kt_comparison"):
    """kT_BB vs episode, joint vs GBM-only, for episodes where both fits keep a BB model.

    Colour encodes fit type (Joint vs GBM-only); marker shape encodes episode,
    taken from EpisodeMarkerResolver (grb_research.grb_utils) -- the same
    per-episode marker convention used throughout this project -- rather than
    being invented locally. The two dimensions get two separate legends
    (CLAUDE.md's "identify burst/episode/model" requirement is met instead by
    naming the model on each x-tick, since burst is fixed and every episode
    uses one model in both fits -- models_agree is True throughout).
    """
    update_style()

    plotted = df.dropna(subset=["kt_bb_keV_joint"])
    fig, ax = plt.subplots()

    x = np.arange(len(plotted))
    width = 0.14

    for i, (_, row) in enumerate(plotted.iterrows()):
        marker = episode_markers[row["episode"]]
        ax.errorbar(
            x[i] - width,
            row["kt_bb_keV_joint"],
            yerr=row["kt_bb_err_keV_joint"],
            fmt=marker,
            mfc="none",
            mec=JOINT_COLOR,
            ecolor=JOINT_COLOR,
            markersize=MARKER_SIZE,
            markeredgewidth=MARKER_EDGE_WIDTH,
            linewidth=LINE_WIDTH,
        )
        if np.isfinite(row["kt_bb_keV_gbm"]):
            ax.errorbar(
                x[i] + width,
                row["kt_bb_keV_gbm"],
                yerr=row["kt_bb_err_keV_gbm"],
                fmt=marker,
                mfc="none",
                mec=GBM_COLOR,
                ecolor=GBM_COLOR,
                markersize=MARKER_SIZE,
                markeredgewidth=MARKER_EDGE_WIDTH,
                linewidth=LINE_WIDTH,
            )

    ax.set_xticks(x)
    _ff = lambda x: x.replace("_", "+")
    ax.set_xticklabels([f"{row['episode']}\n({_ff(row['model_name_joint'])})"
                        for _, row in plotted.iterrows()])
    ax.set_xlabel("Episode", fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel(r"$kT_\mathrm{BB}$ [keV]", fontsize=LABEL_FONT_SIZE)
    ax.tick_params(labelsize=TICK_FONT_SIZE)

    fit_handles = [
        mlines.Line2D(
            [],
            [],
            marker="o",
            mfc="none",
            mec=JOINT_COLOR,
            color=JOINT_COLOR,
            linestyle="None",
            markeredgewidth=MARKER_EDGE_WIDTH,
            markersize=MARKER_SIZE,
            label="Joint (GBM+LAT)",
        ),
        mlines.Line2D(
            [],
            [],
            marker="o",
            mfc="none",
            mec=GBM_COLOR,
            color=GBM_COLOR,
            linestyle="None",
            markeredgewidth=MARKER_EDGE_WIDTH,
            markersize=MARKER_SIZE,
            label="GBM-only",
        ),
    ]
    fit_legend = ax.legend(
        handles=fit_handles,
        title="Fit type",
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_FONT_SIZE,
        loc="best",
    )
    ax.add_artist(fit_legend)
    save_fig(fig, path_stub)


def make_delta_cstat_plot(df, episode_markers, path_stub="gbm_only_refit_delta_cstat"):
    """Delta-C-stat (BASE -> BASE+BB improvement) vs episode, joint vs GBM-only.

    This is the paper's own statistical currency for whether BB is required
    at all (section-1-introduction.tex's Delta-C-stat >= 28.74 rule), plotted
    directly rather than left in the CSV -- it shows not just that kT_BB
    matches between the two fits, but that BB is strongly preferred in both,
    by a wide margin over threshold in every episode. Log-scaled y-axis: the
    values span about an order of magnitude (30-380), and the threshold line
    needs to stay legible at the low end.
    """
    update_style()

    plotted = df.dropna(subset=["delta_cstat_joint", "delta_cstat_gbm"])
    fig, ax = plt.subplots()

    x = np.arange(len(plotted))
    width = 0.14

    ax.axhline(DELTA_CSTAT_THRESHOLD, color="black", linewidth=LINE_WIDTH, linestyle="--")

    for i, (_, row) in enumerate(plotted.iterrows()):
        marker = episode_markers[row["episode"]]
        ax.plot(
            x[i] - width,
            row["delta_cstat_joint"],
            marker=marker,
            mfc="none",
            mec=JOINT_COLOR,
            color=JOINT_COLOR,
            mew=MARKER_EDGE_WIDTH,
            markersize=MARKER_SIZE,
            linestyle="None",
        )
        ax.plot(
            x[i] + width,
            row["delta_cstat_gbm"],
            marker=marker,
            mfc="none",
            mec=GBM_COLOR,
            mew=MARKER_EDGE_WIDTH,
            color=GBM_COLOR,
            markersize=MARKER_SIZE,
            linestyle="None",
        )

    ax.set_yscale("log")
    ax.set_xticks(x)
    _ff = lambda x: x.replace("_", "+")
    ax.set_xticklabels([f"{row['episode']}\n({_ff(row['model_name_joint'])})"
                        for _, row in plotted.iterrows()])
    ax.set_ylabel(r"$\Delta C\mathrm{-stat}$ (BASE $\to$ BASE+BB)", fontsize=LABEL_FONT_SIZE)
    ax.tick_params(labelsize=TICK_FONT_SIZE)

    fit_handles = [
        mlines.Line2D(
            [],
            [],
            marker="o",
            mfc="none",
            mec=JOINT_COLOR,
            mew=MARKER_EDGE_WIDTH,
            color=JOINT_COLOR,
            linestyle="None",
            markersize=MARKER_SIZE,
            label="Joint (GBM+LAT)",
        ),
        mlines.Line2D(
            [],
            [],
            marker="o",
            mfc="none",
            mec=GBM_COLOR,
            mew=MARKER_EDGE_WIDTH,
            color=GBM_COLOR,
            linestyle="None",
            markersize=MARKER_SIZE,
            label="GBM-only",
        ),
        mlines.Line2D(
            [], [], color="black", linewidth=LINE_WIDTH, linestyle="--", label=f"threshold = {DELTA_CSTAT_THRESHOLD}"
        ),
    ]
    fit_legend = ax.legend(
        handles=fit_handles,
        title="Fit type",
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_FONT_SIZE,
        loc="best",
    )
    ax.add_artist(fit_legend)
    save_fig(fig, path_stub)


def make_fractional_diff_plot(df, episode_markers, path_stub="gbm_only_refit_fractional_diff"):
    """Fractional difference, (GBM-only - Joint) / Joint, in kT_BB and f_BB per episode.

    Quantifies directly how much the GBM-only refit shifts each parameter
    relative to the joint fit -- the number the kT comparison plot leaves the
    reader to eyeball from two overlapping error bars. Errors come from
    paired_ratio_pct()'s own MC ratio distribution (percentiles of
    (GBM-only draw / Joint draw - 1), not a linearized combination of two
    marginal errors), so they are asymmetric where the underlying
    distribution is.
    """
    update_style()

    plotted = df.dropna(subset=["kt_bb_keV_joint", "kt_bb_keV_gbm"])
    fig, axes = plt.subplots(2, 1, sharex=True)

    panels = [
        (
            axes[0],
            "kt_bb_diff_pct",
            "kt_bb_diff_pct_err_lower",
            "kt_bb_diff_pct_err_upper",
            r"$\Delta kT_\mathrm{BB}$ [%]",
        ),
        (axes[1], "f_bb_diff_pct", "f_bb_diff_pct_err_lower", "f_bb_diff_pct_err_upper", r"$\Delta f_\mathrm{BB}$ [%]"),
    ]
    x = np.arange(len(plotted))
    for ax, col, err_lo_col, err_hi_col, ylabel in panels:
        ax.axhline(0.0, color="black", linewidth=LINE_WIDTH, linestyle="--")
        for i, (_, row) in enumerate(plotted.iterrows()):
            marker = episode_markers[row["episode"]]
            ax.errorbar(
                x[i],
                row[col],
                yerr=[[row[err_lo_col]], [row[err_hi_col]]],
                fmt=marker,
                mfc="none",
                mec=GBM_COLOR,
                ecolor=GBM_COLOR,
                capsize=CAP_SIZE,
                markersize=MARKER_SIZE,
                markeredgewidth=MARKER_EDGE_WIDTH,
                linewidth=LINE_WIDTH,
            )
        ax.set_ylabel(ylabel, fontsize=LABEL_FONT_SIZE)
        ax.tick_params(labelsize=TICK_FONT_SIZE)

    axes[1].set_xticks(x)
    _ff = lambda x: x.replace("_", "+")
    axes[1].set_xticklabels([f"{row['episode']}\n({_ff(row['model_name_joint'])})"
                             for _, row in plotted.iterrows()])

    save_fig(fig, path_stub)


def main():
    root = find_project_root()
    with open(root / "results.json") as f:
        data = json.load(f)

    joint = GRB.from_dictionary(name=JOINT_KEY, grb_data=data[JOINT_KEY])
    gbm = GRB.from_dictionary(name=GBM_KEY, grb_data=data[GBM_KEY])

    joint_eps = {episode_label(i): i for i in joint.intervals}
    gbm_eps = {episode_label(i): i for i in gbm.intervals}

    rows = []
    for ep in EPISODE_ORDER:
        if ep not in joint_eps or ep not in gbm_eps:
            continue
        rows.append(build_row(ep, joint_eps[ep], gbm_eps[ep]))

    df = pd.DataFrame([r.__dict__ for r in rows])
    df.to_csv("gbm_only_refit_comparison.csv", index=False)

    resolver = EpisodeMarkerResolver(t90_marker=T90_MARKER)
    episode_markers = {ep: resolver.resolve(joint_eps[ep]) for ep in EPISODE_ORDER if ep in joint_eps}

    make_plot(df, episode_markers)
    make_delta_cstat_plot(df, episode_markers)
    make_fractional_diff_plot(df, episode_markers)

    pd.set_option("display.width", 200)
    print(df[["episode", "model_name_joint", "model_name_gbm", "models_agree", "kt_bb_keV_joint", "kt_bb_keV_gbm",
              "kt_consistent_1sigma", "any_bb_step_clears_threshold_gbm", "recorded_winner_is_bb_gbm",
              "status_consistent_with_threshold_gbm"]])
    print(df[["episode", "kt_bb_diff_pct", "kt_bb_diff_pct_err_lower", "kt_bb_diff_pct_err_upper",
              "f_bb_diff_pct", "f_bb_diff_pct_err_lower", "f_bb_diff_pct_err_upper"]])


if __name__ == "__main__":
    main()
