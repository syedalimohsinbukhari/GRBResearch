"""Created on Apr 01 14:02:04 2026"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grb_research import (
    EpisodeTypes,
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
    MARKER_SIZE,
    SAVE_DPI,
    find_project_root,
    get_rng,
    prepare_grbs,
    seed_from_name,
    update_style,
)
from grb_research.grb_constants import CAP_SIZE, N_SAMPLES
from utils import extract_kt_epeak_from_models, fit_and_plot_odr

# -- Load data ----------------------------------------------------------------

update_style()

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"
grb_name = ["080916C", "140206B", "131014A", "231129C"]

_, _, grb_objs = prepare_grbs(grb_name, result_file, get_best=False)
grb_080916C, grb_140206B, grb_131014A, grb_231129C = grb_objs

# The SBPL -> Band E_peak conversion inside extract_kt_epeak_from_models is Monte-Carlo;
# BAND/CPL rows read e_peak directly from the fit and carry no MC provenance.
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)

# -- Fetch models (curated selection) -----------------------------------------
#
# Each list picks a SPECIFIC model per episode via get_model(name, interval, tr_index)
# -- not necessarily that episode's BEST-fit model -- and deliberately omits some
# episodes entirely (e.g. GRB080916C: only T90/EX0/TR1/TR3 of its 8 episodes). This is
# a curated selection, preserved exactly as-is; see epeak_vs_kt.md for what motivated it.
#
# bb_status is NOT a "has BB" flag (every model here already has BB, or kT wouldn't
# exist): "Y" marks the points used in the second, restricted ODR fit per panel; "X"
# points still appear on the plot but are excluded from that fit.

models_080916C = [
    grb_080916C.get_model("SBPL_BB", interval=EpisodeTypes.T90),
    grb_080916C.get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb_080916C.get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb_080916C.get_model("SBPL_BB", interval=EpisodeTypes.TR, tr_index=3),
]
bb_status_080916C = ["Y", "Y", "Y", "X"]

models_140206B = [
    grb_140206B.get_model("BAND_BB", interval=EpisodeTypes.T90),
    grb_140206B.get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb_140206B.get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb_140206B.get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=2),
    grb_140206B.get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=3),
]
bb_status_140206B = ["X", "Y", "Y", "X", "X"]

models_131014A = [
    grb_131014A.get_model("BAND_BB", interval=EpisodeTypes.T90),
    grb_131014A.get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb_131014A.get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb_131014A.get_model("SBPL_BB", interval=EpisodeTypes.TR, tr_index=2),
    grb_131014A.get_model("BAND_BB", interval=EpisodeTypes.EX1),
]
bb_status_131014A = ["Y", "Y", "Y", "Y", "Y"]

models_231129C = [
    grb_231129C.get_model("SBPL_BB", interval=EpisodeTypes.T90),
    grb_231129C.get_model("SBPL_BB", interval=EpisodeTypes.EX0),
    grb_231129C.get_model("SBPL_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb_231129C.get_model("SBPL_BB", interval=EpisodeTypes.EX1),
]
bb_status_231129C = ["Y", "Y", "Y", "X"]

# -- Extract kT / E_peak -----------------------------------------------------

kt_080916C, ep_080916C, mkr_080916C, clr_080916C, lbl_080916C, mc_080916C = extract_kt_epeak_from_models(
    models_080916C, rng=rng
)
kt_140206B, ep_140206B, mkr_140206B, clr_140206B, lbl_140206B, mc_140206B = extract_kt_epeak_from_models(
    models_140206B, rng=rng
)
kt_131014A, ep_131014A, mkr_131014A, clr_131014A, lbl_131014A, mc_131014A = extract_kt_epeak_from_models(
    models_131014A, rng=rng
)
kt_231129C, ep_231129C, mkr_231129C, clr_231129C, lbl_231129C, mc_231129C = extract_kt_epeak_from_models(
    models_231129C, rng=rng
)

# -- Plot ---------------------------------------------------------------------

f, ax = plt.subplots(2, 2, figsize=(10.5, 8))
ax = np.array(ax).flatten()

grb_panels = [
    (ax[0], kt_080916C, ep_080916C, mkr_080916C, clr_080916C, lbl_080916C, bb_status_080916C),
    (ax[1], kt_140206B, ep_140206B, mkr_140206B, clr_140206B, lbl_140206B, bb_status_140206B),
    (ax[2], kt_131014A, ep_131014A, mkr_131014A, clr_131014A, lbl_131014A, bb_status_131014A),
    (ax[3], kt_231129C, ep_231129C, mkr_231129C, clr_231129C, lbl_231129C, bb_status_231129C),
]

for a, kt, ep, mkrs, clrs, lbls, status_list in grb_panels:
    for kt_i, ep_i, mkr, clr, lbl, status in zip(kt, ep, mkrs, clrs, lbls, status_list):
        a.errorbar(
            kt_i[1],
            ep_i[1],
            xerr=[[kt_i[0]], [kt_i[2]]],
            yerr=[[ep_i[0]], [ep_i[2]]],
            fmt=mkr,
            mfc="w" if status == "X" else None,
            ms=MARKER_SIZE,
            capsize=CAP_SIZE,
            color=clr,
            linestyle="--" if status == "X" else "-",
            label=lbl,
        )

# -- ODR fits -----------------------------------------------------------------

# GRB 080916C: "all" points, then "restricted" to the Y-marked full set.
odr_all_080916C = fit_and_plot_odr(kt_080916C, ep_080916C, ax[0], color="teal", annotation_xy=(0.05, 0.12))
full_mask_080916C = [s == "Y" for s in bb_status_080916C]
odr_full_080916C = fit_and_plot_odr(kt_080916C, ep_080916C, ax[0], mask=full_mask_080916C)

odr_all_140206B = fit_and_plot_odr(kt_140206B, ep_140206B, ax[1], color="teal", annotation_xy=(0.05, 0.82))
full_mask_140206B = [s == "Y" for s in bb_status_140206B]
odr_full_140206B = fit_and_plot_odr(kt_140206B, ep_140206B, ax[1], mask=full_mask_140206B)

# GRB 131014A: every point is already "Y" -- a single unmasked fit covers both cases.
odr_all_131014A = fit_and_plot_odr(kt_131014A, ep_131014A, ax[2])

odr_all_231129C = fit_and_plot_odr(kt_231129C, ep_231129C, ax[3], color="teal", annotation_xy=(0.05, 0.12))
full_mask_231129C = [s == "Y" for s in bb_status_231129C]
odr_full_231129C = fit_and_plot_odr(kt_231129C, ep_231129C, ax[3], mask=full_mask_231129C, annotation_xy=(0.05, 0.22))

ax[-1].set_xlabel("kT [keV]", fontsize=LABEL_FONT_SIZE)
[a.set_ylabel(r"$E_\text{peak}$ [keV]", fontsize=LABEL_FONT_SIZE) for i, a in enumerate(ax) if i % 2 == 0]
[
    a.legend(fontsize=LEGEND_FONT_SIZE, title=f"GRB{grb_name[i]}", title_fontsize=LEGEND_TITLE_FONT_SIZE)
    for i, a in enumerate(ax)
]
# [a.grid(True, which="both", alpha=0.5, ls="--") for a in ax]
f.tight_layout()

for ext in ["png", "pdf"]:
    plt.savefig(f"epeak_vs_kt.{ext}", dpi=SAVE_DPI)

# ─── Save the values ──────────────────────────────────────────────────────

grb_sets = [
    ("GRB080916C", models_080916C, bb_status_080916C, kt_080916C, ep_080916C, mc_080916C),
    ("GRB140206B", models_140206B, bb_status_140206B, kt_140206B, ep_140206B, mc_140206B),
    ("GRB131014A", models_131014A, bb_status_131014A, kt_131014A, ep_131014A, mc_131014A),
    ("GRB231129C", models_231129C, bb_status_231129C, kt_231129C, ep_231129C, mc_231129C),
]

points_rows = []
for full_name, models, status_list, kt_arr, ep_arr, mc_arr in grb_sets:
    for model, status, kt_i, ep_i, mc in zip(models, status_list, kt_arr, ep_arr, mc_arr):
        idx = "" if model.interval.index is None else str(model.interval.index)
        points_rows.append(
            {
                "grb_name": full_name,
                "episode": f"{model.interval.kind.value}{idx}",
                "model_name": model.name,
                "kt_keV": kt_i[1],
                "kt_err_lower_keV": kt_i[0],
                "kt_err_upper_keV": kt_i[2],
                "e_peak_keV": ep_i[1],
                "e_peak_err_lower_keV": ep_i[0],
                "e_peak_err_upper_keV": ep_i[2],
                "used_in_restricted_fit": status == "Y",
                "n_samples": N_SAMPLES if mc else np.nan,
                "seed": SEED if mc else np.nan,
            }
        )

pd.DataFrame(points_rows).to_csv("epeak_vs_kt_points.csv", index=False)
print(f"Saved: epeak_vs_kt_points.csv  ({len(points_rows)} rows)")

odr_rows = [
    {
        "grb_name": "GRB080916C",
        "fit_type": "all",
        "n_points": len(kt_080916C),
        "slope": odr_all_080916C.beta[0],
        "slope_err": odr_all_080916C.sd_beta[0],
        "intercept": odr_all_080916C.beta[1],
        "intercept_err": odr_all_080916C.sd_beta[1],
    },
    {
        "grb_name": "GRB080916C",
        "fit_type": "restricted",
        "n_points": int(sum(full_mask_080916C)),
        "slope": odr_full_080916C.beta[0],
        "slope_err": odr_full_080916C.sd_beta[0],
        "intercept": odr_full_080916C.beta[1],
        "intercept_err": odr_full_080916C.sd_beta[1],
    },
    {
        "grb_name": "GRB140206B",
        "fit_type": "all",
        "n_points": len(kt_140206B),
        "slope": odr_all_140206B.beta[0],
        "slope_err": odr_all_140206B.sd_beta[0],
        "intercept": odr_all_140206B.beta[1],
        "intercept_err": odr_all_140206B.sd_beta[1],
    },
    {
        "grb_name": "GRB140206B",
        "fit_type": "restricted",
        "n_points": int(sum(full_mask_140206B)),
        "slope": odr_full_140206B.beta[0],
        "slope_err": odr_full_140206B.sd_beta[0],
        "intercept": odr_full_140206B.beta[1],
        "intercept_err": odr_full_140206B.sd_beta[1],
    },
    {
        "grb_name": "GRB131014A",
        "fit_type": "all",
        "n_points": len(kt_131014A),
        "slope": odr_all_131014A.beta[0],
        "slope_err": odr_all_131014A.sd_beta[0],
        "intercept": odr_all_131014A.beta[1],
        "intercept_err": odr_all_131014A.sd_beta[1],
    },
    {
        "grb_name": "GRB231129C",
        "fit_type": "all",
        "n_points": len(kt_231129C),
        "slope": odr_all_231129C.beta[0],
        "slope_err": odr_all_231129C.sd_beta[0],
        "intercept": odr_all_231129C.beta[1],
        "intercept_err": odr_all_231129C.sd_beta[1],
    },
    {
        "grb_name": "GRB231129C",
        "fit_type": "restricted",
        "n_points": int(sum(full_mask_231129C)),
        "slope": odr_full_231129C.beta[0],
        "slope_err": odr_full_231129C.sd_beta[0],
        "intercept": odr_full_231129C.beta[1],
        "intercept_err": odr_full_231129C.sd_beta[1],
    },
]
pd.DataFrame(odr_rows).to_csv("epeak_vs_kt_odr_fits.csv", index=False)
print(f"Saved: epeak_vs_kt_odr_fits.csv  ({len(odr_rows)} rows)")
