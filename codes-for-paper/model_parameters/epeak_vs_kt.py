"""Created on Apr 01 14:02:04 2026"""

import json

import matplotlib.pyplot as plt
import numpy as np

from grb_research import update_style
from utils import (
    extract_kt_epeak_from_models,
    fit_and_plot_odr,
    find_project_root,
    GRBCatalog,
    short_to_long,
    EpisodeTypes,
    LEGEND_FONT_SIZE,
    LABEL_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
)

# -- Load data ----------------------------------------------------------------

update_style()

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"
grb_name = ["080916C", "140206B", "131014A", "231129C"]

with open(result_file, "r") as f:
    data = json.load(f)

grb = GRBCatalog.from_iterable(grb_list=grb_name, data=data, name_mapping=short_to_long)
grb_full_name = [short_to_long[i] for i in grb_name]

# -- Fetch models (manual specification) --------------------------------------

models_080916C = [
    grb[grb_full_name[0]].get_model("SBPL_BB", interval=EpisodeTypes.T90),
    grb[grb_full_name[0]].get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb[grb_full_name[0]].get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb[grb_full_name[0]].get_model("SBPL_BB", interval=EpisodeTypes.TR, tr_index=3),
]

bb_status_080916C = ["Y", "Y", "Y", "X"]

models_140206B = [
    grb[grb_full_name[1]].get_model("BAND_BB", EpisodeTypes.T90),
    grb[grb_full_name[1]].get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb[grb_full_name[1]].get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb[grb_full_name[1]].get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=2),
    grb[grb_full_name[1]].get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=3),
]

bb_status_140206B = ["X", "Y", "Y", "X", "X"]

models_131014A = [
    grb[grb_full_name[2]].get_model("BAND_BB", EpisodeTypes.T90),
    grb[grb_full_name[2]].get_model("BAND_BB", interval=EpisodeTypes.EX0),
    grb[grb_full_name[2]].get_model("BAND_BB", interval=EpisodeTypes.TR, tr_index=1),
    grb[grb_full_name[2]].get_model("SBPL_BB", interval=EpisodeTypes.TR, tr_index=2),
    grb[grb_full_name[2]].get_model("BAND_BB", interval=EpisodeTypes.EX1),
]
bb_status_131014A = ["Y", "Y", "Y", "Y", "Y"]

models_231129C = [
    grb[grb_full_name[3]].get_model("SBPL_BB", EpisodeTypes.T90),
    grb[grb_full_name[3]].get_model("SBPL_BB", EpisodeTypes.EX0),
    grb[grb_full_name[3]].get_model("SBPL_BB", EpisodeTypes.TR, tr_index=1),
    grb[grb_full_name[3]].get_model("SBPL_BB", EpisodeTypes.EX1)
]

bb_status_231129C = ["Y", "Y", "Y", "X"]

# -- Extract kT / E_peak -----------------------------------------------------

kt_080916C, ep_080916C, mkr_080916C, clr_080916C, lbl_080916C = extract_kt_epeak_from_models(models_080916C)
kt_140206B, ep_140206B, mkr_140206B, clr_140206B, lbl_140206B = extract_kt_epeak_from_models(models_140206B)
kt_131014A, ep_131014A, mkr_131014A, clr_131014A, lbl_131014A = extract_kt_epeak_from_models(models_131014A)
kt_231129C, ep_231129C, mkr_231129C, clr_231129C, lbl_231129C = extract_kt_epeak_from_models(models_231129C)
# -- Plot ---------------------------------------------------------------------

f, ax = plt.subplots(4, 1, figsize=(6, 8))

ax = np.array(ax).flatten()

grb_panels = [
    (ax[0], kt_080916C, ep_080916C, mkr_080916C, clr_080916C, lbl_080916C, bb_status_080916C),
    (ax[1], kt_140206B, ep_140206B, mkr_140206B, clr_140206B, lbl_140206B, bb_status_140206B),
    (ax[2], kt_131014A, ep_131014A, mkr_131014A, clr_131014A, lbl_131014A, bb_status_131014A),
    (ax[3], kt_231129C, ep_231129C, mkr_231129C, clr_231129C, lbl_231129C, bb_status_231129C)
]

# ax[2].get_xlim()

for a, kt, ep, mkrs, clrs, lbls, status_list in grb_panels:
    for kt_i, ep_i, mkr, clr, lbl, status in zip(kt, ep, mkrs, clrs, lbls, status_list):
        a.errorbar(
            kt_i[1],
            ep_i[1],
            xerr=[[kt_i[0]], [kt_i[2]]],
            yerr=[[ep_i[0]], [ep_i[2]]],
            fmt=mkr,
            mfc="w" if status == "X" else None,
            ms=8,
            capsize=5,
            color=clr,
            linestyle="--" if status == "X" else "-",
            label=lbl,
        )

# -- ODR fits -----------------------------------------------------------------

# GRB 080916C: fit "full" points only
fit_and_plot_odr(kt_080916C, ep_080916C, ax[0], color='teal', annotation_xy=(0.05, 0.12))
full_mask_080916C = [s == "Y" for s in bb_status_080916C]
fit_and_plot_odr(kt_080916C, ep_080916C, ax[0], mask=full_mask_080916C)
fit_and_plot_odr(kt_140206B, ep_140206B, ax[1], color='teal', annotation_xy=(0.05, 0.82))
full_mask_140206B = [s == "Y" for s in bb_status_140206B]
fit_and_plot_odr(kt_140206B, ep_140206B, ax[1], mask=full_mask_140206B)
fit_and_plot_odr(kt_131014A, ep_131014A, ax[2])
fit_and_plot_odr(kt_231129C, ep_231129C, ax[3], color='teal', annotation_xy=(0.05, 0.12))
full_mask_231129C = [s == "Y" for s in bb_status_231129C]
fit_and_plot_odr(kt_231129C, ep_231129C, ax[3], mask=full_mask_231129C,
                 annotation_xy=(0.05, 0.22))

print(grb_name[0])
print(kt_080916C)
print(ep_080916C)
print(full_mask_080916C)
print(grb_name[2])
print(kt_131014A)
print(ep_131014A)
print(grb_name[1])
print(kt_140206B)
print(ep_140206B)
print(full_mask_140206B)
print(grb_name[3])
print(kt_231129C)
print(ep_231129C)
print(full_mask_231129C)


ax[-1].set_xlabel("kT [keV]", fontsize=LABEL_FONT_SIZE)
[a.set_ylabel(r"$E_\text{peak}$ [keV]", fontsize=LABEL_FONT_SIZE) for a in ax]
[
    a.legend(fontsize=LEGEND_FONT_SIZE, title=f"GRB{grb_name[i]}", title_fontsize=LEGEND_TITLE_FONT_SIZE)
    for i, a in enumerate(ax)
]
[a.grid(True, which="both", alpha=0.5, ls="--") for a in ax]
f.tight_layout()

plt.show()
# for ext in ["png", "pdf"]:
#     plt.savefig(f"epeak_vs_kt.{ext}", dpi=300)
# plt.close()
