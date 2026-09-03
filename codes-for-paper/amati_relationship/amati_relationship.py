"""Created on Jan 24 02:32:18 2026 — refactored Mar 27 2026"""

from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from amati_helpers import (
    amati_relationship_dirirsa2019,
    plot_unknown_redshift_grb,
    plot_grbs_over_amati_relationship,
    EP_NORM,
    EI_NORM,
)
from grb_research import ModelSet, find_project_root, get_rng, prepare_grbs, seed_from_name, update_style, \
    LABEL_FONT_SIZE, LEGEND_FONT_SIZE, LEGEND_TITLE_FONT_SIZE
from grb_research.grb_utils import save_fig

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

update_style()

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

grb_list = ["080916C", "131014A", "140206B", "231129C"]
gc, grb_list_long, grb_objects, grb_best = prepare_grbs(grb_list=grb_list, result_file=result_file, get_best=True)

grb_best = [ModelSet([i for i in j if i.name != "PL"]) for j in grb_best]

redshifts = [4.35]
t90_markers = ["o", "s", "X", "D"]

# ---------------------------------------------------------------------------
# Sampling config
# ---------------------------------------------------------------------------

n_sample = 5_000
n_grid = 1000
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

f, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 10), sharex=True, sharey=True)
ax = ax.flatten()

amati_kw = dict(x_lim=(180, 7e4), y_lim=(6.25e50, 4.75e55), num_points=n_grid, use_average=True)
for a in ax:
    amati_relationship_dirirsa2019(axis=a, **amati_kw)

# ---------------------------------------------------------------------------
# Known-redshift GRBs — one per subplot
# ---------------------------------------------------------------------------

ep_total, ei_total = [], []
ep_err_total, ei_err_total = [], []
ep_label, g_name = [], []
model_list = []

for i, a in enumerate([ax[0]]):
    _ = plot_grbs_over_amati_relationship(
        best_model_list=[grb_best[i]],
        redshift_list=[redshifts[i]],
        t90_marker_list=[t90_markers[i]],
        n_grid=n_grid,
        n_sample=n_sample,
        rng=rng,
        axis=a,
    )

    a.legend(loc="best", title=f"GRB{grb_list[i]}",
             fontsize=LEGEND_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE)

    ep_total.append(_[0])
    ei_total.append(_[1])
    ep_label.append(_[2])
    model_list.append(_[3])
    ep_err_total.append(_[4])
    ei_err_total.append(_[5])
    g_name.append([f"GRB{grb_list[i]}"] * len(_[0]))

# ---------------------------------------------------------------------------
# Unknown-redshift GRB (GRB150210A) — redshift locus across z = 1, 3, 5, 7
# ---------------------------------------------------------------------------

for idx, m_ in enumerate(grb_best[1:]):
    for m in m_:
        plot_unknown_redshift_grb(
            models=[m],
            t90_marker=t90_markers[idx + 1],
            z_values=(1, 3, 5, 7),
            n_grid=n_grid,
            n_sample=n_sample,
            rng=rng,
            axis=ax[idx + 1],
        )

    ax[idx + 1].legend(loc="best", title=f"GRB{grb_list[idx + 1]}",
             fontsize=LEGEND_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE)

ep_total = list(chain.from_iterable(ep_total))
ei_total = list(chain.from_iterable(ei_total))
ep_label = list(chain.from_iterable(ep_label))
model_list = list(chain.from_iterable(model_list))
g_name = list(chain.from_iterable(g_name))
ep_err_total = list(chain.from_iterable(ep_err_total))
ei_err_total = list(chain.from_iterable(ei_err_total))

# Convert to arrays
ep_total, ei_total = np.array(ep_total), np.array(ei_total)
ep_label = np.array(ep_label)
g_name = np.array(g_name)
model_list = np.array(model_list)

# Extract asymmetric errors from (2, 1) arrays
# Known-redshift GRBs have errors, unknown-redshift ones will get NaNs
ep_err_lower = np.array([err[0, 0] for err in ep_err_total])
ep_err_upper = np.array([err[1, 0] for err in ep_err_total])
ei_err_lower = np.array([err[0, 0] for err in ei_err_total])
ei_err_upper = np.array([err[1, 0] for err in ei_err_total])

# Pad error arrays with NaNs for unknown-redshift GRB entries
n_unknown = len(ep_total) - len(ep_err_lower)
if n_unknown > 0:
    ep_err_lower = np.concatenate([ep_err_lower, np.full(shape=n_unknown, fill_value=np.nan)])
    ep_err_upper = np.concatenate([ep_err_upper, np.full(shape=n_unknown, fill_value=np.nan)])
    ei_err_lower = np.concatenate([ei_err_lower, np.full(shape=n_unknown, fill_value=np.nan)])
    ei_err_upper = np.concatenate([ei_err_upper, np.full(shape=n_unknown, fill_value=np.nan)])

n_samples_col = np.full(len(g_name), n_sample)
seed_col = np.full(len(g_name), SEED)

q = pd.DataFrame(
    [
        g_name,
        model_list,
        ep_label,
        ep_total,  # already in units of EP_NORM (10³ keV)
        ep_err_lower,
        ep_err_upper,
        ei_total,  # already in units of EI_NORM (10⁵² erg)
        ei_err_lower,
        ei_err_upper,
        n_samples_col,
        seed_col,
    ]
).T
q.columns = [
    "GRBName",
    "Model",
    "EpisodeName",
    f"E_i_peak__{EP_NORM:.0e}_keV",
    f"E_i_peak_err_lower__{EP_NORM:.0e}_keV",
    f"E_i_peak_err_upper__{EP_NORM:.0e}_keV",
    f"E_0_iso__{EI_NORM:.0e}_erg",
    f"E_0_iso_err_lower__{EI_NORM:.0e}_erg",
    f"E_0_iso_err_upper__{EI_NORM:.0e}_erg",
    "n_samples",
    "seed",
]
q.to_csv("amati_relationship.csv", index=False)

# ---------------------------------------------------------------------------
# Shared axis labels and export
# ---------------------------------------------------------------------------

for a in ax[2:]:
    a.set_xlabel(r"$E_{i,\mathrm{peak}}^{2}$ [keV]", fontsize=LABEL_FONT_SIZE)
    # a.tick_params(axis="both")
for a in ax[::2]:
    a.set_ylabel(r"$E_\mathrm{iso}^{52}$ [erg]", fontsize=LABEL_FONT_SIZE)
#     a.tick_params(axis="both")

save_fig(f, "amati_relationship")
