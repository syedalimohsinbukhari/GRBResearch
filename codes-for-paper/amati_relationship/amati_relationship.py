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
from grb_research import (
    ModelSet,
    episode_order,
    find_project_root,
    get_rng,
    prepare_grbs,
    seed_from_name,
    update_style,
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
)
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

    a.legend(loc="best", title=f"GRB{grb_list[i]}", fontsize=LEGEND_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE)

    ep_total.append(_[0])
    ei_total.append(_[1])
    ep_label.append(_[2])
    model_list.append(_[3])
    ep_err_total.append(_[4])
    ei_err_total.append(_[5])
    g_name.append([f"GRB{grb_list[i]}"] * len(_[0]))

# ---------------------------------------------------------------------------
# Unknown-redshift GRBs — redshift locus across z = 1, 3, 5, 7
# ---------------------------------------------------------------------------

Z_VALUES_UNKNOWN = (1, 3, 5, 7)

# Previously discarded: plot_unknown_redshift_grb() already computes (Ep, Eiso) at every z in
# Z_VALUES_UNKNOWN for these three bursts' episodes, purely to draw the dashed locus -- its return
# value was never captured. Both quantities increase monotonically with z (Eiso via d_L(z), Ep,rest
# via the (1+z) factor), so the z=1/z=7 endpoints of each episode's track are its full range; no
# need to re-derive anything, just keep what was already computed and take min/max over it.
unknown_z_rows = []

for idx, m_ in enumerate(grb_best[1:]):
    for m in m_:
        # `uz_` prefix deliberately distinct from ep_total/ei_total/model_list/ep_label above --
        # this is plain top-level script code with no function scope, so a same-named unpacking
        # target here would silently clobber those variables before the CSV write further down
        # reads them (exactly what happened here originally: `model_list` collided and the
        # known-redshift GRB080916C rows in amati_relationship.csv came out with a garbled
        # Model column, sourced from whatever unknown-z episode ran last, not from GRB080916C).
        uz_ep_all, uz_ei_all, uz_ep_name, uz_model_list, uz_ep_err_all, uz_ei_err_all = plot_unknown_redshift_grb(
            models=[m],
            t90_marker=t90_markers[idx + 1],
            z_values=Z_VALUES_UNKNOWN,
            n_grid=n_grid,
            n_sample=n_sample,
            rng=rng,
            axis=ax[idx + 1],
        )

        # E_i,peak = (1+z) * E_p,obs is an exact identity (ep_samples in _compute_ep_eiso is drawn
        # once per z-call and merely rescaled by (1+z), never refit) -- showing it swept over z
        # carries no information beyond E_p,obs itself plus the caption stating that scaling
        # (2026-09-30 review). E_p,obs is recovered from z=1 (multiplier (1+1)=2); independently
        # cross-checking against z=7 (multiplier 8) below is a free MC-noise diagnostic, since
        # both estimate the same quantity from independent draws -- not stored in the CSV, printed
        # for review.
        ep_obs = uz_ep_all[0] / (1 + Z_VALUES_UNKNOWN[0])
        ep_obs_err_lower = uz_ep_err_all[0][0, 0] / (1 + Z_VALUES_UNKNOWN[0])
        ep_obs_err_upper = uz_ep_err_all[0][1, 0] / (1 + Z_VALUES_UNKNOWN[0])
        ep_obs_crosscheck = uz_ep_all[-1] / (1 + Z_VALUES_UNKNOWN[-1])
        print(
            f"  E_p,obs cross-check {uz_ep_name[0].split('$')[0]:<5s}: "
            f"from z={Z_VALUES_UNKNOWN[0]}: {ep_obs:.4f}, from z={Z_VALUES_UNKNOWN[-1]}: {ep_obs_crosscheck:.4f} "
            f"({100 * abs(ep_obs - ep_obs_crosscheck) / ep_obs:.2f}% apart -- MC noise floor, not a real disagreement)"
        )

        row = {
            "GRBName": f"GRB{grb_list[idx + 1]}",
            "Model": uz_model_list[0],
            # ep_name carries _episode_label()'s LaTeX model suffix (e.g. "TR1$_\text{band+bb}$");
            # stripped the same way the known-redshift ep_label list is stripped above, so
            # EpisodeName is a bare label ("TR1") every downstream table-sorter can match on.
            "EpisodeName": uz_ep_name[0].split("$")[0],
            f"E_p_obs__{EP_NORM:.0e}_keV": ep_obs,
            f"E_p_obs_err_lower__{EP_NORM:.0e}_keV": ep_obs_err_lower,
            f"E_p_obs_err_upper__{EP_NORM:.0e}_keV": ep_obs_err_upper,
        }
        # E_iso does NOT scale as a fixed multiple of z the way E_i,peak does -- it runs through
        # d_L(z)^2 and a per-episode K-correction, so each z gets its own column (2026-09-30 review).
        for z_idx, z in enumerate(Z_VALUES_UNKNOWN):
            row[f"E_0_iso_z{z}__{EI_NORM:.0e}_erg"] = uz_ei_all[z_idx]
            row[f"E_0_iso_z{z}_err_lower__{EI_NORM:.0e}_erg"] = uz_ei_err_all[z_idx][0, 0]
            row[f"E_0_iso_z{z}_err_upper__{EI_NORM:.0e}_erg"] = uz_ei_err_all[z_idx][1, 0]
        unknown_z_rows.append(row)

    ax[idx + 1].legend(
        loc="best", title=f"GRB{grb_list[idx + 1]}", fontsize=LEGEND_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE
    )

# Sorted by temporal episode order (T90, EX0, TR1..TRn, EX1), not source-data order -- see
# grb_research.episode_order's docstring for why that matters (EX1 vs TR5 ordering).
unknown_z_rows.sort(key=lambda r: (r["GRBName"], episode_order(r["EpisodeName"])))

unknown_z_df = pd.DataFrame(unknown_z_rows)
unknown_z_df["n_samples"] = n_sample
unknown_z_df["seed"] = SEED
unknown_z_df.to_csv("amati_relationship_unknown_z.csv", index=False)

ep_total = list(chain.from_iterable(ep_total))
ei_total = list(chain.from_iterable(ei_total))
ep_label = list(chain.from_iterable(ep_label))
ep_label = [i.split("$")[0] for i in ep_label]
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
# Sorted by temporal episode order, not results.json's own row order -- see
# grb_research.episode_order's docstring (this is the table where EX1 used to print before TR5).
q = q.sort_values(by="EpisodeName", key=lambda col: col.map(episode_order), kind="stable")
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
