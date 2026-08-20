"""Created on Dec 17 13:22:15 2025"""

from typing import Tuple

import numpy as np
from matplotlib import pyplot as plt, ticker

from grb_research import update_style, GRID_LINESTYLE, GRID_ALPHA, SAVE_DPI
from utils import (
    extract_parameter,
    convert_sbpl_to_band,
    find_project_root,
    prepare_grbs,
    LABEL_FONT_SIZE,
    TICK_FONT_SIZE,
    LEGEND_FONT_SIZE,
    ModelSet,
    plot_per_episode,
    save_value_error_as_parquet,
)

update_style()


def extract_peak_energy(
    best_model: ModelSet, seed: int = 42
) -> Tuple[np.typing.ArrayLike, np.typing.ArrayLike, np.typing.ArrayLike]:
    """Extract peak energy (Ep) values and asymmetric errors from a model set."""
    values, errors_hi, errors_lo = [], [], []

    for model in best_model:
        if "SBPL" in model.name:
            err_lo, median, err_hi = convert_sbpl_to_band(model, seed=seed)
        else:
            result = extract_parameter(model, "e_peak")
            median, err_hi = result if result is not None else (np.nan, np.nan)
            err_lo = err_hi

        values.append(median)
        errors_hi.append(err_hi)
        errors_lo.append(err_lo)

    return np.asarray(values), np.asarray(errors_hi), np.asarray(errors_lo)


grb_list = ["080916C", "131014A", "140206B", "231129C"]

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

_, grb_list_long, grb_objs, grb_best = prepare_grbs(grb_list, result_file, get_best=True)

has_BB = [['BB' in i.name for i in t_.get_all_best_models()] for t_ in grb_objs]

starts, ends, diffs, midpoints = zip(
    *(g.intervals.extract_interval_arrays(return_include=("diff", "midpoint")) for g in grb_objs)
)
ep_values, ep_errors_hi, ep_errors_lo = zip(*(extract_peak_energy(best) for best in grb_best))

_, ax = plt.subplots(2, 2, figsize=(10, 7.5))
ax = ax.flatten()

for i, axis in enumerate(ax):
    plot_per_episode(values=ep_values[i], errors=[ep_errors_hi[i], ep_errors_lo[i]], m_name=grb_list[i],
                     start=starts[i], end=ends[i], difference=diffs[i], midpoints=midpoints[i], axes=axis, has_BB=has_BB[i])

# --- Inset zoom on the trailing CPL episode of 140206B ---
axins = ax[2].inset_axes([100, 900, 55, 1300], transform=ax[2].transData)
axins.errorbar(midpoints[2][-1:], ep_values[2][-1:],
               xerr=diffs[2][-1:] / 2,
               yerr=[ep_errors_lo[2][-1:], ep_errors_hi[2][-1:]],
               fmt="o", capsize=5, color='g')
axins.text(midpoints[2][-1] + 2.5, ep_values[2][-1] + 250, "CPL", fontsize=LABEL_FONT_SIZE)
axins.set_xlim(100, 155)
axins.set_ylim(5500, 8500)
axins.xaxis.set_major_locator(ticker.MultipleLocator(25))
axins.spines['top'].set_visible(False)
axins.spines['right'].set_visible(False)
axins.grid(True, which="both", alpha=GRID_ALPHA, ls=GRID_LINESTYLE)
axins.tick_params(labelsize=TICK_FONT_SIZE - 2)

[i.grid(True, which="both", alpha=GRID_ALPHA, ls=GRID_LINESTYLE) for i in ax]
[v.set_xlabel("Time [s]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i >= 2]
[v.set_ylabel("Energy [keV]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i % 2 == 0]
ax[2].set_ylim(top=3000)
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
[i.legend(loc="upper right", frameon=False, fontsize=LEGEND_FONT_SIZE) for i in ax]
plt.tight_layout()
# plt.show()
[plt.savefig(f"./peak_energy_best__all.{i}", dpi=SAVE_DPI) for i in ["png", "pdf"]]
plt.close()

######################################################################################################################
# SAVE THE VALUES
######################################################################################################################

list_of_eps = []
for grb in grb_objs:
    episode_labels = []
    for interval in grb.intervals:
        if interval.index is None:
            episode_labels.append(interval.kind.value)
        else:
            episode_labels.append(f"{interval.kind.value}{interval.index}")
    list_of_eps.append(episode_labels)

list_of_names = [[model.name for model in best] for best in grb_best]

save_value_error_as_parquet(
    list_of_ep=list_of_eps,
    grb_names=grb_list_long,
    list_of_values=list(ep_values),
    list_of_errors=(list(ep_errors_hi), list(ep_errors_lo)),
    list_of_names=list_of_names,
    filename="peak_energy.parquet",
    asym_errs=True,
)
