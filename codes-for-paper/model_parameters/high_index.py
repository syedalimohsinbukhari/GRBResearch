"""Created on Dec 17 13:22:15 2025"""

from typing import Tuple

import numpy as np
from matplotlib import pyplot as plt

from grb_research import (update_style, plot_per_episode, LABEL_FONT_SIZE, TICK_FONT_SIZE, LEGEND_FONT_SIZE,
                          SAVE_DPI, save_value_error_as_parquet)
from utils import (
    extract_parameter,
    find_project_root,
    prepare_grbs,
    ModelSet,
)

update_style()


def extract_high_index(model_collection: ModelSet) -> Tuple[np.typing.ArrayLike, np.typing.ArrayLike]:
    """
    Extract high (beta) spectral index values and errors based on model type.

    Rules:
    - BAND / SBPL and their derivatives:
        extract parameter where 'index2' is in the name
    - CPL_PL / CPL_PL_BB:
        extract parameter where 'add_index_pl' is in the name
    - CPL only:
        append np.nan (no high-energy index)
    """

    values = []
    errors = []

    sbpl_band_models = {"band", "band_bb", "band_pl_bb", "sbpl", "sbpl_bb", "sbpl_pl_bb"}
    cpl_pl_models = {"cpl_pl_bb"}

    for model in model_collection:
        m_name = model.name.lower()

        # --- BAND / SBPL family ---
        if m_name in sbpl_band_models:
            result = extract_parameter(model=model, param_pattern="index2")
            if result is not None:
                values.append(result[0])
                errors.append(result[1])

        # --- CPL + PL family ---
        elif m_name in cpl_pl_models:
            result = extract_parameter(model=model, param_pattern="add_index_pl")
            if result is not None:
                values.append(result[0])
                errors.append(result[1])

        # --- Pure CPL ---
        else:
            values.append(np.nan)
            errors.append(np.nan)

    return np.asarray(values), np.asarray(errors)


grb_list = ["080916C", "131014A", "140206B", "231129C"]

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

_, grb_list_long, grb_objs, grb_best = prepare_grbs(grb_list, result_file, get_best=True)

has_BB = [['BB' in i.name for i in t_.get_all_best_models()] for t_ in grb_objs]

start_080916, end_080916, diff_080916, midpoint_080916 = grb_objs[0].intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_131014, end_131014, diff_131014, midpoint_131014 = grb_objs[1].intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_140206, end_140206, diff_140206, midpoint_140206 = grb_objs[2].intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_231129, end_231129, diff_231129, midpoint_231129 = grb_objs[3].intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)

start_140206 = start_140206[:-1]
end_140206 = end_140206[:-1]
diff_140206 = diff_140206[:-1]
midpoint_140206 = midpoint_140206[:-1]

beta_value_080916c, beta_error_080916c = extract_high_index(grb_best[0])
beta_value_131014a, beta_error_131014a = extract_high_index(grb_best[1])
beta_value_140206b, beta_error_140206b = extract_high_index(grb_best[2])
beta_value_231129c, beta_error_231129c = extract_high_index(grb_best[3])

_, ax = plt.subplots(2, 2, figsize=(10, 7.5))
ax = ax.flatten()

plot_per_episode(values=beta_value_080916c, errors=beta_error_080916c, m_name=grb_list[0], start=start_080916,
                 end=end_080916, difference=diff_080916, midpoints=midpoint_080916, axes=ax[0], has_BB=has_BB[0])

plot_per_episode(values=beta_value_131014a, errors=beta_error_131014a, m_name=grb_list[1], start=start_131014,
                 end=end_131014, difference=diff_131014, midpoints=midpoint_131014, axes=ax[1], has_BB=has_BB[1])

plot_per_episode(values=beta_value_140206b, errors=beta_error_140206b, m_name=grb_list[2], start=start_140206,
                 end=end_140206, difference=diff_140206, midpoints=midpoint_140206, axes=ax[2], has_BB=has_BB[2])

plot_per_episode(values=beta_value_231129c, errors=beta_error_231129c, m_name=grb_list[3], start=start_231129,
                 end=end_231129, difference=diff_231129, midpoints=midpoint_231129, axes=ax[3], has_BB=has_BB[3])

# [i.grid(True, which="both", alpha=0.5, ls="--") for i in ax]
[v.set_xlabel("Time [s]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i >= 2]
[v.set_ylabel("Energy [keV]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i % 2 == 0]
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
[i.legend(loc="lower right", frameon=False, fontsize=LEGEND_FONT_SIZE) for i in ax]
plt.tight_layout()
# plt.show()
[plt.savefig(f"./high_index_best__all.{i}", dpi=SAVE_DPI) for i in ["png", "pdf"]]
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

list_of_values = [beta_value_080916c, beta_value_131014a, beta_value_140206b, beta_value_231129c]
list_of_errors = [beta_error_080916c, beta_error_131014a, beta_error_140206b, beta_error_231129c]
list_of_names = [[i.name for i in j] for j in grb_best]

save_value_error_as_parquet(
    list_of_ep=list_of_eps,
    grb_names=grb_list_long,
    list_of_values=list_of_values,
    list_of_errors=list_of_errors,
    list_of_names=list_of_names,
    filename="high_index.parquet",
)
