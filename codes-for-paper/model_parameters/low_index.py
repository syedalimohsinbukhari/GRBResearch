"""Created on Dec 17 13:22:15 2025"""

import json
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from grb_research import update_style
from grb_research.grb_core import prepare_grbs
from utils import (
    extract_parameter,
    find_project_root,
    short_to_long,
    LABEL_FONT_SIZE,
    TICK_FONT_SIZE,
    LEGEND_FONT_SIZE,
    GRBCatalog,
    ModelSet,
    plot_per_episode,
    save_value_error_as_parquet,
)

update_style()


def extract_low_index(best_model: ModelSet) -> Tuple[np.ndarray, np.ndarray]:
    """Extract the low (alpha) spectral index values and errors from a model set."""
    value, error = [], []
    for model in best_model:
        result = extract_parameter(model, "index1")
        if result is not None:
            value.append(result[0])
            error.append(result[1])

    return np.array(value), np.array(error)


SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

with open(result_file, "r") as f:
    example_data = json.load(f)

grb_list = ["080916C", "131014A", "140206B", "231129C"]
grb_list_long = [short_to_long[i] for i in grb_list]

_, _, grb_objs, _ = prepare_grbs(grb_list, result_file, get_best=True)

gc = GRBCatalog.from_iterable(grb_list=grb_list, data=example_data, name_mapping=short_to_long)

grb080916c = gc.get_grb(grb_list_long[0])
grb131014a = gc.get_grb(grb_list_long[1])
grb140206b = gc.get_grb(grb_list_long[2])
grb231129c = gc.get_grb(grb_list_long[3])

grb080916c_best = grb080916c.get_all_best_models()
grb131014a_best = grb131014a.get_all_best_models()
grb140206b_best = grb140206b.get_all_best_models()
grb231129c_best = grb231129c.get_all_best_models()

start_080916, end_080916, diff_080916, midpoint_080916 = grb080916c.intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_131014, end_131014, diff_131014, midpoint_131014 = grb131014a.intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_140206, end_140206, diff_140206, midpoint_140206 = grb140206b.intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)
start_231129, end_231129, diff_231129, midpoint_231129 = grb231129c.intervals.extract_interval_arrays(
    return_include=("diff", "midpoint")
)

alpha_value_080916c, alpha_error_080916c = extract_low_index(grb080916c_best)
alpha_value_131014a, alpha_error_131014a = extract_low_index(grb131014a_best)
alpha_value_140206b, alpha_error_140206b = extract_low_index(grb140206b_best)
alpha_value_231129c, alpha_error_231129c = extract_low_index(grb231129c_best)

_, ax = plt.subplots(4, 1, figsize=(5.5, 12))

plot_per_episode(values=alpha_value_080916c, errors=alpha_error_080916c, m_name=grb_list[0], start=start_080916,
                 end=end_080916, difference=diff_080916, midpoints=midpoint_080916, axes=ax[0])

plot_per_episode(values=alpha_value_131014a, errors=alpha_error_131014a, m_name=grb_list[1], start=start_131014,
                 end=end_131014, difference=diff_131014, midpoints=midpoint_131014, axes=ax[1])

plot_per_episode(values=alpha_value_140206b, errors=alpha_error_140206b, m_name=grb_list[2], start=start_140206,
                 end=end_140206, difference=diff_140206, midpoints=midpoint_140206, axes=ax[2])

plot_per_episode(values=alpha_value_231129c, errors=alpha_error_231129c, m_name=grb_list[3], start=start_231129,
                 end=end_231129, difference=diff_231129, midpoints=midpoint_231129, axes=ax[3])

[i.grid(True, which="both", alpha=0.5, ls="--") for i in ax]
[i.set_xlabel("Time [s]", fontsize=LABEL_FONT_SIZE) for i in ax]
[i.set_ylabel(r"Lower Index [$\alpha$]", fontsize=LABEL_FONT_SIZE) for i in ax]
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
[i.legend(loc="best", frameon=False, fontsize=LEGEND_FONT_SIZE) for i in ax]
plt.tight_layout()
# plt.show()
[plt.savefig(f"./low_index_best_all.{i}", dpi=600) for i in ["png", "pdf"]]
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

list_of_values = [alpha_value_080916c, alpha_value_131014a, alpha_value_140206b, alpha_value_231129c]
list_of_errors = [alpha_error_080916c, alpha_error_131014a, alpha_error_140206b, alpha_error_231129c]
list_of_names = [[i.name for i in j] for j in [grb080916c_best, grb131014a_best, grb140206b_best, grb231129c_best]]

save_value_error_as_parquet(
    list_of_ep=list_of_eps,
    grb_names=grb_list_long,
    list_of_values=list_of_values,
    list_of_errors=list_of_errors,
    list_of_names=list_of_names,
    filename="low_index.parquet",
)
