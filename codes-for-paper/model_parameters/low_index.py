"""Created on Dec 17 13:22:15 2025"""

from typing import Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from grb_research import (
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    SAVE_DPI,
    TICK_FONT_SIZE,
    find_project_root,
    plot_per_episode,
    prepare_grbs,
    update_style,
)
from utils import extract_parameter, prepare_panel_data

update_style()


def extract_low_index(model_collection) -> Tuple[np.ndarray, np.ndarray]:
    """Extract the low (alpha) spectral index values and errors from a model set."""
    values, errors = [], []
    for model in model_collection:
        result = extract_parameter(model, "index1")
        if result is not None:
            values.append(result[0])
            errors.append(result[1])

    return np.array(values), np.array(errors)


grb_list = ["080916C", "131014A", "140206B", "231129C"]

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

_, _, grb_objs, grb_best = prepare_grbs(grb_list, result_file, get_best=True)
panels = prepare_panel_data(grb_objs)

alpha_results = [extract_low_index(best) for best in grb_best]

_, ax = plt.subplots(2, 2, figsize=(11.5, 8.5))
ax = ax.flatten()

for i, (panel, (values, errors)) in enumerate(zip(panels, alpha_results)):
    plot_per_episode(
        values=values,
        errors=errors,
        m_name=grb_list[i],
        start=panel.start,
        end=panel.end,
        difference=panel.diff,
        midpoints=panel.midpoint,
        axes=ax[i],
        has_BB=panel.has_bb,
        episode_labels=panel.episode_labels,
        model_names=panel.model_names,
        markers=panel.markers,
    )

[v.set_xlabel("Time [s]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i >= 2]
[v.set_ylabel("Energy [keV]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i % 2 == 0]
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
# GRB080916C (ax[0]) and GRB140206B (ax[2]) have the most episodes -- widen their
# legend into columns so the box stays compact instead of running down the panel.
# GRB140206B's points span its full y-range (TR2 near the top, EX0 near the bottom),
# so no horizontal band inside the axes is free of data for a near-full-width 3-column
# box -- reserve headroom below the data instead of fighting for a spot inside it.
ylim = ax[2].get_ylim()
ax[2].set_ylim(ylim[0] - 0.3 * (ylim[1] - ylim[0]), ylim[1])

legend_ncols = {0: 2, 2: 1}
loc_ = {0: "upper center", 1: "upper right", 2: "lower left", 3: "upper right"}
[
    a.legend(
        loc=loc_.get(i, None),
        fontsize=LEGEND_FONT_SIZE,
        title=f"GRB{grb_list[i]}",
        title_fontsize=LEGEND_FONT_SIZE,
        ncols=legend_ncols.get(i, 1),
    )
    for i, a in enumerate(ax)
]
plt.tight_layout()
for extension in ("png", "pdf"):
    plt.savefig(f"./low_index_best__all.{extension}", dpi=SAVE_DPI)
plt.close()

# ─── Save the values ──────────────────────────────────────────────────────

frame = pd.DataFrame(
    [
        {
            "grb_name": f"GRB{short_name}",
            "episode": panel.episode_labels[i],
            "model_name": panel.model_names[i],
            "alpha_index": values[i],
            "alpha_index_err": errors[i],
        }
        for short_name, panel, (values, errors) in zip(grb_list, panels, alpha_results)
        for i in range(len(panel.episode_labels))
    ]
)
frame.to_csv("low_index.csv", index=False)
print(f"Saved: low_index.csv  ({len(frame)} rows)")
