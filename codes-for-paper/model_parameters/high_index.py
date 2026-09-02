"""Created on Dec 17 13:22:15 2025"""

from typing import Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from plotez.typing import ArrayLike

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


def extract_high_index(model_collection) -> Tuple[ArrayLike, ArrayLike]:
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

_, _, grb_objs, grb_best = prepare_grbs(grb_list, result_file, get_best=True)
panels = prepare_panel_data(grb_objs)

beta_results = [extract_high_index(best) for best in grb_best]

_, ax = plt.subplots(2, 2, figsize=(11.5, 8.5))
ax = ax.flatten()

for i, (panel, (values, errors)) in enumerate(zip(panels, beta_results)):
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
# GRB140206B's wide box collides with its EX0/TR1 points under "lower right"; its data
# clusters at low x/high y, so "upper right" is the clear corner instead.
legend_ncols = {0: 2, 2: 2}
[
    a.legend(
        loc="lower right",
        fontsize=LEGEND_FONT_SIZE,
        title=f"GRB{grb_list[i]}",
        title_fontsize=LEGEND_FONT_SIZE,
        ncols=legend_ncols.get(i, 1),
    )
    for i, a in enumerate(ax)
]
plt.tight_layout()
for extension in ("png", "pdf"):
    plt.savefig(f"./high_index_best__all.{extension}", dpi=SAVE_DPI)
plt.close()

# ─── Save the values ──────────────────────────────────────────────────────

frame = pd.DataFrame(
    [
        {
            "grb_name": f"GRB{short_name}",
            "episode": panel.episode_labels[i],
            "model_name": panel.model_names[i],
            "beta_index": values[i],
            "beta_index_err": errors[i],
        }
        for short_name, panel, (values, errors) in zip(grb_list, panels, beta_results)
        for i in range(len(panel.episode_labels))
    ]
)
frame.to_csv("high_index.csv", index=False)
print(f"Saved: high_index.csv  ({len(frame)} rows)")
