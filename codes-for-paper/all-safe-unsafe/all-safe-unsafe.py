"""Created on Feb 09 16:18:36 2026"""

import json

import matplotlib.pyplot as plt

from grb_research import (
    GRBCatalog,
    GRBPlotStyle as grbStyle,
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
    MODEL_ORDER,
    TICK_FONT_SIZE,
    find_project_root,
    short_to_long,
    update_style,
)

update_style()

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"
grb_name = ["080916C", "131014A", "140206B", "231129C"]

with open(result_file, "r") as f:
    data = json.load(f)

grb = GRBCatalog.from_iterable(grb_list=grb_name, data=data, name_mapping=short_to_long)

f, ax = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
ax = ax.flatten()

colors = grbStyle.GRB_COLORS_ITERABLE
MODEL_ORDER_R = [i.replace("_", "+") for i in MODEL_ORDER]

for index, value in enumerate(grb.grb_list):
    dict_ = value.get_model_count(separate=True)
    safe = dict_["SAFE"]
    unsafe = dict_["UNSAFE"]
    marginal = dict_["MARGINAL"]
    best_ = dict_["BEST"]

    # Create data for the bar chart
    safe_counts = [safe.get(model, 0) for model in MODEL_ORDER]
    unsafe_counts = [unsafe.get(model, 0) for model in MODEL_ORDER]
    marginal_counts = [marginal.get(model, 0) for model in MODEL_ORDER]
    best_counts = [best_.get(model, 0) for model in MODEL_ORDER]

    # X positions for bars
    x_positions = list(range(len(MODEL_ORDER)))

    # Plot bars with different hatches and colors
    # Plot UNSAFE data first (keep the alpha value low)
    ax[index].bar(
        x_positions,
        unsafe_counts,
        label="UNSAFE",
        color=colors[index],
        alpha=0.25,
        edgecolor="black",
        linewidth=0.5,
        hatch="///",
    )

    # Plot MARGINAL data on top (with hatch)
    ax[index].bar(
        x_positions,
        marginal_counts,
        label="MARGINAL",
        color=colors[index],
        hatch="O",
        alpha=0.5,
        edgecolor="black",
        linewidth=0.5,
        bottom=unsafe_counts,
    )

    # Plot SAFE data on top (with hatch)
    ax[index].bar(
        x_positions,
        safe_counts,
        label="SAFE",
        color=colors[index],
        alpha=0.75,
        hatch="\\",
        edgecolor="black",
        linewidth=0.5,
        bottom=[a + b for a, b in zip(unsafe_counts, marginal_counts)],
    )

    # Plot the BEST data on top (without hatch style)
    ax[index].bar(
        x_positions,
        best_counts,
        label="BEST",
        color=colors[index],
        edgecolor="black",
        linewidth=0.5,
        bottom=[a + b + c for a, b, c in zip(safe_counts, marginal_counts, unsafe_counts)],
    )

    # Set tick positions and labels
    ax[index].set_xticks(x_positions)
    ax[index].set_xticklabels(MODEL_ORDER_R)
    ax[index].legend(
        loc="best",
        title=f"GRB{grb_name[index]}",
        ncol=2,
        shadow=True,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_FONT_SIZE,
    )

[v.tick_params("x", rotation=90, labelsize=TICK_FONT_SIZE) for v in ax]
[v.tick_params("y", labelsize=TICK_FONT_SIZE) for v in ax]
ax[0].set_ylim(top=11.5)
[i.set_ylabel("No. of models", fontsize=LABEL_FONT_SIZE) for i in [ax[0], ax[2]]]
plt.tight_layout()
# plt.show()
[plt.savefig(f"./all-safe-unsafe.{i}", dpi=600) for i in ["png", "pdf"]]
plt.close()
