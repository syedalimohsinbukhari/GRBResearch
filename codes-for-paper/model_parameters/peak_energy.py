"""Created on Dec 17 13:22:15 2025"""

from typing import Tuple

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker

from grb_research import (
    GRID_ALPHA,
    GRID_LINESTYLE,
    LABEL_FONT_SIZE,
    LEGEND_FONT_SIZE,
    SAVE_DPI,
    TICK_FONT_SIZE,
    find_project_root,
    get_rng,
    plot_per_episode,
    prepare_grbs,
    seed_from_name,
    update_style,
    EpisodeMarkerResolver,
)
from grb_research.grb_constants import N_SAMPLES
from utils import convert_sbpl_to_band, extract_parameter, prepare_panel_data

update_style()

# SBPL -> Band E_peak conversion is Monte-Carlo (utils.convert_sbpl_to_band); BAND/CPL
# rows read e_peak directly from the fit and carry no MC provenance.
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)


def extract_peak_energy(
    model_collection, rng=None, n_samples: int = N_SAMPLES
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract peak energy (Ep) values and asymmetric errors from a model set.

    Returns
    -------
    values, errors_lo, errors_hi : np.ndarray
    is_mc : np.ndarray of bool
        True for rows derived via the SBPL -> Band Monte-Carlo conversion (so their
        n_samples/seed provenance is meaningful); False for BAND/CPL rows read
        directly from a fitted ``e_peak`` parameter (no MC involved).
    """
    values, errors_hi, errors_lo, is_mc = [], [], [], []

    for model in model_collection:
        if "SBPL" in model.name:
            err_lo, median, err_hi = convert_sbpl_to_band(model, n_sample=n_samples, rng=rng)
            is_mc.append(True)
        else:
            result = extract_parameter(model, "e_peak")
            median, err_hi = result if result is not None else (np.nan, np.nan)
            err_lo = err_hi
            is_mc.append(False)

        values.append(median)
        errors_hi.append(err_hi)
        errors_lo.append(err_lo)

    return np.asarray(values), np.asarray(errors_lo), np.asarray(errors_hi), np.asarray(is_mc)


grb_list = ["080916C", "131014A", "140206B", "231129C"]

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

_, _, grb_objs, grb_best = prepare_grbs(grb_list, result_file, get_best=True)
panels = prepare_panel_data(grb_objs)

ep_results = [extract_peak_energy(best, rng=rng) for best in grb_best]

_, ax = plt.subplots(2, 2, figsize=(11.5, 8.5))
ax = ax.flatten()

for i, (panel, (values, errors_lo, errors_hi, _)) in enumerate(zip(panels, ep_results)):
    # errors=[lower, upper]: matplotlib's 2-row yerr convention, matching plot_per_episode's
    # own y_low/y_high computation (errors[0]=lower magnitude, errors[1]=upper magnitude).
    # The previous version passed [hi, lo] here, silently swapping the error-bar direction
    # for every SBPL-modeled (asymmetric-error) episode -- fixed as part of this rewrite.
    plot_per_episode(
        values=values,
        errors=[errors_lo, errors_hi],
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

# --- Inset zoom on the trailing CPL episode of 140206B ---
values_140206, errors_lo_140206, errors_hi_140206, _ = ep_results[2]
panel_140206 = panels[2]
axins = ax[2].inset_axes([100, 1500, 55, 1300], transform=ax[2].transData)
axins.errorbar(
    panel_140206.midpoint[-1:],
    values_140206[-1:],
    xerr=panel_140206.diff[-1:] / 2,
    yerr=[errors_lo_140206[-1:], errors_hi_140206[-1:]],
    fmt=panel_140206.markers[-1],
    capsize=5,
    color="g",
)
axins.text(panel_140206.midpoint[-1] + 2.5, values_140206[-1] + 250, "CPL", fontsize=LABEL_FONT_SIZE)
axins.set_xlim(100, 155)
axins.set_ylim(5500, 8500)
axins.xaxis.set_major_locator(ticker.MultipleLocator(20))
axins.spines["top"].set_visible(False)
axins.spines["right"].set_visible(False)
# axins.grid(True, which="both", alpha=GRID_ALPHA, ls=GRID_LINESTYLE)
axins.tick_params(labelsize=TICK_FONT_SIZE - 2)

[v.set_xlabel("Time [s]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i >= 2]
[v.set_ylabel("Energy [keV]", fontsize=LABEL_FONT_SIZE) for i, v in enumerate(ax) if i % 2 == 0]
ax[2].set_ylim(top=3000)
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
# GRB080916C (ax[0]) and GRB140206B (ax[2]) have the most episodes -- widen their
# legend into columns so the box stays compact instead of running down the panel.
legend_ncols = {0: 2, 2: 2}
loc_ = {0: "upper center", 1: "upper right", 2: "center left", 3: "upper right"}
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
    plt.savefig(f"./peak_energy_best__all.{extension}", dpi=SAVE_DPI)
plt.close()

# ─── Save the values ──────────────────────────────────────────────────────

frame = pd.DataFrame(
    [
        {
            "grb_name": f"GRB{short_name}",
            "episode": panel.episode_labels[i],
            "model_name": panel.model_names[i],
            "e_peak_keV": values[i],
            "e_peak_err_lower_keV": errors_lo[i],
            "e_peak_err_upper_keV": errors_hi[i],
            "n_samples": N_SAMPLES if is_mc[i] else np.nan,
            "seed": SEED if is_mc[i] else np.nan,
        }
        for short_name, panel, (values, errors_lo, errors_hi, is_mc) in zip(grb_list, panels, ep_results)
        for i in range(len(panel.episode_labels))
    ]
)
frame.to_csv("peak_energy.csv", index=False)
print(f"Saved: peak_energy.csv  ({len(frame)} rows)")
