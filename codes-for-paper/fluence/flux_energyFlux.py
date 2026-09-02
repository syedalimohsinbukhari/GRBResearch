"""Created on May 26 20:25:58 2026"""

from grb_research import EpisodeMarkerResolver, EpisodeTypes, TimeInterval, update_style
from grb_research.grb_constants import (
    ANNOTATION_FONT_SIZE,
    CAP_SIZE,
    LEGEND_FONT_SIZE,
    LINE_WIDTH,
    MARKER_SIZE,
)
from grb_research.grb_utils import save_fig

"""
GRB Spectral Properties Visualization Script
Publication-ready figures for research paper
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# 1. LOAD DATA
# ============================================================
df = pd.read_csv("./flux_energy_flux.csv")

# ============================================================
# 2. STYLE CONFIGURATION
# ============================================================

update_style()

# MARKER_SIZE /= 1.4

t90_markers = ["o", "s", "X", "D"]

grbs = df["grb_name"].unique()

# ============================================================
# 5. FIGURE 3: Flux vs Fluence Correlation
# ============================================================
fig, ax = plt.subplots(2, 2, figsize=(10, 8), squeeze=False, sharex=True, sharey=True)
ax_flat = ax.flatten()

[i.set_xscale("log") for i in ax_flat]
[i.set_yscale("log") for i in ax_flat]

kev_to_erg = 1.602e-9
flux_ref = np.logspace(0, 3, 200)  # spans your x-axis range
for idx2, grb in enumerate(grbs):
    emr = EpisodeMarkerResolver(t90_marker=t90_markers[idx2])
    tr_count = 0
    sub = df[df["grb_name"] == grb]
    for _, row in sub.iterrows():
        if "TR" not in row["ep_type"]:
            ep_type = EpisodeTypes[row["ep_type"]]
            mm = emr.resolve(TimeInterval(ep_type))
            col_ = emr.get_color(TimeInterval(ep_type))
        else:
            ep_type = EpisodeTypes.TR
            ep_type = TimeInterval(ep_type, index=tr_count)
            mm = emr.resolve(ep_type)
            col_ = emr.get_color(ep_type)
            tr_count += 1

        # errorbar (not scatter) so the flux/fluence MC uncertainty already in the CSV is shown,
        # per CLAUDE.md's "error bars wherever the underlying value has an MC uncertainty" rule.
        ax_flat[idx2].errorbar(
            row["flux_ph_cm2_s"],
            row["fluence_erg_cm2"],
            xerr=[[row["flux_err_lower_ph_cm2_s"]], [row["flux_err_upper_ph_cm2_s"]]],
            yerr=[[row["fluence_err_lower_erg_cm2"]], [row["fluence_err_upper_erg_cm2"]]],
            marker=mm,
            markersize=MARKER_SIZE,
            color=col_,
            linestyle="none",
            capsize=CAP_SIZE,
            elinewidth=LINE_WIDTH,
            label=row["ep_type"],
            zorder=1,
        )

        # grb_name is already the paper name (e.g., GRB080916C), not the raw directory name
        # long_to_short expects -- no lookup needed.
        ax_flat[idx2].legend(loc="best", title=grb, ncol=2, fontsize=LEGEND_FONT_SIZE)

    for ls, e_kev in zip(["--", ":", "-."], [10, 100, 300]):
        e_erg = e_kev * kev_to_erg
        ax_flat[idx2].plot(flux_ref, e_erg * flux_ref, "k", ls=ls, alpha=0.25, lw=1, zorder=0)
        ax_flat[idx2].text(
            flux_ref[-1],
            e_erg * flux_ref[-1],
            f"⟨E⟩ = {e_kev / 1e3} MeV",
            fontsize=ANNOTATION_FONT_SIZE,
            color="gray",
            va="bottom",
            ha="right",
        )

[i.set_xlim(left=0.75) for i in ax_flat[2:]]
[i.set_xlabel(r"Flux (ph cm$^{-2}$ s$^{-1}$)") for i in ax_flat[2:]]
[i.set_ylabel(r"Energy flux (erg cm$^{-2}$ s$^{-1}$)") for i in ax_flat[::2]]

save_fig(fig, 'flux_vs_energy_flux')
