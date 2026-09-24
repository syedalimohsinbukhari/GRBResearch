"""Plotting for the joint N-pulse fit: data + total joint model + each pulse's own component,
overlaid -- what plot_diagnostics.plot_fit_overlay (one pulse at a time) can't show. Same
generated-output convention as the rest of this project (update_style() + save_fig(), CSV via
pandas) -- see ../plot_diagnostics.py's module docstring for the convention itself.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[4]  # GRB231129C_joint -> new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import update_style  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

from joint_pulse3 import unflatten_params  # noqa: E402

update_style()

_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:cyan", "tab:olive", "tab:pink"]


def plot_joint_overlay(t_window, y_window, sigma_window, model, popt, out_dir, label, pulse_indices=None):
    """pulse_indices: the REAL (1-indexed) pulse number for each position in popt, e.g. [1,2,3,5,6]
    when pulse 4 has been dropped from the model -- without this, legend labels default to
    positional 1..N, which silently mislabels every pulse after a gap (a real bug caught once:
    the first version of this function did exactly that on the 5-pulse variant)."""
    n = model.n_pulses
    if pulse_indices is None:
        pulse_indices = list(range(1, n + 1))
    params = unflatten_params(popt, n)
    y_total = model(t_window, *popt)

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.errorbar(t_window, y_window, yerr=sigma_window, fmt=".", ms=2, alpha=0.25, color="tab:gray", label="data", zorder=1)
    for i, (amplitude, t_peak, t_v) in enumerate(params):
        y_i = model.pulse_model(i, t_window, amplitude, t_peak, t_v)
        ax.plot(t_window, y_i, "--", lw=1.0, color=_COLORS[i % len(_COLORS)], label=f"pulse {pulse_indices[i]}", zorder=2, alpha=0.8)
    ax.plot(t_window, y_total, "-", lw=1.6, color="tab:red", label="joint total", zorder=3)

    ax.set_xlabel("t [s]")
    ax.set_ylabel("rate")
    ax.set_title(f"{label}: joint {n}-pulse 3-param fit")
    ax.legend(fontsize=7, ncol=2)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_stem = out_dir / f"joint_overlay_{label}"
    save_fig(fig, fig_stem)

    df = pd.DataFrame({"t": t_window, "y_data": y_window, "sigma": sigma_window, "y_total_fit": y_total})
    for i, (amplitude, t_peak, t_v) in enumerate(params):
        df[f"y_pulse{pulse_indices[i]}"] = model.pulse_model(i, t_window, amplitude, t_peak, t_v)
    csv_path = out_dir / f"joint_overlay_{label}.csv"
    df.to_csv(csv_path, index=False)

    return {"csv": csv_path, "pdf": Path(f"{fig_stem}.pdf"), "png": Path(f"{fig_stem}.png")}
