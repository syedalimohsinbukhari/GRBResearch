"""Reporting for sections 4, 7 (T3/T4/T5), and 8 diagnostics: CSV + PDF/PNG outputs, following
this project's generated-output convention (same pattern as plot_profile_scan.py, which covers
section 5 -- see that module's docstring). Split out so the tests that exercise sections 4/7/8
can each save their own artifacts without duplicating the matplotlib/grb_research boilerplate.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[3]  # new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import update_style  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

update_style()


def plot_fit_overlay(t_window, y_window, y_fit, out_dir, label, sigma=None, half_max=None, extra_title=""):
    """Data + best-fit pulse3 overlay -- the standard "does the fit look right" check. Used by
    section 4 (finalized pulse) and section 7 T3 (with half_max=(x_left, x_right, half_value)
    marking the brentq crossing points)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    if sigma is not None:
        ax.errorbar(t_window, y_window, yerr=sigma, fmt=".", ms=2, alpha=0.3, color="tab:gray", label="data", zorder=1)
    else:
        ax.plot(t_window, y_window, ".", ms=2, alpha=0.3, color="tab:gray", label="data", zorder=1)
    ax.plot(t_window, y_fit, "-", color="tab:red", lw=1.2, label="3-param fit", zorder=2)
    if half_max is not None:
        x_left, x_right, half_value = half_max
        ax.axhline(half_value, color="tab:green", ls="--", lw=0.8, label="half max")
        ax.axvline(x_left, color="tab:green", ls=":", lw=0.8)
        ax.axvline(x_right, color="tab:green", ls=":", lw=0.8)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("rate")
    ax.set_title(f"{label}{(': ' + extra_title) if extra_title else ''}")
    ax.legend(fontsize=7)

    fig_stem = out_dir / f"fit_overlay_{label}"
    save_fig(fig, fig_stem)

    df = pd.DataFrame({"t": t_window, "y_data": y_window, "y_fit": y_fit, "residual": np.asarray(y_window) - np.asarray(y_fit)})
    if sigma is not None:
        df["sigma"] = sigma
    csv_path = out_dir / f"fit_overlay_{label}.csv"
    df.to_csv(csv_path, index=False)

    return {"csv": csv_path, "pdf": Path(f"{fig_stem}.pdf"), "png": Path(f"{fig_stem}.png")}


def plot_t4_scatter(results, r_true, out_dir, label):
    """r_seed vs. converged r, and converged t_peak vs. r_seed -- visualizes section 7 T4's
    reproducibility/runaway finding (box-edge pin, or a wrong-but-stable interior minimum)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    r_seed = np.array([r["r_seed"] for r in results])
    r_fit = np.array([r["r"] for r in results])
    t_peak_fit = np.array([r["t_peak"] for r in results])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))
    ax1.plot(r_seed, r_fit, "o-", color="tab:blue")
    ax1.axhline(r_true, color="tab:green", ls="--", lw=0.8, label=f"r_true={r_true:.0f}")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("seed r")
    ax1.set_ylabel("fitted r")
    ax1.set_title("4-param fit: seed r -> converged r")
    ax1.legend(fontsize=7)

    # Same fix as plot_window_widening's y-axis: when t_peak is genuinely stable across seeds
    # (to ~1e-5 s, as expected once r has run away/converged -- see PROGRESS.md), matplotlib's
    # autoscale-plus-offset-notation zooms into that noise floor and makes it look like a real
    # trend (e.g. "1e-5+1.49e1"). Pin the range to a fixed absolute band around the mean instead,
    # so a flat line reads as flat and only a real, tolerance-breaking deviation shows visibly.
    t_peak_mean = float(np.mean(t_peak_fit))
    t_peak_band = max(0.05, 0.01 * abs(t_peak_mean))
    ax2.plot(r_seed, t_peak_fit, "o-", color="tab:red")
    ax2.axhspan(t_peak_mean - t_peak_band, t_peak_mean + t_peak_band, color="tab:red", alpha=0.08)
    ax2.set_ylim(t_peak_mean - t_peak_band, t_peak_mean + t_peak_band)
    ax2.ticklabel_format(useOffset=False, style="plain", axis="y")
    ax2.set_xscale("log")
    ax2.set_xlabel("seed r")
    ax2.set_ylabel("fitted t_peak [s]")
    ax2.set_title(f"t_peak stability across seeds (+/-{t_peak_band:.3g} s band)")

    fig.suptitle(label)
    fig_stem = out_dir / f"t4_scatter_{label}"
    save_fig(fig, fig_stem)

    df = pd.DataFrame({"r_seed": r_seed, "r_fit": r_fit, "t_peak_fit": t_peak_fit})
    csv_path = out_dir / f"t4_scatter_{label}.csv"
    df.to_csv(csv_path, index=False)

    return {"csv": csv_path, "pdf": Path(f"{fig_stem}.pdf"), "png": Path(f"{fig_stem}.png")}


def plot_window_widening(rows, out_dir, label, rel_tol=0.05):
    """t_peak/t_v vs. window-widening factor -- section 8's window-widening invariance check.
    rows: iterable of (name, t_min, t_max, t_peak, t_v).

    y-limits are pinned to +/- rel_tol around the base (first row) value, NOT autoscaled: with
    near-machine-precision invariance (as expected for a well-posed fit -- see PROGRESS.md),
    matplotlib's default autoscale-plus-offset-notation zooms into the ~1e-9 s noise floor and
    makes it visually look like a large swing. Pinning the range to the actual invariance
    tolerance this check enforces shows a flat line when the check truly passes, and only shows
    a visible trend when the deviation is large enough to actually matter.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(rows, columns=["name", "t_min", "t_max", "t_peak", "t_v"])
    x = np.arange(len(df))
    t_peak_base, t_v_base = df["t_peak"].iloc[0], df["t_v"].iloc[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))
    ax1.plot(x, df["t_peak"], "o-", color="tab:blue")
    ax1.axhspan(t_peak_base * (1 - rel_tol), t_peak_base * (1 + rel_tol), color="tab:blue", alpha=0.08)
    ax1.set_ylim(t_peak_base * (1 - rel_tol), t_peak_base * (1 + rel_tol))
    ax1.ticklabel_format(useOffset=False, style="plain", axis="y")
    ax1.set_xticks(x)
    ax1.set_xticklabels(df["name"])
    ax1.set_ylabel("t_peak [s]")
    ax1.set_title(f"t_peak vs. window (+/-{rel_tol:.0%} band)")

    ax2.plot(x, df["t_v"], "o-", color="tab:red")
    ax2.axhspan(t_v_base * (1 - rel_tol), t_v_base * (1 + rel_tol), color="tab:red", alpha=0.08)
    ax2.set_ylim(t_v_base * (1 - rel_tol), t_v_base * (1 + rel_tol))
    ax2.ticklabel_format(useOffset=False, style="plain", axis="y")
    ax2.set_xticks(x)
    ax2.set_xticklabels(df["name"])
    ax2.set_ylabel("t_v [s]")
    ax2.set_title(f"t_v vs. window (+/-{rel_tol:.0%} band)")

    fig.suptitle(label)
    fig_stem = out_dir / f"window_widening_{label}"
    save_fig(fig, fig_stem)

    csv_path = out_dir / f"window_widening_{label}.csv"
    df.to_csv(csv_path, index=False)

    return {"csv": csv_path, "pdf": Path(f"{fig_stem}.pdf"), "png": Path(f"{fig_stem}.png")}
