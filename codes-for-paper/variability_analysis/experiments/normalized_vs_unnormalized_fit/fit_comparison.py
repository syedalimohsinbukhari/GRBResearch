"""Experiment: does peak-normalizing the light curve before fitting (this project's existing
convention, `y /= Y_MAX_CTS_PER_S`, used by every fitter_*.py in the parent folder) actually change
the fitted Norris-pulse decomposition, or is it just a numerical convenience?

Two methods, same data, same P0 (up to a unit conversion), same bounds:

1. NORMALIZED (existing production convention) -- fit `NorrisFitter` (../../norris_fit.py, wraps
   pymultifit's `BaseFitter`, which calls plain `scipy.optimize.curve_fit(..., bounds=...)` with no
   scaling option) against y peak-normalized to 1.
2. UNNORMALIZED, explicit x_scale -- fit the same pulses directly against scipy.optimize.least_squares
   on the RAW (not normalized) counts/s data, bypassing pymultifit entirely (its BaseFitter.fit() does
   not expose x_scale/method, so this cannot be done through NorrisFitter). The optimizer is told each
   parameter's natural scale explicitly: amplitude ~ the data's own peak (Y_MAX_CTS_PER_S), t_s/tau1
   ~5s, tau2 ~1s -- the same order-of-magnitude split that peak-normalization encodes implicitly.

`x_scale='jac'` (scipy's automatic per-parameter scaling) was tried first, for GRB231129C, and
rejected: it does NOT recover a good fit -- SSE ~1.03e8 vs the reference ~5.71e7, i.e. it lands in the
same bad basin as no scaling at all (x_scale=1.0). Automatic scaling from the Jacobian apparently isn't
enough here; only an explicit, domain-informed x_scale array works. That's why this script never calls
'jac' and goes straight to the explicit array.

This script does NOT modify anything outside this folder. It imports light_curves.py and norris_fit.py
from the parent folder read-only (sys.path insert, not a package install) and reads each burst's P0 /
fit window from its own fitter_*.py file's *committed* values (copied here as literals, not imported --
this script must not import or execute any fitter_*.py, since those are separate scripts with side
effects of their own).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
VARIABILITY_DIR = HERE.parent.parent  # codes-for-paper/variability_analysis/
PROJECT_ROOT = VARIABILITY_DIR.parent.parent  # GRBResearchWork/
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from light_curves import lightcurve_data  # noqa: E402
from norris_fit import NorrisFitter, norris_pulse, t_peak, tv_value  # noqa: E402
from grb_research import update_style, LINE_WIDTH  # noqa: E402
from grb_research.grb_utils import save_fig  # noqa: E402

LC_DIR = PROJECT_ROOT / "light_curves"
ENERGY_LOW, ENERGY_HIGH = 10, 400

# Amplitude left free (scaled by each burst's own Y_MAX at run time); t_s/tau1 given a slightly wider scale than tau2
# since onset times and rise timescales tend to range wider than decay timescales in these fits -- the exact split
# doesn't need to be precise, x_scale only needs to be the right *order of magnitude*, unlike a seed value which needs
# to be close to the true optimum.
X_SCALE_TS_TAU1 = 5.0
X_SCALE_TAU2 = 1.0

# Per-burst config, copied (not imported) from each burst's own fitter_*.py as of 2026-09-23 -- read
# only, this script does not modify or execute those files.
BURSTS = {
    "GRB080916009": dict(
        paper_name="GRB080916C",
        # Widened to the light curve's own full x.min()/x.max() range, for now (2026-09-24) --
        # summed_nai_curve()'s mask (t > window[0]) & (t < window[1]) with (-inf, inf) bounds keeps every
        # point without needing to look up t.min()/t.max() first. Was (-1, 70), the production window.
        window=(-np.inf, np.inf),
        # 6-pulse P0 (was 7) -- updated 2026-09-23 to match the production fix in ../../fitter.py: the
        # dropped pulse (was A=0.1, t_s=20, tau1=9, tau2=1) converged to a degenerate near-delta-function
        # spike (tau2~0.009) on this burst's summed-detector curve in BOTH methods tested here (see the
        # now-stale "GRB080916C, pulse 5" divergence this experiment originally flagged in comparison.md --
        # that finding is what confirmed the pulse was unreliable, not a normalization artifact, and is why
        # the user asked for it to be dropped here too, for consistency with production).
        p0=[
            (0.4, -0.7, 2.34, 0.471),
            (0.7, 0.3, 0.88, 6),
            (0.2, 5.3, 0.6, 6),
            (0.3, 1.3, 43, 14),
            (0.3, 52, 18, 0.7),
            (0.3, 61, 0.3, 3),
        ],
        max_nfev=20000,
    ),
    # "GRB131014215": dict(
    #     paper_name="GRB131014A",
    #     window=(-1, 10),
    #     p0=[
    #         (0.115, -0.9, 1, 1),
    #         (0.2, -0.8, 1, 1),
    #         (0.9, 1.19, 1, 1),
    #         (0.4, 2.4, 1, 1),
    #         (0.3, 2.71, 1, 1),
    #     ],
    #     # Known issue (see ../../variability_analysis.md, BUG-23 follow-up): this burst's
    #     # summed-3-detector curve needs a larger iteration budget than the NorrisFitter/pymultifit
    #     # default (5000) to converge at all -- confirmed independently in the parent session. Applied
    #     # here to both methods so a convergence failure isn't mistaken for a real method difference.
    #     max_nfev=20000,
    # ),
    # "GRB140206275": dict(
    #     paper_name="GRB140206B",
    #     # Widened from (-1, 160) to (-20, 300), 2026-09-23 (user call): pulse 7 (the broad, low-amplitude
    #     # "pedestal" pulse, tau1~2000) was landing with t_s pinned against the narrow window's own lower
    #     # bound (t_min=-1) in the unnormalized fit -- NorrisFitter.fit_boundaries() sets t_s's lower bound
    #     # to the fit window's own left edge, so a pinned t_s means the box constraint, not the data, was
    #     # deciding that parameter. Matches the "wide_-20_300" window already used as this burst's own
    #     # robustness check in ../window_sensitivity_GRB140206275/window_sensitivity.py.
    #     window=(-20, 300),
    #     # COMPLEX_P0 from fitter_GRB140206275.py -- the user-preferred, finer-grained decomposition
    #     # used for that burst's main per-episode results (see that file's own comments).
    #     p0=[
    #         (0.13, -0.3, 0.12, 1.68),
    #         (0.35, 4.0, 6, 14),
    #         (0.5, 11, 5, 1.4),
    #         (0.5, 28, 2, 1),
    #         (0.23, 24, 0.3, 1.434),
    #         (0.23, 23, 83, 1.1),
    #         (0.05, -0.9, 3400, 4.4),
    #     ],
    #     max_nfev=20000,
    # ),
    # "GRB231129779": dict(
    #     paper_name="GRB231129C",
    #     window=(-1, 10),
    #     p0=[
    #         (0.6, -0.2, 1, 1),
    #         (0.3, 0.08, 1, 1),
    #         (0.6, 2, 0.5, 0.5),
    #         (0.6, 4, 0.5, 0.5),
    #         (0.3, 4.2, 2, 1),
    #     ],
    #     max_nfev=20000,
    # ),
}


def summed_nai_curve(grb_dir: Path, window):
    dat_nai = sorted(f.stem for f in grb_dir.glob("*.dat") if "n" in f.stem)
    nai_data = [lightcurve_data(str(grb_dir / f"{d}.dat"), ENERGY_LOW, ENERGY_HIGH) for d in dat_nai]
    t = nai_data[0][0]
    for det, (t_i, _, _) in zip(dat_nai[1:], nai_data[1:]):
        assert np.array_equal(t, t_i), f"{det}'s time grid differs from {dat_nai[0]}'s"
    r = np.sum([r for _, r, _ in nai_data], axis=0)
    b = np.sum([b for _, _, b in nai_data], axis=0)
    mask = np.logical_and(t > window[0], t < window[1])
    return t[mask], (r - b)[mask], dat_nai


def total_model(t, params, n_pulses):
    return sum(norris_pulse(t, params[i * 4:(i + 1) * 4]) for i in range(n_pulses))


def fit_normalized(t, y_raw, y_max, p0, max_nfev):
    """Method 1: existing production convention -- peak-normalize, fit with NorrisFitter."""
    y_norm = y_raw / y_max
    nf = NorrisFitter(t, y_norm, max_iterations=max_nfev)
    print(p0)
    nf.fit(p0=p0)
    return nf.params  # normalized-scale params; amplitude must be rescaled by y_max by the caller


def fit_unnormalized_xscale(t, y_raw, y_max, p0, n_pulses, max_nfev):
    """Method 2: raw data, explicit x_scale telling the optimizer each parameter's natural
    magnitude (amplitude ~ y_max, t_s/tau1 ~5, tau2 ~1) instead of transforming the data."""
    p0_raw = np.array([(a * y_max, ts, tau1, tau2) for (a, ts, tau1, tau2) in p0]).flatten()
    lb = np.tile([0.0, t.min(), 1e-4, 1e-4], n_pulses)
    ub = np.tile([np.inf, t.max(), np.inf, np.inf], n_pulses)
    x_scale = np.tile([y_max, X_SCALE_TS_TAU1, X_SCALE_TS_TAU1, X_SCALE_TAU2], n_pulses)

    def residuals(params):
        return total_model(t, params, n_pulses) - y_raw

    res = least_squares(residuals, p0_raw, bounds=(lb, ub), x_scale=x_scale, max_nfev=max_nfev)
    return res  # raw-scale params in res.x; res.success/res.status carry convergence info


def run_burst(grb_dir_name: str, cfg: dict):
    paper_name = cfg["paper_name"]
    grb_dir = LC_DIR / grb_dir_name
    t, y_raw, dat_nai = summed_nai_curve(grb_dir, cfg["window"])
    y_max = float(np.max(y_raw))
    n_pulses = len(cfg["p0"])

    rows = []
    convergence = {}

    # --- Method 1: normalized ---
    try:
        params_norm = fit_normalized(t, y_raw, y_max, cfg["p0"], cfg["max_nfev"])
        model_norm_raw = total_model(t, params_norm, n_pulses) * y_max
        sse_norm = float(np.sum((model_norm_raw - y_raw) ** 2))
        convergence["normalized"] = "converged"
    except RuntimeError as exc:
        params_norm, model_norm_raw, sse_norm = None, None, None
        convergence["normalized"] = f"FAILED: {exc}"

    # --- Method 2: unnormalized, explicit x_scale ---
    res_raw = fit_unnormalized_xscale(t, y_raw, y_max, cfg["p0"], n_pulses, cfg["max_nfev"])
    if res_raw.success:
        params_raw = res_raw.x
        model_raw = total_model(t, params_raw, n_pulses)
        sse_raw = float(np.sum((model_raw - y_raw) ** 2))
        convergence["unnormalized_xscale"] = "converged"
    else:
        params_raw, model_raw, sse_raw = None, None, None
        convergence["unnormalized_xscale"] = f"FAILED: {res_raw.message}"

    for i in range(n_pulses):
        row = {"grb_dir": grb_dir_name, "grb_name": paper_name, "pulse_index": i + 1,
               "n_pulses": n_pulses, "y_max_cts_per_s": y_max, "detectors": "+".join(dat_nai)}
        if params_norm is not None:
            A, ts, tau1, tau2 = params_norm[i * 4:(i + 1) * 4]
            row.update({
                # A_norm is the raw fitted amplitude in peak-normalized units (y/y_max scale, in [0,1] by
                # construction of the fit bounds) -- the value NorrisFitter actually optimizes over, before
                # rescaling by y_max back to physical counts/s.
                "normalized_A_norm": A, "normalized_A_cts_per_s": A * y_max, "normalized_t_s": ts,
                "normalized_tau1": tau1, "normalized_tau2": tau2,
                "normalized_t_peak_s": t_peak(ts, tau1, tau2), "normalized_t_v_s": tv_value(tau1, tau2),
            })
        if params_raw is not None:
            A, ts, tau1, tau2 = params_raw[i * 4:(i + 1) * 4]
            row.update({
                "unnorm_xscale_A_cts_per_s": A, "unnorm_xscale_t_s": ts,
                "unnorm_xscale_tau1": tau1, "unnorm_xscale_tau2": tau2,
                "unnorm_xscale_t_peak_s": t_peak(ts, tau1, tau2), "unnorm_xscale_t_v_s": tv_value(tau1, tau2),
            })
        row["sse_normalized_rescaled"] = sse_norm
        row["sse_unnorm_xscale"] = sse_raw
        row["sse_ratio_unnorm_over_norm"] = (sse_raw / sse_norm) if (sse_norm and sse_raw) else None
        rows.append(row)

    # --- Plot: data + both fits overlaid, whichever converged ---
    update_style()
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.plot(t, y_raw, color="0.6", lw=LINE_WIDTH * 0.6,
            label=f"10-400 keV NaI ({'+'.join(dat_nai)}, summed)\nBackground subtracted")
    if model_norm_raw is not None:
        ax.plot(t, model_norm_raw, color="tab:blue", lw=LINE_WIDTH, label="Normalized-fit total (rescaled)")
    if model_raw is not None:
        ax.plot(t, model_raw, color="tab:red", ls="--", lw=LINE_WIDTH, label="Unnormalized fit (explicit x_scale)")

    # --- Individual pulses, both methods -- same show_individuals=True convention every production
    # fitter_*.py gets from NorrisFitter.plot_fit() (pymultifit's _plot_individual_fitter: dashed line,
    # one shade per pulse, parameters in the label), reproduced by hand here since this script fits one
    # method (fit_normalized) through NorrisFitter and the other (fit_unnormalized_xscale) through a bare
    # scipy.optimize.least_squares call, neither of which hands back a fitter object with its own
    # plot_fit(). Normalized-method pulses shaded blue, unnormalized-method pulses shaded red/orange --
    # same two-color convention as the two total-fit lines above -- so a method is identifiable at a
    # glance even with both sets of individuals overlaid.
    if params_norm is not None:
        norm_shades = plt.cm.Purples(np.linspace(0.4, 0.85, n_pulses))
        for i in range(n_pulses):
            par = params_norm[i * 4:(i + 1) * 4]
            A, ts, tau1, tau2 = par
            y_i = norris_pulse(t, par) * y_max
            ax.plot(t, y_i, ls=":", lw=LINE_WIDTH * 0.7, color=norm_shades[i],
                    label=f"N pulse {i + 1} (A={A:.3f}, t_s={ts:.2f}, tau1={tau1:.2f}, tau2={tau2:.2f})")
    if params_raw is not None:
        raw_shades = plt.cm.Oranges(np.linspace(0.4, 0.85, n_pulses))
        for i in range(n_pulses):
            par = params_raw[i * 4:(i + 1) * 4]
            A, ts, tau1, tau2 = par
            y_i = norris_pulse(t, par)
            ax.plot(t, y_i, ls="-.", lw=LINE_WIDTH * 0.7, color=raw_shades[i],
                    label=f"U pulse {i + 1} (A={A:.1f}, t_s={ts:.2f}, tau1={tau1:.2f}, tau2={tau2:.2f})")

    ax.set_xlabel("Time since trigger [s]")
    ax.set_ylabel("Count rate [counts/s]")
    ax.set_title(f"{paper_name}: normalized vs. unnormalized (x_scale) Norris fit")
    legend = ax.legend(fontsize="x-small", loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig_path = HERE / f"fit_comparison_{grb_dir_name}"
    # plt.show()
    save_fig(fig, fig_path, bbox_extra_artists=(legend,))

    return pd.DataFrame(rows), convergence


def main():
    all_rows = []
    all_convergence = {}
    for grb_dir_name, cfg in BURSTS.items():
        print(f"=== {cfg['paper_name']} ({grb_dir_name}) ===")
        df, convergence = run_burst(grb_dir_name, cfg)
        all_rows.append(df)
        all_convergence[cfg["paper_name"]] = convergence
        for method, status in convergence.items():
            print(f"  {method}: {status}")
        if not df.empty and df["sse_ratio_unnorm_over_norm"].notna().any():
            print(f"  SSE ratio (unnorm_xscale / normalized): {df['sse_ratio_unnorm_over_norm'].iloc[0]:.4f}")

    combined = pd.concat(all_rows, ignore_index=True)
    csv_path = HERE / "fit_comparison_results.csv"
    combined.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path} ({len(combined)} rows)")
    return combined, all_convergence


if __name__ == "__main__":
    main()
