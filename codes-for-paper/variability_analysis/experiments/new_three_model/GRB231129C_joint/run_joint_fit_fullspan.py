"""GRB231129C, joint 6-pulse fit over the FULL light-curve x-span -- the direct empirical test of
whether this 3-param joint fitter is affected by the fit window's x.min/x.max the way NorrisFitter
is (fitter_GRB231129779.py: START1=-np.inf, END1=np.inf -- the whole light curve, no window choice
at all). Uses the data's own t.min()/t.max() explicitly (not literal -inf/inf passed into a mask --
same effect once compared against real finite timestamps, but explicit about what it resolves to,
per user request) rather than run_joint_fit.py's WINDOW=(-2, 20).

Why this matters (see PROGRESS.md's "is it affected by x.min/x.max" discussion): bounds_seeding.py
ties t_peak's box to [x_min, x_max] and t_v's to [2*dt, window_width] -- structurally the same
pattern as NorrisFitter.fit_boundaries() tying t_s's lower bound to x_values.min(). NorrisFitter's
window-sensitivity comes specifically from the t_s<->tau1 degenerate direction that widening gives
more room to drift along; this model has no such direction (r0 is fixed), so the claim is that the
box-tied-to-window-edge pattern doesn't cause the same drift here. Section 8's window-widening
check already supports that claim indirectly (small widenings don't move the answer); fitting over
the ENTIRE light curve (the most extreme possible "widening") is the direct version of the same
test, not an extrapolation from it.

Plot is boxed to [T05-PLOT_PAD_S, T95+PLOT_PAD_S] (the production convention -- see
../GRB231129C/load_data.py's T05/T95/PLOT_PAD_S, copied from ../../../GRB231129C/_common.py) --
plot-axis bounding only, matching that file's own documented meaning: it changes what's shown, not
what's fitted, since what's fitted here is already the entire light curve.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_joint_fit import (  # noqa: E402
    check_edge_pinning,
    check_multistart_reproducibility,
    fit_joint,
    load_prior_single_pulse_results,
)
from bounds_seeding import median_dt  # noqa: E402
from joint_pulse3 import flatten_seed, joint_bounds, unflatten_params  # noqa: E402
from load_data import PLOT_PAD_S, T05, T95, load_light_curve  # noqa: E402
from plot_joint_fit import plot_joint_overlay  # noqa: E402
from grb_research import get_rng, seed_from_name  # noqa: E402

SEED = seed_from_name(__file__)


def main():
    print("=== GRB231129C: JOINT 6-pulse fit over the FULL light-curve x-span ===")
    t_full, y_full, sigma_full = load_light_curve()
    prior = load_prior_single_pulse_results()
    n_pulses = len(prior)

    # The window IS the data's own extent -- t_full.min()/t_full.max() explicitly, not a literal
    # -inf/inf mask (same effect against real finite timestamps, but explicit about what it
    # resolves to).
    t_window, y_window, sigma_window = t_full, y_full, sigma_full
    dt = median_dt(t_window)
    print(f"fit window: [{t_window.min():.3f}, {t_window.max():.3f}] s ({t_window.size} bins) -- "
          f"the data's own full extent, matching fitter_GRB231129779.py's START1=-inf, END1=inf")

    r0_list = prior["r0_chosen"].tolist()
    p0_flat = flatten_seed(list(zip(prior["amplitude"], prior["t_peak_s"], prior["t_v_s"])))
    lb, ub = joint_bounds(t_window, dt, n_pulses)
    print(f"t_peak bounds now: [{lb[1]:.3f}, {ub[1]:.3f}]  (was [-2.0, 20.0] in run_joint_fit.py)")
    print(f"t_v upper bound now: {ub[2]:.3f}  (was 22.0 in run_joint_fit.py)")

    print("\n--- Joint fit over the full span ---")
    model, popt, pcov, chi2_full = fit_joint(t_window, y_window, sigma_window, r0_list, p0_flat, dt, maxfev=60000)
    params_full = unflatten_params(popt, n_pulses)
    print(f"chi2 (full span, {t_window.size} bins) = {chi2_full:.4f}")

    print("\n--- Compare against the windowed (-2, 20) joint fit ---")
    windowed = pd.read_csv(HERE / "GRB231129C_joint_results.csv")
    print(f"{'pulse':>6} {'t_peak(windowed)':>18} {'t_peak(full-span)':>18} {'delta':>10}   "
          f"{'t_v(windowed)':>15} {'t_v(full-span)':>15} {'delta':>10}")
    max_dt_peak_rel = 0.0
    max_dt_v_rel = 0.0
    for i, pulse_index in enumerate(prior["pulse_index"].tolist()):
        row = windowed[windowed["pulse_index"] == pulse_index].iloc[0]
        a_f, tp_f, tv_f = params_full[i]
        d_tp = tp_f - row["t_peak_s"]
        d_tv = tv_f - row["t_v_s"]
        max_dt_peak_rel = max(max_dt_peak_rel, abs(d_tp) / max(abs(row["t_peak_s"]), 1.0))
        max_dt_v_rel = max(max_dt_v_rel, abs(d_tv) / row["t_v_s"])
        print(f"{pulse_index:>6} {row['t_peak_s']:>18.5f} {tp_f:>18.5f} {d_tp:>+10.5f}   "
              f"{row['t_v_s']:>15.5f} {tv_f:>15.5f} {d_tv:>+10.5f}")
    print(f"\nmax relative t_peak difference: {100*max_dt_peak_rel:.4f}%")
    print(f"max relative t_v difference: {100*max_dt_v_rel:.4f}%")

    print("\n--- Audit ---")
    pinned = check_edge_pinning(popt, lb, ub, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")

    rng = get_rng(seed=SEED)
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(
        t_window, y_window, sigma_window, r0_list, p0_flat, dt, rng
    )
    print(f"multi-start reproducibility: {'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")
    print("window-widening: not applicable -- this IS the full data span, nothing left to widen into "
          "(that's the whole point of this test).")

    print("\n--- Plot (fit: full span; display: boxed to [T05-pad, T95+pad], production convention) ---")
    plot_lo, plot_hi = T05 - PLOT_PAD_S, T95 + PLOT_PAD_S
    box_mask = (t_full >= plot_lo) & (t_full <= plot_hi)
    t_box, y_box, sigma_box = t_full[box_mask], y_full[box_mask], sigma_full[box_mask]
    paths = plot_joint_overlay(
        t_box, y_box, sigma_box, model, popt, out_dir=str(HERE), label="GRB231129C_joint_fullspan",
        pulse_indices=prior["pulse_index"].tolist(),
    )
    print(f"plot boxed to [{plot_lo:.3f}, {plot_hi:.3f}] s (T05={T05}, T95={T95}, pad={PLOT_PAD_S})")
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print("\n=== Verdict ===")
    if max_dt_peak_rel < 0.01 and max_dt_v_rel < 0.01 and not pinned and repro_ok:
        print("Fitting over the ENTIRE light curve (vs. the chosen (-2,20) window) changes every pulse's")
        print(f"t_peak/t_v by well under 1% ({100*max_dt_peak_rel:.4f}% / {100*max_dt_v_rel:.4f}% max), with no edge")
        print("pinning and exact multi-start reproducibility. Confirmed directly, not by extrapolation:")
        print("this fitter is NOT sensitive to the window's x.min/x.max the way NorrisFitter is.")
    else:
        print("Meaningful movement or audit failure found -- window choice DOES matter here after all;")
        print("see the per-pulse table above for which pulse(s) and by how much.")


if __name__ == "__main__":
    main()
