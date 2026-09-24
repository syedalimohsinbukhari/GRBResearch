"""GRB231129C, joint N-pulse 3-param fit: pulses' (A, t_peak, t_v) fit SIMULTANEOUSLY against the
raw light curve in one curve_fit call, each pulse keeping the r0 already chosen by
../GRB231129C/run_pipeline.py's per-pulse residual method (pass 2) as a FIXED input.

Why this exists as a separate run from ../GRB231129C/ (kept in its own folder per user request):
that folder's method fits one pulse at a time against a residual built by subtracting an
*approximation* of the other 5 pulses (archived, then pass-1's own results) -- pulse 6 there
showed a visible neighbor-contamination artifact traceable to that approximation. This script
removes the approximation entirely: every pulse sees every other pulse's actual current shape
during optimization, not a frozen estimate of it, because they are all being fit together.

Runs TWO variants, both kept (user request, 2026-09-25): the full 6-pulse decomposition, and a
5-pulse variant dropping pulse 4. Pulse 4 comes from a LATER, full-range 6-pulse fit
(GRB231129C/_common.py) that added a pulse not present in the ORIGINAL production 5-pulse fit
(fitter_GRB231129779.py, window (-1,10)) -- its own P0 entry is marked "# replaceable" in that
source, and the user independently confirmed it's negligible. Our own numbers agreed before being
told that: pulse 4 showed the largest joint-vs-single-pulse movement of any pulse in the 6-pulse
run, and the largest window-widening drift before the window was corrected. See PROGRESS.md for
the full writeup of this discovery.

Seeds for the joint fit come from ../GRB231129C/GRB231129C_3param_results.csv (pass 2), which is
read-only reuse of an already-computed, already-validated near-answer -- not a re-derivation, and
not a dependency this script writes to.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
GRB231129C_DIR = NEW_THREE_MODEL_DIR / "GRB231129C"
PROJECT_ROOT = NEW_THREE_MODEL_DIR.parents[3]
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(GRB231129C_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from bounds_seeding import median_dt  # noqa: E402
from joint_pulse3 import flatten_seed, joint_bounds, make_joint_pulse3, pulse_block_cov, residual_excluding, unflatten_params  # noqa: E402
from load_data import EPISODE_BOUNDS, load_archived_pulses, load_light_curve  # noqa: E402
from plot_joint_fit import plot_joint_overlay  # noqa: E402
from profile_scan import chi_square, refine_r0_zone  # noqa: E402
from pulse3 import norris_raw, width_function  # noqa: E402
from report_pulse import _zone_spread  # noqa: E402

WINDOW = (-2.0, 20.0)  # widened from ../GRB231129C/run_pipeline.py's (-2, 10): that window's own
# window-widening audit failed here (pulse 4/6 t_v drifted up to ~10%) because [10, 16]s has real
# excess flux (~1.9 sigma/bin, checked directly against the raw light curve) that (-2, 10) was
# cutting off -- almost certainly pulse 6's own decay tail (it's the broadest pulse, tau2~1.8-2.0s).
# Widened to comfortably include where the excess decays back to noise (~t=16-20s).
WIDEN_FACTORS = (1.5, 2.0)
N_JITTER_SEEDS = 3  # joint 18-dim refits are pricier than single-pulse ones; 3 is enough to catch drift
EDGE_REL_TOL = 1e-3
REPRO_REL_TOL = 1e-3
SEED = seed_from_name(__file__)


def load_prior_single_pulse_results() -> pd.DataFrame:
    """pass-2 rows from ../GRB231129C/'s already-computed run: r0_chosen (fixed here) and
    (amplitude, t_peak_s, t_v_s) as the joint fit's seed."""
    path = GRB231129C_DIR / "GRB231129C_3param_results.csv"
    df = pd.read_csv(path)
    df = df[df["pass"] == "pass2"].sort_values("pulse_index").reset_index(drop=True)
    return df


def fit_joint(t_window, y_window, sigma_window, r0_list, p0_flat, dt, maxfev=40000):
    n = len(r0_list)
    model = make_joint_pulse3(r0_list)
    lb, ub = joint_bounds(t_window, dt, n)
    popt, pcov = curve_fit(
        model, t_window, y_window, p0=p0_flat, bounds=(lb, ub),
        sigma=sigma_window, absolute_sigma=True, maxfev=maxfev,
    )
    y_fit = model(t_window, *popt)
    chi2 = chi_square(y_window, y_fit, sigma=sigma_window)
    return model, popt, pcov, chi2


def check_edge_pinning(popt, lb, ub, n_pulses):
    """Bug found and fixed (2026-09-25, during the full-x-span test): the original version used a
    tolerance relative to the bound SPAN (hi - lo), matching audit_checklist.py's single-pulse
    `_is_pinned` helper. That breaks down once the window (and so the span) becomes very large
    relative to the fitted value itself -- fitting over GRB231129C's full 614s light curve made
    t_v's upper bound 614.3, so 1e-3*span ~= 0.6s of absolute slop, enough to falsely flag
    ordinary ~0.5s t_v values as "pinned" when they were nowhere near either edge. Fixed to match
    profile_scan.py's already-correct `_is_pinned_r`: tolerance relative to the bound VALUE
    itself (np.isclose-style), not the span between bounds.
    """
    pinned = []
    for i in range(n_pulses):
        for j, name in enumerate(("A", "t_peak", "t_v")):
            idx = 3 * i + j
            value, lo, hi = popt[idx], lb[idx], ub[idx]
            near_lo = np.isclose(value, lo, rtol=EDGE_REL_TOL, atol=EDGE_REL_TOL * max(abs(lo), 1e-3))
            near_hi = np.isfinite(hi) and np.isclose(value, hi, rtol=EDGE_REL_TOL, atol=EDGE_REL_TOL * max(abs(hi), 1e-3))
            if near_lo or near_hi:
                pinned.append((i + 1, name, value, lo, hi))
    return pinned


def check_multistart_reproducibility(t_window, y_window, sigma_window, r0_list, p0_flat, dt, rng):
    """N_JITTER_SEEDS joint refits from perturbed seeds (+/-10% on A/t_v, +/-3 bins on t_peak,
    smaller than bounds_seeding.jittered_seeds' single-pulse defaults since 18 simultaneous
    parameters have more room to drift if perturbed too hard) -- checks the joint optimum is the
    same basin regardless of small seed perturbations, the joint analog of section 7 T4a.

    Jitter is clamped into bounds (bug found 2026-09-25, GRB131014215_joint's r0-iteration run):
    a relative jitter on a seed value already sitting near its own bound (e.g. t_v at the 2*dt
    resolution floor) can land just outside it, which curve_fit rejects outright
    ("Initial guess is outside of provided bounds") instead of a graceful clip.
    """
    n = len(r0_list)
    lb, ub = joint_bounds(t_window, dt, n)
    solutions = []
    for _ in range(N_JITTER_SEEDS):
        jittered = []
        for i in range(n):
            a0, tp0, tv0 = p0_flat[3 * i : 3 * i + 3]
            a = float(np.clip(a0 * (1 + rng.uniform(-0.1, 0.1)), lb[3 * i], ub[3 * i]))
            tp = float(np.clip(tp0 + rng.uniform(-3, 3) * dt, lb[3 * i + 1], ub[3 * i + 1]))
            tv = float(np.clip(tv0 * (1 + rng.uniform(-0.1, 0.1)), lb[3 * i + 2], ub[3 * i + 2]))
            jittered.extend((a, tp, tv))
        _, popt_j, _, chi2_j = fit_joint(t_window, y_window, sigma_window, r0_list, tuple(jittered), dt)
        solutions.append((popt_j, chi2_j))
    arr = np.array([s[0] for s in solutions])
    spread = arr.max(axis=0) - arr.min(axis=0)
    scale = np.abs(arr).mean(axis=0)
    ok = True
    for i in range(n):
        a_ok = spread[3 * i] < REPRO_REL_TOL * scale[3 * i]
        tp_ok = spread[3 * i + 1] < REPRO_REL_TOL
        tv_ok = spread[3 * i + 2] < REPRO_REL_TOL * scale[3 * i + 2]
        ok &= a_ok and tp_ok and tv_ok
    chi2_spread = max(s[1] for s in solutions) - min(s[1] for s in solutions)
    return ok, spread, chi2_spread


def check_window_widening(t_full, y_full, sigma_full, r0_list, popt_base, dt, widen_factors=WIDEN_FACTORS):
    n = len(r0_list)
    base_params = unflatten_params(popt_base, n)
    results = {factor: None for factor in widen_factors}
    t_min0, t_max0 = WINDOW
    ok = True
    for factor in widen_factors:
        half_extra = 0.5 * (factor - 1.0) * (t_max0 - t_min0)
        t_min, t_max = t_min0 - half_extra, t_max0 + half_extra
        mask = (t_full >= t_min) & (t_full <= t_max)
        t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
        _, popt_w, _, _ = fit_joint(t_w, y_w, sigma_w, r0_list, tuple(popt_base), dt)
        widened_params = unflatten_params(popt_w, n)
        results[factor] = widened_params
        for i in range(n):
            t_peak_ok = abs(widened_params[i][1] - base_params[i][1]) <= 0.05 * max(abs(base_params[i][1]), 1.0)
            t_v_ok = abs(widened_params[i][2] - base_params[i][2]) <= 0.05 * base_params[i][2]
            ok &= t_peak_ok and t_v_ok
    return ok, results


def run_variant(pulse_indices, label, title):
    """One full joint-fit run (fit -> per-pulse extraction -> audit -> plot -> CSV) restricted to
    `pulse_indices` (a list of 1-indexed pulse numbers, e.g. [1,2,3,4,5,6] or [1,2,3,5,6]).
    Outputs are named by `label`, so different variants never overwrite each other."""
    print(f"\n{'=' * 70}\n=== {title} ===\n{'=' * 70}")
    t_full, y_full, sigma_full = load_light_curve()
    archived_all = load_archived_pulses()
    prior_all = load_prior_single_pulse_results()

    archived = archived_all[archived_all["pulse_index"].isin(pulse_indices)].reset_index(drop=True)
    prior = prior_all[prior_all["pulse_index"].isin(pulse_indices)].sort_values("pulse_index").reset_index(drop=True)
    n_pulses = len(prior)
    print(f"pulses included: {pulse_indices} (n={n_pulses})")

    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_window, y_window, sigma_window = t_full[mask], y_full[mask], sigma_full[mask]
    dt = median_dt(t_window)

    r0_list = prior["r0_chosen"].tolist()
    p0_per_pulse = list(zip(prior["amplitude"], prior["t_peak_s"], prior["t_v_s"]))
    p0_flat = flatten_seed(p0_per_pulse)
    lb, ub = joint_bounds(t_window, dt, n_pulses)

    print("\n--- Joint fit ---")
    model, popt, pcov, chi2_joint = fit_joint(t_window, y_window, sigma_window, r0_list, p0_flat, dt)
    joint_params = unflatten_params(popt, n_pulses)
    print(f"chi2_joint ({3 * n_pulses} free params) = {chi2_joint:.4f}")

    y_archived = np.zeros_like(y_window)
    for row in archived.itertuples():
        y_archived += norris_raw(t_window, row.A_cts_per_s, row.t_s, row.tau1, row.tau2)
    chi2_archived_joint = chi_square(y_window, y_archived, sigma=sigma_window)
    print(f"chi2_archived ({4 * n_pulses} free params, scored not refit on this window) = {chi2_archived_joint:.4f}")
    print(f"delta_chi2 (archived - joint 3-param) = {chi2_archived_joint - chi2_joint:.4f}")

    print("\n--- Per-pulse extraction + r0-systematic (via refine_r0_zone on joint-consistent residuals) ---")
    rows = []
    pulse_index_list = prior["pulse_index"].tolist()
    for i, pulse_index in enumerate(pulse_index_list):
        amplitude, t_peak, t_v = joint_params[i]
        r0 = r0_list[i]
        block = pulse_block_cov(pcov, i)
        a_err, t_peak_err_stat, t_v_err_stat = np.sqrt(np.diag(block))

        residual_i = residual_excluding(model, t_window, y_window, popt, exclude_index=i)
        refined_scan, refined_classification = refine_r0_zone(t_window, residual_i, r0, sigma=sigma_window, dt=dt)
        t_peak_err_sys, t_v_err_sys = _zone_spread(refined_scan, *refined_classification.flat_zone)
        t_peak_err = float(np.hypot(t_peak_err_stat, t_peak_err_sys))
        t_v_err = float(np.hypot(t_v_err_stat, t_v_err_sys))

        mapping = model.mappings[i]
        _, t_s, tau1, tau2 = mapping.to_raw(amplitude, t_peak, t_v)

        arch = archived[archived["pulse_index"] == pulse_index].iloc[0]
        t_v_archived_point = arch["tau2"] * width_function(arch["r"])
        matched_episodes = [name for name, (lo, hi) in EPISODE_BOUNDS.items() if lo <= t_peak <= hi]

        prior_row = prior[prior["pulse_index"] == pulse_index].iloc[0]

        print(
            f"  pulse {pulse_index}: r0={r0:.4g}  A={amplitude:.2f}+/-{a_err:.2f}  "
            f"t_peak={t_peak:.4f}+/-{t_peak_err:.4f} (stat={t_peak_err_stat:.4f},sys={t_peak_err_sys:.4f})  "
            f"t_v={t_v:.4f}+/-{t_v_err:.4f} (stat={t_v_err_stat:.4f},sys={t_v_err_sys:.4f})"
        )

        rows.append(
            {
                "pulse_index": pulse_index,
                "episode": "+".join(matched_episodes) or None,
                "r0": r0,
                "amplitude": amplitude,
                "amplitude_err": a_err,
                "t_peak_s": t_peak,
                "t_peak_err_stat_s": t_peak_err_stat,
                "t_peak_err_sys_s": t_peak_err_sys,
                "t_peak_err_s": t_peak_err,
                "t_v_s": t_v,
                "t_v_err_stat_s": t_v_err_stat,
                "t_v_err_sys_s": t_v_err_sys,
                "t_v_err_s": t_v_err,
                "t_s": t_s,
                "tau1": tau1,
                "tau2": tau2,
                "archived_t_peak_s": arch["t_peak_s"],
                "archived_t_v_point_s": t_v_archived_point,
                "delta_t_peak_vs_archived_s": t_peak - arch["t_peak_s"],
                "delta_t_v_vs_archived_s": t_v - t_v_archived_point,
                "single_pulse_t_peak_s": prior_row["t_peak_s"],
                "single_pulse_t_v_s": prior_row["t_v_s"],
                "delta_t_peak_vs_single_pulse_s": t_peak - prior_row["t_peak_s"],
                "delta_t_v_vs_single_pulse_s": t_v - prior_row["t_v_s"],
            }
        )

    df = pd.DataFrame(rows)
    csv_path = HERE / f"{label}_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path.name}")

    print("\n--- Audit (joint analogs of section 8) ---")
    pinned = check_edge_pinning(popt, lb, ub, n_pulses)
    print(f"edge pinning: {'none' if not pinned else pinned}")

    rng = get_rng(seed=seed_from_name(f"{__file__}::{label}"))
    repro_ok, repro_spread, repro_chi2_spread = check_multistart_reproducibility(
        t_window, y_window, sigma_window, r0_list, p0_flat, dt, rng
    )
    print(f"multi-start reproducibility ({N_JITTER_SEEDS} jittered joint refits): "
          f"{'PASS' if repro_ok else 'FAIL'}, chi2 spread={repro_chi2_spread:.6f}")

    widen_ok, widen_results = check_window_widening(t_full, y_full, sigma_full, r0_list, popt, dt)
    print(f"window-widening invariance (x1.5, x2.0): {'PASS' if widen_ok else 'FAIL'}")

    audit_all_passed = (not pinned) and repro_ok and widen_ok
    print(f"\naudit overall: {'ALL PASS' if audit_all_passed else 'SOME CHECKS FAILED'}")

    print("\n--- Plot ---")
    paths = plot_joint_overlay(
        t_window, y_window, sigma_window, model, popt, out_dir=str(HERE), label=label, pulse_indices=pulse_index_list
    )
    print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print("\n--- Joint vs. single-pulse-residual method: how much did switching methods move things ---")
    for row in rows:
        print(
            f"  pulse {row['pulse_index']}: delta_t_peak={row['delta_t_peak_vs_single_pulse_s']:+.5f} s, "
            f"delta_t_v={row['delta_t_v_vs_single_pulse_s']:+.5f} s"
        )

    return {
        "df": df,
        "chi2_joint": chi2_joint,
        "chi2_archived": chi2_archived_joint,
        "n_pulses": n_pulses,
        "audit_all_passed": audit_all_passed,
    }


def main():
    result_6pulse = run_variant(
        pulse_indices=[1, 2, 3, 4, 5, 6],
        label="GRB231129C_joint",
        title="GRB231129C: JOINT fit, full 6-pulse decomposition",
    )
    result_5pulse = run_variant(
        pulse_indices=[1, 2, 3, 5, 6],
        label="GRB231129C_joint_5pulse",
        title="GRB231129C: JOINT fit, 5-pulse decomposition (pulse 4 dropped -- see module docstring)",
    )

    print(f"\n{'=' * 70}\n=== 6-pulse vs. 5-pulse summary ===\n{'=' * 70}")
    print(f"6-pulse: chi2_joint={result_6pulse['chi2_joint']:.4f} (18 params), audit={result_6pulse['audit_all_passed']}")
    print(f"5-pulse: chi2_joint={result_5pulse['chi2_joint']:.4f} (15 params), audit={result_5pulse['audit_all_passed']}")
    print("(these chi2 are over different data -- 5-pulse's residual is NOT pulse 4's flux -- so this is not a")
    print(" direct nested-model comparison; see PROGRESS.md for how the two variants are actually compared.)")


if __name__ == "__main__":
    main()
