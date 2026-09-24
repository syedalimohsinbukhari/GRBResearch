"""GRB231129C end-to-end run: the first real multi-pulse light curve through the full section
1-8 pipeline (profile scan -> finalize -> audit), per norris_3param_reduction.md's "Known
limitations" and the external review's sequencing note (review-on-progress.md, ../PROGRESS.md).

Method (review-on-progress.md: "run the scan on the current-best joint model's residuals for
that pulse (subtract other pulses, scan the target on what remains), and iterate once if the
joint fit moves things"):

  Pass 1: for each of the archived 6 pulses, build a residual light curve = full data minus the
          OTHER 5 pulses' archived 4-param model, then run select_r0 -> finalize_pulse ->
          run_audit on that residual within a shared window.
  Pass 2: rebuild each pulse's residual using pass 1's NEW 3-param results (converted to raw
          Norris params) for the other 5 pulses instead of the archived ones, and refit -- checks
          how much iterating once actually moves things.

Comparison baseline (review point 5): new point-estimate t_v is compared against the archived
POINT estimate tv_value(tau1, tau2) at the archived best-fit params, NOT the archived t_v_s
column (which section 6 already established is an MC-propagated median, not a point estimate --
comparing against it would manufacture disagreements on exactly the degenerate pulses this run
cares about).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
PROJECT_ROOT = NEW_THREE_MODEL_DIR.parents[3]
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from audit_checklist import run_audit  # noqa: E402
from bounds_seeding import median_dt  # noqa: E402
from load_data import EPISODE_BOUNDS, load_archived_pulses, load_light_curve, load_project_r_population  # noqa: E402
from plot_diagnostics import plot_fit_overlay  # noqa: E402
from profile_scan import DEFAULT_R0_GRID, select_r0  # noqa: E402
from pulse3 import make_pulse3, norris_raw, width_function  # noqa: E402
from report_pulse import finalize_pulse  # noqa: E402

WINDOW = (-2.0, 10.0)  # shared fit window: covers all 6 archived t_peak (0.65-5.5s) with margin
WIDEN_FACTORS = (1.5, 2.0)
SEED = seed_from_name(__file__)


def build_residual(t_full, y_full, other_pulses_raw):
    """y_full minus the sum of norris_raw(t; *params) for every (A, t_s, tau1, tau2) in
    other_pulses_raw -- the "residual after subtracting already-fitted pulses" section 3 already
    names for multi-pulse seeding, made explicit here for the whole per-pulse fit, not just the seed."""
    model_sum = np.zeros_like(y_full)
    for amplitude, t_s, tau1, tau2 in other_pulses_raw:
        model_sum += norris_raw(t_full, amplitude, t_s, tau1, tau2)
    return y_full - model_sum


def make_window_data_factory(t_full, residual_full, sigma_full):
    def make_window_data(t_min, t_max):
        mask = (t_full >= t_min) & (t_full <= t_max)
        return t_full[mask], residual_full[mask], sigma_full[mask]

    return make_window_data


def select_r0_with_edge_check(t_window, y_window, sigma, dt, archived_r_median, pulse_index):
    """select_r0, but if the flat zone touches the DEFAULT_R0_GRID edge, the spec's own
    instruction ("extend if the profile is still falling at the edges") applies -- retry once
    with a grid extended by 2 more decades in the direction that was touching."""
    selection = select_r0(t_window, y_window, sigma=sigma, dt=dt, archived_r_median=archived_r_median)
    zone_lo, zone_hi = selection.classification.flat_zone
    grid_lo, grid_hi = DEFAULT_R0_GRID[0], DEFAULT_R0_GRID[-1]
    touches_lo = np.isclose(zone_lo, grid_lo)
    touches_hi = np.isclose(zone_hi, grid_hi)
    if not (touches_lo or touches_hi):
        return selection, False

    new_lo_exp = np.log10(grid_lo) - (2 if touches_lo else 0)
    new_hi_exp = np.log10(grid_hi) + (2 if touches_hi else 0)
    n_points = int(round((new_hi_exp - new_lo_exp) / (np.log10(grid_hi) - np.log10(grid_lo)) * len(DEFAULT_R0_GRID)))
    widened_grid = np.logspace(new_lo_exp, new_hi_exp, n_points)
    print(
        f"     pulse {pulse_index}: flat zone touched the default r0 grid edge "
        f"({selection.classification.flat_zone}) -- retrying with grid logspace({new_lo_exp:.0f}, {new_hi_exp:.0f}, {n_points})"
    )
    selection2 = select_r0(t_window, y_window, r0_grid=widened_grid, sigma=sigma, dt=dt, archived_r_median=archived_r_median)
    return selection2, True


def fit_one_pass(t_full, y_full, sigma_full, pulses_raw: dict, archived_r_median: float, rng, pass_label: str):
    """One full pass over all pulses: build each pulse's residual from the CURRENT estimate of
    the other 5 pulses (pulses_raw), then select_r0 -> finalize_pulse -> run_audit on it."""
    results = {}
    for i, params_i in pulses_raw.items():
        others = [p for j, p in pulses_raw.items() if j != i]
        residual_full = build_residual(t_full, y_full, others)
        make_window_data = make_window_data_factory(t_full, residual_full, sigma_full)

        t_window, y_window, sigma_window = make_window_data(*WINDOW)
        dt = median_dt(t_window)

        selection, grid_widened = select_r0_with_edge_check(
            t_window, y_window, sigma_window, dt, archived_r_median, pulse_index=i
        )
        reported = finalize_pulse(t_window, y_window, selection, sigma=sigma_window, dt=dt)
        audit = run_audit(
            t_window, y_window, selection, reported, rng,
            sigma=sigma_window, dt=dt, make_window_data=make_window_data, widen_factors=WIDEN_FACTORS,
        )

        print(
            f"  [{pass_label}] pulse {i}: r0={reported.r0:.4g} ({selection.classification.shape}, "
            f"grid_widened={grid_widened})  A={reported.amplitude:.2f}+/-{reported.amplitude_err:.2f}  "
            f"t_peak={reported.t_peak:.4f}+/-{reported.t_peak_err:.4f}  t_v={reported.t_v:.4f}+/-{reported.t_v_err:.4f}  "
            f"audit_pass={audit.all_passed}"
        )

        results[i] = {
            "selection": selection,
            "reported": reported,
            "audit": audit,
            "t_window": t_window,
            "y_window": y_window,
            "sigma_window": sigma_window,
            "grid_widened": grid_widened,
        }
    return results


def raw_params_from_reported(reported) -> tuple:
    return reported.amplitude, reported.t_s, reported.tau1, reported.tau2


def results_dataframe(archived: pd.DataFrame, pass_results: dict, pass_label: str) -> pd.DataFrame:
    rows = []
    for i, r in pass_results.items():
        rep = r["reported"]
        sel = r["selection"]
        arch = archived[archived["pulse_index"] == i].iloc[0]
        t_v_archived_point = arch["tau2"] * width_function(arch["r"])
        matched_episodes = [name for name, (lo, hi) in EPISODE_BOUNDS.items() if lo <= rep.t_peak <= hi]
        rows.append(
            {
                "pass": pass_label,
                "pulse_index": i,
                "episode": "+".join(matched_episodes) or None,
                "r0_chosen": rep.r0,
                "profile_shape": sel.classification.shape,
                "grid_widened": r["grid_widened"],
                "amplitude": rep.amplitude,
                "amplitude_err": rep.amplitude_err,
                "t_peak_s": rep.t_peak,
                "t_peak_err_s": rep.t_peak_err,
                "t_v_s": rep.t_v,
                "t_v_err_s": rep.t_v_err,
                "t_s": rep.t_s,
                "tau1": rep.tau1,
                "tau1_err": rep.tau1_err,
                "tau2": rep.tau2,
                "tau2_err": rep.tau2_err,
                "delta_chi2": sel.delta_chi2["delta_chi2"],
                "delta_chi2_verdict": sel.delta_chi2["verdict"],
                "r_pinned_4param": sel.delta_chi2["r_pinned"],
                "audit_all_passed": r["audit"].all_passed,
                "archived_t_peak_s": arch["t_peak_s"],
                "archived_t_v_point_s": t_v_archived_point,
                "archived_r": arch["r"],
                "delta_t_peak_vs_archived_s": rep.t_peak - arch["t_peak_s"],
                "delta_t_v_vs_archived_s": rep.t_v - t_v_archived_point,
            }
        )
    return pd.DataFrame(rows)


def main():
    print("=== GRB231129C: 3-param pipeline, real multi-pulse light curve ===")
    t_full, y_full, sigma_full = load_light_curve()
    archived = load_archived_pulses()
    r_population = load_project_r_population()
    archived_r_median = float(np.median(r_population))
    print(f"loaded light curve: {t_full.size} bins, t in [{t_full.min():.3f}, {t_full.max():.3f}]")
    print(f"archived pulses: {len(archived)}")
    print(f"project-wide archived r population: n={len(r_population)}, median={archived_r_median:.4g}, "
          f"p16/p84=[{np.percentile(r_population, 16):.4g}, {np.percentile(r_population, 84):.4g}]")

    pulses_raw_archived = {
        int(row.pulse_index): (row.A_cts_per_s, row.t_s, row.tau1, row.tau2) for row in archived.itertuples()
    }

    rng = get_rng(seed=SEED)

    print("\n--- Pass 1: residuals built from ARCHIVED other-pulse params ---")
    pass1 = fit_one_pass(t_full, y_full, sigma_full, pulses_raw_archived, archived_r_median, rng, "pass1")

    pulses_raw_pass1 = {i: raw_params_from_reported(r["reported"]) for i, r in pass1.items()}

    print("\n--- Pass 2: residuals rebuilt from PASS-1 3-param other-pulse params ---")
    pass2 = fit_one_pass(t_full, y_full, sigma_full, pulses_raw_pass1, archived_r_median, rng, "pass2")

    df1 = results_dataframe(archived, pass1, "pass1")
    df2 = results_dataframe(archived, pass2, "pass2")
    df = pd.concat([df1, df2], ignore_index=True)
    csv_path = HERE / "GRB231129C_3param_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path.name} ({len(df)} rows)")

    print("\n--- Pass1 -> Pass2 movement (how much iterating once changed things) ---")
    for i in pulses_raw_archived:
        t_peak_move = pass2[i]["reported"].t_peak - pass1[i]["reported"].t_peak
        t_v_move = pass2[i]["reported"].t_v - pass1[i]["reported"].t_v
        print(f"  pulse {i}: delta_t_peak={t_peak_move:+.5f} s, delta_t_v={t_v_move:+.5f} s")

    print("\n--- Fit overlay plots (pass 2, final) ---")
    for i, r in pass2.items():
        model = make_pulse3(r["reported"].r0)
        y_fit = model(r["t_window"], r["reported"].amplitude, r["reported"].t_peak, r["reported"].t_v)
        paths = plot_fit_overlay(
            r["t_window"], r["y_window"], y_fit, out_dir=str(HERE), label=f"pulse{i}_pass2",
            sigma=r["sigma_window"],
            extra_title=f"r0={r['reported'].r0:.3g}, t_peak={r['reported'].t_peak:.3f}, t_v={r['reported'].t_v:.3f}",
        )
        print(f"  wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    n_audit_pass = sum(1 for r in pass2.values() if r["audit"].all_passed)
    print(f"\naudit: {n_audit_pass}/{len(pass2)} pulses fully passed all section-8 checks in pass 2")
    for i, r in pass2.items():
        if not r["audit"].all_passed:
            print(f"  pulse {i} audit detail:\n" + "\n".join(f"    {item.name}: {'PASS' if item.passed else 'FAIL'} -- {item.detail}" for item in r["audit"].items))


if __name__ == "__main__":
    main()
