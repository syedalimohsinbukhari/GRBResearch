"""Validation for audit_checklist.py against norris_3param_spec.md section 8, on a clean
well-resolved synthetic pulse (the sharp scenario used throughout this folder) -- every checklist
item is expected to PASS here, since this is the "good fit" reference case; a real audit that
fails one of these on real data is the signal something upstream needs attention, not this test.

Dataset from synthetic_datasets.audit_source, seeded via this project's seed_from_name
convention (not a literal int).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[3]  # new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from audit_checklist import run_audit
from plot_diagnostics import plot_window_widening
from profile_scan import select_r0
from report_pulse import finalize_pulse
from synthetic_datasets import audit_source


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def main():
    print("\n[section 8] audit checklist on a clean, well-resolved synthetic pulse")
    dt = 0.03

    make_window_data = audit_source(t_full_min=-20.0, t_full_max=60.0, noise_frac=0.03)
    t_window, y_window, sigma = make_window_data(0.0, 30.0)  # base window: full pulse, some baseline padding

    selection = select_r0(t_window, y_window, sigma=sigma, dt=dt)
    reported = finalize_pulse(t_window, y_window, selection, sigma=sigma, dt=dt)
    print(f"     r0_chosen={reported.r0:.4g}  A={reported.amplitude:.4f}  t_peak={reported.t_peak:.4f}  t_v={reported.t_v:.4f}")

    report = run_audit(
        t_window,
        y_window,
        selection,
        reported,
        sigma=sigma,
        dt=dt,
        make_window_data=make_window_data,
        widen_factors=(1.5, 2.0),
        rng=get_rng(seed=seed_from_name(__file__)),
    )
    print("\n" + str(report))

    window_item = next(item for item in report.items if item.name == "window-widening invariance")
    if window_item.data is not None:
        paths = plot_window_widening(window_item.data, out_dir=".", label="section8_window_widening")
        print(f"\n     wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    results = [check(item.name, item.passed) for item in report.items]
    return all(results)


if __name__ == "__main__":
    ok = main()
    print()
    if ok:
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
