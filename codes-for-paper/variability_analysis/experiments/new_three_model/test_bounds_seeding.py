"""Validation for bounds_seeding.py against norris_3param_spec.md section 3, and an end-to-end
smoke test: bounds + seed -> scipy.optimize.curve_fit on a synthetic pulse actually converges.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[3]  # new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from bounds_seeding import median_dt, pulse3_bounds, pulse3_bounds_and_seed, pulse3_seed
from pulse3 import make_pulse3

SEED = seed_from_name(__file__)  # project RNG convention -- src/grb_research/SEEDING.md


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def bounds_basic():
    print("\n[bounds] window=[0, 50], dt=0.5")
    t_window = np.array([0.0, 50.0])
    lb, ub = pulse3_bounds(t_window, dt=0.5)
    print(f"     lb={lb} ub={ub}")
    ok = lb == (0.0, 0.0, 1.0) and ub == (np.inf, 50.0, 50.0)
    return check("bounds match section-3 formula (A>=0, t_peak in window, t_v in [2dt, window])", ok)


def bounds_narrow_window_raises():
    print("\n[bounds] window narrower than 2*dt must raise")
    t_window = np.array([0.0, 0.5])  # window=0.5, dt inferred small but forced dt=1 below
    try:
        pulse3_bounds(t_window, dt=1.0)
        ok = False
    except ValueError as e:
        print(f"     raised ValueError: {e}")
        ok = True
    return check("ValueError raised when window <= 2*dt", ok)


def median_dt_check():
    print("\n[median_dt] evenly spaced grid")
    t = np.arange(0, 10, 0.25)
    dt = median_dt(t)
    ok = np.isclose(dt, 0.25)
    print(f"     median_dt = {dt!r}")
    return check("median_dt recovers the true bin width", ok)


def seed_recovers_synthetic_peak():
    print("\n[seed] neutral seed on a synthetic pulse (no noise)")
    r0_true = 190.0
    a_true, t_peak_true, t_v_true = 12.0, 20.0, 1.5
    model = make_pulse3(r0_true)
    dt = 0.05
    t_window = np.arange(0.0, 40.0, dt)
    y_window = model(t_window, a_true, t_peak_true, t_v_true)

    lb, ub, p0 = pulse3_bounds_and_seed(t_window, y_window, dt=dt)
    a0, t_peak0, t_v0 = p0
    print(f"     true:  A={a_true} t_peak={t_peak_true} t_v={t_v_true}")
    print(f"     seed:  A0={a0!r} t_peak0={t_peak0!r} t_v0={t_v0!r}")
    print(f"     bounds: lb={lb} ub={ub}")

    ok1 = check("t_peak0 within dt of true t_peak (argmax on a fine grid)", abs(t_peak0 - t_peak_true) <= dt)
    ok2 = check("A0 within 1% of true peak amplitude", abs(a0 - a_true) / a_true < 0.01)
    ok3 = check("t_v0 respects section-3 bounds [2dt, window]", lb[2] <= t_v0 <= ub[2])
    return ok1 and ok2 and ok3


def seed_clamped_in_narrow_window():
    print("\n[seed] t_v0 clamped when window is too narrow for the default 2.5*dt")
    dt = 1.0
    t_window = np.arange(0.0, 4.0, dt)  # window=3, 2.5*dt=2.5 (fine), but 3*dt=3 would need clamp
    y_window = np.array([0.0, 1.0, 3.0, 2.0])
    lb, ub, p0 = pulse3_bounds_and_seed(t_window, y_window, dt=dt)
    a0, t_peak0, t_v0 = p0
    print(f"     bounds t_v in [{lb[2]}, {ub[2]}], seed t_v0={t_v0!r}")
    ok = lb[2] <= t_v0 <= ub[2]
    return check("t_v0 stays within bounds even in a narrow window", ok)


def end_to_end_curve_fit():
    print("\n[end-to-end] bounds+seed -> curve_fit converges on a noisy synthetic pulse")
    rng = get_rng(seed=SEED)
    r0_true = 190.0
    a_true, t_peak_true, t_v_true = 12.0, 20.0, 1.5
    model = make_pulse3(r0_true)
    dt = 0.05
    t_window = np.arange(0.0, 40.0, dt)
    y_clean = model(t_window, a_true, t_peak_true, t_v_true)
    y_window = y_clean + rng.normal(scale=0.05 * a_true, size=y_clean.shape)

    lb, ub, p0 = pulse3_bounds_and_seed(t_window, y_window, dt=dt)
    popt, pcov = curve_fit(model, t_window, y_window, p0=p0, bounds=(lb, ub), maxfev=20000)
    a_fit, t_peak_fit, t_v_fit = popt
    print(f"     true: A={a_true} t_peak={t_peak_true} t_v={t_v_true}")
    print(f"     fit:  A={a_fit!r} t_peak={t_peak_fit!r} t_v={t_v_fit!r}")

    ok1 = check("fitted A within 5% of truth", abs(a_fit - a_true) / a_true < 0.05)
    ok2 = check("fitted t_peak within 0.1 of truth", abs(t_peak_fit - t_peak_true) < 0.1)
    ok3 = check("fitted t_v within 10% of truth", abs(t_v_fit - t_v_true) / t_v_true < 0.10)
    return ok1 and ok2 and ok3


if __name__ == "__main__":
    results = [
        bounds_basic(),
        bounds_narrow_window_raises(),
        median_dt_check(),
        seed_recovers_synthetic_peak(),
        seed_clamped_in_narrow_window(),
        end_to_end_curve_fit(),
    ]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
