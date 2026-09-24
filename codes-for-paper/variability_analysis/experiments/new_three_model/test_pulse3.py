"""Validation for pulse3.py against norris_3param_spec.md, before this wrapper touches real data.

Covers the section-2/section-1 pieces relevant to what pulse3.py implements: T1 (round-trip) and
T2 (f(r) properties) from section 7, plus the overflow-safety and API-consistency checks section 2
calls for. T3-T5 (FWHM/2 via brentq, synthetic degeneracy fit, chi^2 profile) need an optimizer and
belong with the section-3/5 fitting code, not this wrapper.
"""
import numpy as np

from pulse3 import (
    LN2,
    Pulse3Mapping,
    make_pulse3,
    norris_raw,
    pulse3,
    pulse4,
    reduced_from_raw,
    width_function,
    width_function_inverse,
)


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def t1_round_trip():
    print("\n[T1] round-trip: (tau1, tau2) = (85.6, 0.453)")
    tau1, tau2, t_s0 = 85.6, 0.453, 3.0
    t_peak, t_v, r = reduced_from_raw(t_s0, tau1, tau2)
    print(f"     t_peak={t_peak!r} t_v={t_v!r} r={r!r}")

    mapping = Pulse3Mapping(r0=r)
    _, t_s, tau1_rec, tau2_rec = mapping.to_raw(amplitude=1.0, t_peak=t_peak, t_v=t_v)

    ok = (
        np.isclose(tau1_rec, tau1, atol=0, rtol=1e-12)
        and np.isclose(tau2_rec, tau2, atol=0, rtol=1e-12)
        and np.isclose(t_s, t_s0, atol=1e-10, rtol=1e-12)
    )
    print(f"     recovered tau1={tau1_rec!r} tau2={tau2_rec!r} t_s={t_s!r}")
    return check("round-trip recovers (tau1, tau2, t_s) to machine precision", ok)


def t2_width_function():
    print("\n[T2] width_function(r) properties")
    results = []

    f0 = width_function(0.0)
    results.append(check(f"f(0) = {f0!r} == ln2/2 = {LN2/2!r}", np.isclose(f0, LN2 / 2)))

    r_grid = np.logspace(-6, 8, 2000)
    f_grid = width_function(r_grid)
    results.append(check("monotone increasing on logspace(-6, 8)", np.all(np.diff(f_grid) > 0)))

    f_big = width_function(1e8)
    results.append(check(f"f(1e8) = {f_big!r} ~ 83.26", np.isclose(f_big, 83.26, atol=0.01)))

    f_sample = width_function(r_grid)
    r_back = width_function_inverse(f_sample)
    results.append(check("r(f) inverse round-trips over the same grid", np.allclose(r_back, r_grid, rtol=1e-8)))

    return all(results)


def overflow_safety():
    print("\n[extra] overflow safety for mu > ~355 (factored form would overflow here)")
    # r0 = 2e6 -> mu = sqrt(r0) ~ 1414, exp(2*mu) alone would overflow float64 (max exp arg ~709).
    mapping = Pulse3Mapping(r0=2e6)
    t_v = 1.0
    _, t_s, tau1, tau2 = mapping.to_raw(amplitude=1.0, t_peak=10.0, t_v=t_v)
    mu = np.sqrt(tau1 / tau2)
    print(f"     r0=2e6 -> mu={mu!r} (2*mu={2*mu!r}, exp(2*mu) alone would overflow)")

    t = np.linspace(0, 20, 2000)
    y = norris_raw(t, amplitude=1.0, t_s=t_s, tau1=tau1, tau2=tau2)
    ok = np.all(np.isfinite(y)) and not np.any(y < 0)
    return check("norris_raw stays finite and non-negative across the whole grid", ok)


def api_consistency():
    print("\n[extra] make_pulse3(r0) matches pulse3(..., r0) and pulse4 at r=r0")
    t = np.linspace(-5, 30, 500)
    r0 = 190.0
    a, t_peak, t_v = 7.5, 12.0, 1.3

    y_a = pulse3(t, a, t_peak, t_v, r0)
    model = make_pulse3(r0)
    y_b = model(t, a, t_peak, t_v)
    y_c = pulse4(t, a, t_peak, t_v, r0)

    ok1 = check("pulse3(...) == make_pulse3(r0)(...)", np.array_equal(y_a, y_b))
    ok2 = check("pulse3(...) == pulse4(..., r=r0) (same point, both parametrizations)", np.allclose(y_a, y_c, rtol=1e-10))
    ok3 = check("make_pulse3(r0).mapping.r0 == r0", model.mapping.r0 == r0)
    return ok1 and ok2 and ok3


if __name__ == "__main__":
    results = [t1_round_trip(), t2_width_function(), overflow_safety(), api_consistency()]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
