"""Section 7 T3 (norris_3param_spec.md): numerically locate the two half-max points of I(t) via
brentq and confirm (x_right - x_left)/2 == the width_function()-based t_v -- an independent,
root-finding check on norris_raw() that doesn't go through width_function() itself at all, so it
can't share a transcription bug with pulse3.py's own formula.
"""
import numpy as np
from scipy.optimize import brentq

from plot_diagnostics import plot_fit_overlay
from pulse3 import norris_raw, reduced_from_raw, width_function


def half_max_points(tau1: float, tau2: float, t_s: float = 0.0):
    """(x_left, x_right, half_value) via brentq on I(t) - I(t_peak)/2, independent of width_function()."""
    amplitude = 1.0
    t_peak = t_s + np.sqrt(tau1 * tau2)
    i_peak = norris_raw(np.array([t_peak]), amplitude, t_s, tau1, tau2)[0]
    half = i_peak / 2

    def f(t):
        return norris_raw(np.array([t]), amplitude, t_s, tau1, tau2)[0] - half

    # bracket the rise-side root in (t_s, t_peak), the decay-side root in (t_peak, far right)
    x_left = brentq(f, t_s + 1e-9, t_peak)
    x_right_hi = t_peak
    step = max(tau2, 1.0)
    while f(x_right_hi) > 0:
        x_right_hi += step
    x_right = brentq(f, t_peak, x_right_hi)

    return x_left, x_right, half


def half_max_width(tau1: float, tau2: float, t_s: float = 0.0) -> float:
    """(x_right - x_left) / 2, via half_max_points."""
    x_left, x_right, _ = half_max_points(tau1, tau2, t_s)
    return (x_right - x_left) / 2


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def t3_spec_worked_example():
    print("\n[T3] spec worked example: tau1=85.6, tau2=0.453 -> t_v = 1.40710035884524")
    expected = 1.40710035884524
    tau1, tau2, t_s = 85.6, 0.453, 3.0
    x_left, x_right, half = half_max_points(tau1, tau2, t_s)
    got = (x_right - x_left) / 2
    print(f"     brentq FWHM/2 = {got!r}")

    t_peak = t_s + np.sqrt(tau1 * tau2)
    t_plot = np.linspace(t_s - 0.05, t_peak + 6 * tau2, 2000)
    y_plot = norris_raw(t_plot, 1.0, t_s, tau1, tau2)
    paths = plot_fit_overlay(
        t_plot, y_plot, y_plot, out_dir=".", label="section7_T3_fwhm",
        half_max=(x_left, x_right, half), extra_title=f"tau1={tau1}, tau2={tau2}, FWHM/2={got:.6f}",
    )
    print(f"     wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    return check("matches spec's expected value to 1e-6", np.isclose(got, expected, atol=1e-6))


def t3_matches_width_function_across_r():
    print("\n[T3] brentq FWHM/2 vs width_function()-based t_v, across a range of r=tau1/tau2")
    all_ok = True
    for tau2, r in [(0.453, 189.0), (1.0, 0.01), (1.0, 1.0), (1.0, 1e4), (5.0, 50.0)]:
        tau1 = r * tau2
        got = half_max_width(tau1, tau2, t_s=-2.0)
        expected = tau2 * width_function(r)
        ok = np.isclose(got, expected, rtol=1e-6)
        print(f"     tau1={tau1:.4g} tau2={tau2:.4g} (r={r:.4g}): brentq={got!r} width_function={expected!r} -> {'PASS' if ok else 'FAIL'}")
        all_ok &= ok
    return check("brentq FWHM/2 matches width_function()-based t_v for all tested (tau1, tau2)", all_ok)


def t3_consistent_with_reduced_from_raw():
    print("\n[T3] cross-check against pulse3.reduced_from_raw's t_v for the same worked example")
    tau1, tau2, t_s = 85.6, 0.453, 3.0
    _, t_v, _ = reduced_from_raw(t_s, tau1, tau2)
    got = half_max_width(tau1, tau2, t_s=t_s)
    return check(f"brentq FWHM/2 ({got!r}) == reduced_from_raw's t_v ({t_v!r})", np.isclose(got, t_v, rtol=1e-8))


if __name__ == "__main__":
    results = [t3_spec_worked_example(), t3_matches_width_function_across_r(), t3_consistent_with_reduced_from_raw()]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
