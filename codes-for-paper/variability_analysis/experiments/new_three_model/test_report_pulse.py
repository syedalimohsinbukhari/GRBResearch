"""Validation for report_pulse.py against norris_3param_spec.md section 4, reusing the two
synthetic scenarios from test_profile_scan.py: a sharp-minimum pulse (small but non-zero
r0-choice systematic expected, from the section-4 Delta-chi2<=1 local refinement around r0_chosen
-- NOT exactly zero; see report_pulse.py's module docstring for why reusing the coarse section-5
flat zone verbatim was wrong for a sharp profile) and a flat-bottom pulse (large systematic
expected, dominating the statistical term, since the flat zone spans the whole grid and
r0_chosen sits far from the profile minimum's own grid neighborhood).
"""
from plot_diagnostics import plot_fit_overlay
from profile_scan import select_r0
from pulse3 import make_pulse3
from report_pulse import finalize_pulse
from synthetic_datasets import flat_scenario, sharp_scenario


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def sharp_case():
    print("\n[sharp] r_true=190, full window, low noise")
    r_true = 190.0
    a_true, t_peak_true, t_v_true = 20.0, 15.0, 1.2
    t_window, y_window, sigma, dt = sharp_scenario()  # shared dataset -- synthetic_datasets.py

    selection = select_r0(t_window, y_window, sigma=sigma, dt=dt)
    reported = finalize_pulse(t_window, y_window, selection, sigma=sigma, dt=dt)

    print(f"     r0_chosen={reported.r0:.3g}  flat_zone={reported.r0_flat_zone}")
    print(f"     A={reported.amplitude:.4f} +/- {reported.amplitude_err:.4f} (true {a_true})")
    print(
        f"     t_peak={reported.t_peak:.4f} +/- {reported.t_peak_err:.4f} "
        f"(stat={reported.t_peak_err_stat:.4f}, sys={reported.t_peak_err_sys:.4f}) (true {t_peak_true})"
    )
    print(
        f"     t_v={reported.t_v:.4f} +/- {reported.t_v_err:.4f} "
        f"(stat={reported.t_v_err_stat:.4f}, sys={reported.t_v_err_sys:.4f}) (true {t_v_true})"
    )
    print(f"     tau1={reported.tau1:.4f} +/- {reported.tau1_err:.4f}")
    print(f"     tau2={reported.tau2:.4f} +/- {reported.tau2_err:.4f}")
    print(f"     t_s={reported.t_s:.4f} +/- {reported.t_s_err:.4f}")

    model = make_pulse3(reported.r0)
    y_fit = model(t_window, reported.amplitude, reported.t_peak, reported.t_v)
    paths = plot_fit_overlay(
        t_window, y_window, y_fit, out_dir=".", label="section4_sharp", sigma=sigma,
        extra_title=f"r0={reported.r0:.3g}, t_peak={reported.t_peak:.3f}+/-{reported.t_peak_err:.3f}, t_v={reported.t_v:.3f}+/-{reported.t_v_err:.3f}",
    )
    print(f"     wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    ok1 = check("A recovered within 3-sigma of truth", abs(reported.amplitude - a_true) < 3 * reported.amplitude_err)
    ok2 = check("t_peak recovered within 3-sigma of truth", abs(reported.t_peak - t_peak_true) < 3 * reported.t_peak_err)
    ok3 = check(
        "refined r0-choice systematic is small (but non-zero) relative to the statistical term on a sharp profile",
        0.0 <= reported.t_peak_err_sys < 3 * reported.t_peak_err_stat and 0.0 <= reported.t_v_err_sys < 3 * reported.t_v_err_stat,
    )
    ok4 = check(
        "combined t_peak_err is at least as large as the statistical term alone (quadrature never shrinks it)",
        reported.t_peak_err >= reported.t_peak_err_stat,
    )
    return ok1 and ok2 and ok3 and ok4


def flat_case():
    print("\n[flat] r_true=190, window starts partway up the rise, more noise")
    r_true = 190.0
    a_true, t_peak_true, t_v_true = 20.0, 15.0, 1.2
    t_window, y_window, sigma, dt = flat_scenario()  # shared dataset -- synthetic_datasets.py

    selection = select_r0(t_window, y_window, sigma=sigma, dt=dt, archived_r_median=200.0)
    reported = finalize_pulse(t_window, y_window, selection, sigma=sigma, dt=dt)

    print(f"     r0_chosen={reported.r0:.3g}  flat_zone={reported.r0_flat_zone}")
    print(f"     A={reported.amplitude:.4f} +/- {reported.amplitude_err:.4f} (true {a_true})")
    print(
        f"     t_peak={reported.t_peak:.4f} +/- {reported.t_peak_err:.4f} "
        f"(stat={reported.t_peak_err_stat:.4f}, sys={reported.t_peak_err_sys:.4f}) (true {t_peak_true})"
    )
    print(
        f"     t_v={reported.t_v:.4f} +/- {reported.t_v_err:.4f} "
        f"(stat={reported.t_v_err_stat:.4f}, sys={reported.t_v_err_sys:.4f}) (true {t_v_true})"
    )
    print(f"     tau1={reported.tau1:.4f} +/- {reported.tau1_err:.4f}")
    print(f"     tau2={reported.tau2:.4f} +/- {reported.tau2_err:.4f}")

    model = make_pulse3(reported.r0)
    y_fit = model(t_window, reported.amplitude, reported.t_peak, reported.t_v)
    paths = plot_fit_overlay(
        t_window, y_window, y_fit, out_dir=".", label="section4_flat", sigma=sigma,
        extra_title=f"r0={reported.r0:.3g}, t_peak={reported.t_peak:.3f}+/-{reported.t_peak_err:.3f}, t_v={reported.t_v:.3f}+/-{reported.t_v_err:.3f}",
    )
    print(f"     wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    ok1 = check(
        "r0-choice systematic dominates over statistical term for a flat (whole-grid) profile",
        reported.t_v_err_sys > reported.t_v_err_stat,
    )
    ok2 = check("t_v_err (combined) is strictly larger than the statistical term alone", reported.t_v_err > reported.t_v_err_stat)
    ok3 = check(
        "t_peak still recovered reasonably despite the missing rise (informative from the decay alone)",
        abs(reported.t_peak - t_peak_true) < 2.0,
    )
    return ok1 and ok2 and ok3


if __name__ == "__main__":
    results = [sharp_case(), flat_case()]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
