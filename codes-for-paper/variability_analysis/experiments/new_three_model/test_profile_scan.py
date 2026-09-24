"""Validation for profile_scan.py against norris_3param_spec.md section 5, using the two
synthetic scenarios the spec itself describes (section 7 T4/T5): a well-resolved pulse (clean
rise+decay both visible -> data constrains r -> SHARP minimum) and a pulse whose window starts
partway up the rise (spec's T4 setup -> asymmetry unconstrained -> FLAT bottom).

Also exercises plot_profile_scan.py so this run leaves real CSV/PDF/PNG outputs in this folder,
per this project's generated-output convention.
"""
from plot_profile_scan import save_profile_scan_outputs
from profile_scan import select_r0
from synthetic_datasets import flat_scenario, sharp_scenario


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def sharp_minimum_scenario():
    print("\n[sharp] full pulse window, low noise, r_true=190 -> expect a sharp chi2(r0) minimum")
    r_true = 190.0
    t_window, y_window, sigma, dt = sharp_scenario()  # shared dataset -- synthetic_datasets.py

    result = select_r0(t_window, y_window, sigma=sigma, dt=dt)
    c = result.classification
    print(f"     shape={c.shape}  r0_at_min={c.r0_at_min:.3g} (true r={r_true})")
    print(f"     flat_zone={c.flat_zone}  span=x{c.flat_zone_span:.2f}")
    print(f"     chi2_min_3param={c.chi2_min:.3g}")
    if result.fit4 is not None:
        print(f"     best 4-param: chi2={result.fit4['chi2']:.3g} r={result.fit4['r']:.3g}")
    print(f"     delta_chi2={result.delta_chi2['delta_chi2']} -> {result.delta_chi2['verdict']}")
    print(f"     r0_chosen={result.r0_chosen:.3g} ({result.r0_choice_reason})")

    paths = save_profile_scan_outputs(result, label="synthetic_sharp")
    for k, v in paths.items():
        print(f"     wrote {k}: {v.name}")

    ok1 = check("classified as sharp minimum", c.shape == "sharp")
    ok2 = check("r0_at_min within a factor of 2 of the true r", 0.5 * r_true <= c.r0_at_min <= 2 * r_true)
    return ok1 and ok2


def flat_bottom_scenario():
    print("\n[flat] window starts partway up the rise, more noise, r_true=190 -> expect a flat chi2(r0) bottom")
    r_true = 190.0
    # window starts at t=14.5, just before the peak (t_peak=15) -- most of the rise (which happens
    # on the fast tau1/tau2=r_true=190 timescale) is cut off, matching spec T4's setup.
    t_window, y_window, sigma, dt = flat_scenario()  # shared dataset -- synthetic_datasets.py

    result = select_r0(t_window, y_window, sigma=sigma, dt=dt, archived_r_median=200.0)
    c = result.classification
    print(f"     shape={c.shape}  r0_at_min={c.r0_at_min:.3g} (true r={r_true})")
    print(f"     flat_zone={c.flat_zone}  span=x{c.flat_zone_span:.2f}")
    print(f"     chi2_min_3param={c.chi2_min:.3g}")
    if result.fit4 is not None:
        print(f"     best 4-param: chi2={result.fit4['chi2']:.3g} r={result.fit4['r']:.3g}")
    print(f"     delta_chi2={result.delta_chi2['delta_chi2']} -> {result.delta_chi2['verdict']}")
    print(f"     r0_chosen={result.r0_chosen:.3g} ({result.r0_choice_reason})")

    paths = save_profile_scan_outputs(result, label="synthetic_flat")
    for k, v in paths.items():
        print(f"     wrote {k}: {v.name}")

    ok1 = check("classified as flat bottom", c.shape == "flat")
    ok2 = check("flat zone spans a factor of >= 3 in r0", c.flat_zone_span >= 3.0)
    ok3 = check("r0_chosen uses the archived median (flat-bottom action)", result.r0_chosen == 200.0)
    return ok1 and ok2 and ok3


if __name__ == "__main__":
    results = [sharp_minimum_scenario(), flat_bottom_scenario()]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
