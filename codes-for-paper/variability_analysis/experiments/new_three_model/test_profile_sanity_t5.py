"""Section 7 T5 (norris_3param_spec.md): "the chi2 profile from section 5 on the T4 simulation
must be flat-bottomed (best 3-param chi2 within ~1 of the best multi-seed 4-param chi2)."

Runs the full section-5 select_r0() decision rule (profile_scan.py) on both T4 datasets used in
test_multistart_t4.py -- the moderately-truncated one (T4a) and the more severely truncated one
(T4b) -- and checks both the FLAT classification and the tight Delta-chi2. Both datasets come
from synthetic_datasets.py, seeded via this project's seed_from_name convention (not a literal
int), specifically so T4 and T5 share the identical noise realization by construction.
"""
from plot_profile_scan import save_profile_scan_outputs
from profile_scan import select_r0
from synthetic_datasets import moderate_truncation_scenario, severe_truncation_scenario


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    return ok


def t5_on_dataset(name, t_window, y_window, sigma, dt, delta_chi2_tol=1.0):
    print(f"\n[T5: {name}]")
    result = select_r0(t_window, y_window, sigma=sigma, dt=dt)
    c = result.classification
    d = result.delta_chi2
    print(f"     shape={c.shape}  flat_zone={c.flat_zone}  span=x{c.flat_zone_span:.2f}")
    print(f"     chi2_3param(best r0)={d['chi2_3param']:.4f}  chi2_4param(best)={d['chi2_4param']:.4f}  delta_chi2={d['delta_chi2']:.4f}  r_pinned={d['r_pinned']}")

    paths = save_profile_scan_outputs(result, out_dir=".", label=f"section7_T5_{name.split()[0]}")
    print(f"     wrote {paths['grid_csv'].name}, {paths['summary_csv'].name}, {paths['fig_pdf'].name}, {paths['fig_png'].name}")

    ok1 = check(f"[{name}] classified as flat-bottomed", c.shape == "flat")
    ok2 = check(f"[{name}] |delta_chi2| within ~{delta_chi2_tol}", abs(d["delta_chi2"]) <= delta_chi2_tol)
    return ok1 and ok2


if __name__ == "__main__":
    results = [
        t5_on_dataset("moderate truncation (T4a dataset)", *moderate_truncation_scenario()),
        t5_on_dataset("severe truncation (T4b dataset)", *severe_truncation_scenario()),
    ]
    print()
    if all(results):
        print("ALL CHECKS PASSED")
    else:
        print("CHECKS FAILED -- see above")
        raise SystemExit(1)
