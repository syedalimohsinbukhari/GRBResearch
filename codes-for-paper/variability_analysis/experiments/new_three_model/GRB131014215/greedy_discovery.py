"""GRB131014215 blind multi-pulse discovery -- "blank slate test for the pipeline" (user,
2026-09-25): no archived P0, fitted parameters, or episode bounds for this burst were read to
build this. Only: the raw light curve (load_data.py) and T90 = [0.960, 4.160]s, a catalog/
observational duration (results.json's top-level "T90" key), not a fitted Norris-pulse parameter.
User: "stay inside T90 duration" (for now) -- WINDOW below is exactly T90, no padding.

Method: greedy iterative decomposition, the "from-scratch" strategy discussed but not built until
now (norris_3param_reduction.md's "Known limitations" -- this pipeline previously always needed a
prior decomposition to bootstrap a residual). Each iteration:
  1. Neutral-seed (argmax, no P0) a candidate pulse on the CURRENT residual.
  2. Run the full section 5/4/8 pipeline on it (profile scan -> finalize -> audit).
  3. Accept only if it clears real significance/reliability bars (Delta-chi2, amplitude sigma,
     multi-start reproducibility, no resolution-floor pinning) -- the same bars T4b's box-edge
     runaway and the GRB231129C hidden-peak test used to reject an unreliable candidate.
  4. If accepted, subtract it from the residual and repeat; if rejected, stop.

Caveat (flagged honestly, not hidden): WINDOW = T90 exactly means the window's own right edge
(4.160s) cuts into real, still-significant decay flux (~37 sigma at t=4.16 in the raw light
curve) -- the same truncation risk GRB231129C's pulse 6 hit with too narrow a window. Kept as
instructed ("for now"); the true pulse count/shapes here should be re-checked with a properly
padded window before treating this as final.
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
from load_data import load_light_curve  # noqa: E402
from plot_diagnostics import plot_fit_overlay  # noqa: E402
from profile_scan import chi_square, select_r0  # noqa: E402
from pulse3 import make_pulse3  # noqa: E402
from report_pulse import finalize_pulse  # noqa: E402

T90 = (0.960, 4.160)  # results.json, catalog value -- not a fitted parameter
WINDOW = T90  # user: "stay inside T90 duration" (for now)
MAX_CANDIDATES = 10
MIN_DELTA_CHI2 = 10.0  # this project's usual "the parameters earned their keep" bar
MIN_AMPLITUDE_SIGMA = 3.0
SEED = seed_from_name(__file__)


def main():
    print("=== GRB131014215: blind greedy multi-pulse discovery (T90-only window, no prior P0) ===")
    t_full, y_full, sigma_full = load_light_curve()
    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
    dt = median_dt(t_w)
    print(f"window: {WINDOW} ({t_w.size} bins, dt={dt:.4f})")
    print(f"raw peak significance in window: {(y_w/sigma_w).max():.1f} sigma")

    residual = y_w.copy()
    accepted = []
    rng = get_rng(seed=SEED)

    for k in range(1, MAX_CANDIDATES + 1):
        print(f"\n--- candidate {k}: neutral seed on current residual ---")
        idx_peak = int(np.argmax(residual))
        print(f"argmax(residual) at t={t_w[idx_peak]:.4f}, value={residual[idx_peak]:.1f} "
              f"({residual[idx_peak]/sigma_w[idx_peak]:.2f} sigma)")

        selection = select_r0(t_w, residual, sigma=sigma_w, dt=dt)
        reported = finalize_pulse(t_w, residual, selection, sigma=sigma_w, dt=dt)

        model = make_pulse3(reported.r0)
        y_candidate = model(t_w, reported.amplitude, reported.t_peak, reported.t_v)
        chi2_before = chi_square(residual, np.zeros_like(residual), sigma=sigma_w)
        chi2_after = chi_square(residual, y_candidate, sigma=sigma_w)
        delta_chi2 = chi2_before - chi2_after
        a_sigma = reported.amplitude / reported.amplitude_err if reported.amplitude_err > 0 else 0.0

        report = run_audit(t_w, residual, selection, reported, rng, sigma=sigma_w, dt=dt)

        print(f"fit: A={reported.amplitude:.1f}+/-{reported.amplitude_err:.1f}  "
              f"t_peak={reported.t_peak:.4f}+/-{reported.t_peak_err:.4f}  "
              f"t_v={reported.t_v:.4f}+/-{reported.t_v_err:.4f}  r0={reported.r0:.4g} ({selection.classification.shape})")
        print(f"delta_chi2={delta_chi2:.2f}  amplitude_sigma={a_sigma:.2f}")
        for item in report.items:
            if item.name in ("no parameter on a box edge", "multi-start reproducibility"):
                print(f"  [{'PASS' if item.passed else 'FAIL'}] {item.name}: {item.detail}")

        reliability_ok = report.items[0].passed and report.items[1].passed  # edge-pinning, reproducibility
        significance_ok = delta_chi2 > MIN_DELTA_CHI2 and a_sigma > MIN_AMPLITUDE_SIGMA

        if significance_ok and reliability_ok:
            print(f"-> ACCEPTED as pulse {len(accepted) + 1}")
            accepted.append(
                {
                    "candidate_order": k,
                    "amplitude": reported.amplitude,
                    "amplitude_err": reported.amplitude_err,
                    "t_peak_s": reported.t_peak,
                    "t_peak_err_s": reported.t_peak_err,
                    "t_v_s": reported.t_v,
                    "t_v_err_s": reported.t_v_err,
                    "r0": reported.r0,
                    "shape": selection.classification.shape,
                    "delta_chi2": delta_chi2,
                    "amplitude_sigma": a_sigma,
                }
            )
            residual = residual - y_candidate
        else:
            reasons = []
            if not significance_ok:
                reasons.append(f"significance (delta_chi2={delta_chi2:.2f}<{MIN_DELTA_CHI2} or sigma={a_sigma:.2f}<{MIN_AMPLITUDE_SIGMA})")
            if not reliability_ok:
                reasons.append("reliability (edge-pinned or non-reproducible)")
            print(f"-> REJECTED: fails {' and '.join(reasons)}. Stopping search.")
            break

    print(f"\n=== Result: {len(accepted)} pulses accepted (blind, no prior P0) ===")
    df = pd.DataFrame(accepted)
    if len(df):
        df.insert(0, "pulse_index", range(1, len(df) + 1))
        df = df.sort_values("t_peak_s").reset_index(drop=True)
        df["pulse_index"] = range(1, len(df) + 1)
        print(df[["pulse_index", "t_peak_s", "t_peak_err_s", "t_v_s", "t_v_err_s", "amplitude", "shape", "delta_chi2", "amplitude_sigma"]].to_string(index=False))
        csv_path = HERE / "GRB131014215_blind_discovery.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nwrote {csv_path.name}")

        y_total_fit = np.zeros_like(y_w)
        for row in accepted:
            model = make_pulse3(row["r0"])
            y_total_fit += model(t_w, row["amplitude"], row["t_peak_s"], row["t_v_s"])
        paths = plot_fit_overlay(
            t_w, y_w, y_total_fit, out_dir=str(HERE), label="GRB131014215_blind_discovery", sigma=sigma_w,
            extra_title=f"{len(accepted)} pulses found blind, T90-only window",
        )
        print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print(f"\n=== Compared to the hint (5 currently fitted, potentially 7 total, NOT told the values) ===")
    print(f"blind pipeline found: {len(accepted)}")


if __name__ == "__main__":
    main()
