"""GRB131014215 blind multi-pulse discovery, v2: LOCALIZED greedy search.

greedy_discovery.py (v1, global-window) found only 1 pulse: with no local windowing, the
unconstrained single-pulse fit on this densely-packed (5-7 pulses in 3.2s) residual converged to
a broad ~1.5s-wide "average" envelope spanning several real sub-peaks at once, rather than
isolating any one of them -- visible directly in fit_overlay_GRB131014215_blind_discovery.png
(the fit completely misses the t=1.79s spike and the t~2.4s dip). Still fully blind (no archived
P0/fitted-parameter file read for this burst) -- the fix is algorithmic, not information leakage:

  1. scipy.signal.find_peaks on the raw light curve locates candidate LOCATIONS (indices), with a
     minimum height/prominence and minimum separation -- purely a property of the data, found by
     an algorithm, not eyeballed from a plot.
  2. Candidates are processed tallest-first. Each gets a LOCAL sub-window sized to the midpoint
     between it and its nearest still-unprocessed neighboring candidate (clipped to T90 and to a
     sane minimum/maximum width) -- so the fit can no longer "see" past a neighboring peak and
     smear across it.
  3. Same acceptance bars as v1 (Delta-chi2, amplitude sigma, multi-start reproducibility, no
     resolution-floor pinning) run on that LOCAL window's residual.
  4. If accepted, the fitted pulse is evaluated and subtracted over the FULL T90 window (not just
     the local sub-window -- Norris tails extend beyond it), then the next candidate is processed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

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

T90 = (0.960, 4.160)
WINDOW = T90
MIN_LOCAL_WIDTH = 0.35  # s, floor on the local sub-window half-width (a few x typical dt=0.064)
MAX_LOCAL_WIDTH = 1.2  # s, cap so an isolated candidate doesn't just re-grab the whole window
MIN_DELTA_CHI2 = 10.0
MIN_AMPLITUDE_SIGMA = 3.0
FIND_PEAKS_MIN_HEIGHT_SIGMA = 5.0
FIND_PEAKS_MIN_DISTANCE_S = 0.15
SEED = seed_from_name(__file__)


def local_window(t_candidate, other_candidate_times, dt):
    """Local sub-window around t_candidate, extending toward the nearest unresolved neighbor on
    each side (halfway), clipped to [MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH] half-width and to WINDOW."""
    left_neighbors = [t for t in other_candidate_times if t < t_candidate]
    right_neighbors = [t for t in other_candidate_times if t > t_candidate]
    half_left = (t_candidate - max(left_neighbors)) / 2 if left_neighbors else MAX_LOCAL_WIDTH
    half_right = (min(right_neighbors) - t_candidate) / 2 if right_neighbors else MAX_LOCAL_WIDTH
    half_left = float(np.clip(half_left, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    half_right = float(np.clip(half_right, MIN_LOCAL_WIDTH, MAX_LOCAL_WIDTH))
    lo = max(WINDOW[0], t_candidate - half_left)
    hi = min(WINDOW[1], t_candidate + half_right)
    if hi - lo < 2 * dt * 3:  # need enough bins for a real fit
        mid = (lo + hi) / 2
        lo, hi = max(WINDOW[0], mid - 3 * dt * 1.5), min(WINDOW[1], mid + 3 * dt * 1.5)
    return lo, hi


def main():
    print("=== GRB131014215: blind LOCALIZED greedy multi-pulse discovery (T90-only window) ===")
    t_full, y_full, sigma_full = load_light_curve()
    mask = (t_full >= WINDOW[0]) & (t_full <= WINDOW[1])
    t_w, y_w, sigma_w = t_full[mask], y_full[mask], sigma_full[mask]
    dt = median_dt(t_w)
    print(f"window: {WINDOW} ({t_w.size} bins, dt={dt:.4f})")

    peak_idx, props = find_peaks(
        y_w, height=FIND_PEAKS_MIN_HEIGHT_SIGMA * np.median(sigma_w), distance=max(1, int(FIND_PEAKS_MIN_DISTANCE_S / dt))
    )
    candidate_times = t_w[peak_idx].tolist()
    candidate_heights = y_w[peak_idx].tolist()
    order = np.argsort(candidate_heights)[::-1]  # tallest first
    candidate_times_ordered = [candidate_times[i] for i in order]
    print(f"find_peaks located {len(candidate_times)} candidate locations (>{FIND_PEAKS_MIN_HEIGHT_SIGMA} sigma, "
          f">={FIND_PEAKS_MIN_DISTANCE_S}s apart): {[f'{t:.3f}' for t in sorted(candidate_times)]}")

    residual = y_w.copy()
    accepted = []
    remaining = list(candidate_times_ordered)
    rng = get_rng(seed=SEED)

    while remaining:
        t_c = remaining.pop(0)
        others = remaining  # not-yet-processed candidates define this one's local window
        lo, hi = local_window(t_c, others, dt)
        lmask = (t_w >= lo) & (t_w <= hi)
        t_l, res_l, sigma_l = t_w[lmask], residual[lmask], sigma_w[lmask]
        print(f"\n--- candidate at t={t_c:.4f}, local window [{lo:.4f}, {hi:.4f}] ({t_l.size} bins) ---")

        if t_l.size < 6:
            print("-> SKIPPED: too few bins in local window")
            continue

        try:
            selection = select_r0(t_l, res_l, sigma=sigma_l, dt=dt)
            reported = finalize_pulse(t_l, res_l, selection, sigma=sigma_l, dt=dt)
        except Exception as e:  # noqa: BLE001
            print(f"-> SKIPPED: fit failed ({e})")
            continue

        model = make_pulse3(reported.r0)
        y_candidate_local = model(t_l, reported.amplitude, reported.t_peak, reported.t_v)
        chi2_before = chi_square(res_l, np.zeros_like(res_l), sigma=sigma_l)
        chi2_after = chi_square(res_l, y_candidate_local, sigma=sigma_l)
        delta_chi2 = chi2_before - chi2_after
        a_sigma = reported.amplitude / reported.amplitude_err if reported.amplitude_err > 0 else 0.0

        try:
            report = run_audit(t_l, res_l, selection, reported, rng, sigma=sigma_l, dt=dt)
            reliability_ok = report.items[0].passed and report.items[1].passed
        except Exception as e:  # noqa: BLE001
            print(f"-> audit failed ({e}); treating as unreliable")
            reliability_ok = False

        print(f"fit: A={reported.amplitude:.1f}+/-{reported.amplitude_err:.1f}  "
              f"t_peak={reported.t_peak:.4f}+/-{reported.t_peak_err:.4f}  "
              f"t_v={reported.t_v:.4f}+/-{reported.t_v_err:.4f}  r0={reported.r0:.4g} ({selection.classification.shape})")
        print(f"delta_chi2={delta_chi2:.2f}  amplitude_sigma={a_sigma:.2f}  reliability_ok={reliability_ok}")

        significance_ok = delta_chi2 > MIN_DELTA_CHI2 and a_sigma > MIN_AMPLITUDE_SIGMA
        if significance_ok and reliability_ok:
            print("-> ACCEPTED")
            accepted.append(
                {
                    "seed_t": t_c,
                    "local_window": (lo, hi),
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
            y_candidate_full = model(t_w, reported.amplitude, reported.t_peak, reported.t_v)
            residual = residual - y_candidate_full
        else:
            print("-> REJECTED (not added back to the queue)")

    print(f"\n=== Result: {len(accepted)} pulses accepted (blind, localized, no prior P0) ===")
    df = pd.DataFrame(accepted)
    if len(df):
        df = df.sort_values("t_peak_s").reset_index(drop=True)
        df.insert(0, "pulse_index", range(1, len(df) + 1))
        print(df[["pulse_index", "t_peak_s", "t_peak_err_s", "t_v_s", "t_v_err_s", "amplitude", "shape", "delta_chi2", "amplitude_sigma"]].to_string(index=False))
        csv_path = HERE / "GRB131014215_blind_discovery_local.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nwrote {csv_path.name}")

        y_total_fit = np.zeros_like(y_w)
        for row in accepted:
            model = make_pulse3(row["r0"])
            y_total_fit += model(t_w, row["amplitude"], row["t_peak_s"], row["t_v_s"])
        paths = plot_fit_overlay(
            t_w, y_w, y_total_fit, out_dir=str(HERE), label="GRB131014215_blind_discovery_local", sigma=sigma_w,
            extra_title=f"{len(accepted)} pulses found blind (localized), T90-only window",
        )
        print(f"wrote {paths['csv'].name}, {paths['pdf'].name}, {paths['png'].name}")

    print("\n=== Compared to the hint (5 currently fitted, potentially 7 total, NOT told the values) ===")
    print(f"blind pipeline found: {len(accepted)}")


if __name__ == "__main__":
    main()
