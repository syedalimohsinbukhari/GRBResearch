"""Paper-facing check, 2026-09-23: for all four bursts, does fitting the NORMALIZED (production-method)
light curve at its own full range (x.min(), x.max()) change t_peak/t_v by more than 5% relative to the
current production fit (each burst's own narrower window)? If not, per the user's explicit instruction,
the full range becomes the recommended window and its parameters are what should feed
lorentz_factor.py's Gamma_min pipeline -- not yet wired in, this script only checks and reports.

Normalized method only (NorrisFitter, peak-normalized y) -- the paper uses the normalized spectrum, so
the unnormalized/explicit-x_scale comparison in fit_comparison.py is not needed here; it already
established (comparison.md) that the two methods agree to <0.15% SSE and the normalized method is what
every fitter_*.py actually runs in production.

GRB140206B uses the pulse-6-reseeded P0 (A=0.2, t_s=23, tau1=1, tau2=1) at full range -- the plain
production P0 is already known (grb140206b_window_scan.py) to collapse pulse 6 there; testing the
unreseeded P0 again would just reproduce a known failure, not answer the paper-facing question.

Does not modify any fitter_*.py. Reads each burst's *current* production norris_fit_results_*.csv as the
comparison baseline (not re-derived from memory).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fit_comparison import LC_DIR, HERE, VARIABILITY_DIR, summed_nai_curve, fit_normalized, total_model
from norris_fit import t_peak, tv_value
from grb_research import update_style, LINE_WIDTH
from grb_research.grb_utils import save_fig

MAX_NFEV = 20000

# T05/T95 (T90 start/end) per burst, read from results.json's interval strings via the same convention
# already used in each fitter_*.py's own EPISODE_BOUNDS comments (T05 = the first TR interval's start,
# T95 = the last TR interval's end -- not re-derived, copied from those files' own sourced comments).
# Plot-only: bounds the x-axis to T05-25s .. T95+25s so the (much wider) full fit window doesn't make
# the actual burst illegible. The FIT itself still runs on the full x.min()/x.max() range regardless.
T90 = {
    "GRB080916009": (1.280, 64.256),
    "GRB131014215": (0.960, 4.160),
    "GRB140206275": (7.488, 154.240),
    "GRB231129779": (0.384, 7.296),
}
PLOT_PAD_S = 25.0

BURSTS = {
    "GRB080916009": dict(
        paper_name="GRB080916C",
        prod_window=(-1, 70),
        p0=[
            (0.4, -0.7, 2.34, 0.471),
            (0.7, 0.3, 0.88, 6),
            (0.2, 5.3, 0.6, 6),
            (0.3, 1.3, 43, 14),
            (0.3, 52, 18, 0.7),
            (0.3, 61, 0.3, 3),
        ],
        prod_csv="norris_fit_results_GRB080916009.csv",
    ),
    "GRB131014215": dict(
        paper_name="GRB131014A",
        prod_window=(-1, 10),
        p0=[
            (0.115, -0.9, 1, 1),
            (0.2, -0.8, 1, 1),
            (0.9, 1.19, 1, 1),
            (0.4, 2.4, 1, 1),
            (0.3, 2.71, 1, 1),
        ],
        prod_csv="norris_fit_results_GRB131014215.csv",
    ),
    "GRB140206275": dict(
        paper_name="GRB140206B",
        prod_window=(-1, 160),
        # Pulse 6 re-seeded (A=0.2, t_s=23, tau1=1, tau2=1) -- see grb140206b_fullrange_fitter.py /
        # comparison.md "Recommended configuration". The plain production P0 collapses pulse 6 at full
        # range (grb140206b_window_scan.py); this is the seed that doesn't.
        p0=[
            (0.13, -0.3, 0.12, 1.68),
            (0.35, 4.0, 6, 14),
            (0.5, 11, 5, 1.4),
            (0.5, 28, 2, 1),
            (0.23, 24, 0.3, 1.434),
            (0.2, 23, 1, 1),
            (0.05, -0.9, 3400, 4.4),
        ],
        prod_csv="norris_fit_results_GRB140206275.csv",  # COMPLEX model -- the one production uses for episode results
    ),
    "GRB231129779": dict(
        paper_name="GRB231129C",
        prod_window=(-1, 10),
        p0=[
            (0.6, -0.2, 1, 1),
            (0.3, 0.08, 1, 1),
            (0.6, 2, 0.5, 0.5),
            (0.6, 4, 0.5, 0.5),
            (0.3, 4.2, 2, 1),
        ],
        prod_csv="norris_fit_results_GRB231129779.csv",
    ),
}

THRESHOLD_PCT = 5.0


def pct_diff(new, old):
    return 100.0 * abs(new - old) / abs(old) if old else float("nan")


def main():
    all_rows = []
    for grb_dir, cfg in BURSTS.items():
        grb_dir_path = LC_DIR / grb_dir
        t_full, _, _ = summed_nai_curve(grb_dir_path, (-np.inf, np.inf))
        full_window = (float(t_full.min()), float(t_full.max()))
        n_pulses = len(cfg["p0"])

        t, y_raw, dat_nai = summed_nai_curve(grb_dir_path, full_window)
        y_max = float(np.max(y_raw))
        params = fit_normalized(t, y_raw, y_max, cfg["p0"], MAX_NFEV)
        model = total_model(t, params, n_pulses) * y_max
        sse = float(np.sum((model - y_raw) ** 2))

        prod_df = pd.read_csv(VARIABILITY_DIR / cfg["prod_csv"]).drop_duplicates(subset=["pulse_index"]).set_index("pulse_index")

        print(f"\n=== {cfg['paper_name']} ({grb_dir}) -- full range {full_window} vs production {cfg['prod_window']} ===")
        print(f"{'pulse':>5}  {'t_peak(full)':>12}  {'t_peak(prod)':>12}  {'%diff':>7}  {'t_v(full)':>10}  {'t_v(prod)':>10}  {'%diff':>7}  {'flag':>6}")
        for i in range(n_pulses):
            A, ts, tau1, tau2 = params[i * 4:(i + 1) * 4]
            tp_full, tv_full = t_peak(ts, tau1, tau2), tv_value(tau1, tau2)
            pulse_idx = i + 1
            if pulse_idx not in prod_df.index:
                print(f"{pulse_idx:5d}  (no matching production pulse index -- skipped)")
                continue
            tp_prod = prod_df.loc[pulse_idx, "t_peak_s"]
            tv_prod = prod_df.loc[pulse_idx, "t_v_s"]
            tp_pct = pct_diff(tp_full, tp_prod)
            tv_pct = pct_diff(tv_full, tv_prod)
            flag = "FLAG" if max(tp_pct, tv_pct) > THRESHOLD_PCT else "ok"
            print(f"{pulse_idx:5d}  {tp_full:12.4f}  {tp_prod:12.4f}  {tp_pct:6.2f}%  {tv_full:10.4f}  {tv_prod:10.4f}  {tv_pct:6.2f}%  {flag:>6}")
            all_rows.append({
                "grb_name": cfg["paper_name"], "grb_dir": grb_dir, "pulse_index": pulse_idx,
                "full_window_min_s": full_window[0], "full_window_max_s": full_window[1],
                "prod_window_min_s": cfg["prod_window"][0], "prod_window_max_s": cfg["prod_window"][1],
                "full_A_cts_per_s": A * y_max, "full_t_s": ts, "full_tau1": tau1, "full_tau2": tau2,
                "full_t_peak_s": tp_full, "full_t_v_s": tv_full,
                "prod_t_peak_s": tp_prod, "prod_t_v_s": tv_prod,
                "t_peak_pct_diff": tp_pct, "t_v_pct_diff": tv_pct,
                "flag_gt_5pct": max(tp_pct, tv_pct) > THRESHOLD_PCT,
                "sse_full_range": sse,
            })

        # --- Plot: fit computed on the full x.min()/x.max() range (above), but the AXIS is bounded to
        # T05-25s .. T95+25s -- plot legibility only, does not change what was fitted.
        t05, t95 = T90[grb_dir]
        xlim = (t05 - PLOT_PAD_S, t95 + PLOT_PAD_S)
        update_style()
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(t, y_raw, color="0.6", lw=LINE_WIDTH * 0.6,
                label=f"10-400 keV NaI ({'+'.join(dat_nai)}, summed)\nBackground subtracted")
        ax.plot(t, model, color="tab:blue", lw=LINE_WIDTH,
                label=f"Normalized fit total (fitted on full range {full_window[0]:.1f}, {full_window[1]:.1f}s)")
        ax.axvline(t05, color="0.3", ls=":", lw=1, label=f"T05={t05:.3f}s")
        ax.axvline(t95, color="0.3", ls="--", lw=1, label=f"T95={t95:.3f}s")
        ax.set_xlim(*xlim)
        ax.set_xlabel("Time since trigger [s]")
        ax.set_ylabel("Count rate [counts/s]")
        ax.set_title(f"{cfg['paper_name']}: full-range normalized fit (plot bounded to T05-25s, T95+25s)")
        ax.legend(fontsize="small")
        fig_path = HERE / f"full_range_check_{grb_dir}"
        save_fig(fig, fig_path)

    combined = pd.DataFrame(all_rows)
    csv_path = HERE / "full_range_all_bursts_check_results.csv"
    combined.to_csv(csv_path, index=False)
    print(f"\nwrote {csv_path} ({len(combined)} rows)")
    n_flagged = combined["flag_gt_5pct"].sum()
    print(f"\n{n_flagged} of {len(combined)} pulses flagged (>5% t_peak or t_v difference from production)")
    return combined


if __name__ == "__main__":
    main()
