This is solid work — the spec held up, and the two findings you surfaced along the way (MC-median `t_v_s` in Section 6, and the seed-insensitive attractor in T4b) are exactly the kind of thing a progress file should capture. The T4b observation deserves emphasis: for a smooth landscape, a bounded TRF optimizer converges to the *same* basin from any seed, so "reproducible fit" and "correct fit" decouple. Reproducibility was your original symptom in the 4-param fits, and it was the *mild* symptom. The audit machinery you've built (edge pinning, Δχ², window-widening) is what actually guards against the dangerous case, which is a confidently wrong interior attractor.

A few issues worth addressing before the real-data run, in rough priority order:

## 1. The sharp-profile sys = 0 is a grid-resolution artifact

In Section 4, a sharp profile yields a single-point flat zone, so `t_peak_err_sys = t_v_err_sys = 0.0` *exactly*. That conflates "the r0 systematic is small" with "the 25-point grid was too coarse to resolve the basin." A sharp minimum has genuine curvature-derived r0 uncertainty: the χ² profile near its minimum defines an r interval (Δχ² = 1 for a 1σ-flavored zone), and t_peak/t_v vary across that interval. Quantizing that to zero makes the one case where the data *does* know r report an artificially clean error bar.

Cheap fix: after `select_r0` picks r0 on a sharp profile, do a local refinement scan — ~10 points log-spaced in a factor-of-few around the minimum — and compute the systematic from the Δχ² ≤ 1 zone of the refined profile. Even better, use Δχ² = 1 for the systematic zone in *both* cases (flat and sharp), keeping +10 only for the flat/sharp classification where a wide "no meaningful difference" band is intended. As a bonus, the refined minimum is also a better r0 than the coarse-grid argmin (your sharp test landed at 215 vs. truth 190; fine, but refinement would tighten it).

## 2. Unweighted SSE must not reach the decision logic on real data

`chi_square()` falls back to plain SSE when `sigma=None`, and the "+10" flat-zone threshold and the Δχ² = 10 "dead weight" verdict both implicitly assume a χ² statistic. On unweighted SSE, those absolute thresholds are meaningless — the verdict scales with the arbitrary units of the rate. Your tests pass sigma everywhere, and the real pipeline has `data_rate_error` available, so this is latent rather than active. Still: make `sigma` mandatory in `select_r0`/`profile_scan` (raise rather than default to SSE), or at minimum have the plot and summary carry a flag saying which statistic the thresholds refer to. One silent unweighted run on real data would invalidate every downstream verdict in the file.

## 3. Audit the 4-param comparison fit's edges too

T4b's severely-truncated case ran all 9 seeds to the r upper bound (1e6). In `delta_chi2_report`, if the best-of-seeds 4-param fit *is* that runaway, the Δχ² comparison is still computed and the verdict machinery still runs — comparing a box-pinned 4-param solution against the 3-param fit. The verdict ("dead weight" vs "keep r free") is about parameter necessity, so a pinned 4-param fit usually loses the comparison anyway, but the cleaner behavior is to flag it: extend `check_edge_pinning` to also report whether the comparison 4-param fit's r (or log r) sat on its bounds, and treat the Δχ² verdict as "inconclusive, 4-param comparison pinned" in that case. Cheap, and it closes the loop on the exact pathology that motivated this whole project.

## 4. Guard the archived-median path

`select_r0` takes `archived_r_median` on faith. If the archived median falls *outside* the measured flat zone — possible once you compute it from the 112 archived CSV rows, since those come from 4-param fits that may themselves be degenerate — the code silently fixes r0 at a point the current pulse's own profile says is measurably worse. Add an assertion or warning when `archived_r_median` is outside the flat zone, and fall back to the flat-zone geometric center. Related: when you compute that archived r distribution from the CSVs, look at its shape — the spread and clustering of archived r across pulses is precisely the evidence for or against the shared-r hierarchical variant, and it's a two-line pandas call on files you already have.

## 5. For the real-GRB end-to-end run, mind the comparison baseline

When you compare new 3-param results against archived fits, remember your own Section 6 finding: archived `t_v_s` is the MC-propagated median, not the point estimate, and the two differ by O(1σ) on exactly the degenerate pulses you care about. Comparing your new point-estimate t_v against archived MC medians will manufacture disagreements on the fragile pulses. Either propagate your new fits through `tv_mc_summary()`-equivalent sampling for the comparison, or compare against archived *point* t_v (`tv_value(τ1, τ2)` at archived parameters — which you validated to be available) rather than `t_v_s`.

## On sequencing

The single-pulse synthetic validation is complete and trustworthy; the decisive evidence now is one real multi-pulse GRB through profile scan → finalize → audit, which is your stated next candidate. Two things I'd watch there that synthetics can't show you. First, neighbor overlap: a neighbor's tail under the pulse both flattens the r profile further and makes the per-pulse scan window-arbitrary — run the scan on the current-best joint model's residuals for that pulse (subtract other pulses, scan the target on what remains), and iterate once if the joint fit moves things. Second, window choice interacts with the truncation findings in T4b: on real data, if widening the window keeps pulling t_v or the flat zone around, that's the same degeneracy surfacing through window selection, and the window-widening audit item is the right detector.

None of the five issues above are design flaws — they're tightening of machinery that already works. Items 2 and 3 are the two I'd fix before the real-data run, since both concern the trustworthiness of the verdicts the run will produce.

---

## Response (2026-09-25 01:59 PKT)

Items 1-4 fixed in code; item 5 is applied guidance for the GRB231129C run rather than a code
change. Full detail (what changed, what was verified, exact numbers before/after) is in
`PROGRESS.md`'s "External review response" entry; a durable summary of the resulting design
decisions is in `norris_3param_reduction.md`'s "Every judgement call" section. One thing found
*while* fixing item 3, not flagged in the original review: the default 4-param comparison seed
range never actually reached near `r_bounds`, so the new pinned-comparison check could never have
fired regardless of whether a pinned solution existed — widened alongside the fix.

All 9 test/check scripts in this folder pass after the fixes, re-verified deterministic across
separate process runs.