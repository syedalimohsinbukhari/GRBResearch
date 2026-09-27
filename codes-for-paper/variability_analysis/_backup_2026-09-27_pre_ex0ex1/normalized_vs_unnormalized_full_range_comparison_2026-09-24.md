# Normalized-vs-unnormalized full-range Norris fit: current per-burst comparison (2026-09-24)

Snapshot of the `GRB{080916C,131014215,140206275,231129C}/norris_fit_results_*_{normalized,unnormalized}.csv`
comparison as of the latest P0 fixes (commit `[main-minor-98]` plus subsequent working-tree iteration on
GRB140206B). Companion to `variability_analysis.md`'s "Session update, 2026-09-24" section and
`t_peak_tv_reparametrization_idea.md` — this file is the numeric record; those two carry the narrative and
the theory, respectively.

**Re-verified 2026-09-26 01:53 PKT — GRB140206B's numbers below are confirmed reproducible, after fixing a
regression.** Commit `a57d164` (2026-09-25, after this doc was first written) re-enabled one of
`GRB140206275/_common.py`'s two `# replaceable` P0 pulses and reran only `fitter_normalized.py`, leaving
the normalized CSV at 8 pulses against the unnormalized CSV's stale 7 — no valid comparison was sitting on
disk. Re-commenting that pulse reproduced this doc's own "Stage 2" 7-pulse numbers (pulse 5/TR3's `t_v`
gap: 31.5% here vs. 31.4% below — same fit, confirmed independently). The *second* `# replaceable` pulse
(`(0.2, 15, 1, 1)`) turned out to have never actually been committed as commented-out anywhere in this
file's git history, despite this doc describing that 6-pulse state as "current" — the Stage 3 P0 reduction
had apparently only ever been applied on the machine that produced these numbers, not synced back here.
Commenting out both `# replaceable` pulses and rerunning both fitters reproduced the GRB140206B table
below bit-for-bit. `_common.py` now carries a comment recording that both must stay commented together.

**CLOSED, 2026-09-26 02:14 PKT — user decision.** Following the re-verification above, the user asked
for a third configuration to be tried: re-enable the near-t=0 pulse instead of the t=15s one (7 pulses
either way). That configuration ("Q0", `t_v` worst gap 2.15%) beat both the original 7-pulse
configuration ("P0", 31.5%) and the in-between 6-pulse one (5.58%) — see the GRB140206B section below
for the full numeric comparison. `_common.py` now defines `P0` and `Q0` as two permanent, fully-explicit
arrays (no more comment-toggling one shared list, which is what caused both prior desyncs). Both
configurations' fits are kept on disk under distinct filenames for future reference: `Q0` writes the
canonical (unsuffixed) files via `fitter_{normalized,unnormalized}.py`; `P0` writes `*_p0`-suffixed
files via `fitter_{normalized,unnormalized}_p0.py`. **Q0 is the definitive result. This closes the
normalized-vs-unnormalized comparison task for all four bursts — no open items remain anywhere in this
file.**

**Method reminder**: "normalized" = production convention, peak-normalize the light curve then fit with
`NorrisFitter`/pymultifit. "unnormalized" = fit the raw counts/s directly via
`scipy.optimize.least_squares` with an explicit `x_scale` array, bypassing pymultifit. Both use the same
data, same P0 (up to a unit conversion on amplitude), same bounds, same full `x.min()/x.max()` fit window.
Covariance for `t_v`'s MC propagation comes from the same kind of estimator for both methods (the
unnormalized fit's covariance is derived from its Jacobian the way `scipy.optimize.curve_fit` does
internally). Full derivation and rationale: `GRB080916C/fitter_unnormalized.py`'s docstring.

**How to read `kept`**: `mc_kept_fraction` from `tv_mc_summary()` — the fraction of MC-drawn covariance
samples that landed on physically valid (`tau1>0, tau2>0`) draws. Lower means the pulse is less tightly
constrained by the fit; the normalized/unnormalized pair should be read together with this, not `t_peak`/
`t_v` percent-differences alone, since a low-`kept` pulse is expected to show more spread between methods.

---

## GRB080916C — 6 pulses, clean

| pulse | episode | t_peak (norm / unnorm) | Δ | t_v (norm / unnorm) | Δ | kept (norm / unnorm) |
|---|---|---|---|---|---|---|
| 1 | — | 0.3342 / 0.3337 | 0.13% | 0.5606 / 0.5569 | 0.68% | 0.764 / 0.767 |
| 2 | TR1 | 2.6093 / 2.6087 | 0.02% | 4.0349 / 4.0351 | 0.01% | 0.999 / 0.998 |
| 3 | TR2 | 5.8496 / 5.8498 | 0.00% | 0.4337 / 0.4336 | 0.03% | 0.696 / 0.688 |
| 4 | TR3 | 25.9076 / 25.9073 | 0.00% | 15.3921 / 15.4009 | 0.06% | 1.000 / 1.000 |
| 5 | TR4 | 57.7057 / 57.7424 | 0.06% | 1.1973 / 1.1766 | 1.73% | 0.542 / 0.480 |
| 6 | TR5 | 62.6540 / 62.6521 | 0.00% | 1.6007 / 1.5991 | 0.10% | 0.870 / 0.867 |

Pulse 5/TR4 is the only pulse with `kept < 0.7` in either method, and correspondingly the largest `t_v`
spread (1.73%) — still small in absolute terms. Everything else is sub-0.1% on `t_v`. No open items here.

## GRB131014A — 5 pulses, clean (previously the "unresolved pulses 1/2" burst — now fixed)

| pulse | episode | t_peak (norm / unnorm) | Δ | t_v (norm / unnorm) | Δ | kept (norm / unnorm) |
|---|---|---|---|---|---|---|
| 1 | EX0 | 0.5099 / 0.5117 | 0.35% | 0.4575 / 0.4650 | 1.66% | 0.876 / 0.885 |
| 2 | EX0 | 1.2422 / 1.2416 | 0.05% | 0.1923 / 0.1932 | 0.46% | 0.602 / 0.621 |
| 3 | EX0 | 1.7596 / 1.7599 | 0.02% | 0.3745 / 0.3746 | 0.02% | 1.000 / 1.000 |
| 4 | TR2 | 2.6988 / 2.6990 | 0.01% | 0.5564 / 0.5564 | 0.00% | 1.000 / 1.000 |
| 5 | TR2 | 3.5926 / 3.5923 | 0.01% | 0.2890 / 0.2888 | 0.05% | 1.000 / 1.000 |

**Before the P0 fix**, pulses 1 and 2 (then labelled differently — no episode / TR1) disagreed badly
between methods: `t_peak` off by 128%/10%, `t_v` off by 4%/76%. That reproduced this burst's own
already-documented "unresolved" finding elsewhere in `variability_analysis.md` (two different seedings
already gave two different answers within the single-method production fit, `kept` 0.10–0.21 there).

**How it was actually fixed, per the user (2026-09-24):** the unnormalized (`least_squares`, explicit
`x_scale`) fits were consistently landing on better-constrained solutions than the normalized/pymultifit
path for this burst — plausibly because of how the `x_scale` fix conditions the optimizer for this
particular data, though the exact mechanism isn't nailed down. The fix was to take the unnormalized fit's
converged `t_s`/`tau1`/`tau2` (not `A` — amplitude isn't transferable between a peak-normalized and a
physical-units fit) as the new `P0` for the normalized fit, then hand-tune from there. Done on a Xeon
processor for the extra raw compute, which cut the iteration time for this burst from the ~6 hours spent
on it in the shared environment down to roughly 2 hours.

The P0 fix brought `kept` up to 0.60–0.88 for these two pulses and both methods now agree to ≤1.7% on
everything. Worth noting in case the original "unresolved" characterization needs revisiting upstream —
this P0 may have resolved a limitation that was previously written up as fundamental to the data.
**If GRB140206B's pulse 4/TR3 (still the one open item below) needs the same treatment, this is the
playbook: seed the normalized fit from its own unnormalized fit's converged shape parameters rather than
hand-guessing a fresh P0.**

## GRB140206B — RESOLVED 2026-09-26 02:14 PKT: Q0 (7-pulse-2) is the definitive result

**Closed, user decision.** Everything below this point (the 6-pulse "Stage 3" table and its pulse
4/TR3 open item) is superseded — kept as historical narrative, not the current status. The actual
resolution: `GRB140206275/_common.py` now defines two fully-explicit, non-overlapping P0 arrays
instead of toggling comment lines in one shared list (that toggling caused two real desyncs between
the normalized and unnormalized CSVs — commit `a57d164`, then a repeat on 2026-09-26 — see
`_common.py`'s own comments for the full incident history):

- **P0** — 7 pulses, near-t=0 pulse OUT, t=15s pulse IN. Historical record only
  (`fitter_{normalized,unnormalized}_p0.py`, `*_p0`-suffixed output files). Worst t_v gap: **31.5%**
  (pulse 5/TR3) — this is what the "Stage 2" table further down actually describes.
- **Q0** — 7 pulses, near-t=0 pulse IN, t=15s pulse OUT. **DEFINITIVE**
  (`fitter_{normalized,unnormalized}.py`, canonical unsuffixed output files). Worst t_v gap:
  **2.15%** (pulse 5/TR3) — better than P0 (31.5%) *and* better than the 6-pulse configuration
  below (5.58%), because dropping the t=15s pulse alone, with the near-t=0 pulse still available to
  absorb whatever flux sits at the window's start, lets pulse 5/TR3 settle into a shape both fitting
  methods agree on far more tightly than either P0 or the 6-pulse compromise:

| pulse | episode | t_peak (norm / unnorm) | Δ | t_v (norm / unnorm) | Δ | kept (norm / unnorm) |
|---|---|---|---|---|---|---|
| 1 | — | 0.0116 / 0.0116 | 0.01% | 0.5909 / 0.5908 | 0.02% | 0.506 / 0.504 |
| 2 | TR2 | 13.2038 / 13.2093 | 0.04% | 7.0849 / 7.1008 | 0.22% | 1.000 / 1.000 |
| 3 | TR2 | 13.9390 / 13.9394 | 0.00% | 1.8832 / 1.8811 | 0.11% | 1.000 / 1.000 |
| 4 | TR4 | 29.1286 / 29.1156 | 0.04% | 14.4947 / 14.5375 | 0.29% | 1.000 / 1.000 |
| 5 | TR3 | 23.8157 / 23.8460 | 0.13% | 1.0210 / 0.9990 | 2.15% | 0.800 / 0.800 |
| 6 | TR4 | 29.9013 / 29.9006 | 0.00% | 3.2328 / 3.2298 | 0.09% | 1.000 / 1.000 |
| 7 | TR6 | 122.0541 / 121.9736 | 0.07% | 22.6489 / 22.6272 | 0.10% | 0.997 / 0.998 |

Pulse 1 (the near-t=0 pulse) has a low `mc_kept_fraction` (~0.50, still weakly constrained in an
absolute sense) but the tightest normalized/unnormalized agreement of any pulse in any of the three
configurations tried (0.02% on `t_v`) — it is not adding noise to the fit, it's what lets pulse
5/TR3 resolve cleanly. **No open items remain for this burst.**

---

## GRB231129C — 6 pulses, clean (first full run — previously incomplete)

| pulse | episode | t_peak (norm / unnorm) | Δ | t_v (norm / unnorm) | Δ | kept (norm / unnorm) |
|---|---|---|---|---|---|---|
| 1 | EX0 | 0.6518 / 0.6516 | 0.03% | 0.8094 / 0.8088 | 0.08% | 0.996 / 0.996 |
| 2 | EX0 | 2.4657 / 2.4657 | 0.00% | 0.4816 / 0.4813 | 0.06% | 0.640 / 0.639 |
| 3 | EX0 | 1.3195 / 1.3197 | 0.01% | 0.6068 / 0.6072 | 0.06% | 0.989 / 0.985 |
| 4 | TR2 | 3.4902 / 3.4903 | 0.00% | 0.4331 / 0.4319 | 0.27% | 0.811 / 0.807 |
| 5 | TR2 | 4.2874 / 4.2874 | 0.00% | 0.8245 / 0.8242 | 0.04% | 1.000 / 1.000 |
| 6 | TR2 | 5.5075 / 5.5075 | 0.00% | 1.3471 / 1.3471 | 0.00% | 1.000 / 1.000 |

Cleanest of the four bursts — every pulse agrees to ≤0.3% on both `t_peak` and `t_v`. This is the burst
whose light curve/detector set changed mid-session (more NaI detectors than previously accounted for,
per `variability_analysis.md`); these numbers are from the resulting refit, now complete in both methods.

---

## Cross-burst summary

| burst | pulses | worst t_peak Δ | worst t_v Δ | open items |
|---|---|---|---|---|
| GRB080916C | 6 | 0.13% | 1.73% (pulse 5/TR4) | none |
| GRB131014A | 5 | 0.35% | 1.66% (pulse 1/EX0) | none — previously "unresolved" pulses 1/2 now fixed |
| GRB140206B | 7 (Q0) | 0.13% | 2.15% (pulse 5/TR3) | none — resolved 2026-09-26, see Q0 above |
| GRB231129C | 6 | 0.03% | 0.27% (pulse 4/TR2) | none |

**All four bursts are now clean by any reasonable threshold, with no open items remaining.**
GRB140206B needed the most iteration (three P0 configurations tried — see the historical section
below) before landing on Q0.

---

## Historical narrative, superseded by the above — kept for the record only

### GRB140206B — 6 pulses (COMPLEX model), one open item on pulse 4/TR3

| pulse | episode | t_peak (norm / unnorm) | Δ | t_v (norm / unnorm) | Δ | kept (norm / unnorm) |
|---|---|---|---|---|---|---|
| 1 | TR2 | 13.3713 / 13.4077 | 0.27% | 6.9376 / 6.9317 | 0.09% | 1.000 / 1.000 |
| 2 | TR2 | 13.9335 / 13.9240 | 0.07% | 1.8587 / 1.8484 | 0.55% | 1.000 / 1.000 |
| 3 | TR4 | 28.5912 / 28.8588 | 0.94% | 14.4642 / 14.4429 | 0.15% | 1.000 / 1.000 |
| 4 | TR3 | 23.8196 / 23.8486 | 0.12% | 1.2762 / 1.2049 | **5.58%** | 0.775 / 0.798 |
| 5 | TR4 | 29.9180 / 29.9170 | 0.00% | 3.3527 / 3.3261 | 0.79% | 1.000 / 1.000 |
| 6 | TR6 | 122.1343 / 121.9946 | 0.11% | 22.5193 / 22.5748 | 0.25% | 0.995 / 0.997 |

This is the only burst that needed real iteration this session — full history below, since the path here
matters for anyone re-deriving this or picking similar P0s for other bursts.

**Stage 1 (stale, do not use): an 8th pulse degenerated to a spike.** Before this burst's `P0` was pared
down, the unnormalized fit's CSV had 8 pulse rows against the normalized fit's 7 — a leftover from an
older `P0` version. The extra unnormalized pulse: `t_s=-0.0007, tau1=0.0001 (pinned at its 1e-4 lower
bound), t_peak=0.012s, kept=0.504` — a degenerate near-delta-function, the same failure mode already seen
and dropped for GRB080916C's old 7th pulse. Matching the *other* 7 unnormalized pulses to the normalized
7 by physical identity (not index) showed they agreed fine (≤2.5%) — the mismatch was entirely this one
spurious extra pulse, not a real disagreement. Resolved by simply re-running the unnormalized fit against
the then-current `P0`.

**Stage 2 (resynced, 7 pulses): pulse 5/TR3 showed a real 31% `t_v` gap.** Once both methods used the
same 7-pulse `P0`, every pulse's `t_peak` agreed to ≤1.1% — but pulse 5 (TR3, at the time) showed:

| param | normalized | unnormalized | Δ |
|---|---|---|---|
| `t_s` | 23.287 | 23.271 | 0.07% |
| `tau1` | 0.0429 | 0.0840 | **95.9%** |
| `tau2` | 8.402 | 5.175 | **38.4%** |
| `t_peak` | 23.887 | 23.930 | 0.18% |
| `t_v` | 3.321 | 2.276 | **31.4%** |
| `kept` | 0.839 | 0.899 | — |

Classic `tau1`↔`tau2` sloppy-parameter behavior (`tau1` almost doubles, `tau2` drops by over a third,
`t_s`/`t_peak` barely move) — except here, unusually, the degeneracy leaked into `t_v` itself rather than
staying confined to the raw parameters. `tau1*tau2` (which sets `t_peak`) was 0.360 vs 0.435 between the
two methods — close but not exactly conserved, meaning this pulse sits *near* but not *on* a pure
product-preserving ridge, with enough slack to drag `t_v` along too. **This is a real counterexample to
`t_peak_tv_reparametrization_idea.md`'s working assumption that `t_v` is reliably stable** — worth citing
there if that write-up gets revisited.

**Stage 3 (current, 6 pulses): dropping a pulse fixed most of it.** `P0` was reduced from 7 to 6 pulses
(one pulse removed — the exact one isn't recorded in this file, check `GRB140206275/_common.py`'s git
history/current `P0` for the specific parameters dropped). Pulse indices shifted, but matching by physical
identity: the pulse that inherited the old TR3 role (now pulse 4) improved from a 31% to a **5.58%** `t_v`
gap — still the single largest disagreement in this burst's table, and `kept` (0.775/0.798) is still the
only sub-1.0 value in the whole set, but this now reads as an ordinary residual rather than an open
problem. **This is the current, in-sync state** (both CSVs postdate `_common.py`) and the numbers in the
table at the top of this section are Stage 3, not Stage 1 or 2.

**Takeaway for pulse 4/TR3**: if this pulse's `t_v` ever needs to be trusted precisely (e.g. feeding
`Gamma_min`), treat 5.58% as its current normalized-vs-unnormalized systematic, and consider a
window-sensitivity check (same diagnostic already used elsewhere in this project) before relying on either
method's point value in isolation.

