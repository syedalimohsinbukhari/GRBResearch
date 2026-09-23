# Normalized vs. unnormalized Norris fitting — does peak-normalization change the answer?

Standalone diagnostic, 2026-09-23. Prompted by a question raised in the parent session while testing
the BUG-23 detector-summing fix (`../../variability_analysis.md`): every `fitter_*.py` in the parent
folder peak-normalizes the light curve (`y /= Y_MAX_CTS_PER_S`) before fitting. Is that step load-bearing
for the fit result, or just a numerical convenience that could be dropped?

**Everything in this folder is self-contained.** No file outside `experiments/normalized_vs_unnormalized_fit/`
was created or modified to produce it. `fit_comparison.py` reads `light_curves.py`/`norris_fit.py` from
the parent folder (read-only, `sys.path.insert`) and copies each burst's `P0`/fit-window as literals —
it does not import or execute any `fitter_*.py`.

## What the code computes

Same Norris (2005) single-pulse model as the parent folder (`../../norris_fit.py`):

$$I(t) = A\,\exp\!\left[2\sqrt{\tau_1/\tau_2}\right]\exp\!\left[-\tau_1/(t-t_s) - (t-t_s)/\tau_2\right],\quad t>t_s$$

linear in amplitude $A$. Two ways of fitting a sum of $N$ such pulses to the same background-subtracted,
summed-NaI light curve $y(t)$ (all four bursts, using the BUG-23-fixed all-detectors-summed curve, raw,
no per-detector normalization — see the parent folder's writeup):

1. **Normalized (existing production method).** `NorrisFitter` (wraps `pymultifit.BaseFitter`, which
   calls plain `scipy.optimize.curve_fit(f, x, y, p0, bounds=...)` with no scaling option exposed) fit
   against $y/\max(y)$, amplitude seeds in $[0,1]$.
2. **Unnormalized, explicit `x_scale`.** The same pulses fit directly against raw $y(t)$ via
   `scipy.optimize.least_squares`, called directly (bypassing `pymultifit`, whose `BaseFitter.fit()` —
   `.venv/.../pymultifit/fitters/backend/baseFitter.py:262-291` — hardcodes the `curve_fit` call with no
   `x_scale`/`method` passthrough, so this cannot be reached through `NorrisFitter`). `x_scale` is an
   explicit array, one 4-tuple per pulse: `(Y_MAX_CTS_PER_S, 5.0, 5.0, 1.0)` — amplitude scaled to the
   data's own peak, `t_s`/`tau1` given an O(5s) scale, `tau2` O(1s). Amplitude seed is the normalized
   seed times `Y_MAX_CTS_PER_S`; `t_s`/`tau1`/`tau2` seeds unchanged.

Both methods use identical bounds ($A\ge0$, $t_s\in[t_\min,t_\max]$, $\tau_1,\tau_2\ge10^{-4}$) and,
where a burst needed it, an identical raised iteration budget (`max_nfev=20000`).

## Every judgement call, and why

- **`x_scale='jac'` (scipy's automatic per-parameter scaling) was tried first and rejected.** Tested on
  GRB231129C in the parent session before this folder existed: it converges, but to a *worse* local
  optimum (SSE $\approx1.03\times10^8$) — no better than no scaling at all (`x_scale=1.0`, SSE
  $\approx1.02\times10^8$), both far worse than the normalized method's $5.71\times10^7$. Automatic
  Jacobian-based scaling apparently isn't enough for this problem; only an explicit, domain-informed
  `x_scale` recovers the normalized method's result. This is why the script here goes straight to the
  explicit array and never calls `'jac'`.
- **`x_scale` values (5.0, 5.0, 1.0) are a deliberate order-of-magnitude choice, not a tuned fit.** They
  only need to be roughly right — `x_scale` sets step-size conditioning for the optimizer's trust region,
  not the objective function itself, so precision here doesn't matter the way a seed value's precision
  does. Chosen once and reused across all four bursts without per-burst tuning, and it worked for all
  four — the point of the experiment is that this generalizes without hand-tuning per burst.
- **GRB131014A's `max_nfev=20000` (raised from `NorrisFitter`'s default 5000)** mirrors a fix applied
  independently to the *production* `fitter_CLAUDE_GRB131014215.py` in the parent conversation while this
  experiment was running (that burst's summed 3-detector curve, peak $\approx1.3\times10^5$ cts/s, needs
  more iterations than the single-detector curve it used to be fit against, peak $\approx4.4\times10^4$
  cts/s — an iteration-budget issue, not a bad seed). Applied here to both methods for a fair comparison,
  so a convergence failure isn't mistaken for a genuine method difference.
- **GRB140206B uses `COMPLEX_P0`** (the 7-pulse decomposition), matching the paper's preferred model in
  `fitter_GRB140206275.py`, not the alternative `SIMPLE_P0` — reusing an existing decision, not remaking it.

## Results

**Overall fit quality is essentially identical between methods for every burst** (SSE computed in the
same raw counts/s² units both ways — the normalized method's model is rescaled back to physical units
before comparing):

| GRB | SSE, normalized (rescaled) | SSE, unnormalized (x_scale) | ratio |
|---|---|---|---|
| GRB080916C | 7.2019e7 | 7.2012e7 | 0.9999 |
| GRB131014A | 2.8928e9 | 2.8903e9 | 0.9991 |
| GRB140206B | 2.2587e8 | 2.2556e8 | 0.9986 |
| GRB231129C | 5.7087e7 | 5.7086e7 | 1.0000 |

(Per-pulse SSE values are identical within a row for a given burst — SSE is a whole-fit quantity — see
`fit_comparison_results.csv`'s `sse_normalized_rescaled`/`sse_unnorm_xscale` columns for exact figures.)
Every ratio sits within 0.15% of 1 — **normalization changes essentially nothing about fit quality**, for
all four bursts, confirming this generalizes beyond the one burst (GRB231129C) checked in the parent
session.

**Most individual pulse parameters agree closely between methods** — `t_peak` typically matches to
better than 0.01 s, `t_v` to a few percent, across the bulk of the 23 fitted pulses (see the CSV).

**GRB080916C updated 2026-09-23: pulse 5 dropped, matching a production fix.** This experiment originally
ran GRB080916C at 7 pulses and found its pulse 5 (`t_s≈20`, inside TR3) diverged sharply between methods
($t_v=0.046$ vs. $0.153$ s, a ~3.3x disagreement) despite near-identical overall SSE. That divergence
was the signal, not a coincidence: in the parent `fitter.py`, the same pulse independently converged to a
degenerate near-delta-function fit (`tau2` collapsing to 0.009, `mc_kept_fraction` dropping from 0.52 to
0.05) on this burst's summed-detector curve, confirmed against the raw data to be overfitting two noisy
bins rather than resolving a real feature. The user decided to drop it from the production 7-pulse model,
and this experiment was updated to match (now 6 pulses for GRB080916C, `p0` in `fit_comparison.py`). With
it gone, GRB080916C's SSE ratio improved from 1.0031 to 0.9999, and every remaining pulse agrees closely
between methods (`t_peak` within 0.02 s, `t_v` within ~5%) — no more weakly-constrained outliers for this
burst.

**Two GRB140206B pulses still disagree substantially between methods despite near-identical overall SSE,
and this is reported honestly rather than smoothed over:**

- **Pulse 6**: $t_s=22.6$ vs. $14.6$ s (nearly 8 s apart), $\tau_1=164$ vs. $550$ — though $t_v$ still
  agrees closely (3.95 vs. 4.06 s) despite the very different $t_s$/$\tau_1$. Classic $t_s\leftrightarrow
  \tau_1$ degeneracy: two very different-looking pulses produce almost the same peak time and width.
- **Pulse 7** (the broadest, lowest-amplitude "pedestal" pulse, $\tau_1\sim2000$): the unnormalized fit's
  $t_s$ sits at $-0.96$, essentially pinned against the fit window's own lower bound ($t_\min\approx-1$) —
  a sign this particular pulse is only weakly constrained by the data at all.

**Interpretation:** these are low-amplitude, broad, or otherwise weakly-constrained pulses whose exact
parameters barely affect the total model's SSE — the loss surface is genuinely flat in their direction, so
small differences in optimizer path (here, normalized vs. unnormalized) land in different but
similarly-good local optima. This is a property of *those specific pulses* being weakly constrained, not
evidence that one fitting method is more correct than the other, and not a normalization artifact — it
shows up in a pulse-by-pulse comparison precisely because both methods otherwise agree so closely
everywhere else. GRB080916C's now-resolved pulse 5 is the cautionary case: here, the divergence *did* turn
out to flag a real problem (a degenerate fit in production), so this kind of disagreement is worth treating
as a diagnostic signal, not dismissed by default.

## GRB140206B follow-up: window scan + pulse-6 reseed, and the recommended configuration

Both GRB140206B pulses flagged above were investigated further, in two more standalone scripts in this
same folder. **Neither production `fitter_GRB140206275.py`/`_simple.py` has been changed** — both still
use their original window `(-1, 160)` and converge cleanly there; this section documents a diagnosed
*alternative* configuration, not a correction to a broken one.

**Step 1 — `grb140206b_window_scan.py`: does widening the fit window to the light curve's own
`t.min()`/`t.max()` resolve the pinning?** Three windows tested, same `P0` otherwise: `narrow (-1, 160)`
(production), `wide (-20, 300)` (this burst's existing robustness-check window from
`../window_sensitivity_GRB140206275/window_sensitivity.py`), and `full (-133.184, 481.152)` (the data's
own bounds, read directly, not hardcoded).

| window | pulse 7 $t_s$ (norm / unnorm) | pulse 6 $\tau_2$ (norm / unnorm) |
|---|---|---|
| narrow | $-0.04$ / $-0.96$ (pinned at $t_\min$) | $1.440$ / $1.004$ (fine) |
| wide | $2.12$ / $-19.97$ (pinned at $t_\min$) | $6.108$ / $6.165$ (fine) |
| full | $-10.84$ / $-10.35$ (**resolved** — both methods agree to 0.5s) | $1.497$ / $\mathbf{0.00045}$ (**collapsed**) |

Widening to the full range lets pulse 7 finally converge on both methods independently agreeing on
$t_s\approx-10$ to $-11$s (a real precursor onset ~10s before trigger, invisible to either narrower
window) — but it destabilizes pulse 6 instead, which collapses to the exact same near-delta-function
pathology (`tau2`→0) that got GRB080916C's pulse 5 dropped from production. Net: no single window is
uniformly correct; each resolves one weak pulse while unsettling another. Full results:
`grb140206b_window_scan_results.csv`, `grb140206b_window_scan.png/.pdf` (all three windows' normalized-fit
total models overlaid on the full light curve — visually indistinguishable across the actual burst despite
the underlying per-pulse differences).

**Step 2 — `grb140206b_fullrange_fitter.py`: reseed pulse 6 instead of accepting the collapse.** At the
full-range window, pulse 6's seed was changed from `COMPLEX_P0`'s original $(A, t_s, \tau_1, \tau_2) =
(0.23, 23, 83, 1.1)$ — every other pulse's seed unchanged. Two amplitude values were tried for the
reseeded $\tau_1,\tau_2\to(1,1)$ neutral seed:

- **$A=0.23$ (unchanged):** pulse 6 does converge (no collapse), but drifts entirely out of its intended
  region — $t_s$ lands at $5.5$s (normalized) / $17.5$s (unnormalized), inside the *early* pulse complex
  rather than the $t\approx23$–$28$s region it's meant to model. Overall SSE drops (~5% lower than the
  original seed), but two *other* pulses (1 and 5) pick up new $\tau_1\to0.0001$ degeneracies as a side
  effect of the joint fit redistributing flexibility. Rejected: numerically better, physically worse.
- **$A=0.2$ (recommended):** pulse 6 stays anchored at $t_s\approx23$s in *both* methods and converges
  cleanly, with no knock-on degeneracy in any other pulse. This is the configuration below.

### Recommended configuration for GRB140206B (diagnosed, not yet applied to production)

Full range $(-133.184, 481.152)$s, `COMPLEX_P0` with pulse 6 re-seeded to $(0.2, 23, 1, 1)$:

| pulse | $A_\text{norm}$ | $t_s$ (norm) | $t_s$ (unnorm) | $\tau_1$ (norm) | $\tau_1$ (unnorm) | $\tau_2$ (norm) | $\tau_2$ (unnorm) | $t_\text{peak}$ (norm) | $t_\text{peak}$ (unnorm) | $t_v$ (norm) | $t_v$ (unnorm) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.144 | $-0.292$ | $-0.292$ | 0.087 | 0.089 | 1.799 | 1.774 | 0.105 | 0.106 | 0.940 | 0.931 |
| 2 | 0.367 | 4.515 | 4.491 | 4.495 | 4.550 | 16.479 | 16.485 | 13.121 | 13.152 | 11.442 | 11.472 |
| 3 | 0.567 | 10.890 | 10.882 | 6.513 | 6.593 | 1.488 | 1.478 | 14.004 | 14.004 | 1.865 | 1.861 |
| 4 | 0.241 | 26.559 | 26.478 | 3.467 | 3.674 | 3.212 | 3.158 | 29.895 | 29.885 | 2.944 | 2.942 |
| 5 | 0.023 | 24.191 | 24.191 | 0.0088 | 0.0072 | 4.053 | 1.737 | 24.380 | 24.303 | 1.583 | 0.705 |
| **6** | 0.051 | **23.114** | **23.184** | **1.000** | 0.580 | **0.447** | 0.528 | 23.782 | 23.737 | **0.481** | **0.486** |
| **7** | 0.070 | $-27.650$ | $-28.422$ | 3978.8 | 4048.2 | 5.675 | 5.637 | 122.622 | 122.645 | 24.393 | 24.375 |

SSE: $5.838\times10^8$ (normalized) vs. $5.839\times10^8$ (unnormalized) — within 0.02% of each other, and
close to the original-seed full-range fit's $5.843\times10^8$ (a comparably good optimum, not an
artificially lower one bought with a worse decomposition, unlike the rejected $A=0.23$ attempt).

**Both target pulses are fixed simultaneously, with no new side effects:**
- Pulse 6 stays at its intended $t_s\approx23$s, no longer collapses, and the two methods now agree on it
  to $\sim$1% ($t_v=0.481$ vs. $0.486$).
- Pulse 7 converges consistently ($t_s=-27.65$ vs. $-28.42$s, within 0.8s of each other, no longer pinned
  to any window edge; $t_v=24.39$ vs. $24.38$, agreeing to $\sim$0.1%).
- Pulse 1 is confirmed unaffected — its $\tau_1\approx0.087$–$0.089$ matches the original-seed full-range
  fit almost exactly, so the earlier $A=0.23$ attempt's pulse-1 degeneracy really was a side effect of that
  specific seed choice, not something inherent to widening the window.
- **Pulse 5 remains the one open loose end** — it was already inconsistent between methods in the
  *original*-seed full-range fit ($t_v=0.049$ vs. $1.000$, a ~20$\times$ disagreement, predating this
  reseed and unrelated to pulse 6). This configuration narrows that gap ($t_v=1.583$ vs. $0.705$, ~2$\times$)
  but does not fully resolve it — flagged for follow-up, not silently left out.

Full results: `grb140206b_fullrange_pulse6_reseed_results.csv`,
`grb140206b_fullrange_pulse6_reseed.png/.pdf`.

**Why this stays a diagnosed alternative, not a production change (yet):** this required a hand-picked
amplitude seed found by trial (`0.23` failed, `0.2` worked) at a specific window, on top of the already
hand-tuned `COMPLEX_P0`/`SIMPLE_P0` this burst's production fit carries. That is a real, deliberate
judgement call requiring the same kind of scrutiny this project gives every pulse-decomposition choice
(per `../../variability_analysis.md`'s existing conventions) — recorded here as the best configuration
found, available to promote to production if/when that decision is made, not applied unasked.

## Conclusion

**Peak-normalization is not load-bearing for the physics — it's a solver convenience, and a correct
one.** An equally-valid unnormalized method (explicit `x_scale`) reproduces the same overall fit quality
and the same well-constrained pulse parameters, for all four bursts, without touching the data. The two
weakly-constrained-pulse disagreements above are a property of those pulses, not of the normalization
choice — the same ambiguity would likely surface as fit-to-fit instability (different random seeds,
slightly different `P0`, etc.) even within a single method. This confirms the parent session's original
conclusion (based on GRB231129C alone) generalizes: normalizing before fitting is the right call, not
because raw-scale fitting is impossible, but because it's the simplest implementation of a fix that's
otherwise available (explicit `x_scale`) and doesn't need a burst-by-burst tuned scale array.

## Files in this folder

- `fit_comparison.py` — the main script, all four bursts at their production window; run with
  `PYTHONPATH=src:codes-for-paper/variability_analysis` from `GRBResearchWork/` (add
  `:codes-for-paper/variability_analysis/experiments/normalized_vs_unnormalized_fit` to the `PYTHONPATH`
  for the two follow-up scripts below, since they import from this one).
- `fit_comparison_results.csv` — 23 rows (one per fitted pulse across all 4 bursts; GRB080916C is 6 pulses
  as of the 2026-09-23 update above), both methods' full
  parameter sets, `t_peak`, `t_v`, and per-burst SSE side by side.
- `fit_comparison_<GRB_dir>.png/.pdf` — one figure per burst: the summed light curve with both fitted
  total models overlaid (they are visually indistinguishable at this scale for all four bursts).
- `grb140206b_window_scan.py` / `..._results.csv` / `..._scan.png/.pdf` — GRB140206B only, three fit
  windows (narrow/wide/full) compared, both methods. See "GRB140206B follow-up" above.
- `grb140206b_fullrange_fitter.py` / `..._pulse6_reseed_results.csv` / `..._pulse6_reseed.png/.pdf` —
  GRB140206B at the full-range window with pulse 6 re-seeded; produces the recommended configuration
  documented above.
- `why_not_unnormalized_fitting.md` — standalone reviewer-facing justification for why this project's
  production pipeline fits the peak-normalized light curve rather than the raw one, referencing the full
  investigation in this file. Written to stand alone if a referee asks the question directly.
