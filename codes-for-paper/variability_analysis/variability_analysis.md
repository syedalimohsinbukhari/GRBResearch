# Variability timescale (t_v) via manually-seeded joint Norris fit — GRB131014A

Companion to `PHASE5_TV_PLAN.md` (project root), which documented the earlier automated
(MEPSA/scipy peak-detection + per-window local fit) pipeline for this same Phase 5 goal. That pipeline's
source files (`variability_timescale/norris_fit.py`, `variability_timescale.py`, etc.) have since been
deleted — see `PHASE5_TV_PLAN.md`'s 2026-09-16 status update — so this folder is now the only live track.
This folder is a **separate, manually-driven track**, worked one GRB at a time starting with GRB131014A on
2026-09-07 after the automated pipeline stalled on GRB080916C's TR2 (zero MEPSA detections even with
padding, because TR2's true shape is one broad, smoothly-declining pulse with no local excess for a spike
detector to find — see the MEPSA-fallback discussion in this session's history). GRB140206B
(`fitter_GRB140206275.py` + a `_simple` variant) and GRB231129C (`fitter_GRB231129779.py`) were fit in
this same manual style in a later session (2026-09-15/16), and GRB080916C (`fitter.py`, 8-pulse) followed
in a session on 2026-09-16/17 — **none of the three are yet written up below** — the "What the code
computes" / "Every judgement call" / "Results" sections that follow describe the GRB131014A fit only; see
each script's own docstring/comments and the sections further down for the other three bursts until this
note is extended. Not yet wired into `lorentz_factor.py`'s `Gamma_min` pipeline — same deliverable-boundary
stance as the abandoned automated pipeline.

**Open cross-burst issue, flagged 2026-09-18, not yet addressed:** every fit in this folder (all four
bursts) uses only a single NaI detector (`nai_data[0]`, whichever sorts first alphabetically — `n3` for
GRB080916C, `na` for GRB131014A, `n3` for GRB140206B, `n7` for GRB231129C) rather than summing across the
GRB's full NaI detector set, unlike the convention typical of published GRB analyses. This wasn't a
deliberate choice recorded anywhere — it's just what every `fitter*.py`/`window_sensitivity.py` script has
done since the very first one. Revisiting this (which detectors, how to combine background/rate/errors
across them — see BUG-22 below for the right way to combine per-channel errors) is planned as follow-up
work; every `t_v`/pulse-decomposition result in this file predates that rework and should be treated as
single-detector until it's redone.

## What the code computes

Same Norris (2005) pulse and $t_v$ definition as `variability_timescale/norris_fit.py` (Bukhari et al. 2022,
Adv. Space Res., eq. 9/10) — this folder's `fitter.py`/`fitter_CLAUDE_GRB131014215.py` import that module directly rather
than redefining it. $t_v = (\tau_2/2)\sqrt{(\ln2+2\sqrt{\tau_1/\tau_2})^2-4\tau_1/\tau_2}$, the FWHM/2 of the
fitted pulse.

The fit itself: a single joint `NorrisFitter` (pymultifit `BaseFitter`, 5 pulses, 20 free parameters) over
one fixed window, `[-1, 10]` s since trigger, of the 10–400 keV NaI background-subtracted light curve
(GRB131014215, detectors `na`/`nb`/`n9`), rather than the automated pipeline's per-episode local windows.

## Every judgement call, and who made it

- **Manual joint fit over automated per-window detection.** After the automated MEPSA-based pipeline
  (`variability_timescale/`) hit real, GRB080916C-specific problems (TR2 needing a direct-max fallback,
  TR4 landing on a degenerate solution), the user instead hand-fit GRB131014A directly with carefully
  chosen $p_0$ seeds spanning the whole EX0/TR1/TR2/EX1 region at once. **User decision**, motivated by
  wanting a clean, verifiable reference fit for this burst before deciding whether/how to generalize the
  automated pipeline further.
- **Amplitude/shape decoupling.** The fit runs on the light curve normalized to peak=1 (for numerical
  conditioning — feeding raw counts/s directly into $p_0$ made the optimizer behave badly), and $A$ is
  rescaled back to physical counts/s *after* fitting via the recorded `Y_MAX_CTS_PER_S`, exact since $A$ is
  purely multiplicative in the Norris formula. **Claude decision**, prompted directly by the user's own
  observation that pre-scaling $p_0$'s amplitude by $y_\text{max}$ "screws the fitter over" — confirmed the
  root cause is optimizer conditioning, not a flaw in the rescale itself.
- **Guess-sensitivity used as the actual reliability criterion, not `mc_kept_fraction` alone.** Refitting
  with every $\tau_1=\tau_2=1$ (a neutral, uninformative shape seed, only $A$/$t_s$ roughly hinted) reproduced
  pulses 3, 4, 5 to within 0.01 s in $t_\text{peak}$ and comparable `mc_kept_fraction` (0.77–1.00 either way),
  but landed pulses 1 and 2 on genuinely different $t_\text{peak}$/$t_v$ each time, with `mc_kept_fraction`
  staying poor (0.10–0.21) under both seedings. **Jointly established**: the user proposed the neutral-seed
  test as an argument that careful guessing wasn't load-bearing; the result confirmed it only for pulses
  3–5, and additionally showed pulses 1/2 are unreliable *because* they sit in the earliest, heavily
  overlapping part of the light curve (Hakkila 2026's "overlapping pulses are not unique" failure mode), not
  because of any fitter/seeding problem that better guessing would fix.
- **$(\tau,\xi)$ reparameterization tried and found not to help, for this fit specifically.** Converting the
  same $p_0$ via `norris_reparam.tau_xi()` and refitting with `NorrisReparamFitter` landed pulses 1–3 on
  different local minima entirely (pulse 2 converged to a *worse* degenerate solution, $t_v=252$ s in an
  11 s window), and where it did land near the same point (pulses 4, 5), the $t_s$–$\tau$/$t_s$–$\xi$
  correlations stayed at 0.94–0.99 — no better than the original $t_s$–$\tau_1$ correlation it was meant to
  fix. **Negative result, recorded so it is not re-tried blindly**: Recipe C's premise
  (`variability_timescale/norris_shenanigans.md`) doesn't hold for this joint 5-pulse fit.
- **Photon-to-pulse assignment: nearest preceding pulse onset.** Each LAT photon is assigned to the pulse
  with the largest $t_s \le$ the photon's arrival time — i.e. whichever pulse has "turned on" and not yet
  been superseded by the next pulse's onset. **Claude-proposed, verified against the user's own manual
  identification** of two specific photons (445.28 MeV at $t=2.533$ s → pulse 4; 552.79 MeV at $t=4.111$ s →
  pulse 5) before being generalized to all 10 photons in the window. Not the only possible rule (nearest by
  absolute time, or nearest to $t_\text{peak}$, would be alternatives) but this one is physically motivated
  and reproduced the verified cases exactly.
- **Episode boundaries from `results.json` directly**, not eyeballed or taken from the automated pipeline's
  photon-centered heuristic windows: `EX0 -0.192_2.432`, `TR1 0.960_2.432`, `TR2 2.432_4.160`,
  `EX1 2.432_6.976`. A pulse can match more than one episode (EX0 fully contains TR1's window, EX1 fully
  contains TR2's) — matching is "does this pulse's $t_\text{peak}$ fall inside this episode's boundary", not
  a 1:1 assignment.

## Conventions that could plausibly have gone another way

- The fixed `[-1, 10]` s window and 5-pulse count were chosen by inspection for this burst specifically —
  not derived from a general rule, and not expected to transfer unchanged to another GRB.
- $p_0$ itself (in `fitter.py`) is the user's own hand-tuned values, arrived at by trial and error against
  the normalized light curve; left untouched throughout this session's work per explicit instruction.
- The twin-axis photon overlay's y-limits are set to `(0, data.max())` on both axes (no autoscale padding)
  specifically so 0 and each series' own peak align visually — a display choice, not a data transform.

## Validation performed

- **Neutral-seed reproducibility test** (above): the actual mechanism used to separate "real, data-determined"
  pulses (3, 4, 5) from "seed-dependent, non-unique" ones (1, 2) — not just asserted, run twice with visibly
  different $\tau_1,\tau_2$ starting points and compared.
- **Photon table cross-checked against `LAT_analysis/014__GRB131014215/` and `lorentz_results.csv`**: the
  1232.92 MeV (TR1/EX0), 1021.24 MeV (TR2), 1193.12 MeV (EX1) "defining" photons already used elsewhere in
  the project appear at the same arrival times in this window's full photon list, confirming
  `GRB131014215_lat.fits` (copied here from `light_curves/GRB131014215/lat.fits`, not re-derived) and the
  hardcoded `T0_MET_S = 403420143.2` (read from
  `LAT_analysis/014__GRB131014215/Ep1__0.960_2.432/GRB131014A_fit_results_0.96_2.432.txt`'s own `T_0` line)
  are correct: the FITS file's own GTI (403420144.16–403420147.36) equals `T0_MET_S+0.96` to `T0_MET_S+4.16`,
  exactly the T90 window.
- $A_\text{cts/s} = A_\text{norm}\times Y_\text{MAX\_CTS\_PER\_S}$ is an exact linear identity, not an
  approximation — verified by re-plotting both data and every fitted curve at physical scale and confirming
  the total fit still tracks the data (a scale error would show as a visible mismatch; none appeared).

## Window-widening as a pulse-count diagnostic (found via GRB231129C, applies generally)

Found while investigating a *different* burst (GRB231129C, `experiments/window_sensitivity_GRB231129779/`)
but the conclusion applies to every joint Norris fit in this project, including the GRB131014A one
documented above — recorded here since this is the method note for that fit family.

**The check:** refit the same pulses, from the same $p_0$, over two different window widths (e.g.
$[-1,10]$s vs $[-10,20]$s), and compare not just each pulse's $t_v$ but its `mc_kept_fraction`. This is
possible because of a specific, verified mechanism: `NorrisFitter.fit_boundaries()`
(`codes-for-paper/variability_analysis/norris_fit.py:49-50`) sets $t_s$'s lower bound to `x_values.min()`
— the fit window's own left edge is $t_s$'s box constraint — so widening the window directly widens how
far the optimizer can push $t_s$ along the already-established near-total $t_s$↔$\tau_1$ degeneracy (this
session's correlation matrices, TR4, the $(\tau,\xi)$ reparam test above).

**What it found for GRB231129C, concretely:** a 4-pulse fit's pulse 3 kept its $t_v$ point estimate stable
under the window change (0.505→0.490, ~3%) but its `mc_kept_fraction` collapsed (0.451→0.148) — the median
survived, the reliability didn't, and checking $t_v$ alone would have missed this entirely. Adding a 5th
pulse (there was real structure near $t\approx4.6$–$5.5$s the 4-pulse model had nowhere to put) made the
instability disappear outright: every pulse's $t_v$ shift dropped to <10%, and `mc_kept_fraction` *improved*
for all five pulses going from narrow to wide, instead of collapsing for one.

**The methodological conclusion:** window-sensitivity is a **necessary, not sufficient**, diagnostic for an
under-specified pulse count. If parameters (especially `mc_kept_fraction`, not just $t_v$) are stable
under a window change, that's real evidence the count is right. If they're not, the confirming step is
adding a pulse and checking whether the instability actually goes away — the instability alone is a red
flag, not proof; a genuinely isolated pulse could in principle be unstable for its own unrelated reasons.

**GRB131014A's pulses 1 and 2, checked against this diagnostic (2026-09-16) — no collapse found, and a
third window strengthens the result.** `experiments/window_sensitivity_GRB131014215/window_sensitivity.py`
ran the same comparison on the full 5-pulse fit across **three** window widths — narrow `[-1,10]`s (the
fitter's own default), wide `[-10,20]`s, and widest `[-135.9,478.4]`s (the light curve's own
`np.min(t)`/`np.max(t)`, the most extreme box `fit_boundaries()` can ever present) — using each pulse's
own narrow-window converged parameters as $p_0$ for all three (so only the window bound changes, isolating
it from a seed-choice effect). Result:

| pulse | $t_s$: narrow → wide → widest | `mc_kept_fraction`: narrow → wide → widest |
|---|---|---|
| 1 | −0.680 → −0.683 → −0.682 | 0.206 → 0.342 → **0.860** |
| 2 | −0.028 → −0.028 → −0.014 | 0.127 → 0.209 → **0.650** |
| 3 | 1.185 → 1.185 → 1.184 | 0.998 → 1.000 → 1.000 |
| 4 | 2.471 → 2.471 → 2.471 | 0.999 → 1.000 → 1.000 |
| 5 | 2.659 → 2.658 → 2.658 | 0.774 → 0.908 → 1.000 |

Every pulse's $t_s$ stays essentially anchored even at the widest window (≤0.014 s drift) — nothing here
rides the box-constraint edge the way GRB140206B's SIMPLE pulse 5 does (below) — and `mc_kept_fraction`
climbs *monotonically* for every pulse, including 1 and 2, whose values roughly triple to quadruple by the
widest window. This is the same signature GRB231129C showed only *after* its missing 5th pulse was added,
not the collapse signature that flagged that gap in the first place. By this diagnostic's own logic, that's
real evidence the 5-pulse count is right here — **no hint that a 6th pulse is needed** — and the monotonic
*improvement* specifically argues against a missing-pulse explanation: if pulses 1/2's poor constraint came
from an absent 6th pulse, widening the window should have made things worse or unstable, not steadily
better. The improvement itself is attributable to the extra off-burst baseline in the wider windows
tightening the fit's overall covariance, not to any change in what the pulses themselves resolve to.

This resolves the open question above: pulses 1 and 2's poor *absolute* `mc_kept_fraction` at the narrow
window does not collapse under widening — it improves steadily and their $t_s$ stays fixed — so it reads as
a genuine, irreducible degeneracy in that early overlapping structure (consistent with the Hakkila 2026
"overlapping pulses are not unique" framing already invoked above) rather than an under-fitting artifact.
The "unresolved" call in Known Limitations below now rests on two independent checks (seed-sensitivity and
window-widening) agreeing, not on seed-sensitivity alone. Diagnostic plots (fitted light curve + individual
pulses, one per window) are saved alongside the CSV in that same experiment folder.

**GRB231129C's own fixed 5-pulse model, re-checked at the same third window (2026-09-16).** The two-window
result above was GRB231129C's *post-fix* state; extending it to the widest `[-138.5,475.8]`s window
confirms the same monotonic pattern holds all the way out: every pulse's `mc_kept_fraction` keeps climbing
(e.g. pulse 3 — the one whose collapse originally flagged the missing 5th pulse — goes 0.372→0.479→0.637;
pulse 5 goes 0.896→0.942→0.991), $t_s$ shifts stay small (≤0.24 s, largest for pulse 3), and $t_v$ drifts up
to ~12% by the widest window for pulses 4/5 but with no reversal or instability. No evidence the fix was
incomplete.

**GRB140206B, both candidate models, checked at all three windows (2026-09-16) — the one case where a
pulse's $t_s$ genuinely rides the box edge.** `experiments/window_sensitivity_GRB140206275/window_sensitivity.py`
ran narrow `[-1,160]`s (the fitter's own default, covering all 7 episodes), wide `[-20,300]`s, and widest
`[-133.2,481.2]`s (this burst's own `np.min(t)`/`np.max(t)`) on both SIMPLE (5-pulse) and COMPLEX (7-pulse)
decompositions, using each model's own narrow-window converged parameters as $p_0$ for all three windows.

Both models show the same "no collapse" pattern as the other two bursts for their compact, well-localized
pulses — $t_v$ shifts stay under 2% and `mc_kept_fraction` holds flat or improves. But **SIMPLE's pulse
5** (the broad ~3440 s-$\tau_1$ pedestal underlying the whole burst) is different: its $t_s$ does not
anchor, it drifts continuously along the $t_s$↔$\tau_1$ degenerate ridge as the window opens:

| window | left edge | $t_s$ | $\tau_1$ | $t_v$ | kept |
|---|---|---|---|---|---|
| narrow | −1.0 | −0.960 | 3441 | 17.22 | 0.613 |
| wide | −20.0 | −19.968 | 5509 | 17.09 | 0.647 |
| widest | −133.2 | −60.654 | 11585 | 16.88 | 0.569 |

At the wide window this looked pinned to the box edge (−19.97 vs. a −20 boundary), but the widest window
shows that was coincidental, not a true pin — given more room, $t_s$ moved to −60.65, well short of the
−133.2 boundary, with $\tau_1$ roughly tripling again. $t_s$ and $\tau_1$ are trading off against each
other with the data unable to pin either one; $t_v$ (which depends on both through the eq. 10 combination)
stays comparatively steady, drifting only 2% across the whole three-window span, and `mc_kept_fraction`
wobbles (0.61→0.65→0.57) without a clean collapsing trend, reading as MC noise around a genuinely
poorly-constrained pulse rather than a widening-driven failure.

**COMPLEX's pulse 5 is a different physical object** — one of the narrow 23–28 s-cluster pulses
($t_s\approx23.08$, $\tau_1\approx1.1$), not the pedestal — and it is essentially rock-solid across all
three windows ($t_s$ moves 0.008 s total; kept even improves slightly, 0.607→0.630→0.642). **COMPLEX's
actual pedestal is pulse 7**, and it behaves very differently from SIMPLE's pulse 5 despite similar
$\tau_1$: $t_s$ only creeps from −0.95 to −1.65 (not −60), and $\tau_1$ barely moves (3306→3523 vs.
SIMPLE's 3441→11585). The extra narrow pulses COMPLEX carves out of the 23–28 s region (pulses 4/5/6)
appear to anchor the pedestal in place — splitting that region leaves the optimizer less room to let the
broad component absorb ambiguous flux, so pulse 7's own $t_s$/$\tau_1$ stay far better constrained than
SIMPLE's equivalent pulse. This is a concrete, quantified reason to prefer COMPLEX's decomposition if the
pedestal's $t_s$/$t_v$ ever matters downstream, beyond the photon-anchoring argument already in
`fitter_GRB140206275.py`'s own docstring — and it directly answers that docstring's own open question:
COMPLEX pulses 4/5/6's already-mediocre `mc_kept_fraction` (0.47–0.65) does not get worse under widening,
so the photon-anchored 3-way split is a legitimate resolution-vs-stability tradeoff, not a symptom of a
wrong pulse count or a window artifact.

**Cross-burst comparison.** None of the other 8 pulses checked across GRB131014A and GRB231129C show
anything like SIMPLE pulse 5's runaway drift — everywhere else $t_s$ stays anchored to within a fraction of
a second regardless of window width. GRB140206B's SIMPLE pedestal is a genuine, burst-and-model-specific
finding, not a general property of this fit family.

**GRB080916C, 6-vs-7-pulse question (2026-09-16) — inconclusive via `mc_kept_fraction`, resolved via
residuals instead.** `experiments/window_sensitivity_GRB080916009/window_sensitivity.py` checked whether a
6-pulse model (`fitter.py`'s original decomposition) or a 7-pulse model (an extra pulse at $t_s\approx20$s,
inside TR3's window) is better supported, across narrow `[-1,70]`s, wide `[-20,150]`s, and widest
`[-25.9,300.6]`s. Unlike GRB231129C's original case, **neither model shows a `mc_kept_fraction` collapse**
at any window — the candidate 7th pulse's own reliability is mediocre but stable (0.517→0.604→0.604), so
this diagnostic alone doesn't decide it here. What did: comparing fit residuals directly. The 7-pulse
model's total SSE only improves ~3.2% globally (the amount you'd expect just from 4 more free parameters
soaking up noise everywhere) but ~12% specifically in the `[15,30]`s region the new pulse targets —
disproportionate local improvement is the signature of a real feature, not an overfitting artifact. Plots:
`window_sensitivity_{six,seven}_{narrow_-1_70,wide_-20_150,widest_full_range}.png/.pdf` in that folder.

This same residual-inspection approach, applied next to a *different* region of GRB080916C (the
$t\approx0$–$3$s span, where a single broad pulse was overshadowing two more real peaks), is what led to
the 8-pulse decomposition now in `fitter.py` — see the dedicated section below.

## GRB080916C: 8-pulse decomposition via residual-seeding, and a CERN ROOT cross-check (2026-09-16/18)

Not yet written up as its own "What the code computes" / "Every judgement call" / "Results" set (see the
note at the top of this file) — this section covers two specific, self-contained pieces of that burst's
work: how the final 8-pulse `p0` in `fitter.py` was actually found, and a follow-up cross-check of the
fit in CERN ROOT that surfaced a real, now-fixed bug in `light_curves.py`.

**The problem.** Starting from a 6-pulse fit (`norris1`–`norris6`), `norris2` (a broad envelope/tail,
$\tau_2\approx6.7$) visibly overshadows the pulses after it — the user's own read: "norris2 is overshadowing
[norris3] and is even bleeding into norris3/4/5", with two more real peaks suspected near $t\approx1.5$s and
$t\approx2.4$s that repeated attempts to fit directly (jointly, from hand-typed `p0`) could not recover:
the optimizer kept collapsing the two new pulses to degenerate near-zero-width spikes
($\tau_2\sim10^{-3}$, `kept_fraction=0`), because `norris2`'s tail claims the flux there first.

**The fix: seed from the residual, not the raw light curve.** Subtracting the converged 6-pulse fit from
the data (`fitter_EXPERIMENT.py`) shows two coherent, smoothed residual bumps — `+0.065`..`+0.079` near
$t=1.41$–$1.54$s and `+0.072`..`+0.086` near $t=2.43$–$2.62$s, separated by a `-0.08`..`-0.10` dip — real
structure, not noise (confirmed by 5-bin-moving-average smoothing, same technique used to originally spot
GRB131014A's/GRB231129C's residual bumps). Seeding the two new pulses (`norris2-1`, `norris2-2`) at the
*residual's own amplitude* (~0.08) instead of the raw light curve's amplitude (~0.9), while seeding every
other pulse from the already-converged 6-pulse fit rather than fresh hand-guesses, converges cleanly —
both land close to their target peaks with sane, non-degenerate parameters.

**Neutral-seed reproducibility, the same criterion as GRB131014A's pulses 1/2 (above).** Refitting
`norris2-1`/`norris2-2` from a different, uninformative $\tau_1,\tau_2=0.3,0.3$ seed instead of the
residual-informed one: `norris2-1` reproduces to $t_\text{peak}=1.480$s, $t_v=0.049$s (vs. 1.466s/0.051s) —
**reads as real**, same tight agreement standard as GRB131014A's reliable pulses. `norris2-2` reproduces to
$t_\text{peak}=2.462$s, $t_v=0.177$s (vs. 2.309s/0.148s) — a looser match, **probably real but not yet as
tight**; one more seed trial before trusting its exact numbers. `fitter.py`'s `P0` comment documents both
results and the full recipe directly at the parameter list, so this reasoning isn't lost again (an earlier
hand-edit of that same `p0` had already lost the original 6-pulse values once).

**CERN ROOT cross-check.** To let the fit be reproduced/verified outside this project's Python stack, the
GRB080916C 10–400 keV light curve (NaI `n3+n4`, background-subtracted, summed across detectors — see the
"Open cross-burst issue" note above on why this is still single-burst-only, not yet the project-wide
detector-summing convention) was exported to `GRB080916009_lightcurve_10-400keV.csv`, and
`GRB080916009_norris_fit.C` reproduces the same 8-pulse `TF1` fit in ROOT/Minuit (same `norris_pulse` form,
same `fit_boundaries()`-equivalent parameter limits, seeded from the Python fit's converged values rescaled
to the CSV's own physical count-rate scale). Both live in
`experiments/root_fit_GRB080916009/`, moved there 2026-09-18 (previously sat directly in this folder).

Running the errors-weighted fit (`TGraphErrors`, `use_errors=true`) first came back with
**`chi2/ndf = 33.4/1077 ≈ 0.031`** — far below the ~1 a correctly-weighted fit should give, meaning the
error bars were too large relative to the data's real scatter. An unweighted fit (`TGraph`, ignoring
`net_err_cts_per_s` entirely) tracked the data just as well, with an implied RMS residual of only
`sqrt(chi2/n) ≈ 252` cts/s — much smaller than what the (inflated) `net_err_cts_per_s` column claimed
(~830–1000 cts/s). That discrepancy traced to a real bug in `light_curves.py::lightcurve_data()`:
combining ~89 independent per-channel errors (10–400 keV spans that many energy channels) with a plain
linear `+=` instead of in quadrature, inflating the combined error by close to $\sqrt{89}\approx9.4\times$
— see `BUGS.md` BUG-22 for the full writeup, confirmation, and fix. **Checked before fixing**: every
`lightcurve_data()` call site in the repo (`make_lightcurve.py`, every `fitter*.py`, every
`window_sensitivity.py`) uses `errors=False`, the default — `errors=True` was never actually invoked
anywhere except this ROOT-export exercise, so the bug had zero effect on anything published. After the
fix, the same weighted ROOT fit came back at **`chi2/ndf = 1565.5/1077 ≈ 1.454`** — a physically credible
value, consistent with the unweighted fit's own RMS-residual estimate. `norris2-2` still shows the same
Minuit-flagged degeneracy (`τ₁` pinned at its lower limit, `ERR MATRIX NOT POS-DEF`) in every version of
this fit, weighted or not — expected, not a new problem from any of this.

## Results (this session, 2026-09-07)

| pulse | episode(s) | $t_\text{peak}$ [s] | $t_v$ [s] | kept | anchoring photon |
|---|---|---|---|---|---|
| 3 | TR1 (+EX0) | 1.766–1.774 | $0.375$–$0.381$ | 0.91–0.999 | 1232.92 MeV @ 1.984 s (+6 more, all in decay) |
| 4 | TR2 (+EX1) | 2.697–2.699 | $0.552$–$0.561$ | 0.999–1.00 | 445.28 MeV @ 2.533 s |
| 5 | TR2 (+EX1) | 3.589–3.590 | $0.266$–$0.269$ | 0.77–0.85 | **1021.24 MeV @ 3.184 s** (TR2's existing $\Gamma_\text{min}$ photon) + 552.79 MeV @ 4.111 s |

## Known limitations and open questions

- **Pulses 1 and 2 are unresolved**, not merely unmeasured: two different, reasonable seedings give two
  different answers, and at the narrow `[-1,10]`s window both stay poorly MC-constrained
  (`mc_kept_fraction` 0.10–0.21 under seed-sensitivity, 0.13–0.21 in the window-widening baseline). No LAT
  photon falls in either pulse's active window (first photon arrives at 1.852 s), so this doesn't block
  feeding TR1/TR2's $\Gamma_\text{min}$ — but it does mean the earliest part of this burst's variability
  structure is not characterized here. **Checked against the window-widening pulse-count diagnostic** (see
  that section above), including a third, most-extreme window pushed to the light curve's own
  `np.min(t)`/`np.max(t)` (2026-09-16): neither pulse's `mc_kept_fraction` collapses under widening — both
  climb monotonically instead, reaching 0.860 (pulse 1) and 0.650 (pulse 2) at the widest window, while
  $t_s$ stays fixed to within 0.014 s — so the diagnostic finds no evidence a missing 6th pulse is the
  cause. The "unresolved" call now rests on two independent checks (seed-sensitivity and window-widening)
  agreeing, not on seed-sensitivity alone.
- **TR2 now has two candidate anchor pulses**, not one: pulse 4 (445 MeV photon, extremely well-constrained,
  $t_v\approx0.56$ s) and pulse 5 (the *already-published* 1021 MeV defining photon, less tightly
  constrained, $t_v\approx0.27$ s). Which one should replace TR2's current duration-based $t_v$ in
  `lorentz_factor.py`, if either, is an open decision — not resolved here, and not yet fed back into that
  pipeline (same as the rest of Phase 5).
- **This is one GRB, hand-tuned.** Whether the neutral-seed reliability test, the photon-assignment rule, or
  the 5-pulse/[-1,10]s window choice generalize to the other three bursts is untested; the user's plan is to
  work through them manually one at a time.

## Files here

GRB131014A (documented above):
- `fitter.py` — the user's own working file, reused/repurposed across bursts as this track progressed
  (currently set up for GRB140206B, not GRB131014A — see below); not touched by Claude past the two
  changes explicitly requested for the GRB131014A stage (episode-boundary `axvline`s, fixing a broken
  `savefig` path).
- `fitter_CLAUDE_GRB131014215.py` — Claude's copy, used for the amplitude-rescale, physical-unit plot, LAT-photon overlay,
  and photon-to-pulse assignment work documented above. Diverges from `fitter.py` from that point on.
- `GRB131014215_lat.fits` — copied from `light_curves/GRB131014215/lat.fits`; the FITS-reading convention
  (`fits.open(...)[1].data`, per-source probability column named after the source) is copied from
  `light_curves/make_lightcurve.py`, not re-derived.
- `norris_fit_results_GRB131014215.csv` — one row per (pulse, episode) match, produced by `fitter_CLAUDE_GRB131014215.py`.
- `norris_fitted_GRB131014215.png/.pdf` — the decorated plot (physical units, episode boundaries, LAT
  photon overlay on twin axis).
- `experiments/window_sensitivity_GRB131014215/` — the window-widening pulse-count diagnostic applied to
  this burst (see that section above), added 2026-09-16: `window_sensitivity.py` (+ its own local
  `light_curves.py`/`norris_fit.py` copies), `window_sensitivity_results.csv` (3 windows × 5 pulses), and
  one fitted-light-curve plot per window width (`window_sensitivity_narrow_-1_10.png/.pdf`,
  `window_sensitivity_wide_-10_20.png/.pdf`, `window_sensitivity_widest_full_range.png/.pdf`).

GRB140206B and GRB231129C fits themselves are still **not yet written up** as their own "What the code
computes" / "Every judgement call" / "Results" sections — see the note at the top of this file. The file
inventory below is complete, though, so it's clear what each artifact is and how it was produced.

GRB140206B — two live candidate decompositions, kept side by side rather than one being picked as final
(see `fitter_GRB140206275.py`'s own docstring, and the window-widening section above for the diagnostic
that independently supports COMPLEX for the pedestal pulse specifically):
- `fitter_GRB140206275.py` — the **COMPLEX** model: 7 pulses, splitting the $t\approx23$–$28$s region into
  three (pulses 4/5/6) specifically to give the 753.11 MeV photon at $t=23.998$s (TR3's own $\Gamma_\text{min}$-defining
  photon, already in `lorentz_results.csv`) its own dedicated peak (pulse 5) instead of leaving it buried in
  a single broad pulse's tail. Fits over `[-1, 160]`s (all seven episodes, through TR6's 154.240s end).
  Uses **dominant-flux photon assignment** (`assign_pulse` picks whichever pulse has the largest
  model-predicted flux at the photon's arrival time), not GRB131014A's nearest-preceding-onset rule — see
  `fitter_GRB231129779.py` for why that rule was replaced. Also fits the **SIMPLE** model inline (for the
  diagnostic comparison below) but only writes COMPLEX's results as this file's main results CSV/plot.
- `fitter_GRB140206275_simple.py` — the **SIMPLE** model: 5 pulses, treating that same $t\approx23$–$28$s
  region as one broad pulse (pulse 4) instead of splitting it. More numerically stable there (`kept≈1.0` vs
  COMPLEX's 0.54–0.67 across pulses 4/5/6) but doesn't resolve TR3's photon into its own peak — the
  stability-vs-photon-resolution tradeoff is why both files are kept rather than one being deleted. Same
  dominant-flux photon assignment as the COMPLEX file.
- `GRB140206275_lat.fits` — copied from
  `light_curves/GRB140206275/GRB140206Bfiltered_gti_gtsrcprob_7.488_154.176.fits`; `T0_MET_S = 413361375.84`
  read from `LAT_analysis/007__GRB140206275/Ep1__7.488_11.072/*_fit_results_*.txt`'s own `T_0` line.
- `norris_fit_results_GRB140206275.csv` — COMPLEX model's one-row-per-(pulse, episode) results, produced by
  `fitter_GRB140206275.py`.
- `norris_fit_results_GRB140206275_simple.csv` — SIMPLE model's equivalent, produced by
  `fitter_GRB140206275_simple.py`.
- `norris_fit_diagnostics_GRB140206275.csv` — SIMPLE-vs-COMPLEX comparison, one row per SIMPLE pulse matched
  to its nearest-$t_\text{peak}$ COMPLEX counterpart (`t_peak_offset_s` column), produced by
  `fitter_GRB140206275.py`. This is where the stability tradeoff above is quantified: SIMPLE pulse 4
  (`kept=0.9998`) matches COMPLEX pulse 4 (`kept=0.6727`, `t_peak_offset≈0.34`s), and SIMPLE pulse 5 (the
  broad pedestal, `kept=0.6093`) matches COMPLEX pulse 7 (`kept=0.623`, `t_peak_offset≈0.19`s) — the same
  SIMPLE-pulse-5 / COMPLEX-pulse-7 pedestal pairing independently identified by the window-widening
  diagnostic above.
- `norris_fit_GRB140206275.png` — an earlier, intermediate COMPLEX-model plot (2026-09-07, ~05:42, roughly
  18 minutes before the final `norris_fit_results_GRB140206275.csv`/`norris_fitted_GRB140206275.png` at
  ~06:00). Not reproducible from the current script as written (it only ever saves to
  `norris_fitted_{name}`) — an intermediate/draft output from that session, exact iteration not otherwise
  documented, kept as-is rather than deleted.
- `norris_fitted_GRB140206275.png/.pdf` — COMPLEX model's decorated final plot (physical units, LAT-photon
  overlay on twin axis, same layout as GRB131014A's).
- `norris_fitted_GRB140206275_simple.png/.pdf` — SIMPLE model's equivalent final plot.

GRB231129C:
- `fitter_GRB231129779.py` — the 5-pulse fit (this is already the *corrected* pulse count — see the
  "What it found for GRB231129C" paragraph above; a 4-pulse version was tried first and its `mc_kept_fraction`
  collapse under window-widening is what led to adding this 5th pulse). Fits over `[-1, 10]`s. Introduces
  **dominant-flux photon assignment** (later reused by both GRB140206B scripts above): the earlier
  nearest-preceding-`t_s` rule (used for GRB131014A) broke down specifically on this burst, where a fast
  pulse turning on just before a slower, still-dominant earlier pulse peaks would steal that earlier pulse's
  photons purely because its `t_s` was more recent, even though it barely contributed flux yet.
- `GRB231129779_lat.fits` — copied from
  `light_curves/GRB231129779/GRB231129C_filtered_gti_gtsrcprob_0.384_7.296.fits`; `T0_MET_S = 722977823.114`
  read from `LAT_analysis/GRB231129C/Ep1__0.384_3.136/*_fit_results_*.txt`'s own `T_0` line, cross-checked
  against the FITS file's own GTI the same way as GRB131014A's `T0_MET_S` was (see Validation above).
- `norris_fit_results_GRB231129779.csv` — one row per (pulse, episode) match, produced by
  `fitter_GRB231129779.py`.
- `norris_fit_GRB231129779.png` (narrow window, `[-1,10]`s) and `norris_fit_GRB231129779__bkp.png` (wide
  window, `[-10,20]`s) — **the user's own original manual before/after comparison** (2026-09-07, 03:02–03:13,
  predating the final `norris_fit_results_GRB231129779.csv`/`norris_fitted_GRB231129779.png` at ~05:14 by
  roughly two hours). This is the actual origin of the whole window-widening diagnostic above: the user
  noticed the fit moved visibly between these two plots, which motivated building
  `experiments/window_sensitivity_GRB231129779/window_sensitivity.py` as a controlled, quantified version of
  the same comparison (per that script's own docstring). Not reproducible from the current
  `fitter_GRB231129779.py` (which only saves to `norris_fitted_{name}`) — kept as the original artifact, not
  regenerated.
- `norris_fitted_GRB231129779.png/.pdf` — the decorated final plot.

GRB080916C (see the dedicated section above for the 8-pulse decomposition and ROOT cross-check):
- `fitter.py` — the live 8-pulse fit (`norris1`–`norris6` + `norris2-1`/`norris2-2`), the user's own
  working file, reused/repurposed across bursts as this track progressed (was GRB140206B's working file
  earlier, per GRB131014A's own `fitter.py` entry above). Its `P0` comment documents the full
  residual-seeding recipe and neutral-seed reproducibility results (see section above). Produces
  `norris_fit_results_GRB080916009.csv`.
- `fitter_EXPERIMENT.py` — Claude's diagnostic copy, added 2026-09-17 specifically to plot the 6-pulse
  data-minus-fit residual (raw + 5-bin-smoothed overlay) that motivated `norris2-1`/`norris2-2`; not the
  main fit, a supporting visualization for the section above. Produces
  `fitter_EXPERIMENT_residual_GRB080916009.png/.pdf`.
- `norris_fit_results_GRB080916009.csv` — one row per (pulse, episode) match, produced by `fitter.py`.
- `experiments/window_sensitivity_GRB080916009/` — the window-widening 6-vs-7-pulse diagnostic (see
  section above), added 2026-09-16: `window_sensitivity.py` (+ its own local `light_curves.py`/
  `norris_fit.py` copies), `window_sensitivity_results.csv`, and one fitted-light-curve plot per model per
  window width (`window_sensitivity_{six,seven}_{narrow_-1_70,wide_-20_150,widest_full_range}.png/.pdf`,
  6 plots total).
- `experiments/root_fit_GRB080916009/` — the CERN ROOT cross-check (see section above), added
  2026-09-17, moved into its own experiment folder 2026-09-18: `GRB080916009_lightcurve_10-400keV.csv`
  (NaI `n3+n4` summed, background-subtracted, 10–400 keV, quadrature-combined errors post-BUG-22-fix),
  `GRB080916009_norris_fit.C` (ROOT/Minuit macro, `use_errors` toggle for weighted vs. unweighted fitting),
  and its weighted/unweighted output plots (`GRB080916009_norris_fit_ROOT_{weighted,unweighted}.png/.pdf`;
  the no-suffix `GRB080916009_norris_fit_ROOT.png/.pdf` is a pre-`use_errors`-toggle duplicate of the
  weighted run, superseded but not deleted).

Shared infrastructure (all four bursts above):
- `norris_fit.py` — the local `NorrisFitter`/`norris_pulse`/`tv_value`/`tv_mc_summary` implementation every
  `fitter*.py` script now imports (`from norris_fit import ...`). Fixed 2026-09-16 (BUG-21, `BUGS.md`):
  this used to be the stray-dot `norris..py` (pre-Phase-5, `main-minor-75`), unimported by anything, while
  the scripts imported the broken `from variability_timescale.norris_fit import NorrisFitter` instead — now
  renamed and wired up, per `CLAUDE.md`'s "copy rather than fight `sys.path`" convention. `light_curves.py`
  (copied in alongside it) resolves the matching `variability_timescale.light_curves` import the same way.
  Fixed again 2026-09-18 (BUG-22, `BUGS.md`): `light_curves.py`'s `lightcurve_data()` combined per-channel
  errors linearly instead of in quadrature; patched here and in every copy of this file across this folder
  (the four `experiments/*/` copies plus the original `light_curves/light_curves.py`) to stay in sync.
- `dry_run.png` — a dry-run diagnostic plot (2026-09-07, ~02:57, the earliest timestamp of any file in this
  folder); which script/burst produced it is not documented here.
- `experiments/window_sensitivity_GRB231129779/` — the window-widening pulse-count diagnostic (see
  section above), `window_sensitivity.py` + its `window_sensitivity_results.csv` (3 windows × 5 pulses); a
  fitted-light-curve plot per window width (`window_sensitivity_narrow_-1_10.png/.pdf`,
  `window_sensitivity_wide_-10_20.png/.pdf`, `window_sensitivity_widest_full_range.png/.pdf`) was added
  2026-09-16, same pattern as the GRB131014A experiment above.
- `experiments/window_sensitivity_GRB140206275/` — the window-widening pulse-count diagnostic applied to
  both of this burst's candidate models (see section above), added 2026-09-16: `window_sensitivity.py`
  (+ its own local `light_curves.py`/`norris_fit.py` copies), `window_sensitivity_results.csv` (2 models ×
  3 windows × 5 or 7 pulses = 36 rows), and one fitted-light-curve plot per model per window width
  (`window_sensitivity_{simple,complex}_{narrow_-1_160,wide_-20_300,widest_full_range}.png/.pdf`, 6 plots
  total).
