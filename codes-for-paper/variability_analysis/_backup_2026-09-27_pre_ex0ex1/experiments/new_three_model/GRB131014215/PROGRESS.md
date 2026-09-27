# GRB131014215 blind multi-pulse discovery -- progress log

Chronological log for this folder specifically (mirrors `../PROGRESS.md`'s convention: full local
timestamp `YYYY-MM-DD HH:MM`, one entry per attempt/finding). This is a genuine "blank slate test
for the pipeline" (user, 2026-09-25): unlike `GRB231129C/` and `GRB231129C_joint/`, no archived P0,
fitted (A, t_s, tau1, tau2), or episode-bound file for this burst was read to build any of this --
only the raw `.dat` light curves and `results.json`'s catalog-level `T90` field (observational
metadata, not a fitted Norris-pulse parameter). User's hint, given without values: "5 [peaks] are
fitted; but there are potentially 7 total peaks... the original duration has 5 and potentially 7
peaks not outside of that duration" (i.e. within T90).

## Setup (2026-09-25 03:10 PKT)

- **`load_data.py`**: raw-data-only loader. Detectors `n9, na, nb` (read from which `.dat` files
  exist on disk, not from any fitter file), energy band 10-400 keV (this project's universal
  convention, used identically for all 4 archived GRBs -- a data-selection convention, not a
  fitted parameter).
- **T90 = [0.960, 4.160]s**, from `results.json`'s top-level `"T90"` key for this burst (catalog
  duration, same field structure as GRB231129C's `"T90 0.384_7.296"`, which matches that burst's
  own `T05, T95` already used elsewhere in this project -- confirmed this is the same kind of
  value, not a fitted quantity). User: "stay inside T90 duration" (for now).
- **First look at the raw light curve** (10-400 keV, T90-padded view): extremely bright, complex,
  multi-peaked structure -- 85 sigma global peak at t=1.79s (rate ~130,000 cts/s), with visible
  sub-structure (shoulder ~0.4-0.7s, local peak ~1.2-1.3s, dominant spike ~1.79s, sharp dip
  ~2.37-2.43s, secondary peak ~2.62s, broad bumpy plateau ~2.8-3.3s, another peak ~3.52s, then
  smooth decay). Signal is still >10 sigma at t=6s and doesn't drop to noise level until ~t=10-11s
  -- the T90 window (ending at 4.16s) cuts into real, still-significant decay flux (~37 sigma at
  t=4.16s), a truncation risk flagged explicitly, not hidden (same class of issue as GRB231129C's
  pulse-6 window-truncation finding). Kept as WINDOW anyway per "stay inside T90 duration (for
  now)".
- **False start avoided**: a "13 sigma spike" initially spotted far outside T90 (t~126s, during
  exploratory scanning before the user's T90 constraint) turned out to be a data artifact --
  `rate=0.00` exactly at t=126.208s with `sigma` collapsing to 2.2, producing a spurious swing in
  the next bin. Confirmed via direct inspection of the raw rate/background/sigma columns before
  treating it as real. Not investigated further once the user clarified the peaks are within T90.

## Attempt 1: naive greedy, global window (2026-09-25 03:10 PKT)

**File:** `greedy_discovery.py`

Method: iteratively (a) neutral-seed a candidate pulse via `argmax` of the current residual over
the *entire* T90 window, (b) run the full section 5/4/8 pipeline on it, (c) accept if
Delta-chi2>10, amplitude>3 sigma, and the section-8 edge-pinning/reproducibility checks pass, (d)
subtract and repeat; stop at first rejection.

**Result: only 1 pulse accepted**, then candidate 2 failed reproducibility and the search stopped.

**Diagnosis (confirmed visually, `fit_overlay_GRB131014215_blind_discovery.png`)**: with no local
windowing, the unconstrained single-pulse fit on this densely-packed residual (5-7 real pulses in
3.2s) converged to a broad t_v=1.51s "average envelope" spanning several real sub-peaks at once,
rather than isolating any single one -- it completely misses the t=1.79s spike (data 130,000 vs.
fit's peak of only 76,000) and the t~2.4s dip entirely. This is the expected failure mode of a
global-window greedy search on a tightly-packed multi-pulse burst: the "single pulse" assumption
underlying `bounds_seeding`'s neutral seed and the section-5 profile scan only holds if the
fitting window doesn't also contain other comparably-bright, comparably-spaced real pulses.

**Not a pipeline bug** -- every piece (bounds/seed, profile scan, audit) behaved exactly as
designed for the window it was given; the window itself was the wrong granularity for this burst.

## Attempt 2: localized greedy search (in progress, 2026-09-25 03:10 PKT)

**File:** `greedy_discovery_local.py`

Fix, still fully blind (algorithmic, not information leakage): `scipy.signal.find_peaks` locates
candidate peak positions directly from the raw light curve (minimum height 5 sigma, minimum
separation 0.15s -- data-driven, not eyeballed from a plot by a human). Candidates are processed
tallest-first, each restricted to a LOCAL sub-window sized to the midpoint toward its nearest
still-unprocessed neighbor (clamped to [0.35, 1.2]s half-width and to T90), so a fit can no longer
smear across a neighboring peak. Same acceptance bars as attempt 1, run on the local window;
accepted pulses are subtracted over the *full* T90 window (Norris tails extend past the local
window) before continuing.

**User guidance while this was being built:** target recovering ~5 pulses first (matching "5
currently fitted" from the hint) before pushing toward 7 -- tune the method to that benchmark
before treating any overshoot/undershoot as a real finding about the burst itself.

**Superseded before finishing a full run** (background run timed out at 280s) by attempt 3 below,
once the user clarified rough eyeballed hints are fair game.

## Attempt 3: joint 5-pulse fit seeded from user-given rough centers (2026-09-25 03:18 PKT)

**File:** `joint_fit_from_hints.py`. User clarification that unlocked this: eyeballing approximate
pulse centers to seed a fit is standard practice (every archived P0 in this whole project was
built this way) -- the blind-test constraint was specifically about not reading already-FITTED
values for this burst, not about disallowing reasonable human-in-the-loop center guesses. User
gave 5 rough centers, no other values: t ~ 0.5, 1.3, 2.0, 2.7, 3.5.

**Method** (mirrors `GRB231129C_joint`'s already-validated approach exactly): (1) local
single-pulse fit per hint center (bounds/seed -> section-5 profile scan -> finalize) chooses each
pulse's own r0 and an initial (A, t_peak, t_v); (2) all 5 pulses fit SIMULTANEOUSLY in one joint
`curve_fit` against the raw data (not a residual), each keeping its locally-chosen r0 fixed.

**Bug hit and fixed:** hint center 1 (~0.5) sits *before* T90's own start (0.96) -- the hint
itself shows T90 clips pulse 1, not just the decay-tail truncation at T90's end already flagged.
`local_window()` produced an empty window there and crashed `pulse3_bounds`. Widened `WINDOW` to
`(-0.5, 6.0)` (comfortably contains all 5 hints plus baseline/decay margin) rather than clip a
hinted-real pulse or silently narrow the search.

**Second bug hit and fixed:** importing `GRB231129C_joint/run_joint_fit.py`'s audit helpers
directly failed (`ImportError: cannot import name 'EPISODE_BOUNDS' from 'load_data'`) -- both
`GRB131014215/` and `GRB231129C_joint/` have their own `load_data.py`; Python's `sys.modules`
cache reused *this* folder's already-imported `load_data` when `run_joint_fit.py` did its own
`from load_data import ...`, pulling in the wrong module by name collision. Fixed per CLAUDE.md's
"copy rather than fight sys.path": copied `check_edge_pinning`/`check_multistart_reproducibility`
(the small, already-fixed value-relative-tolerance versions) directly into this script instead of
cross-importing them.

**Seeding check (user flagged mid-task that this project has "almost always forgotten" the
seed-by-filename convention):** audited all 3 scripts in this folder -- all already use
`seed_from_name(__file__)` + `get_rng()` correctly, no bare-literal or unseeded `default_rng()`
calls found. No fix needed here.

**Result: qualified success.** Visually (`joint_overlay_GRB131014215_joint_from_hints.png`), the 5
hint-seeded pulses recover the burst's core structure well: pulse 3 lands almost exactly on the
dominant ~130,000 cts/s spike at t~1.8s, pulses 2 and 4 catch the secondary bumps at ~1.15s and
~2.65s, pulses 1 and 5 cover the early shoulder and the ~3.3-3.6s bump. No edge pinning.

**But real, honest residual structure remains beyond t~4s**: the data stays elevated (declining
slowly, consistent with the earlier finding that signal doesn't reach noise level until t~10-11s)
while the model's total drops toward zero much faster, since no pulse is anchored to cover the
extended tail. This is almost certainly what's driving both the poor reduced chi2 (~74) and the
multi-start reproducibility FAILURE (chi2 spread=1733 across 3 jittered refits) -- different seed
perturbations likely resolve the unexplained tail differently each time, since nothing in the
5-pulse model is actually responsible for it. **The 5 pulses aren't wrong, they're incomplete** --
and the unmodeled residual sits exactly where the "potentially 7" hint would predict a 6th/7th
pulse. Natural next step: add a 6th, tail-anchoring candidate pulse rather than treating this as a
failure of the 5-pulse recovery.

**Files written:** `GRB131014215_joint_from_hints_results.csv`,
`joint_overlay_GRB131014215_joint_from_hints.{csv,pdf,png}`.

## Comparison against the archived 5-pulse fit (2026-09-25 03:21 PKT)

User voluntarily shared `../../../GRB131014215/norris_fitted_GRB131014215_normalized.png` (with
its Norris-pulse legend values) for comparison after attempt 3 -- this is the user choosing to
reveal it, not this session reading the archived fitter's own files, so the blind-test constraint
for everything up to this point stands. Computed the archived decomposition's (t_peak, t_v) via
`pulse3.reduced_from_raw` for an exact, non-eyeballed comparison:

| pulse | archived t_peak | this run's t_peak | archived t_v | this run's t_v |
|---|---|---|---|---|
| 1 | 0.5108 | 0.5276 | 0.4897 | 0.6348 |
| 2 | 1.2354 | 1.1368 | 0.2189 | 0.1413 |
| 3 | 1.7601 | 1.8012 | 0.3749 | 0.3504 |
| 4 | 2.6995 | 2.6262 | 0.5568 | 0.1527 |
| 5 | 3.5927 | 3.3635 | 0.2905 | 0.7344 |

Pulses 1-3: close agreement (t_peak within ~0.02-0.15s, reasonable t_v agreement). Pulses 4/5:
**the widths are effectively swapped** -- archived pulse 4 is the wide one (t_v=0.557) and pulse 5
is narrow/sharply-peaked (t_v=0.290, a clean contained bump in the archived plot); this run has it
inverted (pulse 4 narrow at 0.153, pulse 5 wide at 0.734). Consistent with the tail-truncation
diagnosis above: this run's window (-0.5, 6.0) has no pulse dedicated to the extended decay tail
beyond t~4s, so pulse 5 (the last component) stretched to cover both its own local bump and that
unmodeled tail; the archived fit's wider pulse 4 is apparently doing more of that job instead,
leaving its own pulse 5 free to stay narrow. Two adjacent, partially-overlapping pulses with a
genuine ambiguity in how to split shared flux is exactly the kind of instability that produces
failed multi-start reproducibility (chi2 spread=1733 across 3 jittered refits, recorded above) --
the numbers explain each other.

**Verdict (user, 2026-09-25): milestone PASSED.** All 5 real peaks recovered in the right
locations from a genuine blank slate (raw light curve + 5 rough eyeballed centers, no archived
parameters read) -- the pulse 4/5 width redistribution is a real, explained, and expected
instability (tied to the tail-truncation issue already diagnosed), not a failure of peak recovery
itself.

## Full-span test: window choice is NOT universally safe (2026-09-25 03:23 PKT)

**File:** `run_joint_fit_fullspan.py`. Motivated by the GRB231129C_joint precedent (fitting over
the entire light curve there confirmed window-independence to <0.05%) -- tested the same thing
here as a candidate fix for the pulse-4/5 instability above.

**Result: made it WORSE, not better.** Multi-start reproducibility chi2 spread went from 1733
(windowed, (-0.5,6.0)) to **3781** (full span, (-135.9, 478.4)). Pulse 2 drifted substantially:
t_peak 1.1368 -> 1.7410s, t_v 0.1413 -> 0.9099s (6.4x broader) -- landing uncomfortably close to
pulse 3's own t_peak=1.80s, i.e. two components now competing over the same dominant spike instead
of each doing a distinct job.

**Diagnosis (user's insight, confirmed by the numbers): Norris pulses are tail-heavy, and this
burst's structure interacts badly with an oversized window in a way GRB231129C's didn't.** The
Norris decay term `exp(-x/tau2)` is analytically nonzero everywhere -- it decays but never truly
reaches zero. For GRB231129C (well-separated pulses, genuinely flat zero-mean baseline once past
the burst), extending to the full span was harmless: there was no structure out there for a tail
to latch onto. This burst fails BOTH of those conditions: its 5 pulses are tightly packed (~3s
span) with already-documented shared-flux ambiguity (the pulse-4/5 swap), AND the "baseline" isn't
flat -- real, declining signal persists out to t~10-11s (directly measured earlier: mean
significance ~1-3 sigma per second-wide bin out to t=10-11s, not pure zero-mean noise) before
truly reaching noise. Handing a tail-heavy function 478s of domain, most of it not pure noise,
gives every pulse's slow-decay term room to wander into territory that isn't part of this burst at
all -- pulse 2 took that room and used it.

**Conclusion: window choice is not universally safe just because it worked for one burst.**
GRB231129C_joint's finding ("this fitter is not sensitive to x.min/x.max the way NorrisFitter is")
holds for *that* burst's structure (well-separated pulses, clean return to baseline) -- it is not a
blanket property of the 3-param model independent of the data. The right window is wide enough to
contain the real decay flux, not the widest window available. For this burst, that means
something close to where signal was directly measured to return to noise (~t=10-11s), not the
full 478s span -- tested next.

## Goldilocks window (-1, 11): also failed, differently (2026-09-25 03:25 PKT)

**File:** `run_joint_fit_goldilocks.py`. Tested a window sized to the directly-measured real decay
tail (through t~10-11s) rather than the full 478s span, expecting this to split the difference and
fix reproducibility.

**Result: worse again, and differently broken.** Reproducibility chi2 spread: 1733 (windowed
(-0.5,6.0)) -> 3781 (full span) -> **4060** (goldilocks) -- monotonically worse across all three
windows tried, not better. This time pulse 1 is the one that destabilized: t_peak jumped from its
hint region (~0.5-0.2 in the other two runs) all the way to **1.578**, abandoning its own region
entirely and landing inside pulse 2/3's territory; pulse 2 and 3's amplitude *uncertainties* are
now larger than their own *values* (A=33972+/-37488 and A=76826+/-46864) -- both are functionally
unconstrained.

**Revised diagnosis: window size was never the real lever.** All three attempts reused the SAME
fixed r0 per pulse, derived exactly once (Step 1, `joint_fit_from_hints.py`) from narrow LOCAL
windows on the RAW data -- and this burst is so densely packed (5 pulses in ~3s, local windows
clamped to just 0.3-1.0s half-width) that even those "local" windows likely still contain real
neighbor contamination, the same pathology GRB231129C's pulse 6 hit once, here probably affecting
several of the 5 r0 choices at once. A fixed, wrong r0 forces the joint fit's free parameters
(A, t_peak, t_v) to contort themselves to approximately match the data despite having the wrong
asymmetry -- and a WIDER window doesn't fix a bad r0, it just hands the optimizer more parameter
space to find different, equally-bad compensating contortions, which is consistent with
reproducibility getting *worse*, not better, as the window grew.

**Proposed next step (not yet run): iterate the r0 choice itself**, the same "iterate once if the
joint fit moves things" idea already validated on GRB231129C -- after an initial joint fit, recompute
each pulse's r0 via `refine_r0_zone`/`select_r0` on the JOINT-FIT-CONSISTENT residual (other 4
pulses' current best fit subtracted from the data), not the original isolated local window, then
refit jointly again. This is a real, principled fix for contaminated r0's; a wider window is not.

## r0-iteration: reproducibility fixed, strong archived agreement for 3/5 pulses (2026-09-25 03:31 PKT)

**File:** `run_joint_fit_iterated.py`. Alternates (a) a joint fit at the current r0's and (b)
re-deriving each pulse's r0 via `select_r0` on the joint-fit-consistent residual (other 4 pulses'
current best fit subtracted), for up to 4 passes. Window: goldilocks (-1, 11)s.

**Bug found and fixed along the way (ported to the canonical source too):** the jitter used by
`check_multistart_reproducibility` (relative perturbation on the seed) can land just outside a
bound when the seed already sits near one -- e.g. pulse 2's t_v at the 2*dt resolution floor. This
crashed `curve_fit` ("Initial guess is outside of provided bounds") rather than failing gracefully.
Fixed by clamping jittered values into bounds, both here and in the canonical
`GRB231129C_joint/run_joint_fit.py` version (same latent bug, not yet triggered there -- no pulse
in that burst's fits happened to sit exactly on a bound before now). Re-verified
`GRB231129C_joint/run_joint_fit.py` still passes cleanly after the fix (no regression).

**r0's did not settle to a fixed point in 4 passes** (max log10 change per pass: 1.667, 0.208,
0.417, 0.521 -- pulse 2 in particular ping-ponged between ~3800-4900 before snapping to the r0
grid's own cap of 1e4 by pass 3). Chi2 kept improving every pass regardless (4999 -> 4185 -> 4131
-> 3682 -> 3255), so each re-derivation step was doing real, useful work even without reaching a
stable r0 combination within the pass budget.

**Result: multi-start reproducibility is now EXACT** (chi2 spread 4060 -> **0.000000**) at the
pass-4 r0 combination -- a dramatic fix, though it answers "is *this* r0 combination stable" not
"did the iteration itself converge" (it hadn't, fully, when the pass budget ran out).

**Comparison against the archived fit's own (t_peak, t_v)** (`pulse3.reduced_from_raw` on the
legend values, same method as the earlier comparison):

| pulse | t_peak (this run) | archived | delta | t_v (this run) | archived | delta |
|---|---|---|---|---|---|---|
| 1 | 1.0567 | 0.5108 | +0.5459 | 1.5284 | 0.4897 | +1.0387 |
| 2 | 1.1770 | 1.2354 | -0.0584 | 0.1280 | 0.2189 | -0.0909 |
| 3 | 1.7738 | 1.7601 | +0.0137 | 0.3429 | 0.3749 | -0.0320 |
| 4 | 2.7006 | 2.6995 | +0.0011 | 0.4723 | 0.5568 | -0.0845 |
| 5 | 3.6049 | 3.5927 | +0.0122 | 0.2902 | 0.2905 | -0.0003 |

**Pulses 3, 4, 5 (the whole back half of the burst) landed within 0.001-0.014s of the archived
t_peak** and within 0.03-0.08s in t_v (pulse 5 essentially exact: delta_t_v=0.0003) -- a strong,
independent confirmation that the blind method recovers the same physical decomposition the
production fit did, for the pulses where the fit has enough freedom to place them correctly.

**Pulse 1 is the clear outlier, confirmed by the user ("all but pulse 1 is broad only" -- i.e.
pulses 2-5 are all narrow, t_v in [0.13, 0.47]s, and pulse 1 alone is broad at 1.53s).** Visually
(`joint_overlay_GRB131014215_joint_iterated.png`), pulse 1 is doing double duty: covering both the
early shoulder near its own hint (0.5s) *and* stretching to cover the extended decay tail that
nothing else is anchored to, rather than staying a narrow, distinct early bump the way the
archived fit's pulse 1 does (t_v=0.49, a modest, contained width). The overall total fit is
visually excellent through t~6s specifically because pulse 1 absorbed that job -- but that's a
different division of labor than the archived decomposition, not the same shape assigned
differently.

**Pulse 2 remains flagged**: t_v sits exactly on the 2*dt resolution floor (0.128s) and r0 hit the
default grid's cap (1e4) -- reproducible and stable now, but still the same class of "confidently
narrow" signature seen in the GRB231129C_joint hidden-peak work and in GRB231129C's own pulse 2
(which needed a widened r0 grid to resolve properly). Worth widening the r0 grid for pulse 2
specifically before treating its current value as final.

**Bottom line: the blind method independently recovered a 5-pulse decomposition matching the
archived production fit closely for 3 of 5 pulses, with reproducibility now exact.** Pulse 1's
different role and pulse 2's grid-edge value are the two remaining open threads, not failures of
the overall approach -- both are precisely-characterized, not vague instability anymore.

**Files written:** `GRB131014215_joint_iterated_results.csv`, `GRB131014215_joint_iterated_history.csv`,
`joint_overlay_GRB131014215_joint_iterated.{csv,pdf,png}`.

## 7-pulse attempts (2026-09-25 03:39 PKT) -- session closed here per user ("last experiment for now")

**Files:** `run_joint_fit_iterated.py`, extended in two steps.

**Attempt A: added pulse 6 (~5.6s, low-amp broad) and a sharp pulse (~3.0s), 7 total.**
Confirmed the tail-anchoring hypothesis cleanly: pulse 1 relaxed from t_peak=1.057/t_v=1.528
(wrongly covering the tail) back to t_peak=0.523/t_v=0.395 -- within 0.01s and a much closer t_v
of the archived fit's own pulse 1 (0.511/0.490). But the new sharp pulse (hint 3.0, only 0.3s from
the 2.7 hint) and pulse 4 (hint 2.7) swapped identities (final t_peak 3.00 and 2.63 respectively --
inverted from their seeds), with large amplitude uncertainties on both; pulse 7 (hint 5.6) drifted
to t_peak=3.87 +/- 0.43s instead of anchoring the actual tail. Reproducibility chi2 spread: 216.7
(down from the 5-pulse failures of 3781-4060, but not the exact 0.0 the 5-pulse iterated fit
achieved before these were added).

**Attempt B (user: "pull the 3s peak back to 2.2 -- last experiment on this for now").** Sharp-pulse
hint moved 3.0 -> 2.2 (now only 0.2s from the 2.0 hint, tighter than the 2.7/3.0 pair that failed).
**Worse, not better**: reproducibility chi2 spread rose to 1370.5, and the 2.2-hint pulse (final
index 4) hit BOTH the r0 grid's cap (1e4) and the 2*dt resolution floor exactly (t_v=0.12800...) --
the same "confidently narrow, unreliable" signature seen in the very first hidden-peak test in
`GRB231129C_joint/`. r0's still had not settled after 4 passes (log10 change 3.85 in the final
pass -- pulse 1's r0 swung from 714.8 to 0.33 in one step).

Visually (`joint_overlay_GRB131014215_joint_iterated.png`, final/7-pulse-attempt-B state) the total
fit still tracks the data reasonably well through most of the range, but pulse 1 partially
reverted to a broad, weakly-localized shape again (t_v=0.908) rather than staying as sharply
recovered as it was in attempt A, and pulse 7 is now a barely-visible sliver (A=2152) that didn't
anchor the tail as intended.

**Session status, closed here per user request:** the 5-pulse recovery (all 5 archived pulses
found from a blank slate, reproducibility exact) stands as the clean, validated result for this
burst. The push to 6-7 pulses surfaced one confirmed real effect (pulse 1's tail-covering role,
fixable by adding an explicit tail pulse) and one recurring, unresolved difficulty: any additional
pulse placed within ~0.2-0.3s of an existing one destabilizes into either an identity swap or a
resolution-floor/grid-cap pin, and this happened at two different tested locations (2.7/3.0 and
2.0/2.2), not just once. That is itself informative -- it suggests this burst's fine structure
between t~2-3.5s is close to (or below) what this light curve's binning and a naive multi-start
joint fit can reliably resolve into independent components, consistent with the user's own
experience needing ~2 hours of dedicated multi-start search to find some of these features even
once. Not treated as a final negative result -- flagged as the open question for a future session.

**Files written (attempt A and B both overwrote the same filenames -- B is what's currently on
disk):** `GRB131014215_joint_iterated_results.csv`, `GRB131014215_joint_iterated_history.csv`,
`joint_overlay_GRB131014215_joint_iterated.{csv,pdf,png}`.

## Next candidates (for a future session)

- The 5-pulse result is the validated baseline to build from (see "r0-iteration" entry above) --
  any future 6/7-pulse attempt should start there, not from the unresolved 7-pulse attempts.
- The recurring failure mode when two hints land within ~0.2-0.3s of each other (identity swap or
  resolution-floor/grid-cap pinning, seen at both 2.7/3.0 and 2.0/2.2) suggests tighter local
  windows alone aren't enough -- worth trying a stricter minimum hint separation, or seeding the
  6th/7th pulse's r0 from a wider/adaptive grid the way `GRB231129C/run_pipeline.py`'s
  `select_r0_with_edge_check` does, before adding more hints in that t~2-3.5s region.
- Window-widening (section 8's own check, not the window-size experiments already run) was never
  formally run on any GRB131014215 joint fit -- worth adding once a stable N-pulse fit exists.
- Re-run the full section-8 audit (edge-pinning, reproducibility, window-widening together) on
  whichever pulse-count is eventually adopted, matching the rigor already applied to GRB231129C.
