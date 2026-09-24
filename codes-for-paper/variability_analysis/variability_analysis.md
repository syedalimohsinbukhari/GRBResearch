# Variability timescale (t_v) via manually-seeded joint Norris fit — GRB131014A

Companion to `GRBResearchWork/PHASE5_TV_PLAN.md`, which documented the earlier automated
(MEPSA/scipy peak-detection + per-window local fit) pipeline for this same Phase 5 goal. That pipeline's
source files (`variability_timescale/norris_fit.py`, `variability_timescale.py`, etc.) have since been
deleted — see `PHASE5_TV_PLAN.md`'s 2026-09-16 status update — so this folder is now the only live track.
This folder is a **separate, manually-driven track**, worked one GRB at a time starting with GRB131014A on
2026-09-07 after the automated pipeline stalled on GRB080916C's TR2 (zero MEPSA detections even with
padding, because TR2's true shape is one broad, smoothly-declining pulse with no local excess for a spike
detector to find — see the MEPSA-fallback discussion in this session's history). GRB140206B
(`fitter_GRB140206275.py` + a `_simple` variant) and GRB231129C (`fitter_GRB231129779.py`) were fit in
this same manual style in a later session (2026-09-15/16), and GRB080916C (`fitter.py`) followed
in a session on 2026-09-16/17, originally as an 8-pulse fit (`norris1`-`norris6` plus a `norris2-1`/
`norris2-2` split of `norris2`'s tail) — **reverted 2026-09-22, user decision**, back to the standard
7-pulse model (`norris1`-`norris6` plus one extra pulse at $t_s\approx20$s inside TR3) already
independently validated via the window-widening/residual-SSE diagnostic in
`experiments/window_sensitivity_GRB080916009/`; see the "Reversion" note atop the dedicated section below
for the (now-superseded) 8-pulse work — **none of the three are yet written up below** — the "What the code
computes" / "Every judgement call" / "Results" sections that follow describe the GRB131014A fit only; see
each script's own docstring/comments and the sections further down for the other three bursts until this
note is extended. Not yet wired into `lorentz_factor.py`'s `Gamma_min` pipeline — same deliverable-boundary
stance as the abandoned automated pipeline.

**Open cross-burst issue, flagged 2026-09-18 — resolved 2026-09-23.** Every fit in this folder (all four
bursts) used to use only a single NaI detector (`nai_data[0]`) — documented at the time as `n3` for GRB080916C,
`na` for GRB131014A, `n3` for GRB140206B, `n7` for GRB231129C — rather than summing across the GRB's full
NaI detector set, unlike the convention typical of published GRB analyses. This wasn't a deliberate choice
recorded anywhere — it's just what every `fitter*.py`/`window_sensitivity.py` script had done since the
very first one. **The "whichever sorts first alphabetically" part of that original note was itself
wrong — see BUG-23's cross-burst extension below**, found while verifying whether the other three bursts
share GRB080916C's detector-selection bug.

**Fix (2026-09-23, user decision): sum every NaI detector, raw, no per-detector normalization.** Rather than
patch the single-detector pick (which would have meant hardcoding a per-burst detector name, per BUG-23's
original writeup below), every `fitter*.py`/`window_sensitivity.py` script now sums all of a burst's NaI
detectors' background-subtracted count rates directly — each detector's own effective area/viewing angle is
trusted to weight its own contribution, rather than peak- or integral-normalizing each detector to equal
weight first, or applying an angle/effective-area-weighted combination. This was an explicit choice among
those alternatives, not a default — the user picked raw summing specifically because it matches the one
existing precedent already in this file (the GRB080916C ROOT cross-check's "n3+n4 summed" light curve,
described further down). It also closes BUG-23 structurally rather than patching it: there is no longer a
single detector to mis-pick, so `os.listdir()` order can never matter for any burst, present or future.
Each script asserts every detector's time grid matches the first before summing (confirmed identical for
all four bursts' full detector sets by direct check, not assumed, across all 4 files' data). Applied to all
9 affected scripts (5 `fitter*.py` + 4 `window_sensitivity.py`); all compile.

**All four bursts re-run, 2026-09-23 — all converge; every `norris_fit_results_*.csv`/plot in this folder
now reflects the summed-detector fit, not the pre-fix single-detector one.**

- **GRB231129C**: converged cleanly, existing 5-pulse seed unchanged. `mc_kept_fraction` improved for every
  pulse (e.g. 0.35→0.72, 0.94→0.99) — a genuine S/N gain from combining detectors here, not dilution.
- **GRB140206B**: both the SIMPLE/COMPLEX pair (`fitter_GRB140206275.py`) and the `_simple.py` variant
  converged cleanly, existing seeds unchanged.
- **GRB131014A**: initially failed (`RuntimeError: Optimal parameters not found`) with both the original
  seed and one reconverged from the old single-detector fit. Diagnosed as a pure iteration-budget issue —
  both seeds land on the identical solution once `NorrisFitter`'s default `max_iterations=5000` is raised to
  20000 (the summed 3-detector curve peaks at ~130,030 cts/s vs. the single-detector curve's ~44,084, needing
  more optimizer steps to settle). Fixed with `max_iterations=20000` in `fitter_CLAUDE_GRB131014215.py`
  (also applied to `fitter.py` as a precaution). Converged parameters match the old single-detector fit
  closely (`t_peak`s agree to <0.01s), confirming this was purely a budget issue, not a different optimum.
- **GRB080916C**: converged, but pulse 7 (`t_s≈20`, inside TR3, the pulse added in the 2026-09-22 7-pulse
  reversion) landed on a genuinely different, degenerate optimum on the summed curve: `tau2` collapsed from
  1.02 to 0.0087 and `mc_kept_fraction` from 0.52 to 0.05 — a visible near-delta-function spike in the plot
  at t≈20.5s. Checked against the raw summed data directly: that region sits at ~1650-1700 cts/s against a
  ~1000-1800 cts/s noisy floor throughout the whole t=17-23s window, i.e. no real standout feature — the fit
  was overfitting two noisy bins. Confirmed this wasn't a seed or iteration-budget artifact: the original
  seed and one reconverged from the old single-detector fit land on the same degenerate solution, at every
  iteration budget tried up to 50000. **User decision (2026-09-23): drop this pulse**, reverting to the
  6-pulse model (`P0_6`, already an independent baseline in
  `experiments/window_sensitivity_GRB080916009/window_sensitivity.py`'s own 6-vs-7-pulse check) rather than
  bound `tau2` or reseed via the residual-subtraction technique used elsewhere in this file. Verified before
  applying: dropping the pulse leaves the other six essentially undisturbed — every other pulse's
  `t_peak`/`t_v`/`A` matches the (still-degenerate) 7-pulse fit's corresponding pulse within ordinary
  `t_s`<->`tau1`-degeneracy noise (largest `t_peak` shift 0.03s), and total SSE rises only ~0.4% (7.202e7 vs
  7.173e7, physical-units counts/s²) — the expected signature of removing a pulse that was fitting noise
  rather than signal. Applied to `fitter.py`; the regenerated plot confirms the spike is gone and the other
  six pulses track the data the same way as before.

**Checked, not assumed: GRB080916C's LAT-photon-to-pulse assignment is unaffected by dropping pulse 7.**
Compared the printed photon-arrival/pulse table between the 7-pulse and 6-pulse runs directly — identical,
all 14 photons assigned to pulse 3 or pulse 4 in both. The removed pulse never had any photons assigned to
it to begin with (its `t_s≈20` onset always lost the nearest-preceding-onset/active-threshold competition to
pulse 4's broad, already-active TR3 pedestal), so dropping it changes nothing about the photon table.

**Not yet done:** none of these four re-fits have been fed into `lorentz_factor.py`'s `Gamma_min` pipeline —
same deliverable-boundary stance as always (§ above).

**Standalone follow-up, 2026-09-23: `experiments/normalized_vs_unnormalized_fit/`.** Prompted by the
question of whether this folder's `y /= Y_MAX_CTS_PER_S` normalization step (used by every `fitter*.py`)
changes the fit result. Full writeup, method, and results in that folder's `comparison.md`; short version:

- Tested a from-scratch unnormalized alternative (raw counts/s, `scipy.optimize.least_squares` with an
  explicit `x_scale` array, bypassing `pymultifit` entirely) against the production normalized method, for
  all four bursts. **Normalization changes nothing about the physics** — SSE agrees to within 0.15% for
  every burst, most pulse parameters agree to <1%. `comparison.md`'s companion file,
  `why_not_unnormalized_fitting.md`, is a standalone reviewer-facing answer to "why wasn't the raw light
  curve fit directly", kept as a backup in case this is ever asked about the paper's methodology.
- The comparison also independently corroborated the GRB080916C pulse-5 drop above (that pulse diverged
  sharply between methods before it was known to be degenerate in production) and led to a GRB140206B
  follow-up: widening the fit window to the light curve's own `t.min()`/`t.max()` resolves pulse 7's
  window-edge pinning but destabilizes pulse 6 instead (same `tau2`-collapse pathology as GRB080916C's
  dropped pulse); re-seeding pulse 6's `(tau1, tau2)` from `(83, 1.1)` to a neutral `(1, 1)` at amplitude
  `A=0.2` (found by trial — `A=0.23` converges too, but lets the pulse drift out of its intended
  $t\approx23$s region and destabilizes two other pulses instead) fixes both pulses 6 and 7 simultaneously,
  with SSE matching the original-seed fit (not an artificially lower one bought with a worse
  decomposition). This is recorded in `comparison.md` as a **diagnosed, recommended configuration for
  GRB140206B — not yet applied to `fitter_GRB140206275.py`/`_simple.py`**, both of which still use their
  original `(-1, 160)` window and converge cleanly there. Pulse 5 remains a known, pre-existing (not newly
  introduced) weak point at the full-range window, unresolved.
- **Paper-facing follow-up, superseded 2026-09-24 — see the new section below.** `full_range_all_bursts_check.py`
  (same folder) was split into `full_range_single_burst.py` (shared logic) + one driver per burst
  (`full_range_check_GRB*.py`), specifically so each burst could run as its own process instead of one
  script blocking on all four — the combined script had been killed twice before finishing. Two of the
  four (GRB080916C, GRB140206B) ran to completion this way and **failed** the user's ≤5% `t_peak`/`t_v`
  tolerance (see below); GRB131014A/GRB231129C were killed mid-run once the session's focus moved to the
  per-burst `GRB*/` folders described next. Full detail in the new section below.

## Session update, 2026-09-24: per-burst normalized-vs-unnormalized full-range treatment

Builds on the two threads above (the `≤5%` full-range tolerance check and
`experiments/normalized_vs_unnormalized_fit/`), but promotes the comparison out of the `experiments/`
sandbox into a proper per-burst treatment with the full production machinery — episode matching,
LAT-photon-to-pulse assignment using each burst's own already-established rule, MC-propagated `t_v`, and
CSV+PNG+PDF outputs — rather than `fit_comparison.py`'s lighter point-estimate-only comparison. **The
user is planning to continue this work on a heavier machine** — this section is written so a fresh
session (here or there) can pick it up without re-deriving anything. **All of this — the split scripts,
the four new `GRB*/` folders, the standalone `fitter_GRB080916C.py`, the `window_sensitivity_GRB140206275`
split, and `t_peak_tv_reparametrization_idea.md` — is exploratory/testing work, and the user intends to
stage it for a git commit**, an explicit exception to this project's usual working-tree-only rule (same
kind of exception already made for §13's RNG/seeding overhaul and the 2026-09-23 Phase 5/BUG-23 commit —
see those sections for precedent). None of it has fed `lorentz_factor.py`'s `Gamma_min` pipeline yet, and
two of four bursts already failed the ≤5% full-range tolerance check below, so "staged for commit" here
means "this testing work is being checkpointed," not "these results are settled or paper-ready."

**`full_range_check_GRB*.py` (split from `full_range_all_bursts_check.py`), run to completion for two of
four bursts:**

- **GRB080916C**: `t_peak` stable across all 6 pulses (≤0.3%), but `t_v` exceeded the 5% tolerance on 4
  of 6 pulses (up to ~17%). **Fails** the ≤5% gate.
- **GRB140206B**: pulses 1/2/3/7 fine; pulses 4/5/6 (the already-fragile 23–28s cluster) blow far past
  tolerance — pulse 6's `t_peak` shifts 37%, `t_v` shifts 86%, even with the pulse-6 reseed from
  `experiments/normalized_vs_unnormalized_fit/comparison.md`. **Fails** the ≤5% gate.
- **GRB131014A, GRB231129C**: killed mid-run (`TaskStop`) once the session's focus moved to the per-burst
  `GRB*/` folders below — no result recorded, per this project's own "an interrupted run proves nothing"
  convention.
- **Net effect**: two of four bursts have now failed the full-range tolerance check outright, so a
  blanket switch of `lorentz_factor.py`'s `Gamma_min` input to the full range is **not** viable as a
  single decision — this looks like it will need to be a burst-by-burst (or even pulse-by-pulse) call,
  not a project-wide toggle.

**New per-burst folders, `GRB080916C/`, `GRB131014215/`, `GRB140206275/`, `GRB231129C/`** (siblings of
this file, each with `_common.py` + `fitter_normalized.py` + `fitter_unnormalized.py`). Each `_common.py`
loads the burst's full-range summed light curve once, reuses that burst's own established P0 and
photon-assignment rule (nearest-preceding-active-onset for GRB080916C, plain nearest-preceding-onset for
GRB131014A, dominant-flux for GRB140206B/GRB231129C — copied from each burst's own production
`fitter*.py`, not standardized across bursts), and exposes a `FitResult` duck-type class so
`norris_fit.tv_mc_summary()` works identically whether the fit came from a real `NorrisFitter` (normalized
method) or a bare `scipy.optimize.least_squares` result (unnormalized method — covariance derived from the
Jacobian the same way `scipy.optimize.curve_fit` does internally, `absolute_sigma=False` convention, so
both methods' `t_v` uncertainties come from the same kind of estimator). `AutoMinorLocator` is set
explicitly (same `n` on both twinned y-axes) so the count-rate/photon-energy axes' minor gridlines
actually align — `update_style()`'s default minor-tick behavior picks a different subdivision count per
axis depending on each axis's own major-tick step, so without this fix the two axes' minor ticks drift
out of alignment even when the (explicitly-set) major ticks match.

- **GRB080916C**: both methods run and agree closely — all pulses within ~2% on `t_peak`/`t_v` except the
  already-known-fragile pulse 5/TR4, whose raw `tau1` differs by 52% between two different P0 seeds tried
  in the standalone `fitter_GRB080916C.py` (see below) while `t_peak` moved only 0.05% and `t_v` only
  ~1–1.7% — a clean illustration of the `t_s`/`tau1`/`tau2` sloppy-parameter behavior discussed in
  `t_peak_tv_reparametrization_idea.md`.
  - Also promoted to a **standalone, self-contained top-level script**, `../fitter_GRB080916C.py` (no
    `_common.py` import, matching every other `fitter_*.py`'s pattern per `CLAUDE.md`'s "copy rather than
    fight sys.path") — same 6-pulse `P0`, full range, outputs suffixed `_fullrange` so they never collide
    with `fitter.py`'s own `(-1, 70)`-window results. The user has been iterating `P0` directly in this
    file since; treat its current on-disk `P0` as a live, possibly-mid-edit value, not a settled one.
- **GRB131014A**: both methods run (slow — the largest dataset of the four, ~21 min/~13 min CPU time for
  normalized/unnormalized). Pulses 3/4/5 (TR1/EX1/EX1) agree well; **pulses 1 and 2 disagree
  substantially between methods** (`t_peak` off by 128%/10%, `t_v` off by 4%/76%) — but this reproduces,
  not contradicts, this file's own already-documented finding above that pulses 1/2 are genuinely
  **unresolved** for this burst (two different seedings already gave two different answers within the
  single-method, single-window production fit; see the GRB131014A "Known limitations" section further
  down). Two different *methods* landing on two different answers here is the same underlying degeneracy
  showing up again, not a new discrepancy introduced by the normalized/unnormalized split.
- **GRB140206B**: COMPLEX (7-pulse) model only — the decomposition `fitter_GRB140206275.py` itself uses
  for its main per-episode results. Fit at the **true** full `x.min()/x.max()` range, deliberately, not
  the `(-20, 300)` partial widening `experiments/normalized_vs_unnormalized_fit/fit_comparison.py` used
  for this burst specifically (that partial widening existed only to avoid pulse 7's window-edge pinning
  at `(-1, 160)` — this folder's `_common.py` flags the true full range as a known risk for exactly that
  reason before it was run).
  - **Pulses 2 and 7 initially looked catastrophically discrepant (821%/89% on `t_peak`) — this is a
    pulse-index label swap between the two fits, not a real disagreement.** Confirmed by matching
    physical identity instead of raw index: normalized-pulse-2 (`t_peak=13.19, t_v=7.10`) matches
    unnormalized-**pulse-7** (`t_peak=13.25, t_v=7.09`) to <0.5%; normalized-pulse-7 (`t_peak=121.68,
    t_v=22.62`) matches unnormalized-**pulse-2** (`t_peak=121.50, t_v=22.80`) to <1%. The unnormalized
    optimizer converged with the ~13s pulse and the broad pedestal swapped in array position relative to
    `COMPLEX_P0`'s original order — worth remembering before trusting any raw index-to-index diff on a
    joint multi-pulse fit across two independently-converged methods.
  - Once corrected for the swap, pulses 2/3/4/6/7 all agree well (`t_peak`/`t_v` within a few percent).
    **Pulse 5 (TR3) has a real, non-swap 21% `t_v` disagreement** (1.089s normalized vs. 0.864s
    unnormalized) despite similar `mc_kept_fraction` (~0.78) in both methods — genuinely open, not yet
    diagnosed. Worth checking first on the heavier machine.
  - Pulse 1 shows a 101% `t_peak` difference but is trivial in absolute terms (0.045s vs 0.091s, both
    essentially at the trigger) — checked against every other pulse and ruled out as a swap.
- **GRB231129C**: folder built (`_common.py` + both drivers), **not run to completion in this track**.
  Superseded mid-session: the user found GRB231129C's light curve/detector set had changed significantly
  (more NaI detectors than previously accounted for) and is handling that burst's refit manually, directly
  in production `fitter_GRB231129779.py` (working-tree edits as of this writeup: `START1/END1` widened to
  `(-inf, inf)`; `P0` replaced with converged parameters read off the `window_sensitivity_GRB231129779`
  widest-window run) rather than through the `GRB231129C/` folder. **As of this writeup that file is
  mid-interactive-edit** (`nf.dry_run(); plt.show()`, no `nf.fit()` call active) — not a finished result.
  The `GRB231129C/fitter_normalized.py` driver in the new folder was independently left mid-edit too
  (debug prints, a temporary 4-pulse `P0`, `plt.show()`) by the same exploration — neither should be
  trusted as settled until re-run cleanly.

**`experiments/window_sensitivity_GRB140206275/window_sensitivity.py` also split** the same session, into
`window_sensitivity_common.py` (shared light-curve loading + a `run_model(model_name, p0, also_track=...)`
function) plus `window_sensitivity_simple.py`/`window_sensitivity_complex.py` — same reasoning as the
`full_range_check_GRB*.py` split (independent, individually re-runnable). The original combined script and
its pre-split outputs are preserved verbatim in `_backup_2026-09-24_pre_split/` before removal. Re-running
COMPLEX reproduced the pre-split findings essentially exactly (pulse 6 the fragile one, `mc_kept_fraction`
0.41→0.37→0.58 non-monotonically across narrow/wide/widest; pulse 7 the actual pedestal, healthy monotonic
improvement 0.90→0.97→0.98) — confirms the split didn't change anything, as intended.

**Speculative, not scheduled:** `t_peak_tv_reparametrization_idea.md` (this folder) sketches refitting the
Norris pulse directly in `(t_peak, t_v, r=tau1/tau2)` instead of `(t_s, tau1, tau2)`, with a full
closed-form derivation of the change of variables, and explains why it targets a different (more central)
degeneracy than this project's earlier `(tau, xi)` reparametrization attempt (which left `t_s` as a direct
fit parameter and never actually decoupled the dominant `t_s`<->`tau1` degeneracy — see that section
further down). Written to solicit outside review before deciding whether it's worth a real trial; not
implemented, not on any to-do list.

**BUG-23 extended to all four bursts, 2026-09-22 (verification, no fix applied yet — see `BUGS.md`).**
Checked every burst's *current* `os.listdir()[0]`-picked detector against its *committed* CSV's
`y_max_cts_per_s` (which pins the detector actually used when those results were generated, cross-checked
against each candidate `.dat` file's own peak count rate):

| GRB | documented/committed detector | current `os.listdir()[0]` | live-broken right now? | would `sorted()` fix it? |
|---|---|---|---|---|
| GRB080916C | `n3` | `n4` (pre-fix) | was, now fixed | yes — but only because `n3 < n4`, a coincidence, not a principled fix |
| GRB131014A | `na` | `nb` | **yes, unfixed** | no — `sorted(['n9','na','nb'])[0] = 'n9'`, still wrong |
| GRB140206B | `n3` | `n0` | **yes, unfixed** | no — `sorted(['n0','n1','n3'])[0] = 'n0'`, still wrong |
| GRB231129C | `n7` | `n7` | no, currently correct | **no — would break it**: `sorted(['n3','n6','n7'])[0] = 'n3'` |

Same root cause as GRB080916C's original BUG-23 (`os.listdir()` order is filesystem-dependent, not
alphabetical), but the `sorted()` fix applied to GRB080916C does **not** generalize: for GRB131014A and
GRB140206B, the documented/correct detector is not the alphabetically-first one among that burst's actual
detector set, so `sorted()` would keep them broken (131014A) or leave them broken (140206B) without
actually fixing anything, and applying it to GRB231129C would break a burst that is currently fine. The
correct fix is to hardcode each burst's documented detector name explicitly, not to sort. **Affects 6
files, not yet touched:** `fitter_GRB140206275.py`, `fitter_GRB140206275_simple.py`,
`fitter_GRB231129779.py`, and the three matching `experiments/window_sensitivity_GRB*/window_sensitivity.py`
copies (`fitter_CLAUDE_GRB131014215.py` and GRB131014A's own `experiments/window_sensitivity_GRB131014215/`
copy have the same unsorted pattern too, confirmed by inspection, though GRB131014A's case is already
covered by the table above). GRB080916C's own `sorted()` fix should also be revisited to an explicit pin,
since it is currently correct by coincidence rather than by design.

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
the 8-pulse decomposition (later reverted — see the dedicated section below).

**Re-verified 2026-09-22** after the BUG-23 detector fix (see the dedicated section below): rerunning
`window_sensitivity.py` reproduces this section's numbers essentially exactly (SEVEN pulse 5's
`mc_kept_fraction` 0.5172→0.6037→0.6039 vs. the 0.517→0.604→0.604 documented above) — the findings in
this section stand, this was purely a detector mix-up introduced by an unrelated data resync, not a
change to the underlying result.

## GRB080916C: 8-pulse decomposition via residual-seeding, and a CERN ROOT cross-check (2026-09-16/18)

**Reversion, 2026-09-22 (user decision): `fitter.py` no longer uses this 8-pulse decomposition.** The user
judged the `norris2-1`/`norris2-2` split (below) not sound, and `fitter.py` now runs the standard 7-pulse
model instead — identical to the 6-pulse fit below plus one extra pulse at $t_s\approx20$s inside TR3's
window, the same model already independently validated (not via this residual-seeding argument, but via the
window-widening/residual-SSE diagnostic) in `experiments/window_sensitivity_GRB080916009/` (see the
"6-vs-7-pulse question" section above). `light_curves/GRB080916009/`'s NaI `.dat` files were missing on this
machine (`iqra-siddique`) earlier in this same session — see the note above — but were synced in
mid-session; `fitter.py` has now been rerun and `norris_fit_results_GRB080916009.csv` (7 rows, one per
pulse) reflects the 7-pulse fit as of 2026-09-22. Also added in the same session: a LAT-photon overlay
(twin y-axis, red hollow circles, only $E>1$ GeV) on `fitter.py`'s final plot, per user request, matching
`fitter_GRB231129779.py`'s convention; the LAT FITS file was copied locally
(`GRB080916009_lat.fits`, from `light_curves/GRB080916009/lat.fits`) rather than referenced live, same
"copy rather than fight `sys.path`" convention as that file. `T0_MET_S = 243216766.62`, read from
`LAT_analysis/018__GRB080916009/Ep1A__m0.128_4.864/GRB080916C_fit_results_-0.128_4.864.txt`'s own `T_0`
line and cross-checked against the FITS file's own GTI (1.280–64.257s matches the T90 window exactly),
same validation pattern as the other three bursts.

**BUG-23 (see `BUGS.md`), found and fixed 2026-09-22 while re-running the window-sensitivity check below
at the user's request:** the same NaI `.dat` resync that unblocked the above silently made `fitter.py`'s
and `window_sensitivity.py`'s `dat_NaI[0]` pick `n4` instead of the documented `n3` — `os.listdir()`
order is filesystem-dependent, not alphabetical, and the resync happened to change it. This affected the
CSV/plot regeneration described in the paragraph above (both had to be rerun a second time after the
fix); `y_max_cts_per_s` in the current `norris_fit_results_GRB080916009.csv` (1879.6628) confirms `n3` is
now the detector in use. Fixed via `sorted()` in both scripts. **Not yet checked for the other three
bursts** — same scope boundary as the cross-burst single-detector issue above.

**Photon-to-pulse assignment: nearest-preceding-onset (temporal-proximity), not dominant-flux — user
decision, 2026-09-22.** `fitter.py` now assigns each $E>1$ GeV photon to whichever pulse has the largest
$t_s \le$ the photon's arrival time, same rule and same `assign_pulse()` implementation as
`fitter_CLAUDE_GRB131014215.py` (GRB131014A), not the dominant-flux rule used for GRB231129C/GRB140206B.
**Why, for this burst specifically:** the two nearest $E>1$ GeV photons to norris3's peak (5.869s) arrive
at 6.072s and 6.857s — only 0.20s/0.99s after it, and inside TR2 (norris3's own episode). The same two
photons sit 3.43s/4.21s after norris2's peak (2.646s, in TR1), yet dominant-flux would still assign them
to norris2, because norris2's broad tail ($\tau_2\approx6.17$) out-predicts norris3's small, narrow pulse
at that time — the same overshadowing behavior already documented for norris2 against the (now-reverted)
8-pulse split's `norris2-1`/`norris2-2`. Temporal proximity assigns both photons to norris3 instead,
matching the episode boundary and the physical expectation that a high-energy photon arrives during or
shortly after the pulse that produced it.

**Found and fixed the same session, user decision:** applied burst-wide with pure onset-recency, the rule
above assigned *every* photon from 6.072s up to pulse 5's onset (19.729s) to pulse 3 — six photons total,
including three (10.215s, 16.538s, 16.798s) that sit well past norris3's own decay and deep inside pulse
4's broad TR3 pedestal, because pulse 4's earlier onset ($t_s=1.010$s) can never "supersede" pulse 3's
later one under pure onset-recency. User caught this by inspection ("6 photons kaise; pulse 3 k kareeb to
sirf 3 photons hain") and asked for a refinement.

**Fix: `ACTIVE_THRESHOLD_FRAC` restricts candidacy to pulses that are still "active."** A pulse is only a
temporal-proximity candidate at time $t$ if its own predicted value is $\ge1\%$ of its own peak value at
that instant; among the pulses that pass this filter, the same largest-preceding-$t_s$ tie-break as before
picks the winner. 1% was picked for margin, not tuned: norris3's own value falls from 11.4% of its own
peak at 7.445s to 0.076% at 10.215s — a $>2$-orders-of-magnitude drop — so any threshold roughly between
0.1% and 10% gives the identical result. Verified directly against each pulse's own predicted flux at
every candidate photon time (not asserted from the formula alone) before adopting it.

**Result: pulse 3 now gets exactly the 3 genuinely-nearby photons** (6.072s, 6.857s, 7.445s — 0.2–1.6s
after its own peak), and its "defining" (max-energy) photon is now 2110.10 MeV @ 6.857s, one of the two
photons that actually motivated the temporal-proximity argument in the first place — resolving the
mismatch the unrestricted rule had introduced. The three displaced photons (10.215s, 16.538s, 16.798s) now
go to pulse 4 (physically sensible — pulse 4's broad pedestal is what's actually producing flux there),
along with the two latest photons (40.502s, 43.992s) that the unrestricted rule had also misassigned to
pulse 5 well after pulse 5's own decay. Pulse 5 keeps the 5 photons genuinely inside its own active window
(22.209–28.203s).

New CSV columns (`norris_fit_results_GRB080916009.csv`): `n_lat_photons_assigned`, `photon_e_max_MeV`,
`photon_t_arr_at_e_max_s`, one row per pulse — pulse 3: 3 photons, max 2110.10 MeV @ 6.857s; pulse 4: 6
photons, max 27428.80 MeV @ 40.502s; pulse 5: 5 photons, max 6721.31 MeV @ 28.203s; every other pulse: 0.

**TR2's candidate for downstream analysis: pulse 3 (norris3), decided 2026-09-22 — user decision.** TR2
has only one pulse matching its boundary (4.864–15.040s) in the first place — pulse 3, $t_\text{peak}=5.869$s
— so there was no TR2-internal ambiguity to resolve (unlike TR3, below, which has two). The decision
recorded here is to use it, on the strength of its now-corrected photon assignment: its defining photon
(2110.10 MeV @ 6.857s) arrives only 0.99s after its own emission peak, the closest peak-to-photon gap of
any pulse/episode pairing worked out in this file so far. **Not yet acted on**: TR2's $t_v$ in
`lorentz_factor.py` is still duration-sourced (see below) — this only fixes *which* pulse TR2 will draw
from once that wiring happens.

**TR3's candidate is still an open decision, not resolved by the above.** TR3 (15.040–55.296s) has two
matching pulses — pulse 4 (broad pedestal, $t_\text{peak}=26.273$s, `kept=0.982`, $t_v=16.08$s) and pulse 5
($t_\text{peak}=22.991$s, `kept=0.517`, $t_v=1.31$s) — and which one (if either) should feed `Gamma_min` for
TR3 has not been decided.

**Known naming bug, pre-existing, not introduced by this revert:** `fitter.py`'s `fig_path`
(`fig_path = Path(__file__).parent / f".norris_fitted_{GRB_080916C.name}"`) has a stray leading dot, so
the plot actually saves as the hidden files `.norris_fitted_GRB080916009.png/.pdf`, not
`norris_fitted_GRB080916009.png/.pdf` as named throughout this document and in every other burst's
equivalent line — not fixed here since it wasn't part of what was asked; flagging so it isn't mistaken for
a missing file.

The CERN ROOT cross-check below was also run against the now-superseded 8-pulse fit — the BUG-22
error-combining fix it surfaced is independent of pulse count and remains valid, but the ROOT macro's own
fitted values are not representative of the current (7-pulse) model. The rest of this section is kept
as-is, unedited, as the historical record of how the 8-pulse decomposition was reached.

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

GRB080916C (see the dedicated section above for the 8-pulse decomposition and ROOT cross-check, and its
"Reversion" note for the 2026-09-22 change described here):
- `fitter.py` — **as of 2026-09-22, the live 7-pulse fit** (`norris1`–`norris6` plus one extra pulse at
  $t_s\approx20$s inside TR3), reverted from the 8-pulse `norris2-1`/`norris2-2` split per user decision.
  The user's own working file, reused/repurposed across bursts as this track progressed (was GRB140206B's
  working file earlier, per GRB131014A's own `fitter.py` entry above). Produces
  `norris_fit_results_GRB080916009.csv`, regenerated at the 7-pulse model 2026-09-22. Its plot output
  actually lands at the hidden filenames `.norris_fitted_GRB080916009.png/.pdf` — see the naming-bug note in
  the dedicated section above.
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
