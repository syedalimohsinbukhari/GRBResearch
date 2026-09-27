# 3-parameter Norris pulse fit: method note

**Status: closed for now, 2026-09-26 01:41 PKT -- user decision.** Not a negative result -- the method is
validated (see "Validation performed" below) but not wired into production, and several
production-relevant decisions (pulse-count calls for GRB231129C/GRB131014A, GRB080916C/GRB140206B
never run) are left open. Full closure note and state-to-resume-from: `PROGRESS.md`'s final entry.

Per-folder reference note (CLAUDE.md's "Per-folder method notes" convention). Written for the
user's own understanding of *how and why* this folder's code works, as a standalone reference --
not a chronological log (that's `PROGRESS.md`) and not a replacement for the spec itself
(`norris_3param_spec.md`, the source of truth for every formula below).

## What the code computes

Norris et al. (2005) pulse, standard 4-parameter form (amplitude $A$, start time $t_s$,
rise/decay timescales $\tau_1,\tau_2$):

$$I(t) = A\exp\left(2\sqrt{\tau_1/\tau_2} - \frac{\tau_1}{t-t_s} - \frac{t-t_s}{\tau_2}\right), \quad t > t_s$$

This project's existing pipeline (`norris_fit.py`) fits $(A, t_s, \tau_1, \tau_2)$ directly. This
folder implements a *reparametrization*, fitting $(A, t_\text{peak}, t_v)$ at a fixed asymmetry
ratio $r_0 = \tau_1/\tau_2$ instead, via the exact, closed-form linear map (`pulse3.py`,
`Pulse3Mapping`):

$$\tau_1 = \frac{r_0}{f(r_0)}\,t_v, \qquad \tau_2 = \frac{t_v}{f(r_0)}, \qquad t_s = t_\text{peak} - \frac{\sqrt{r_0}}{f(r_0)}\,t_v$$

where $f(r) = \frac{\sqrt{\ln 2}}{2}\sqrt{\ln 2 + 4\sqrt{r}}$ is the pulse's own width function
(FWHM/2 $= \tau_2 f(\tau_1/\tau_2)$), matching Bukhari et al. (2022, Adv. Space Res., eq. 10),
itself attributed there to Norris et al. (2005) -- validated bit-for-bit against the pipeline's
own `tv_value()` in `check_tv_consistency.py` (section 6, see "Validation performed" below).

Fixing $r_0$ removes the $\tau_1 \leftrightarrow \tau_2$ degenerate direction that trades the two
off at constant $(t_\text{peak}, t_v)$ -- the reduction's whole purpose, motivated by that
degeneracy showing up as unreliable free-$r$ fits on real pulses (see "Known limitations" for how
that manifests concretely).

## Every judgement call, and who made it

The governing spec (`norris_3param_spec.md`) was supplied by the user, fully specifying the model,
the bounds/seeding formulas, and the section-5 decision-rule *shape* (grid range, flat vs. sharp
distinction, cost-report requirement). Everything below is a concrete choice Claude made filling in
numbers or mechanism the spec left open, proposed and applied without a separate check-in unless
noted -- all open to being revisited by the user:

- **Flat-zone classification threshold: span $\geq 3\times$.** Spec says "spanning factors of
  several or more"; 3x was picked as the concrete cutoff (`FLAT_ZONE_MIN_SPAN`, `profile_scan.py`).
- **r0-systematic refinement window: $\pm$ a factor of 3 around $r_{0,\text{chosen}}$, 10
  log-spaced points, Delta-chi2 $\leq 1$.** Not in the spec at all -- added in response to external
  review (`review-on-progress.md` point 1): the coarse 25-point/5-decade section-5 grid
  under-resolves a sharp profile's own r0 uncertainty. Deliberately a *different* number from
  section 5's own $+10$ classification threshold; see `refine_r0_zone`'s docstring for why the two
  questions ("does this pulse constrain r at all" vs. "what is r0's own local uncertainty") need
  different tolerances.
- **4-param comparison fit: default 7 seeds, `geomspace` spanning within a decade of
  `DEFAULT_R_BOUNDS = (1e-3, 1e6)`.** Also review-driven (point 3) -- the original 5-seed default
  (`geomspace(1, 1000, 5)`) never reached near the bound, so a real box-edge runaway could go
  undetected by construction, not because none existed.
- **`sigma` mandatory, not optional, at every decision-driving entry point.** Review point 2:
  unweighted SSE must never silently reach the flat/sharp or Delta-chi2 verdicts, since both
  thresholds are only meaningful in chi^2 units.
- **Archived-median guard: warn + fall back to flat-zone center if the archived value falls
  outside this pulse's own measured flat zone.** Review point 4.
- **Neutral-seed jitter scale (`bounds_seeding.jittered_seeds`): relative to the seed's own value
  (A by $\pm30\%$, t_peak by $\pm5$ bins, t_v by $\pm50\%$), not a fraction of the full bound box.**
  Claude's own fix after the first version (jittering by bound-range fraction) scattered seeds into
  uninformative regions on a wide, mostly-empty fit window and produced a spurious
  reproducibility failure -- see `PROGRESS.md`'s section-8 entry.
- **RNG seeding: `seed_from_name(f"{__file__}::<dataset>")` per synthetic dataset, never a literal
  int.** Following this project's existing convention (`src/grb_research/SEEDING.md`); the initial
  implementation used literal seeds per test file and was corrected after the user flagged it --
  see `PROGRESS.md`'s "Seeding migration" entry.

## Conventions that could plausibly have gone another way

- **Default r0 grid: `logspace(-1, 4, 25)`** (0.1 to 10000), per the spec's own section 5 step 1.
  Not independently chosen, but worth flagging as a real, tunable input: pulses whose true $r$ or
  whose degenerate flat zone extends outside this range would need it widened (the spec itself
  says "extend if the profile is still falling at the edges" -- not yet automated, currently a
  manual judgement call per pulse).
- **`chi_square()` keeps an optional `sigma=None` -> unweighted-SSE fallback**, deliberately, even
  though every *decision-driving* caller requires `sigma` -- it remains available as a generic
  residual-reporting utility (e.g. `plot_diagnostics.py`'s fit overlays), just walled off from the
  verdict logic.
- **Bounds edge-pinning / box-edge-pinning tolerance: relative, not absolute** (`EDGE_REL_TOL`,
  `_is_pinned_r`'s `rtol`) -- appropriate given `r_bounds` and the fit window can each span many
  orders of magnitude, where an absolute tolerance would be meaningless at one end and useless at
  the other.
- **Window-widening invariance tolerance: 5% relative** (`WINDOW_WIDEN_REL_TOL`) -- arbitrary but
  generous; every synthetic test so far has passed to $<10^{-8}$ relative, so this has not yet been
  stress-tested against a real deviation.

## Validation performed

Full detail and every numeric result: `PROGRESS.md` (chronological, one entry per section/fix,
each timestamped). Summary:

- **Section 6** (`check_tv_consistency.py`): pipeline's `tv_value()` validated against the spec's
  width formula to machine precision, and against all 112 rows of archived 4-param fit results
  (within 1.31 sigma worst case, accounting for the archived column being an MC median, not a
  point estimate -- see PROGRESS.md for why that's expected, not a bug).
- **Sections 1-2, 3** (`pulse3.py`, `bounds_seeding.py`): round-trip, width-function, overflow-safety,
  and bounds/seeding tests -- `test_pulse3.py`, `test_bounds_seeding.py`, both all-pass.
- **Section 5** (`profile_scan.py`): validated on the spec's own T4-style synthetic scenarios (a
  well-resolved sharp case and a rise-truncated flat case) -- `test_profile_scan.py`.
- **Section 4** (`report_pulse.py`): uncertainty propagation validated on the same two scenarios,
  confirming the r0-choice systematic dominates on a flat profile and is small-but-nonzero on a
  sharp one (post-review fix) -- `test_report_pulse.py`.
- **Section 7 T1-T5**: T1/T2 in `test_pulse3.py`; T3 (FWHM/2 via independent `brentq` root-finding,
  not sharing any code path with `width_function()`) in `test_fwhm_brentq.py`; T4 (multi-start
  reproducibility of the 3-param fit, contrasted against the 4-param fit's degenerate behavior) in
  `test_multistart_t4.py`; T5 (profile-sanity Delta-chi2 check on the T4 datasets) in
  `test_profile_sanity_t5.py`.
- **Section 8** (`audit_checklist.py`): all 5 checklist items, including a new window-widening
  invariance diagnostic (not in any earlier section) -- `test_audit_checklist.py`.
- **External review** (`review-on-progress.md`, 2026-09-25): 5 points raised, 4 addressed in code
  (see "Every judgement call" above), 1 (comparison baseline) is applied guidance for the real-data
  run below. Full response: `PROGRESS.md`'s "External review response" entry.
- **GRB231129C real-data run** (`GRB231129C/`, 2026-09-25): the full pipeline (residual
  construction -> select_r0 -> finalize_pulse -> run_audit, two passes) run on a real 6-pulse GBM
  light curve for the first time. All 6 pulses converged and passed the full section-8 audit; new
  point-estimate (t_peak, t_v) agreed with the archived 4-param point estimate to <0.5% for 5/6
  pulses. Full detail, including the one pulse (6) with a larger (~3%) deviation traced to a
  neighbor-subtraction artifact, not a pipeline bug: `PROGRESS.md`'s "GRB231129C real-data run"
  entry.

All 9 synthetic-data test/check scripts pass as of the last full rerun (`PROGRESS.md`,
2026-09-25 01:59 PKT), verified deterministic (byte-identical output across separate process
invocations); the GRB231129C run (below) is the first real-data validation.

## Known limitations and open questions

- **Neighbor-pulse contamination in the shared-window residual method.** GRB231129C's pulse 6
  (the last of 6 in a tightly-packed episode pair) showed visible residual structure -- an
  undershot peak, a decay-shape mismatch, and a not-quite-flat pre-pulse baseline -- consistent
  with imperfect subtraction of the preceding pulse's decaying tail. The section-8 audit passed
  cleanly regardless: it checks internal consistency of the fit *given the residual it was
  handed*, not freedom from neighbor contamination or absolute goodness-of-fit -- a real scope
  boundary, not a bug. A genuine joint multi-pulse fit (spec section 5's "3n parameters" default)
  or a tighter per-pulse window are the two candidate fixes, neither implemented yet.
  (`PROGRESS.md`'s GRB231129C entry has the full writeup.)
- **Single-pulse-per-fit only; no joint multi-pulse fit implemented.** Every pulse (synthetic or
  real) has so far been fit one at a time, either alone or against a residual with other pulses
  subtracted -- never all pulses of a burst fit simultaneously. This is the spec's own stated
  default (section 5: "independent 3-param pulses (3n parameters)"), but it's also the most
  direct fix for the neighbor-contamination limitation above.
- **The T4b degeneracy symptom is noise-draw-dependent.** The spec names "large r scatter
  (possibly a runaway pinned to a box edge)" as the expected symptom of a degenerate 4-param fit.
  Two different (both properly-seeded) noise realizations of the same synthetic scenario produced
  two different concrete symptoms -- one a genuine box-edge runaway, one a confident-but-wrong
  interior attractor -- both are the same underlying pathology (an unreliable free-r estimate that
  *looks* numerically stable), but neither is guaranteed for a given pulse. Don't expect a specific
  symptom on real data; expect the pathology, and rely on the profile-scan/audit machinery (not
  visual inspection of one multi-start run) to catch it.
- **Not wired into production.** This is a standalone validation folder; `norris_fit.py` (the
  pipeline every `fitter_*.py` script actually uses) is untouched. Promoting the 3-param reduction
  into production is a separate, not-yet-scoped decision.
- **Hierarchical (one r shared per episode) variant not implemented.** Spec section 5's last
  paragraph mentions this as optional, contingent on archived per-pulse r values clustering by
  episode. A whole-project population check ran during the GRB231129C run (`load_project_r_population`,
  83 archived r values across all 4 GRBs): median 3.31, spanning ~9 orders of magnitude (6e-5 to
  1.9e4) end to end -- clearly not clustered *project-wide*. That's not yet the real test the spec
  asks for (well-determined, non-degenerate r values *within one episode*, filtered by fit
  quality) -- still open.
- **Default r0 grid range is not adaptive at the library level.** `profile_scan.select_r0` itself
  still always uses the fixed `logspace(-1, 4, 25)` grid unless the caller overrides it. The
  GRB231129C run added an adaptive retry *in that folder's own driver script*
  (`GRB231129C/run_pipeline.py`'s `select_r0_with_edge_check`, widening by 2 decades and retrying
  once if the flat zone touches the grid edge -- fired for one of that burst's 6 pulses), not in
  `profile_scan.py` itself; promoting that retry logic into the shared library is still open.
- ~~`audit_checklist.py`'s `_is_pinned` has a latent span-relative-tolerance bug.~~ **Fixed
  2026-09-25.** Found first in `GRB231129C_joint/run_joint_fit.py`'s own `check_edge_pinning`
  (during the full-x-span window-sensitivity test): a tolerance relative to the *bound span*
  rather than the *bound value* falsely flags ordinary values as "pinned" once the span becomes
  large relative to the fitted value (e.g. a very wide fit window). Ported the same value-relative
  fix (matching `profile_scan.py`'s `_is_pinned_r`) to `audit_checklist.py`'s single-pulse version;
  full synthetic suite + the real GRB231129C single-pulse run re-verified with no regressions, and
  `test_hidden_peak.py`'s genuine pinning case (t_v at the real 2*dt floor) still correctly flags.
  `_is_pinned`/`_is_pinned_r` now share the same value-relative form everywhere in this folder.
