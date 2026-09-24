# 3-param Norris fit -- implementation progress

Tracks work against `norris_3param_spec.md`, section by section. Each entry: what was built, how
it was tested, and the result. Update this file whenever a section's code or tests change. This is
the chronological log; `norris_3param_reduction.md` (this folder) is the standalone, non-chronological
method-note reference (CLAUDE.md's per-folder method note convention) -- defining equations,
judgement-call attribution, conventions, and known limitations, pointing back here for detail.

Convention: every dated entry uses full local timestamp `YYYY-MM-DD HH:MM` (not date-only), taken
from `date "+%Y-%m-%d %H:%M %Z"` at the time the test was run.

## Status at a glance

| Spec section | What | File | Status |
|---|---|---|---|
| 6 | tv_value() consistency check | `check_tv_consistency.py` | done, all pass |
| 1, 2 | pulse3 model (f(r), linear map, safe eval) | `pulse3.py` | done, all pass |
| 3 | bounds + neutral seeding | `bounds_seeding.py` | done, all pass |
| 4 | forward conversion + uncertainty propagation | `report_pulse.py` | done, all pass |
| 5 | chi^2 profile scan over r0 | `profile_scan.py`, `plot_profile_scan.py` | done, all pass |
| 7 | T3/T4/T5 (FWHM brentq, synthetic degeneracy, profile sanity) | `test_fwhm_brentq.py`, `test_multistart_t4.py`, `test_profile_sanity_t5.py` | done, all pass (T1/T2 in `test_pulse3.py`) |
| 8 | post-fit audit checklist | `audit_checklist.py` | done, all pass |

## Section 6 -- tv_value() consistency check

**File:** `check_tv_consistency.py`

Validated the pipeline's `tv_value()` (`../../norris_fit.py`) against the spec's width formula
before writing anything downstream, per the spec's explicit "READ BEFORE IMPLEMENTING" instruction.

**Results (2026-09-25 01:17 PKT):**
- Worked example (tau1=85.6, tau2=0.453): both `tv_value()` and an independently-transcribed spec
  formula give `1.4071003588452418`, matching the spec's `1.40710035884524` to machine precision.
- 2000 random (tau1, tau2) pairs in [1e-3, 200]: 0.0 abs/rel diff between `tv_value()` and the spec
  formula -- confirmed to be literally the same expression.
- All 112 rows across the 14 archived `norris_fit_results*.csv` files: `tv_value()` at the archived
  point-estimate (tau1, tau2) lands within 1.31 sigma of the archived `t_v_s` (worst case), well
  inside the 5-sigma tolerance.

**Finding:** archived `t_v_s` is *not* `tv_value(tau1, tau2)` -- it's the MC-propagated median over
the fit covariance (`tv_mc_summary()` in `norris_fit.py`), which differs from the point estimate by
O(1 sigma) on skewed/degenerate pulses. Largest deviations coincided with the lowest
`mc_kept_fraction` rows (0.48-0.90) -- exactly the degenerate pulses the 3-param reduction targets.
Not a bug; documented in the script.

**Conclusion:** f(r) in the spec matches the pipeline exactly. No rebuild needed.

## Sections 1-2 -- pulse3 model

**File:** `pulse3.py` (self-contained, no import of `norris_fit.py`)

Implements: `width_function(r)`, `width_function_inverse(f)`, `norris_raw(t, A, t_s, tau1, tau2)`
(numerically-safe combined-exponent form, clipped at arg<=50), `Pulse3Mapping(r0)` (precomputes
c_tau1, c_tau2, c_ts), `make_pulse3(r0)` (curve_fit-ready closure with `.mapping` exposed),
`pulse3(t, A, t_peak, t_v, r0)` (one-off convenience form), `reduced_from_raw(t_s, tau1, tau2)`
(inverse direction), `pulse4(t, A, t_peak, t_v, r)` (optional 4-param comparison wrapper).

**Tests:** `test_pulse3.py` -- covers spec section 7 T1 and T2, plus overflow-safety and
API-consistency checks.

**Results (2026-09-25 01:17 PKT): ALL PASS**
- T1 round-trip (tau1=85.6, tau2=0.453, t_s=3.0): recovered to machine precision
  (tau1=85.60000000000001, tau2=0.45300000000000007, t_s=3.0).
- T2 width_function: f(0)=0.3465735902799726 (== ln2/2), monotone increasing on
  logspace(-6,8), f(1e8)=83.2561824662465 (~83.26 as spec requires), r(f) inverse round-trips.
- Overflow safety: r0=2e6 -> mu~1414.2 (2*mu alone would overflow float64, max exp arg ~709);
  `norris_raw` stays finite and non-negative across the full evaluation grid.
- API consistency: `pulse3(...)` == `make_pulse3(r0)(...)`; `pulse3` matches `pulse4` at r=r0;
  `make_pulse3(r0).mapping.r0 == r0`.

**Not yet covered here:** T3 (FWHM/2 via brentq) -- deferred to section 5/7 work since it needs a
root-finder harness, not just the model definition.

## Section 3 -- bounds and seeding

**File:** `bounds_seeding.py` (self-contained, matches `pulse3.py`'s style)

Implements: `median_dt(t)`, `pulse3_bounds(t_window, dt)` -> `(lb, ub)` in curve_fit's `(lb, ub)`
3-tuple form, `pulse3_seed(t_window, y_window, dt, background_rate=0)` -> `(A0, t_peak0, t_v0)`
(t_v0 = 2.5*dt, clamped into `[2*dt, window]`), `pulse3_bounds_and_seed(...)` convenience wrapper
that also asserts the seed falls inside the bounds it returns.

Design notes:
- This pipeline's light curves are already background-subtracted count *rates* (per
  `NorrisFitter`'s docstring in `norris_fit.py`), so `y_window` is used directly as a rate and
  `background_rate` defaults to 0 -- the spec's more general `counts_per_bin/dt - background_rate`
  form is available via the `background_rate` argument for count-valued data.
- `t_v0` is explicitly clamped to `[2*dt, window]` so the seed is always bounds-valid even in an
  unusually narrow fit window (spec only says "~2-3 dt"; clamping was an addition beyond the literal
  spec text, needed for `pulse3_bounds_and_seed`'s in-bounds assertion to always hold).

**Tests:** `test_bounds_seeding.py`

**Results (2026-09-25 01:17 PKT): ALL PASS**
- Bounds formula check: window=[0,50], dt=0.5 -> lb=(0.0, 0.0, 1.0), ub=(inf, 50.0, 50.0).
- Narrow-window guard: window <= 2*dt raises `ValueError` (would otherwise produce an empty t_v
  bound).
- `median_dt` recovers the true bin width on an evenly-spaced grid.
- Neutral seed on a clean synthetic pulse (A=12, t_peak=20, t_v=1.5, r0=190): t_peak0 within 1 dt
  of truth, A0 within 1% of truth, t_v0 respects bounds.
- Clamping check: narrow window (dt=1, window=3) still yields an in-bounds seed (t_v0=2.5, bounds
  [2,3]).
- End-to-end smoke test: bounds+seed fed straight into `scipy.optimize.curve_fit` on a noisy
  synthetic pulse (5% Gaussian noise on A) converges to A within 5%, t_peak within 0.1, t_v within
  10% of the injected truth.

## Section 5 -- chi^2 profile scan over r0

**Files:** `profile_scan.py` (algorithm, self-contained -- no matplotlib/grb_research import),
`plot_profile_scan.py` (reporting: CSV + PDF/PNG, split out so the algorithm stays dependency-light)

Implements all 4 steps of the spec's decision rule:
- `chi_square(y_data, y_model, sigma=None)` -- weighted if sigma given, else plain SSE.
- `fit_single_r0(...)` / `profile_scan(...)` -- step 1-2: curve_fit at each r0 in
  `DEFAULT_R0_GRID = logspace(-1, 4, 25)` from the section-3 neutral seed, recording chi2/A/t_peak/t_v.
- `classify_profile(...)` -- step 3: FLAT (chi2 within `+10` of its min over a contiguous r0 run
  spanning `>=3x`, both thresholds named constants) vs. SHARP; returns the flat zone and a
  recommended r0 (flat-zone geometric center, or the minimum for sharp).
- `fit_4param_multistart(...)` -- multi-seed (5 seeds by default) free-r `pulse4` fit for comparison.
- `delta_chi2_report(...)` -- step 4: chi2_3param(best r0) - chi2_4param(best multi-seed); verdict
  "dead weight" if `<=10`, "keep r free" if `>10`.
- `select_r0(...)` -- orchestrates all of the above into one call; on a flat-bottom profile, fixes
  r0 to `archived_r_median` if supplied (spec's preferred action) else the flat-zone geometric center;
  on sharp, fixes r0 at the profile minimum. Returns an `R0SelectionResult` with the full audit trail.

`plot_profile_scan.py` adds, per this project's generated-output convention (CSV via
`df.to_csv`, figures via `grb_research.update_style()` + `grb_research.grb_utils.save_fig`, both
written into this folder -- matching `experiments/window_sensitivity_GRB140206275/`):
- `save_profile_scan_outputs(result, out_dir, label)` writes `r0_profile_scan_<label>.csv` (full
  grid), `r0_profile_scan_summary_<label>.csv` (one-row decision outcome), and
  `r0_profile_scan_<label>.pdf`/`.png` (chi2 vs r0, log x-axis, flat-zone shaded, chosen r0 marked,
  best 4-param chi2 overlaid) -- the plot the spec's step 3 asks to make, saved rather than shown.

**Tests:** `test_profile_scan.py` -- the two synthetic scenarios the spec itself describes (T4/T5):
a well-resolved pulse (sharp minimum expected) and a pulse whose window starts partway up the rise
(spec's T4 setup; flat bottom expected). Also exercises `save_profile_scan_outputs` so the run
leaves real files in this folder.

**Results (2026-09-25 01:22 PKT): ALL PASS**
- Sharp scenario (r_true=190, full window, 2% noise): classified `sharp`, r0_at_min=215 (within a
  factor of 2 of truth), flat zone span x1.0 (single grid point), chi2_3param=971 vs. best
  4-param chi2=971 (r_fit=219) -> Delta-chi2=0.014 -> "4th parameter is dead weight" (expected: at
  the correct r0, the 3-param and 4-param fits should coincide). Wrote
  `r0_profile_scan_synthetic_sharp.{csv,pdf,png}` + summary CSV.
- Flat scenario (same r_true=190, window starts at t=14.5 just before t_peak=15, cutting off most
  of the rise; 8% noise): classified `flat`, flat zone spans the *entire* grid ((12.1, 10000),
  x825) -- the missing rise leaves essentially no asymmetry leverage, exactly the spec's
  degenerate case. r0_chosen correctly took the supplied `archived_r_median=200` rather than the
  flat-zone center. Delta-chi2=0.015 -> also "dead weight" here, consistent with the flat profile
  (any r0 fits about as well). Wrote `r0_profile_scan_synthetic_flat.{csv,pdf,png}` + summary CSV.

**Design choices beyond the literal spec text** (documented since they weren't fully specified):
- Flat-zone span threshold for flat vs. sharp: 3x (spec only says "factors of several or more").
- `chi_square()` sigma-weighting: unweighted SSE if no `sigma` passed, chi^2-weighted if it is --
  matches this pipeline's real per-bin uncertainties being available (`light_curves.py`'s
  `data_rate_error`, quadrature-summed across NaI channels) but not required for a scan run without them.
- 4-param multi-start seeds: 5 seeds geometrically spaced in r over `[1, 1000]` by default,
  matching section 7 T4's "5 neutral seeds" reproducibility framing.

## Section 4 -- forward conversion + uncertainty propagation

**File:** `report_pulse.py` (builds on `pulse3.py`'s `Pulse3Mapping` and `profile_scan.py`'s
`R0SelectionResult`/`fit_single_r0`)

Built *after* section 5 rather than in spec order, deliberately: the r0-choice systematic this
section needs (spread of best-fit t_peak/t_v across the flat zone) only exists once the section-5
profile scan has been run, so completing section 4 was blocked until section 5 existed. See "build
order" note at the end of this file.

Implements `finalize_pulse(t_window, y_window, selection, ...)`: refits at
`selection.r0_chosen` directly if that r0 isn't exactly a profile-scan grid point (so the
reported values are always for the r0 actually used downstream), forward-converts
(A, t_peak, t_v) -> raw (A, t_s, tau1, tau2) via `Pulse3Mapping.to_raw`, and reports uncertainty
as statistical (fit covariance) and systematic (r0-choice, i.e. half the flat-zone spread of
t_peak(r0)/t_v(r0)) combined in quadrature -- exactly the spec's "add in quadrature with the
statistical term" instruction. `tau1_err`/`tau2_err` propagate from `t_v_err` through the (exact,
linear) map; `t_s_err` additionally uses the fit's t_peak-t_v covariance for its statistical part.
Returns a `ReportedPulse` dataclass carrying every raw/reduced value and its error, plus
`r0_flat_zone` (the spec's prescribed "uncertainty" on the fixed r).

**Tests:** `test_report_pulse.py` -- reuses the sharp/flat synthetic scenarios from
`test_profile_scan.py`, since section 4's systematic term should behave oppositely in each.

**Results (2026-09-25 01:25 PKT): ALL PASS**
- Sharp case: flat zone is a single grid point, so `t_peak_err_sys = t_v_err_sys = 0.0` exactly
  and the combined error equals the statistical term alone -- A, t_peak, t_v all recovered within
  3-sigma of the injected truth (A=20.0021+/-0.0631 vs 20.0, t_peak=14.9999+/-0.0037 vs 15.0,
  t_v=1.2047+/-0.0044 vs 1.2).
- Flat case: t_v's systematic term (0.1398) dominates its statistical term (0.0347) by ~4x --
  the flat (whole-grid) profile correctly demotes the naive covariance-only precision, which would
  otherwise understate the real r0-choice uncertainty by about that factor. t_v=1.1861+/-0.1440
  vs true 1.2 (now correctly inside the combined error bar).

## Section 7 -- T3, T4, T5 (T1/T2 were already covered in `test_pulse3.py`)

**Files:** `test_fwhm_brentq.py` (T3), `bounds_seeding.jittered_seeds` (new, section-3 file) +
`profile_scan.fit_4param_all_seeds` (new, section-5 file) + `test_multistart_t4.py` (T4),
`test_profile_sanity_t5.py` (T5)

**T3 -- FWHM/2 via brentq.** `half_max_width()` locates the two half-max points of `norris_raw`
by root-finding (`scipy.optimize.brentq`), independent of `width_function()` entirely (can't share
a transcription bug with it). Checked against the spec's worked example, against
`width_function()` across 5 different r spanning 0.01-1e4, and against `reduced_from_raw()`.

**Results (2026-09-25 01:29 PKT): ALL PASS**
- Worked example: brentq FWHM/2 = 1.4071003588452453 vs. spec's 1.40710035884524.
- All 5 (tau1, tau2) pairs (r = 0.01, 1, 50, 189, 1e4): brentq matches `width_function()`-derived
  t_v to rtol 1e-6.

**T4 -- multi-start reproducibility / degeneracy contrast.** Added `jittered_seeds()` to
`bounds_seeding.py` (neutral seed + 4 perturbed variants, clamped into bounds) and
`fit_4param_all_seeds()` to `profile_scan.py` (like `fit_4param_multistart` but returns every
seed's result, not just the best) specifically to support this test.

**Results (2026-09-25 01:29 PKT): ALL PASS, with one finding worth flagging:**
- T4a (3-param reproducibility): on a moderately rise-truncated synthetic pulse, 5 jittered
  seeds at 5 different r0 across the grid all converge to the same solution (chi2 agrees to 6
  decimal places; A/t_v agree to 1e-3 relative, t_peak to 1e-3 s -- limited only by optimizer
  floating-point precision along the shallow direction, not by seed choice).
- T4b (4-param degeneracy): the spec's literal expectation is "large r scatter (possibly a
  runaway pinned to a box edge)". On the *same* moderate truncation used for T4a, the free-r fit
  instead converged reproducibly to a single wrong interior value (r=11.06, seed-independent) --
  not scatter, and not a box-edge runaway. A *more* severely truncated window (starting almost at
  the peak, cutting off essentially all rise and most of the decay) did reproduce the box-edge
  case exactly: every one of 9 widely-spaced seeds ran away to the upper r-bound (1e6, vs.
  r_true=190, a >5000x miss) while t_peak stayed stable to 3e-5 s. **Finding:** for a smooth,
  single-valley chi2(r) landscape, a bounded local optimizer (scipy's trust-region-reflective)
  converges to the *same* attractor from any start -- "scatter" isn't guaranteed even when the
  fit is badly degenerate; what's guaranteed is that the free-r estimate can be confidently
  wrong (interior false-minimum or box-edge runaway) while looking numerically stable, which is
  arguably the more dangerous failure mode for downstream code that doesn't itself check for it.

**T5 -- profile sanity.** Confirmed on both T4 datasets: `select_r0()` classifies both as `flat`,
and Delta-chi2 (3-param best vs. 4-param best) is 0.006 and 0.510 respectively -- both comfortably
inside the spec's "within ~1" tolerance.

## Section 8 -- post-fit audit checklist

**File:** `audit_checklist.py` (builds on `bounds_seeding.py`, `profile_scan.py`, `report_pulse.py`)

Implements all 5 checklist items as independent `AuditItem` checks, assembled by `run_audit()`
into an `AuditReport`:
- `check_edge_pinning` -- t_peak/t_v against the section-3 bounds.
- `check_multistart_reproducibility` -- reuses `jittered_seeds` (see below) + `fit_single_r0`.
- `check_delta_chi2_recorded` -- reads `selection.delta_chi2` (section 5.4).
- `check_flat_zone_recorded` -- confirms `reported.r0_flat_zone` is valid and the section-4
  quadrature combination was actually applied (combined error >= statistical alone).
- `check_window_widening` -- new diagnostic, not built anywhere else: takes a
  `make_window_data(t_min, t_max)` callable (data-source-agnostic -- regenerated synthetic data in
  tests, a slice of the same real light curve in production) and confirms t_peak/t_v stay within
  `WINDOW_WIDEN_REL_TOL` (5%) as the window widens, with the result also stashed on the
  `AuditItem.data` field so callers can plot it without duplicating the fit loop.

**Bug found and fixed while testing this:** `jittered_seeds()` (added in the T4 work, section 3
file) originally jittered each parameter by a fraction of its *entire bound range*. That's fine
for a tight fit window, but for a wide, mostly-empty window (this section's own audit test: a
30s baseline around a 1.2s pulse) it dragged t_peak seeds ~9s away and t_v seeds into the noise
floor, landing genuinely different local optima and failing the reproducibility check --
correctly, in the sense that the *test* was wrong, not the 3-param fit. Fixed by jittering
relative to the neutral seed's *own* scale instead (A by +/-30%, t_peak by +/-5 bins, t_v by
+/-50%), which is what "a slightly different but still reasonable guess at this pulse" should
mean. Re-verified T4a still passes with the new scaling.

**Results (2026-09-25 01:37 PKT): ALL PASS**
- On a clean, well-resolved synthetic pulse (30s window, full pulse + baseline, 3% noise):
  no parameter pinned to a box edge; 5 jittered seeds agree to ~1e-8 (A, t_v relative) / ~1e-8 s
  (t_peak); Delta-chi2 recorded (=62.05, correctly flagged "keep r free" here since this pulse
  *is* well-resolved, unlike the T4/T5 degenerate cases); flat-zone recorded (a single grid point
  here, sys=0, matching the sharp-profile case from section 4); window-widening invariant to
  ~1e-9 s under x1.5/x2.0 widening (baseline padding, same noise realization).

## Plotting: retrofitted CSV/PDF/PNG outputs to sections 4, 7, 8

Section 5 (`plot_profile_scan.py`) already wrote CSV/PDF/PNG outputs per this project's
convention; sections 4/7/8 initially didn't. Added `plot_diagnostics.py` (same
matplotlib-Agg + `grb_research.update_style()`/`save_fig()` pattern) and wired it into every test:
`plot_fit_overlay` (section 4's finalized-pulse overlay; section 7 T3's half-max-point overlay),
`plot_t4_scatter` (T4's seed-r -> fitted-r and t_peak-stability scatter), `plot_window_widening`
(section 8's new diagnostic). `test_profile_sanity_t5.py` reuses section 5's own
`save_profile_scan_outputs` rather than a new function, since T5 literally is the profile scan
applied to the T4 datasets.

**Axis-scaling bug found and fixed:** the first versions of the T4 t_peak-stability panel and the
window-widening panel autoscaled y-axis limits to the data's actual (tiny) range. When a result is
genuinely invariant to ~1e-5-1e-9 s -- the *correct* outcome here -- matplotlib's autoscale
zooms into that noise floor and its offset notation (e.g. axis labeled "1e-5+1.49e1") reads, at a
glance, as a real trend across a huge swing, when the swing is actually floating-point-level
noise. Fixed by pinning y-limits to a fixed, physically meaningful tolerance band around the mean
(+/-0.05 s / +/-1% for t_peak-stability, +/-5% for window-widening, matching each check's own pass
threshold) instead of autoscaling, with `ticklabel_format(useOffset=False)` so labels always read
as plain values. A genuinely flat result now draws as a flat line inside a visible band; only a
real, tolerance-breaking deviation would show up as visible movement.

## Spec coverage: complete

All 8 sections of `norris_3param_spec.md` now have working, tested code in this folder. Remaining
optional items from the spec text itself (not required by any section's checklist): the
hierarchical multi-pulse r-sharing variant mentioned in section 5's last paragraph, and wiring this
3-param fitter into the production pipeline (`norris_fit.py`) in place of / alongside the existing
4-param `NorrisFitter` -- both are new scope beyond what the spec asked to be validated first.

## External review response (2026-09-25 01:59 PKT)

`review-on-progress.md` (this folder) reviewed the section 1-8 work above and flagged 5 issues
before the real-GRB run. All 4 code issues fixed; issue 5 (comparison baseline for the real run)
is guidance applied when that run happens, not a code change here.

**1. Sharp-profile `sys=0` was a grid-resolution artifact, not a real zero.** The coarse
section-5 grid (25 points, 5 decades, flat_tol=10) collapses to a single point on a sharp profile,
so `report_pulse.py` was reading an exact 0.0 systematic there -- correct-looking, not
correct: a sharp minimum still has real curvature-based r0 uncertainty, just narrower than the
coarse grid can resolve. Fixed: `profile_scan.refine_r0_zone()` (new) always runs a dedicated,
finer local scan (10 points, log-spaced over r0_center/3 to r0_center*3) with its own
Delta-chi2<=1 zone -- a real ~1-sigma-flavored interval, deliberately different from section 5's
own +10 classification threshold (which answers "does this pulse constrain r at all", not "what
is r0's local uncertainty"). `report_pulse.finalize_pulse` now uses this refined zone for the
t_peak/t_v systematic on *every* pulse, flat or sharp. Verified: the sharp scenario's sys term
went from exactly (0.0, 0.0) to (0.0043, 0.0001) -- small, genuinely non-zero, and still
correctly dominated by the statistical term (test_report_pulse.py's ok3/ok4 updated to check
"small but non-zero" and "combined >= stat" rather than exact equality to zero).

**2. Unweighted SSE must never reach the decision logic.** `chi_square()`'s `sigma=None` ->
plain-SSE fallback is fine as a generic utility, but the flat/sharp +10 threshold and the
Delta-chi2=10 verdict are only meaningful in chi^2 units; on raw SSE both scale with the
arbitrary magnitude of the rate, silently invalidating every verdict. Fixed: `_require_sigma()`
(new, `profile_scan.py`) raises `ValueError` at `fit_single_r0`, `profile_scan`, and
`fit_4param_all_seeds` if `sigma is None` -- every decision-driving entry point now refuses to run
unweighted rather than silently producing meaningless thresholds. Every existing test already
passed `sigma` explicitly, so this was a pure tightening with zero call-site breakage (confirmed
by the full rerun below).

**3. The 4-param comparison fit's own edges weren't audited.** If the best-of-seeds 4-param
comparison fit itself ran away and pinned to an r-bound, `delta_chi2_report` was still emitting a
normal "dead weight"/"keep r free" verdict -- comparing the 3-param fit against a comparison that
is itself untrustworthy. Fixed: `delta_chi2_report` now checks `fit4["r"]` against `r_bounds`
(`_is_pinned_r`, new) and emits an explicit `"inconclusive: 4-param comparison fit ran away and
pinned at r=..."` verdict plus a `r_pinned` key instead, so a pinned comparison is flagged rather
than silently trusted either way. `select_r0` and `fit_4param_multistart`/`fit_4param_all_seeds`
now share one `DEFAULT_R_BOUNDS = (1e-3, 1e6)` constant, threaded through consistently.

Found *while* implementing this fix: the old default 4-param seed range
(`geomspace(1, 1000, 5)`) never reached anywhere near `r_bounds`'s 1e6 edge, so the pinned-check
could never fire in the first place regardless of whether a pinned solution existed -- fixed
alongside, default seeds now span `geomspace(max(lo*10, 1e-2), min(hi/10, 1e5), 7)`, reaching
within a decade of both bounds by default.

**4. The archived-median path took `archived_r_median` on faith.** If it falls outside this
pulse's own measured flat zone (possible, since it comes from potentially-degenerate archived
4-param fits), `select_r0` was fixing r0 there anyway -- silently worse than the flat-zone center
by this pulse's own profile. Fixed: `select_r0` now checks `zone_lo <= archived_r_median <=
zone_hi`; outside the zone, it `warnings.warn()`s and falls back to the flat-zone geometric
center instead. Verified with a deliberately-bad median (0.01, far outside a (0.1, 10000) flat
zone): warning fired, r0_chosen correctly fell back to the geometric center (31.6) instead of the
bad value.

**5. Real-GRB comparison baseline (guidance for the upcoming run, not a code change here).**
Section 6 already found archived `t_v_s` is the MC-propagated median, not
`tv_value(tau1, tau2)` at the point estimate, differing by O(1 sigma) on exactly the degenerate
pulses this project cares about. When comparing new 3-param point-estimate fits against archived
values for GRB231129C, compare against archived *point* `tv_value(tau1, tau2)` (validated
available and correct in section 6), not the archived `t_v_s` MC-median column, or the comparison
will manufacture disagreements on the fragile pulses that aren't real.

**Sequencing note applied:** for the real multi-pulse run, the profile scan for pulse *i* will run
on the residual after subtracting the other pulses' current-best model (not the raw summed light
curve), iterating once if the joint fit moves things -- this is what section 3's own
"residual light curve" language already implies for multi-pulse seeding, now explicit for the
real-data case. Window-widening (section 8) is the designated detector if truncation-driven
degeneracy (T4b's finding) resurfaces through window choice on real data.

**Full rerun after all 4 fixes: all 9 scripts still ALL PASS** (`test_report_pulse.py`'s two
assertions that depended on the old exact-zero sys behavior were updated to match the corrected,
more honest behavior -- see point 1 above; everything else needed no changes).

## Seeding migration (2026-09-25 01:59 PKT)

Every noisy synthetic dataset in this folder previously used a hardcoded literal RNG seed
(`np.random.default_rng(1)`, `(7)`, etc.) instead of this project's convention
(`src/grb_research/SEEDING.md`): `SEED = seed_from_name(__file__)` + `get_rng(seed=SEED)`, never a
bare literal. Worse, `test_multistart_t4.py` (T4) and `test_profile_sanity_t5.py` (T5, "the chi2
profile ... on the T4 simulation") independently retyped the same literal seeds (3, 7) in two
files to get "the same dataset" -- exactly the fragile drift-prone pattern the project's real
seeding overhaul (2026-09-03) eliminated everywhere else.

Fixed: new `synthetic_datasets.py` centralizes every shared scenario (`sharp_scenario`,
`flat_scenario`, `moderate_truncation_scenario`, `severe_truncation_scenario`, `audit_source`),
each seeded via `seed_from_name(f"{__file__}::<dataset-name>")` -- a name-derived seed, not a
literal, so any caller in any process gets a bit-identical dataset by construction. Every test file
that duplicated one of these constructions (`test_profile_scan.py`, `test_report_pulse.py`,
`test_multistart_t4.py`, `test_profile_sanity_t5.py`, `test_audit_checklist.py`) now imports the
shared function instead. Files with their own one-off RNG use (`check_tv_consistency.py`,
`test_bounds_seeding.py`) now use `seed_from_name(__file__)` directly. `bounds_seeding.jittered_seeds`,
`audit_checklist.check_multistart_reproducibility`, and `audit_checklist.run_audit` had their
`rng` parameter's silent `None` -> unseeded-`default_rng()` fallback removed and made a required
argument, matching this project's stated tightening for library helpers ("`rng` made required,
since every real caller already passed rng=" -- SEEDING.md).

**Consequence worth flagging:** migrating to name-derived seeds changed the *specific* noise
realizations in `test_multistart_t4.py`'s T4b scenario -- the old literal-seeded draw happened to
show the free-r fit running away to a box edge (r=1e6) from every seed; the new, properly-seeded
draw instead converges reproducibly to a different, still-wrong interior value (r~32.4 vs true
190). Both are genuine instances of the same underlying pathology the spec's 3-param reduction
exists to guard against (a confident, seed-independent, badly wrong free-r estimate) -- just not
the literal "pinned to a box edge" phrasing from the spec text for this particular draw. Test
assertions were generalized to check for that pathology (reproducible convergence to a value far
from truth) rather than the specific box-edge symptom; T4b's docstring/prints now say so
explicitly. All numeric values below and in earlier sections of this file reflect the
seed_from_name-derived datasets, not the original literal-seeded ones.

**Verified deterministic:** `test_multistart_t4.py` run twice as separate process invocations,
output diffed byte-for-byte identical.

## GRB231129C real-data run (2026-09-25 02:07 PKT)

**Files:** `GRB231129C/load_data.py`, `GRB231129C/run_pipeline.py` (new subfolder, per user request)

The first real multi-pulse light curve through the full section 1-8 pipeline: 6 pulses, real
GBM data (10-400 keV, NaI n3+n6+n7 summed, background-subtracted, dt=0.064s), real per-bin
uncertainties (quadrature-combined across detectors -- required, since `sigma` is now mandatory
per the review fix). Data/episode bounds/archived 6-pulse P0 copied from
`../../../GRB231129C/_common.py` (CLAUDE.md's copy-not-cross-import convention: that folder is a
sibling under `variability_analysis/`, not an ancestor already on this folder's own import chain).

**Method:** per-pulse residual construction, exactly the "subtract other pulses, scan the target
on what remains" the review's sequencing note called for, run twice (`fit_one_pass`, pass 1 using
the archived 4-param model for the other 5 pulses, pass 2 rebuilding residuals from pass 1's own
3-param results instead, to check how much iterating once moves things). Shared fit window
(-2, 10) s across all 6 pulses (covers every archived t_peak, 0.65-5.5s, with margin).
`archived_r_median` = project-wide population median of archived 4-param r = tau1/tau2 (83 values
across all 4 GRBs' archived CSVs, deduplicated per pulse; `load_project_r_population()`) --
matching the spec's literal "population median ... for this project", not just this burst's own 6
pulses. New: `select_r0_with_edge_check` retries with a 2-decades-wider r0 grid if the flat zone
touches the default grid's edge (spec section 5: "extend if the profile is still falling at the
edges") -- fired for pulse 2 (archived r=18442, an extreme-asymmetry pulse), which needed the grid
extended to `logspace(-1, 6, 35)` before its flat zone stopped touching the boundary.

**Results: all 6 pulses converged, all 6 fully passed the section-8 audit, both passes.**
Pass1->pass2 movement was millisecond-level for every pulse (largest: pulse 1's t_v moved 5.5ms) --
confirms one iteration is enough to reach self-consistency here, no further iteration needed.
`archived_r_median` (3.31) fell *outside* the pulse's own flat zone for 4 of 6 pulses (1, 2, 3, 5)
-- the guard fired correctly each time (warned, fell back to the flat-zone geometric center)
rather than trusting a population value this pulse's own data said was measurably worse. Only
pulse 4's flat zone happened to contain the population median, so it was the only pulse where the
archived value was actually used as-is.

New point-estimate (t_peak, t_v) vs. the archived **point** estimate (`tv_value`/`t_peak` at the
archived (tau1, tau2), NOT the archived `t_v_s` MC-median column -- review point 5):

| pulse | episode(s) | shape | r0_chosen | t_peak (new) | t_peak (archived) | delta t_peak | t_v (new) | t_v (archived, point) | delta t_v |
|---|---|---|---|---|---|---|---|---|---|
| 1 | EX0+TR1 | sharp | 1.10 | 0.6334 | 0.6516 | -0.0182 | 0.8147 | 0.8099 | +0.0048 |
| 2 | EX0+TR1 | flat | 5.82e4 | 2.4666 | 2.4657 | +0.0009 | 0.5431 | 0.5428 | +0.0003 |
| 3 | EX0+TR1 | sharp | 4.64 | 1.3162 | 1.3197 | -0.0034 | 0.6186 | 0.6237 | -0.0051 |
| 4 | TR2+EX1 | flat | 3.31 | 3.5034 | 3.4903 | +0.0131 | 0.4492 | 0.4511 | -0.0019 |
| 5 | TR2+EX1 | sharp | 1.78 | 4.2962 | 4.2874 | +0.0088 | 0.8284 | 0.8269 | +0.0015 |
| 6 | TR2+EX1 | flat | 0.332 | 5.5264 | 5.5075 | +0.0188 | 1.3173 | 1.3576 | -0.0403 |

Delta-chi2 (3-param vs. best 4-param) is small for every pulse (0.0005-0.78, all "4th parameter is
dead weight") -- the reduction costs essentially nothing statistically anywhere in this burst.

**Honest finding: pulse 6's fit has visible residual structure the audit didn't catch.**
`fit_overlay_pulse6_pass2.png` shows the model undershooting the data's peak and mismatching the
decay shape, sitting on a residual baseline that isn't quite flat before the pulse starts (t=0-4.5s)
-- consistent with imperfect subtraction of pulse 5's decaying tail (pulse 5 sits right before
pulse 6 in the same TR2+EX1 episode pair). This is pulse 6's own largest deviation from the
archived point estimate (t_v: -0.040s, ~3% relative, vs. <0.5% for every other pulse) and is
exactly the "neighbor overlap...makes the per-pulse scan window-arbitrary" risk the review's
sequencing note flagged in advance. The section-8 audit (edge-pinning, reproducibility, Delta-chi2,
flat-zone, window-widening) passed cleanly regardless: it verifies internal
consistency/reliability of the fit *given the residual it was handed*, not freedom from
neighbor-subtraction contamination or absolute goodness-of-fit -- a real scope boundary of the
audit machinery, not a bug in it. A tighter per-pulse window (rather than the shared (-2, 10) used
here) or a proper joint multi-pulse 3-param fit (spec section 5's "3n parameters" default, not yet
implemented) would be the next step to check whether this specific artifact resolves.

**Archived r population** (`load_project_r_population`, 83 values across all 4 GRBs): median 3.31,
p16/p84 = [0.10, 87.9] -- spans nearly 9 orders of magnitude (6e-5 to 1.9e4) end to end. This is a
very broad, clearly non-clustered distribution at the whole-project level, which is why the
archived-median fallback fired for most of this burst's own pulses above; it also argues against
the spec's optional hierarchical (one r per episode) variant being well-supported *project-wide*
-- though the real test per spec section 5's last paragraph is whether *well-determined* r values
within a single episode cluster, which needs a quality filter (e.g. sharp-profile pulses only)
this quick population check didn't apply. Left as an open question in
`norris_3param_reduction.md`.

**Cosmetic fix found along the way:** `plot_diagnostics.plot_fit_overlay`'s fit-curve legend label
was the literal string `"pulse3 fit"` (naming the `pulse3.py` *model*), which reads as "pulse #3"
once real numbered pulses are on the same plots -- renamed to `"3-param fit"`.

**Files written:** `GRB231129C_3param_results.csv` (12 rows, both passes), 6x
`fit_overlay_pulse<N>_pass2.{csv,pdf,png}`.

## GRB231129C joint N-pulse fit (2026-09-25 02:49 PKT)

**Files (new folder `GRB231129C_joint/`, kept separate from `GRB231129C/` per user request):**
`joint_pulse3.py` (N-pulse joint model, bounds, seed-flattening), `run_joint_fit.py` (driver),
`plot_joint_fit.py` (overlay plot with each pulse's own component shown).

Removes the approximation in `GRB231129C/run_pipeline.py`'s method (each pulse fit one at a time
against a residual built from an *estimate* of the other 5): every pulse's (A, t_peak, t_v) is now
fit **simultaneously** in one `curve_fit` call, each keeping the r0 already chosen by the
single-pulse pipeline. Seeded from `GRB231129C/GRB231129C_3param_results.csv` (pass 2) -- read-only
reuse, not a re-derivation. Joint analogs of the section-8 audit were written specifically for this
script (edge-pinning, multi-start reproducibility via small joint-seed jitters, window-widening) --
reusing the exact single-pulse `audit_checklist.py` functions wasn't possible since they assume a
3-parameter, not 3N-parameter, fit.

**Window had to be widened from (-2, 10) to (-2, 20).** The narrower window's own window-widening
audit failed (pulse 4/6 t_v drifted up to ~10%) -- checked directly against the raw light curve
and confirmed real: [10, 16]s has genuine excess flux (~1.9 sigma/bin, decaying to noise by
~t=16-20s), almost certainly pulse 6's own decay tail (broadest pulse, tau2~1.8-2.0s) being cut
off. Widened once, re-verified window-widening PASS. This is the audit doing exactly its job --
catching a real window-truncation problem, not a method flaw.

**Discovery: GRB231129C has two different archived decompositions, not one.** The ORIGINAL
production fit (`fitter_GRB231129779.py`, root level) is 5 pulses, window (-1, 10)s. A LATER,
full-range fit (`GRB231129C/_common.py`, ~5 hours after) added a 6th pulse -- its own P0 entry is
marked `# replaceable` in that source. `GRB231129C/norris_fit_results_GRB231129779_unnormalized.csv`
(the 6-pulse one) is what this project's entire GRB231129C/ and GRB231129C_joint/ analysis has been
built on. Also found: `_common.py`'s docstring claim that its P0 was "copied verbatim from
`fitter_GRB231129779.py`" is false -- the arrays don't match, not even in pulse count. Worth a
`BUGS.md` entry (not filed yet). User confirmed independently (without being told the above) that
pulse 4 is the one in question -- and our own numbers had already flagged it: pulse 4 showed the
largest pass1->pass2 movement in `GRB231129C/run_pipeline.py` and the largest pre-widening
window-sensitivity, both before this discovery.

**5-pulse vs. 6-pulse, tested directly (user request: keep both variants, not just one).**
`run_joint_fit.py` now runs both: `GRB231129C_joint_results.csv` (6-pulse) and
`GRB231129C_joint_5pulse_results.csv` (pulses 1,2,3,5,6 -- pulse 4 dropped, remaining 5 refit
jointly against the SAME raw data, not a residual with pulse 4 pre-subtracted -- the correct,
unbiased test of "does the burst need a 6th pulse").

| | 6-pulse | 5-pulse (pulse 4 dropped) |
|---|---|---|
| chi2 | 349.08 (18 params) | 378.89 (15 params) -- **worse**, despite 3 fewer params |
| audit | edge-pinning none, reproducibility PASS (chi2 spread 0.0 exactly), window-widening PASS | same 3 checks all PASS -- internally stable, just a worse fit |
| pulse 5 vs. its own single-pulse-residual answer | A=11893 (+0.6%), t_v=0.825 (~0%) -- barely moves | A=14841 (**+25%**), t_v=1.046 (**+27%**), t_peak shifts -0.29s |

Visually (compare `joint_overlay_GRB231129C_joint.png` vs. `..._5pulse.png`): the 6-pulse fit's
pulse 4 is a small, narrow bump sitting in a real notch in the data at t~3.2-4.0s (the data itself
dips to ~14,300 then rises to a local peak ~16,500 there, and the 6-pulse total line tracks that
wiggle); the 5-pulse fit's pulse 5 alone smooths over the same region rather than resolving it,
and has to inflate to do so.

**Verdict on 5 vs. 6, as put to the user:** pulse 4 is small in amplitude (9,462, vs. 14,900-20,700
for its 6-pulse neighbors) and has the largest *relative* uncertainty of any of the 6 pulses
(~9.5% on A, ~20% on t_v) -- but unlike the hidden-peak candidate below, it is NOT pinned to any
bound, does NOT fail multi-start reproducibility, and does NOT cause catastrophic neighbor
distortion. Dropping it costs a real, non-marginal Delta-chi2 (~30, well past the ~7-10
usually-justifies-the-parameters threshold used elsewhere in this project) and forces pulse 5 to
unphysically stretch. Recommendation given to the user: keep 6 pulses -- the evidence (clean
audit pass, real chi2 cost, neighbor distortion when dropped) favors it over the historical
"replaceable" judgment, which predated this audit machinery. Final call is the user's; open item
below.

**Caveat:** the 5-pulse variant tested here is "6-pulse minus pulse 4, refit" -- it reuses pulses
1/2/3/5/6's own 6-pulse-fit shapes as the seed, not the independently-optimized parameter values
from the *original* production 5-pulse fit (different window, different P0, likely different local
optimum for the remaining 5). This is the directly relevant test for "does the burst need a 6th
pulse" (same data, controlled comparison), but is not literally a comparison against that original
5-pulse fit's own numbers.

## Hidden peak near t=2s (2026-09-25 02:49 PKT)

User reported a candidate peak near t=2s from a separate ~2-hour multi-start search on another
machine, converging only once or twice -- no known initial parameters for it. Directly suited to
this pipeline's neutral seed (`bounds_seeding.py`), which is argmax-based, not P0-based.

**Step 1 -- checked the residual directly.** Real signal: two *adjacent* bins in the 6-pulse
fit's residual, t=1.856s and t=1.984s, sit at +2.15 and +3.43 sigma above the model, immediately
followed by a compensating dip (down to -2.49 sigma) as pulse 2's rise gets pulled to partly
absorb it. Adjacent-bin correlated excess is much less likely to be pure noise than an isolated
outlier.

**Step 2 -- isolated single-pulse test** (`test_hidden_peak.py`): fit a candidate pulse against
the 6-pulse residual in window (0.5, 2.3)s, using `bounds_seeding`'s neutral seed (no P0).
Result: chi2_null (no pulse) = 42.87, chi2_with_pulse (best r0) = 35.28, Delta-chi2 = 7.59 for 3
params -- marginal, right at the usual justification threshold. Amplitude significance 2.12 sigma.
**Audit FAILED**: t_v landed exactly on the 2*dt resolution floor (0.128s = two bin-widths) and
multi-start reproducibility failed. The fitted pulse (peak ~1049 counts/s) visibly undershoots the
actual spike (~2472 counts/s at t=1.984s) in `fit_overlay_candidate_hidden_peak.png` -- the model
is structurally prevented from chasing a feature narrower than its own resolution floor allows.
Finding: this data's time resolution (64ms bins, confirmed by the user to be the finest available
for this burst -- no finer product exists to switch to) is the limiting factor, not the search
strategy; a smooth Norris pulse shape is fundamentally fighting a spike-like feature at or below
that resolution.

**Step 3 -- joint 7-pulse test** (`run_joint_fit_7pulse.py`, per user request: "one test; if it
works it works, otherwise it may just be a statistical fluctuation caught in a rare minima"):
added the candidate as a genuine 7th pulse, fit simultaneously with the other 6 (not against a
frozen residual). From a well-chosen seed (built from step 2's own near-answer), this looks
*much* stronger in isolation: chi2_7pulse=320.00 vs chi2_6pulse=349.08, Delta-chi2=**29.1** (vs.
7.59 isolated), pulse 7 at 3.57 sigma, nothing pinned to an edge, stable under window-widening.
**But multi-start reproducibility failed badly**: 3 small seed perturbations produced a chi2
spread of 12.9 (vs. 0.000000 for the real 6-pulse fit) and pulse 7's own t_peak spread was **10.5
seconds** across those 3 nearby starts -- the signature of a narrow, fragile local optimum that a
good seed happens to find, not a robust global one. Also: adding pulse 7 forced pulse 2 to lose
9% amplitude and narrow 27% in t_v, and pulse 3 to gain 12% amplitude and widen 17% -- two
neighbors reshaping by double-digit percentages is consistent with real parameter degeneracy
between "add a 7th pulse" and "let 2 and 3 cover the same region differently," not a clean,
independent detection.

**Verdict, in the user's own framing: this is the rare-minima case, not "it works."** Both
independent tests (isolated: pinned at the resolution floor; joint: unreproducible across nearby
seeds) point the same direction. Matches the user's 2-hour, converged-once-or-twice experience
exactly -- there is a narrow basin a well-chosen search occasionally finds, not a stable answer a
thorough search should reliably reach. **Decision (user, 2026-09-25): do not include the hidden
peak in the joint pipeline.** Real Delta-chi2 exists at both the isolated and joint level (7.6 and
29.1 respectively) and the same t~2s location was independently flagged before either test ran, so
this is not dismissed as pure noise -- it is left out on reliability grounds (fails reproducibility
and/or the resolution floor in both tests run), consistent with the section 7/8 machinery's whole
purpose: distinguishing "a real bump exists" from "this specific fit reliably describes it."

**Files written:** `GRB231129C_joint_results.csv`, `GRB231129C_joint_5pulse_results.csv`,
`joint_overlay_GRB231129C_joint{,_5pulse,_7pulse}.{csv,pdf,png}`,
`fit_overlay_candidate_hidden_peak.{csv,pdf,png}`.

**Bug found and fixed along the way:** `plot_joint_fit.plot_joint_overlay`'s legend originally
labeled pulses by their *position* in the params list (1, 2, 3, ...), not their real pulse index --
once pulse 4 was dropped for the 5-pulse variant, the legend's "pulse 4" was silently actually real
pulse 5, and "pulse 5" was real pulse 6. Fixed with an explicit `pulse_indices` argument; verified
by rerunning and confirming the corrected legend against known amplitudes.

## Window x.min/x.max sensitivity, tested directly (2026-09-25 02:57 PKT)

User question: is this joint fitter affected by the fit window's x.min/x.max the way NorrisFitter
is? (`NorrisFitter.fit_boundaries()` ties `t_s`'s lower bound to `x_values.min()`; widening the
window gives the `t_s<->tau1` degenerate direction more room to drift -- the mechanism the
GRB140206275 window-sensitivity diagnostic exists to catch.) `bounds_seeding.py` ties `t_peak` and
`t_v`'s bounds to the window edges the same structural way, so the question is legitimate.

**Answer, tested directly rather than argued from the section-8 window-widening check alone**
(`run_joint_fit_fullspan.py`, new): refit the 6-pulse joint model over the light curve's ENTIRE
x-span (`t.min()`/`t.max()` explicitly -- not a literal `-inf`/`inf` mask, same effect against
real finite timestamps but explicit about what it resolves to, per user preference), i.e. the
data's own full extent, matching `fitter_GRB231129779.py`'s own `START1=-np.inf, END1=np.inf`
convention (the most extreme "widening" possible, not an extrapolation from smaller ones). Plot
boxed to `[T05-PLOT_PAD_S, T95+PLOT_PAD_S]` = `[-4.616, 12.296]`s (production convention, copied
into `GRB231129C/load_data.py` from `../../../GRB231129C/_common.py`) -- plot-axis bounding only,
same documented meaning as there.

| pulse | t_peak (windowed, -2..20) | t_peak (full span, -138.5..475.8) | delta | t_v (windowed) | t_v (full span) | delta |
|---|---|---|---|---|---|---|
| 1 | 0.62619 | 0.62622 | +0.00003 | 0.80978 | 0.80981 | +0.00004 |
| 2 | 2.46630 | 2.46629 | -0.00001 | 0.53367 | 0.53368 | +0.00001 |
| 3 | 1.31753 | 1.31750 | -0.00003 | 0.63321 | 0.63315 | -0.00006 |
| 4 | 3.46383 | 3.46370 | -0.00014 | 0.40968 | 0.40954 | -0.00014 |
| 5 | 4.27884 | 4.27886 | +0.00002 | 0.82542 | 0.82564 | +0.00022 |
| 6 | 5.56213 | 5.56243 | +0.00031 | 1.37573 | 1.37618 | +0.00045 |

Max relative difference: t_peak 0.0055%, t_v 0.0341%. Edge pinning: none. Multi-start
reproducibility: PASS (chi2 spread 4e-6, vs. 0.0 for the windowed fit -- both effectively exact).
**Confirmed directly: this fitter is not sensitive to the window's x.min/x.max the way NorrisFitter
is**, even under the most extreme widening possible (22s window -> the full 614s light curve).

**Bug found and fixed along the way:** `run_joint_fit.py`'s `check_edge_pinning` used a tolerance
relative to the bound *span* (`hi - lo`), matching `audit_checklist.py`'s single-pulse `_is_pinned`
helper. That's fine for a normal-sized window but breaks down once the span becomes huge relative
to the fitted value -- the full-span fit's `t_v` upper bound is 614.3, so `1e-3*span ~= 0.6s` of
absolute slop, enough to falsely flag ordinary ~0.5s `t_v` values (pulses 2/3/4) as "pinned" when
they were nowhere near either edge. Fixed to match `profile_scan.py`'s already-correct
`_is_pinned_r` (tolerance relative to the bound *value*, `np.isclose`-style, not the span).
Re-verified the normal-window 6- and 5-pulse runs still show "edge pinning: none" after the fix
(no regression).

**Follow-up (2026-09-25 03:01 PKT): ported the same fix to `audit_checklist.py`'s single-pulse
`_is_pinned`**, which had the identical latent span-relative-tolerance bug (per user request).
Verified: (1) `test_audit_checklist.py` -- ALL PASS, no regression; (2) `GRB231129C/run_pipeline.py`
rerun on real data -- bit-for-bit identical pass1->pass2 deltas to the pre-fix run, still 6/6
pulses pass audit; (3) `test_hidden_peak.py` rerun -- still correctly flags `t_v` as pinned (it
genuinely sits at the 2*dt floor in that normal-sized window), confirming the fix discriminates
real pinning from the false-positive case rather than just suppressing all pinning detection; (4)
full 9-script synthetic suite -- ALL PASS. `_is_pinned` now matches `profile_scan._is_pinned_r`'s
value-relative form everywhere it's used in this folder.

**Files written:** `run_joint_fit_fullspan.py`, `joint_overlay_GRB231129C_joint_fullspan.{csv,pdf,png}`.
`GRB231129C/load_data.py` gained `T05`, `T95`, `PLOT_PAD_S` (copied from `_common.py`, same
convention as that module's existing copied constants).

## Open item: 5 vs. 6 pulses for GRB231129C

Recommendation given (keep 6, see above); awaiting the user's final call before treating either as
"the" GRB231129C decomposition for downstream work (the other-3-GRB run below, and any eventual
paper section) -- this determines what those next steps actually run against.

## GRB131014215: blind multi-pulse discovery (2026-09-25 03:42 PKT)

**Folder:** `GRB131014215/` -- own dedicated chronological log at `GRB131014215/PROGRESS.md`
(mirrors this file's convention); this is a summary, not the full detail.

Genuinely different scope from GRB231129C: a **blank-slate test** (user's framing) of whether this
pipeline can discover a burst's pulse decomposition with NO archived parameters read at all --
unlike GRB231129C/GRB231129C_joint, which built on that burst's existing 6-pulse production fit.
For GRB131014215 the user explicitly withheld the archived fit and only gave a hint without values
("5 are fitted; but there are potentially 7 total peaks... not outside of [T90] duration"), later
sharing 5 rough eyeballed peak-center guesses only (t ~ 0.5, 1.3, 2.0, 2.7, 3.5) -- clarified as
fair game since eyeballing centers is how every archived P0 in this project was built; the
constraint was specifically about not reading already-fitted values.

**Only legitimate inputs used:** the raw `.dat` light curves, `results.json`'s catalog-level `T90`
field (observational metadata, not a fitted parameter), and the user's 5 (then 7) rough centers.

**What worked, validated:** a naive global-window greedy peak-finder failed (one broad "average"
pulse swallowed several real peaks); a localized version and, more successfully, a joint 5-pulse
fit seeded from the 5 hint centers (mirroring GRB231129C_joint's method) recovered **all 5 real
peaks in the right locations**, reproducibility exact after fixing the r0 choice (see below) --
**3 of 5 pulses landed within 0.001-0.014s of the archived fit's own t_peak** once the user
voluntarily shared it for comparison (t_v agreement within 0.03-0.08s, one pulse essentially
exact). User-confirmed milestone: 5-pulse recovery **passed**.

**Real methodological finding, not specific to this burst:** window size is not a universal safety
property of this 3-param model the way GRB231129C_joint's full-x-span test suggested -- fitting
this burst's full 478s span made a real instability *worse* (reproducibility chi2 spread
1733->3781), because Norris pulses are tail-heavy (`exp(-x/tau2)` never truly reaches zero) and
this burst's own tightly-packed structure gives that tail real room to wander once handed
uninformative-but-not-flat baseline. The actual fix was re-deriving each pulse's r0 from the
joint-fit-consistent residual across several alternating passes (not window size at all) --
reproducibility chi2 spread went from 4060 (worst window attempt) to **0.000000** (exact).

**Two real bugs found and fixed, both ported back to the canonical GRB231129C_joint source too:**
a module-name collision (`load_data.py` exists in both `GRB131014215/` and `GRB231129C_joint/`;
fixed per CLAUDE.md's copy-not-cross-import convention) and a jitter-outside-bounds crash in the
joint multi-start reproducibility check (fixed by clamping; re-verified no regression on
GRB231129C_joint's own tests).

**Open, not resolved this session:** pushing to 6-7 pulses confirmed one real effect (an explicit
tail-anchoring 6th pulse lets pulse 1 relax back to its correct, narrow, archived-matching shape)
but surfaced a recurring difficulty -- any extra pulse placed within ~0.2-0.3s of an existing one
destabilizes into an identity swap or a resolution-floor/grid-cap pin, at two different tested
locations. Session closed at the user's request with the 5-pulse result as the validated baseline;
6/7-pulse work is logged as open for a future session, not as a failure of the method.

## Next candidates

- Resolve the 5-vs-6-pulse open item for GRB231129C (above).
- Resolve the 6-vs-7-pulse open item for GRB131014215 (`GRB131014215/PROGRESS.md`'s own "Next
  candidates" -- start from the validated 5-pulse baseline, try a wider/adaptive r0 grid before
  adding more tightly-spaced hints in the t~2-3.5s region).
- Investigate pulse 6's neighbor-contamination artifact from the single-pulse (non-joint) method --
  likely already resolved by the joint fit (pulse 6 no longer relies on an approximate residual),
  worth explicitly confirming.
- The hierarchical r-sharing variant's real test (well-determined, non-degenerate archived r values
  within one episode, filtered by fit quality, not the whole-project population checked earlier).
- Run the same joint pipeline on the other 2 archived GRBs (GRB080916C, GRB140206B) -- now that
  both GRB231129C and GRB131014215 have working joint 3-param fitters.
- File the `_common.py` P0-docstring/pulse-count inconsistency discovered during the GRB231129C
  work as a `BUGS.md` entry (separate from this project's own scope, but a real issue in the
  production code).

## Build order note (2026-09-25)

Sections were implemented out of the spec's own 1-8 order, driven by dependencies rather than
document order: 6 (consistency, explicitly "read before implementing" in the spec) -> 1-2 (model)
-> 3 (bounds/seed) -> 5 (profile scan) -> 4 (uncertainty propagation, blocked on 5's flat-zone
output). Remaining work (7's T3-T5, then 8) follows the same rule: T3 has no dependencies and could
run anytime; T4/T5 and section 8 depend on each other in that order.
