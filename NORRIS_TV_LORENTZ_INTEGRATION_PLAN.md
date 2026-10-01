# Wiring measured Norris-profile t_v into lorentz_factor.py's Gamma_min (working plan)

Companion to `PLAN.md` Phase 5, `PHASE5_TV_PLAN.md`, and
`codes-for-paper/lorentz_factor/lorentz_factor.md` §6b/§11. Written 2026-09-27 — this is the plan
for the "feed Phase 5's Norris measurements back into `Gamma_min`" step that `PLAN.md`/
`PHASE5_TV_PLAN.md` both flag as deliberately out of scope for the fitting work itself.

**Status: implemented 2026-09-27, same day.** All four decisions below were carried out exactly as
written; the resulting numbers, verification results, and table/documentation changes are recorded
in `lorentz_factor.md` §11 (not duplicated here — this file stays as the pre-implementation plan
and rationale). `GRBResearchPaper` remains untouched, per the scope boundary below.

## Context

`lorentz_factor.py` (Limit A) and `lorentz_factor_limit_b.py` (Limit B) currently take the
variability timescale $t_v$ as each episode's **duration** — a deliberate upper bound adopted
2026-08-21 specifically because no measured $t_v$ existed yet (`lorentz_factor.md` §4). Phase 5
has since produced an actual measurement: a manual, per-burst joint Norris-pulse fit to the
10–400 keV summed light curve, living in `codes-for-paper/variability_analysis/`, with per-pulse
results (including MC-propagated $t_v$ and its uncertainty) in
`norris_fit_results_GRB<name>.csv` for all four bursts. That work is now stable — BUG-25 (bad
seeding) and BUG-26 (NaN overflow) are both fixed, and the four canonical fitter scripts were
restructured 2026-09-27 — but per `PLAN.md`'s own stated deliverable boundary, none of it has been
fed back into `Gamma_min` yet. This plan is that feed-back step.

**Why now:** the Norris measurement replaces a known-conservative upper bound with an actual
measurement, which is the single largest remaining systematic in $\Gamma_\text{min}$
(`lorentz_factor.md` §5: a factor of ~2 spread between candidate $t_v$ conventions, versus ~1-8%
statistical error). It only changes a *published* number for GRB080916C (the only burst with a
redshift), but the CSV computes informational rows for all four bursts today, so the change
threads through all of them for consistency.

**Explicit scope boundary, matching precedent (BUG-18, the RNG overhaul):** this plan covers
`GRBResearchWork` only — regenerating `lorentz_results.csv`/`lorentz_results_limit_b.csv`/both
`.tex` tables/`gamma_comparison` figure, and updating `lorentz_factor.md`. It does **not** touch
`GRBResearchPaper`. Paper integration (new prose, `VERIFY.md` entries, copying the regenerated
tables) is a separate follow-up once the user has reviewed the new numbers — every past
$\Gamma_\text{min}$ change (197→258→134→507) was verified this way before touching the paper.

## Decisions locked in (confirmed with user, 2026-09-27)

1. **Multi-pulse → episode selection rule** — **Superseded 2026-10-02 (user decision):** the pulse is now the one the Norris fitters themselves assign the episode's photon to (nearest preceding onset among still-active pulses), not the nearest `t_peak`; see `lorentz_factor.md` §15. The text below describes the original rule.

   Original:  when more than one fitted pulse's `t_peak_s` falls
   in an episode's window, pick the pulse whose `t_peak_s` is closest to that episode's own
   $\Gamma_\text{min}$-defining LAT photon arrival time (`t_arr_s`, from `LAT_PHOTONS`/
   `lat_photons.csv`, already loaded in `lorentz_factor.py`). This is the exact rule already
   recorded in `PHASE5_TV_PLAN.md` for the (abandoned) automated pipeline, applied here to the
   actual manual-track CSVs instead. It resolves GRB131014A's TR2 ambiguity (pulse 5, the 1021 MeV
   photon's own pulse, over pulse 4) and generalizes cleanly, since every episode this rule is
   ever invoked for is already restricted to `episode in photons` (i.e. `t_arr_s` always exists).

2. **GRB080916C's missing EX0/EX1 bounds are a prerequisite fix**: `shared_utilities.py`'s
   `grb080916C_bounds()` is the only one of the four bursts' bounds functions without EX0/EX1
   entries, so today only TR1–TR5 are eligible for a measured $t_v$ for this burst (the one that
   matters). Add `"EX0": (-0.128, 4.864), "EX1": (59.520, 67.904)` (values to be re-confirmed
   against `results.json`'s own interval strings before editing, not hand-copied from memory),
   matching the other three bursts' pattern and `CLAUDE.md`'s EX0/EX1 definition. Checked in
   advance: EX0 maps unambiguously to pulse 1, EX1 to pulse 6 — no new multi-pulse conflict is
   introduced by this addition.

3. **Quality gate**: a Norris measurement is only used if the selected pulse's `mc_kept_fraction`
   is at or above a threshold (default **0.5**); below it, fall back to duration, exactly as
   today, but record why. This matches `variability_analysis.md`'s own flagged low-confidence
   cases (e.g. GRB131014A pulses 1–2, `mc_kept_fraction` 0.10–0.21) — the gate should reject
   exactly those and no others; that's a concrete check in Verification below.

4. **Error propagation**: once $t_v$ is a real measurement (not a fixed bound), its own MC
   uncertainty (`t_v_err_lower_s`/`upper_s`) should show up in $\Gamma_\text{min}$'s error bars,
   per `CLAUDE.md`'s "every MC-derived quantity carries an error" rule. The underlying per-draw
   $t_v$ samples aren't persisted (only 16/50/84 percentiles are, in the CSV), and re-deriving them
   exactly would mean importing the fitted `NorrisFitter` object across folders — against
   `CLAUDE.md`'s "copy rather than fight `sys.path`" convention. **Approximation, to be stated
   explicitly in `lorentz_factor.md`**: resample $t_v$ per MC draw from a two-piece (split-normal)
   distribution built from `t_v_s` and its asymmetric lower/upper errors. This is a documented
   approximation, not exact reproduction of the Norris fit's own MC — acceptable given
   $\Gamma_\text{min}$'s proven weak sensitivity to $t_v$ (exponent ≈0.15, `lorentz_factor.md` §5).
   Duration-sourced rows are unaffected (no resampling — a duration is still a fixed bound, not a
   distribution).

## Implementation

### Step 0 — Backups (before any edit)

- `codes-for-paper/lorentz_factor/_backup_2026-09-27_pre_norris_tv/`: copy every existing file in
  that folder (`lorentz_factor.py`, `lorentz_factor_limit_b.py`, `lorentz_results*.csv`,
  `lorentz_table*.tex`, `lorentz_factor.md`, `gamma_comparison_plot.py`, `gamma_comparison.png/pdf`)
  — same pattern as the existing `experiments/window_sensitivity_GRB140206275/_backup_2026-09-24_pre_split/`
  precedent.
- `codes-for-paper/variability_analysis/_backup_2026-09-27_pre_ex0ex1/`: copy `shared_utilities.py`
  and `norris_fit_results_GRB080916C.csv` before the bounds edit + rerun.

### Step 1 — `variability_analysis/shared_utilities.py`

- Confirm GRB080916C's exact EX0/EX1 bounds against `results.json` (`GRB080916009` entry), then
  add them to `grb080916C_bounds()`'s `episode_bounds` dict.
- Rerun the canonical GRB080916C fitter script to regenerate `norris_fit_results_GRB080916C.csv`.
  Diff old vs. new CSV: expect exactly two new rows (EX0/pulse1, EX1/pulse6) and every existing
  row byte-identical — confirm this, don't assume it.

### Step 2 — `lorentz_factor.py`

- New constants: path to the four `norris_fit_results_GRB*.csv` files (via `find_project_root()`),
  `MC_KEPT_FRACTION_MIN = 0.5`.
- New function `load_norris_tv(photons)` returning `{(short_name, episode): {t_v_s, t_v_err_lower_s,
  t_v_err_upper_s, pulse_index, mc_kept_fraction, rejected_reason}}`, implementing decision 1's
  selection rule and decision 3's gate. `rejected_reason` is set (e.g. `"mc_kept_fraction<0.5"`)
  whenever a candidate existed but didn't pass the gate, so the fallback to duration is visible,
  not silent.
- `variability_timescale()`: extend to accept `t_arr` and the norris lookup; precedence stays
  `VARIABILITY_TIMESCALE` override → Norris (if present and gated in) → duration. Returns the
  value, its source (`"literature"`/`"norris"`/`"duration"`), and error bars (`None` for
  non-Norris sources).
- `compute_gamma_min`'s MC loop in `main()`: when source is `"norris"`, draw `delta_t` per
  iteration from the split-normal approximation (decision 4); otherwise keep today's fixed-scalar
  behavior unchanged.
- CSV output: add `t_v_err_lower_s`/`t_v_err_upper_s` columns (blank/`None` for non-Norris rows),
  per `CLAUDE.md`'s uncertainty rule.
- `build_latex_table()`: replace the current `assert t_v_sources == {"duration"}` /
  single-header-dagger design (already anticipated in its own comment) with a per-cell marker —
  keep `\dagger` for duration-sourced, add a new, currently-unused marker (`\S`) for
  Norris-sourced, both explained in the table's `tablenotes`.
- Update the module docstring / inline comments describing the $t_v$ convention.

### Step 3 — `lorentz_factor_limit_b.py`

- Same `variability_timescale()` call-site update (it currently imports and calls this function
  identically to Limit A) and the same MC delta_T resampling in `compute_gamma_min_limit_b`'s
  loop, for consistency — `lorentz_factor.md` §8.4 already states both limits share every $t_v$
  caveat.
- Same table marker fix as Step 2.

### Step 4 — `lorentz_factor.md`

- New section documenting: the four locked-in decisions above with `user, 2026-09-27`
  attribution, the split-normal approximation caveat, the quality-gate threshold and rationale,
  the EX0/EX1 bounds fix, and a rejected-candidates list (which (grb, episode) fell back to
  duration and why).
- Refresh §3/§8.3's results tables with the new numbers, keeping the old ones marked **STALE**
  rather than deleted — matching the file's existing convention.
- Update §6b's status note to mark the wiring done, cross-referencing this file.

## Verification

1. `py_compile` both scripts after edits.
2. Rerun both scripts (`PYTHONPATH=src .venv/bin/python codes-for-paper/lorentz_factor/lorentz_factor.py`
   and `..._limit_b.py`) and confirm they complete without error.
3. Diff regenerated `lorentz_results.csv`/`lorentz_results_limit_b.csv` against the Step-0 backups:
   confirm T90 rows are unchanged (Norris deliberately never covers T90), and TR1–TR5/EX0/EX1 rows
   for GRB080916C change from `duration` to `norris` source wherever the quality gate passes.
4. Hand-recompute one changed GRB080916C episode's $\Gamma_\text{min}$ from the new CSV's own
   `f_1`/`alpha_LS`/`E_max_MeV`/`z`/`t_v_s` and confirm it matches the pipeline's own value —
   the project's standing "a formula change gets a numeric check, not a read-through" rule.
5. Confirm the quality gate rejected exactly the already-known low-confidence cases (e.g. GRB131014A
   pulses 1–2) and no others — cross-check against `variability_analysis.md`'s own flagged list.
6. Re-verify the Limit A/Limit B shared `tau_hat` bit-identity still holds (same check as
   `lorentz_factor.md` §8.2) now that both scripts read the same Norris-derived delta_T.
7. Regenerate `gamma_comparison_plot.py`'s figure and visually inspect the PNG (this is a
   GRBResearchWork-side plot regeneration, not a paper build, so it's fine for Claude to run and
   check directly, unlike `pdflatex`/`pdftotext`).
8. `git status` in `GRBResearchPaper` to confirm it is untouched.
9. No commits in either repo — leave both working trees for the user's review, and report exactly
   which episodes' $\Gamma_\text{min}$ changed and by how much before any further step (paper
   integration, or a commit) is taken.