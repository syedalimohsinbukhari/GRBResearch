# Plan: Take the BB-component / fractional-ratio / photospheric-physics work to a finished paper

## Context

The paper (`GRBResearchPaper`) and its analysis code (`GRBResearchWork`) already cover spectral fitting (Band/CPL/SBPL/PL with optional BB add-on), model selection (AIC/BIC, BEST/SAFE/UNSAFE classification), and Amati-relation energetics for four core GRBs: **GRB080916C, GRB131014A, GRB140206B, GRB231129C** (`best_names.txt`, 26 intervals, 13 of them BB-inclusive). Sections 1–4 and the appendices are content-complete — they're just commented out of `main.tex` to keep local compile times down while iterating.

What's still missing is the paper's core physical payoff, and it's exactly what you flagged: **the BB component's fractional contribution to total flux, and the physical quantities that follow from it** (photospheric radius, Lorentz factor via the Pe'er 2007 method). This is called out as unchecked in `GRBResearchWork/ToDo.md` under "Blackbody Component Diagnostics" and is not implemented anywhere in `src/grb_research/` or `codes-for-paper/` — confirmed by grep (no `bb_fraction`/`flux_ratio`/`energy_fraction` hits). Section 5 of the paper ("Physics Calculations") is correspondingly thin: only peak-energy evolution, E_iso, and Amati are live; a γγ-opacity Lorentz-factor subsection and a hardness-ratio subsection are drafted but commented out.

Decisions already made with you:
- **Redshift handling**: run the new BB-fraction + Pe'er-method analysis on all 4 core GRBs; for the 3 without a confirmed z (131014A, 140206B, 231129C), sweep a plausible z range rather than assuming one value (per `ToDo.md`'s own suggestion).
- **Sequencing**: get the paper's already-finished content re-enabled and internally consistent first (cheap, low-risk), then build the new science on top of a clean base.
- **Chapters 1–4 and the appendices are already final** — they were commented out purely to save local compile time, not because they're incomplete. Phase 0 below is a re-enable-and-verify pass on them, not a content pass.

**Research scope**: `GRBResearchWork` holds raw-interval directories and `results.json` entries for eight GRBs, but this paper's sample is only four: GRB080916C (`GRB080916009`), GRB131014A (`GRB131014215`), GRB140206B (`GRB140206275`), GRB231129C (`GRB231129779`). The other four directories present in the repo — `GRB110721200`, `GRB110731465`, `GRB150210935`, `GRB190114873` — are not part of this paper and must not be pulled into any of the phases below (figures, tables, sample statistics, new analysis code). This is also now recorded in `CLAUDE.md` for future sessions. It's the reason Phase 0.4 treats `lorentz_factor.py`'s old-sample entries as a bug to fix, not data to keep.

No commits to either repo — everything below is working-tree only, per your standing rule.

**Output standard**: `update_style()` (`src/grb_research/__init__.py:58-105`, not `grb_styles.py` — that file only holds the `GRBPlotStyle` color/marker-convention class) is the existing consistency convention and every new plot must call it. But consistency isn't enough — every new plot and CSV from Phases 1–2 must also be *maximally defensible*: self-contained enough that a reader can verify or reproduce the number without opening the script. See the "Output standard" section below for the concrete column/figure requirements this implies, and note this is now also recorded in `CLAUDE.md` as a standing project convention beyond just this plan.

---

## Output standard for new plots & CSVs (applies to Phases 1–2)

Grounded in an audit of existing outputs: `codes-for-paper/amati_relationship/amati_relationship.csv` is the best existing example (units baked into column names, asymmetric MC errors, GRB/episode/model identity all present) — follow it, not `flux_fluence.csv` (unitless columns) or `codes-for-paper/evolution_of_kt/kt_evolution_080916c.py` (no `update_style()`, no axis labels, no legend, no `savefig` — ends in a bare `plt.show()`, so it must **not** be used as a template despite being the closest existing analog).

**CSVs** — every row self-contained:
- Identity: `grb_name` (paper name, e.g. `GRB080916C`), `episode` (e.g. `T90`, `TR1`), `model_name` (e.g. `BAND_BB`).
- Assumptions used for that row: `z` (or the swept value, for the 3 GRBs without confirmed redshift), `H0`, `Om0` (the single reconciled cosmology from Phase 0.3 — don't let a third value creep in).
- Value + uncertainty as `<name>_<unit>`, `<name>_err_lower_<unit>`, `<name>_err_upper_<unit>` for every MC-derived quantity (`f_BB`, `R0_cm`, `Gamma`, etc.) — never emit a bare median when the resampler already gives percentiles.
- MC provenance: `n_samples`, `seed` as constant-valued columns, so the CSV doesn't depend on a companion file to be trusted. **Updated by the RNG/seeding overhaul, Phase 4c below:** don't hardcode a literal seed or thread `seed=` through function boundaries — each script derives its own via `SEED = seed_from_name(__file__)` then builds one `rng = get_rng(seed=SEED)` and threads `rng=` (only) through every downstream call, per `src/grb_research/SEEDING.md`.
- Nothing beyond that — no intermediate arrays or debug columns; stop at what's needed to recompute or sanity-check the final number.

**Plots** — every new figure:
- Calls `update_style()`.
- Axis labels with explicit units.
- Uses `GRBPlotStyle` (`src/grb_research/grb_styles.py`) for per-GRB color / per-episode marker / BB-vs-non-BB fill, matching the paper's existing figure conventions (solid = no BB, hollow = BB-augmented).
- Error bars/bands wherever the underlying value has an MC uncertainty.
- Explicit `plt.savefig(..., dpi=SAVE_DPI)` to both `.png` and `.pdf` (follow `amati_relationship.py`'s save loop) — never a bare `plt.show()`.

## Phase 0 — Re-enable finished content & fix known inconsistencies (`GRBResearchPaper`)

1. In `main.tex`, uncomment the `\input`s for `tex_files/abstract`, `tex_files/section-6-discussion`, `tex_files/section-7-conclusion`, and the appendix block (`appendices/table_*_large.tex`, `tex_files/generated/lat_info_table.tex` — `\input` directly, no wrapper file as of 2026-09-06). These are done content, per you — just wire them back in and recompile to catch any stale `\ref`/`\cite`/figure-path breakage from being disconnected.
2. Recompile (`pdflatex` → `bibtex` → `pdflatex` ×2, matching the existing `out/` artifacts) and fix any compile errors that surface from the re-enabled content (undefined refs, missing figures, etc.) — content edits only where needed to make it build, not a rewrite.
3. **Cosmology mismatch**: `tex_files/section-5-data-analysis.tex` (Isotropic Energy, line ~77) uses $H_0=69.6$, $\Omega_m=0.286$ (Fana Dirirsa 2019), while `codes-for-paper/lorentz_factor/lorentz_factor.py:36` uses $H_0=67.4$, $\Omega_m=0.315$. Pick one (recommend keeping Fana Dirirsa 2019 since it's already the Amati-relation reference) and update `lorentz_factor.py` to match.
4. **Stale GRB sample in `lorentz_factor.py`**: `grb_inputs` (lines 125–179) mixes the old sample (GRB110721A, GRB110731A, GRB150210A) with the current one (GRB080916C). Replace the three old-sample entries with the current sample's other three GRBs (131014A, 140206B, 231129C; `z: None` since unconfirmed), so the regenerated `lorentz_table.tex` lines up with the rest of the paper's 4-GRB framing — same pattern already used in the Amati/E_iso table.
5. Leave the abstract's headline numbers as a placeholder for now (marked clearly, e.g. `% TODO: fill after Phase 3`) — it needs the new BB-fraction/photospheric results before it can be finalized, so don't rewrite it twice.
6. Deduplicate `images/section4` vs `images/section5` (near-identical `amati_relationship`, `butterfly_all`, `peak_energy_best__all` PNGs) — check which set the re-enabled sections actually reference and drop the stale copies.
7. Leave the "Hardness-to-Softness Evolution" subsection commented out — no hardness-ratio code exists in `GRBResearchWork`, and it's orthogonal to the BB-fraction direction. Flag as future work, not in scope here.

## Phase 1 — BB fractional flux/energy ratio (`GRBResearchWork`)

New module, e.g. `codes-for-paper/bb_fraction/bb_flux_fraction.py`, built on existing infrastructure:

- **Component separation** is already there: `SpectralModels._evaluate_components` (`src/grb_research/grb_sed.py:96-110`) returns per-component spectra + total for any BB-composite model via `MODEL_MAP` (`grb_sed.py:14-25`). Use `model_type="energy"` (`grb_sed.py:128-129`) so each component is energy flux, then integrate over the same 1 keV–10 MeV bolometric band already used for fluence (`tex_files/section-5-data-analysis.tex:73`, matching `Bloom2001`/`Nava2012` convention) via `np.trapz` over the log-spaced grid `_compute_model` already builds.
- `f_BB(t) = ∫BB dE / ∫total dE` per interval, computed for all 13 BB-inclusive intervals across the 4 core GRBs (params pulled from `results.json`, `amp_bb`/`kt_bb` keys confirmed present).
- **Uncertainty propagation**: reuse the existing per-model MC resamplers in `grb_calculations.py:259-286` (`_pl_bb_resampler`, `_cpl_bb_resampler`, `_band_bb_resampler`, `_sbpl_bb_resampler`, `_cpl_pl_bb_resampler`, `_band_pl_bb_resampler`, `_sbpl_pl_bb_resampler`) to draw parameter samples, recompute `f_BB` per draw, and take percentiles — same MC convention already used for `E_iso` (`mc_e_iso_sampler`, `grb_calculations.py:333`).
- **Output**: per-GRB time-resolved `f_BB(t)` plot meeting the Output standard above (style/legend conventions from `codes-for-paper/amati_relationship/amati_relationship.py` and `src/grb_research/grb_styles.py` — **not** `kt_evolution_080916c.py`, which lacks styling, labels, and even a `savefig` call); plus a results table (CSV + LaTeX), following the `lorentz_results.csv`/`lorentz_table.tex` generation pattern in `lorentz_factor.py:221-288`, upgraded to the new column standard (units, asymmetric errors, MC provenance).

## Phase 2 — Photospheric radius & Lorentz factor, Pe'er (2007) method (`GRBResearchWork`)

New module, e.g. `codes-for-paper/photospheric_radius/pe_er_photosphere.py`:

- Before coding, read the Pe'er 2006/2007 PDFs in `GRBResearchPaper/Literature Review/` to pin down the exact R0/Γ equations and the `ℛ` (blackbody-radius) parameter definition — don't hand-derive from memory.
- Inputs per interval: `kt_bb`, `amp_bb` (from `results.json`), the BB fraction / $Y = F_\text{total}/F_\text{BB}$ ratio from Phase 1, and $d_L(z)$ via `astropy.cosmology.FlatLambdaCDM` using the cosmology settled in Phase 0.3.
- **GRB080916C**: confirmed $z=4.35$ → single physical $R_0$, $\Gamma$ per interval.
- **131014A, 140206B, 231129C**: per your call above, sweep a plausible z range and report $R_0(z)$, $\Gamma(z)$ sensitivity curves/tables rather than a single assumed value. Confirm the actual z-grid bounds with domain judgment before finalizing (e.g. typical long-GRB range) — flag this as a value judgment, not something to silently invent.
- **Optional/stretch** (only if time remains after the core R0/Γ derivation): fireball baryon-loading $\eta$ and magnetization $\sigma_0$ from $\Gamma$, $R_0$ — third unchecked `ToDo.md` bullet under "Blackbody Component Diagnostics." Superseded/scoped by Phase 6 below.
- Output: results table (CSV + LaTeX) + plot meeting the Output standard above — same generation pattern as `lorentz_table.tex`/`lorentz_results.csv`, upgraded to the new column standard (units, asymmetric errors, MC provenance, z/cosmology used per row — critical here since 3 of the 4 GRBs use a swept z range rather than one value).

## Phase 3 — Paper integration (`GRBResearchPaper`)

- Add two new subsections to `tex_files/section-5-data-analysis.tex`, after "Amati relationship" (~line 143): **"Blackbody Fractional Flux Contribution"** (Phase 1 figure + table) and **"Photospheric Radius and Bulk Lorentz Factor"** (Phase 2 figure + table) — presented alongside the existing (now re-enabled per Phase 0) γγ-opacity $\Gamma_\min$ subsection as a complementary constraint (opacity = lower bound, Pe'er = model-dependent estimate).
- Add a few sentences to `tex_files/section-6-discussion.tex` and `tex_files/section-7-conclusion.tex` incorporating the new findings — these are otherwise finished, so this is addition, not rewrite.
- Finalize `tex_files/abstract.tex`: replace the stale old-sample sentences (GRB110721A/110731A/150210A leftovers) with the current 4-GRB sample, and fill in the Phase 0.5 placeholder with real headline numbers from Phases 1–2.
- Update `GRBResearchWork/ToDo.md`: check off "Multi-component decomposition" and "Photospheric radius & Lorentz factor" (and "Fireball parameter estimation" if the stretch goal lands).

## Phase 4 — Final QA

- ✅ Full recompile of `main.tex` — zero errors, zero undefined refs/citations, 23 pages. Note a clean build proves nothing about content: verification is by `pdftotext` + grep, after a whole subsection once sat commented out while the build reported success.
- 🟨 `ruff` is **not installed** despite `ruff.toml` existing. Substitutes (`py_compile`, `awk 'length>120'`, an AST unused-import scan) have been run over `LAT_analysis/*.py` and `lorentz_factor.py`; **`bb_flux_fraction.py` and `pe_er_photosphere.py` have not been linted.**
- ⬜ Spot-check new table values against `results.json` by hand — still outstanding for the Phase 1/2 tables, and the largest untested surface. The LAT table is now independently verified (see Phase 4b).
- ⬜ No commits in either repo — leave the working trees as-is for your review.

## Phase 4b — Unplanned: LAT data layer and the MeV/keV units fix

*Not in the original plan; done 2026-08-21 after Phases 0–3. Full account in `HANDOFF.md` §9, `BUGS.md` BUG-18 and OBS-08, and `LAT_analysis/LAT_analysis.md`.*

Checking whether `LAT_analysis/` could replace `lorentz_factor.py`'s hand-transcribed `LAT_PHOTONS` (OBS-08) turned up a units error: the LAT photon energies are **MeV**, but `compute_gamma_min` normalises by $m_ec^2$ in **keV**. $\Gamma_\text{min}$ was low by 3.4–4.5×.

- **Fixed** (BUG-18), and **BUG-16 reverted** — it had "corrected" the appendix header from MeV to keV on an arithmetic slip.
- **The narrative moved more than the numbers**: the thermal $\Gamma$ exceeds the opacity floor by **1.5–2.4**, not "five to ten". Reworded in `abstract.tex`, §5 (×2), §7 (×2).
- **OBS-08 closed.** `LAT_analysis/txt_to_csv.py` now builds `lat_photons.csv` from the raw gtlike output; `csv_to_latex.py` renders the appendix table from it; `lorentz_factor.py` reads it. No hand-typed LAT numbers remain, except the `\pFlux` upper limits, which are genuinely not derivable from the files in the repo and are guarded by an assertion instead.

## Phase 4c — Unplanned: RNG/seeding overhaul

*Not in the original plan; reviewer-driven, done 2026-09-02/03 after Phase 4b. Full account in `SEED_PLAN.md` (the plan) and `SEED_PLAN-implementation.md` (what actually landed, item by item, with verification results); the scheme itself is documented in `src/grb_research/SEEDING.md`.*

A paper reviewer flagged that Monte Carlo draws across the codebase reused the same literal seed (`12345`, sometimes `1234`/`42`) independently in many scripts. Investigation found something worse: several scripts rebuilt an identically-seeded generator on *every iteration* of a loop over models/GRBs/episodes — a correctness bug (silently correlated/duplicated draws across iterations meant to be independent), not just a style concern.

- **Fixed** by deriving a deterministic per-script seed from `(script filename, MASTER_SEED)` via SHA-256 (`seed_from_name`), replacing every hardcoded literal; each script builds one `rng` and threads it through every downstream call instead of re-passing `seed=`. `MASTER_SEED = 2828702241`, drawn via `secrets.randbits(32)` with its exact provenance recorded in `SEEDING.md`, so it's auditable rather than hand-picked.
- **Two real latent bugs found and fixed along the way** (BUGS.md BUG-19, BUG-20): `convert_sbpl_to_band` silently overwrote a caller-supplied `rng` whenever `seed` was also given, and `FluxFluenceCalculator.__init__`'s falsy `if seed:` check silently dropped `seed=0`. Neither was ever triggered by an actual caller, so no published number was affected.
- **A follow-up cleanup pass** stripped the now-vestigial `seed=` fallback parameter from 11 downstream functions across `src/grb_research/` and `codes-for-paper/`, making `rng` required everywhere — verified behavior-preserving by re-running every affected script and diffing output against the pre-cleanup CSVs (byte-identical in all 8 cases fast enough to check that way; the two slow scripts were smoke-tested at a reduced sample count instead, per the "Determinism checks on slow MC scripts" convention this added to `CLAUDE.md`).
- **Integrated into the paper**: a two-sentence summary in `section-5-data-analysis.tex`, the per-script seed cited in every genuinely Monte Carlo-derived figure's caption, and a new generated appendix table (`tab:seed_table`) listing every seeded script — built from `codes-for-paper/seed_registry.yaml`/`figure_registry.yaml`/`table_registry.yaml`, the same registry pattern used to drive the new `sync_paper_assets.py`, which now copies every figure/table this paper uses from `GRBResearchWork` into `GRBResearchPaper` in one pass instead of a fully manual copy.

## Phase 5 (optional) — Measured variability timescale via Norris-profile fits

*Added 2026-08-21 after assessing Bukhari et al. (2022) §6.5. Feasible; improves a systematic rather than fixing an error.*

**Status update (2026-09-18): all four bursts now have a manual fit; still not fed back into the pipeline.** `GRBResearchWork/PHASE5_TV_PLAN.md` records the judgement calls this section calls for (pulse identification per episode, the rule for which pulse defines $t_\text{v}$ in a multi-pulse episode, EX0/EX1 window handling) and is the current working plan — read it, not just this section, before touching Phase 5. Implementation has moved through two approaches: an automated MEPSA/peak-detection pipeline (attempted, stalled on GRB080916C's TR2, since abandoned — its source files are deleted) and a manual, GRB-by-GRB joint-Norris-fit track in `codes-for-paper/variability_analysis/` (`variability_analysis.md`), which as of this date covers all four bursts, including GRB080916C (redone manually 2026-09-16/17, the last one to be done, at the burst where the automated approach originally stalled). None of this has been fed back into `lorentz_factor.py`'s $\Gamma_\text{min}$ pipeline yet — still an explicit deliverable boundary, per `PHASE5_TV_PLAN.md`. **Open issue, not yet addressed:** every fit so far uses only a single NaI detector rather than summing across detectors, unlike typical published GRB analyses — flagged 2026-09-18, see `PHASE5_TV_PLAN.md`.

**Status update (2026-09-22):** GRB080916C's `fitter.py` reverted from the 8-pulse `norris2-1`/`norris2-2` decomposition back to the standard 7-pulse model (user decision — the split wasn't judged sound); a LAT-photon overlay ($E>1$ GeV) and a temporal-proximity photon-to-pulse assignment (refined with an active-window flux threshold so a decayed pulse stops claiming later photons) were added to that same script. **BUG-23, found this session:** the single-NaI-detector selection flagged above isn't just a missing-summing gap — it's non-deterministic. Every `fitter*.py`/`window_sensitivity.py` picks its detector via `os.listdir()[0]`, which is filesystem order, not alphabetical; a routine data resync silently flipped GRB080916C's detector from the documented `n3` to `n4` mid-session. Fixed for GRB080916C (`sorted()`, verified against the previously-committed results). **Verified still live and unfixed for GRB131014A and GRB140206B** (their documented detector isn't alphabetically first among their actual detector set, so `sorted()` won't fix them — needs an explicit per-burst hardcoded detector name instead); GRB231129C is currently correct by luck, not design. Full writeup: `variability_analysis.md`'s GRB080916C section and `BUGS.md` BUG-23. Not yet fixed for the two exposed bursts, and TR3 still has two unresolved candidate pulses for whichever burst's decomposition eventually feeds `Gamma_min`.

**Status update (2026-09-23): BUG-23 closed for all four bursts, with a structural fix rather than the per-burst hardcode planned above.** Instead of hardcoding each burst's documented detector, every `fitter*.py`/`window_sensitivity.py` now sums *all* of a burst's NaI detectors (raw background-subtracted count rates, no per-detector normalization — user's explicit choice, matching the one existing precedent already in the repo). This closes both the single-detector-summing gap flagged 2026-09-18 above and BUG-23 itself in one fix, since there's no longer a single detector to mis-pick. Re-running all four fits with the summed curves surfaced two further, separate fixes: **GRB131014A** initially failed to converge at all (`NorrisFitter`'s default `max_iterations=5000` wasn't enough for the larger summed-curve scale; fixed with `max_iterations=20000`, confirmed to land on the same solution as before, not a different one). **GRB080916C's 7th pulse** (the one added in the 2026-09-22 revert) converged to a degenerate near-delta-function fit on the summed curve — confirmed against the raw data to be overfitting noise, not a real feature — and was dropped, reverting to a 6-pulse model (verified the other six pulses are undisturbed). GRB140206B and GRB231129C converged cleanly with no changes.

A new standalone folder, `codes-for-paper/variability_analysis/experiments/normalized_vs_unnormalized_fit/`, then asked whether this project's peak-normalization convention (`y /= Y_MAX_CTS_PER_S` before fitting) changes any of this: tested against a from-scratch unnormalized method for all four bursts and found no — SSE agrees to <0.15% either way. It also diagnosed (but has not applied to production) a recommended full-range fit window for GRB140206B, whose two weakest pulses turned out to need the light curve's own full `x.min()`/`x.max()` range plus a specific pulse re-seed to resolve cleanly. `why_not_unnormalized_fitting.md` in that folder is a standalone reviewer-facing backup for "why wasn't the raw light curve fit directly", in case this is asked about the paper's methodology.

**Still open, paper-facing:** whether the *other three* bursts (not just GRB140206B) should also move to the full-range window before their parameters feed `Gamma_min` — a script for this (`full_range_all_bursts_check.py`, same folder) is written but has been killed twice mid-run before finishing all four bursts; no result yet. Full writeup of all of the above: `variability_analysis.md` and `BUGS.md` BUG-23.

`lorentz_factor.py` currently adopts each episode's **duration** as $t_\text{v}$ — an upper bound on the true variability timescale, which makes $\Gamma_\text{min}$ conservative by construction. Bukhari et al. (2022) instead *fit* the variability timescale from the light curve using the Norris et al. (2005) pulse profile,

$$I(t) = F \exp\!\left[2(\tau_1/\tau_2)^{1/2}\right]\exp\!\left[-\tau_1/(t-t_s) - (t-t_s)/\tau_2\right],\qquad t>t_s$$

fitted to a single pulse of the 10–400 keV light curve (their eq. 9; they obtain $t_\text{v}=0.7\pm0.02$ s for GRB110721A).

- **Data is in hand:** `light_curves/` has background-subtracted, energy-resolved RMFIT `.dat` files for all four bursts, a `lightcurve_data()` reader and per-GRB drivers. All are 64 ms binned with peak rates of 2 900–44 300 s$^{-1}$ — ample to constrain a sub-second pulse.
- **Expected effect:** $\Gamma_\text{min}\propto t_\text{v}^{-1/(2\alpha+2)}$, so replacing a 40 s duration with a $\sim$0.7 s fitted width raises TR3's limit from 137 to $\sim$256, closing part of the gap to Abdo et al.'s 887.
- **Prerequisite, and the real work:** the Norris profile describes *one* pulse. Several episodes are Bayesian-block intervals spanning multiple pulses, and T90 covers the entire burst. Pulse identification per episode, plus a stated rule for which pulse defines $t_\text{v}$ in a multi-pulse episode, must be settled first — these are judgement calls, not mechanical steps, and should be recorded in `lorentz_factor.md` when made.
- **Do not import numbers from Bukhari et al.** — they use Planck 2020 cosmology ($H_0=67.4$, $\Omega_M=0.315$), not this project's 69.6/0.286.
- A second use of the same fit, if pursued: fitting in several energy bands tests whether $R_\text{GeV}\simeq R_\text{MeV}$, as they do.

## Phase 6 (optional) — Outflow composition ($\eta$, $\sigma_0$) via hybrid-outflow inversion

*Discussed and scoped 2026-08-22, not started. Motivated by the user's "what are we building towards" question about the paper's throughline (see `HANDOFF.md` §11 for the short version) — this is the candidate headline result that would make the Lorentz-factor apparatus (Limit A, Limit B) the paper's actual point rather than a validation footnote to the BB-fraction census, matching the shape of Bukhari et al. (2022).*

**Headline method: Gao & Zhang (2015)** (`GRBResearchPaper/Literature Review/GaoHeZhang2015.pdf`), a "top-down" inversion that takes an episode's detected $(kT_\text{BB}, F_\text{BB})$ plus one more constraint and returns the central-engine pair $(\eta, \sigma_0)$ directly — composition, not just a bound. Their own worked example (GRB110721A) supplied the missing constraint as an assumed $r_0$; ours doesn't have to, because **Limit A/B's $\Gamma_\text{min}$ is already sitting in `lorentz_results.csv`/`lorentz_results_limit_b.csv`** and can serve as that second constraint instead. This reframes Limit A/B's narrative role in the paper: currently written up (§5) as an *independent cross-check* on the thermal $\Gamma$, it becomes instead a *required input* to the composition result — revisit that framing in `section-5-data-analysis.tex` when this phase lands, since the current wording ("an independent cross-check") will read as stale once Limit A/B are doing load-bearing work.

**Non-BB episodes (the majority, ~13 of ~26 across the sample): Zhang & Pe'er (2009) suppression bound** (`GRBResearchPaper/Literature Review/Zhang2009 - Evidence of an Initially Magnetically Dominated Outflow in GRB 080916C.pdf`), folded in as the *same* physical argument applied where Gao & Zhang cannot run — no detected BB means no $(kT,F_\text{BB})$ to invert, but the predicted pure-fireball photosphere flux vs. our actual non-detection threshold still yields a $\sigma$ lower limit. One paragraph in the paper, not a rival subsection: this is a regime of the same question, not a second method.

**Explicitly dropped from the computed-results path: Hascoët et al. (2013)**'s $Q = f_\text{th}/f_\text{IS}$ efficiency-ratio framework — a third, independent theoretical route to interpreting $f_\text{BB}$ in terms of $\sigma$ that would duplicate what Gao & Zhang (detected case) and Zhang & Pe'er (non-detected case) already cover between them. Decided explicitly 2026-08-22, before any code was written, because five candidate $\Gamma$/$\sigma$ frameworks (Limit A, Limit B, Gao & Zhang, Zhang & Pe'er, Hascoët) in one paper reads as a methods survey, not a paper with a result — the user's own concern, raised and acted on proactively. Hascoët may still earn a one-line citation as supporting context; it gets no equations, table, or figure of its own.

**Generalizing beyond GRB080916C:** both routes need $d_L$, same as every other $\Gamma$-dependent quantity in this paper — so the honest generalization is not "redshift-free," it's reusing the z-sweep-at-fiducial-$z=2$ convention `pe_er_photosphere.py` already applies to the three bursts without a spectroscopic redshift. `lorentz_factor.py`/`lorentz_factor_limit_b.py` currently gate on measured z (GRB080916C only); extending them to the same sweep would unlock this phase for all four bursts. Direct precedent: Bukhari et al. (2022) (`GRBResearchPaper/Literature Review/2022_06_06_Dr_Saeeda_Sajjad_Syed_Ali_Mohsin_Bukhari_Urooj_Murtaza.pdf`) computes $\Gamma_{\gamma\gamma}$ via the Hascoët (2012) multi-zone formalism at swept $z=1,2,3$ for GRB120709A, which has no measured redshift at all.

**Real complexity, not yet resolved — flag before implementing:**
- Gao & Zhang's inversion has six $r_\text{ph}$ regimes (their Table 3), selected by where $r_\text{ph}$ falls relative to $r_\text{ra}$ and $r_\text{c}$; regime selection needs an assumed acceleration power-law index $\delta$ — a judgement call to make and record, not a default to reach for.
- $\Gamma_\text{min}$ is a *lower bound*, not the true $\Gamma$ their equations expect as an input. Plugging it in directly needs care about which direction the resulting $\sigma_0$ bound points (upper vs. lower limit) — not yet worked out, and wrong without it.
- Yassine et al. (2017) (`GRBResearchPaper/Literature Review/Yassine2017 - Time Evolution of the Spectral Break in the High Energy Extra Component of GRB 090926A.pdf`) is precedent for discussing $\Gamma$ (and by extension $\eta,\sigma_0$) evolution *across episodes within a burst*, not just a single burst-averaged number — worth keeping in mind for how results get presented once computed.

**Status:** scoped only. No code, CSV, or `.tex` written for this phase yet.

Full assessment in `GRBResearchWork/codes-for-paper/lorentz_factor/lorentz_factor.md` §6b.

---

## Critical files

- `GRBResearchPaper/main.tex` — include toggles
- `GRBResearchPaper/tex_files/{abstract,section-5-data-analysis,section-6-discussion,section-7-conclusion}.tex`
- `GRBResearchWork/src/grb_research/grb_sed.py` (`SpectralModels._evaluate_components`) — component separation, reuse as-is
- `GRBResearchWork/src/grb_research/grb_calculations.py` — existing BB MC resamplers, reuse as-is; also `get_rng()` (lines 23-44) for RNG seeding
- `GRBResearchWork/src/grb_research/__init__.py` (`update_style()`, lines 58-105) and `grb_styles.py` (`GRBPlotStyle`) — required plotting conventions
- `GRBResearchWork/codes-for-paper/amati_relationship/amati_relationship.csv` / `.py` — template for the CSV column standard and a correct `update_style()`-based plotting script
- `GRBResearchWork/codes-for-paper/lorentz_factor/lorentz_factor.py` — pattern to follow for new modules' table/CSV generation; also needs its own GRB-sample/cosmology fix (Phase 0)
- `GRBResearchWork/LAT_analysis/` — the LAT source of truth: raw gtlike output per episode, `txt_to_csv.py` → `lat_photons.csv` → `csv_to_latex.py` → the appendix table, plus `LAT_analysis.md`
- `GRBResearchWork/results.json` — source of `amp_bb`/`kt_bb` per interval, and of the canonical episode labels
- `GRBResearchWork/codes-for-paper/best_names.txt` — defines the 13 BB-inclusive intervals to process
- `GRBResearchWork/ToDo.md` — checklist to update as phases complete
- `GRBResearchPaper/Literature Review/` — Pe'er 2006/2007 references for Phase 2 formulas