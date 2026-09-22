# Handoff — GRB photospheric-emission paper

**Written 2026-08-21, last substantively extended 2026-09-03 (through §13).** Picks up after Phases 0–3 of `PLAN.md`. Read `CLAUDE.md` first (scope, conventions, hard boundary), then this file.

**Staleness note (added 2026-09-16, commit refs updated 2026-09-22):** this file's narrative sections (§1–§13) stop at the RNG/seeding overhaul (2026-09-03). Both repos have moved on since — `GRBResearchWork` is now at `[main-minor-93]` (was `[main-minor-75]`/the `seed-correction` branch at the time §13 was written; `[main-minor-88]` as of when this note was first added, `[main-minor-91]` as of the previous update) and `GRBResearchPaper` is still at `[main-minor-21]` (was `[main-minor-15]`, unchanged since this note was added) — including at least a GBM-only 3-way refit exercise (`codes-for-paper/gbm_only_refit/`) and the Phase 5 variability-timescale work below, neither reconstructed here. Treat `git log` in each repo, `GRBResearchWork/PHASE5_TV_PLAN.md`, and `BUGS.md` (now 23 bugs + 11 observations + 1 plan note — all closed/resolved except OBS-11, DEFERRED, and BUG-23, fixed for GRB080916C only and confirmed live/unfixed for GRB131014A and GRB140206B as of 2026-09-22) as more current than this file's specifics; the status table in §1 has been spot-corrected for Phase 5 and the weakness-review pass (§12) only.

Everything below is **working-tree only. Nothing has been committed in either repo**, per the standing rule — **except §13** (RNG/seeding overhaul), where the user explicitly asked for commits partway through that session; see §13 for the exact commits.

---

## 1. Where the work stands

| Phase | Status |
|---|---|
| **Phase 0** — re-enable finished content, fix known inconsistencies | ✅ Complete |
| **Phase 1** — BB fractional flux $f_\text{BB}$ | ✅ Complete |
| **Phase 2** — photospheric radius & Lorentz factor (Pe'er 2007) | ✅ Complete |
| **Phase 3** — paper integration | ✅ Complete |
| **Phase 4** — final QA | 🟨 **In progress — start here** (see §6) |
| Unplanned — LAT data layer, units fix (BUG-18, OBS-08) | ✅ Complete |
| Unplanned — Lithwick & Sari Limit B cross-check (see §10) | ✅ Complete |
| Unplanned — external weakness-review pass (see §12) | ✅ Complete (all 10 items + Priority 3 pass, see §12) |
| Unplanned — RNG/seeding overhaul + paper integration (see §13) | ✅ Complete |
| Stretch: fireball magnetisation $\sigma_0$ | 🟨 Scoped as Phase 6 (see §11), not started |
| Optional: Norris-profile $t_v$ (see §6b) | 🟨 In progress since 2026-09-07 — see `GRBResearchWork/PHASE5_TV_PLAN.md` and `codes-for-paper/variability_analysis/variability_analysis.md`; not yet fed into `lorentz_factor.py`. GRB080916C reverted to a standard 7-pulse fit 2026-09-22 (was an 8-pulse `norris2-1`/`norris2-2` split, judged unsound); BUG-23 (wrong NaI detector silently picked after a resync) found and fixed there, confirmed live/unfixed for GRB131014A and GRB140206B |

The paper compiles clean: **0 errors, 0 undefined references, 24 pages** (was 18 at the start).

**Build command changed 2026-08-21, at the user's request: `pdflatex`+`bibtex` directly, not `latexmk`.** `latexmk` was this project's build tool up to that point (chosen for its automatic rerun-until-stable behavior); the user asked to switch after prior bad experience with it elsewhere. Manual sequence, run from `GRBResearchPaper/`:

```bash
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=out main.tex
BIBINPUTS=$(pwd): BSTINPUTS=$(pwd): bibtex out/main
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=out main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=out main.tex
```

**The `BIBINPUTS`/`BSTINPUTS` are not optional.** `bibtex` resolves `ref.bib` and `bibtex/aa.bst` relative to its own working directory, not `-output-directory`; running plain `bibtex out/main` from `GRBResearchPaper/` fails with "I couldn't open database file ref.bib" and then "I couldn't open style file bibtex/aa.bst", because both files live at the repo root while `bibtex` looks in `out/`. `latexmk` handled this path resolution automatically, which is part of why it was silently working before. Two `pdflatex` passes are needed after `bibtex` (not one) to settle both the bibliography and cross-references/page numbers.

---

## 2. How to run things (none of this is obvious)

- **`grb_research` is NOT installed** into `GRBResearchWork/.venv`. Every script needs `PYTHONPATH=src`:
  ```bash
  cd GRBResearchWork
  PYTHONPATH=src .venv/bin/python codes-for-paper/bb_fraction/bb_flux_fraction.py
  ```
- `codes-for-paper/amati_relationship/*` used to additionally import `from src.grb_research import ...`, needing `PYTHONPATH=.:src`. **Fixed 2026-08-21 (OBS-01)** — now uses `from grb_research import ...` like everything else, so plain `PYTHONPATH=src` is sufficient project-wide.
- **`ruff` is not installed** despite `ruff.toml` existing. Substitutes used so far: `python -m py_compile`, `awk 'length>120'` for the 120-char limit, and an AST scan for unused imports.
- `main.tex` uses `\documentclass[twocolumn, draft, a4paper]`. **`draft` stubs out figures**, so a clean draft build does *not* prove figure paths resolve. Strip `draft` to check them.
- **A clean build proves nothing about content.** A whole subsection sat commented out for most of this session while the build reported zero errors. Always `pdftotext out/main.pdf` and grep for the content you expect.

---

## 3. Results now in the paper

$f_\text{BB} = F_\text{BB}/F_\text{total}$, observer-frame 1 keV – 10 MeV, all 13 BB-inclusive episodes:

| GRB | $f_\text{BB}$ | note |
|---|---|---|
| GRB080916C | 0.050 – 0.068 | early episodes only |
| GRB131014A | 0.072 – 0.127 | all five episodes; TR2 is the outlier and the one `SBPL_BB` |
| GRB140206B | 0.058 – 0.063 | EX0 and TR1 only |
| GRB231129C | 0.206 – 0.247 | **most thermally dominated in the sample** |

Pe'er (2007) fireball parameters at $Y=1$ (GRB080916C at $z=4.35$, others at fiducial $z=2$):

| GRB | $\Gamma$ | $r_0$ [cm] | $r_\text{ph}$ [cm] |
|---|---|---|---|
| GRB080916C | 752 – 861 | $0.9$–$1.4\times10^{7}$ | $6.1\times10^{11}$ – $1.4\times10^{12}$ |
| GRB131014A | 458 – 714 | $0.7$–$4.7\times10^{8}$ | $3.4$–$7.6\times10^{12}$ |
| GRB140206B | 597 – 600 | $0.8$–$1.3\times10^{7}$ | $5.6$–$7.5\times10^{11}$ |
| GRB231129C | 385 – 426 | $6.2$–$8.9\times10^{8}$ | $3.5$–$4.8\times10^{12}$ |

$\Gamma_\text{min}$ (γγ opacity, GRB080916C only — the others lack redshifts): 260–507 across episodes, conservative by construction.

---

## 4. Three corrections that changed published numbers

All were found this session and all are worth understanding before touching the analysis.

**$E_\text{iso}$ was wrong twice, in opposite directions.**
`mc_e_iso_sampler` (a) integrated $d_L$ over redshift instead of evaluating it (×3.88 too large), and (b) evaluated the spectrum on the rest-frame grid while integrating over observed energies (×3.70 too small). These nearly cancelled, which is why the error survived so long. Both fixed. GRB080916C T90 is now $3.61\times10^{54}$ erg, and over Abdo's wider band the corrected value agrees with their published $8.8\times10^{54}$ to ~14% — it was an order of magnitude off before. All eight episodes now sit within $2\sigma$ of the Amati relation (was 3σ, and clustered low).

**$\Gamma_\text{min}$ went 197 → 258 → 134 → 507.** First because $f_1$ for SBPL was extrapolated from the 100 keV pivot with the *steep* index (wrong: SBPL follows the shallow index below its break). Then because $t_v$ was unified to each episode's own duration, replacing a literature 0.9 s that matched no episode and that §5's text did not describe.

**Then $\Gamma_\text{min}$ moved once more, by 3.4–4.5×, because the LAT photon energies are MeV and the code read them as keV** (BUG-18). This is the correction with the largest effect on what the paper *says*, not just what it tabulates: the thermal $\Gamma$ now exceeds the opacity floor by 1.5–2.4 rather than "five to ten". Full account in §9.

---

## 5. Open items

Full detail in `BUGS.md` (`GRBResearchWork/md_file_refs/`) — as of 2026-08-31 (when this section was last written), 18 bugs, 10 observations, 1 plan note, all fixed/closed/resolved per the user's decisions on 2026-08-21 (OBS-01–09) and 2026-08-31 (OBS-10, see §12). **Stale count as of 2026-09-22: `BUGS.md` now has 22 bugs, 11 observations, 1 plan note — all closed/resolved except OBS-11 (DEFERRED, user's call, added 2026-09-22).** The bugs/observations added since 2026-08-31 (BUG-15–22, OBS-11) are not summarized below; read `BUGS.md` directly.

- **OBS-05** — $S_\text{obs}$ in the k-correction is computed from the fitted model, not a catalogue fluence, so the correction cancels to a single integral. **User confirmed this is the deliberate design**, not an open question — already stated in §5. Our model over GBM's band gives $1.67\times10^{-4}$ vs Abdo's published $\sim2.4\times10^{-4}$ erg cm$^{-2}$.
- **OBS-02** — re-verified and found already resolved by BUG-11: `s_obs = detector_flux * model.interval.duration` (line 435) restores the duration multiply, and the k-correction (436-437) divides same-unit fluxes before converting to erg once. The index in `BUGS.md` was simply never updated after BUG-11 landed.
- **OBS-01** — fixed. Every `from src.grb_research import ...` switched to `from grb_research import ...` across `butterfly_all.py`, `best_names.py`, `amati_helpers.py`, `grb_fluence.py`, `model_parameters/utils.py`, `kt_evolution_080916c.py`. All now run with `PYTHONPATH=src` alone (verified by symbol-presence check, not full script execution — these MC pipelines take minutes). `pip install -e .` still not done, considered out of scope for this fix.

**Known limitation, not a bug:** our $\Gamma_\text{min}$ is weaker than published treatments (our TR3 gives 504 where Abdo et al. get 887 for an overlapping interval). The residual gap is methodological — Lithwick & Sari Limit A is a simplified analytic form — and the paper presents ours as conservative limits. Note this gap was a factor of 7 until BUG-18 was fixed on 2026-08-21; most of it was a MeV/keV units error, not method.

---

## 6. What Phase 4 should do

1. ✅ **Full recompile, and grep the PDF for expected content.** Done repeatedly while fixing BUG-18: 0 errors, 0 undefined refs, 23 pages, with the corrected $\Gamma_\text{min}$ values, the `(MeV)` header and the reworded ratio sentences all confirmed in `pdftotext` output rather than assumed from a clean build.
2. ✅ **Lint the new modules.** `ruff` is still not installed; the §2 substitutes were run over `LAT_analysis/{txt_to_csv,csv_to_latex}.py` and `lorentz_factor.py` — all compile, 0 lines over 120, no unused imports. Two real findings were fixed: dead `dataclass`/`EpisodeTypes` imports in `lorentz_factor.py`, and a LaTeX caption that could not satisfy both the 120-char Python limit and the one-sentence-per-line LaTeX rule until it was refactored into a sentence list. **Still to do: the Phase 1/2 modules** (`bb_flux_fraction.py`, `pe_er_photosphere.py`) have not been linted.
3. ✅ **Spot-check table values against `results.json` by hand.** Done 2026-08-21, via four parallel per-GRB subagents (one per burst), each independently pulling raw fit parameters from `results.json`, checking transcription into `bb_flux_fraction.csv`/`pe_er_photosphere.csv` and both generated LaTeX tables, and hand-recomputing $f_\text{BB}$ and Pe'er (2007) eqs. (1)/(4)/(5) against the CSVs. All 13 BB-inclusive intervals across all four bursts checked; every result agreed within ≤3% (well inside MC-median-vs-point-estimate slop) with no transcription errors, sign flips, or unit errors. One non-bug finding worth knowing: `kt_bb_keV` means "raw fit value" in `bb_flux_fraction.csv` but "MC median" in `pe_er_photosphere.csv` — documented in both `.md` files and `BUGS.md` OBS-09. The LAT table was already checked separately (§9).
4. ⬜ Re-read §5–§7 end-to-end for narrative coherence; they were edited in pieces, and BUG-18 has since reworded four sentences across §5, §7 and the abstract.
5. ⬜ Leave both working trees uncommitted for review.

For the **$\sigma_0$ stretch goal**: Pe'er's pure-fireball equations cannot give magnetisation. Use `GRBResearchPaper/Literature Review/GaoHeZhang2015.pdf` (hybrid outflow with arbitrary entropy and magnetisation). `GRBResearchPaper/Literature Review/Zhang2009 - Evidence of an Initially Magnetically Dominated Outflow in GRB 080916C.pdf` is directly about our brightest burst. Note $\eta$ (baryon loading) and $r_0$ are **already delivered** — $\eta = \Gamma$ in the coasting phase — so that ToDo bullet is only partly outstanding.

---

## 6b. Assessed but not done: Norris-profile fit for $t_v$

Raised by the user 2026-08-21 from §6.5 of `GRBResearchPaper/Literature Review/2022_06_06_Dr_Saeeda_Sajjad_Syed_Ali_Mohsin_Bukhari_Urooj_Murtaza.pdf` (Bukhari et al. 2022 — the user is an author).

**Verdict: feasible, worthwhile, not required.** Full write-up in `codes-for-paper/lorentz_factor/lorentz_factor.md` §6b.

They derive $t_v$ by fitting one pulse of the 10–400 keV light curve with the Norris et al. (2005) profile, rather than adopting an interval duration. We currently use the episode duration — deliberately an *upper* bound, hence a conservative $\Gamma_\text{min}$.

**The data is already in the repo and is good.** `GRBResearchWork/light_curves/` holds background-subtracted, energy-resolved RMFIT `.dat` files for all four bursts, plus a reader (`lightcurve_data()`) and per-GRB drivers. All four are binned at **64 ms**, spanning $\sim$600 s, with peak rates 2 900–44 300 s$^{-1}$ — a $\sim$0.7 s pulse is resolved by $\sim$10 bins across the rise alone.

**Why it matters:** $t_v$ is the largest controllable systematic in $\Gamma_\text{min}$ — the three candidate conventions span a factor of two (134 / 209 / 258 for T90). Replacing an upper bound with a measurement would raise $\Gamma_\text{min}$ predictably; for TR3, $t_v$ going from 40.256 s to $\sim$0.7 s gives $137 \to \sim256$, closing part of the gap to Abdo et al.'s 887.

**Blockers to settle first:** the Norris profile describes a *single* pulse, but several of our episodes are Bayesian-block intervals containing several, and T90 spans the whole burst. Pulse identification per episode, and a stated rule for which pulse defines $t_v$, are prerequisites and are judgement calls. Also note Bukhari et al. use Planck cosmology (67.4/0.315), not this project's 69.6/0.286 — borrow the method, not the numbers.

---

## 7. Files created or changed

**New at the project root** (`~/Pictures/grb_research/`, outside both repos, so untracked by either):
- `BUGS.md` — the bug log; read this before trusting any number.
- `PLAN.md`, `HANDOFF.md` — the phased plan and this file.

**New in `GRBResearchWork`:**
- `codes-for-paper/bb_fraction/` — `bb_flux_fraction.py`, `csv_to_latex.py`, CSV, PNG/PDF, table, `bb_fraction.md`.
- `codes-for-paper/photospheric_radius/` — `pe_er_photosphere.py`, `csv_to_latex.py`, CSV, PNG/PDF, table, `photospheric_radius.md`.
- `codes-for-paper/lorentz_factor/lorentz_factor.md` + regenerated CSV/table.
- `codes-for-paper/lorentz_factor/lorentz_factor_limit_b.py` + `lorentz_results_limit_b.csv` + `lorentz_table_limit_b.tex`. Closes the Limit B cross-check, §10.
- `LAT_analysis/` — `lat_photons.csv` (the LAT source of truth, built by the rewritten `txt_to_csv.py`), `csv_to_latex.py`, `lat_info_table.tex`, `LAT_analysis.md`. Closes OBS-08.

**Modified in `GRBResearchWork`:** `src/grb_research/{__init__,grb_calculations,grb_seds,grb_styles}.py`, `codes-for-paper/{amati_relationship,lorentz_factor}/*`, `ToDo.md`.

Two helpers were promoted into the package so Phases 1 and 2 share one implementation: `component_energy_fluxes()` and `draw_model_samples()` in `grb_calculations.py`. A third, `compute_tau_hat()` in `lorentz_factor.py`, is shared between Limit A and Limit B for the same reason (§10).

**Modified in `GRBResearchPaper`:** `main.tex`, `abstract.tex`, `section-{4,5,6,7}`, `ref.bib` (added `Peer2007`), plus new `tex_files/generated/` (six auto-generated tables, including `lat_info_table.tex` and `lorentz_table_limit_b.tex`) and two new figures in `images/section5/`. (At the time this was written, `main.tex` reached `lat_info_table.tex` via a two-line wrapper, `appendices/appendix_LAT_info.tex` — removed 2026-09-06 once every generated table's `dest` was already a stable path in its own right; `main.tex` now `\input{tex_files/generated/lat_info_table}` directly. See `CLAUDE.md`'s "Generated LaTeX tables".)

**Seven unreferenced images were deleted** from `GRBResearchPaper/images/` (all git-tracked, so recoverable).

---

## 8. Method notes — read these before changing any calculation

Each folder of analysis code carries a `<topic>.md` recording *why* each choice was made, who made it, what was rejected, and what was validated. They are the reasoning record the paper cannot hold:

- `codes-for-paper/bb_fraction/bb_fraction.md`
- `codes-for-paper/photospheric_radius/photospheric_radius.md`
- `codes-for-paper/lorentz_factor/lorentz_factor.md` — §4 (the $t_v$ decision) and §5 (statistical vs systematic budget) matter most.
- `LAT_analysis/LAT_analysis.md` — the data-plumbing layer, not physics: units, the episode-label join, and the one input that cannot be derived.

This is a required convention, not a nicety — see `CLAUDE.md`.

---

## 9. The LAT data layer, and a units bug that changed $\Gamma_\text{min}$

Unplanned work, done after Phases 0–3. Two `BUGS.md` entries: **BUG-18** (the units error) and **OBS-08** (now closed).

**The bug.** `compute_gamma_min` normalises by `E_max / 511.0`, where 511 keV is $m_ec^2$ — so it requires keV. The values it was given are **MeV**. Three independent confirmations: the source field is named `P > 0.9 Max (E) MeV`; the smallest values in the sample (117.7, 122.9, 194.5) sit just above the LAT selection floor of 100 MeV, and as keV would fall three orders of magnitude below the LAT band; and the user confirmed their LAT analysis is thresholded at >100 MeV. $\Gamma_\text{min}$ was low by $1000^{(\alpha-1)/(2\alpha+2)}$ = **3.4–4.5×**.

This overturned **BUG-16**, which had "corrected" the appendix header from MeV to keV on the strength of an arithmetic slip (27428.80 MeV is 27.43 GeV; 27428.80 *GeV* would be 27 TeV). That entry is now marked REVERTED.

**What it changed in the paper.** The headline consequence is not the numbers but the narrative: the thermal $\Gamma$ exceeded the opacity floor by "factors of five to ten", and now exceeds it by **1.5–2.4**. The two independent methods agree far more closely than the paper claimed — a stronger consistency check, not a weaker one. Reworded in `abstract.tex`, §5 (×2) and §7 (×2). The gap to Abdo et al. also narrows from a factor of 7 to 1.8, so most of what §5 attributed to Lithwick & Sari Limit A being a simplified analytic form was in fact this bug.

**The data layer (OBS-08).** `LAT_PHOTONS` was a hand-transcribed dict; both consumers now read one generated CSV:

```
LAT_analysis/<grb>/Ep*/  *_analysis_result_*.txt + *_fit_results_*.txt
        │  txt_to_csv.py          (rewritten — it now really does emit a CSV)
        ▼
    lat_photons.csv   (26 rows, the single source of truth)
        ├── csv_to_latex.py ─→ lat_info_table.tex ─→ tex_files/generated/lat_info_table.tex ─→ main.tex \input (direct, since 2026-09-06)
        └── lorentz_factor.py ─→ Gamma_min
```

Episode labels are resolved by matching interval bounds against `results.json`, never from directory names — `Ep5A` is `EX1`, not `EX4`. The hardcoded `LOW_SIGNIFICANCE` set is gone too; weak detections derive from `TS < 25`, which reproduced the old list exactly.

**Verified.** All 26 episodes cross-checked against source twice (by value, then by interval-bound mapping): zero mismatches, so the old transcription had been faithful — but faithful to the *rounded* appendix table. The CSV carries full precision (`301.204` vs `301.20`), which moved $\Gamma_\text{min}$ by at most $2.9\times10^{-6}$ relative and changed no rounded value.

**One input is still hand-supplied and cannot be derived:** the `\pFlux` footnote upper limits are 95% profile-likelihood limits from gtlike's `UpperLimits`, whose output is not in `LAT_analysis/`. They sit 1–8% above (Flux + 2σ), so they cannot be reconstructed from the fit results either. They live in `FLUX_UPPER_LIMITS` in `csv_to_latex.py`, and the generator **raises** if that dict disagrees with the derived `TS < 25` set — the failure mode is a crash, not a silent wrong number.

**Left for you to decide:** the four per-GRB `LAT_analysis/*/lat_analysis_<start>_<stop>.tex` files are superseded by `lat_info_table.tex` and nothing generates them any more. They are git-tracked, so deletion is recoverable; they were left in place rather than removed unasked.

---

## 10. Lithwick & Sari Limit B — an independent cross-check on $\Gamma_\text{min}$

Unplanned, done 2026-08-21 after cross-checking `lorentz_factor.py`'s equations against the Lithwick & Sari (2001) PDF itself (now in `GRBResearchPaper/Literature Review/`, at the user's request). Full account in `codes-for-paper/lorentz_factor/lorentz_factor.md` §8.

**What it is.** The paper derives two independent lower bounds on $\Gamma$: Limit A (photon annihilation — what the paper already used) and Limit B (Compton scattering off the $e^\pm$ pairs that annihilation itself creates). Limit B shares the same $\hat\tau$ as Limit A but does **not** depend on the LAT photon energy $E_\text{max}$ at all, so it remains usable for an episode whose highest-energy photon association is insecure — something Limit A cannot offer. The paper's own convention (Table 3 of Lithwick & Sari) is to report $\max(A,B)$.

**Implemented as a fully separate pipeline**, per your instruction: `codes-for-paper/lorentz_factor/lorentz_factor_limit_b.py`, its own `lorentz_results_limit_b.csv` and `lorentz_table_limit_b.tex` — nothing shared with Limit A's output files. `compute_tau_hat()` was factored out of `lorentz_factor.py` and imported by both, since $\hat\tau$ is algebraically identical between the two limits and duplicating it would risk exactly the kind of drift `BUGS.md` already logs for other formulas in this project.

**Verified against the source paper, not just cross-checked internally.** Every exponent and the numeric coefficient in both $\hat\tau$ (eq. 9) and Limit A (Table 2) were confirmed line-for-line against the Lithwick & Sari PDF before Limit B was implemented. $\hat\tau$ is bit-identical (max abs diff `0.0`) between the two CSVs across all 26 rows, as it must be for two functions computing the same closed form from the same seeded MC draws.

**Result: Limit A dominates in every episode of GRB080916C** (the only burst with a redshift, same restriction as Limit A) — so $\max(A,B)=A$ throughout, and the paper's adopted $\Gamma_\text{min}$ values are unchanged. This is a real finding, not a wasted computation: it confirms Limit A's dominance for this sample rather than assuming it.

**Integrated into the paper**, `tex_files/section-5-data-analysis.tex`, immediately after the existing Limit A paragraph: a numbered equation (`eq:gamma_min_limitB`, renders as eq. 8), citing the shared $\hat\tau$ via `\cref`, and a new generated table (Table 4, `\label{tab:lorentz_limit_b}`) via the same `% Generated by ... — do not hand-edit` \input convention as every other generated table. Rebuilt clean (0 errors, 0 undefined refs, 23→24 pages) and verified by `pdftotext`, not just a successful build: the equation renders correctly and Table 4's $\Gamma_{\min,B}$ values (112, 190, 213, 190, 115, 154, 131, 141) match the CSV exactly.

**Follow-up, same day: a comparison figure.** The user noticed `codes-for-paper/lorentz_factor/` was the only Phase 1/2 topic folder without a plot (`bb_fraction/` and `photospheric_radius/` both have one). Built `gamma_comparison_plot.py`: a single figure putting Limit A, Limit B, and the Phase 2 thermal $\Gamma$ (from `photospheric_radius/pe_er_photosphere.csv`) on one log-scaled axis, per episode of GRB080916C, with full MC error bars. It reads the three already-computed CSVs rather than recomputing anything — a fourth independent implementation of any of the three formulas was exactly the risk flagged in `CLAUDE.md`'s new "Verify, don't assume" section.

**Two real bugs caught by the user during review, not by the script running successfully:**
1. The in-axes "upper left" `Method` legend sat directly on top of the thermal-$\Gamma$ points — the highest values on the log axis — the same class of defect as `BUGS.md` BUG-15. Fixed by moving both legends outside the axes.
2. After that fix, the legend text was still clipped at the figure's right edge across three separate resize attempts, because `axis.add_artist(legend1)`-added legends aren't reliably included in matplotlib's automatic tight-bbox artist search. Fixed with an explicit `bbox_extra_artists=(legend1, legend2)` at `savefig` — verified by cropping and re-inspecting the actual saved PNG pixels, not by re-running the script and assuming success.

The user then asked to shrink the figure back down and shorten the "Method" legend's long labels with embedded newlines instead of growing the canvas further — both applied; final figure is `figsize=(8, 5.5)` with `subplots_adjust(right=0.62)`.

**Integrated into the paper** as Figure 9 (`fig:gamma_comparison`) in `section-5-data-analysis.tex`, with a cross-reference added to §6's existing thermal-vs-opacity comparison sentence, now also quoting the verified thermal-vs-Limit-B ratio (4.0–6.7×, computed from the CSVs, not eyeballed). **Verified without `draft` mode** — `draft` stubs figures, so a normal draft rebuild proves nothing about whether the image path resolves. Built to a throwaway `out_nodraft/` (removed after), confirmed 0 missing-figure warnings, `pdftotext`-grepped the actual caption text, and rendered page 12 to a PNG to visually confirm the figure and both legends display correctly at column width. `main.tex` and `out/` were left exactly as found (draft mode, normal `out/` rebuilt).

Full write-up, including the judgement calls (why episode markers/order come from live `TimeInterval` objects rather than re-parsing CSV strings, why the CSVs are read rather than recomputed), in `lorentz_factor.md` §8.6.

## 11. Paper throughline discussion, and Phase 6 scoping (2026-08-22)

User raised a standing concern across their own prior work-with-wife and this session: "what are we building towards" — is this paper a numbers dump (BB detection, photospheric radius, opacity limits, each reported independently) or does it build toward one physical claim, the way Bukhari et al. (2022) centered on the outflow/opacity result with everything else as support material?

Working conclusion: the candidate headline claim is **outflow composition** — whether GRB080916C's (and possibly the sample's) ejecta is a pure baryonic fireball or carries a significant Poynting-flux component — using the already-computed $\Gamma$/$f_\text{BB}$/$\Gamma_\text{min}$ machinery as ingredients rather than parallel results. Full literature synthesis (5 papers: Zhang & Pe'er 2009, Gao & Zhang 2015, Yassine et al. 2017, Hascoët et al. 2013, Bukhari et al. 2022) and the resulting plan are in `PLAN.md` Phase 6 — not duplicated here.

**Explicit scope decision, made before any code was written:** Hascoët et al. (2013) is dropped from the computed-results path (no equations/table/figure), kept at most as a citation. Reason: with Limit A, Limit B, Gao & Zhang, and Zhang & Pe'er already covering the detected/non-detected regimes between them, a fifth framework would make the paper read as a methods survey rather than a paper with one result — the user's own objection, raised proactively before scoping went further, not a retrospective fix.

Phase 6 is scoped, not started. Next session should read `PLAN.md` Phase 6 in full before touching this — it records real unresolved complexity (Gao & Zhang's six $r_\text{ph}$ regimes, and the direction $\Gamma_\text{min}$-as-lower-bound pushes the resulting $\sigma_0$ bound) that must be settled before writing code, not defaulted past.

## 12. External weakness-review pass (2026-08-31) — all 10 items + Priority 3 pass now CLOSED

A separate document, `grb_paper_weaknesses_and_fixes.md` (project root, outside both repos — not `PLAN.md`/`BUGS.md`), landed a 10-item referee-style review of the paper (plus a Priority 3 editorial-consistency list). Tracked issue-by-issue in a new file, `review-resolution.md` (project root) — **read that file, not this section, for the full status of all 10 items and the Priority 3 list.** This section writes up #2 and #4 in full (the two closed in the original session that changed a real published number); #6–#10 and Priority 3 closed later (2026-09-01, except #7 on 2026-09-06) and are only summarized briefly below — this section previously said they were still open, which was stale as of at least 2026-09-06.

**#1 (missing fitted-parameter appendix) was explicitly marked IGNORED at the user's instruction** — not evaluated, not fixed, still commented out in `main.tex`. Do not assume it's been looked at.

**#4 — statistical validity of the $\Delta$C-stat threshold.** The likelihood-ratio-test $\Delta\text{C-stat}\approx\chi^2(\Delta k)$ approximation used for BB-detection thresholds (28.74, 36.86) is formally invalid when a component's normalization sits at a parameter-space boundary under the null — the Protassov et al. (2002) result. **User decided:** cite the caveat, don't run a Monte Carlo null-distribution calibration. Added `Protassov2002` to `ref.bib` and two sentences to `section-4-joint-analysis-results.tex` §4.2 stating the caveat and leaning on Guiriec (2011)/Fana Dirirsa (2019) as precedent for the same thresholds. No numbers changed.

**#2 — hidden redshift-dependence in the cross-burst $f_\text{BB}$ comparison.** The paper's headline claim ("GRB231129C is 3-4x more thermally dominated than GRB080916C") compared $f_\text{BB}$ across bursts using a fixed **observer**-frame band (1 keV–10 MeV) — true redshift-independent *within* an episode, but not across bursts at different redshifts, since a fixed observer-frame band is a different rest-frame band per burst. GRB080916C has a measured $z=4.35$; GRB231129C (like the other two) has none.

**Fix:** recomputed $f_\text{BB}$ in a fixed **rest-frame** band (1 keV–10 MeV, matching the existing $S_\text{bol}$/$E_\text{iso}$ convention) alongside the existing observer-frame value. $z=4.35$ for GRB080916C; for the other three, **user chose a full sweep** ($z=0.5$–$5.0$, $z=2$ fiducial marked) over a single fixed fiducial value, mirroring the exact precedent `photospheric_radius/pe_er_photosphere.py` already set for the same problem (three of four bursts lacking a spectroscopic $z$). `bb_flux_fraction.py` rerun: 263 rows (was 13), ~10 minutes.

**Result — the actual finding, not just a re-derivation:**

| | Observer-frame $f_\text{BB}$ | Rest-frame $f_\text{BB}^\text{rest}$ |
|---|---|---|
| GRB080916C | 0.050–0.068 | 0.074–0.099 (at $z=4.35$) |
| GRB231129C | 0.206–0.247 | 0.209–0.248 (at fiducial $z=2$) |

The cross-burst **ranking is robust** — GRB231129C stays the most thermally dominated burst in both bands — but the headline **multiplier is not**: roughly 3–4× observer-frame, roughly 2–3× rest-frame. GRB080916C's higher redshift shifts its rest-frame band further from the fixed observer-frame band than GRB231129C's fiducial $z=2$ does, so its rest-frame $f_\text{BB}$ rises more. Capture-fraction sanity check (`bb_fraction_captured_in_band_rest`) stayed ≥99.9999% at every $z$ in the sweep, including $z=5$, so this isn't a band-truncation artifact — a genuine k-correction effect, not a computational one.

**User decided** (over keeping "three to four times" unchanged, or softening to one blended range) to report both bands' numbers explicitly everywhere the claim appears: `abstract.tex`, `section-5-data-analysis.tex` §5.2 (new definition paragraph, new figure `fig:bbfraction-rest-vs-z`, revised 3-sentence headline paragraph explaining the $z$-asymmetry mechanism), `section-6-discussion.tex` (shorter robustness-confirmation sentence — it never had the explicit "3-4x" wording to begin with), `section-7-conclusion.tex`. Table 4 (`tab:bbfraction`) gained $z$ and $f_\text{BB}^\text{rest}$ columns and bumped to `table*`.

**Also fixed in passing:** `bb_flux_fraction.py`'s RNG bypassed the shared `get_rng` helper despite `bb_fraction.md` documenting that convention — not a correctness bug (see `BUGS.md` OBS-10), but fixed while the file was open anyway.

**#3 — missing synchrotron/multi-break alternative model.** The paper's own description of the BB signature in GRB131014A and GRB231129C ("shallow shoulder on the low-energy wing," low-energy indices "close to the line of death") is exactly the signature Oganesyan et al. (2018), Ravasio et al. (2018), and Burgess et al. (2020) attribute to a marginally fast-cooling synchrotron spectrum with an added low-energy break, not necessarily a photosphere. The doc's fuller resolution asks for actually fitting that alternative (a 2SBPL model — the name traces to Ravasio et al. 2018) against 1-2 strong-BB episodes and comparing C-stat/BIC. **Scoping check: this is not something Claude can do in this environment.** `GRBResearchWork` has no raw GBM/LAT count-spectrum data (PHA/RSP) and no independent fitting pipeline — the whole codebase is downstream of RMFIT, which produces `results.json` and the per-episode `.fit` files (confirmed by inspecting one: FITS header says `FILETYPE = 'SPECTRAL FITS'`, "time-sequenced spectral fit parameters" — RMFIT's *output*, not a re-fittable input), and RMFIT itself (external IDL GUI software) isn't available to Claude. **User decided:** citation-only, same tradeoff as #4 — no new fits. Added a 5-sentence paragraph to the Introduction and a one-sentence caveat at the exact "shallow shoulder"/"line of death" spot in §6, citing all three papers (bibliographic details verified against arXiv/DOI via web search, not taken from memory, given the stakes of a wrong citation in a submitted paper). If the user later wants the real test, the division of labor would mirror everything else in this project: they run the RMFIT fit and drop the `.fit` file into the matching episode directory, Claude adds the model to `grb_enums.py`, parses it, and does the comparison.

**Verification:** bare `pdflatex` compile-error check passed (19 pages, was 18) — this is not content verification. The corresponding `VERIFY.md` entry ("Synchrotron/multi-break alternative caveat (weakness #3)", 2026-08-31) is now **CONFIRMED** (2026-09-03, source + build-artifact check) — updated here 2026-09-16 after finding this section still described it as open.

**#5 — asymmetric model grid (PL extension only tested with BB).** Checked `results.json` across all four bursts (26 episodes): confirmed zero bare BASE+PL (non-BB) fits exist anywhere, and `grb_enums.py` has `BAND_PL`/`CPL_PL`/`SBPL_PL` commented out. Initially treated this as blocked the same way as #3 (needs new RMFIT fits) and drafted a caveat along those lines — **the user corrected this twice, worth reading both corrections in full in `review-resolution.md` since they change the substance:**

1. First correction: a fully symmetric grid (BASE, BASE+PL, BASE+BB, BASE+PL+BB) *was* run once, in a much earlier iteration predating this repo, and dropped for complexity — not an oversight, and not something Claude needed to newly discover as "blocked."
2. Second, more important correction: the real reason BASE+PL isn't tested as its own branch in the current pipeline is that **BASE+PL and BASE+BB are not nested with respect to each other** — both add 2 dof to a BASE model, but along different, non-hierarchical branches, so a direct likelihood-ratio comparison between them isn't statistically well defined. BASE+BB→BASE+BB+PL, by contrast, *is* a valid nested comparison (2 more dof on top of an already-accepted BASE+BB), which is exactly the stepwise procedure §4.2 of the paper already describes (BASE→BASE+BB at Δk=2, BASE→BASE+BB+PL at Δk=4, then stepwise BASE+BB→BASE+BB+PL at Δk=2) — so the paper's existing model-selection procedure was already correct; only the one ambiguous framing sentence needed fixing.
3. Separately, the supplementary PL's physical origin traces to Fermitools' own power-law characterization of the LAT-band data (matches `LAT_analysis.md`'s documented `gtlike` fields `Index`/`Flux (0.1-100.0 GeV)`), not an independently RMFIT-tested continuum alternative.

**User explicitly asked to keep the paper text simple** — it states the Fermitools origin and the non-nested-comparison reason in two sentences, with none of the project history (the dropped earlier grid, etc.) in the manuscript. That history lives in `review-resolution.md` only. Implementation: `section-4-joint-analysis-results.tex` §4 opening paragraph, one sentence expanded to two. Verified with a bare `pdflatex` pass (19 pages, unchanged from before #5).

**Committed (#2, #4 only):** `GRBResearchWork` (`2fa38b3`, `70b8654`), `GRBResearchPaper` (`1f1c65f`) — see those repos' `git log` for exact diffs. `70b8654` is a fix-up for a mistake made mid-session: 14 pre-existing staged `.yaml~` files (unrelated to this work, present in the index before this session started) rode along into the first commit because `git add <my files> && git commit` was run without checking what else was already staged. Caught immediately via `git show --stat HEAD` and corrected with a second commit rather than an amend, per the project's git-safety rules. **#3's changes committed separately:** `GRBResearchPaper@875d1b1`. Note this commit was staged carefully (`ref.bib`, `section-1-introduction.tex`, `section-6-discussion.tex`, `out/main.pdf` only) after finding untracked root-level LaTeX build artifacts (`main.aux`, `main.pdf`, etc., not gitignored) that predated this commit and look like the user's own manual `VERIFY.md` check — left untouched, not committed, not investigated further.

**#6–#10 and the Priority 3 editorial-consistency list are also CLOSED** (2026-09-01, except #7 on 2026-09-06) — summarized briefly here since this section previously said they were still open; full detail in `review-resolution.md`, not repeated here:

- **#6** (GRB131014A off-axis angle) — GBM-only refit, all 5 episodes confirm the same BB-augmented model with $kT_\text{BB}$ agreeing within $1\sigma$ of the joint fit; written up in `section-6-discussion.tex` with two new figures (`fig:gbmonly-kt`, `fig:gbmonly-cstat`).
- **#7** (instrumental-artifact/low-energy threshold robustness) — a second GBM-only refit at a raised 40 keV NaI threshold confirmed BB detection is not a low-energy-band artefact; **user decided to keep this as an internal-only check, no paper text added** (a real but modest $kT_\text{BB}$ systematic was found alongside it, judged not to warrant a manuscript caveat).
- **#8** (GRB080916C vs. prior literature) — new paragraph in §6 reconciling this paper's BB finding and $\Gamma_\text{min}$ with Abdo (2009)/Guiriec (2015).
- **#9** (Amati-relation interpretation overstated) — softened wording plus new time-resolved-Amati-validity literature (Basak & Rao 2012, Ghirlanda et al. 2010) in §5.2/abstract.
- **#10** (N=4 sample overclaiming) — softened "universal"/"confirms...population-level" framing in three places (intro, abstract, discussion/conclusion).
- **Priority 3** (8-item editorial pass, including the "25 units" vs. 28.74 $\Delta$C-stat inconsistency and the other items originally logged as pending here) — all CLOSED 2026-09-01, backups in `_backups/2026-09-01_priority3_editorial/` in each repo.

A citation bug found while working #8 (wrong `Guiriec2015` cite for the GRB120323A sentence, should have been `Guiriec2013`) was also fixed — see `review-resolution.md`.

## 13. RNG/seeding overhaul, and paper integration (2026-09-02/03)

Unplanned, separate session from §1–12. Full detail in two new project-root files, `SEED_PLAN.md` (the plan) and `SEED_PLAN-implementation.md` (the implementation log, updated as each piece landed, not written from any agent's self-report) — this section summarizes, doesn't duplicate. `GRBResearchWork/src/grb_research/SEEDING.md` is the method-note equivalent for this thread, same convention as `bb_fraction.md` etc.

**The problem.** A paper reviewer flagged that Monte Carlo draws across the codebase reuse the same literal seed (`12345`, sometimes `1234`/`42`) hardcoded independently in many scripts. Investigation found it was worse than a style concern: several scripts reseeded an identical `SEED` constant on *every iteration* of a loop over models/GRBs/episodes (`lorentz_factor.py`, `amati_relationship.py` among them) — a fresh, identically-seeded generator per call means those iterations drew literally identical underlying random numbers, a real correctness bug (silently correlated/duplicated MC draws across iterations meant to be independent), already logged as `BUGS.md` OBS-10 for the narrower `bb_flux_fraction.py` case before this overhaul generalized the fix everywhere.

**The fix.** `MASTER_SEED = 2828702241`, drawn via `secrets.randbits(32)` (OS CSPRNG, not `np.random`, so it's independently auditable) on 2026-09-02T15:40:46 UTC — provenance reproduced verbatim in `SEEDING.md`. `seed_from_name(name, master_seed)` derives a deterministic per-script seed via SHA-256 of the script's basename + `MASTER_SEED`, so every script gets a distinct, reproducible seed without hardcoding one. Each script now builds **one** `rng` at the top (`SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`) and threads it through every downstream call — the actual fix for the reseed-per-iteration bug, since one `Generator`'s state advances across calls instead of resetting.

**Two dispatch phases, both independently re-verified — not trusted from either background agent's self-report.** Phase A (`src/grb_research/`): added `MASTER_SEED`/`seed_from_name`, routed 4 bypass sites (`ModelResampler`, `FluxFluenceCalculator`, `ParameterSet.get_populated_values`, `plot_covariance_corner`) through `get_rng`. A predicted circular import (`grb_atomic → grb_calculations → grb_model → grb_atomic`) was confirmed real by testing a module-level import first (it failed), fixed with a local import instead. Phase B (`codes-for-paper/`): replaced literal seeds with `seed_from_name(__file__)` across 10 scripts (`grb_fluence.py`, `lorentz_factor.py`, `lorentz_factor_limit_b.py` — deliberately given its *own* independent seed, decorrelated from Limit A — `gbm_only_refit.py`, `bb_flux_fraction.py`, `pe_er_photosphere.py`, `epeak_vs_kt.py`, `peak_energy.py`, `amati_relationship.py`, `butterfly_all.py`; `variability_analysis/norris..py` updated too but dormant, no caller). Verification: every script rerun twice from a clean state, output CSVs diffed byte-for-byte, `seed` column confirmed single-valued — all deterministic. The `tau_hat` bit-identical cross-check between Limit A and Limit B (§10 above) was re-confirmed to still hold after decorrelating their seeds (max abs diff `0.0` across all 26 rows), since that identity depends on shared closed-form math on point values, not a shared seed.

**Follow-up cleanup, same threads.** Once every real call site was confirmed to pass `rng=` rather than `seed=`, the now-vestigial `seed=` parameter was stripped entirely from 11 downstream functions (`draw_model_samples`, `mc_e_iso_sampler`, `mc_spectra_sampler`, `plot_all_models`, `ModelResampler.__init__`, `FluxFluenceCalculator.__init__`, `ParameterSet.get_populated_values`, `plot_covariance_corner`, `convert_sbpl_to_band`, `plot_grbs_over_amati_relationship`, `plot_unknown_redshift_grb`), making `rng` required everywhere. Proven a no-op by construction (`get_rng` already returned `rng` unchanged whenever given one) and confirmed empirically: 8 of the 10 Phase B scripts came back byte-identical against their pre-cleanup CSVs on rerun; the 2 slow ones (`bb_flux_fraction.py`, `pe_er_photosphere.py`, ~10 min/run at full sample count) got a clean smoke test at a temporarily reduced `N_SAMPLES=10` instead — the new `CLAUDE.md` "Determinism checks on slow MC scripts" convention, added at the user's request specifically because of this.

**Paper integration (`GRBResearchPaper`).** New appendix section "Monte Carlo seeds" (`tex_files/generated/seed_table.tex`, reached via a wrapper file `appendices/appendix_seeding_table.tex` at the time this was written — that wrapper was removed 2026-09-06, `main.tex` now `\input{tex_files/generated/seed_table}` directly), a 2-sentence summary in `section-5-data-analysis.tex` after the existing cosmology sentence, and a "Monte Carlo seed: N" caption sentence on the 7 figures whose plotted quantities are actually MC-derived (`peak_energy_best__all`, `amati_relationship`, `bb_flux_fraction`, `bb_flux_fraction_rest_vs_z`, `pe_er_photosphere`, `butterfly_all`, `gamma_comparison`). The two `gbm_only_refit` figures were deliberately excluded after checking the plotting code, not the caption wording: `kt_bb_err_keV` is a raw fit-covariance error and $\Delta$C-stat is a deterministic fit-statistic difference, so neither actually plots an MC quantity even though the producing script does MC elsewhere. User explicitly declined converting the mixed-provenance `peak_energy_best__all` (partly MC, partly direct fit error) or `low_index`/`high_index` (no MC at all) to "fully MC" for consistency — would trade an exact number for a seed-dependent approximation with no accuracy gain, a real methodology change, not a caption fix.

**Registry + sync tooling, added at the user's request.** `codes-for-paper/seed_registry.yaml`, `figure_registry.yaml`, `table_registry.yaml` are the source-of-truth mapping from script/asset to seed/destination, each entry categorized `active`/`unused`/`dormant` (or `no_known_source`) so nothing is ever silently dropped from view — new general convention, documented in `SEEDING.md`. `seed_table_to_latex.py` generates the appendix table from the registry, computing every seed live rather than hardcoding it. `sync_paper_assets.py` copies every registered figure/table from this repo into `GRBResearchPaper` in one pass; its first real run found several files (`amati_relationship.png`, `bb_flux_fraction.png`, `lorentz_table.tex`, and others) that had been regenerated during this session but never actually copied into the paper repo — a real gap the tool caught, not just a validation exercise.

**Committed, at the user's explicit request** (the one exception to this project's working-tree-only rule): `GRBResearchWork` on branch `seed-correction`, three commits — `ab854ad [seed-correction-minor-1]` (Phase B script/data changes), `1eaad7b [seed-correction-major-1]` (registry/sync tooling — this branch's first major, its counter separate from `main`'s major-1..5), `375fea6 [seed-correction-minor-2]` (`src/grb_research/` core changes + `SEEDING.md`, the latter a deliberate tracked exception to the usual untracked-`.md` convention since it's a code-root reference doc, not a working document). `GRBResearchPaper` on `main`, one commit, `f80338e [main-minor-15]`, bundling the appendix/caption/paragraph changes and the resync of previously-stale figures/tables.

**Verification status.** All code-level claims above (byte-identical reruns, the `tau_hat` cross-check, the no-op proof) were checked directly, not assumed. Paper *content* is different: a bare `pdflatex` compile check passed (0 errors, 0 undefined references, 28 pages) after fixing two real formatting bugs (a `table` environment overflowing single-column width, then overflowing worse at `table*` width because plain columns don't wrap long text — fixed with `\resizebox`). That proves the LaTeX is syntactically correct, not that the new sentences read well or the table is legible at its rendered size — the corresponding `VERIFY.md` entry ("RNG-seeding paper integration: paragraph, caption seeds, appendix table", 2026-09-03) is now **CONFIRMED** (2026-09-04, source + user visual confirmation) — updated here 2026-09-16 after finding this section still described it as PENDING.

**Supersedes `BUGS.md` OBS-10** (`bb_flux_fraction.py`'s RNG bypassing `get_rng`) — that fix is now folded into the general scheme; `compute_fraction` no longer has a `seed` parameter to bypass anything with.
