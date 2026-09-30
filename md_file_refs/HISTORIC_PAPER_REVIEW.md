# Historic archive — paper review & weaknesses

Consolidated 2026-10-01 from six retired standalone files, each preserved **verbatim**
(including raw, unedited reviewer notes) under its own section below, in reading/genealogy
order: a first raw review pass → the systematized task list it fed → the resolution log
tracking that list → a later separate review round → two incremental closure checklists.
This is archival/historic reference material, not living project memory — for current
conventions and status see `CLAUDE.md`, `BUGS.md`, `PLAN.md`, `HANDOFF.md`, `VERIFY.md` in
this same directory. See `HISTORIC_INDEX.md` for a one-line manifest of where every retired
file's content now lives.

---

## quick-fixes-high-priority.md

1. Stepwise selection GRB231129C ka headline claim kamzor kar sakta hai. Ye sab se serious hai. Aap base model pehle select karte hain, phir BB add karte hain. 231129C ke liye base SBPL (895.36) < Band (905.96), to SBPL+BB chuna gaya = 837.19. Lekin table 12 mein Band+BB = 818.24 hai — 19 C-stat units behtar, same 6 parameters, aur parameters bhi SAFE hain. Yehi TR1 (679.30 vs 688.35) aur EX0 (658.54 vs 666.15) mein bhi hai.

Masla ye hai ke BB normalization dono mein alag hai (T90: 190.4 vs 314.1), yaani f_BB badal jayega — aur f_BB ≈ 0.21–0.25 hi wo number hai jispe "3–4× more thermal" wala headline khada hai. Referee ye pakar lega. Do options:

Band+BB ke saath f_BB recompute karke dikha dein ke ranking survive karti hai (best case), ya
§5.2 mein stepwise ordering ko explicitly justify karein aur maan lein ke global minimum nahi choona ja raha.

Filhaal §5.2 mein "proceeds stepwise" sirf BB→PL step ke liye illustrate hua hai, BASE→BASE+BB ke liye nahi.

2. Unphysical kT fits tables mein hain, rejection criterion likha nahi. Table 11 Episode IV: SBPL+BB kT = 5.11 keV, amplitude 8782; Band+BB kT = 4.97. Table 10 EX–B: SBPL+BB kT = 9.67 keV, amplitude 5382. Ye ΔC-stat ~30 clear karte hain aur errors bhi chhote hain — yaani aapke likhe hue criteria ke hisaab se ye BEST hone chahiye the. Aap ne (sahi tor par) reject kiya, lekin kis basis par, ye paper mein kahin nahi hai. Ek line chahiye: "Fits with kT below ~15 keV, i.e. within a factor of a few of the NaI lower bound, are rejected as unphysical runaways at the edge of the band."

3. NaI low-energy artifact test abhi bhi missing. GBM-only refit ne LAT wala concern band kar diya, lekin asal instrumental risk NaI side hai — iodine K-edge 33.2 keV par, aur aapki kT range 27–59 keV bilkul wahin baith rahi hai. Aur point 2 ke runaway fits isi ki taraf ishara kar rahe hain. Sasta test: NaI lower bound 10 keV se badha kar ~40 keV karke refit, aur dikhayein ke kT aur ΔC survive karte hain. Ye figure 13/14 wale format mein hi ho jayega.

4. r_ph ka order-of-magnitude offset quantify nahi hua. §6.6 kehta hai z-dependence "mild" hai, phir 3–8×10¹² vs 6×10¹¹ ke offset ko "assumed redshift" par daal deta hai. Ye do baatein aapas mein nahi milti — agar dependence mild hai to wo ek order of magnitude explain nahi kar sakti. Figure 12 se numbers already available hain: z = 0.5 aur z = 5 par 131014A/231129C ka r_ph likh dein. Agar offset z=0.5 par bhi rehta hai to seedha likhein ke ye redshift ka artifact nahi hai — wo zyada honest aur zyada interesting hai.

5. Synchrotron alternative Conclusions ki future-work list mein nahi hai. Aap ne intro aur §7.1 mein acknowledge kiya, lekin §8 ki teen-item list (Y, 140206B, redshifts) mein 2SBPL nahi hai. Referee ke liye ye ajeeb lagega — jo cheez aap khud "genuine competing explanation" keh rahe hain, wo unresolved-issues list se ghayab hai. Item 4 bana dein.

---

## grb_paper_weaknesses_and_fixes.md

# GRB Paper — Weak Aspects & Resolution Plan

**Paper:** Thermal and Non-Thermal Emission in Bright Fermi Gamma-ray Bursts
**Authors:** Iqra Siddique, Syed Ali Mohsin Bukhari
**Purpose of this file:** Structured task list for an AI agent (or human collaborator) to work through the identified weaknesses, in priority order. Each item has: the problem, why it matters, and a concrete resolution path. Where relevant, exact source lines/sections are cited so the agent can locate them without re-reading the whole paper.

---

## Status note (already discovered)

The appendix containing the fitted-parameter tables (`table_080916C_large.tex`, `table_131014A_large.tex`, `table_140206B_large.tex`, `table_231129C_large.tex`, `appendix_LAT_info.tex`) and the `\label{sec:grb-parameters}` anchor are **present in the LaTeX source but commented out** in `main.tex`. This is the direct cause of the unresolved `table ??` / `appendix ??` references in the compiled PDF. This explains Weakness #1 below and should be the first thing checked/fixed, since it may already resolve most of that item with no new writing required.

**Immediate action:** uncomment the `\begin{appendix}...\end{appendix}` block in `main.tex`, verify the four `table_*_large.tex` files are actually populated (not placeholders), recompile (2 passes), and confirm every in-text `\ref{}` now resolves.

---

## Priority 1 — Blocking issues (must fix before any submission)

### 1. Missing fitted-parameter appendix (partially a compile bug, partially a content gap)
- **Status: IGNORED (2026-08-31)** — out of scope for this pass per user instruction. See `review-resolution.md` (§review-resolution.md in this archive).
- **Problem:** Body text references `appendix ??` and `table ??` in multiple places (Sections 4, 5, 6). No table anywhere in the compiled paper gives per-episode fitted values (α, β, E_p, kT, amplitude, C-stat, dof, ΔC-stat) for each of the four continuum models.
- **Why it matters:** A spectral-fitting paper is unverifiable without these numbers. This is the first thing a referee will ask for.
- **Resolution:**
  1. Uncomment the appendix block in `main.tex` (see status note above).
  2. Verify each `table_*_large.tex` contains: model name, episode, all free parameters with 1σ errors, C-stat, dof, ΔC-stat vs. the relevant baseline model, and the SAFE/MARGINAL/UNSAFE flag (the red-highlight-with-% overshoot scheme described in the commented intro paragraph).
  3. Recompile twice to resolve cross-references.
  4. Spot-check that values in the appendix match the summary numbers already quoted in Table 5 and Table 6 of the main text (e.g. GRB080916C T90 kT = 44.45 ± 2.05 in Table 5 vs. 44.46⁺²·⁰¹₋₂.₀₂ in Table 6 — resolve this discrepancy, they should be identical or the difference should be explained, e.g. different rounding pipelines).

---

### 2. Hidden redshift-dependence in the cross-burst f_BB comparison
- **Status: CLOSED (2026-08-31)** — rest-frame f_BB computed (z-swept for the 3 unmeasured-z bursts), both bands reported everywhere the claim appears; ranking robust, multiplier revised 3-4x → also reports ~2-3x rest-frame. See `review-resolution.md` (§review-resolution.md in this archive).
- **Problem:** f_BB is computed in a fixed **observer-frame** band (1 keV–10 MeV, Eq. 9), and described as "independent of redshift." That is true for a single episode, but it does **not** make f_BB comparable *across* bursts at different redshifts, because a fixed observer-frame band corresponds to a different rest-frame band for each z. GRB080916C (z = 4.35) vs. the other three bursts (assumed z = 2) are being compared in genuinely different rest-frame energy windows relative to where the BB peak (≈3.92 kT) sits.
- **Why it matters:** This directly undermines the headline claim that "GRB231129C is 3–4× more thermally dominated than GRB080916C" — some or all of that difference could be a k-correction/band-selection artifact rather than a physical difference in thermal strength.
- **Resolution:**
  1. Recompute f_BB in a fixed **rest-frame** band (e.g. 1 keV–10 MeV rest-frame) for all four bursts, using each burst's actual/assumed z.
  2. Report both the observer-frame and rest-frame f_BB values side by side.
  3. Quantify how much the burst ranking (which is "most thermal") changes between the two band choices.
  4. If the ranking is robust to this choice, state that explicitly as a robustness check — it strengthens the claim. If it isn't robust, revise the discussion/conclusion language accordingly.

---

### 3. Missing synchrotron / multi-break alternative model
- **Status: CLOSED (2026-08-31)** — citation-only. Fitting a 2SBPL/synchrotron-break model against real data requires RMFIT, which Claude cannot run in this environment (no raw PHA/RSP data, no independent fitting pipeline in this repo) — user chose citation-only over waiting for new RMFIT fits. See `review-resolution.md` (§review-resolution.md in this archive).
- **Problem:** The paper's own description of the BB signature — "shallow shoulder on the low-energy wing," "low-energy indices close to the line of death" — is exactly the signature attributed in recent literature (Oganesyan et al. 2018–2019; Ravasio et al. 2018–2019; Burgess et al. 2020) to a **marginally fast-cooling synchrotron spectrum with a second low-energy break**, not necessarily a photosphere. This is currently not mentioned or tested against.
- **Why it matters:** This is the most likely first substantive objection from any referee working in GRB prompt-emission spectroscopy. Leaving it unaddressed risks the paper being seen as not engaging with the current state of the field.
- **Resolution:**
  1. Add a citation and one paragraph in the Introduction and/or Discussion acknowledging the synchrotron-break alternative explanation for low-energy spectral softening.
  2. Select at least one or two representative episodes with strong BB significance (e.g. GRB131014A TR2, or GRB231129C T90) and fit a double-smoothly-broken-power-law (2SBPL) or equivalent synchrotron-motivated model with an extra break as a competing model.
  3. Compare C-stat / BIC between BASE+BB and the synchrotron-break alternative for these episodes.
  4. Report the outcome honestly — even a partial test ("BB preferred by ΔC-stat = X in the cases checked") substantially strengthens the paper's credibility, regardless of which model wins.

---

### 4. Statistical validity of the ΔC-stat threshold for BB detection
- **Status: CLOSED (2026-08-31)** — resolved via citation only (step 4 of the resolution below), no MC calibration. See `review-resolution.md` (§review-resolution.md in this archive).
- **Problem:** The likelihood-ratio-test approximation ΔC-stat ≈ χ²(Δk) is used to justify the BB-detection thresholds (28.74, 36.86). This approximation is formally invalid when a new model component's normalization sits at a boundary of parameter space (BB flux → 0 under the null) and an additional parameter (kT) is unconstrained under the null — the classic Protassov et al. (2002) result. Since this test underlies every BB detection claim in the paper, it is a foundational statistical issue, not a minor one.
- **Why it matters:** If the actual null distribution of ΔC-stat is broader than χ²(2)/χ²(4), the true significance of some "BEST" BB detections could be lower than claimed.
- **Resolution:**
  1. Cite Protassov et al. (2002) and acknowledge the caveat explicitly in Section 5.2.
  2. Ideally, run a Monte Carlo calibration for at least one representative burst/episode: simulate ~1000 synthetic spectra from the best-fit non-BB model, refit both BASE and BASE+BB to each, and build the empirical ΔC-stat null distribution.
  3. Compare the empirical 5σ-equivalent threshold to the analytically assumed 28.74/36.86 and report whether the adopted thresholds are conservative, matched, or need revision.
  4. If a full simulation is infeasible before submission, at minimum state the assumption is conservative and cite prior work (e.g. Guiriec et al. 2011, Fana Dirirsa et al. 2019) that used similar thresholds successfully, to justify not deviating from established practice.

---

## Priority 2 — Substantive but not blocking

### 5. Asymmetric model grid (PL extension only tested with BB)
- **Status: CLOSED (2026-08-31, re-confirmed 2026-09-01)** — clarification only, not actually a gap: BASE+PL and BASE+BB are non-nested (not directly LRT-comparable), so the paper tests the valid nested chain BASE→BASE+BB→BASE+BB+PL instead; the supplementary PL also traces to Fermitools' LAT-band characterization, not an independently RMFIT-tested model. See `review-resolution.md` (§review-resolution.md in this archive) for full detail (kept out of the paper text at the user's request).
- **Problem:** Section 4 states the supplementary power-law (for high-energy LAT emission) was added "only for BASE+BB models," not for the plain BASE models.
- **Why it matters:** Without testing BASE+PL as its own model, the analysis cannot cleanly distinguish "this episode needs a BB" from "this episode needs *any* extra spectral component." This could inflate apparent BB significance.
- **Resolution:**
  1. Re-run the model grid symmetrically: BASE, BASE+PL, BASE+BB, BASE+PL+BB for every episode.
  2. Re-derive the BEST-model classification with the full grid.
  3. Report whether any episode previously classified as BASE+BB is actually better described by BASE+PL alone.

### 6. GRB131014A off-axis angle concern
- **Problem:** Section 4.2.2 states the source off-axis angle was 70°, exceeding the nominal LAT field-of-view boundary for the entire T90 interval. This is the same burst for which every single time-resolved episode shows a BB-augmented BEST model.
- **Why it matters:** At large off-axis angles, LAT effective area drops sharply and response-matrix systematics increase, which could produce spurious low-energy curvature that mimics a BB component.
- **Resolution:**
  1. Add explicit discussion connecting the off-axis geometry to the BB detections for this burst.
  2. Re-fit GRB131014A using GBM data only (drop LAT) and check whether the BB component survives without the LAT contribution.
  3. Report the GBM-only vs. joint-fit comparison as a robustness check in the results or discussion section.

### 7. No instrumental-artifact / low-energy threshold robustness test
- **Problem:** All four bursts return BB temperatures in a narrow, suspiciously consistent range (kT ≈ 27–59 keV). NaI response systematics near the low-energy end (10–30 keV) and the iodine K-edge (~33 keV) are known sources of spurious curvature in this exact energy range.
- **Why it matters:** A cross-burst consistency this tight, if instrumental, would undermine the entire thermal-component narrative.
- **Resolution:**
  1. Re-fit a representative subset of episodes with the NaI low-energy threshold raised from 10 keV to 20–25 keV.
  2. Report whether BB detections (both significance and kT) survive this change.
  3. Include this as a short "systematics" paragraph or subsection.

### 8. GRB080916C thermal detection vs. prior literature (Abdo et al. 2009, Guiriec et al. 2015)
- **Problem:** GRB080916C is a well-known case in the literature for being well-described by a pure Band function, with prior work finding at most a marginal thermal component. This paper claims clear BB detections in T90, EX0, and TR1 for this burst, with no direct discussion of why the result differs.
- **Why it matters:** Silently disagreeing with established literature on a well-studied burst invites referee skepticism unless the discrepancy is explained.
- **Resolution:**
  1. Add a paragraph directly comparing this paper's GRB080916C fits to Abdo et al. (2009) and Guiriec et al. (2015): different response files (GBMRSP v2.0), different episode boundaries, inclusion of EX windows, or different ΔC-stat threshold could all explain the difference — identify which.
  2. Separately, verify the Abdo et al. (2009) published Γ_min value for this burst (paper implies ~890 but this should be confirmed from the source) and directly compare it to this paper's thermal Γ ≈ 750–860. If the thermal Γ is *below* a previously published Γ_min, this is a real tension that needs discussion — e.g. via Y > 1 (since Γ ∝ Y^(1/4), Y ~ 1.5–2 could resolve it) or measurement differences in the two Γ_min derivations.

### 9. Amati-relation interpretation is overstated
- **Problem:** Section 6.3 states that episodes lying within the 2σ Amati band "confirms the physical self-consistency of the spectral fits." Given the intrinsic scatter of the Amati relation, landing within 2σ is a fairly weak test. Additionally, applying the (originally time-integrated) Amati relation to time-resolved episodes is itself a contested practice in the literature.
- **Why it matters:** Overstating what this check demonstrates weakens the paper's overall rigor in referees' eyes.
- **Resolution:**
  1. Soften the language — describe this as "consistent with" rather than "confirms."
  2. Cite prior work specifically validating (or debating) time-resolved Amati relations (e.g. Ghirlanda et al., Basak & Rao) to justify the approach.

### 10. Overclaiming population-level generality from N = 4 non-complete sample
- **Problem:** The Introduction explicitly frames the goal as testing "whether photospheric emission is a universal or selective feature," and the Conclusion uses similarly broad language, despite the sample being explicitly stated elsewhere as a non-flux-limited, brightness-selected case study.
- **Why it matters:** Internal inconsistency in framing; referees will flag the mismatch between the modest N=4 case-study design and the universalizing claims.
- **Resolution:** Do a consistency pass on Abstract, Introduction, and Conclusion to uniformly frame the paper as a case study of spectral complexity in bright LAT-detected GRBs, removing "universal" framing language, consistent with the caveat already present in the Introduction ("not a complete or flux-limited sample").

---

## Priority 3 — Editorial / consistency fixes (quick pass, low effort)

- [ ] Introduction lists "three criteria simultaneously" for sample selection but only (i) and (ii) are given — add missing criterion (iii) or fix the count.
- [ ] Introduction (~line 161) states "ΔC-stat criterion of 25 units" — this contradicts the 28.74 value used consistently elsewhere (Abstract, Section 5.2). Fix to one consistent number.
- [ ] Section 4.2 states ROI defined at 12°, but 4.2.1–4.2.4 repeatedly refer to "the 10° ROI." Make consistent.
- [ ] Section 6.4 states a burst-averaged spectrum yields "a systematically weaker" Γ_min bound, but Table 3 shows the T90 (burst-averaged) Γ_min = 507 is actually the *highest* value in that column for GRB080916C. Resolve the contradiction (either fix the claim or explain why this burst is an exception).
- [ ] Section 6.5/6.6 report GRB080916C T90 kT as 44.45 ± 2.05 (Table 5) vs. 44.46⁺²·⁰¹₋₂.₀₂ (Table 6) — same quantity, should match exactly or the pipeline difference should be footnoted.
- [ ] Terminology: "posterior samples," "credible interval," and "MCMC" are used throughout (Section 5.3, Table 2 caption, Section 6.5), but the actual method (per Section 5) is RMFIT with Cash statistics — a frequentist fit — with Monte Carlo error propagation from the fit covariance matrix. This is not Bayesian MCMC sampling. Replace with accurate language (e.g. "Monte Carlo–propagated 1σ intervals") throughout, or clarify if an actual MCMC/Bayesian step was performed somewhere not currently described.
- [ ] Title/Abstract calls the sample "bright" GRBs (selection criterion ii: "High GBM fluence"), but no fluence values are actually reported anywhere. Add a small table of fluences (or peak flux) for the four bursts to substantiate the "bright" framing.
- [ ] SBPL model is cited to Ryde (1999); double-check this is the correct/most standard citation for this functional form, or whether Preece et al. (1994)/Kaneko et al. (2006) should be cited alongside or instead.

---

## Suggested order of operations for the processing agent

1. Fix the LaTeX compile issue (uncomment appendix, verify tables, recompile) — resolves #1 and exposes any further gaps in the appendix content.
2. Run the Priority 3 consistency pass (fast, no new analysis needed).
3. Address #2 (rest-frame f_BB) and #5 (symmetric model grid) — both are re-analysis of existing fits, no new data collection needed.
4. Address #6 (GBM-only refit for GRB131014A) and #7 (threshold robustness test) — both are additional fitting runs on existing data.
5. Address #3 (synchrotron alternative) and #4 (Protassov calibration) — these are the most work-intensive but also the most referee-critical; prioritize accordingly if time-constrained.
6. Address #8, #9, #10 — writing/framing fixes, can be done in parallel with the analysis items above.

---

## review-resolution.md

# Review Resolution Log

Tracks each issue from `grb_paper_weaknesses_and_fixes.md` (§grb_paper_weaknesses_and_fixes.md in this archive) against what was actually decided/done. One entry per issue, updated as each is resolved. Status values: **OPEN** (not yet discussed), **DECIDED** (approach agreed, not yet implemented), **CLOSED** (implemented), **IGNORED** (explicitly out of scope per user instruction).

---

## Priority 1

### #1 — Missing fitted-parameter appendix
**Status:** CLOSED (2026-09-01, superseding the 2026-08-31 IGNORED status)
The appendix block in `GRBResearchPaper/main.tex` was re-activated by the user directly (commit `464a747` `[main-minor-12]`, 2026-09-01, no Claude attribution) and has been live since — `table ??`/`appendix ??` references now resolve. It has since grown: a third "Monte Carlo seeds" section was added on top of it (commit `f80338e`, 2026-09-03). Caught 2026-09-03 while investigating a VERIFY.md discrepancy note — the old IGNORED status had gone stale for 3 days without being re-checked against current `main.tex`.

### #2 — Hidden redshift-dependence in cross-burst f_BB comparison
**Status:** CLOSED (2026-08-31)
**Decision:** recompute $f_\text{BB}$ in a fixed rest-frame band (1 keV–10 MeV) alongside the existing observer-frame value, using $z=4.35$ for GRB080916C and a **swept** $z=0.5$–$5.0$ (fiducial $z=2$ marked) for the other three bursts — user chose the sweep over a single-fiducial value, matching the exact precedent this paper already set for the same problem in `photospheric_radius/pe_er_photosphere.py`. Report both bands' numbers explicitly wherever the "3-4x" headline claim appears, rather than picking one.
**Result:** the cross-burst *ranking* is robust (GRB231129C remains most thermally dominated in both bands), but the *multiplier* is not — observer-frame gives ~3-4x GRB080916C's thermal fraction, rest-frame gives ~2-3x, because GRB080916C's higher redshift shifts its rest-frame band further from the fixed observer-frame band than GRB231129C's fiducial z does. Capture-fraction sanity check passed at every z in the sweep (≥99.9999%), so this isn't a band-truncation artifact.
**Implementation:**
- `GRBResearchWork/codes-for-paper/bb_fraction/bb_flux_fraction.py`: added rest-frame band + per-burst z-sweep (mirroring `pe_er_photosphere.py`'s `REDSHIFTS`/`Z_MIN`/`Z_MAX`/`Z_POINTS`/`Z_FIDUCIAL`), extended `FractionResult`/CSV with `z`, `z_source`, `H0`, `Om0`, `e_min_keV_rest`, `e_max_keV_rest`, `f_bb_rest` (+errors), `bb_fraction_captured_in_band_rest`; added `bb_flux_fraction_rest_vs_z.png/.pdf`. Also fixed a bare RNG call to route through `get_rng` (matching what `bb_fraction.md` already claimed but the code didn't do).
- `bb_flux_fraction.csv` regenerated: 263 rows (was 13) — one row per BB-inclusive interval per redshift.
- `csv_to_latex.py`: added `select_rows()` (fiducial/spectroscopic collapse, mirroring `photospheric_radius/csv_to_latex.py`), new `$z$` and `$f_\text{BB}^\text{rest}$` columns with the dagger-footnote convention; table bumped to `table*`.
- `GRBResearchPaper/tex_files/section-5-data-analysis.tex` (§bb-fraction): added rest-frame definition paragraph, new figure (`fig:bbfraction-rest-vs-z`), and the revised 3-sentence headline paragraph (both bands' numbers + robustness statement + the z-asymmetry mechanism).
- `abstract.tex`, `section-6-discussion.tex`, `section-7-conclusion.tex`: revised to state both bands' multipliers (discussion.tex gets a shorter robustness-confirmation sentence rather than restating numbers, since it didn't have the explicit "3-4x" claim to begin with).
- `bb_fraction.md`: new §2.2b documenting the judgement calls (rest-frame band choice, sweep-vs-fiducial decision, RNG fix), §3.1b validation results, updated §4 results table.
- Backups of every touched file taken before editing, in `_backups/2026-08-31_restframe_fbb/` in each repo.
**Verification:** logged in `VERIFY.md` — bare `pdflatex` compile-error check passed (19 pages, was 18); full PDF content check needs the user's manual pass.
**Committed:** `GRBResearchWork@2fa38b3` (+ fix-up `70b8654`, unrelated stray files accidentally swept into the first commit, see `HANDOFF.md` §12), `GRBResearchPaper@1f1c65f` (bundled with #4 below — same commit).

### #3 — Missing synchrotron / multi-break alternative model
**Status:** CLOSED (2026-08-31) — citation-only
**Scoping finding:** the doc's fuller resolution (fit a 2SBPL/synchrotron-break model to 1-2 strong-BB episodes, compare C-stat/BIC) requires RMFIT, which is not available to Claude in this environment — `GRBResearchWork` has no raw GBM/LAT count-spectrum data (PHA/RSP) and no independent Python fitting pipeline; the per-episode `.fit` files are RMFIT's own *output*, not re-fittable inputs. **User decided:** citation-only, same tradeoff as #4 — no new fits, no numeric comparison.
**Implementation:**
- Added three verified bib entries to `GRBResearchPaper/ref.bib`: `Oganesyan2018` (A&A 616, A138), `Ravasio2018` (A&A 613, A16 — the paper that names the "2SBPL" model the weakness doc itself suggested fitting), `Burgess2020` (Nature Astronomy 4, 174-179). All three verified against arXiv/DOI via web search before adding, not taken from memory.
- `GRBResearchPaper/tex_files/section-1-introduction.tex`: new 5-sentence paragraph after the existing photospheric-evidence literature review, presenting the synchrotron-break alternative as a genuine competing explanation and stating explicitly that it was not tested here.
- `GRBResearchPaper/tex_files/section-6-discussion.tex`: one-sentence caveat at the exact spot (the "shallow shoulder" / "line of death" description of GRB131014A's and GRB231129C's BB signature) where the paper's own language matches what this literature attributes to synchrotron — the most direct engagement point, per the weakness doc's own framing.
**Verification:** bare `pdflatex` compile-error check passed (19 pages, unchanged from the #2/#4 pass). New citations not yet resolved (needs `bibtex`) — logged in `VERIFY.md`.
**Committed:** `GRBResearchPaper@875d1b1`.

### #4 — Statistical validity of the ΔC-stat threshold (Protassov et al. 2002)
**Status:** CLOSED (2026-08-31)
**Decision:** no Monte Carlo null-distribution calibration — cite Protassov et al. (2002) as the caveat and lean on precedent (Guiriec et al. 2011, Fana Dirirsa et al. 2019) having used the same analytic thresholds for the same BB-detection problem, per the doc's own fallback resolution path (step 4).
**Implementation:**
- Added `Protassov2002` bib entry to `GRBResearchPaper/ref.bib` (ApJ 571, 545).
- Added two sentences to `GRBResearchPaper/tex_files/section-4-joint-analysis-results.tex`, end of §4.2 (`subsec:best-models`, after the $\Delta$C-stat threshold paragraph): states the boundary/unconstrained-parameter condition that breaks the asymptotic $\chi^2$ approximation, cites Protassov et al. (2002), and states the threshold is adopted as established practice per Guiriec (2011) / Fana Dirirsa (2019) rather than recalibrated.
- Verified with a bare `pdflatex` pass: compiles clean, no new errors (citations not yet resolved since `bibtex` wasn't re-run — see VERIFY.md).
**Verification:** logged in `VERIFY.md` — full build + citation resolution needs the user's manual check.
**Committed:** `GRBResearchPaper@1f1c65f` (bundled with #2 above — same commit; `GRBResearchWork` has no code changes for this item).

---

## Priority 2

### #5 — Asymmetric model grid (PL extension only tested with BB)
**Status:** CLOSED (2026-08-31, re-confirmed 2026-09-01) — clarification, no new fits, not actually a gap; paper-text wording finalized at user's discretion, see below
**Finding:** checked `results.json` directly across all four bursts (26 episodes) — confirmed zero bare BASE+PL (non-BB) fits exist anywhere; `grb_enums.py` has `BAND_PL`/`CPL_PL`/`SBPL_PL` commented out.
**Full history (project context, not in the paper — user asked to keep the paper text simple):**
1. A fully symmetric grid (BASE, BASE+PL, BASE+BB, BASE+PL+BB) was tried in a much earlier iteration, well before this repo existed, and dropped for complexity.
2. For the actual paper work, only BASE, BASE+BB, and BASE+BB+PL (PL added on top of an already-accepted BASE+BB) were used — bare BASE+PL was never tested as its own comparison point.
3. **The real reason, per the user:** BASE+PL and BASE+BB are not nested/incremental relative to each other — both add 2 dof to BASE, but along different, non-hierarchical branches (one models an extra continuum component, one a thermal one), so a direct likelihood-ratio comparison between them isn't statistically well defined the same way BASE→BASE+BB or BASE+BB→BASE+BB+PL are. Since BASE+BB→BASE+BB+PL only adds 2 dof on top of an accepted BASE+BB and *is* a valid nested comparison, that's the branch that was kept.
4. Separately, the supplementary PL component's physical origin traces to Fermitools' own power-law characterization of the LAT-band data (consistent with `LAT_analysis.md`'s documented `gtlike` output fields: `Index`, `Flux (0.1-100.0 GeV)`) — not an independently RMFIT-tested continuum alternative.
So the reviewer's underlying premise (PL tested asymmetrically as an arbitrary choice) doesn't hold — it's a deliberate, statistically-motivated design (non-nested models aren't directly comparable via LRT), not an oversight.
**Explicitly not claimed anywhere:** whether the very first dropped symmetric-grid exploration (item 1) changed any BEST-model classification — no record of its outcome survives (predates this repo).
**Implementation, revised 2026-09-01:** the two-sentence Fermitools-origin/non-nested-comparison expansion described above as the original implementation was checked on 2026-09-01 and found **not present** in the current file — neither the 2026-08-31 pre-edit backup nor the 2026-09-01 backup taken at the start of that day's session contain it. The current sentence in `section-4-joint-analysis-results.tex` reads: "Additionally, a supplementary power law was used to account for high-energy LAT emission as a further, nested extension of an already-accepted BASE+BB model." — a single sentence, narrower than either the original "only for BASE+BB models" wording or the claimed two-sentence expansion. Cause not determined (possibly overwritten by a manual edit at some point; not investigated further).
**User's decision (2026-09-01), left at user's discretion:** the fuller two-sentence explanation made the paragraph "quite complex" — close #5 as-is with the current single-sentence wording. The user will handle any reviewer pushback on this point directly themselves rather than have the paper text carry the full non-nested-comparison reasoning. The reasoning itself (items 1–4 above) remains valid project context in this file regardless of what the paper text says.
**Verification:** bare `pdflatex` compile-error check passed previously (19 pages, matching the pre-#5 baseline); current wording confirmed present via direct file read, 2026-09-01.
**Committed:** the current single-sentence wording is already committed as part of `GRBResearchPaper@cdaf8bb` ([main-minor-13], 2026-09-01) since that commit covered the whole file. No further paper-text change needed for this item.

### #6 — GRB131014A off-axis angle concern
**Status:** CLOSED (2026-09-01) — written up in the paper, pending the user's PDF read-through (`VERIFY.md`).
**Paper text:** new paragraphs in `GRBResearchPaper/tex_files/section-6-discussion.tex`, inserted after the existing shallow-shoulder/synchrotron-alternative paragraph for GRB131014A+GRB231129C, before the GRB140206B paragraph. States the off-axis coincidence plainly, rules out the direct effective-area mechanism (LAT sensitive only above 100 MeV, cited to the already-used `Atwood2009`, vs. this burst's kT_BB ≈ 27–44 keV — four orders of magnitude apart), then presents the GBM-only refit as the test for the remaining indirect joint-likelihood pathway, with the two primary figures embedded (`fig:gbmonly-kt`, `fig:gbmonly-cstat`) and the per-episode numbers quoted in text. No new citations needed. Two new images copied to `GRBResearchPaper/images/section6/` (new folder, first Discussion-section figures in this paper): `gbm_only_refit_kt_comparison.pdf/png`, `gbm_only_refit_delta_cstat.pdf/png`.
**Verification:** bare `pdflatex` ×2 (per CLAUDE.md's compile-error-catching exception) — clean, 0 errors, 0 undefined references after the second pass, both new images load correctly (main.pdf now 27 pages, up from 26). Content itself (does the new text read correctly, do the figures render legibly at column width) is logged in `VERIFY.md` for the user's manual check, not confirmed here.
**Superseded below:** the "Analysis CLOSED... paper-text framing still to be decided" status and the two entries under it are the pre-write-up state, kept for the record.
**Physics rebuttal established (2026-09-01):** WebSearch-confirmed that Fermi-LAT's effective area is ~0 below 100 MeV regardless of off-axis angle, so the weakness doc's proposed mechanism (off-axis LAT systematics producing spurious low-energy curvature that mimics a BB) cannot directly apply — this burst's BB signature (kT ≈ 27–59 keV) lives entirely in GBM's band, not LAT's. This leaves a softer, indirect concern: LAT data could still bias other continuum parameters through the shared joint likelihood in a way that happens to favour a spurious BB.
**GBM-only RMFIT refit (2026-09-01):** the user ran a GBM-only refit (LAT dropped) for all 5 episodes of GRB131014A and added it to `results.json` as `GRB131014215GBM`. New analysis in `GRBResearchWork/codes-for-paper/gbm_only_refit/` (script `gbm_only_refit.py`, method note `gbm_only_refit.md`, output `gbm_only_refit_comparison.csv` + plot) compares it against the joint fit `GRB131014215`.
**Result — all 5 episodes, clean positive confirmation:** T90, EX0, TR1, EX1 (BAND_BB) and TR2 (SBPL_BB) all pick the *same* BB-augmented model in both fits, kT_BB agrees within 1σ in every case (e.g. TR1: 44.05±1.62 keV joint vs 44.49±1.54 keV GBM-only; T90: 34.92±0.80 keV joint vs 35.36±0.76 keV GBM-only), and the BASE→BASE+BB step clears the paper's own Δ-C-stat ≥ 28.74 threshold by a wide margin in every GBM-only episode (Δ-C-stat = 46–382). This directly closes the #6 concern: the BB detections do not depend on LAT being in the fit, closing the indirect joint-likelihood pathway as well as the effective-area mechanism originally named.
**T90 anomaly caught and fixed in this pass:** the script's independent `any_bb_step_clears_threshold()` check initially flagged T90 alone as inconsistent (`results.json` had no `BAND_BB` fit for GBM-only T90; the recorded winner was a non-BB `BAND` despite `SBPL_BB` clearing the threshold by Δ-C-stat = 491). Traced to a user-side labelling mistake — the T90 GBM-only BB fit had been saved as `BAND_PL` instead of `BAND_BB`. User renamed it and regenerated `results.json`; re-run confirms `BAND_BB` is now present and correctly recorded as BEST, matching the joint fit. All 5 episodes now consistent, no remaining exceptions.
**Not yet decided:** how to frame this in the paper — citation-only physics rebuttal, vs. a fuller paragraph/appendix table citing the GBM-only robustness check (now that it's a clean 5-for-5 result, not a partial one). Deferred to the next pass on this item.

### #7 — No instrumental-artifact / low-energy threshold robustness test
**Status:** CLOSED (2026-09-06) — no paper-text change.
**User decision:** the raised-NaI-threshold (40 keV) GBM-only refit stays a project-internal robustness check, not written into the paper. The core reviewer concern (are BB detections an artefact of the low-energy NaI band?) is already answered — BB is still overwhelmingly required at every episode — and the kT_BB systematic found alongside it doesn't warrant a caveat in the manuscript. No sentence added anywhere; `gbm_only_refit_3way.py`/`gbm_only_refit.md` §6-9 remain the record of the experiment for anyone revisiting this later.

**GBM-only RMFIT refit at a raised NaI threshold (2026-09-06):** the user re-ran GRB131014A's GBM-only refit (already used for weakness #6) a second time, with the NaI lower-energy bound raised to 40 keV (band 40-900 keV, `results.json` key `GRB131014215GBM40keV`), motivated by this burst's extreme off-axis geometry (ROI = 12 deg, zenith = 100 deg per the user). New analysis in `GRBResearchWork/codes-for-paper/gbm_only_refit/gbm_only_refit_3way.py` (method note: `gbm_only_refit.md` Sec 6-9), a three-way comparison against the existing joint (`GRB131014215`) and standard-range GBM-only (`GRB131014215GBM`) fits, kept as a separate script/output set from the existing two-way #6 analysis so as not to change the figures already embedded in `section-6-discussion.tex`.

**Result — mixed, not a clean confirmation like #6:** BB is still overwhelmingly required at every episode at the raised threshold (delta-C-stat = 49-241 against the 28.74 threshold, dropping relative to the standard-range fits as expected from fewer energy channels but nowhere near the boundary), which directly answers the weakness doc's core concern — the BB detections are not an artefact of the low-energy NaI band. But kT_BB shifts measurably higher under the raised threshold in every episode (+2% to +18%), not consistent with zero within 1-sigma in 3 of 5 episodes (T90, TR2, EX1), and EX1's best-fit model itself changes (BAND_BB -> SBPL_BB) — a real, if modest, systematic tied to losing low-energy leverage on the Wien tail, concentrated in the two weakest-BB episodes (TR2, EX1).

**Not yet decided:** how (or whether) to state this in the paper — options include a systematics caveat quoting the size of the shift alongside the existing kT_BB table, a footnote scoped to the two affected episodes, or leaving it out given BB detection itself is unaffected. Deferred to the user, same as #6's post-analysis framing decision.

### #8 — GRB080916C thermal detection vs. prior literature
**Status:** CLOSED (2026-09-01)
**Implementation:** New paragraph in `section-6-discussion.tex`, after the existing Peer2007/photospheric-radius sentence for GRB080916C. Two parts:
1. **Why the BB detection differs:** `Abdo2009FermiObservations080916C` fit only Band functions to five fine time bins, never testing a BB-augmented model, and didn't cover the EX0/EX1 excess windows this paper analyses; `Guiriec2015`'s later 3-component re-fit found only a weak thermal preference in some time-resolved spectra, consistent with (not contradicting) this paper's finding that the BB is confined to the earliest episodes.
2. **The Γ_min tension, stated plainly (user-confirmed framing):** Abdo et al. (2009) reported opacity Γ_min = 887±21 and 608±15 from two fine time bins; this paper's thermal Γ ≈ 750–860 sits below the higher value. Traced the likely mechanism: Abdo's finer bins imply a shorter variability timescale $t_v$ than this paper's conservative episode-duration-based $t_v$, and since $\Gamma_{\min}\propto t_v^{-1/(2\alpha+2)}$ (already stated in §5), a shorter $t_v$ inflates the opacity bound. This paper's own opacity Limit A for GRB080916C tops out at 507 — well below the thermal Γ — so once the different $t_v$ conventions are accounted for there's no indication of an internal contradiction. Stated as the likely explanation, not a closed case.
**Verification:** bibtex + 2×pdflatex clean (26 pages, 0 errors, 0 undefined citations/refs); `pdftotext`-confirmed the new paragraph renders with both citations and the `table 4` cross-reference resolving correctly.

### #9 — Amati-relation interpretation overstated
**Status:** CLOSED (2026-09-01)
**Decision:** soften "confirms...physically self-consistent" to "consistent with...physically self-consistent" (`abstract.tex`, `section-5-data-analysis.tex`), and add literature on time-resolved Amati validity rather than staying silent on it.
**Literature verified via ADS/arXiv abstract fetch (not memory) before citing:** Basak & Rao (2012, ApJ 749, 132) report Pearson $r$ dropping from 0.80 (time-integrated) to 0.37 for fine "intensity-guided" time bins that ignore pulse structure, but recovering to 0.89 for pulse-wise (physically-motivated) binning — directly relevant since this paper's TR/EX episodes are Bayesian-block intervals tied to light-curve structure, not arbitrary fine bins. Ghirlanda et al. (2010, A&A 511, A43) is cited as a complementary point (the Yonetoku $E_p$-$L_\text{iso}$ relation also recovers time-resolved within individual bursts) — confirmed via direct abstract fetch that this paper does *not* address Ep-Eiso/Amati specifically, so it's used only for the Yonetoku point, not misattributed to Amati.
**Implementation:** `section-5-data-analysis.tex` (end of §5.2 Amati subsection): softened sentence + 2 new sentences with the r=0.80/0.37/0.89 statistic and the TR/EX-episodes-are-pulse-scale argument. `abstract.tex`: one-sentence softening only (length constraint). New `Ghirlanda2010`, `BasakRao2012` bib entries added to `ref.bib`.
**Verification:** same bibtex+pdflatex pass as #8, same file — confirmed clean.

### #10 — Overclaiming population-level generality from N=4 sample
**Status:** CLOSED (2026-09-01)
**Implementation:**
- `section-1-introduction.tex`: "we test whether photospheric emission is a universal or selective feature of bright LAT-detected GRBs" → "we examine the incidence and strength of photospheric emission across this sample as a case study of spectral complexity in bright LAT-detected GRBs" — removes "universal" framing, consistent with the existing "not a complete or flux-limited sample" caveat a few lines later.
- `abstract.tex`: "...component of prompt GRB spectra" → "...component of prompt emission in this sample" (unqualified generalization → sample-scoped).
- **Two more instances of the same pattern found while reading, not in the original doc's literal wording but the same overclaiming shape:** `section-6-discussion.tex` and `section-7-conclusion.tex` both had "This confirms that population-level studies... will systematically underestimate..." — a strong "confirms...population-level" claim resting on one burst (GRB080916C). Softened both to "This illustrates how population-level studies... can underestimate...".
**Verification:** same bibtex+pdflatex pass, clean.

### Citation bug found during #8 (not one of the 10 tracked items)
**Status:** CLOSED (2026-09-01) — **user confirmed fix in this pass**
**Problem:** `section-1-introduction.tex`'s GRB120323A sentence cited `\citet{Guiriec2015}`, but `Guiriec2015` in `ref.bib` is the *GRB080916C* paper (ApJ 807, 148) — the one #8 itself needed. The real GRB120323A paper, Guiriec et al. 2013 (ApJ 770, 32), had no bib entry at all. Verified via ADS (`2013ApJ...770...32G`) before fixing.
**Implementation:** added `Guiriec2013` to `ref.bib`; changed the intro's citation to `\citet{Guiriec2013}`. Both PDFs (`Guiriec2013 - ...pdf`, `Guiriec2015 - ...pdf`) were already present in `Literature Review/`, just the bib/citation link was wrong.
**Also fetched and saved to `Literature Review/` (per user request):** `Basak2012 - ...pdf` (arXiv:1202.3089) and `Ghirlanda2010 - ...pdf` (arXiv:0908.2807) — neither existed there before; confirmed valid (correct title/author on page 1) after download. `Abdo, A. A. et al. (2009).pdf` (already present, generically named) was confirmed via `pdftotext` to be the GRB080916C paper, not a duplicate of the folder's other, differently-titled GRB090902B Abdo 2009 paper — no new fetch needed for it.
**Backups:** all touched files backed up first, in `_backups/2026-09-01_priority2_items_8_9_10/`.
**Committed:** not yet — pending, per the no-auto-commit rule.

---

## Priority 3 — editorial/consistency pass

**Status: CLOSED (2026-09-01)**, all 8 items. Backups of every touched file taken first, in `_backups/2026-09-01_priority3_editorial/` in each repo.

1. **"Three criteria" count mismatch — CLOSED.** Checked git history: the enumerate block (`section-1-introduction.tex:47`) has only ever had 2 items; no third criterion was ever dropped. Fixed "three criteria" → "two criteria". No third criterion exists anywhere else in the paper or `PLAN.md` either, so nothing was added.

2. **"25 units" vs. 28.74 ΔC-stat — CLOSED.** `section-1-introduction.tex:54` said "25 units"; `section-4-joint-analysis-results.tex:26-28` establishes 28.74 (Δk=2)/36.86 (Δk=4) as the actual thresholds. Changed "25" → "28.74" (matches the Δk=2 case the intro sentence is describing generically).

3. **ROI 12° vs. 10° — CLOSED.** `section-3-data-preparation-and-analysis.tex:72` defines the ROI at 12°; the four per-burst subsubsections (lines 88, 96, 106, 115) each said "10° ROI". Verified against ground truth — every `LAT_analysis/*/fit_results_*.txt` in `GRBResearchWork` reports `ROI, 12` — and changed all four "10" → "12".

4. **§6.4 Γ_min "systematically weaker" claim — CLOSED.** Verified against `lorentz_results.csv`: GRB080916C's T90 Γ_min (507.48) is actually the *highest* value in the column, not weaker than any time-resolved episode (TR3: 503.996, the next closest). Traced the mechanism precisely rather than just softening blind: T90's `tau_hat` (7.10e9) is lower than TR3's (7.24e9) as expected from its longer `t_v` (62.976 vs. 40.256 s), but Γ_min also depends on the fitted spectral index (`alpha_LS` = 2.2507 for T90 vs. 2.2156 for TR3) through the `e1,e2,e3` exponents in `compute_gamma_min`, and this difference nearly cancels the `t_v` effect. `section-5-data-analysis.tex`: removed "systematically", added a clause on the spectral-index dependence, citing GRB080916C's T90-vs-TR3 near-tie as the concrete example already in `tab:lorentz`. **User confirmed this approach (soften + brief mechanism note) over a minimal deletion or leaving it open**, via `AskUserQuestion`.

5. **kT mismatch (`tab:bbfraction` vs. `tab:photospheric`) — CLOSED.** Root cause traced precisely: `bb_flux_fraction.py` takes `kt_bb = model.get_parameter_value("kt_bb")` — the literal RMFIT fit value and its native symmetric error, **not** an MC draw — while `pe_er_photosphere.py`'s `interval_samples()` draws kT via MC and reports the median/16th/84th percentile. This was already documented precisely in both folders' method notes (`bb_fraction.md` §4.5, `photospheric_radius.md` §4.5, `BUGS.md` OBS-09) — the paper-facing table captions just didn't say so. Fixed both `csv_to_latex.py` scripts' captions to state this explicitly and cross-reference each other; added a pointer in both `.md` files back to this resolution. Regenerated both tables — only caption text changed, all data values identical to before.

6. **MCMC/posterior/credible-interval language — CLOSED.** Actual method is RMFIT + Cash-stat + Monte Carlo error propagation, not Bayesian MCMC. Fixed in `section-4-joint-analysis-results.tex` (lines 56, 65, 66 — "credible interval"/"posterior samples" → "Monte Carlo-propagated interval"/"Monte Carlo parameter samples") and in `codes-for-paper/amati_relationship/csv_to_latex.py`'s caption string ("MCMC samples" → "Monte Carlo samples" — confirmed via `amati_helpers.py` that it's plain `np.random.default_rng`, not MCMC). Regenerated `amati_relationship_table.tex`; only the caption changed.

7. **SBPL citation — CLOSED.** Verified via web search: Ryde (1999) is the correct originating citation for the SBPL functional form; Kaneko et al. (2006) is the standard re-parameterization cited alongside it in the literature (Preece 1994 is not the relevant citation for this functional form). `section-4-joint-analysis-results.tex:8`: `\citep[\sbpl]{Ryde1999}` → `\citep[\sbpl]{Ryde1999, Kaneko2006}` (`Kaneko2006` already in `ref.bib`).

8. **"Bright" framing with no fluence table — CLOSED.** Found the existing `codes-for-paper/fluence/grb_fluence.py` computed fluence over `FluxFluenceCalculator`'s default 10 keV–1 MeV band, not the paper's stated 8 keV–40 MeV GBM band (`section-1-introduction.tex:50`) — so it couldn't be cited as-is. **User confirmed the full-scope fix** (recompute + compliant CSV + new generator + method note) over a lighter caveat-only fix, via `AskUserQuestion`. Fixed `grb_fluence.py` to pass the correct `log_energy_range`, regenerated `flux_fluence.csv`/`flux_energy_flux.csv` with full CLAUDE.md-convention columns (paper-name `grb_name`, `e_min_keV`/`e_max_keV`, `n_samples`, `seed`, units baked into column names); wrote `codes-for-paper/fluence/csv_to_latex.py` (new, T90-only rows) and `fluence.md` (new per-folder method note — this folder previously had none). New `tex_files/generated/fluence_table.tex` (`tab:fluence`) referenced from `section-1-introduction.tex` right after the two selection criteria.

**Verification:** logged in `VERIFY.md` — bare `pdflatex` compile-error check pending; full PDF content check (new fluence table, revised wording throughout, both kT-caption footnotes, Γ_min caveat) needs the user's manual pass.

---

## GRBResearch_Issues_List.md

# Issues List for Revision

*"Thermal and Non-Thermal Emission in Bright Fermi Gamma-ray Bursts"* — I. Siddique, S. A. M. Bukhari
Reviewed draft: `GRBResearch.pdf` (Sept 4, 2026)

The four items below are outstanding from the latest reviewed draft. Three were flagged in an earlier round and remain unresolved; the fourth (Issue 4) is newly found in this pass. All other previously flagged Priority-2/3 items (abstract attribution, the σ_T gloss, the Γ_min reconciliation paragraph, the TS<25 footnote, the Protassov caveat, etc.) are confirmed present and correct in this draft and are not repeated here.

---

## Issue 1 — Table 2 episode boundaries for GRB231129C contradict Tables 8 and 12

**Status:** RESOLVED (confirmed by user, 2026-09-06)

The root cause was a transcription error in the (then hand-typed) Table 2 source: GRB231129C's `EX0`/`TR1` cells wrongly ended at `1.792` instead of `3.136`, which fabricated the extra "bin 2" and shifted the real `TR2` into a mislabeled "bin 3". Table 2 is now generated directly from `results.json` (`GRBResearchWork/codes-for-paper/time_integrated_table/time_integrated_table.py`), which fixed the value and collapsed the row back to the correct 4-bin scheme (`EX0`, `1`≡TR1, `2`≡TR2, `EX1`) matching Tables 8/12. Confirmed against the rendered PDF by the user.

**Where:** Table 2 (main text, Sec. 4.1) vs. Table 8 (Appendix A) and Table 12 (Appendix B.4).

**Problem:**
Table 2 lists five time-resolved bins for GRB231129C:

```
EX0 (−0.192–1.792 s), 1 (0.384–1.792 s), 2 (1.792–3.136 s), 3 (3.136–7.296 s), EX1 (3.136–10.048 s)
```

but Tables 8 and 12 (which carry the actual LAT and joint-fit results) both use a different, four-bin scheme for the same burst:

```
EX0 (−0.192–3.136 s), TR1 (0.384–3.136 s), TR2 (3.136–7.296 s), EX1 (3.136–10.048 s)
```

Table 2's bin "2" (1.792–3.136 s) does not appear anywhere in the analysis tables at all. Table 2's bin "3" (3.136–7.296 s) is the same time range as Table 8/12's "TR2," but carries a different label. Every downstream reference to this burst's episode numbering (e.g. Fig. 2(b), Sec. 5.3.4, Sec. 6.5) follows the Table 8/12 scheme, not Table 2. The other three bursts (080916C, 131014A, 140206B) were checked and are internally consistent — only the GRB231129C row of Table 2 is affected.

**Suggested fix:**
Replace the GRB231129C column of Table 2 with the four-bin scheme already used in Tables 8 and 12 (EX0, TR1, TR2, EX1), dropping the extra "1"/"2" split.

---

## Issue 2 — Equation (4): Amati relation exponent is inverted

**Status:** Unresolved (carried over from prior review round)

**Where:** Section 6.3, Equation (4) and surrounding text.

**Problem:**
The text defines y = log₁₀(E_iso/E_0,iso), x = log₁₀(E_i,p/E_0,p), and states the slope m = 0.5. As written, this gives E_iso ∝ E_p^0.5. The Amati relation is E_iso ∝ E_p² (m ≈ 2, not 0.5) — this is also the standard form cited from Amati (2006) and Nava et al. (2012) in the same sentence.

Table 3's own numbers confirm the data follow the steep (m≈2) relation, not the shallow one written in the text: going from E_i,p ≈ 1.3–1.4 MeV (TR4/TR5/EX1) up to ≈7 MeV (EX0/TR1), E_iso rises roughly seven-fold in log(E_p) but far more than the square root of that in log(E_iso) — consistent with slope ≈1–2, not 0.5. Figure 8 itself appears correct, since its best-fit line is taken directly from Fana Dirirsa et al. (2019) rather than recomputed from Equation (4); only the text/equation is wrong.

**Suggested fix:**
Correct m to its intended value (≈2, or whatever value is actually cited/used for the plotted best-fit line) so the equation matches both the cited literature form and Figure 8/Table 3.

---

## Issue 3 — Fitted β is inconsistent with the LAT photon index (GRB131014A, GRB231129C)

**Status:** CLOSED — discussed at length, no text change made (2026-09-06)

Investigated in depth: the numbers check out exactly as the reviewer reported, and a diagnostic
comparison (`GRBResearchWork/codes-for-paper/beta_lat_consistency/beta_lat_consistency.md`) confirms
the underlying cause is as suspected — β describes the GBM-dominated continuum curvature, not a
genuinely LAT-constrained slope, for these low-LAT-photon-count episodes. A short in-text caveat was
drafted for the two affected subsubsections in `section-4-joint-analysis-results.tex`, but on review
the author judged it unnecessary: this is expected behavior to anyone in the field (Band/SBPL's β
has no physical mandate to describe a distinct high-energy component two-plus decades above where
the fit is actually constrained) and doesn't need spelling out in the paper itself. No text change
was kept — closing per author's discretion, not because the reviewer's diagnosis was wrong.

**Where:** Table 9/10/12 (fitted β) vs. Table 8 (LAT photon index) — specifically GRB131014A T90 and GRB231129C T90/EX0/TR1.

**Problem:**

| Burst / episode | Fitted β (±) | LAT photon index (±) | Discrepancy |
|---|---|---|---|
| GRB131014A, T90 | −2.904 ± 0.035 | −1.99 ± 0.31 | ≈2.9σ |
| GRB231129C, T90 | −4.438 ± 0.275 | −2.46 ± 0.41 | ≈4.0σ |
| GRB231129C, EX0 / TR1 | steep (SBPL+BB fit) | −6.00 (fixed, TS<25) | large, low-significance |

The high-energy spectral index β from the joint fit should, in principle, be constrained by the LAT photons at those episodes, but it disagrees with the LAT-band index at ~3–4σ for these two bursts (and by even more for GRB231129C's EX0/TR1, where the LAT index is pinned at −6.00 with TS<25 — a low-significance association). A fitted β ≈ −4.4 (GRB231129C) is also not a value seen elsewhere in the GRB spectral-index literature. The likely cause is that these episodes have very few LAT photons (~10 or fewer), so β is effectively set by the GBM extrapolation rather than by real high-energy data — but this is not stated anywhere in the draft.

**Suggested fix:**
Add a short caveat (e.g. in Sec. 5.3.2/5.3.4 or the Discussion) noting that for these low-LAT-count episodes, β is effectively a GBM-only extrapolation and should not be read as an independently LAT-constrained index — similar in spirit to the off-axis-systematics paragraph already written for GRB131014A in Sec. 7.1.

---

## Issue 4 — "...C-Stat units" value is wrong in Section 5.3.2 (GRB131014A, TR2) — NEW

**Status:** RESOLVED — no source fix needed (confirmed 2026-09-06)

The current source (`GRBResearchWork/GRBResearchPaper/tex_files/section-4-joint-analysis-results.tex:95`) already reads "74 C-Stat units", not "82". Checked against Table 10's actual C-stat values (SBPL `778.9741` − SBPL+BB `704.6270` = `74.35` → 74; Band `786.98` − SBPL `778.97` = `8.01` → matches "eight" too). Git history shows exactly one commit changed "a further 82 C-Stat units" → "74 C-Stat units", with no later revert — confirming this note's own suspicion: the Sept 4 reviewed PDF was built from a stale/older draft, and the fix was never actually lost.

**Recommendation:** rebuild `main.pdf` fresh from the current source (`pdflatex` → `bibtex` → `pdflatex` ×2) before the next review round, to avoid re-flagging fixes that are already in place.

**Where:** Section 5.3.2, second paragraph, vs. Table 10 (Episode II: 2.432–4.160 s).

**Problem:**
The text reads: "For TR2, the SBPL base model was favored over Band by eight C-Stat units. Adding a blackbody component to SBPL improved the fit by a further **82** C-Stat units..." The "eight" figure checks out (Band 786.98 vs. SBPL 778.97 → difference = 8.01). The "82" figure does not: Table 10 gives SBPL = 778.9741 and SBPL+BB = 704.6270, so the C-stat improvement is 74.35, not 82.

Note: this exact wording — "82" → "74" — was one of the items already confirmed fixed in the previously-reviewed `main.pdf`, so its reappearance here suggests either an older draft was uploaded by mistake, or the fix was inadvertently reverted. Worth double-checking which source file is authoritative before re-applying the fix.

**Suggested fix:**
Change "82 C-Stat units" to "74 C-Stat units" (or the precisely rounded value of 74.35) in Section 5.3.2.

---

*Prepared from a review of GRBResearch.pdf (draft dated September 4, 2026). All page/section/table references are to that draft.*

---

## quick-fixes-mid-priority.md

# Quick fixes — mid priority (content/argument strength)

## Done

- [x] **Protassov marginal case** — your ΔC-stat values are mostly 46–382, comfortably clear of threshold; EX1 (GRB131014A, joint fit) at 30.4 is the one genuinely marginal case, confirmed against the exact same number already independently stated in `section-6-discussion.tex`'s GBM-only-refit robustness check. Added to the Protassov caveat paragraph in `section-4-joint-analysis-results.tex`: "The concern is material only for the single episode near threshold (EX1, $\Delta C\text{-stat} = 30.4$); every other detection exceeds it by a wide margin." *(2026-09-04)*
- [x] **Asymmetric model grid** — same question as `quick-fixes.md`'s (§quick-fixes.md in this archive) supplementary-PL item; already settled in `review-resolution.md` (§review-resolution.md in this archive) #5 and reconfirmed by the user directly in this same session: BASE+PL and BASE+BB aren't nested (both add 2 dof along different, non-hierarchical branches), so a direct likelihood-ratio comparison between them isn't statistically well-defined the way BASE→BASE+BB or BASE+BB→BASE+BB+PL are — that's the branch actually tested. No re-litigation, no new paper text. *(2026-09-04)*
- [x] **TS<25 episodes** — verified: "Table 4" = `tab:lorentz` (covers GRB080916C only, hence never flags other bursts); "table 8" = `lat_info_table.tex` (`tab:burst_table`), already flags GRB231129C's EX0 (TS=8.640) and TR1 (TS=9.522) via footnotes — matches the user's numbers exactly. Added to `section-4-joint-analysis-results.tex`'s GRB231129C paragraph, right after the existing BB-significance sentence: "For TR1 and EX0, the \ac{LAT} association is not statistically secure ($\mathrm{TS} < 25$; \cref{tab:burst_table}), so the joint fit's BB detection in these two episodes is effectively \ac{GBM}-dominated rather than resting on the LAT data." *(2026-09-04)*
- [x] **Γ_min vs thermal Γ reconciliation** — verified the math (t_v ≈ 1.669 s to reconcile our Γ_min=507 with Abdo's 887, using this episode's actual α=2.2507). Added to `section-6-discussion.tex` right after the existing Γ_min=507 sentence, reconciliation first then a reworded (not removed) caveat about $f_1$ not being reconciled. Also added a small reusable script, `codes-for-paper/gamma_min_reconciliation/reconcile_tv.py` (+ its required method note), since this kind of literature-Γ_min check could recur — it reads $(\Gamma_{\min}, t_\text{v}, \alpha)$ straight from `lorentz_factor/lorentz_results.csv` rather than hardcoding them, and is now registered in `runner_registry.yaml` so it reruns automatically alongside the rest of the pipeline. *(2026-09-04)*
- [x] **Table 8 pinned indices** — confirmed real: GRB140206B TR2 and GRB231129C EX0/TR1's raw fitted photon indices are -5.999986/-5.999931/-5.999390 — all within 0.0007 of -6.00, a fit-boundary pin rather than three independent measurements. Fixed at the source (`LAT_analysis/csv_to_latex.py`, not hand-edited): added a tolerance-based pin detector (`PINNED_INDEX_VALUE`/`PINNED_INDEX_TOL`) and a shared footnote letter applied to all three rows' index cells (GRB231129C's two rows correctly carry *two* footnotes now — their existing TS<25 flux-limit note plus the new pin note). Regenerated and copied into the paper repo. *(2026-09-04)*
- [x] **Table 8 caption: explain the TS<25 upper-limit footnotes** — added a fourth caption sentence (sourced off `TS_SECURE_DETECTION` so the "25" can't drift out of sync with the actual threshold): "For episodes with $\mathrm{TS} < 25$, the \ac{LAT} association is not statistically secure; footnotes give the corresponding 95\% photon flux upper limit from Fermitools' `gtlike` `UpperLimits` tool in place of a point measurement." Regenerated and copied into the paper repo. *(2026-09-04)*

## Deferred

- **GRB231129C zenith cut violation** — §4.2.4 kehta hai 12° ROI ka outer boundary ~100.2–100.3° tak jata hai, jabke zenith cut 100° hai. Ye maana to gaya hai lekin bataya nahi ke kya kiya gaya — exposure loss correct hua ya nahi? *(2026-09-04: checked the repo — the actual gtmktime/gtselect commands were run outside it, not checked in, so I can't verify from source what was actually done. Confirmed with the user this was not actually corrected. Left as-is for now — properly addressing it would mean going back to LAT data regeneration, a much bigger task than a text fix. Revisit when that's in scope.)*

---

## quick-fixes.md

# Quick fixes — Priority 3 (editorial, sab quick)

## Done

- [x] **EX–A/EX–B vs EX0/EX1** — unified to EX0/EX1 everywhere (per-burst appendix tables, LAT-info appendix table, `time-integrated-table.tex`). Root cause: `log_to_latex_parser.py`/`LAT_analysis/csv_to_latex.py` derived the label from a directory letter suffix that was always "A" for *both* excess episodes; now derived from actual interval bounds against T90. Along the way: fixed a `...GBM` refit directory silently merging into its base burst's episode list, and a regex that dropped whichever model printed last in a fit run (recovered GRB231129C's missing EX0 `\band` row). Per-burst tables and the LAT-info table are now regenerated end-to-end (`generate_table_from_log.py`) instead of hand-maintained; BEST-row highlighting and citation macros are auto-applied. *(2026-09-04, commits `main-minor-80/81` in GRBResearchWork, `main-minor-18/19` in GRBResearchPaper)*
- [x] **§B.3 heading "GRB140602B" → GRB140206B** — auto-fixed by the regeneration above; the typo only existed in the hand-maintained copy. *(2026-09-04)*
- [x] **§4.2.2 unescaped "<"** — already fixed by the user directly ("reached < 90°" → "was less than 90°"); confirmed no other text-mode `<` remains anywhere (remaining `<` occurrences are all safely inside `$...$` math mode). *(2026-09-04)*
- [x] **Table 13 caption "(not itself MC! (MC!)-derived)"** — root cause was `\ac{MC}` used with no acronym defined for "MC"; added `\acro{MC}{Monte Carlo}` to `glossary.tex`, and the generator's cref text was separately changed to `\acl{MC}` (long form). Both together render correctly now. *(2026-09-04)*
- [x] **Line 433 "given in appendix appendix B"** — `\cref` already auto-prepends "Appendix"; removed the redundant manual "appendix" from the sentence in `section-4-joint-analysis-results.tex`. *(2026-09-04)*
- [x] **§3.1 "Detected by both the GBM) and the LAT"** — stray `)` removed after `\ac{GBM}` in `section-2-grb-observations.tex`. *(2026-09-04)*
- [x] **§4.2 "an additional are of 25°"** — fixed to "area", and reworded the surrounding double-verb clause ("...was checked" → "...checked") in `section-3-data-preparation-and-analysis.tex`. *(2026-09-04)*
- [x] **§5.3.3 "models is" + undefined peak energy ratio** — "model is" (singular); defined the ratio as $E_\text{p}/E_\text{peak,BB}$, the continuum-to-thermal peak energy ratio using the $E_\text{peak,BB}\approx3.92\,kT_\text{BB}$ relation already established earlier in the same section — verified against GRB140206B's actual TR1/EX0 fit values (ratios 9.9 and 10.4, consistent with "$\geq 9$"). *(2026-09-04)*
- [x] **Abstract SAFE/BEST attribution** — split into two clauses: SAFE/UNSAFE error thresholds cite `\citet{Kaneko2006}` as before; MARGINAL is flagged as "extended here" (confirmed by the user: MARGINAL is not standard literature, it's this paper's own one-parameter-40–50%-error relaxation); BEST is now attributed to its actual stepwise ΔC-stat procedure (`\citep{Nava2012, Gruber2014}`), not Kaneko2006. *(2026-09-04)*
- [x] **§5.3.2 "further 82 C-Stat units"** — verified against `table_131014A_large.tex`'s actual TR2 cstat values: SBPL→SBPL_BB is 74 units (778.9741 → 704.6270); 82.35 is BAND→SBPL_BB, the cumulative figure, not the direct step. Changed "a further 82" → "74". *(2026-09-04)*
- [x] **§5.3 "energy flux in units of erg cm⁻²" missing s⁻¹** — root cause: `\eFlnc` (correctly a *fluence* macro, erg/cm², matching `\pFlnc`) was misused to describe *flux* (a rate). Added a proper `\eFlux` macro (erg/cm²/s, matching the existing `\pFlux` pattern) and switched the one usage site to it. *(2026-09-04)*
- [x] **Eq. (10)/(11) σ vs σ_T collision** — verified both equations against Pe'er et al. 2007 (arXiv:astro-ph/0703734) directly: `eq:script_r`, `eq:peer_gamma`, `eq:peer_r0` match the paper's equations 1, 4, 5 exactly, every constant included (1.06, 1.48, 4^{3/2}, all exponents) — no transcription error. The source paper itself uses σ (Stefan-Boltzmann) and σ_T (Thomson) side by side too, but glosses σ_T in text; ours didn't. Added the missing gloss and switched `\sigma_T` → `\sigma_\mathrm{T}` (upright subscript, so "T" reads as a label rather than a variable). *(2026-09-04)*
- [x] **§5.3.4 "could not be constrained in the second episode"** — root cause: the sentence was ambiguous, not simply wrong. `SBPL_BB` (the extension of TR2's *winning* base model) genuinely never converged (confirmed against both the old and newly-regenerated `table_231129C_large.tex` — identical row set, not something my earlier fixes touched), but three non-winning BB variants (BAND_BB, CPL_BB, PL_BB) *did* converge with real κT values, which the original wording could easily be misread as denying. Reworded to name `\sbplbb` and `TR2` explicitly, and named EX1 explicitly in the following sentence for parallelism. *(2026-09-04)*
- [x] **Table 1 burst order (080916C, 140206B, 131014A, 231129C) differs from the rest of the paper** — root cause: `codes-for-paper/fluence/grb_fluence.py`'s `GRB_LIST`, fixed to `["080916C", "131014A", "140206B", "231129C"]` matching every other table's convention. Regenerated `flux_fluence.csv`/`flux_energy_flux.csv` (real N_SAMPLES=10000 run, user's own machine) and `fluence_table.tex`, synced into the paper repo via `sync_paper_assets.py` — order now correct everywhere. *(2026-09-04)*
- [x] **§5.3.3 "that was had"** — fixed to "that had" in `section-4-joint-analysis-results.tex`. *(2026-09-04)*
- [x] **§6.1.1 "less then"** — fixed to "less than" in `section-5-data-analysis.tex`. *(2026-09-04)*
- [x] **§4.2.2 "was 70°, exceeded"** — fixed to "exceeding" (participial clause) in `section-3-data-preparation-and-analysis.tex`. *(2026-09-04)*
- [x] **§6.5 "roughly 70% above the ≈0.072"** — verified against `bb_fraction_table.tex`'s actual GRB131014A values: the "other four episodes" span $f_\mathrm{BB}=0.0715$–$0.0784$ (T90 alone is 0.0784, notably higher than EX0/TR1/EX1's tight 0.0715–0.0720 cluster), so a single "≈0.072" understates the range. Changed to "≈0.072–0.078". *(2026-09-04)*
- [x] **Seed citations gap: tables missing seed numbers, Table 13 mapping only figures** — confirmed real: `amati_relationship_table.tex`, `lorentz_table.tex`, `lorentz_table_limit_b.tex` had no seed in caption (unlike `fluence_table`/`bb_fraction_table`/`photospheric_table`, which already did). Added `(seed $N$)` to all three captions, reading the seed live from the CSV/module constant rather than hardcoding (matching the existing `bb_fraction`/`photospheric` pattern). Also extended `seed_registry.yaml`'s `amati_relation`/`lorentz_limit_a`/`lorentz_limit_b` entries so the appendix seed table's Fig./Table column now crefs `tab:eiso`/`tab:lorentz`/`tab:lorentz_limit_b` alongside their figures, not just the figures. Regenerated all four affected tables + the seed table, synced into the paper repo. *(2026-09-04)*
- [x] **Eq. (12) orphaned + F^ob undefined** — root cause: my earlier σ_T gloss (added while fixing the σ/σ_T collision item above) landed *between* `eq:peer_gamma` and `eq:peer_r0`, breaking the "Combining R... yields... [eq][eq]" lead-in that originally spanned both equations as one unit. Moved the gloss to a single trailing "where" clause after both equations (restoring the shared lead-in), and used the same clause to define $F^\mathrm{ob}$ (never previously defined despite being used in `eq:peer_gamma`) — verified its meaning directly against Pe'er et al. (2007)'s own text: "the total (thermal + non-thermal) observed flux," matching our own lead-in's "the total flux" that was never tied to the symbol. *(2026-09-04)*
