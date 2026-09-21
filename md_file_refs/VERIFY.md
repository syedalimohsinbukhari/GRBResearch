# PDF build verification queue

Entries Claude adds when a change needs confirming against the rendered paper PDF, per `CLAUDE.md`'s "PDF build
verification" convention. Claude does not run `pdflatex`/`pdftotext` to check these — build (`pdflatex`→`bibtex`→
`pdflatex`×2 from `GRBResearchPaper/`, per `HANDOFF.md` §1) and check manually, then mark resolved (or report back what
you found).

Status key: **PENDING** (needs a manual check), **CONFIRMED** (checked, matches), **FAILED** (checked, doesn't match —
needs a fix).

---

### 2026-09-06 — Time-integrated/time-resolved duration table (`table:time-integrated-duration`) automated; one numeric error fixed

**What changed:**

- New generator `GRBResearchWork/codes-for-paper/time_integrated_table/time_integrated_table.py`: reads `results.json`
  interval strings directly (no CSV needed — every cell is a raw interval boundary, nothing fit or MC'd) and emits the
  table, replacing a table that had been hand-typed with no generator anywhere in the repo.
- `GRBResearchPaper/appendices/time-integrated-table.tex` is now a two-line
  `\input{tex_files/generated/time_integrated_table}` wrapper, same pattern as `appendix_LAT_info.tex`/
  `appendix_seeding_table.tex`. `codes-for-paper/table_registry.yaml` moved this entry from `no_known_source` to
  `active`; `sync_paper_assets.py` now copies it automatically.
- Diffing the generator's output against the previously-published table caught a genuine transcription error:
  **GRB231129C's `EX0` and `TR1` rows both showed end time `1.792`; correct value is `3.136`** (per
  `results.json["GRB231129779"]`: `"EX0 -0.192_3.136"`, `"TR1 0.384_3.136"` — and by CLAUDE.md's own definition, `EX0`
  shares `TR1`'s end time, so the two must agree). Fixed automatically by the regenerated table; no other cell changed.
- Header row now uses the paper's `\grbxxxxxxx` macros (`tex_files/preamble.tex:128-131`) instead of literal
  `GRB080916C` text — those macros expand to exactly that text, so this is a zero-visual-diff consistency change, not a
  content change.
- Dropped the dead, already-commented-out "TR breakdown episodes" (`BR--A/B/C`) rows for GRB131014A — their values don't
  appear in `results.json` under either `GRB131014215` or `GRB131014215GBM`, and were never rendered (commented out) in
  the old table either.
- Full method notes in `GRBResearchWork/codes-for-paper/time_integrated_table/time_integrated_table.md`.
- Bare `pdflatex` compile check: 0 errors, 29 pages, unchanged page count from before this table's change — confirms
  LaTeX validity only, not the content below.

**What to check:**

1. Table `table:time-integrated-duration` (§3, "Time-integrated and time-resolved duration for the selected GRBs"):
   GRB231129C's `EX0` row and its `1` (TR1) row should now both show `-0.192` to `3.136` and `0.384` to `3.136`
   respectively (previously both ended at `1.792`).
2. No other cell in this table changed — every other GRB/episode combination should render identically to before.
3. GRB name header row renders as `GRB080916C`, `GRB131014A`, `GRB140206B`, `GRB231129C` (unchanged in appearance, now
   macro-driven).

**Why:** The old table's `1.792` value for GRB231129C was self-inconsistent with the paper's own stated `EX0`/`TR1`
relationship and with `results.json`; if the fix didn't actually render, that specific duration would still be wrong in
a published table.

**CONFIRMED by user (2026-09-06)** — checked against the rendered PDF. This also closes Issue 1 in
`GRBResearch_Issues_List.md` (the reviewer-flagged phantom "bin 2" for GRB231129C was exactly this transcription error).

Status: **CONFIRMED**

---

### 2026-09-06 — Amati relation constants (Eq. 4/§5.2) corrected to match Fig. 8 and cited source

**What changed:** `GRBResearchPaper/tex_files/section-5-data-analysis.tex`, the paragraph right after Equation~\ref{eq:
amati_linear} (§"Amati relationship").

- Old text: $E_\text{0,p} = 300$ keV, $m = 0.5$, $k = 0$, cited to `\citep{Amati2006, Nava2012}` — self-inconsistent
  with the immediately preceding sentence ("Following Fana Dirirsa (2019)...") and matching no real source. This was
  Issue 2 in `GRBResearch_Issues_List.md`.
- New text: $E_\text{0,p} = 950$ keV, $m = 1.16 \pm 0.37$, $k = 1.67 \pm 0.16$, cited to `\citep{FanaDirirsa2019}` —
  verified directly against Fana Dirirsa et al. (2019), Table 3 "F10" row (their own 25-GRB *Fermi*-LAT
  sample, $E_\text{iso}$ computed over 1 keV–$10^4$ keV, matching this paper's own $S_\text{bol}$ integration band), and
  against `codes-for-paper/amati_relationship/amati_helpers.py:167-172`'s `amati_relationship_dirirsa2019()`, which
  hardcodes these exact values (`e_i_peak_norm=950.0, k=1.67, m=1.16`) to draw Figure 8's best-fit line. So the fix
  makes the text match both the literature source it already claims to follow and the actual plotted figure — no change
  to Figure 8, Table 3 ($E_\text{iso}$/$E_\text{i,p}$ values, computed independently via direct flux integration, never
  touch these constants), or any other section (no other passage restates the old wrong numbers).
- Added one more sentence stating $\sigma_\text{ext} = 0.47 \pm 0.12$ (also from Fana Dirirsa 2019, Table 3 F10 row) —
  this is the third parameter (alongside the $m$/$k$ errors) that the same code function uses to set the width of the
  1σ/2σ/3σ confidence bands in Figure 8's caption, which the text previously never explained.

**What to check:** §5.2 "Amati relationship", the sentence (s) right after Eq. 4: should now read $E_\text{0,p}=950$
keV, $m=1.16\pm0.37$, $k=1.67\pm0.16$ (citing Fana Dirirsa 2019), plus the new $\sigma_\text{ext}=0.47\pm0.12$
sentence — not the old $300$ keV/$0.5$/$0$ values.

**Why:** The old numbers didn't match Figure 8's actual best-fit line, didn't match the source the paragraph claims to
follow, and matched no real Amati relation fit in the literature — a reader checking the equation against the figure
would find them inconsistent.

Status: **CONFIRMED** (2026-09-06, user checked rendered PDF)

---

### 2026-09-06 — Removed per-table wrapper files; main.tex/section-3 now `\input` generated tables directly

**What changed:** Pure plumbing, no table content touched.

- `appendices/appendix_LAT_info.tex`, `appendices/appendix_seeding_table.tex`, `appendices/time-integrated-table.tex`
  deleted — each was only a comment plus one `\input{tex_files/generated/...}` line.
- `main.tex`: `\input{appendices/appendix_LAT_info}` → `\input{tex_files/generated/lat_info_table}`;
  `\input{appendices/appendix_seeding_table}` → `\input{tex_files/generated/seed_table}`.
  `section-3-data-preparation-and-analysis.tex` already pointed at `tex_files/generated/time_integrated_table` directly
  (user's own edit, prompting this cleanup).
- `codes-for-paper/table_registry.yaml` and `CLAUDE.md`'s "Generated LaTeX tables" convention updated to describe the
  direct-`\input` pattern instead of citing the now-deleted wrapper files.

**What to check:** Nothing content-wise — this only changes which file does the `\input`, not what renders. A bare
`pdflatex` (no `bibtex`, so citation/reference warnings are expected and not meaningful here) compiled with 0 errors and
all three `\input` targets resolved. Worth a glance during your next normal full build that `tab:burst_table` (LAT
appendix), `tab:seed_table`, and `table:time-integrated-duration` all still render exactly as before — they should,
since none of their generated `.tex` content changed.

**Why:** Confirms the file-path change didn't silently orphan a table (e.g., a stale `\input` pointing at a deleted file
would be a hard LaTeX error, not a silent gap — so this is lower-risk than a typical content check, logged mainly for
completeness).

Status: **CONFIRMED** (2026-09-06, user checked rendered PDF)

---

### 2026-09-04 — Third panel ($r_\mathrm{ph}$ vs $z$) added to Figure 10 (`fig:photospheric`)

**What changed:**

- `GRBResearchWork/codes-for-paper/photospheric_radius/pe_er_photosphere.py`: `make_plot()` now
  plots $r_0 (z)$, $\Gamma (z)$, **and $r_\mathrm{ph} (z)$** (was 2 panels, now 3) — same per-episode
  curve/marker/error-band treatment as the existing two panels, `r_ph` was already computed and in the CSV, just not
  plotted. `figsize` widened `(12.5,5.0)`→`(18.0,5.0)`; legend-column layout fraction rescaled proportionally so its
  absolute width is unchanged.
- Regenerated `pe_er_photosphere.png/.pdf` — **CSV is byte-identical to before** (`git diff --stat` on the `.csv` shows
  no changes), confirming this is a plot-only change, no new/different computation.
- Synced into `GRBResearchPaper/images/section5/pe_er_photosphere.png` via `codes-for-paper/sync_paper_assets.py` (the
  project's registered-asset sync tool, not a hand-copy).
- `section-5-data-analysis.tex`, `fig:photospheric` caption: "Base radius $r_0$ (left) and bulk Lorentz factor $\Gamma$
  (right)" → names all three panels, left/middle/right.
- **Also synced in the same pass** (pre-existing, already-made changes in `GRBResearchWork` from a prior session, just
  not yet copied over — not something this session authored): seed-citation additions to `tab:eiso`, `tab:lorentz`,
  `tab:lorentz_limit_b` captions and matching `tab:seed_table` cross-reference updates. Confirmed benign by diff before
  accepting (each is a 1-line addition of "(seed N)" to an existing sentence, or a `\cref` list extension) — flagged
  here for transparency since `sync_paper_assets.py` syncs everything registered, not just what this session touched.

**What to check:**

1. Figure 10 (`fig:photospheric`) now shows three panels side by side (r_0, Γ, r_ph), same visual style/legend as the
   existing two, at readable size/resolution across the full page width.
2. Caption reads "Base radius $r_0$ (left), bulk Lorentz factor $\Gamma$ (middle), and photospheric
   radius $r_\mathrm{ph}$ (right)" — grammatically clean, no leftover "left/right" pairing.
3. The new $r_\mathrm{ph}$ panel visually matches the numbers already stated in prose (the item-4 rewrite immediately
   following this figure): GRB080916C's points should sit near the bottom of the swept bursts' curves at low $z$ and
   clearly below them at $z=2$, visually demonstrating the "offset shrinks toward overlap at low $z$" claim.
4. The three unrelated 1-line seed-citation additions (`tab:eiso`, `tab:lorentz`, `tab:lorentz_limit_b`) render as clean
   trailing additions, not awkward insertions.

**Why:** This closes a real gap the user identified directly — the paper argues extensively about $r_\mathrm{ph} (z)$
behavior (item 4 above) without ever showing $r_\mathrm{ph} (z)$ itself, only $r_0 (z)$ and $\Gamma (z)$. A bare
`pdflatex` pass compiled clean (0 errors, 29 pages, unchanged from the item-4 edit) — confirms LaTeX/image-path validity
only; the figure's visual legibility at 3-panel width needs the user's own look, hence PENDING.

Status: **CONFIRMED** (2026-09-06, user checked rendered PDF)

---

### 2026-09-04 — r_ph "mild" z-dependence contradiction fixed with explicit z=0.5/z=5 numbers (quick-fixes item 4)

**What changed:** `section-5-data-analysis.tex`, §5.2 (`subsec:photospheric`):

- Line ~320: dropped "mild" (which contradicted the very next paragraph's order-of-magnitude offset claim); added a
  sentence quantifying the actual swing — $r_\mathrm{ph}$ varies by a uniform factor of $\approx4.2$ across the full
  swept range ($z=0.5$–$5.0$) for every episode, verified directly from `pe_er_photosphere.csv` (identical ratio,
  4.2032, for all 8 swept episodes across GRB131014A/GRB140206B/GRB231129C — not an estimate).
- The offset paragraph (previously 3 sentences, now 4): added explicit combined $r_\mathrm{ph}$ ranges for
  GRB131014A+GRB231129C at $z=0.5$ ($1.1$–$2.6\times10^{12}$~cm) and $z=5$ ($4.8$–$10.7\times10^{12}$~cm), computed from
  the same CSV (min/max across all thermally-detected episodes of both bursts at each z, not just T90). At $z=0.5$ this
  range **overlaps** GRB080916C's ($6.1$–$14\times10^{11}$~cm) rather than sitting an order of magnitude above it —
  turning the previous unsupported assertion ("should not be over-interpreted... scales with the assumed redshift") into
  a demonstrated result.
- Clarified that GRB140206B ($\approx7\times10^{11}$~cm) is *also* at an assumed fiducial $z=2$, same as
  GRB131014A/GRB231129C — it's grouped with GRB080916C by where its r_ph *lands* (in the
  literature's $10^{11}$–$10^{12}$~cm range), not because its redshift is measured. This was a point the user
  specifically asked to have clarified before applying the fix.

**What to check:**

1. §5.2, redshift-dependence paragraph: confirm it no longer says "mild," and instead states the $\approx4.2\times$
   swing factor.
2. §5.2, offset paragraph: confirm it now reads as 4 sentences with the two new numeric ranges
   ($1.1$–$2.6\times10^{12}$~cm at $z=0.5$; $4.8$–$10.7\times10^{12}$~cm at $z=5$), and that the conclusion ("should not
   be read as evidence of an intrinsically larger photosphere... a genuine comparison awaits a measured redshift") reads
   as a demonstrated finding rather than a bare assertion.
3. No other passage in the paper still calls this z-dependence "mild" (grep confirmed zero remaining occurrences of
   "mild" project-wide as of this edit).

**Why:** This resolves a direct self-contradiction the reviewer flagged (§6.6 — actually §5.2 — calls the dependence
"mild" one sentence before attributing an order-of-magnitude offset to that same dependence). A bare `pdflatex` pass
compiled clean (0 errors, 29 pages, unchanged from the item-1 edit) — confirms LaTeX validity only, not that the numbers
read correctly in context, hence this entry.

Status: **CONFIRMED** (2026-09-06, user checked rendered PDF)

---

### 2026-09-04 — GRB231129C stepwise-selection sensitivity check (quick-fixes item 1)

**What changed:**

- `section-4-joint-analysis-results.tex`, §4.2 (`subsec:best-models`), after the stepwise sentence: three new sentences
  stating that the BASE-model choice is never revisited after \bb\ is added — only the preferred BASE's own \bb\
  extension is tested — framed as a deliberate two-stage working-group convention (fix the continuum family first, then
  ask only whether \bb\ improves that continuum), not an oversight. Cross-references the sensitivity check in §5.2.
- `section-5-data-analysis.tex`, §5.2 (`subsec:bb-fraction`), after the existing band-robustness sentence and before the
  GRB131014A paragraph: new paragraph reporting the diagnostic recomputation of $f_\mathrm{BB}$ using `\bandbb` instead
  of the selected `\sbplbb` for GRB231129C's T90/TR1/EX0 (the only three episodes behind this burst's headline range).
  `\bandbb` has a *lower* C-stat than the selected `\sbplbb` in all three (by 19/9/8 units — visible in
  `tab:GRB231129C-large`, unchanged), even though `\band` wasn't the preferred BASE model in any of them. Reports the
  resulting alternative range ($f_\mathrm{BB}\approx0.14$–$0.19$ obs, $\approx0.14$–$0.20$ rest) and states plainly that
  it overlaps GRB131014A's TR2 outlier ($f_\mathrm{BB}=0.127\pm0.013$) within $1\sigma$ at its lowest point
  (T90, $0.142^{+0.016}_{-0.014}$) — the ranking holds at the median, not sharply.
- **No data, CSV, table, or figure changes.** `\sbplbb` remains the officially reported BEST model and every number
  elsewhere in the paper (abstract, headline "3–4×" claim, `tab:bbfraction`, `fig:bbfraction`) is unchanged — this is a
  caveat/sensitivity-check addition only, not a numbers swap. The BAND_BB alternative numbers quoted above were computed
  directly from `results.json`'s existing `BAND_BB` fit for these three episodes via `bb_flux_fraction.py`'s own
  `compute_fraction()` (same MC machinery, same seed) in a one-off diagnostic run — not written to any script or CSV in
  the repo.

**What to check:**

1. §4.2: confirm the three new sentences read as a deliberate methodological choice, not a hedge, and that the
   cross-reference to §5.2 resolves to the correct section/subsection number.
2. §5.2: confirm the new paragraph sits between the "ranking is robust to comparison band" sentence and the GRB131014A
   paragraph, reads clearly, and the two numeric ranges (0.14–0.19 obs, 0.14–0.20 rest) don't visually clash with the
   primary 0.21–0.25 range stated one paragraph earlier — i.e., it should read as an explicit caveat, not a
   contradiction.
3. `\cref{tab:GRB231129C-large}` and `\grbthirteentenfourteenA's TR2` cross-reference resolve correctly (no "??").

**Why:** This directly addresses a referee-facing weakness in the paper's headline claim — the appendix table already
shows Band+BB beating the selected SBPL+BB by up to 19 C-stat units for the exact three episodes the "most thermally
dominated burst" claim rests on. A bare `pdflatex` pass compiled clean (0 errors, 29 pages, was 28) — that confirms
LaTeX validity, not that the new paragraphs read well or land in the intended spot, hence this entry.

Status: **CONFIRMED** (2026-09-06, user checked rendered PDF)

---

### 2026-09-04 — kT rejection criterion (§4.2) + 2SBPL added to Conclusions future work (quick-fixes items 2, 5)

**What changed:**

- `section-4-joint-analysis-results.tex`, end of §4.2 (`subsec:best-models`), after the Protassov caveat paragraph: new
  sentence stating the physical rejection criterion for unphysically low blackbody temperatures ($kT \lesssim 15$\,keV,
  within a factor of a few of the \SI{10}{\kilo\electronvolt} NaI lower bound), cross-referencing
  `sec:data-preparation`. Motivated by two examples already in the appendix tables that satisfy $\Delta C\text{-stat}$
  /SAFE but were correctly excluded with no stated basis: GRB140206B Episode IV (`\sbplbb` $kT=5.11$\,keV,
  `\bandbb` $kT=4.97$\,keV) and GRB131014A EX1 (`\sbplbb` $kT=9.67$\,keV).
- `section-7-conclusion.tex`, future-work `enumerate` list: new 4th item citing the 2SBPL/synchrotron-break alternative
  (`Oganesyan2018`, `Ravasio2018`, `Burgess2020` — already cited elsewhere in the paper, no new bib entries needed) as
  unresolved, matching the Intro/§6 acknowledgment that it's a genuine untested competing explanation.

**What to check:**

1. §4.2, after "...following their use for the same BB-detection problem in Guiriec (2011) and Fana Dirirsa (2019)": one
   new sentence about the $kT \lesssim 15$\,keV rejection floor. Confirm it reads clearly as an additional, independent
   criterion (not folded into the ΔC-stat/SAFE discussion).
2. Conclusions (§7), future-work list: now 4 items, not 3. Confirm the new 4th item (2SBPL/synchrotron alternative)
   reads consistently with the existing Intro §1 and Discussion §6 mentions of the same citations, and doesn't duplicate
   their wording verbatim.

**Why:** Both are prose-only additions with no data/figure/table change. A bare `pdflatex` pass compiled clean (0
errors, 28 pages, only a pre-existing "multiply-defined labels" warning unrelated to these edits) — that confirms LaTeX
validity, not that the sentences read well or land in the right place, hence this entry.

**Source-checked 2026-09-06 (not a PDF read):** both additions confirmed present verbatim.
`section-4-joint-analysis-results.tex:37` has the $kT \lesssim 15$\,keV rejection sentence, reading clearly as an
independent criterion after the Protassov caveat paragraph, not folded into the ΔC-stat/SAFE discussion.
`section-7-conclusion.tex:39-40` has the new 4th future-work item citing `Oganesyan2018`/`Ravasio2018`/`Burgess2020`,
consistent with (not duplicating) the Intro/§6 mentions of the same citations.

Status: **CONFIRMED** (2026-09-06, source read; visual PDF layout not separately checked but content matches exactly)

---

### 2026-09-03 — RNG-seeding paper integration: paragraph, caption seeds, appendix table

**What changed:**

- `tex_files/section-5-data-analysis.tex`: two new sentences after the existing cosmology sentence (right before the $k$
  -correction equation), introducing the per-script deterministic seeding scheme and pointing to the new appendix table.
- Seven figure captions gained a trailing "Monte Carlo seed: N (\cref{tab:seed_table})" sentence:
  `fig:peakenergybestall` (2933599809), `fig:amati_relationship` (3271127181), `fig:bbfraction` and
  `fig:bbfraction-rest-vs-z` (668559939, same script), `fig:photospheric` (3315768413), `fig:butterfly` in
  `section-4-joint-analysis-results.tex` (2000318690). `fig:gamma_comparison` instead got a sentence noting its three
  constituent scripts each have an independent seed, listed in the table (it plots pre-computed CSVs from three other
  scripts, so no single seed applies to it).
- Both `gbm_only_refit` figures (`fig:gbmonly-kt`, `fig:gbmonly-cstat` in §6) were deliberately **not** given a seed
  citation — checked the actual plotting code, not just the caption: `kt_bb_err_keV` is a raw RMFIT fit-covariance error
  and $\Delta$C-stat is a deterministic fit-statistic difference, so neither figure actually plots an MC-derived
  quantity, even though the producing script does MC elsewhere (for content not shown as a figure in this paper).
- New generator `GRBResearchWork/codes-for-paper/seed_table_to_latex.py` (root-level, not inside a topic folder, at the
  user's direction — the table spans every topic folder's scripts) computes every seed live via `seed_from_name` rather
  than hardcoding, and writes `seed_table.csv` + `seed_table.tex`.
- New appendix: `main.tex` gained a third appendix section, "Monte Carlo seeds" (`\label{sec:seed-table}`), `\input`-ing
  a new thin wrapper `appendices/appendix_seeding_table.tex`, which in turn `\input`s the generated
  `tex_files/generated/seed_table.tex`.

**Correction (2026-09-03, source-checked):** the table is **9 rows**, not 11 as originally logged here.
`codes-for-paper/seed_registry.yaml` has 11 total script entries, but 2 sit under separate `unused:`/`dormant:`
top-level keys (`kt_epeak_correlation` — "not referenced by any figure or table in the current draft";
`norris_pulse_fit_tv` — the in-progress Phase 5 variability-timescale work, "no caller anywhere in codes-for-paper/"),
not the main mapping the generator iterates. Since neither has a real `cref` to print, the generator correctly excludes
them rather than emitting rows with blank Fig./Table cells. The "11-row, for completeness" framing in the original
bullet above was never actually implemented this way — this is the current, correct behavior, not a bug.

**What to check:**

1. The new appendix section ("Monte Carlo seeds") renders as a third appendix, after "GRB parameters", with a visible
   9-row table spanning the full page width. **CONFIRMED by user (2026-09-04)**.
2. Every `\cref{tab:seed_table}` in the body/captions resolves to the correct appendix table number, not "??" or a wrong
   section. **CONFIRMED** via `out/main.log` (0 undefined references/citations after the final pass) and
   `out/main.aux` — source-checkable, no PDF read needed.
3. The seven caption additions read naturally as a trailing sentence, not awkwardly tacked on. **CONFIRMED** by reading
   the source directly (`section-5-data-analysis.tex` lines 30, 147, 249, 259, 310;
   `section-4-joint-analysis-results.tex` line 58) — all seven read as clean trailing sentences.
4. `fig:gamma_comparison`'s new sentence reads clearly against the existing caption text above it. **CONFIRMED** by
   source read (`section-5-data-analysis.tex` line 221).
5. The two `\ac{GBM}` and multi-label `\cref{fig:bbfraction,fig:bbfraction-rest-vs-z}`-style entries render as expected.
   **CONFIRMED by user (2026-09-04)**, on top of the source-level reasoning above (GBM already expanded earlier in the
   body; multi-figure cref is standard, already-used syntax).

**Why:** All items are now confirmed — text/reference checks from source plus `out/main.log`/`out/main.aux`, and the
remaining layout items are confirmed directly by the user against their completed rebuild.

Status: **CONFIRMED** (2026-09-04, source + user visual confirmation)

---

### 2026-08-31 — Protassov caveat for ΔC-stat threshold (weakness #4)

**What changed:** `GRBResearchPaper/tex_files/section-4-joint-analysis-results.tex`, end of §4.2
(`subsec:best-models`) — two new sentences after the existing ΔC-stat threshold paragraph. New bib entry `Protassov2002`
added to `GRBResearchPaper/ref.bib`.

**What to check:** §4.2 "BEST models", after the paragraph ending "...requiring ΔC-stat ≥ 28.74." — confirm two new
sentences citing Protassov et al. (2002) and Guiriec (2011)/Fana Dirirsa (2019) as precedent. Confirm
`\citet{Protassov2002}` resolves.

**CONFIRMED (2026-09-03, source + build-artifact check):** `section-4-joint-analysis-results.tex` lines 30–31 have
exactly these two sentences, immediately after the Δk=2/Δk=4 threshold paragraph (line 28) — "...requiring ΔC-stat ≥
28.74" appears at line 28, matching the described anchor point closely (one paragraph break earlier than described, same
location in practice). `Protassov2002` appears as a proper `\bibitem[{{Protassov} {et~al.}(2002)...}]{Protassov2002}`
entry in the rebuilt `out/main.bbl`, and `out/main.log` shows 0 undefined citations after the final pass — resolves
correctly. No visual-only component to this entry (no figure/table).

Status: **CONFIRMED** (2026-09-03)

---

### 2026-08-31 — Rest-frame f_BB, revised "3-4x" claim, new table/figure (weakness #2)

**What changed:**

- `GRBResearchPaper/tex_files/section-5-data-analysis.tex`, §5.2 (`subsec:bb-fraction`): rewrote the $f_\mathrm{BB}$
  definition paragraph to add the rest-frame definition and z-sweep description; added a new figure block
  (`fig:bbfraction-rest-vs-z`, `images/section5/bb_flux_fraction_rest_vs_z.png`); rewrote the "most thermally dominated"
  paragraph (previously one sentence, now three) to report both observer- and rest-frame multipliers and the robustness
  statement.
- `GRBResearchPaper/tex_files/generated/bb_fraction_table.tex`: regenerated — now `table*` (was `table`), with new $z$
  and $f_\mathrm{BB}^\mathrm{rest}$ columns and a dagger footnote for fiducial-redshift rows.
- `GRBResearchPaper/images/section5/bb_flux_fraction.png/.pdf`: regenerated (same figure, now deduplicated to one point
  per episode rather than one per swept z internally — visually should be unchanged from before).
- `GRBResearchPaper/images/section5/bb_flux_fraction_rest_vs_z.png/.pdf`: new figure.
- `GRBResearchPaper/tex_files/abstract.tex`, `section-6-discussion.tex`, `section-7-conclusion.tex`: revised sentences
  reporting both bands' numbers (abstract/conclusion) or a robustness-confirmation sentence (discussion).

**What to check:**

1. §5.2 renders two figures now (`fig:bbfraction`, unchanged in appearance, and the new `fig:bbfraction-rest-vs-z` — a
   redshift-swept curve plot with three burst panels/color groups, dotted vertical line at z=2). **CONFIRMED by user (
   2026-09-04)**.
2. Table 4 (`tab:bbfraction`) now has eight columns (Model, Episode, z, kT, F_BB, F_total, f_BB^obs, f_BB^rest) spanning
   the full-text width (`table*`), with a `†` marking fiducial-z rows and a footnote explaining it. **CONFIRMED** by
   reading `tex_files/generated/bb_fraction_table.tex` directly — exactly 8 columns, `table*`, `\dagger` footnote
   present on every fiducial-z row.
3. §5.2's "most thermally dominated" paragraph: confirm it reads as three sentences ending "...raising its rest-frame
   thermal fraction more than GRB231129C's." **CONFIRMED** — `section-5-data-analysis.tex` lines 269–272, exactly three
   sentences, correct ending.
4. Abstract: confirm it now says "...three to four times... in a fixed observer-frame band... and roughly two to three
   times in a fixed rest-frame band..." **CONFIRMED** — `abstract.tex` line 9, exact wording present.
5. Conclusion (§7): same check as the abstract. **CONFIRMED** — `section-7-conclusion.tex` line 9, same wording.
6. Discussion (§6, end of the "strong/persistent" paragraph): confirm one new sentence about the ranking being unchanged
   under a rest-frame band. **CONFIRMED** — `section-6-discussion.tex` line 12.
7. No leftover reference to the old single-sentence "3-4x" wording anywhere. **CONFIRMED** by grep across `tex_files/` —
   every occurrence of "three to four times the thermal fraction of" is immediately followed by "in a fixed
   observer-frame band".

**Also found while checking this pass (2026-09-03):** the `f_BB` range quoted in three places (`abstract.tex` line 8,
`section-5-data-analysis.tex` line 266, `section-7-conclusion.tex` line 8) said "$0.050$ to $0.247$", but the
regenerated `bb_fraction_table.tex` (after the deterministic-seed rerun) actually tops out at $0.2462$, i.e., $0.246$ — a
0.001 staleness from the reseeding not being back-propagated into the hardcoded prose range. **Fixed**: all three now
read "$0.050$ to $0.246$".

**Why:** This revises the paper's headline cross-burst thermal-dominance claim across four files based on a new Monte
Carlo computation — exactly the kind of content change a clean compile doesn't verify, but one that source-reading plus
checking the generated table against the prose numbers *can* verify without opening the PDF, backed by the user's direct
visual confirmation for the one item that couldn't be.

Status: **CONFIRMED** (2026-09-04, source and user visual confirmation — one stale number found and fixed along the way)

---

### 2026-08-31 — Synchrotron/multi-break alternative caveat (weakness #3)

**What changed:** `GRBResearchPaper/tex_files/section-1-introduction.tex` (new paragraph after the photospheric-evidence
literature review, before the Amati-relation paragraph); `GRBResearchPaper/tex_files/section-6-discussion.tex` (one new
sentence after the "line of death" sentence in the GRB131014A/GRB231129C discussion). New bib entries `Oganesyan2018`,
`Ravasio2018`, `Burgess2020` in `ref.bib`.

**What to check:**

1. Introduction: confirm the new 5-sentence paragraph appears, citing all three new references, and that it explicitly
   states the synchrotron-break model was *not* fit against this paper's data. **CONFIRMED** —
   `section-1-introduction.tex` lines 35–39, exactly 5 sentences, cites all three, and line 39 explicitly disclaims
   fitting a 2SBPL/synchrotron-break model against this paper's own data.
2. Discussion §6: confirm the one-sentence caveat appears immediately after "...line of death, which further reduces the
   spectral leverage available to separate the two components." **CONFIRMED** — `section-6-discussion.tex` line 27,
   immediately follows.
3. All three new citations resolve correctly (author/year, not "?"). **CONFIRMED** — all three (`Oganesyan2018`,
   `Ravasio2018`, `Burgess2020`) appear as proper `\bibitem` entries in the rebuilt `out/main.bbl`, and `out/main.log`
   shows zero undefined citations after the final pass.

**Why:** Source read plus the rebuilt `.bbl`/`.log` (from your non-draft rebuild, 2026-09-03) confirm both the text and
the citation resolution without needing a PDF read. Nothing purely visual to check in this entry — no new figure or
table.

Status: **CONFIRMED** (2026-09-03, source + build-artifact check)

---

### 2026-08-31 — Supplementary-PL clarification (weakness #5)

**What changed:** `GRBResearchPaper/tex_files/section-4-joint-analysis-results.tex`, §4 opening paragraph — one
ambiguous sentence expanded into two: the Fermitools origin of the supplementary power-law term, and the statistical
reason BASE+PL is not tested as an independent branch (it's not nested with BASE+BB, so a direct likelihood-ratio
comparison between them isn't well-defined; a PL component is instead only added as a nested extension on top of an
accepted BASE+BB). Kept deliberately brief per the user's instruction — no project history (an earlier, pre-repo
symmetric-grid exploration) in the paper text; that stays in `review-resolution.md`/`HANDOFF.md` only.

**What to check:** §4 "Joint Analysis" opening paragraph — confirm the two-sentence version reads clearly and matches
the already-existing stepwise-procedure text a few lines below it (§4.2 "BEST models": BASE→BASE+BB at Δk=2,
BASE→BASE+BB+PL at Δk=4, then the stepwise BASE+BB→BASE+BB+PL at Δk=2) without contradicting or duplicating it.

**Discrepancy found (2026-09-03, source-checked), then resolved same day:** the actual sentence at
`section-4-joint-analysis-results.tex` line 10 is one sentence, not the two-sentence Fermitools-origin expansion this
entry originally described. Checking `review-resolution.md` item #5 (missed on first pass) shows this was already fully
investigated and decided on **2026-09-01**: the two-sentence version was found not to have survived into the file (cause
undetermined — possibly a manual edit), and the user explicitly chose to close it as-is rather than restore the fuller
wording, because the complete explanation made the paragraph "quite complex." The reasoning itself is preserved in
`review-resolution.md` #5, items 1–4, and reconfirmed directly by the user on 2026-09-04: (1) the supplementary PL
term's origin is `fermitools`' own LAT-band power-law characterization, used to make LAT data RMFIT-compatible; (2)
BASE+PL was deliberately never tested as an independent branch — BASE+PL and BASE+BB both add 2 dof to BASE along
different, non-hierarchical branches, so they aren't directly comparable via a likelihood-ratio test the way
BASE→BASE+BB or BASE+BB→BASE+BB+PL are; only the latter nested chain (same root/stem throughout) is statistically
well-defined, so that's the one kept. The user will handle any reviewer pushback on this point directly rather than have
the paper text carry the full reasoning.

**Why:** Content and scope both resolved — this was never a live gap, just a VERIFY.md entry written without checking
the prior decision already on record.

Status: **CONFIRMED** (2026-09-04 — decision already final as of 2026-09-01, per `review-resolution.md` #5; no
paper-text change needed)

---

### 2026-09-01 — Priority 3 editorial/consistency pass (all 8 items)

**What changed:**

- `section-1-introduction.tex`: "three criteria" → "two criteria"; "25 units" → "28.74 units"; new sentence +
  `\input{tex_files/generated/fluence_table}` (new `tab:fluence`) right after the two selection criteria.
- `section-3-data-preparation-and-analysis.tex`: four "10°" ROI mentions (lines 88, 96, 106, 115) → "12°", matching the
  12° already stated at line 72 and every `LAT_analysis/*/fit_results_*.txt`.
- `section-4-joint-analysis-results.tex`: `\sbpl` citation now `\citep[\sbpl]{Ryde1999, Kaneko2006}` (was `Ryde1999`
  alone); "credible interval"/"posterior samples" language (lines 56, 65, 66) replaced with "Monte Carlo-propagated
  interval"/"Monte Carlo parameter samples".
- `section-5-data-analysis.tex`: the Γ_min "systematically weaker" sentence rewritten to note the bound isn't strictly
  monotonic (spectral-index dependence can offset the $t_v$ effect), citing GRB080916C's T90-vs-TR3 near-tie
  in `tab:lorentz`.
- `tex_files/generated/amati_relationship_table.tex`, `bb_fraction_table.tex`, `photospheric_table.tex`: regenerated —
  caption-text-only changes (MCMC→Monte Carlo wording; kT-provenance cross-reference between the latter two), no data
  values changed.
- `tex_files/generated/fluence_table.tex`: **new** — 4-row table (one per burst, T90 only), from a corrected
  `grb_fluence.py` run (8 keV–40 MeV band, was silently 10 keV–1 MeV).

**What to check:**

1. Introduction: "selected to satisfy two criteria simultaneously" (not "three"); "$\Delta C$-stat criterion of 28.74
   units" (not "25"); new `tab:fluence` right after the selection-criteria list. **CONFIRMED** —
   `section-1-introduction.tex` lines 47, 60, 53–57. Fluence values now read 1.6643/1.7374/2.3808/0.7852 (were
   1.6644/1.7382/2.3810/0.7854) — a ~0.01% shift from the deterministic-seed rerun reviewed and committed this session,
   not a content change.
2. §3 (data preparation): all four per-burst "ROI" mentions now read 12°. **CONFIRMED** —
   `section-3-data-preparation-and-analysis.tex` lines 88, 96, 106, 115, all 12°, matching the definition at line 72.
3. §4 (Joint Analysis): SBPL citation now shows two references; no "credible interval"/"posterior samples" anywhere.
   **CONFIRMED** — line 8 has `\citep[\sbpl]{Ryde1999, Kaneko2006}`; project-wide grep across `tex_files/` finds zero
   remaining "credible interval", "posterior samples", or "MCMC".
4. §5 (data analysis), Γ_min paragraph: confirm it doesn't contradict `tab:lorentz`'s printed values. **CONFIRMED** —
   `tex_files/generated/lorentz_table.tex` shows T90 = 507 (highest) and TR3 = 504 (its closest time-resolved
   competitor, within error bars) — the text's "comparable to, rather than weaker than" claim is accurate.
5. `tab:eiso` caption says "Monte Carlo samples", not "MCMC". `tab:bbfraction`/`tab:photospheric` kT-provenance
   cross-reference. **CONFIRMED** by direct read of all three generated `.tex` files — correct wording, cross-reference
   present.
6. `\citet{Kaneko2006}` resolves alongside `Ryde1999`. **CONFIRMED** — proper
   `\bibitem[{{Kaneko}...(2006)...}]{Kaneko2006}` entry present in the rebuilt `out/main.bbl`.

**Why:** All six checks above are text/reference-resolution checks, verifiable from source plus the rebuilt `.bbl`
/generated `.tex` files without reading the PDF — done 2026-09-03. Nothing purely visual in this entry (no new figure,
only a small 4-row table whose content is already source-confirmed).

**Also observed, resolved 2026-09-03:** the appendix-block discrepancy flagged here has been investigated and closed —
see `review-resolution.md` item #1, now updated to CLOSED. The appendix was intentionally re-activated by you directly
in commit `464a747` (2026-09-01), not left commented out; it has since grown (the seeding-table section was added on top
of it 2026-09-03).

Status: **CONFIRMED** (2026-09-03, source + build-artifact check — no visual-only items remain in this entry)

---

### 2026-09-01 — GRB080916C-vs-literature paragraph, Amati softening, N=4 overclaiming pass, one citation fix (weaknesses #8, #9, #10)

**What changed:**

- `section-6-discussion.tex`: new paragraph after the GRB080916C photospheric-radius sentence, comparing this paper's BB
  detections and $\Gamma$ values to `Abdo2009FermiObservations080916C` and `Guiriec2015`; also softened "This confirms
  that population-level studies... will systematically underestimate..." → "This illustrates how population-level
  studies... can underestimate...".
- `section-7-conclusion.tex`: same softening as above, same sentence.
- `section-5-data-analysis.tex` (§5.2 Amati subsection): softened "confirming...physically self-consistent" →
  "consistent with...physically self-consistent"; added two sentences on time-resolved-Amati validity (Basak & Rao
  2012's r=0.80→0.37→0.89 statistic; Ghirlanda et al. 2010 on the Yonetoku relation) and why this paper's Bayesian-block
  TR/EX episodes sit at the pulse scale where the relation is expected to hold.
- `abstract.tex`: same Amati softening (one sentence); "...component of prompt GRB spectra" → "...component of prompt
  emission in this sample".
- `section-1-introduction.tex`: "we test whether photospheric emission is a universal or selective feature..." →
  case-study framing without "universal"; separately, the GRB120323A sentence's citation fixed from
  `\citet{Guiriec2015}` (wrong paper — that key is GRB080916C) to `\citet{Guiriec2013}` (the real GRB120323A paper,
  newly added to `ref.bib`).
- `ref.bib`: three new entries — `Guiriec2013`, `BasakRao2012`, `Ghirlanda2010`.

**What to check:**

1. §6 Discussion, right after "...placing this burst squarely within the range reported for other Fermi GRBs with
   confirmed photospheric components": confirm a new paragraph discussing Abdo (2009) and Guiriec (2015), ending with
   the $\Gamma_{\min}=507$-vs-thermal-$\Gamma$ comparison and a resolving table cross-reference. **CONFIRMED** —
   `section-6-discussion.tex` lines 17–22, immediately follows that sentence, cites both, ends on the Γ_min=507
   comparison with `\cref{tab:lorentz}` (resolves per `out/main.log`, 0 undefined refs).
2. §5.2 Amati subsection: softened sentence plus two new sentences citing Basak & Rao (2012) and Ghirlanda et al.
   (2010). **CONFIRMED** — `section-5-data-analysis.tex` lines 154–155; both citations present as proper `\bibitem`
   entries in `out/main.bbl`.
3. Abstract: "consistent with the physical self-consistency"; "...in this sample". **CONFIRMED** — `abstract.tex` lines
   11, 13.
4. Introduction: case-study framing without "universal"; GRB120323A cites Guiriec (2013). **CONFIRMED** —
   `section-1-introduction.tex` lines 32, 63.
5. §6/§7: both copies read "illustrates how... can underestimate". **CONFIRMED** — `section-6-discussion.tex` line 70,
   `section-7-conclusion.tex` line 18.
6. Reference list: `BasakRao2012`, `Ghirlanda2010`, `Guiriec2013` correctly formatted. **CONFIRMED** — all three present
   as proper `\bibitem` entries in the rebuilt `out/main.bbl` (checked directly, not via `pdftotext`).

**Why:** All six items above are text/citation checks, fully verifiable from source plus the rebuilt `.bbl`/`.log`
without a PDF read — redone 2026-09-03 against the fresh non-draft build (the original `pdftotext` spot-check mentioned
above predates that rebuild, so this source check supersedes it). Nothing purely visual in this entry — no new figure or
table.

Status: **CONFIRMED** (2026-09-03, source + build-artifact check)

---

### 2026-09-01 — GRB131014A off-axis-angle robustness write-up (weakness #6)

**What changed:**

- `section-6-discussion.tex`: five new paragraphs inserted after the existing GRB131014A/GRB231129C
  shallow-shoulder/synchrotron-alternative paragraph, before the GRB140206B paragraph. States the off-axis coincidence
  (70° off-axis for the entire T90, same burst where every episode is BB-augmented) plainly; rules out the direct
  LAT-effective-area mechanism (LAT sensitive only above 100 MeV vs. this burst's kT_BB ≈ 27–44 keV); introduces the
  GBM-only refit as the test for the remaining indirect joint-likelihood pathway; reports the result (same model wins,
  kT_BB consistent within 1σ, Δ-C-stat clears threshold by a wide margin and is actually *higher* without LAT in every
  episode).
- Two new figures embedded: `fig:gbmonly-kt` (kT_BB, joint vs. GBM-only) and `fig:gbmonly-cstat` (Δ-C-stat vs. the 28.74
  threshold, log scale), both `\citep{Atwood2009}`-adjacent text but no new bib entries — `Atwood2009` was already cited
  elsewhere in the paper for the LAT's 100 MeV threshold, reused here.
- New folder `images/section6/` (first Discussion-section figures in this paper) holding
  `gbm_only_refit_kt_comparison.pdf/.png` and `gbm_only_refit_delta_cstat.pdf/.png`, copied from
  `GRBResearchWork/codes-for-paper/gbm_only_refit/`.

**What to check:**

1. §6 Discussion, right after "...a competing explanation we have not tested against our fits": confirm five new
   paragraphs, covering the off-axis coincidence, the energy-band rebuttal, the indirect-pathway caveat, the GBM-only
   test description, and the result — before the GRB140206B paragraph starts. **CONFIRMED** — `section-6-discussion.tex`
   lines 29–30, 31, 32, 34–35, 55–57 are exactly these five paragraphs in this order, immediately after line 27's caveat
   sentence and before line 59's GRB140206B paragraph.
2. Both new figures render legibly at column width, and whether the "Fit type"/"Episode" legends are actually illegible
   or just tight. **CONFIRMED by user (2026-09-04)** — both figures render, legends fine.
3. `\Cref{fig:gbmonly-kt}` and `\cref{fig:gbmonly-cstat}` resolve to actual figure numbers. **CONFIRMED** via
   `out/main.log` (0 undefined references after the final pass).
4. The `\qty{100}{\mega\electronvolt}` and `\SI{70}{\degree}` values render correctly. **CONFIRMED by user (
   2026-09-04)**, plus macro usage confirmed present in source (lines 29, 31).
5. Numbers in text match the figures. **CONFIRMED** — text values verified from source (line 55–56); user confirmed
   2026-09-04 that §7.2's LAT-photon-correction writing is proper against the rendered figures.

**Why:** All items now confirmed — text/reference checks from source plus `out/main.log`, and the remaining visual items
(legend legibility, figure rendering) confirmed directly by the user against their completed rebuild.

Status: **CONFIRMED** (2026-09-04, source + user visual confirmation)
