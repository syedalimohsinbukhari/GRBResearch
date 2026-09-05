# Time-integrated / time-resolved duration table — method notes

Companion to `time_integrated_table.py`. Written 2026-09-06, closing the last `no_known_source`
entry in `table_registry.yaml` for a table that is actually `\input` by the paper
(`tex_files/section-3-data-preparation-and-analysis.tex:33`, which pulls the generated table in
directly -- there is no intermediate wrapper file; see table_registry.yaml's header comment).

Like `LAT_analysis/`, this folder holds no physics — it is bookkeeping. Every cell is a raw
`results.json` interval boundary; nothing is fit, integrated, or run through Monte Carlo.

## 1. What the code computes

For each of the paper's four sample GRBs, `results.json[<long name>]`'s keys are interval strings
(`"TR1 1.280_4.864"`, etc.), parsed by `TimeInterval.from_string` (`src/grb_research/grb_time.py:66`,
documented under CLAUDE.md's "Episode naming"). The table is a straight transcription of every
interval's `(start, end)` into `\sirangeDuration{start}{end}`, laid out with GRBs as columns and
episodes as rows — `T90` in the time-integrated block, `EX0`, `TR1..TRn`, `EX1` in the time-resolved
block. A GRB missing a given row's episode gets `\tabledash`.

## 2. The problem this solves

Before this script, this table was hand-typed directly into `appendices/time-integrated-table.tex`
with no generator anywhere in the repo — flagged in `table_registry.yaml`'s `no_known_source`
section. Regenerating it against `results.json` surfaced a real transcription error (see §4): the
same failure mode as BUG-10 and the reason `LAT_analysis/csv_to_latex.py` exists at all. That
appendix file was later deleted (2026-09-06) once the section file was pointed at the generated
table directly — see §3's wrapper-removal note.

## 3. Judgement calls

- **Row lookup key vs. display label are kept separate** *(Claude, caught in review before this was
  ever synced)*. The time-resolved block displays bare digits (`1`, `2`, …) for TR episodes, not
  `TR1`/`TR2` — that's the pre-existing table's own convention, kept as-is. A first draft of this
  script used the bare digit as the dict lookup key too, so every TR row silently rendered as all
  `\tabledash`; caught by diffing the generator's output against the previously-published table
  before syncing anything, not by inspection.
- **GRB names in the header use the paper's `\grbxxxxxxx` macros**, not literal `GRB080916C` text
  *(Claude, changed from the previous hand-typed table)*. The macros (`tex_files/preamble.tex:128-131`)
  expand to exactly that literal text, so this is a zero-visual-diff change — done for consistency
  with every other generated table (LAT, Amati, lorentz), which all use the macros.
- **The commented-out "TR breakdown episodes" rows (`BR--A/B/C`, GRB131014A) are dropped, not
  reproduced** *(Claude, confirm with user if these should come back)*. Their values
  (`1.344_6.784`, `6.784_10.688`, `10.688_21.824`) don't appear anywhere in `results.json` —
  not under `GRB131014215` nor `GRB131014215GBM` — so they aren't derivable from the data this
  script reads. They were already commented out in the hand-typed table (never rendered), and
  per CLAUDE.md's "Episode naming", `BRn`/`SP` intervals aren't used by any of this paper's four
  bursts anyway, so dropping the dead comment block changes nothing about what compiles.
- **No wrapper file — the section `\input`s the generated table directly** *(user's own edit,
  2026-09-06; extended to the LAT/seed tables and documented in `table_registry.yaml`'s header
  comment)*. `appendices/time-integrated-table.tex` had only ever held a comment plus one
  `\input{tex_files/generated/time_integrated_table}` line — pure indirection, and this table was
  never a real appendix item anyway (it's `\input` from Section 3's body, unlike the LAT/seed
  tables which genuinely live in `main.tex`'s Appendix). Deleted; `section-3-...tex` now `\input`s
  the generated path directly, matching how amati/lorentz/photospheric/fluence already work.

## 4. Validation performed

Generator output was diffed cell-by-cell against the previously-published
`appendices/time-integrated-table.tex` (since deleted, see §3) before syncing. Every cell matched
except one:

- **GRB231129C's `EX0` and `TR1` rows both showed `1.792` as the end time; the correct value is
  `3.136`.** `results.json["GRB231129779"]` has `"EX0 -0.192_3.136"` and `"TR1 0.384_3.136"` — and
  per CLAUDE.md's own definition, `EX0` "shares `TR1`'s end time" by construction, so both cells
  must agree, and neither should be `1.792` (that number doesn't appear in any interval for this
  GRB in `results.json`). This was a genuine transcription error in the previously-published table,
  caught only because the generator's output disagreed with it — not found any other way. Not yet
  confirmed in the rendered PDF; see `VERIFY.md`.

No other discrepancies. Row/column layout, `T90` and `EX1` values, and the `\tabledash` placements
for episodes a GRB doesn't have all matched the previous table exactly.

## 5. Known limitations

- **No MC, no cosmology, no seed** — same reasoning as `LAT_analysis.md` §4: every number is read
  straight from `results.json`, so there is nothing to seed and no `n_samples`/`seed`/`H0`/`Om0`
  columns apply.
- **The copy into `GRBResearchPaper/tex_files/generated/` is manual** (`sync_paper_assets.py`), same
  as every other generated table — the two repos are independent (BUG-10).
