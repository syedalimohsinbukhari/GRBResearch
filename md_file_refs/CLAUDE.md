# Scope — read this before touching the filesystem

This is a two-repo project split across sibling directories under `~/Pictures/grb_research/`:

- **`GRBResearchPaper`** — the LaTeX manuscript. `tex_files/` for sections, `appendices/` for per-GRB tables, `images/` for figures, `bibtex/`, `Literature Review/` (background reading for the paper — `literature review.pdf` and `literature review.docx`).
- **`GRBResearchWork`** (sibling repo, `../GRBResearchWork`) — the Python analysis code and data that produce the paper's numbers and figures. Key locations: `src/grb_research/` (the `grb_research` package: model/fit parsing, flux/fluence Monte Carlo calculators, styling), `codes-for-paper/` (per-figure/table scripts, organized by topic e.g. `fluence/`, `evolution_of_kt/`, `lorentz_factor/`, `amati_relationship/`), `results.json` (the fitted-model database everything reads from), and per-GRB raw-interval directories (e.g. `GRB080916009/`).

Four working documents — `PLAN.md` (the phased plan), `BUGS.md` (the bug/observation log), `HANDOFF.md` (current state — read it first), and `VERIFY.md` (pending PDF-build verification requests — see "PDF build verification" below) — are mirrored at `GRBResearchWork/md_file_refs/`. Bare references to those filenames anywhere in this project mean the `GRBResearchWork/md_file_refs/` copies — check there, not the project-root originals (`~/Pictures/grb_research/`), which are being phased out.

## Research scope — which GRBs are actually in the paper

`GRBResearchWork` holds raw-interval directories and `results.json` entries for **eight** GRBs, but this paper's sample is only **four** of them:

- GRB080916C → `GRB080916009`
- GRB131014A → `GRB131014215`
- GRB140206B → `GRB140206275`
- GRB231129C → `GRB231129779`

The other four directories present in the repo — `GRB110721200`, `GRB110731465`, `GRB150210935`, `GRB190114873` — are **not** part of this paper's sample (leftover/other-project data). Don't pull them into paper figures, tables, sample statistics, or new analysis code unless the user explicitly asks to bring one in. Existing code/scripts that still reference the old sample (e.g. `codes-for-paper/lorentz_factor/lorentz_factor.py`) are a known inconsistency to fix, not a precedent to follow.

## Episode naming — `T90`, `TRn`, `EX0`/`EX1`

Episode labels come from the interval strings in `GRBResearchWork/results.json`, parsed by `TimeInterval.from_string` (`src/grb_research/grb_time.py:66`).
They are **not** derived from the per-GRB directory names, which use a different scheme (`Ep0`, `Ep1`, `Ep1A`, …).

- **`T90`** — the time-integrated interval covering the burst's $T_{90}$.
- **`TRn`** — time-resolved Bayesian-block interval *n*, numbered from 1 in time order. Carries an index.
- **`EX0` / `EX1`** — **excess episodes**: emission falling *outside* the $T_{90}$ window, captured by extending an adjacent `TR` interval past the $T_{90}$ boundary while holding its other edge fixed.
  - `EX0` is the **leading** excess — it shares `TR1`'s end time and starts before $T_{90}$ does. GRB080916C: `EX0 -0.128_4.864` against `TR1 1.280_4.864`, with $T_{90}$ starting at 1.280.
  - `EX1` is the **trailing** excess — it shares the last covered `TR`'s start time and ends after $T_{90}$ does. GRB080916C: `EX1 59.520_67.904` against `TR5 59.520_64.256`, with $T_{90}$ ending at 64.256.
  - A burst may have either, both, or neither. Of this paper's four, only GRB140206B lacks a trailing excess (`EX0` but no `EX1`).
- **`BRn`** (enum name `SP`) — not used by any of this paper's four bursts. Note the canonical string is `BR<n>` while `episode_label()` in `lorentz_factor.py` returns the *enum name*, giving `SP<n>`; the two disagree, which matters only if a `BR` interval is ever brought in.

Three things that are easy to trip over:

- **The digit in `EX0`/`EX1` is part of the type, not an index.** They are separate `EpisodeTypes` enum members (`grb_time.py:17-25`), and `TimeInterval.__post_init__` **raises** if an `index` is passed for an `EX` interval (`grb_time.py:50-56`). Only `TR` and `SP` take an index.
- **Only `EX0` and `EX1` parse.** The regex is `(EX[01])` (`grb_time.py:46`); a third excess episode falls through to `EpisodeTypes.UNKNOWN` *silently*, with no error. Adding one means changing both the enum and the regex.
- **Map directory names to labels by interval bounds, never by name.** `Ep0` → `T90` and `EpN` → `TRn`, but `EpNA` is "the excess extending `TRn`", so GRB080916C's `Ep1A` is `EX0` while its `Ep5A` is `EX1`. Deriving the `EX` number from the directory's `N` gives the wrong answer.

## Output conventions — plots & CSVs

Every plotting script must call `update_style()` (`GRBResearchWork/src/grb_research/__init__.py:58-105` — sets shared rcParams: fonts, dpi, grid, tick style; **not** `grb_styles.py`, which only holds the `GRBPlotStyle` color/marker-convention class). This is the existing project-wide consistency convention, already followed by most of `codes-for-paper/` (good examples: `amati_relationship/amati_relationship.py`, `model_parameters/high_index.py`). Do not copy `codes-for-paper/evolution_of_kt/kt_evolution_080916c.py` as a style template — it's a known broken example (no `update_style()`, no axis labels, no legend, no `savefig`, ends in a bare `plt.show()`).

Beyond style consistency, every new plot or CSV this project produces should be **maximally defensible**: self-contained enough that a reader can verify or reproduce the number without hunting through the script. Minimum parameters needed, nothing missing:

**CSVs** — every row carries:
- Identity: GRB name (paper name, e.g. `GRB080916C`), episode/interval label (e.g. `T90`, `TR1`), model name (e.g. `BAND_BB`).
- Assumptions actually used to produce that row: redshift `z` (or the swept value, when a z-range is used instead of a fixed z), cosmology (`H0`, `Om0`). The project cosmology is **`H0 = 69.6`, `Om0 = 0.286`** (Fana Dirirsa 2019) — reconciled across `tex_files/section-5-data-analysis.tex`, `grb_calculations.py::mc_e_iso_sampler` and `lorentz_factor.py`. Do not introduce a third value.
- Value + uncertainty: derived quantities as `<name>_<unit>`, `<name>_err_lower_<unit>`, `<name>_err_upper_<unit>` — follow `codes-for-paper/amati_relationship/amati_relationship.csv`'s column-naming convention (units baked into the name, asymmetric MC errors), not `flux_fluence.csv`'s (unitless columns).
- MC provenance: `n_samples`, `seed` as constant-valued columns, so the row is self-describing without a companion file that can go missing. Seed via the existing `get_rng(seed=None, rng=None)` helper (`grb_calculations.py:23-44`) and thread `seed=`/`rng=` through function boundaries — never a bare global `np.random.seed()`.
- Nothing beyond that — no intermediate arrays or debug columns; stop at what's needed to recompute or sanity-check the final number.

**Plots** — every new figure:
- Calls `update_style()`.
- Has axis labels with explicit units.
- Uses the shared `GRBPlotStyle` conventions (`src/grb_research/grb_styles.py`) for per-GRB color / per-episode marker / BB-vs-non-BB marker fill, so figures stay visually consistent with the paper's existing ones (solid = no BB, hollow = BB-augmented).
- Shows error bars/bands wherever the underlying value has an MC uncertainty.
- **Legend entries identify the burst, the episode _and_ the fitted model** — e.g. `GRB131014A TR1 (BAND_BB)`, not a bare `TR1`. Which model won for a given episode is a result in its own right and varies within a single burst (GRB131014A alone spans `BAND_BB` and `SBPL_BB`), so a legend that omits it forces the reader back to the table to interpret the figure. Where the burst is already fixed by the panel title, the episode and model still both belong in the entry. `codes-for-paper/photospheric_radius/pe_er_photosphere.py` is the reference. Whatever encodes the episode (marker, line style) must ride on the *labelled* artist, or the legend handle renders without it.
- Is saved explicitly to both `.png` and `.pdf` (per `amati_relationship.py`'s save loop) — never left as a bare `plt.show()`. Use `dpi=SAVE_DPI` rather than a literal `600`.
- **Takes every font size, line width, marker size and DPI from `grb_constants`** (`LABEL_FONT_SIZE`, `LEGEND_FONT_SIZE`, `TICK_FONT_SIZE`, `TITLE_FONT_SIZE`, `ANNOTATION_FONT_SIZE`, `LINE_WIDTH`, `MARKER_SIZE`, `SAVE_DPI`, …), which `update_style()` installs as rcParams. Never hardcode `fontsize=8`, `linewidth=1.3`, `dpi=600` and the like — that silently overrides the shared convention and defeats the single-source-of-truth design stated in `grb_constants`. If a value genuinely isn't covered by a constant (e.g. hollow-marker edge width), define it as a named module-level constant rather than inlining a literal.

## Formatting conventions

**Python** — line length is **120**, not the default 80. This is already set in `GRBResearchWork/ruff.toml` (`line-length = 120`, matching Black's `[format]` block); match it when writing or reformatting code, and don't reflow existing lines to 80.

**LaTeX** — write **one sentence per line**. Never wrap a sentence across multiple lines, and never put two sentences on one line. The existing `tex_files/section-*.tex` already follow this, and it is what keeps `git diff` readable: a reworded sentence shows as a one-line change instead of a re-flowed paragraph. This applies to generated `.tex` too — table captions emitted by helper scripts should break at sentence boundaries.

**Captions describe, they don't interpret.** A figure or table caption should let a reader parse the panel/table standalone — what's plotted, symbols, units, error definition, notation needed for a label or column (e.g. what `TR`/`EX` mean) — but never a claim about what the data show, a comparison between quantities, or a conclusion. Those belong in the body text, where they can also carry more precision than a caption should repeat (an exact ratio, not just "exceeds"). Test: if the sentence would still be true after the underlying numbers changed, it's description and stays; if a rerun of the analysis could make it false, it's interpretation and moves to the text. Caught in Phase 4 QA: `fig:gamma_comparison`'s caption in `tex_files/section-5-data-analysis.tex` originally ended with "The thermal Γ exceeds both opacity limits… and Limit A exceeds Limit B throughout this sample" — a result already stated with actual ratios in the preceding paragraph; trimmed from the caption, left in the text.

## Per-folder method notes (required)

Every folder of analysis code must carry a `<topic>.md` alongside its scripts, written for the user's own understanding rather than for the paper. This covers `codes-for-paper/<topic>/` and also data-preparation folders such as `LAT_analysis/` — see `LAT_analysis/LAT_analysis.md`, which documents a layer that computes no physics but decides units, joins and provenance. The paper reports *what* the result is; this file records *how and why*. It must cover:

- **What the code computes**, stated as the defining equation with its source (author, year, equation number), not just prose.
- **Every judgement call and who made it** — attribute explicitly: decided by the user, proposed by Claude, or settled jointly. Include the options that were rejected and the reason, so a decision is never silently re-litigated later.
- **Conventions that could plausibly have gone another way**: energy band and reference frame, observer vs rest frame, which components are included or excluded, integration limits, MC sample count and seed.
- **Validation performed** — the independent checks that were run and what they returned (analytic cross-checks, literature comparisons, agreement between two code paths). A number nobody checked should be labelled as such.
- **Known limitations and open questions**, cross-referenced to `BUGS.md` entries where relevant.

Write it as the work happens, not retrospectively. `codes-for-paper/bb_fraction/bb_fraction.md` is the reference example.

**Table values carry uncertainties** — every quantity in a generated table that has a Monte Carlo or fit uncertainty must be printed with it, as `$value^{+upper}_{-lower}$`, not as a bare median. This applies to *inputs* as much as to derived results: if $kT$, $f_\text{BB}$ or a flux has an error in the CSV, it has an error in the table. Print a bare number only when the quantity genuinely has no uncertainty (an assumed redshift, a fixed cosmology, an episode label). The CSVs have consistently been complete here — it is the LaTeX projection that keeps losing columns, so check the generated `.tex` against the CSV headers rather than assuming the helper carried everything across. Two rounds of this were caught only by reading the rendered PDF (see `BUGS.md`, BUG-17).

**Normalise table columns by a common power of ten**, carried in the header rather than repeated on every row: `$r_\text{ph}$ [$10^{11}$ cm]` with entries like `6.10`, not `$6.10\times10^{11}$` thirteen times. Normalise only where it helps — a column already of order 0.01–1000 (an index, a temperature in keV, $\Gamma$) reads better unnormalised.

**Derive decimal places numerically, in the generator, never by eye.** The rule used here: enough decimals that the *smallest error in the column* still shows two significant figures, i.e. `decimals = 1 - floor(log10(min_error/norm))`. Both `csv_to_latex.py` scripts implement this as `decimals_for()` and call it on the dataframe at generation time, so the formatting adapts if the data changes instead of silently going stale. Eyeballing produces inconsistent precision across columns and is how a column ends up showing an error of `0.00` — check the rendered numbers, not the intent.

**Generated LaTeX tables** — tables that come from data are produced by a helper script and `\input` into the paper, never retyped into a section file. `codes-for-paper/amati_relationship/csv_to_latex.py` is the reference implementation (reads the CSV, emits the paper's `\grb...`/`\sbplbb` macros); `lorentz_factor.py` emits its table inline. Generated `.tex` lands in `GRBResearchPaper/tex_files/generated/` and carries an `% AUTO-GENERATED by <script>` header; the section or appendix content that uses it `\input`s that path directly, e.g. `main.tex`'s `\input{tex_files/generated/lat_info_table}`. (Before 2026-09-06 this went through an extra hand-authored wrapper file per table, e.g. `appendices/appendix_LAT_info.tex` — just a comment plus one `\input` line; removed as pure indirection once every generated table's `dest` was already a stable, descriptive path in its own right. `codes-for-paper/table_registry.yaml`'s header comment has the details.) Copying the `.tex` between the two repos is manual, since they are independent. Hand-copying numbers into prose or tables is how the $E_\text{iso}$ table silently drifted out of sync with the pipeline (see `BUGS.md`, BUG-10).

**An input that genuinely cannot be derived must fail loudly, not sit quietly.** Where a generated table needs a number that is not reconstructable from the data in the repo — the `\pFlux` upper limits in the LAT appendix are the live example, being profile-likelihood limits from a gtlike step whose output was never saved — keep it in one clearly-named constant in the generator, key it so a stale entry cannot be mismatched silently, and **assert** it against whatever *is* derivable. `LAT_analysis/csv_to_latex.py` raises if its `FLUX_UPPER_LIMITS` dict disagrees with the episodes it derives from `TS < 25`. A crash is a recoverable failure; a silently wrong number in a published table is not (BUG-10, OBS-08).

**Units belong in the name, and in a comment giving the evidence.** A factor-1000 error in $\Gamma_\text{min}$ survived because a variable called `E_max_keV` was fed MeV, and a previous fix argued the wrong way from an arithmetic slip (BUG-16, then BUG-18). Where a physical constant pins the unit — `511` keV versus `0.511` MeV for $m_ec^2$ — say so at the line that uses it, and record *why* the unit is what it is, not just what it is.

## Cross-folder imports — copy rather than fight `sys.path`

If a script needs a function that lives in another folder within the same repo and importing it cleanly would mean `sys.path` manipulation, relative-import gymnastics, or restructuring either folder as a package — don't. Copy the function verbatim into the new folder instead, with a one-line docstring/comment noting where it came from and that it's a copy, not a re-derivation (e.g. `light_curves.py`'s `lightcurve_data()`, copied into `variability_timescale/` on 2026-08-22). Then leave a note in that folder's `<topic>.md` (the per-folder method note already required — see "Per-folder method notes" below) flagging the duplication so the user can rework it into a proper shared import later if they want to. A working duplicate beats a fragile import path — this is a call to keep moving, not a permanent architecture decision.

## PDF build verification — delegated to the user

Running `pdflatex`/`bibtex` and especially `pdftotext | grep` through the Bash tool pulls build logs and rendered-PDF text into context — real token cost for what is fundamentally a sanity check. **Claude does not run the LaTeX build or `pdftotext` to verify content.** When a change needs confirming against the rendered PDF, append an entry to `VERIFY.md` (`GRBResearchWork/md_file_refs/`) instead and stop there.

**Don't log an entry after every individual edit.** Keep editing — through a whole subsection, a whole set of related sentence rewords, a whole multi-file consistency pass — and batch the resulting checks into one `VERIFY.md` entry (or a few, if the edits genuinely touch unrelated claims) once that unit of work is done. Ping the user at that point, not per-edit. The trigger is "a substantial piece of work just landed" or "something crucial needs checking" — a single word fix or a one-line rewording doesn't clear that bar on its own; let it ride until it's part of something worth a ping.

Each `VERIFY.md` entry states:
- **What changed** — files/sections touched, one line.
- **What to check** — the exact section/page and the exact string or number expected (or not expected) to appear in the rendered PDF.
- **Why** — one line on what a failed check would mean.

The user builds (`pdflatex`→`bibtex`→`pdflatex`×2 from `GRBResearchPaper/`, per `HANDOFF.md` §1) and checks manually, then reports back or marks the entry resolved. Treat an unconfirmed `VERIFY.md` entry as an open question, not a passed check — never report a content change as done on the strength of a clean compile alone.

Exception: Claude may still run a bare `pdflatex` (not `pdftotext`) to catch a compile *error* it just introduced while actively editing `.tex` — that's a short log, not content verification, and fixing a broken build is a different task from confirming a claimed number rendered correctly. Once the build is clean, stop there; the content check goes in `VERIFY.md`.

## Verify, don't assume

Every non-trivial claim about the current state of code, data, or the compiled paper must be checked programmatically or by direct inspection of the actual output — never asserted because an operation "should have worked" or "looks right."

- **A clean LaTeX build proves nothing about content.** A whole subsection once sat commented out for most of a session while `latexmk` reported zero errors, and `BUGS.md` BUG-17 (a generated table silently missing uncertainty columns) was caught only by reading the rendered PDF, never by inspecting the generator script. So after any change that's supposed to change what the paper says, don't stop at the exit code or page count — but per "PDF build verification" above, don't run `pdftotext`/grep yourself either. Log the specific check needed in `VERIFY.md` and treat it as unverified until the user reports back.
- **A refactor that "should be behavior-preserving" gets a numeric check, not a read-through.** When a formula or helper is shared between two call sites (e.g. `compute_tau_hat` between Limit A and Limit B in `lorentz_factor.py`), verify the two outputs actually agree — bit-for-bit where the math says they must — rather than trusting that matching source implies matching output.
- **An import or dependency change gets verified against the live package, not assumed from the diff.** After changing how a module is imported, confirm every symbol it uses actually resolves under the new path (e.g. `hasattr` checks against the imported package) rather than assuming a mechanical find-and-replace was safe.
- **A "should be closed" `BUGS.md` entry gets re-run against current line numbers, not re-stated from memory.** OBS-02 sat marked OPEN for an entire session after the fix that actually resolved it (BUG-11) had already landed, because nobody re-checked the current source before re-asserting the old diagnosis.
- **Side effects of a verification step are still side effects.** Checking that code runs by executing it can itself change repo state (stray output files, regenerated plots). Inspect what actually happened (`git status`, file mtimes) afterward rather than assuming a "just testing" run was inert.

This is not a call for exhaustive testing on every change — it is a call to replace "this looks right" with one concrete, reproducible check before reporting something as done.

## Determinism checks on slow MC scripts

Reproducibility (same seed → byte-identical output across two runs) does not depend on the MC sample count — it holds at `n_samples=10` exactly as well as at the script's real `5_000`/`10_000`. When the *only* goal of a run is to confirm determinism (e.g. after a seeding refactor, verifying `rng` is threaded correctly rather than reseeded per iteration), temporarily drop the script's sample-count constant to something small (e.g. `10`) before running it twice, instead of paying the full runtime twice. This is purely a speed optimization for the verification step itself — it changes nothing about what determinism means or what counts as passing.

**If the goal is real work** — generating a number or CSV/plot that will actually be used (in the paper, in a table, in any output the user or a downstream script will read) — always run at the script's real, configured sample count. Never leave a reduced sample count in a script's on-disk default; revert the temporary edit (or restore via `git diff`/`git checkout` on just that line) before finishing, so nothing downstream silently inherits an under-sampled default.

## Using subagents

Spawning Sonnet model subagents (parallel `Explore`/`general-purpose` agents, `fork`s, etc.) is fine when it genuinely helps — e.g. independent per-GRB computation, or splitting research across the two repos. But don't spawn by default: most of this project's work (the phased plan, paper writing) is a tight dependency chain where one agent working serially is simpler and more reliable than coordinating several. When agents are spawned, keep a strict eye on them — verify their actual output/diffs against what they claimed to do before relying on it or reporting it as done (per the "trust but verify" norm for agent results), and always scope their prompts explicitly to `GRBResearchPaper`/`GRBResearchWork` (see the hard boundary below) since a fresh agent won't infer that on its own.

## Hard boundary

**Every file read, edit, write, or shell command in this project must resolve to a path inside `GRBResearchPaper/` or `GRBResearchWork/`.** Nothing else under `~/Pictures/grb_research/` (or elsewhere on the machine) is in scope unless the user explicitly names that other path in their message — and even then, treat it as a one-off, not a standing permission.

Concretely:
- Never run `find`, `grep -r`, `ls`, or launch a search agent rooted at `~/Pictures/grb_research` or `~` — always root them at `GRBResearchPaper` or `GRBResearchWork` explicitly (e.g. `find /home/syedalimohsinbukhari/Pictures/grb_research/GRBResearchWork/... `, not a bare `find ~`).
- Never `cd` out of these two directories. If a command's target isn't already an absolute path inside one of them, check `pwd` first rather than assuming.
- No destructive or global operations (`rm -rf`, `git clean`, `git reset --hard`, package installs/uninstalls, config edits) outside these two repos, full stop — and inside them, follow the normal git-safety rules (confirm before force-push/reset/etc.).
- `GRBResearchWork` has its own `.venv` — use it (not system Python, not another project's venv) for running analysis scripts, since `pyarrow` and other deps live there. Don't install into or modify any other project's environment.
- If a subagent (Explore, general-purpose, etc.) is launched, its prompt must explicitly scope it to these two directories — don't assume it will infer the boundary on its own.
- These two repos are independent git repositories with independent history — don't mix commits/branches across them, and don't assume a `git` command in one affects or should reference the other.
- **No committing to either repo.** Work in the working tree only (edit, run, analyze) — never `git commit`, `git push`, or otherwise write to either repo's history unless the user explicitly asks in a given message.

If a task seems to require touching something outside this boundary, stop and ask rather than proceeding.

## Multiple machines — project root is machine-dependent (added 2026-09-16)

**This project is worked on from more than one machine, not just the one described above.** The project root in the "Scope" section (`~/Pictures/grb_research/`, and the `/home/syedalimohsinbukhari/...` paths in "Hard boundary") describes one such machine and is kept as-is — do not edit those paths to "fix" them, and do not treat them as universally current.

**On this machine** (hostname/user `iqra-siddique`), both repos and the four root working documents (`PLAN.md`, `BUGS.md`, `HANDOFF.md`, `VERIFY.md`) live directly under `/home/iqra-siddique/PycharmProjects/grb_reseach/` — i.e. `/home/iqra-siddique/PycharmProjects/grb_reseach/GRBResearchPaper/` and `/home/iqra-siddique/PycharmProjects/grb_reseach/GRBResearchWork/`, not under any `Pictures/grb_research/` path, and not under any `syedalimohsinbukhari` home directory.

**The rule going forward: the hard boundary always applies, but re-rooted to whichever machine is actively calling for work.** Before treating either root as current, check which machine/environment this session is actually running in (e.g. the working-directory path given in the environment info, or `pwd`/`whoami`) rather than assuming from a stale memory of a prior session. Once that's established:

- Every file read, edit, write, or shell command must resolve to a path inside that machine's `GRBResearchPaper/` or `GRBResearchWork/` — under `~/Pictures/grb_research/` on the `syedalimohsinbukhari` machine, or under `/home/iqra-siddique/PycharmProjects/grb_reseach/` on this one.
- Never run `find`, `grep -r`, `ls`, or launch a search agent rooted at that machine's project-root directory itself (one level up from the two repos) or at `~`/the home directory — always root them at the `GRBResearchPaper` or `GRBResearchWork` subdirectory explicitly, using whichever machine's path is actually live this session.
- The four working `.md` files (`PLAN.md`, `BUGS.md`, `HANDOFF.md`, `VERIFY.md`) are mirrored at `GRBResearchWork/md_file_refs/` on the `syedalimohsinbukhari` machine — check there rather than at that machine's project root, which is being retired for these four (per the "Scope" section above). This mirror path is repo-relative, so it applies the same way regardless of which machine's `GRBResearchWork` checkout is live. Don't assume the same migration has happened on the other machine without checking — its copies may still sit at its own project root per the "On this machine" note above.
- If a subagent is launched, scope its prompt explicitly to the two repo paths for the machine this session is actually running on — don't assume it will infer the layout, and don't hand it the other machine's paths by default.
- Never assume the other machine's copy of either repo is reachable or in sync with the one on this machine — they are separate checkouts unless explicitly stated otherwise.