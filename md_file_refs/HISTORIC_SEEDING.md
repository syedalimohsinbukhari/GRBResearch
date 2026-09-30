# Historic archive — RNG/seeding overhaul

Consolidated 2026-10-01 from two retired standalone files, preserved **verbatim** under
their own sections below: the plan, then its execution/verification log. This is
archival/historic reference material documenting a completed piece of infrastructure work —
for the current seeding scheme and conventions see `src/grb_research/SEEDING.md` (the live,
tracked repo file this work produced) and `CLAUDE.md`'s "Determinism checks on slow MC
scripts" section in this same directory. See `HISTORIC_INDEX.md` for a one-line manifest of
where every retired file's content now lives.

---

## SEED_PLAN.md

# RNG/seeding overhaul — GRBResearchWork

## Context

A paper reviewer flagged that Monte-Carlo draws across the codebase reuse the same literal seed (mostly `12345`, sometimes `1234`/`42`) in many independent scripts — correlated random streams across supposedly-independent analyses. Investigation found the problem is worse than that: several scripts reseed an identical `SEED` constant on *every iteration* of a loop over models/GRBs/episodes (`lorentz_factor.py:308`, `lorentz_factor_limit_b.py:128`, `bb_flux_fraction.py`, `pe_er_photosphere.py`, `gbm_only_refit.py`, `model_parameters/utils.py`+`epeak_vs_kt.py`/`peak_energy.py`, `amati_relationship.py`/`amati_helpers.py`) — a fresh identically-seeded generator per call means those iterations draw *literally identical* underlying random numbers, a correctness bug, not just a style concern.

Fix: a deterministic per-script seed derived from `(script filename, MASTER_SEED)` via SHA-256, replacing every hardcoded literal. Each script builds **one** `rng` object and threads it (`rng=rng`) through every downstream MC call instead of re-passing an integer `seed=` at each call site — this is what actually fixes the reseed-per-iteration bug, since a single `Generator`'s state advances across calls instead of resetting.

`MASTER_SEED` was drawn deliberately rather than hand-picked (this goes into the paper's methods section), via `secrets.randbits(32)` (OS CSPRNG, not `np.random`):
```
$ python3 -c "import secrets, datetime; print(secrets.randbits(32)); print(datetime.datetime.now(datetime.timezone.utc).isoformat())"
2828702241
2026-09-02T15:40:46.344318+00:00
```
**`MASTER_SEED = 2828702241`**, generated 2026-09-02T15:40:46 UTC. This provenance gets written verbatim into a new `src/grb_research/SEEDING.md` (a normal tracked repo file, not a root-level untracked working doc like `PLAN.md`/`HANDOFF.md`) so the value is independently auditable. A paper-side mention (methods section, fixed seed + derivation scheme) is a separate `GRBResearchPaper` edit, out of scope here — flag as a follow-up once this lands.

Execution: **two Sonnet 5 agents, run sequentially** (Phase B imports symbols Phase A creates, so it must run after Phase A lands and is spot-checked) — Phase A scoped strictly to `src/grb_research/`, Phase B scoped strictly to `codes-for-paper/`. Both must follow the parent `CLAUDE.md` hard boundary (nothing outside `GRBResearchWork`/`GRBResearchPaper`), 120-char line length (`ruff.toml`), minimal comments (only non-obvious *why*), use `GRBResearchWork/.venv` for execution, and **must not `git commit`** without an explicit later ask.

All line numbers below were verified directly against current source (both by two independent Explore passes and by direct `Read` in this session) — they are accurate as of now, not from stale memory.

**Status: plan finalized, not yet dispatched — housekeeping pending before implementation. See HANDOFF.md for current session state.** *(Historic note: this status line is frozen at the moment the plan was written; the actual execution outcome is recorded in full in §SEED_PLAN-implementation.md below — both phases completed and were verified 2026-09-02/03.)*

---

## Phase A — core library (`src/grb_research/` only)

### A1. `grb_constants.py` — add `MASTER_SEED`
Insert after `kev_to_erg = 1.6021766208e-09` (line 92), before `GRB_COLORS` (line 94):
```python
MASTER_SEED = 2828702241  # drawn via secrets.randbits(32), 2026-09-02T15:40:46 UTC — see SEEDING.md
```
Add `"MASTER_SEED",` to `__all__` (line 8 block), after `"kev_to_erg",` (line 13).

### A2. `grb_calculations.py` — add `seed_from_name`
Add `import hashlib` to the stdlib block (before `import os`, line 3). Add `MASTER_SEED` to the constants import (line 15): `from .grb_constants import kev_to_erg, MASTER_SEED, model_n_pars`.

Insert immediately after `get_rng` (ends line 46), before `legacy_build_mp` (line 49), numpydoc-style matching `get_rng`'s docstring immediately above:
```python
def seed_from_name(name: str, master_seed: int = MASTER_SEED) -> int:
    """
    Derive a deterministic 32-bit seed from a file name and a master seed.

    Hashes the basename of `name` together with `master_seed` so different scripts get
    different, reproducible seeds without any script hardcoding a literal seed value.
    Uses `os.path.basename` rather than the full path so the derived seed is stable
    across machines and checkout locations.

    Parameters
    ----------
    name : str
        Path or file name to derive the seed from. Typically a script's ``__file__``.
    master_seed : int, optional
        Project-wide master seed mixed into the hash (default: `MASTER_SEED`).

    Returns
    -------
    int
        A deterministic seed in [0, 2**32).
    """
    digest = hashlib.sha256(f"{os.path.basename(name)}-{master_seed}".encode()).hexdigest()
    return int(digest, 16) % (2**32)
```
Do **not** modify `get_rng`'s own signature/body (stays `(seed=None, rng=None)`).

### A3. Route 4 bypass sites through `get_rng`
Confirmed exhaustive (`grep -rn "np\.random\.default_rng\|np\.random\.seed" src/grb_research/*.py`) — exactly these 4 plus `get_rng` itself.

**(a) `ModelResampler.__init__`** (`grb_calculations.py:165-186`) — keep param order `(self, model, samples, rng=None, seed=None, destroy=True)` as-is (no positional callers exist). This is a **DRY refactor, not a bug fix** — current code already gives `rng` correct precedence over `seed`. Replace lines 176-182:
```python
        self._samples = samples if destroy else samples.copy()
        self.rng: np.random.Generator = get_rng(seed=seed, rng=rng)
```
Confirmed callers, all keyword-based, no behavior change: `grb_calculations.py:149`, `grb_calculations.py:923`, `codes-for-paper/amati_relationship/amati_helpers.py:78`.

**(b) `ParameterSet.get_populated_values`** (`grb_atomic.py:36-79`) — file has no `get_rng` import. Add `from .grb_calculations import get_rng`; **first verify no circular import** by running `.venv/bin/python -c "import grb_research"` — `grb_calculations.py` doesn't import `grb_atomic.py` directly but imports `Model` from `grb_model.py`, which imports `grb_atomic.py`, so a cycle is plausible. If the plain module-level import fails, fall back to a local import as the first line of `get_populated_values` (no existing precedent for cross-module runtime cycles in this codebase beyond `TYPE_CHECKING`-gated type imports, so a plain local import is the right minimal fix). Replace lines 69-73:
```python
        # Get or create RNG
        rng_instance = get_rng(seed=seed, rng=rng)
```
Confirmed safe: both current callers (`model_parameters/utils.py:154`, `amati_helpers.py:77`) always pass non-`None` `rng`.

**(c) `plot_covariance_corner`** (`grb_utils.py:117-158`) — **real circular-import risk confirmed**: `grb_calculations.py:22` does `from .grb_utils import save_fig` at module level, so `grb_utils.py` importing `get_rng` from `grb_calculations.py` at module level would deadlock. Use a **local import inside the function body**, placed immediately above its use:
```python
    # Get or create RNG
    from .grb_calculations import get_rng  # local import: grb_calculations imports save_fig from this module

    rng_instance = get_rng(seed=seed, rng=rng)
```
No current callers exist anywhere in the repo — zero behavioral risk, this is a correctness fix for future use.

**(d) `FluxFluenceCalculator.__init__`** (`grb_calculations.py:701-719`) — replace lines 715-719:
```python
        self.rng = get_rng(seed=seed, rng=rng)
```
This also fixes the `if seed:` falsy-check bug (line 717: `seed=0` currently falls through silently). `get_rng` already raises on both-`None`, so the manual raise is now redundant.

### A4. Rename `mc_e_iso_sampler`'s `seed_number` → `seed`
`grb_calculations.py:335-390`. Rename the parameter (line 348: `seed_number=1234` → `seed=1234`), its docstring entry (line 380-381), and the internal call (line 390: `get_rng(seed=seed_number, rng=rng)` → `get_rng(seed=seed, rng=rng)`). Confirmed exhaustive search: the only external caller is `codes-for-paper/amati_relationship/amati_helpers.py:95` — that keyword rename is a **Phase B** edit. No other file calls `mc_e_iso_sampler` directly.

Note (informational, no action needed): because `get_rng` gives `rng` precedence, `amati_helpers.py:282,357`'s `seed_number=seed_number + index2`/`+ ep_idx * ...` offset arithmetic is *already* inert wherever `rng` is also passed non-`None` on those call chains — this was true before this overhaul too. Phase B may delete the now-pointless offset expressions for clarity but this isn't required for correctness.

### A5. Re-export from `__init__.py`
Add `seed_from_name` to the `.grb_calculations` import block (lines 8-15, after `mc_e_iso_sampler`). Add `MASTER_SEED` to the `.grb_constants` import block (lines 16-39).

### A6. New file — `src/grb_research/SEEDING.md`
Tracked repo file (not a root-level untracked working doc). Must cover, following this project's method-note conventions (defining source, judgement calls, provenance, validation):
- The reviewer's original concern and the reseed-per-iteration bug found during investigation (cite `lorentz_factor.py:308` as the concrete example).
- `MASTER_SEED = 2828702241` and its exact generation provenance (the `secrets.randbits(32)` command, raw output, UTC timestamp — reproduced verbatim from this plan's Context section above).
- The `seed_from_name(name, master_seed)` derivation scheme (SHA-256 of basename + master seed, mod 2**32) and why `os.path.basename` is used instead of the full path (portability across machines/checkouts).
- The per-script usage pattern (`SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`, then thread `rng=rng` everywhere) and why threading one `rng` — not re-passing `seed=` — is what fixes the reseed-per-iteration bug.
- A pointer/list of which `codes-for-paper/` scripts consume this (filled in after Phase B; Phase A can leave this list as "see Phase B scripts" or enumerate the known targets from this plan).
- Note that a paper methods-section mention is a separate, not-yet-done step in `GRBResearchPaper`.

### Phase A verification (run from `GRBResearchWork/`, using `.venv`)
```
.venv/bin/python -c "from grb_research import get_rng, seed_from_name, MASTER_SEED; print(MASTER_SEED); print(seed_from_name('foo.py')); print(seed_from_name('foo.py') == seed_from_name('/some/other/path/foo.py'))"
.venv/bin/python -c "import grb_research"
```
Expect: `MASTER_SEED == 2828702241`; `seed_from_name` returns an int in `[0, 2**32)`; path-independence check prints `True`; plain import succeeds (if it fails inside `grb_atomic.py`, that confirms the predicted cycle and the local-import fallback from A3(b) must be applied). No pytest suite exists in this repo (confirmed, no `tests/` dir) — these smoke checks are the verification.

Phase A must not touch `codes-for-paper/`.

---

## Phase B — front-end scripts (`codes-for-paper/` only, runs after Phase A lands and is spot-checked)

General pattern per file: remove the hardcoded seed constant; add `SEED = seed_from_name(__file__)` near the top; build one `rng = get_rng(seed=SEED)`; thread `rng=rng` explicitly through every downstream MC call listed below (explicit at each call site, not relying on a function default silently capturing a module global). CSV writers keep `"seed": SEED` (now a large derived int) and `"n_samples"` unchanged in shape.

**B1. `fluence/grb_fluence.py`** — already structurally correct (spawns 2 independent child streams via `rng.spawn(2)`, threads through the loop) — use as the reference pattern. Only fix the seed source: `RANDOM_SEED = 12345` → `RANDOM_SEED = seed_from_name(__file__)`; `np.random.default_rng(RANDOM_SEED)` (line 84) → `get_rng(seed=RANDOM_SEED)`. Add `get_rng, seed_from_name` to the `grb_research` import.

**B2. `lorentz_factor/lorentz_factor.py`** — `SEED = 12345` (line 127) → `SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`. Line 308 (inside the per-GRB/per-model loop, lines 285/289): `draw_model_samples(model, n_samples=N_SAMPLES, seed=SEED)` → `..., rng=rng)`. CSV (346) unchanged in shape.

**B3. `lorentz_factor/lorentz_factor_limit_b.py`** — derive its **own independent** seed (don't import `SEED` from `lorentz_factor.py` as it does today — `__file__` differs naturally, decorrelating Limit A vs Limit B). Remove `SEED` from its `from lorentz_factor import (...)` list; add `get_rng, seed_from_name` to its `grb_research` import; add `SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`. Line 128 (loop at 106/110): `seed=SEED` → `rng=rng`. CSV (162-163) unchanged in shape. *(Verified independently: the bit-identical `tau_hat` cross-check between this file and `lorentz_factor.py`, documented in `lorentz_factor.md` §8.2, depends on both scripts sharing the same deterministic closed-form `compute_tau_hat` on best-fit point values — not on a shared seed — so giving this file its own seed does NOT break that identity. Only `Gamma_min_err_lower/upper` depend on the seed.)*

**B4. `gbm_only_refit/gbm_only_refit.py`** — `SEED = 12345` (line 96) → derived `SEED`/`rng`. `compute_f_bb(model, n_samples=N_SAMPLES, seed=SEED)` → `(model, n_samples=N_SAMPLES, rng=None)`, forward `rng=rng` internally. `paired_ratio_pct(joint_model, gbm_model, n_samples=N_SAMPLES, seed=SEED)` → `(..., rng=None)`; **delete the `seed`/`seed+1` split** (lines 189-190) — a single shared `rng` threaded sequentially through both `draw_model_samples` calls already gives two independent draws, which is simpler and is what actually fixes cross-episode reuse (the old `+1` only decorrelated joint-vs-gbm *within* one episode, not across the `for ep in EPISODE_ORDER:` loop at line 615). Update the docstring paragraph explaining the old seed/seed+1 mechanism, but **keep** its empirical justification (the `r=0.9996` correlation-collapse finding) — that finding is still the reason a shared, non-reseeded generator matters, only the mechanism changed. Thread `rng` explicitly through `bb_fields` → `build_row(ep, joint_interval, gbm_interval, rng)` → the `main()` loop call site (line 617).

**B5. `bb_fraction/bb_flux_fraction.py`** — `SEED = 12345` (line 119) → derived. `compute_fraction(..., seed=SEED, rng=None)` → default `seed=None` (drop the module-constant default); call site (line 307, in loop at 302/304) passes `rng=rng` explicitly.

**B6. `photospheric_radius/pe_er_photosphere.py`** — `SEED = 12345` (line 123) → derived. `interval_samples(model, n_samples=N_SAMPLES, seed=SEED)` → add `rng=None` param, drop `seed=`; call site (line 276, loop at 263/272) passes `rng=rng`.

**B7. `model_parameters/utils.py`** — add `get_rng` to its `grb_research` import. `convert_sbpl_to_band` (line 137): replace its internal `if seed is not None: rng = np.random.default_rng(seed)` (lines 149-150, which currently **overwrites a caller-supplied `rng` whenever `seed` is also given** — a real bug) with `rng = get_rng(seed=seed, rng=rng)` (confirmed safe: no current caller passes both). `extract_kt_epeak_from_models(models, t90_marker="o", seed=1234)` (line 267) → `(models, t90_marker="o", seed=None, rng=None)`, build `rng = get_rng(seed=seed, rng=rng)` once inside, thread into its internal per-model loop (line 302) as `convert_sbpl_to_band(model_, rng=rng)` instead of `seed=seed`. Note: dropping the `seed=1234` default means a caller passing neither now raises — confirmed both callers (B8, B9 below) are updated in this same phase to always pass `rng=`.

**B8. `model_parameters/epeak_vs_kt.py`** — `SEED = 1234` (line 34) → derived `SEED`/`rng`. The 4 calls to `extract_kt_epeak_from_models(models_XXXX, seed=SEED)` (lines 85,88,91,94, one per GRB) → `rng=rng` (sequential threading of the one `rng` across all 4 calls — simpler than spawning, correct since call order is fixed). CSV (180-181) unchanged in shape.

**B9. `model_parameters/peak_energy.py`** — `SEED = 42` (line 28) → derived `SEED`/`rng`. `extract_peak_energy(model_collection, seed=SEED, n_samples=N_SAMPLES)` → `(model_collection, rng=None, n_samples=N_SAMPLES)`, its internal loop's `convert_sbpl_to_band(model, n_sample=n_samples, seed=seed)` (line 49) → `rng=rng`. Call site (line 72): `extract_peak_energy(best)` → `extract_peak_energy(best, rng=rng)`. CSV (154-155) unchanged in shape.

**B10. `amati_relationship/amati_relationship.py` + `amati_relationship/amati_helpers.py`**:
- `amati_helpers.py`: add `Optional` to its typing import; add `get_rng` to its `grb_research` import. `mc_e_iso_sampler(..., seed_number=seed_number, rng=rng)` (line 95) → `seed=seed_number` (keyword rename only, matches A4; the local var name `seed_number` stays). `plot_grbs_over_amati_relationship` and `plot_unknown_redshift_grb` (lines 214, 298) each add an `rng: Optional[np.random.Generator] = None` param and replace their internal `rng = np.random.default_rng(seed_number)` (lines 265, 332) with `rng = get_rng(seed=seed_number, rng=rng)`. Document `rng` in each docstring.
- `amati_relationship.py`: `n_seed = 12345` (line 43) → `SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`. Add `get_rng, seed_from_name` to its `grb_research` import. The `plot_grbs_over_amati_relationship(...)` call (lines 66-74) drops `seed_number=n_seed`, adds `rng=rng`. The `plot_unknown_redshift_grb(...)` call inside the per-episode loop (lines 91-101) — **this is the actual reseed-per-iteration bug fix**: drops `seed_number=n_seed`, adds `rng=rng`, so the loop no longer rebuilds an identically-seeded fresh generator on every call.
- **Add missing CSV provenance** (confirmed: the `q` DataFrame written to `amati_relationship.csv`, lines 135-159, currently has no `n_samples`/`seed` columns at all — a pre-existing gap). Add `n_samples_col = np.full(len(g_name), n_sample)` and `seed_col = np.full(len(g_name), SEED)` to the DataFrame's row list and `"n_samples"`, `"seed"` to `q.columns`, matching the naming convention used by the other 8 CSV writers in this pass. Check whether `amati_relationship/csv_to_latex.py` reads columns positionally (would break) vs. by name (safe) — adjust if positional.

**B11. `butterfly_plots/butterfly_all.py`** — line 28: `np.random.default_rng(seed=42)` → `get_rng(seed=seed_from_name(__file__))`. Add `get_rng, seed_from_name` to its `grb_research` import. Confirmed: writes no CSV, no provenance-column work needed.

**B12. `variability_analysis/norris..py`** (low priority — dormant, no callers found anywhere in `codes-for-paper/`) — `SEED = 12345` (line 9) → `SEED = seed_from_name(__file__)`; import `get_rng, seed_from_name` from `grb_research` directly (drop the current `from grb_research.grb_calculations import get_rng`). `tv_mc_summary`'s own `rng = get_rng(seed=seed)` stays as-is (correct in isolation). Add a one-line comment near its signature flagging that a future caller looping over `pulse_index` should thread one `rng` rather than reusing `seed=SEED` per iteration, to avoid reintroducing this exact bug class. *(Historic note: this file — `norris..py`, a stray pre-rename duplicate of the later `norris_fit.py` — is itself dead/no-caller as of the 2026-09-27 variability-timescale work; see `BUGS.md`'s norris_fit history and `seed_registry.yaml`'s dormant entry for its current status.)*

**B13. Method-note updates** (per user decision — do this in the same pass): update the 5 existing per-folder notes that document the old literal seed and/or a since-removed judgement call — `lorentz_factor/lorentz_factor.md` (also fix the stale "(same seed, SEED=12345)" parenthetical in its `tau_hat` bit-identical claim — see B3's note; the identity itself still holds, only the wording is stale), `gbm_only_refit/gbm_only_refit.md` (the seed/seed+1 judgement call section — update to describe sequential-rng-threading instead), `bb_fraction/bb_fraction.md`, `fluence/fluence.md`, `photospheric_radius/photospheric_radius.md`. One-line/short-paragraph edits each — no new `.md` files for folders that don't already have one (`amati_relationship/`, `model_parameters/`, `butterfly_plots/`, `variability_analysis/` — pre-existing gaps, out of scope here).

### Phase B verification (per script, run via `.venv/bin/python` from the script's own directory)
1. Run each ported script once — exits 0, no new tracebacks.
2. Check the CSV's `seed` column: now a large derived int (billions-range), constant-valued across all rows of one run (`df["seed"].nunique() == 1`).
3. **Run each script a second time and diff its CSV output against the first run** — must be byte-identical/value-identical (determinism check: `seed_from_name` is a pure function of basename, so reruns reproduce exactly). A script whose two runs differ indicates leftover unseeded randomness or an `rng` being unintentionally shared/mutated — treat as a real bug, not something to wave through.
4. **Do not** compare the new CSV against the pre-port committed CSV and expect a match — every MC-derived column (percentiles, error bars) legitimately shifts because the seed itself changed. State this explicitly per script when reporting results, so it isn't mistaken for a regression.
5. After all scripts run, `git status`/`git diff --stat` across `codes-for-paper/` and report every regenerated file (CSVs, PNGs, PDFs) as a side effect rather than leaving it for the user to discover.
6. Cross-check the `lorentz_factor.md` `tau_hat`-bit-identical claim still holds after B3 (it should, per the note above) — this is a real numeric check, not just a doc edit.

Phase B must not touch `src/grb_research/`.

---

## Execution

1. Launch Phase A as one Sonnet 5 agent, scoped strictly to `src/grb_research/`, given A1–A6 + its verification section verbatim.
2. Review Phase A's actual diff and run its verification commands myself (or confirm the agent ran them and report the output) before proceeding.
3. Launch Phase B as one Sonnet 5 agent, scoped strictly to `codes-for-paper/`, given B1–B13 + its verification section verbatim, only after Phase A is confirmed working.
4. Review Phase B's diff, spot-check a couple of the determinism reruns myself, and report the full `git status` picture (including regenerated data files) back to the user. No commits at any point.

---

## SEED_PLAN-implementation.md

# RNG/seeding overhaul — implementation log

Companion to `SEED_PLAN.md` (§SEED_PLAN.md in this archive). That file is the plan; this file tracks what has
actually landed, item by item, with a short description of what changed. Updated as each phase's
diff is reviewed and verified — not written from an agent's self-report alone.

Status legend: ⬜ not started · 🟨 dispatched / in review · ✅ done and verified

---

## Phase A — core library (`GRBResearchWork/src/grb_research/`)

✅ **Done and verified** 2026-09-02. Implemented by a background agent scoped strictly to `src/grb_research/`; diff reviewed by hand line-by-line against the plan, and both verification commands re-run independently (not just taken from the agent's self-report). No `codes-for-paper/` files touched, nothing committed.

| Item  | Status | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|-------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| A1    | ✅     | Added `MASTER_SEED = 2828702241` to `grb_constants.py`, with provenance comment; added to `__all__`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| A2    | ✅     | Added `seed_from_name(name, master_seed)` to `grb_calculations.py` — SHA-256(basename + master seed) mod 2**32, so each script gets a distinct deterministic seed without a hardcoded literal.                                                                                                                                                                                                                                                                                                                                                                                                          |
| A3(a) | ✅     | `ModelResampler.__init__` — DRY refactor, routed through `get_rng(seed=seed, rng=rng)` instead of duplicating its precedence logic inline. **Superseded 2026-09-03**: `seed` param removed entirely, `rng` now required keyword-only (see "Post-Phase-B cleanup" below) — every real caller already passed `rng=`.                                                                                                                                                                                                                                                                                      |
| A3(b) | ✅     | `ParameterSet.get_populated_values` (`grb_atomic.py`) — same refactor. The predicted circular import (`grb_atomic → grb_calculations → grb_model → grb_atomic`) was confirmed real by testing a module-level import first (it failed with a partially-initialized-module `ImportError`); fixed with a local import inside the method, with a comment explaining why. **Superseded 2026-09-03**: `seed` param and the local `get_rng` import both removed — `rng` is now required keyword-only, and the circular-import workaround is gone entirely since the function no longer needs `get_rng` at all. |
| A3(c) | ✅     | `plot_covariance_corner` (`grb_utils.py`) — same refactor via a **local** import (real circular-import risk: `grb_calculations.py` imports `save_fig` from `grb_utils.py` at module level). No current callers, so zero behavioral risk — a correctness fix for future use. **Superseded 2026-09-03**: `seed` param and the local import both removed the same way as A3(b); `rng` required keyword-only.                                                                                                                                                                                               |
| A3(d) | ✅     | `FluxFluenceCalculator.__init__` — same refactor. Also fixed a real bug: the old `if seed:` check silently ignored `seed=0`. **Superseded 2026-09-03**: `seed` param removed, `rng` required keyword-only — the `seed=0` bug class is now structurally impossible here since there's no `seed` path left to have the bug.                                                                                                                                                                                                                                                                               |
| A4    | ✅     | Renamed `mc_e_iso_sampler`'s `seed_number` parameter → `seed` (signature, docstring, internal call), for naming consistency with the rest of the API. Its one external caller, in `codes-for-paper/amati_relationship/amati_helpers.py`, is untouched — that keyword rename is Phase B. **Superseded 2026-09-03**: `seed` itself removed in the post-Phase-B cleanup — `rng` is now required keyword-only, and the caller's keyword dropped to match.                                                                                                                                                   |
| A5    | ✅     | Re-exported `seed_from_name` and `MASTER_SEED` from `src/grb_research/__init__.py`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| A6    | ✅     | New file `src/grb_research/SEEDING.md` — records the reviewer's original concern, the reseed-per-iteration bug this overhaul fixes (citing `lorentz_factor.py:~308`), `MASTER_SEED`'s exact generation provenance, the `seed_from_name` derivation scheme and why `os.path.basename` is used, the per-script usage pattern (with a code example), a stub pointing to Phase B for the consumer list, and a note that the paper methods-section mention is a separate not-yet-done step in `GRBResearchPaper`.                                                                                            |

**Verification — re-run independently, not just taken from the agent's report:**
```
$ PYTHONPATH=src .venv/bin/python -c "from grb_research import get_rng, seed_from_name, MASTER_SEED; ..."
2828702241
2768334463
True
$ PYTHONPATH=src .venv/bin/python -c "import grb_research"
import OK
```
`py_compile` clean on all 5 touched files. Two lines over 120 chars found by a line-length scan: `grb_atomic.py:98` (pre-existing, not part of this diff — confirmed via `git diff`) and `SEEDING.md:25` (a verbatim shell command inside a markdown code block, not Python source — not a `ruff.toml` violation).

---

## Phase B — front-end scripts (`GRBResearchWork/codes-for-paper/`)

✅ **Done and verified** 2026-09-02/03. Dispatched to a background agent scoped to `codes-for-paper/`; the agent's own verification pass stalled twice (it tried to background long MC scripts and wait across turns, which a subagent can't do — the session that dispatched it also restarted mid-way, orphaning a first round of my own background verification jobs). All code edits landed correctly regardless; every script's determinism check below was run — and re-run after the restart — directly by the main session, not taken from any agent's self-report.

| Item | Status         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|------|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| B1   | ✅             | `fluence/grb_fluence.py` — literal `RANDOM_SEED = 12345` → `seed_from_name(__file__)`. Already threaded `rng` correctly otherwise (reference pattern for the rest of Phase B).                                                                                                                                                                                                                                                                                   |
| B2   | ✅             | `lorentz_factor/lorentz_factor.py` — derived `SEED`/`rng`; fixes the reseed-per-iteration bug in the per-GRB/per-model loop (`draw_model_samples(..., seed=SEED)` → `rng=rng`).                                                                                                                                                                                                                                                                                  |
| B3   | ✅             | `lorentz_factor/lorentz_factor_limit_b.py` — gets its **own** independent derived seed (decorrelated from Limit A) rather than importing Limit A's `SEED`.                                                                                                                                                                                                                                                                                                       |
| B4   | ✅             | `gbm_only_refit/gbm_only_refit.py` — derived `SEED`/`rng`; deleted the old `seed`/`seed+1` split in `paired_ratio_pct` in favor of sequential threading of one `rng` through both `draw_model_samples` calls (keeps the documented `r=0.9996` correlation-collapse finding as justification in the docstring — only the mechanism changed, from independent seeds to sequential draws from one generator).                                                       |
| B5   | ✅             | `bb_fraction/bb_flux_fraction.py` — derived seed, `rng=` threaded explicitly at the call site; `compute_fraction`'s `seed` default dropped.                                                                                                                                                                                                                                                                                                                      |
| B6   | ✅             | `photospheric_radius/pe_er_photosphere.py` — same pattern.                                                                                                                                                                                                                                                                                                                                                                                                       |
| B7   | ✅             | `model_parameters/utils.py` — fixed a real bug in `convert_sbpl_to_band`: it previously overwrote a caller-supplied `rng` whenever `seed` was also given. `extract_kt_epeak_from_models` now threads `rng` instead of a reused `seed=1234` default.                                                                                                                                                                                                              |
| B8   | ✅             | `model_parameters/epeak_vs_kt.py` — derived seed; sequential `rng` threaded across its 4 per-GRB calls.                                                                                                                                                                                                                                                                                                                                                          |
| B9   | ✅             | `model_parameters/peak_energy.py` — derived seed; `extract_peak_energy` takes `rng` instead of `seed`.                                                                                                                                                                                                                                                                                                                                                           |
| B10  | ✅             | `amati_relationship/amati_relationship.py` + `amati_helpers.py` — the actual reseed-per-iteration bug fix for this pair: the per-episode loop calling `plot_unknown_redshift_grb` no longer rebuilds an identically-seeded generator every call. Added the previously-missing `n_samples`/`seed` provenance columns to `amati_relationship.csv`; confirmed `csv_to_latex.py` reads columns by name, not position, so no adjustment needed there. |
| B11  | ✅             | `butterfly_plots/butterfly_all.py` — literal `seed=42` → derived. Writes no CSV; verified by two clean runs instead.                                                                                                                                                                                                                                                                                                                                             |
| B12  | ✅ (code only) | `variability_analysis/norris..py` — low priority, dormant/no callers (confirmed). Literal seed → derived; import changed to pull `get_rng`/`seed_from_name` from the package root; comment added flagging the same bug class for any future caller that loops over `pulse_index`. Not executed — no runnable entrypoint/CSV to verify against.                                                                                                                   |
| B13  | ✅             | Method-note updates — `lorentz_factor.md`, `gbm_only_refit.md`, `bb_fraction.md`, `fluence.md`, `photospheric_radius.md` all updated to describe the new seeding scheme. Three of the five (`bb_fraction.md`, `fluence.md`, `lorentz_factor.md`) are `.gitignore`d, which is why they never appeared in `git status` despite being correctly edited — confirmed by direct content read, not by git.                                                              |

**Verification — every script run twice from a clean state and diffed, not just taken from an agent report:**

| Script                      | Run 1/2 exit | CSV determinism                                                        | seed value |
|-----------------------------|--------------|--------------------------------------------------------------------------|------------|
| `grb_fluence.py`            | 0/0          | seed consistent across both output CSVs (shared `rng.spawn(2)` parent) | 2285298891 |
| `lorentz_factor.py`         | 0/0          | byte-identical                                                         | 3983203661 |
| `lorentz_factor_limit_b.py` | 0/0          | byte-identical                                                         | 4169226666 |
| `gbm_only_refit.py`         | 0/0          | byte-identical                                                         | 3624228654 |
| `bb_flux_fraction.py`       | 0/0          | byte-identical                                                         | 668559939  |
| `pe_er_photosphere.py`      | 0/0          | byte-identical                                                         | 3315768413 |
| `epeak_vs_kt.py`            | 0/0          | byte-identical (both CSVs)                                             | 2177701742 |
| `peak_energy.py`            | 0/0          | byte-identical                                                         | 2933599809 |
| `amati_relationship.py`     | 0/0          | byte-identical                                                         | 3271127181 |
| `butterfly_all.py`          | 0/0          | n/a (no CSV) — two clean runs                                          | 2000318690 |

Every derived seed above is distinct (as expected — each is `seed_from_name` of a different `__file__`), and every determinism check passed. **`tau_hat` bit-identical cross-check (B3's key numeric claim) re-verified after the reseed:** merged `lorentz_results.csv` and `lorentz_results_limit_b.csv` on `GRB`+`episode`, max abs diff in `tau_hat` = `0.0` across all 26 rows.

New CSVs are **not** expected to match the pre-Phase-B committed CSVs — every MC-derived column (percentiles, error bars) legitimately shifted because the seed itself changed; that is not a regression. Full `git status --short`/`git diff --stat codes-for-paper/` after all runs: 48 files changed (15 `.py`, 2 `.md`, 9 CSVs, 22 regenerated PNG/PDF/`.tex`), all consistent with the scripts that were rerun.

**Process note, not a plan defect:** while running Phase B verification, the user asked to add a project convention (now in `CLAUDE.md`, new "Determinism checks on slow MC scripts" section) — when a run's only goal is confirming reproducibility, temporarily drop the script's sample count (e.g. to `10`) rather than paying the full runtime twice, since determinism doesn't depend on sample count. This wasn't applied retroactively to the runs above (already in flight/complete by the time it was added) but should be used for any future re-verification.

---

## Post-Phase-B cleanup: stripping the vestigial `seed=` parameter

✅ **Done and verified** 2026-09-03, at the user's request, after discussing the design with them. Once Phase B landed, every real call site in `codes-for-paper/` passed `rng=`, never `seed=` — so the `seed=`/`get_rng(seed=seed, rng=rng)` fallback on 11 downstream functions was dead weight, kept only for hypothetical ad-hoc/notebook callers. Stripped it, making `rng` a required keyword-only argument instead — `seed` is now gone from every one of these signatures, and `get_rng` itself is called only once per script, at the single root `SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)` line (see "Usage pattern" in `SEEDING.md`), never inside any downstream function:

- `src/grb_research/grb_calculations.py`: `draw_model_samples`, `mc_e_iso_sampler`, `mc_spectra_sampler`, `plot_all_models`, `ModelResampler.__init__`, `FluxFluenceCalculator.__init__`.
- `src/grb_research/grb_atomic.py`: `ParameterSet.get_populated_values` — this also fully removed the local `get_rng` import (the circular-import workaround from Phase A), since the function no longer needs `get_rng` at all.
- `src/grb_research/grb_utils.py`: `plot_covariance_corner` — same, local import removed entirely.
- `codes-for-paper/model_parameters/utils.py`: `convert_sbpl_to_band`.
- `codes-for-paper/bb_fraction/bb_flux_fraction.py`: `compute_fraction` (not in the original 11, but its call to `draw_model_samples(..., seed=seed, ...)` would have broken otherwise, so fixed in the same pass).
- `codes-for-paper/amati_relationship/amati_helpers.py`: `plot_grbs_over_amati_relationship`, `plot_unknown_redshift_grb`, plus their two private helpers `_plot_model_point` and `_compute_ep_eiso` (not separately listed, but part of the same necessary cascade — `seed_number` threaded through all four and is now fully removed, including the now-dead `seed_number + index2` / `seed_number + ep_idx * len(z_values) + z_idx` offset arithmetic the plan had already flagged as inert). Dropped the now-unused `get_rng` import and the now-unused `index`/`ep_idx` loop variables that existed only to feed that arithmetic.

**Closed 2026-09-03, same day:** `codes-for-paper/model_parameters/utils.py::extract_kt_epeak_from_models` was the one function with this pattern left standing after the pass above — flagged by external review as exactly the kind of asymmetry that invites a future caller to copy the stale `seed=`/`rng=` dual-parameter style since it's "the last one still doing it that way." All 4 real call sites (`epeak_vs_kt.py`) already passed `rng=` only, so this was equally safe to strip: `seed` removed, `rng` now required keyword-only, the internal `get_rng` call removed, and the now-fully-unused `get_rng` import dropped from the file. Re-verified `epeak_vs_kt.py` (its only caller) byte-identical to the pre-change CSVs (both `epeak_vs_kt_points.csv` and `epeak_vs_kt_odr_fits.csv`). There is now exactly one place in the whole codebase where a bare `seed` becomes an `rng` — the root `get_rng(seed=SEED)` call in each script — with no competing pattern left anywhere to be mistaken for current practice.

**Verification — the refactor is provably behavior-preserving, and this was checked empirically, not just argued:** every stripped call site already always passed `rng=<non-None>`, and `get_rng`'s own logic is `if rng is not None: return rng` — so removing the wrapper call cannot change any actual random draw. Confirmed empirically by rerunning all 8 fast/moderate scripts and diffing output against the CSV already sitting in the tree from the pre-refactor, already-verified run:

| Script                      | Result                             |
|-----------------------------|-------------------------------------|
| `lorentz_factor.py`         | byte-identical to pre-refactor CSV |
| `lorentz_factor_limit_b.py` | byte-identical                     |
| `epeak_vs_kt.py`            | byte-identical (both CSVs)         |
| `peak_energy.py`            | byte-identical                     |
| `gbm_only_refit.py`         | byte-identical                     |
| `amati_relationship.py`     | byte-identical                     |
| `grb_fluence.py`            | byte-identical (both CSVs)         |
| `butterfly_all.py`          | clean exit (no CSV to diff)        |

For the two slow scripts (`bb_flux_fraction.py`, `pe_er_photosphere.py`, ~10 min/run at full `N_SAMPLES=10_000`), used the new `CLAUDE.md` "Determinism checks on slow MC scripts" convention: temporarily set `N_SAMPLES=10`, ran once as a pure smoke test (real CSV backed up first, restored after), confirmed exit 0 with no `TypeError`, then reverted `N_SAMPLES` back to `10_000` immediately. Not diffed against baseline (reduced sample count changes the output shape), but the mathematical no-op guarantee above covers correctness regardless of sample count — this was purely to catch a missed call site.

**One unrelated thing caught and fixed while checking `git status` after this pass:** the original Phase A subagent had run `git add` on `src/grb_research/SEEDING.md` despite the "working tree only, no commits" instruction — nothing was ever committed, but the file sat staged. Unstaged it (`git restore --staged`) so it's a plain untracked file like every other change in this session.

## Paper integration (`GRBResearchPaper`)

✅ **Done** 2026-09-03, after a scoping discussion with the user. Explicit decisions made before implementation:
- One short paragraph in the methods/data-analysis section (not a full spawn/child-rng writeup) — `section-5-data-analysis.tex`, right after the existing cosmology sentence.
- The per-script seed number goes in every genuinely MC-derived figure's caption, not every figure and not just a representative one.
- No conversion to "fully MC" for figures currently using a direct fit-covariance error (`peak_energy_best__all`'s BAND/CPL points, `low_index`/`high_index`) — user's call after discussing the tradeoff: MC would add seed-dependent noise to an already-exact number for directly-fitted parameters, with no accuracy gain, and is a real methodology change (not a caption fix).
- A canonical seed table in a new short appendix, generator script kept at `codes-for-paper/` root (not a topic subfolder, not `src/grb_research/`) since the table spans every topic folder and isn't analysis code.

**What was found during the figure audit (not assumed):**
- 6 figures cleanly map to one script's seed each: `peak_energy_best__all`, `amati_relationship`, `bb_flux_fraction`, `bb_flux_fraction_rest_vs_z` (same script as the previous), `pe_er_photosphere`, `butterfly_all` (section 4) — found by grepping for every `\includegraphics` in the paper, not just the ones already discussed.
- `gamma_comparison` reads three already-seeded CSVs and does no sampling of its own — doesn't get a bare number, gets a sentence pointing at the table for its three constituent seeds.
- Both `gbm_only_refit` figures actually in the paper turn out to plot non-MC quantities (checked the plotting code, not the caption wording) — the one MC-derived output this script produces (`gbm_only_refit_fractional_diff`) isn't referenced by the paper at all. Neither gets a seed citation in its own caption, but the script itself still gets a row in the appendix table (it's a real, non-dormant consumer of the scheme, just not reflected in those two particular figures).

**Built — revised once, after user review of the first draft:**
- `codes-for-paper/seed_registry.yaml` — the actual source of truth: slug → script path, categorized `active` / `unused` / `dormant` (`epeak_vs_kt.py` unused — not referenced by any figure/table in the current draft; `norris..py` dormant — no caller anywhere) so a script that stops being used is never silently dropped from view, always kept with a one-line reason. New general convention, also written into `SEEDING.md`. Only the script path is stored — never a seed integer — so the registry can't drift from what `seed_from_name` actually computes.
- `codes-for-paper/seed_table_to_latex.py` — reads the registry, computes every seed live via `seed_from_name`, writes `seed_table.csv` (all 11 rows: active + unused + dormant, for a complete audit trail) and `seed_table.tex` (active rows only, 9, ordered figures-first in paper-appearance order then the one table-only row last — not alphabetically, not by topic folder, since a reader arrives here from a caption and scans top-to-bottom the way they read the paper). Columns are `Seed | Item | Fig./Table` (seed first, at the user's request) — no script-filename column in the rendered table (dropped after the first draft; the CSV still has it).
- `GRBResearchPaper/appendices/appendix_seeding_table.tex` — thin wrapper, same pattern as `appendix_LAT_info.tex`.
- `GRBResearchPaper/tex_files/generated/seed_table.tex` — the generated table itself, copied in manually (cross-repo copy is always manual, per the project convention).
- `main.tex` gained a third appendix section, "Monte Carlo seeds."
- Seven figure captions (5 single-seed + `gamma_comparison`'s multi-seed sentence + `butterfly_all`) updated across `section-4-joint-analysis-results.tex` and `section-5-data-analysis.tex`.

**Verification:** Two real formatting bugs were caught and fixed by actually compiling, not just eyeballing the generated `.tex` — a `table` environment overflowed the single-column width (switched to `table*`, full page width), and even at `table*` width, plain columns still overflowed by more (110pt, worse than before) because they don't wrap long text — fixed with `\resizebox{\textwidth}{!}{...}` around the tabular, matching the project's existing pattern for other wide tables (`lat_info_table.tex`). A YAML-ordering mistake (the table-only row landed first instead of last) was caught by inspecting the actual rendered row order after generating, not assumed correct from the source. Final bare `pdflatex` pass (compile-error check only, per `CLAUDE.md`'s exception clause — not content verification), after the registry refactor: 0 errors, 0 undefined references, 0 overfull boxes, 28 pages. Content correctness (does the prose read well, is the table legible at its rendered size) is logged in `VERIFY.md` as **PENDING** for the user's manual check, not claimed as done.

**Side effect:** the bare compile check modified `out/main.pdf`, a tracked file — flagged to the user, who said not to worry about it.

## Asset sync tooling (`sync_paper_assets.py`)

✅ **Done** 2026-09-03, at the user's request — same registry pattern as `seed_registry.yaml`, extended to cover every figure PNG and generated table the paper actually uses, plus a single runner that copies them across the repo boundary instead of the fully-manual process `CLAUDE.md` otherwise documents.

- `codes-for-paper/figure_registry.yaml` and `codes-for-paper/table_registry.yaml` — same `active` / `unused` / `no_known_source` categorization as the seed registry (a script/output is never silently dropped from the file, always kept with a reason). 16 active figures, 2 unused (`epeak_vs_kt.png`, `gbm_only_refit_fractional_diff.png`), 0 no-known-source. 8 active tables, 3 no-known-source (the 8 per-GRB `table_*_large.tex` files, `table_template.tex`, `time-integrated-table.tex` — no generating script found anywhere in `GRBResearchWork` via search, so left un-automated rather than guessed at).
- `codes-for-paper/sync_paper_assets.py` — root-level, reads both registries, copies every `active` entry's `source` (relative to `GRBResearchWork/`) to `dest` (relative to `GRBResearchPaper/`, computed as a sibling directory), creating destination folders as needed. Copies only — never regenerates. Asserts each source exists and reports any missing one by name rather than skipping silently; exits non-zero if anything's missing.
- **First run found a real gap, not just validated the tool**: several images/tables that were regenerated during the earlier Phase B seed-refactor reruns this session (`amati_relationship.png`, `bb_flux_fraction.png`, `bb_flux_fraction_rest_vs_z.png`, `pe_er_photosphere.png`, `peak_energy_best__all.png`, `butterfly_all.png`, both `gbm_only_refit` PNGs, `fluence_table.tex`, `lorentz_table.tex`, `lorentz_table_limit_b.tex`) had never actually been copied into `GRBResearchPaper` — only `seed_table.tex` had been manually copied earlier in this session. The sync brought all of them current in one pass.
- **Unrelated git anomaly found while checking the sync's `git status`, not caused by it**: `out/main.pdf` is staged for deletion in the index (`git status` shows `D  out/main.pdf` and `?? out/`) despite the file existing on disk and being present in `HEAD` — the signature of a `git rm --cached` or similar run outside this session. Not touched; flagged to the user rather than guessed at, per the project's git-safety convention (investigate an unfamiliar state, don't act on it without knowing whose change it is).

## Performance investigation: `bb_flux_fraction.py`, and the `N_GRID` fix

Prompted by the user asking whether `bb_flux_fraction.py` (the slowest script in the pipeline) could be sped up. Two approaches were tried, in order, with the second one actually landing.

**Attempt 1 — batching + multiprocessing (tried, then partly reverted).** Profiling found ~96% of `compute_fraction`'s time went into re-running the full `N_SAMPLES x N_GRID` spectral integration once per swept redshift (25 z-values). Two changes were built and verified correct (`grb_calculations.py::component_energy_fluxes`/`_component_energy_flux_block` extended to accept a 2D `(n_z, n_grid)` energy grid, batching all redshifts into one call — verified bit-identical to looping per z on real model data, and the existing 1D-grid path used by other consumers verified byte-for-byte unchanged; `bb_flux_fraction.py::collect_results()` parallelized across the ~13 independent BB-episodes via `multiprocessing.Pool`, each worker given its own spawned child `rng` via `rng.bit_generator.seed_seq.spawn(n_jobs)` — the exact scenario `SEEDING.md`'s "Parallelism" section had already anticipated). A real memory bug was caught before it could matter: the batched path multiplies per-chunk memory by `n_z`, and at this project's real `N_SAMPLES`/`N_GRID` would have demanded ~50 GB across 8 workers on a 15 GB machine — fixed by scaling `CHUNK_SIZE` down by `n_z` when the grid is 2D, checked against `free -h` before running anything at scale.

Despite all of that being individually correct, the full-scale run projected to ~48 minutes (worse than the ~10 minute baseline) once actually timed at scale, because the batching didn't reduce real floating-point work — it only removed Python-loop overhead, and the `CHUNK_SIZE` shrink needed to keep memory safe put the block-iteration count back to roughly where it started. **User's call: revert the multiprocessing half, keep the batching half.** `collect_results()` is back to the original single-process sequential loop threading one shared `rng`; the 2D-grid batching capability stays in `grb_calculations.py` (harmless, still correct, potentially useful again later) and in `bb_flux_fraction.py`'s `compute_fraction`.

**Attempt 2 — `N_GRID` reduction (the one that actually worked).** A full sweep of every `N_GRID` consumer in the codebase found it collapses to two energy ranges (1 keV-10 MeV, used by `bb_flux_fraction.py`/`pe_er_photosphere.py`/`gbm_only_refit.py`/`mc_e_iso_sampler`'s bolometric integral; 10 keV-10 GeV, used by `mc_spectra_sampler`/`SpectralModels`/`plot_all_models`/`mc_e_iso_sampler`'s detector integral) — and that `gbm_only_refit.py` had its own un-migrated local `N_GRID = 1_000`, live in production the whole time, missed by the user's earlier N_GRID/N_SAMPLES centralization pass. A real convergence test (not assumed) — every model type in this sample (SBPL_BB, BAND_BB, BAND, SBPL, PL, CPL), both bands, 2000 MC-perturbed draws each against a 50,000-point ground truth, 12,000+ evaluations total — found a worst-case relative error at `N_GRID=1000` of `3.67e-8`, six to eight orders of magnitude below any statistical uncertainty this pipeline reports. `N_GRID=5000` only tightened this to `~1e-10`-`1e-11`, i.e. no real gain.

**Landed:** `grb_constants.py`'s shared `N_GRID` reverted `5_000` -> `1_000`, with the convergence numbers recorded in a comment at the constant itself, not just here. `gbm_only_refit.py`'s local override removed, now imports the shared constant (its own value was already `1_000`, so this is a consistency fix with zero numeric change there).

**Full pipeline rerun, all 17 real scripts, verified not just assumed:** `bb_flux_fraction.py` dropped from a projected ~48 min (failed multiprocessing attempt) to **564.3 s (~9.4 min)** with the `N_GRID` fix and batching alone, no parallelism. `lorentz_factor.py`/`lorentz_factor_limit_b.py` dropped to ~3 s each (from minutes). `pe_er_photosphere.py` to 28.1 s. Every other script in the 8-40 s range except `grb_fluence.py` (511.0 s — its bottleneck is `mc_spectra_sampler`'s own repeated `multiprocessing.Pool` spawn overhead, not `N_GRID`, confirmed by it not scaling down proportionally to the other N_GRID-bound scripts) and `amati_relationship.py` (196.7 s). **Total wall time for the entire codebase: 1527.8 s (25.5 min), all 17 scripts exit 0.** Spot-checked 3 output CSVs afterward: seeds and `n_samples` unchanged (expected — seed derivation doesn't depend on `N_GRID`), `f_bb` range (0.050-0.246) and `Gamma_min` range (259.7-507.5) both match the values already documented in `HANDOFF.md` from before this change, confirming the numbers didn't actually move, only the runtime did.

## New tooling: `runner_all.py` / `runner_registry.yaml`

Built at the user's request — "for times like these," i.e. reusable infrastructure for whenever a shared constant changes and the whole pipeline needs rebuilding, rather than re-deriving the script list from memory each time (which is exactly what happened during this session's Phase B work).

- `codes-for-paper/runner_registry.yaml` — the canonical, hand-maintained list of this project's real pipeline scripts (17 entries: the 10 seeded analysis scripts, their 4 `csv_to_latex.py` second-stage table generators, plus `low_index.py`/`high_index.py`/`all-safe-unsafe.py` which have no MC/seed but are still real pipeline scripts), in dependency order (a `csv_to_latex.py` after the script whose CSV it reads). Same transparency convention as `seed_registry.yaml`: what's deliberately excluded (`norris..py`, dormant; `LAT_analysis/`, a different data source entirely, no MC/N_GRID involved) is written down, not silently dropped. *(Historic note: this registry has grown substantially since this log entry — see the current `runner_registry.yaml` for the full, current script list, including the later Norris-fit/lorentz-table generator split.)*
- `codes-for-paper/runner_all.py` — reads the registry, runs each script as a genuine subprocess (`cwd` set to the script's own directory, `PYTHONPATH` set to `src`) — not an in-process import, since several of these scripts execute their logic at bare module level with no `main()`/`if __name__` guard, so importing them in the same process isn't safe. Per-script log files under `runner_all_logs/`, a final summary table with per-script timing, non-zero exit if anything failed. Supports `--only SUBSTRING` (partial pipeline) and `--dry-run` (list without running).

## Notes for future reading

- `MASTER_SEED = 2828702241`, drawn via `secrets.randbits(32)` (OS CSPRNG), 2026-09-02T15:40:46 UTC. Full provenance and the derivation scheme live in `src/grb_research/SEEDING.md` once Phase A lands.
- No commits happen at any point in this work — everything stays in the working tree per this project's standing rule, reviewed by the user before any commit decision is made separately.
- New CSV/plot outputs from Phase B reruns are a real, expected side effect (the seeds change, so MC-derived numbers shift slightly) — not something to be surprised by in `git status` afterward.
