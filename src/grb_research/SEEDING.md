# Seeding scheme for Monte-Carlo draws

## What went wrong

A paper reviewer flagged that Monte-Carlo draws across the codebase reuse the same literal
seed (mostly `12345`, sometimes `1234`/`42`) hardcoded independently in many scripts.
Investigation into that concern found something worse: several scripts reseed an identical
`SEED` constant on *every iteration* of a loop over models/GRBs/episodes. A fresh,
identically-seeded `np.random.Generator` per call means those iterations draw literally
identical underlying random numbers — a correctness bug (silently correlated/duplicated MC
draws across iterations that are supposed to be independent), not just a style concern.
`codes-for-paper/lorentz_factor/lorentz_factor.py` around line 308 was a concrete example of
this pattern (Phase A did not touch it — it was Phase B scope; Phase B has since landed and
fixed it, see "Consumers" below).

## `MASTER_SEED`

```
MASTER_SEED = 2828702241
```

Defined in `grb_constants.py`. Provenance, reproduced verbatim:

```
$ python3 -c "import secrets, datetime; print(secrets.randbits(32)); print(datetime.datetime.now(datetime.timezone.utc).isoformat())"
2828702241
2026-09-02T15:40:46.344318+00:00
```

Drawn via `secrets.randbits(32)` — the OS CSPRNG, not `np.random` — so the master seed itself
does not depend on, or get confused with, any NumPy generator state. This provenance is meant
to also go into the paper's methods section eventually (see "Open items" below).

## Deriving a per-script seed: `seed_from_name`

```python
def seed_from_name(name: str, master_seed: int = MASTER_SEED) -> int:
    digest = hashlib.sha256(f"{os.path.basename(name)}-{master_seed}".encode()).hexdigest()
    return int(digest, 16) % (2**32)
```

Defined in `grb_calculations.py`. Hashes the file's basename together with `MASTER_SEED`
(SHA-256, taken mod `2**32`) to produce a deterministic, script-specific seed, so no script
needs to hardcode its own literal seed value while every script still gets a distinct,
reproducible one.

`os.path.basename` is used rather than the full path so the derived seed is stable across
machines and checkout locations — two different clones of this repo (or the same repo checked
out under a different parent directory) derive the same seed for the same script, since the
full absolute path would differ but the basename would not.

## Usage pattern

Each script builds **one** `rng` object at the top and threads it through every downstream MC
call, rather than re-passing an integer `seed=` at each call site:

```python
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)

# ... later, in a loop over models/GRBs/episodes ...
for item in items:
    result = mc_e_iso_sampler(model, ..., rng=rng)   # NOT seed=SEED
```

This is precisely what fixes the reseed-per-iteration bug: passing `seed=SEED` at each call
site (even with a well-chosen `SEED`) constructs a fresh generator with the same initial state
every time it's called, so iterations that should draw independent samples draw identical ones.
Passing a single already-constructed `rng` instead means the generator's internal state
advances across calls, so each iteration draws a genuinely different sample while the overall
run remains reproducible (same `MASTER_SEED` in, same sequence of draws out).

## Consumers

`seed_from_name` and `MASTER_SEED` are re-exported from `grb_research/__init__.py`. Every script
below builds one `rng` near the top via `SEED = seed_from_name(__file__); rng = get_rng(seed=SEED)`
and threads `rng=rng` through its own MC calls; each derives a distinct seed since `__file__`
differs per script.

- `codes-for-paper/fluence/grb_fluence.py`
- `codes-for-paper/lorentz_factor/lorentz_factor.py`
- `codes-for-paper/lorentz_factor/lorentz_factor_limit_b.py` (independent seed from Limit A —
  deliberately decorrelated; the bit-identical `tau_hat` cross-check between the two does not
  depend on a shared seed, only on both sharing the same closed-form `compute_tau_hat`)
- `codes-for-paper/gbm_only_refit/gbm_only_refit.py`
- `codes-for-paper/bb_fraction/bb_flux_fraction.py`
- `codes-for-paper/photospheric_radius/pe_er_photosphere.py`
- `codes-for-paper/model_parameters/epeak_vs_kt.py`
- `codes-for-paper/model_parameters/peak_energy.py`
- `codes-for-paper/amati_relationship/amati_relationship.py` (`amati_helpers.py` carries the
  `rng`-aware helper functions it calls)
- `codes-for-paper/butterfly_plots/butterfly_all.py`
- `codes-for-paper/variability_analysis/norris..py` — updated for consistency (derived seed
  instead of a literal), but dormant: no caller in `codes-for-paper/` currently invokes it, so
  this one was not run as part of Phase B's verification pass.

All ten runnable scripts above were independently verified deterministic: each run twice from a
clean state, output CSV diffed byte-for-byte between runs, `seed` column confirmed single-valued.
`codes-for-paper/model_parameters/utils.py` is not itself a script, but the shared
`convert_sbpl_to_band`/`extract_kt_epeak_from_models` helpers it exports (consumed by
`epeak_vs_kt.py` and `peak_energy.py`) were fixed in the same pass — `convert_sbpl_to_band`
previously overwrote a caller-supplied `rng` whenever `seed` was also given. Both were later
simplified further (2026-09-03): `seed` was dropped entirely from both signatures and `rng` made
required, since every real caller already passed `rng=` — see "Consumers" above for why `seed` is
kept only at the one root call per script, never on a downstream helper.

## Convention: never silently drop an unused/dormant consumer

When a script that uses this seeding scheme stops being referenced anywhere (e.g. `epeak_vs_kt.py`
isn't cited by any figure or table in the current paper draft) or has no caller at all (e.g.
`norris..py`, dormant), don't just remove it from whatever list is tracking consumers. Keep it,
explicitly categorized (`unused` / `dormant`), with a one-line reason. The concrete instance of
this: `codes-for-paper/seed_registry.yaml` — the source data for the paper's seed-table appendix —
keeps both under those two keys rather than dropping them, so a future reader sees they were
considered and why, not left wondering if they were forgotten. Applies to any future list of this
kind, not just that one file.

## Parallelism does not currently interact with this scheme

`mc_spectra_sampler` (`grb_calculations.py`) is the only place in the codebase that uses
`multiprocessing` (a `Pool` over `n_workers`), and it is safe by construction rather than by
accident: the MC draws (`rng_instance.multivariate_normal(...)`) all happen **before** the pool
is opened, in the parent process, using the single `rng`/`seed` passed in. The pool then
parallelizes `legacy_build_mp` — a deterministic evaluation of the spectral model at each
already-drawn, fixed parameter vector — across workers; no `Generator` object, seed, or further
random draw ever crosses the process boundary. So a shared `rng` is never contended or
duplicated across workers today.

This would stop being true if a future change moved the *sampling* step itself (not just the
post-draw evaluation) inside the `Pool`, e.g. to parallelize MC draws over GRBs/episodes rather
than over already-drawn samples of one model. At that point, do not attempt to share one `rng`
across worker processes (forking/pickling a `Generator` risks each worker starting from an
identical copied state, silently reintroducing correlated draws — the same failure class this
whole scheme exists to fix, one level down). Use `np.random.SeedSequence(MASTER_SEED).spawn(n_workers)`
to hand each worker its own decorrelated child stream instead.

This got a live test 2026-09-03: `bb_fraction/bb_flux_fraction.py`'s episode loop was
parallelized exactly this way (`rng.bit_generator.seed_seq.spawn(n_jobs)`, one child per
worker), and verified deterministic across genuine separate-process reruns. It was reverted
shortly after — not for a correctness problem, but because the win wasn't there in practice
(see `SEED_PLAN-implementation.md`'s "Performance investigation" section): the real fix for
that script's slowness was an unrelated `N_GRID` reduction. The spawn-per-worker pattern itself
worked as designed and remains the right approach if sampling-level parallelism is ever needed
again.

## Paper integration

Done 2026-09-03, in the sibling `GRBResearchPaper` repo: a two-sentence summary of this scheme
(including the `MASTER_SEED` provenance above) in `section-5-data-analysis.tex`, a per-figure
"Monte Carlo seed: N" caption sentence on every genuinely MC-derived figure, and a new "Monte
Carlo seeds" appendix section listing every consumer's seed (`tab:seed_table`, generated by
`codes-for-paper/seed_table_to_latex.py` from `codes-for-paper/seed_registry.yaml` — see that
file for the full script-to-seed mapping, kept independent of this doc so it can't drift).
