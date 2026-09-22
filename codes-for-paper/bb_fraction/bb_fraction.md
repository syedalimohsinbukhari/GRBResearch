# Blackbody fractional flux contribution — method notes

Companion to `bb_flux_fraction.py`. Written for our own understanding: the paper reports the numbers, this file records how they were produced and why each choice was made.

Phase 1 of `PLAN.md`. Produced 2026-08-21.

---

## 1. What is computed

For every interval whose BEST-fitting model includes a blackbody component,

$$f_\text{BB} \;=\; \frac{F^\text{ob}_\text{BB}}{F^\text{ob}_\text{total}}$$

the fraction of the observed energy flux carried by the thermal component.

**Source of the definition.** Pe'er, Ryde, Wijers, Mészáros & Rees (2007), ApJ **664**, L1 — "A new method of determining the initial size and Lorentz factor of GRB fireballs using a thermal emission component". Their §2 defines $F^\text{ob}$ as "the total (thermal + nonthermal) observed γ-ray flux" and $F^\text{ob}_\text{BB}$ as the thermal flux **"integrated over all frequencies"**. Their eq. (5) contains the ratio $F^\text{ob}_\text{BB}/(Y F^\text{ob})$ directly, and their worked example for GRB 970828 quotes

> "the ratio of thermal to total flux $F^\text{ob}_\text{BB}/F^\text{ob} = 0.64 \pm 0.20$ at the break time"

so this is exactly the quantity Phase 2 will need. Computing it in Pe'er's convention now avoids recomputing it later.

**Why this matters for the paper.** Until now the analysis could only say *whether* a BB component was statistically required. $f_\text{BB}$ says *how much* of the burst is thermal, which is the physical payoff, and it is the input to the photospheric radius and Lorentz factor in Phase 2.

---

## 2. Decisions and who made them

### 2.1 Integration band: observer-frame 1 keV – 10 MeV — *Claude, from the Pe'er definition*

Pe'er's $F^\text{ob}_\text{BB}$ is bolometric ("all frequencies"), and his $\mathcal{R} \equiv (F^\text{ob}_\text{BB}/\sigma T^{\text{ob}\,4})^{1/2}$ **requires** it to be bolometric, since $\sigma T^4$ is the all-frequency blackbody flux. A band-limited $F_\text{BB}$ would silently break eq. (1).

Chosen band is 1 keV – 10 MeV in the **observer frame**, matching the bolometric band already used for fluence elsewhere in this work.

*Rejected alternative:* the GBM instrument band (8 keV – 40 MeV). It would avoid extrapolation, but it is an instrument-defined range with no physical meaning, and it is not what Pe'er's equations assume.

*Rejected alternative:* integrating literally over $(0, \infty)$. For some model/parameter combinations $E\,N(E)$ does not converge at both limits (a pure PL diverges), so a finite band is needed for a uniform treatment across all models.

### 2.2 Observer frame, not rest frame — *Claude; consequence of the definition* — **revisited 2026-08-31, see §2.2b**

$f_\text{BB}$ is a ratio of two fluxes over the same band at the same epoch, so **redshift and cosmology cancel exactly and never enter**. This is deliberate and useful: three of the four bursts have no confirmed redshift, and this quantity is still well defined for all of them.

**Consequence for the CSV (superseded).** `CLAUDE.md` requires each row to carry the assumptions used — normally `z`, `H0`, `Om0`. Those were *not* emitted here originally, because emitting them would falsely imply the observer-frame result depends on them. What was emitted instead was the genuine assumption set: `e_min_keV`, `e_max_keV`, `frame`. §2.2b reverses this: the CSV now emits `z`/`H0`/`Om0` as well, since they are genuine inputs to the new rest-frame columns (§2.2b), even though they remain inert for the observer-frame ones.

### 2.2b Cross-burst comparison needs a rest-frame band too — *user, 2026-08-31, resolving weakness #2 of `grb_paper_weaknesses_and_fixes.md`*

§2.2's redshift-cancellation argument is correct **within one episode**, but does not extend to comparing $f_\text{BB}$ **across bursts at different redshifts**: a fixed observer-frame band corresponds to a different rest-frame band per burst, and since the BB peaks near $3.92\,kT$ (Tang 2021, Zhang 2020), this changes how much of the BB's rest-frame shape falls inside vs. outside a shared observer-frame window. This was flagged as a real risk to the paper's headline claim ("GRB231129C is 3-4x more thermally dominated than GRB080916C") since GRB080916C has a measured $z=4.35$ while the other three bursts (including GRB231129C) do not.

**Fix:** also compute $f_\text{BB}^\text{rest}$, integrating both fluxes over a fixed **rest-frame** band 1 keV – 10 MeV (matching the rest-frame band already used for $S_\text{bol}$/$E_\text{iso}$ elsewhere in this work), by evaluating the fitted observer-frame model at the observed energies $E_\text{rest}/(1+z)$. The shift is applied via `redshift_shift = log10(1+z)` on the log-grid bounds — algebraically identical to `mc_e_iso_sampler`'s (1+z) shift (`grb_calculations.py:415-416`) — chosen deliberately to reuse a formula that already has a documented failure mode to avoid: pairing a rest-frame grid with an observed-frame variable underestimated $E_\text{iso}$ by ~3.7x at $z=4.35$ once already (`BUGS.md` BUG-11). Validated here by a z→0 convergence check (the rest-frame and observer-frame grids coincide exactly at z=0) and by confirming `bb_fraction_captured_in_band_rest` stays ≥99.9999% at every z in the sweep (§3.3).

**Redshift for the three unmeasured-z bursts — *user, 2026-08-31*.** Two options were considered: a single fiducial $z=2$ (simpler, matches the weakness doc's literal wording), or a full sweep with $z=2$ marked as fiducial (mirrors the exact precedent this same paper already set for this same problem in `photospheric_radius/pe_er_photosphere.py`: $z$ swept over $[0.5, 5.0]$, 25 log-spaced points, fiducial unioned in, full curve plotted). **User chose the sweep**, since it directly answers "how much does the ranking change" across the full plausible z range rather than at one alternate point, and keeps this exact problem (three of four bursts lacking a spectroscopic z) handled consistently across the two folders that face it. Implemented by copying `REDSHIFTS`/`Z_MIN`/`Z_MAX`/`Z_POINTS`/`Z_FIDUCIAL` verbatim from `pe_er_photosphere.py` (per this project's copy-rather-than-fight-imports convention — these two scripts don't share a package, so the constants are duplicated, not imported).

**Consequence for the CSV.** `z`, `z_source` (`"spectroscopic"`/`"swept"`), `H0`, `Om0`, `e_min_keV_rest`, `e_max_keV_rest`, `f_bb_rest` (+ asymmetric errors), `bb_fraction_captured_in_band_rest` are now emitted, bringing the file into line with the default CLAUDE.md column convention. The CSV now has 263 rows (one per BB-inclusive interval per redshift: 3 spectroscopic-GRB080916C intervals × 1 row + 10 swept intervals × 26 rows = 3 + 260 = 263), up from 13. The observer-frame columns (`f_bb`, `flux_bb_*`, `bb_fraction_captured_in_band`) are unchanged and repeated identically across all rows of a given interval, since they don't depend on z.

**Consequence for the paper claim.** The cross-burst *ranking* (GRB231129C most thermally dominated) is robust to the observer-vs-rest-frame choice, but the exact *multiplier* is not: observer-frame gives roughly 3-4x GRB080916C's thermal fraction, rest-frame gives roughly 2-3x, because GRB080916C's higher redshift ($z=4.35$ vs. the fiducial $z=2$ used for GRB231129C) shifts its rest-frame band further from the fixed observer-frame band, raising its rest-frame $f_\text{BB}$ more than GRB231129C's. **User decided** (over keeping "three to four times" unchanged, or softening to a single "two to four times" range) to report both band's numbers explicitly wherever the claim appears (abstract, §5.2, §6, §7), rather than picking one. See `review-resolution.md` weakness #2 for the exact wording used in each file.

### 2.3 Which intervals — *jointly, from the BEST-model classification*

All intervals whose BEST model name contains `BB`: **13 intervals across all four bursts**.

| GRB | BB-inclusive intervals | of total |
|---|---|---|
| GRB080916C | 3 | 8 |
| GRB131014A | 5 | 5 |
| GRB140206B | 2 | 8 |
| GRB231129C | 3 | 5 |

Note `PLAN.md` originally said "32 intervals, 25 BB-inclusive" — figures that dated from the earlier eight-GRB sample. Logged as `BUGS.md` PLAN-01 and since corrected in `PLAN.md` itself.

### 2.4 Uncertainties: Monte Carlo, 10 000 samples — *existing project convention*

Reuses the project's `n_samples=10_000` convention, the shared `N_SAMPLES` default in `grb_constants.py:103` (imported by `bb_flux_fraction.py`). Note `amati_relationship.py` itself has since diverged to a local `n_sample=5_000` (`amati_relationship.py:41`) — this script was not updated to follow that change, and still uses the shared 10,000 default. **Updated in the RNG/seeding overhaul (Phase B, 2026-09-02):** the seed is no longer the literal `12345` shared by every script, but `seed_from_name(__file__)`, a per-file value deterministically derived from the master seed (`src/grb_research/SEEDING.md`) — this decorrelates each script's draws from every other script's while keeping any single script's own reruns reproducible. **Superseded 2026-09-03**, post-Phase-B cleanup: `compute_fraction`'s `seed` parameter (Phase B first defaulted it to `None`) has since been removed entirely — `rng` is now a required keyword-only argument, since every real call site already passed it and the fallback was dead weight. Same for `draw_model_samples` itself.

Parameters are drawn from the fit covariance and then passed through `ModelResampler.run_resampler()` — the same machinery `mc_spectra_sampler` uses — which enforces the model-specific physical constraints (e.g. the SBPL requirement $\lambda_1 > -2$, $\lambda_2 < -2.05$). Without it, unphysical draws would leak into the percentiles.

RNG is threaded through explicitly, never a global `np.random.seed()`. **Correction, 2026-08-21 → fixed 2026-08-31:** this was true of the *intent* but not the original code — `compute_fraction` built its own `np.random.default_rng(seed)` fallback instead of going through the shared helper, a documentation/implementation mismatch caught while touching this file for §2.2b. Fixed by passing `rng` straight through to `draw_model_samples` instead of constructing a second, redundant generator locally. **As of the 2026-09-03 cleanup**, `draw_model_samples` no longer calls `get_rng` at all — it just uses the `rng` it's given directly, since `get_rng`'s only remaining job project-wide is the one root `seed → rng` conversion at the top of each script (see `SEEDING.md`'s "Usage pattern").

### 2.5 Discarding out-of-range MC draws rather than clipping — *Claude*

$f_\text{BB}$ is physically confined to $[0, 1]$. A small number of draws fall outside because `amp_bb` is drawn from a Gaussian that has support below zero. Those draws are **discarded**, not clipped.

Clipping would pile probability mass onto the boundaries and bias the 16th/84th percentiles inward, making the error bars look artificially tight. Discarding leaves the interval honest. In practice very few draws are affected (the resampler already rejects most).

### 2.6 Vectorised Monte Carlo — *user asked for it after the first version ran too slowly*

The first working version looped over draws in Python, which took several minutes. It now evaluates every draw at once: parameters of shape `(n, 1)` broadcast against energy of shape `(1, n_grid)`, integrated by a single `simpson(..., axis=1)` — the pattern already used in `mc_e_iso_sampler`.

Draws are processed in blocks of `CHUNK_SIZE = 2000` because the un-chunked working array would be $10^4 \times 4000 \times 8$ bytes $\approx 320$ MB *per component*. `N_GRID` was also reduced 4000 → 1000, which is ample for a smooth flux ratio.

Runtime went from several minutes to **28 s**. Measured speed-up was **3.8×** from vectorisation at fixed grid plus ~4× from the grid reduction — worth recording that this is far less than the "one to two orders of magnitude" first estimated, because the 2-D `simpson` call is itself a large part of the cost.

**Results were verified unchanged:** max relative difference in $f_\text{BB}$ of $7\times10^{-9}$, which is Simpson discretisation only and seven orders of magnitude below the statistical error.

This required making `smoothly_broken_power_law` array-safe — it was the only component function with a scalar `if break_energy < 0` guard. That change is bit-identical for scalar callers; see `BUGS.md` BUG-14.

### 2.7 Component decomposition — *Claude; reuses existing machinery*

Parameters are sliced in declaration order exactly as `SpectralModels._evaluate_components` does (`grb_sed.py:96-110`), using `MODEL_MAP` and `model_n_pars`. This guarantees the BB/non-BB split matches the rest of the codebase rather than being a second, independent interpretation of the parameter vector.

---

## 3. Validation

### 3.1 The band really is bolometric — analytic cross-check

For $N(E) = A\,E^2/(e^{E/kT}-1)$ the all-frequency energy flux is analytic:

$$\int_0^\infty E\,N(E)\,\mathrm{d}E \;=\; A\,(kT)^4\,\Gamma(4)\zeta(4) \;=\; A\,(kT)^4\,\frac{\pi^4}{15}$$

Every row reports `bb_fraction_captured_in_band`, the numerically integrated BB flux divided by this analytic value. **Result: 0.999997 – 0.999999 for all 13 intervals** (worst case GRB131014A TR2, $kT = 27.2$ keV).

So the band captures ≥ 99.9997% of the blackbody and $f_\text{BB}$ is bolometric to well below its ~10% statistical uncertainty. This is checked at runtime and written to the CSV rather than assumed once and forgotten.

### 3.1b Rest-frame band capture and shift-formula sanity — added 2026-08-31, §2.2b

`bb_fraction_captured_in_band_rest` (the rest-frame analogue of §3.1) ranges **99.99992% – 99.99999%** across all 263 rows, i.e. every (interval, z) combination in the sweep, including $z=5.0$ (the sweep's upper bound, where the rest-frame band's observer-equivalent upper edge narrows the most, from 10 MeV to $10\,\text{MeV}/6 \approx 1.67$ MeV). The band therefore stays fully bolometric to the BB across the entire redshift range considered, not just at the fiducial $z=2$ — expected, since even at $z=5$ the shifted band's upper edge (~1.67 MeV) sits far above where a $kT\sim27$–$59$ keV blackbody's flux has anything left to contribute.

Two further sanity checks: (1) the rest-frame grid formula (`redshift_shift = log10(1+z)`, applied to the log10 band bounds) is algebraically identical in form to `mc_e_iso_sampler`'s (`grb_calculations.py:415-416`), confirmed by direct code comparison rather than a numeric one, since the two functions integrate different quantities. (2) At $z \to 0$ the rest-frame and observer-frame grids coincide by construction (`redshift_shift \to 0`), so the two bands become identical — not separately re-run, since it follows directly from the formula, but worth stating as the check that would catch a sign error.

### 3.2 Scale sanity against the literature

Pe'er+2007 report $f_\text{BB} = 0.64$ for GRB 970828 — a burst selected *because* it is thermally dominated. Our values, 0.05–0.25, are much smaller, which is the expected direction: our bursts are non-thermally dominated with BB components detected as subdominant additions, described in the paper as "shallow shoulders" rather than resolved thermal peaks. The measured $f_\text{BB}$ is consistent with that description.

---

## 4. Results summary

| GRB | $f_\text{BB}$ range (observer) | $f_\text{BB}^\text{rest}$ range (rest-frame) | note |
|---|---|---|---|
| GRB080916C | 0.050 – 0.068 | 0.074 – 0.099 (at $z=4.35$) | declines from T90 through TR1 |
| GRB131014A | 0.072 – 0.127 | 0.075 – 0.130 (at fiducial $z=2$) | TR2 is the strongest interval |
| GRB140206B | 0.058 – 0.063 | 0.076 – 0.082 (at fiducial $z=2$) | only EX0 and TR1 have a BB |
| GRB231129C | 0.206 – 0.247 | 0.209 – 0.248 (at fiducial $z=2$) | **strongest thermal component in the sample, in both bands** |

Three findings worth carrying into the paper:

1. **GRB231129C is the most thermally dominated burst**, at roughly 3–4× the observer-frame thermal fraction of GRB080916C, or roughly 2–3× in the rest-frame comparison (added 2026-08-31, §2.2b) — the ranking is robust to this choice, the multiplier is not.
2. **GRB140206B's thermal fraction (~6%) is comparable to GRB080916C's (~5–7%)**. This directly contradicts `section-6-discussion.tex:6`, which asserts "the complete absence of thermal components in GRB140206B" — while line 17 of the same subsection correctly says the BB is recovered in TR1 and EX0. See `BUGS.md` BUG-12. The correct statement is that BB is present in **all four** bursts and the distinction is one of degree, not presence.
3. **The rest-frame comparison narrows GRB140206B's gap from GRB131014A's weakest episodes**: observer-frame, GRB140206B (0.058–0.063) sits clearly below GRB131014A's lowest episodes (~0.072); rest-frame, GRB140206B (0.076–0.082) overlaps GRB131014A's lowest episodes (~0.075–0.076) and GRB080916C's T90 (0.099) almost exactly. Not currently stated in the paper text — the existing per-burst strength grouping (§6, "weakly thermal" vs. "intermediate") is unaffected since the grouping is unchanged in both bands, but this specific near-overlap is worth knowing if the discussion is revisited.

---

## 4.5 `kt_bb_keV` here is the raw fit value, not an MC statistic

`bb_flux_fraction.csv`'s `kt_bb_keV`/`kt_bb_err_keV` columns are the point-estimate `kt_bb` and its Gaussian fit error straight from `results.json`, not a Monte Carlo summary. `pe_er_photosphere.csv` (Phase 2) has a same-named `kt_bb_keV` column that **is** an MC median of resampled draws instead, so the two CSVs disagree at the ~0.01–0.02% level despite sharing a column name. Found and confirmed during the Phase 4 spot-check, `BUGS.md` OBS-09 — not a bug, but worth knowing before diffing the two files.

**Surfaced to the reader, 2026-09-01:** `grb_paper_weaknesses_and_fixes.md`'s Priority 3 editorial pass flagged the resulting kT mismatch between `tab:bbfraction` (this table) and `tab:photospheric` (44.45 vs 44.46 keV for GRB080916C T90) as a possible inconsistency. `csv_to_latex.py`'s caption now states explicitly that this table's $kT$ is the raw fit value, not an MC statistic — see `review-resolution.md` Priority 3 item 5.

## 5. Limitations and open questions

- **$f_\text{BB}$ is model-normalised.** Both fluxes come from the fitted model, not from data counts. Same caveat as `BUGS.md` OBS-05 for $E_\text{iso}$; unavoidable for time-resolved episodes where no catalogue flux exists.
- **T90 rows are not independent** of the TR rows that tile them — T90 is the time-integrated spectrum of the same photons. They must not be treated as independent points in any fit or population statistic.
- **The MC assumes a Gaussian parameter posterior** described by the fit covariance. For intervals where `amp_bb` is weakly constrained this is optimistic; those intervals are exactly the ones the SAFE/BEST classification already flags.
- **No Wien-vs-Planck discrimination.** Pe'er+2007 §4 notes that dominant Compton scattering yields a Wien rather than a Planck spectrum, causing a ~5% systematic in $\eta$. The fits assume Planck throughout.

---

---

## Table formatting: normalisation and precision

*Decided by the user, 2026-08-21, after reading the rendered tables.*

Columns are divided by a common power of ten carried in the header rather than repeating `\times 10^{n}` on every row:

| table | column | norm |
|---|---|---|
| Table 4 | $F_\text{BB}$, $F_\text{total}$ | $10^{-6}$ erg cm$^{-2}$ s$^{-1}$ |
| Table 4, Table 5 | $f_\text{BB}$ | $10^{-3}$ |
| Table 5 | $r_0$ | $10^{7}$ cm |
| Table 5 | $r_\text{ph}$ | $10^{11}$ cm |

$kT$ and $\Gamma$ are left unnormalised — both are already of order $10$–$10^{3}$, and normalising them made the numbers *less* readable ($\Gamma = 3.854 \times 10^2$ rather than $385$).

**Decimal places are derived numerically inside `csv_to_latex.py`, not chosen by eye.** The rule, implemented as `decimals_for()`:

$$\text{decimals} = \left\lceil \text{sig} - 1 - \lfloor \log_{10}(\min(\text{error})/\text{norm}) \rfloor \right\rceil, \quad \text{sig} = 2$$

i.e. enough decimals that the *smallest* error in the column still shows two significant figures. It is computed from the dataframe at generation time, so the formatting adapts if the data changes instead of going stale.

This was worth doing numerically: eyeballing had produced 3 decimals where 1 was right ($f_\text{BB}$, whose errors are $\sim5$ in units of $10^{-3}$) and 4 decimals on the Amati table where the errors are of order 7. The same helper now runs in all three `csv_to_latex.py` scripts.

## 6. Files

| file | role |
|---|---|
| `bb_flux_fraction.py` | computation, CSV and figures |
| `bb_flux_fraction.csv` | one self-describing row per BB-inclusive interval per redshift (263 rows) |
| `bb_flux_fraction.png` / `.pdf` | $f_\text{BB}$ (observer-frame) vs time, one panel per GRB |
| `bb_flux_fraction_rest_vs_z.png` / `.pdf` | $f_\text{BB}^\text{rest}$ vs redshift, added 2026-08-31 §2.2b |
| `csv_to_latex.py` | renders `bb_flux_fraction.csv` into the paper's `\input`-ed table |
| `bb_fraction_table.tex` | generated table (also copied to `GRBResearchPaper/tex_files/generated/`) |
| `bb_fraction.md` | this file |
