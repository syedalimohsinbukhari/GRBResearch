# Bug log — found while executing PLAN.md

Maintained by Claude while working through the phased plan. Newest entry at the bottom, so IDs are not in numerical order.

Status key: **OPEN** (needs your decision), **FIXED** (changed in the working tree), **RESOLVED** (superseded by another entry), **NOTED** (recorded, no change required), **DEFERRED** (tracked for a later phase).

Entries that were fixed later carry the resolution first, then **Original diagnosis:** and the reasoning as it stood when found.

---

<!-- INDEX:BEGIN — generated from the ### headings below; regenerate rather than hand-edit. -->
## Index

35 entries: 23 bugs, 11 observations, 1 plan note. **All entries closed or resolved except OBS-11 (DEFERRED, user's call) and BUG-23 (FIXED for GRB080916C only; confirmed present and unfixed for GRB131014A/GRB140206B).**

**Everything else, in discovery order** (the order below is the order found, not ID order):

| ID | Entry | Status |
|---|---|---|
| `BUG-01` | `mc_e_iso_sampler` integrates the luminosity distance instead of evaluating it | FIXED |
| `BUG-02` | `f_1` for SBPL models used the high-energy asymptote, not the actual spectrum | FIXED |
| `BUG-03` | $\Gamma_\text{min}$ truncated rather than rounded in the LaTeX table | FIXED |
| `BUG-04` | wrong burst name in the LAT appendix | FIXED |
| `BUG-05` | `\Gamma_\min` is a LaTeX compile error | FIXED |
| `BUG-06` | cosmology mismatch | FIXED (planned, PLAN.md Phase 0.3) |
| `BUG-07` | stale GRB sample in `lorentz_factor.py` | FIXED (planned, PLAN.md Phase 0.4) |
| `BUG-08` | appendix prose asserts BB findings about bursts not in the sample | FIXED in Phase 3 |
| `BUG-09` | `csv_to_latex.py` crashed on its own CSV | FIXED |
| `BUG-10` | the paper's $E_\text{iso}$ table was hand-copied and had drifted from the pipeline | FIXED |
| `BUG-11` | k-correction evaluated the spectrum on the rest-frame grid but integrated over observed energies | FIXED |
| `OBS-03` | absolute $E_\text{iso}$ scale | RESOLVED by BUG-11 |
| `BUG-12` | §6 contradicts itself and the data on GRB140206B's blackbody | FIXED in Phase 3 |
| `BUG-13` | `GRBPlotStyle` had no colours for three of the paper's four GRBs | FIXED |
| `OBS-06` | `bb_flux_fraction.py`'s Monte Carlo loop is not vectorised | FIXED |
| `BUG-14` | `smoothly_broken_power_law` could not accept array parameters | FIXED |
| `PLAN-01` | PLAN.md's interval counts were wrong for this sample | FIXED |
| `OBS-04` | photospheric-radius claims in the prose are cited, not computed | FIXED in Phase 3 |
| `OBS-05` | $S_\text{obs}$ is model-anchored, not data-anchored, so the k-correction cancels by construction | NOTED, deliberate |
| `OBS-01` | inconsistent package import style | FIXED |
| `OBS-02` | `mc_e_iso_sampler` `method=1` omits the interval duration | RESOLVED by BUG-11 |
| `BUG-17` | generated tables dropped uncertainties on their input columns | FIXED |
| `BUG-16` | LAT appendix photon-energy header changed MeV to keV | REVERTED, the original was correct |
| `OBS-08` | per-episode LAT photon data is hand-transcribed from the paper table | FIXED |
| `OBS-07` | $\Gamma_\text{min}$ is computed only for T90, blocking a clean comparison with the thermal $\Gamma$ | FIXED |
| `BUG-15` | Phase 2 figure did not distinguish episodes; legend handles lacked markers | FIXED |
| `BUG-18` | LAT photon energies are MeV but $\Gamma_\text{min}$ treated them as keV | FIXED |
| `OBS-09` | `kt_bb_keV` means different things in the two Phase 1/2 CSVs | NOTED, deliberate but undocumented until now |
| `OBS-10` | `bb_flux_fraction.py`'s RNG bypassed `get_rng` despite `bb_fraction.md` claiming otherwise | FIXED |
| `BUG-19` | `convert_sbpl_to_band` silently overwrote a caller-supplied `rng` whenever `seed` was also given | FIXED |
| `BUG-20` | `FluxFluenceCalculator.__init__`'s falsy `if seed:` check silently dropped `seed=0` | FIXED |
| `BUG-21` | `codes-for-paper/variability_analysis/fitter*.py` import a `variability_timescale` module that no longer exists | FIXED |
| `BUG-22` | `light_curves.py`'s `lightcurve_data()` combined per-channel errors linearly instead of in quadrature | FIXED |
| `OBS-11` | GRB231129C's zenith-cut boundary (~100.2–100.3°) exceeds the stated 100° cut, exposure-loss correction unverified | DEFERRED, user's call |
| `BUG-23` | GRB080916C's `fitter.py`/`window_sensitivity.py` silently picked `n4` instead of the documented `n3` NaI detector after a data resync; same bug confirmed live and unfixed for GRB131014A/GRB140206B | FIXED (GRB080916C only) |

<!-- INDEX:END -->

---

### BUG-01 — `mc_e_iso_sampler` integrates the luminosity distance instead of evaluating it — **FIXED**

**Where:** `src/grb_research/grb_calculations.py:432-433`

```python
lum_distance = lambda z: FlatLambdaCDM(h0, omega_m).luminosity_distance(z).cgs.value
lum_distance = quad(lum_distance, 0, z)[0]
return 4 * np.pi * lum_distance ** 2 * ... / (1 + z)
```

`astropy`'s `luminosity_distance(z)` **already returns** $d_L(z)$. Wrapping it in `quad(..., 0, z)` computes $\int_0^z d_L(z')\,dz'$, which is not a distance at that redshift — it is an integral of one. The standard relation
$E_\text{iso} = 4\pi d_L^2 F_\text{bol}/(1+z)$ requires $d_L(z)$ itself.

**Impact, measured at $z = 4.35$ (GRB080916C) with $H_0=69.6,\ \Omega_m=0.286$:**

| quantity | value |
|---|---|
| correct $d_L(4.35)$ | $1.2478\times10^{29}$ cm |
| as currently coded | $2.4592\times10^{29}$ cm |
| ratio | 1.971 |
| **error in $E_\text{iso}$** (scales as $d_L^2$) | **×3.88 too large** |

The error is redshift-dependent, so it does not cancel as a constant offset — it distorts the *shape* of the Amati correlation across the swept-$z$ bursts, not just its normalization.

**Blast radius:** `codes-for-paper/amati_relationship/amati_helpers.py:92` is the only caller, but it feeds `amati_relationship.csv`, `amati_relationship.png/pdf`, and `amati_relationship_table.tex` — all of which appear in the paper. Every published $E_\text{iso}$ is affected.

**Fix:** drop the `quad` wrapper and evaluate directly:
```python
lum_distance = FlatLambdaCDM(h0, omega_m).luminosity_distance(z).cgs.value
```

**Applied.** `quad` was then unused and was dropped from the imports. Amati outputs regenerated (`amati_relationship.csv`, `.png`, `.pdf`, `amati_relationship_table.tex`); the figure was copied into `GRBResearchPaper/images/section4`.

**Verified:** the stale `amati_relationship_table.tex` held E_iso values *exactly* $3.8844\times$ the regenerated ones — matching the predicted $(d_{L,\text{old}}/d_{L,\text{new}})^2 = 3.8843$ to five digits.

**§6 claim re-checked:** the discussion states that all GRB080916C episodes lie within the $3\sigma$ Amati band. Inspected the regenerated figure — **still true** after the correction, so no prose change was needed.

---

### BUG-02 — `f_1` for SBPL models used the high-energy asymptote, not the actual spectrum — **FIXED**

**Where:** `GRBResearchWork/codes-for-paper/lorentz_factor/lorentz_factor.py`, former `f1_from_sbpl()`

```python
N_1MeV = amp * (1000.0 / e_piv) ** index2   # WRONG for SBPL
```

This extrapolates from the pivot (100 keV) to 1 MeV using the **steep high-energy index** `index2`. But in the SBPL parameterisation (`grb_seds.py:17-30`) `amp` is the flux *at the pivot*, where the spectrum follows the **shallow** `index1` — the steep index only takes over above `e_break` (702 keV for GRB080916C). Extrapolating with `index2` from 100 keV therefore under-counts the flux at 1 MeV.

Note the same script evaluated the *Band* model correctly (`f1_from_band` called the real Band function), so the two model paths disagreed with each other.

**Impact for GRB080916C (T90, SBPL_BB):**

| quantity | old | corrected |
|---|---|---|
| $f_1$ [ph cm⁻² s⁻¹ MeV⁻¹] | $5.8686\times10^{-2}$ | $3.4483\times10^{-1}$ |
| $\hat\tau$ | $8.460\times10^{10}$ | $4.971\times10^{11}$ |
| $\Gamma_\text{min}$ | **197** | **258** |

Verified three ways: direct evaluation of `smoothly_broken_power_law`, log-log interpolation on a 20 000-point grid, and a by-hand evaluation of the SBPL formula ($3.4497\times10^{-1}$, agreeing to 0.04%).

$\Gamma_\text{min}$ is weakly sensitive to $f_1$ (exponent $1/(2\alpha+2) \approx 0.154$), which is why a 5.9× error in $f_1$ becomes only a 1.31× change in $\Gamma_\text{min}$.

**Consequence for the paper — needs your decision:** `tex_files/section-7-conclusion.tex:19` and the (currently commented) `tex_files/section-5-data-analysis.tex:191` both state $\Gamma_\text{min} = 197$. The regenerated `lorentz_table.tex` now says 258. **These now disagree and must be reconciled before the paper is final.**

---

### BUG-03 — $\Gamma_\text{min}$ truncated rather than rounded in the LaTeX table — **FIXED**

**Where:** `GRBResearchWork/codes-for-paper/lorentz_factor/lorentz_factor.py`, table generation

```python
gamma_str = f"${int(r['Gamma_min'])}$"
```

`int()` truncates: the console printed `197` (via `.0f`, which rounds) while the LaTeX table emitted `196` for the same 196.657 value — so the script silently contradicted itself. Now uses `round()`.

---

### BUG-04 — wrong burst name in the LAT appendix — **FIXED**

**Where:** `GRBResearchPaper/appendices/appendix_LAT_info.tex:30` (that file was hand-typed at the time; superseded by the auto-generated `LAT_analysis/csv_to_latex.py` → `tex_files/generated/lat_info_table.tex`, and the wrapper file itself was deleted 2026-09-06 — see BUG-10/OBS-08 below and `CLAUDE.md`'s "Generated LaTeX tables")

The time-integrated row read `GRB140206A`, while the time-resolved block for the same burst uses `\grbfourteenzerotwozerosixB` (GRB140206**B**) and the interval `7.488–154.240 s` matches `GRB140206275` in `GRBResearchWork/results.json`. A plain typo, but it named a *different real burst*. All four time-integrated rows now use the `\grb...` macros.

---

### BUG-05 — `\Gamma_\min` is a LaTeX compile error — **FIXED**

**Where:** `tex_files/section-7-conclusion.tex:19` (live) and 8 occurrences in `tex_files/section-5-data-analysis.tex` (inside the commented γγ-opacity block)

`\min` is a `\mathop`, so `\Gamma_\min` raises `! Missing { inserted.` Re-enabling section 7 surfaced this. All occurrences are now braced as `\Gamma_{\min}`, including the commented ones so Phase 3 does not re-break the build.

---

### BUG-06 — cosmology mismatch — **FIXED** (planned, PLAN.md Phase 0.3)

`lorentz_factor.py` used $H_0=67.4,\ \Omega_m=0.315$ (Planck 2018) while the paper (`section-5-data-analysis.tex:77`) and `mc_e_iso_sampler` (`grb_calculations.py:342-343`) both use $H_0=69.6,\ \Omega_m=0.286$ (Fana Dirirsa 2019). Reconciled to 69.6/0.286 everywhere; all three sites now agree. $\Gamma_\text{min}$ is insensitive to this change (it did not move the rounded value).

---

### BUG-07 — stale GRB sample in `lorentz_factor.py` — **FIXED** (planned, PLAN.md Phase 0.4)

`grb_inputs` listed GRB110721A, GRB110731A and GRB150210A — three bursts from a previous project that are not in this paper's sample — alongside GRB080916C. Replaced with the current four. The script now derives its spectral parameters from `GRBResearchWork/results.json` via the `grb_research` class API, so this class of drift cannot recur: only the LAT photon properties and redshifts remain hardcoded, because they are not in `GRBResearchWork/results.json`.

Removed with it: the `$^\dagger$` AGN-contamination footnote, which referred solely to GRB110721A.

---

### BUG-08 — appendix prose asserts BB findings about bursts not in the sample — **FIXED in Phase 3**

**Where:** `GRBResearchPaper/main.tex`, the `%%`-commented paragraph above the appendix tables

**Original diagnosis:**

The text claims BB components are most often BEST in early episodes of "GRB 080916C and GRB 110721A", and that "GRB 110731A and GRB 150210A show no episode requiring a BB component". Three of those four bursts are not in this paper.

Left commented out with a `TODO (Phase 3)` marker rather than rewritten, because the corrected statement is an empirical claim that needed the measured $f_\text{BB}$ values from Phase 1.

**Resolved.** With $f_\text{BB}$ in hand the paragraph was rewritten against the real sample and re-enabled: a BB component is now correctly stated to be BEST in at least one episode of *every* burst, with $kT \sim 27$–$59$ keV, persisting throughout GRB131014A and GRB231129C but confined to the earliest intervals of GRB080916C and GRB140206B.

---

### BUG-09 — `csv_to_latex.py` crashed on its own CSV — **FIXED**

**Where:** `codes-for-paper/amati_relationship/csv_to_latex.py:111-113`

Looked up `E_0_iso__1e52_erg`, but `amati_relationship.py` writes `E_0_iso__1e+52_erg` (the column name is built with `f"{EI_NORM:.0e}"`, which emits `1e+52`). Every run raised `KeyError`, so the table could not be regenerated at all — which is *why* the committed `.tex` had drifted so far from the CSV.

Also fixed in the same pass, since they made the output unusable as a paper table:
- `\renewcommand{\arraystretch{1.25}}` — misplaced brace, should be `\renewcommand{\arraystretch}{1.25}`;
- two `\label`s on one table (`tab:eiso` *and* `tab:amati_relationship`) — the second silently won;
- empty `\caption{}`;
- model/burst names emitted as literal text (`SBPL+BB`, `GRB080916C`) instead of the paper's `\sbplbb` / `\grbzeroeightzeroninesixteenC` macros.

---

### BUG-10 — the paper's $E_\text{iso}$ table was hand-copied and had drifted from the pipeline — **FIXED**

**Where:** `GRBResearchPaper/tex_files/section-5-data-analysis.tex`, former hardcoded `tab:eiso`

The table was typed into the `.tex` rather than generated. Its values matched **no** state of the code: they were a uniform $\approx 7.6\%$ below the corrected pipeline output (ratios 1.0717–1.0788 across all eight rows — too systematic for MC noise, and not reproduced by any candidate cosmology: 67.4/0.315 gives 1.0021, 70/0.3 gives 0.9588, 72/0.3 gives 0.9063 against an observed 0.9288).

So the published $E_\text{iso}$ values came from a code state that no longer exists and cannot be reproduced. Replaced with `\input{tex_files/generated/amati_relationship_table}`, generated from the CSV.

**New convention:** generated tables now live in `GRBResearchPaper/tex_files/generated` and are `\input`, never retyped. Currently holds `amati_relationship_table.tex` and `lorentz_table.tex`. Regenerating requires re-copying the `.tex` from `GRBResearchWork` — the two repos are independent, so this copy step is manual.

---

### OBS-02 — `mc_e_iso_sampler` `method=1` omits the interval duration — **RESOLVED by BUG-11**

**This entry is stale.** Its index status was never updated after BUG-11's fix; re-verified 2026-08-21 at the user's prompt and the concern no longer applies.

**Where:** `src/grb_research/grb_calculations.py:418-437` (current line numbers, post BUG-11)

BUG-11's fix pass ("`method=1` completed at the same time … the missing `* model.interval.duration` was restored") already addressed this. Current code:

```python
s_obs = detector_flux * model.interval.duration        # line 435 — the duration multiply is present
k_correction = bolometric_flux / detector_flux          # line 436 — same units (keV/cm^2), so a dimensionless ratio
bolometric_fluence = np.asarray(s_obs * k_correction, dtype=float) * kev_to_erg   # line 437 — converted to erg once
```

`method=1` now agrees with `method=2` to machine precision (BUG-11), so it is a genuine independent cross-check, not a latent luminosity bug.

**Original diagnosis (superseded):**

**Where:** `src/grb_research/grb_calculations.py:418-420`

```python
detector_fluence = simpson(...)  # * model.interval.duration
```

The duration multiplication is commented out, so `method=1` produces a **flux**, not a fluence — making its "E_iso" a luminosity, low by a factor of the interval duration (≈63 s for GRB080916C's T90). `method=2` (line 430) does multiply correctly.

**Not currently affecting any paper number:** `amati_helpers.py:92` passes `method=2`. Flagged so nobody switches to `method=1` assuming the two are equivalent. Left untouched — I could not tell whether `method=1` is intended as a flux diagnostic or is simply broken.

---

### BUG-11 — k-correction evaluated the spectrum on the rest-frame grid but integrated over observed energies — **FIXED**

**Where:** `GRBResearchWork/src/grb_research/grb_calculations.py`, `mc_e_iso_sampler` (was lines 391-430)

This supersedes OBS-03 below, which is now resolved.

```python
energy_bolometric = np.logspace(bol_min, bol_max, n_grid)   # REST-frame band, 1 keV-10 MeV
e_observed = energy_bolometric / (1 + z)
bolometric_samples = mc_spectra_sampler(..., e_range=(bol_min, bol_max), ...)  # evaluated at E_rest
...
bolometric_fluence = simpson(y=bolometric_samples, x=e_observed, ...)          # integrated over E_obs
```

The model is the **observer-frame** fitted spectrum, but it was evaluated at rest-frame energy values while being integrated against observed ones. Substituting $w=(1+z)u$ shows what `method=2` actually computed:

$$\int E_r N(E_r)\,\frac{dE_r}{1+z} = \frac{1}{1+z}\int_{1\,\text{keV}}^{10\,\text{MeV}} E\,N(E)\,dE$$

i.e. the band limits were left unshifted and a spurious $1/(1+z)$ was introduced. Confirmed numerically: the as-coded result matched this closed form to 8 decimal places.

The correct Bloom+2001 bolometric fluence integrates the observed spectrum over the *observed* energies corresponding to the rest-frame band:

$$S_\text{bol}=\int_{E_1/(1+z)}^{E_2/(1+z)} E\,N(E)\,dE ,\qquad E_\text{iso}=\frac{4\pi d_L^2\,S_\text{bol}}{1+z}$$

**Impact at $z=4.35$: $E_\text{iso}$ was $3.696\times$ too small.**

**Fix:** shift the exponents rather than dividing the grid, so the grid the model is evaluated on is the same grid that is integrated over — `logspace(a,b)/(1+z) == logspace(a-s, b-s)` with `s = log10(1+z)`.

**`method=1` completed at the same time** (the author confirmed it had been left unfinished): the missing `* model.interval.duration` was restored, and the spurious extra `* energy` factors removed — with `model_type="energy"` the sampler already returns $E\,N(E)$, so multiplying by $E$ again was integrating $E^2N(E)$, which is not a fluence.

**Both methods now agree to machine precision** (max element-wise relative difference $3.9\times10^{-16}$), which is the expected algebraic identity: the detector term in $S_\text{obs}\times k$ cancels exactly, so the k-correction route *must* reduce to the direct integral. `method=1` is now a genuine independent cross-check rather than a no-op.

**Independent validation.** GRB080916C T90:

| band | $E_\text{iso}$ |
|---|---|
| 1 keV–10 MeV rest (as coded, before) | $9.84\times10^{53}$ erg |
| 1 keV–10 MeV rest (corrected) | $3.61\times10^{54}$ erg |
| 1 keV–10 GeV rest (corrected, Abdo-like band) | $7.58\times10^{54}$ erg |
| **Abdo et al. 2009 published** | $\approx 8.8\times10^{54}$ erg |

Agreement to ~14% over the comparable band, versus roughly an order of magnitude before. On the regenerated Amati figure all eight GRB080916C episodes now straddle the central relation instead of sitting low between $2\sigma$ and $3\sigma$; §6's "within $3\sigma$" statement remains true and is now stronger.

---

### ~~OBS-03~~ — absolute $E_\text{iso}$ scale — **RESOLVED by BUG-11**

Originally raised because the corrected-$d_L$ value ($9.8\times10^{53}$ erg) still sat about an order of magnitude below published estimates. The cause was the frame mismatch in the k-correction; see BUG-11 for the diagnosis, fix and validation.

---

### BUG-12 — §6 contradicts itself and the data on GRB140206B's blackbody — **FIXED in Phase 3**

**Resolved.** §6 now states that a blackbody is required in at least one episode of all four bursts, and separates the sample along two axes — thermal *strength* ($f_\text{BB}$) and thermal *persistence*. GRB140206B is described as weakly thermal but short-lived ($f_\text{BB}\approx0.058$–$0.063$, indistinguishable from GRB080916C's $0.050$–$0.068$), not as thermally absent. The "three of our four bursts" count in §6 and the omission of GRB140206B from §7 were corrected at the same time.

**Original diagnosis:**

`tex_files/section-6-discussion.tex:6` states "the **complete absence** of thermal components in \grbfourteenzerotwozerosixB". Line 17 of the *same subsection* then says "The BB component is recovered only in TR1 and its corresponding excess episode EX0" — which is about that same burst.

The data supports line 17. From `best_names.txt` and `GRBResearchWork/results.json`, GRB140206B has BB-inclusive BEST models in **two** intervals:

| interval | BEST model | $kT$ [keV] |
|---|---|---|
| EX0 (4.288–11.072 s) | `BAND_BB` | 58.51 |
| TR1 (7.488–11.072 s) | `BAND_BB` | 55.84 |

`tex_files/section-4-joint-analysis-results.tex:95` already states this correctly. So line 6 is the outlier.

Knock-on errors from the same mistake:
- `section-6-discussion.tex:20` — "across **three** of our four bursts" should be four; every burst in the sample has at least one BB-BEST interval.
- `section-7-conclusion.tex:8` — lists BB detections in GRB080916C, GRB131014A and GRB231129C only, omitting GRB140206B.

Deferred rather than patched now because the replacement sentence should quantify *how weak* the thermal component is in GRB140206B, which is exactly what the Phase 1 $f_\text{BB}$ values will provide. The likely correct framing is "present in all four bursts but confined to the earliest episodes and contributing only $f_\text{BB}\sim X\%$ in GRB140206B", not a binary present/absent.

---

### BUG-13 — `GRBPlotStyle` had no colours for three of the paper's four GRBs — **FIXED**

`GRBResearchWork/src/grb_research/grb_styles.py`'s `GRB_COLORS` and `GRB_SHORT` listed GRB080916C, GRB110721A, GRB110731A, GRB150210A and GRB190114C — i.e. the *old* sample. GRB131014A, GRB140206B and GRB231129C were absent, so any figure following the mandated `GRBPlotStyle` convention would `KeyError` on three quarters of the current sample.

Added the three missing bursts, assigning the four most mutually distinct Okabe–Ito hues to the current sample (blue / orange / green / purple). The earlier-sample entries are kept and deliberately reuse two of those hues, since no figure mixes the two disjoint samples.

---

### OBS-06 — `bb_flux_fraction.py`'s Monte Carlo loop is not vectorised — **FIXED**

**Outcome.** Whole script went from several minutes to **28 s**, with results unchanged to 1 part in $10^8$ (max relative difference in $f_\text{BB}$ = $7\times10^{-9}$, pure Simpson discretisation — seven orders of magnitude below the ~10% statistical error).

Measured speed-up, benchmarked on GRB231129C T90 with 10 000 draws:

| grid | vectorised | scalar loop | speed-up |
|---|---|---|---|
| `n_grid=1000` | 1.80 s | 6.87 s | **3.8×** |
| `n_grid=4000` | 7.87 s | 11.92 s | 1.5× |

Combined with the grid reduction (4000 → 1000) the overall gain is roughly **15×**. Note this is well short of the "one to two orders of magnitude" estimated below — the 2-D `simpson` call and the large temporaries are themselves a substantial cost, so removing the Python loop recovers less than the loop overhead suggested. Recording the miss rather than quietly restating the estimate.

Enabling this required making `smoothly_broken_power_law` array-safe (see below).

**Original diagnosis:**

The MC loop is scalar over draws:

```python
for i, values in enumerate(samples):
    bb_i, total_i = split_energy_flux(model.name, values, energy)
```

13 intervals × 10 000 samples = 130 000 iterations, each evaluating 2–3 component functions on a 4 000-point grid and making 2 `simpson` calls — roughly 260 000 SciPy calls at Python-loop overhead. This is why the script takes minutes rather than seconds.

The SED functions are elementwise in energy, so parameters shaped `(n_samples, 1)` broadcast against energy shaped `(1, n_grid)` to give `(n_samples, n_grid)` in a single call, followed by one `simpson(..., axis=1)`. That is already the established pattern in `mc_e_iso_sampler` (`axis=1`). Expected speed-up is one to two orders of magnitude.

`N_GRID = 4000` is also generous for a smooth flux ratio; 1000 would be ample. Worth doing before Phase 2, which runs the same machinery across a redshift grid and will multiply the cost.

---

### BUG-14 — `smoothly_broken_power_law` could not accept array parameters — **FIXED**

`src/grb_research/grb_seds.py:17` guarded the break energy with a scalar test:

```python
if break_energy < 0:
    return np.full(shape=energy.shape, fill_value=np.nan)
```

With an array of Monte-Carlo draws this raises `ValueError: The truth value of an array with more than one element is ambiguous`, so SBPL was the *only* one of the five component functions that could not be vectorised (PL, CPL, Band and BB all broadcast correctly). Replaced with a mask, using a placeholder break energy so the logarithm of a non-positive number is never evaluated.

**Verified behaviour-preserving:** for scalar parameters the output is *bit-identical* to the previous implementation (max absolute difference exactly 0.0 over a 2001-point grid), a negative break energy still returns all-NaN with the same shape, and vectorised rows match the corresponding scalar calls to $10^{-15}$.

One deliberate semantic change: the invalid test is now `break_energy <= 0` rather than `< 0`. A break energy of exactly zero previously produced `inf`/`nan` from `log10(E/0)` together with a runtime warning; it now returns NaN cleanly.

---

### PLAN-01 — PLAN.md's interval counts were wrong for this sample — **FIXED**

PLAN.md stated "32 intervals, 25 of them BB-inclusive". The actual BEST-model set for the four-GRB sample is **26 intervals, 13 BB-inclusive**:

| GRB | intervals | BB-inclusive |
|---|---|---|
| GRB080916C | 8 | 3 |
| GRB131014A | 5 | 5 |
| GRB140206B | 8 | 2 |
| GRB231129C | 5 | 3 |
| **total** | **26** | **13** |

The 32/25 figures date from the older eight-GRB set. Phase 1 processes the 13 above.

**Corrected in `PLAN.md` on 2026-08-21**, so the plan and the data now agree; the original wording is preserved here.

---

### OBS-05 — $S_\text{obs}$ is model-anchored, not data-anchored, so the k-correction cancels by construction — **NOTED, deliberate**

**User confirmed 2026-08-21: model-anchored is the intended design**, not an open question. Already stated in §5 (per HANDOFF.md §5); no further action.

Capozziello \& Izzo write the k-correction as

$$S_\text{bol} = S_\text{obs}\,\frac{\int_{1/(1+z)}^{10^4/(1+z)} E\phi\,dE}{\int_{E_\text{min}}^{E_\text{max}} E\phi\,dE}$$

with $S_\text{obs}$ defined as *"the fluence observed for each GRB in a respective detection band"* — an empirical catalog quantity. The denominator is the fitted model over that band. The ratio is therefore a pure **shape** factor: it divides out the model normalisation so the absolute scale is inherited from the data.

In `mc_e_iso_sampler` both terms are computed from the same model, so they cancel identically and `method=1` reduces to `method=2` (verified to $3.9\times10^{-16}$). The identity is real, but it follows from our choice, not from the formula — and it means `method=1` is **not an independent physical check**, only an algebraic one.

Size of the difference for GRB080916C T90:

| band | model-integrated | published |
|---|---|---|
| 8 keV–40 MeV (GBM) | $1.67\times10^{-4}$ erg/cm² | $\approx 2.4\times10^{-4}$ (Abdo+2009) |

A data-anchored $S_\text{obs}$ would raise $E_\text{iso}$ by $\approx 1.4\times$.

**Recommendation: keep model-anchored.** No measured catalog fluence exists per time-resolved episode, so anchoring to data is only possible for T90; mixing the two would make the Amati points internally inconsistent. But this should be stated explicitly in §5 rather than left implicit.

**Latent trap:** `det_min=1.0, det_max=7.0` is 10 keV–10 GeV, which is *not* GBM's 8 keV–40 MeV detection band. Irrelevant while the term cancels; wrong immediately if a measured $S_\text{obs}$ is ever substituted.

---

### OBS-04 — photospheric-radius claims in the prose are cited, not computed — **FIXED in Phase 3**

**Resolved.** §6 and §7 now quote our own $r_\text{ph}$ values, and §5 states plainly that GRB080916C and GRB140206B fall inside the literature $10^{11}$–$10^{12}$ cm range while GRB131014A and GRB231129C come out an order of magnitude higher — with the caveat that the latter two are at an assumed redshift.

**Original diagnosis:**

`tex_files/section-6-discussion.tex:9` and `tex_files/section-7-conclusion.tex:8` both state that the measured $kT \sim 30$–$70$ keV are "consistent with photospheric radii of order $10^{11}$–$10^{12}$ cm", attributed to \citet{Ryde2010, Guiriec2011}.

These are borrowed literature values, not results of this analysis. Phase 2 computes $R_0$ and $R_\text{ph}$ directly via the Pe'er method, so both sentences must be checked against our own numbers and rewritten to cite them. If our derived radii disagree with $10^{11}$–$10^{12}$ cm, that is a finding, not something to paper over.

**Now answered.** Phase 2 computes $r_\text{ph}$ directly. At the measured or fiducial redshift, $Y=1$:

| GRB | $r_\text{ph}$ [cm] | inside the claimed $10^{11}$–$10^{12}$? |
|---|---|---|
| GRB080916C ($z=4.35$) | $6.1\times10^{11}$ – $1.4\times10^{12}$ | yes |
| GRB140206B ($z=2^\dagger$) | $5.6$ – $7.5\times10^{11}$ | yes |
| GRB131014A ($z=2^\dagger$) | $3.4$ – $7.6\times10^{12}$ | **no — an order of magnitude higher** |
| GRB231129C ($z=2^\dagger$) | $3.5$ – $4.8\times10^{12}$ | **no** |

So the borrowed range is right for two bursts and wrong for two. §6:9 and §7:8 must be rewritten against these numbers in Phase 3 rather than continuing to cite Ryde 2010 / Guiriec 2011 for a range our own data does not reproduce. Note the two discrepant bursts are the ones at an *assumed* redshift, so the statement should carry that caveat.

---

### BUG-17 — generated tables dropped uncertainties on their input columns — **FIXED**

Caught by the user reading the rendered Table 5.

In `photospheric_table.tex` the *derived* quantities ($r_0$, $\Gamma$, $r_\text{ph}$) carried asymmetric $1\sigma$ errors, but the two *input* columns did not: $kT$ printed as `44.46` and $f_\text{BB}$ as `0.0684`, both bare medians. Both have Monte Carlo uncertainties — $f_\text{BB}$ was already reported with errors in Table 4 — they simply were not being summarised or emitted. `pe_er_photosphere.py` now percentile-summarises $kT$ and $f_\text{BB}$ from the same draws as everything else, the CSV gained `kt_bb_err_lower/upper_keV` and `f_bb_err_lower/upper`, and every column of the table now carries an interval.

Related, from the same review: **Table 4 reported only the ratio $f_\text{BB}$, never the absolute fluxes.** The user asked whether ratios such as $f_\text{MODEL+BB}$ or $f_\text{BB}/f_\text{MODEL}$ should also appear. They should not — $f_\text{MODEL} = 1 - f_\text{BB}$ and $f_\text{BB}/f_\text{MODEL} = f_\text{BB}/(1-f_\text{BB})$ are exact transforms of a column already present, so they would restate the same number. But the question exposed a genuine gap: from a ratio alone the absolute thermal flux cannot be recovered, and $F_\text{BB}$ spans $1.6\times10^{-7}$ to $7.7\times10^{-6}$ erg cm$^{-2}$ s$^{-1}$ across the sample. `bb_fraction_table.tex` now reports $F_\text{BB}$ and $F_\text{total}$ alongside $f_\text{BB}$, from which any ratio can be formed.

**Lesson worth keeping:** the generated-table scripts were written to showcase the headline quantity and silently narrowed what the CSV already knew. Both CSVs were complete; only the LaTeX projections were lossy.

---

### BUG-16 — LAT appendix photon-energy header changed MeV to keV — **REVERTED, the original was correct**

**This entry was wrong and its "fix" has been undone.** See BUG-18 for the correct diagnosis.

The claim rested on an arithmetic slip: 27428.80 MeV is **27.43 GeV**, not 27 TeV (27428.80 *GeV* would be 27 TeV). The units really are MeV, so the original `(MeV)` header was right and was restored on 2026-08-21.

**Original diagnosis:**

`GRBResearchPaper/appendices/appendix_LAT_info.tex:20` headed the photon-energy column "(MeV)". The tabulated values include `27428.80`, which in MeV would be 27 TeV. Everything else in the project treats them as keV — the paper text and `lorentz_factor.py` both read 27428.80 keV = 27.43 GeV — so the header was simply wrong. Corrected to (keV).

---

### OBS-08 — per-episode LAT photon data is hand-transcribed from the paper table — **FIXED**

**Resolved 2026-08-21.** Both consumers now read one generated CSV, so the table and the code cannot disagree.

- `LAT_analysis/txt_to_csv.py` collects every episode's `*_analysis_result_*.txt` and `*_fit_results_*.txt` into `LAT_analysis/lat_photons.csv` (26 rows). Episode labels are resolved by matching interval bounds against `results.json` via `grb_research`'s own parser, never from directory names — `Ep5A` is `EX1`, not `EX4`.
- `LAT_analysis/csv_to_latex.py` renders the appendix table from that CSV, with an `% AUTO-GENERATED` header. At the time of this fix, `appendices/appendix_LAT_info.tex` was a two-line `\input` of `tex_files/generated/lat_info_table`; that wrapper file was removed 2026-09-06 and `main.tex` now `\input`s the generated table directly.
- `lorentz_factor.py`'s `LAT_PHOTONS` dict is gone, and with it the hardcoded `LOW_SIGNIFICANCE` set — weak detections are now derived from `TS < 25`, which reproduced the old list exactly.

**Verified.** All 26 episodes were cross-checked against the source twice (by value, then by interval-bound mapping): zero mismatches, so the transcription had been faithful. Switching to the CSV moved $\Gamma_\text{min}$ by at most $2.9\times10^{-6}$ relative and changed no rounded value — the CSV carries full source precision (`301.204`) where the hand copy carried the appendix's rounded `301.20`. The generated table reproduces the hand-maintained one except for four deliberate consistency fixes: `40.503`→`40.502` (the same photon was written two ways in two rows), `753.113`/`722.523`→2 dp to match every other row, and `1.28`→`1.280`.

**One input remains hand-supplied:** the `\pFlux` footnote upper limits are 95% profile-likelihood limits from gtlike's `UpperLimits`, whose output is not in `LAT_analysis/`. They sit 1–8% above (Flux + 2σ), so they cannot be reconstructed from the fit results either. They live in `FLUX_UPPER_LIMITS` in `csv_to_latex.py`, and the generator **raises** if that dict disagrees with the derived `TS < 25` set — a crash, not a silent wrong number. Full write-up in `LAT_analysis/LAT_analysis.md`.

**Original diagnosis:**

`lorentz_factor.py`'s `LAT_PHOTONS` dictionary was transcribed by hand from `appendices/appendix_LAT_info.tex` (`tab:burst_table`), because the LAT photon energies, arrival times and TS values exist **only** inside that LaTeX table — there is no machine-readable source in `GRBResearchWork`.

This is exactly the failure mode of BUG-10 (the hand-typed $E_\text{iso}$ table that drifted from the pipeline and matched no code state). The same risk applies here in reverse: edit the LaTeX table and the code silently keeps the old numbers.

**Suggested fix:** put the LAT photon data in a CSV in `GRBResearchWork`, have `lorentz_factor.py` read it, and generate the appendix table from it with a `csv_to_latex.py` — the convention already used for the Amati and photospheric tables. Then the table and the code cannot disagree.

Until then, the transcription is a known single point of drift. It was checked once, on 2026-08-21, against the table as it then stood.

---

### OBS-07 — $\Gamma_\text{min}$ is computed only for T90, blocking a clean comparison with the thermal $\Gamma$ — **FIXED**

**Done.** `lorentz_factor.py` now evaluates every episode with LAT coverage, with MC errors, using each episode's duration as $t_v$ uniformly. **The outcome overturned the earlier suggestion of a $Y$ constraint:** on matching boundaries the opacity floors (78–137) are cleared by the thermal $\Gamma$ (752–861) by factors of 5–10, but are far too weak to bound $Y$ at all — EX0 requires only $Y \geq 7\times10^{-5}$ against $Y \geq 1$ by definition. The apparent tension was entirely an artefact of comparing against Abdo et al.'s differently-binned bound.

**Original diagnosis:**

`GRBResearchWork/codes-for-paper/lorentz_factor/lorentz_factor.py` evaluates the $\gamma\gamma$-opacity bound for the T90 interval alone. Because $\Gamma_\text{min}$ scales with the photon flux $f_1$, a burst-averaged spectrum yields a systematically **weak** bound: ours is $\Gamma_\text{min} = 258$ for GRB080916C, whereas \citet{Abdo2009FermiObservations080916C} obtain $608 \pm 15$ and $887 \pm 21$ from time-resolved bins d and b of the same burst.

This matters because $\Gamma_\text{min}$ is the only independent check on the Phase 2 thermal $\Gamma$, and the two cannot currently be compared on equal footing:

- Abdo et al.'s bin b spans **3.6–7.7 s**.
- Our BB-bearing intervals are EX0 ($-0.128$–4.864 s) and TR1 (1.280–4.864 s) — only $\sim 1.3$ s of overlap.
- Our TR2 (4.864–15.040 s), which covers most of bin b, has **no BB component** (BEST model is plain `BAND`), so no thermal $\Gamma$ exists there.

Since $\Gamma$ varies across episodes (752 → 852 → 861 in our results), comparing across mismatched intervals is not meaningful. Any constraint on $Y$ derived this way — e.g. the naive $Y \gtrsim 1.9$ from requiring $\Gamma \geq 887$ — is indicative at best and **must not be quoted as a result**.

**Fix:** compute $\Gamma_\text{min}$ per episode, on the same boundaries as the thermal analysis. That makes the comparison internal and self-consistent and would turn the $Y$ constraint into a defensible number. Requires per-episode $E_\text{max}$ and $t_\text{v}$, both already available in the LAT appendix table.

---

### BUG-15 — Phase 2 figure did not distinguish episodes; legend handles lacked markers — **FIXED**

Two defects in `pe_er_photosphere.py::make_plot`, both caught by the user:

1. **All episodes of a burst were drawn identically.** GRB131014A's five episodes appeared as five indistinguishable orange curves under a single legend entry, and GRB080916C's three points were identical circles — violating the per-episode-marker convention in `CLAUDE.md`. Fixed by giving each episode its own line style plus its `EpisodeMarkerResolver` marker, and labelling each series with its burst, episode *and* model name.
2. **The marker never reached the legend.** The curve (which carried the label) and the fiducial-redshift marker were drawn as two separate artists, so the legend handle was a bare line — precisely the marker that identifies TR1 from TR2 from EX0 was missing for every swept burst. Fixed by drawing one artist with `markevery=[fiducial]`.

Also fixed in the same pass: the inner `label` variable shadowed the axis-label loop variable, so both y-axes were captioned "GRB231129C TR1 (SBPL_BB)" instead of $r_0$ and $\Gamma$; and the 13-entry legend was moved outside the panels into a single left-hand column, since in-axes placement hid the GRB140206B track completely.

---

### OBS-01 — inconsistent package import style — **FIXED**

**Resolved 2026-08-21, at the user's request** ("the switch to `grb_research` should be made"). Standardized every `from src.grb_research import ...` to `from grb_research import ...`, matching the already-dominant convention (`amati_relationship.py`, `bb_flux_fraction.py`, `pe_er_photosphere.py`, `lorentz_factor.py`, and others).

**Files changed:** `codes-for-paper/butterfly_plots/butterfly_all.py` (also had `from src.grb_research.grb_calculations`/`grb_constants`, both collapsed to the `grb_research` top-level namespace), `codes-for-paper/best_names.py`, `codes-for-paper/amati_relationship/amati_helpers.py`, `codes-for-paper/fluence/grb_fluence.py`, `codes-for-paper/model_parameters/utils.py`, `codes-for-paper/evolution_of_kt/kt_evolution_080916c.py` (the known-broken example per `CLAUDE.md` — import fixed anyway for consistency, since it doesn't affect its other defects).

**Verified without re-running the heavy MC pipelines:** every symbol each file imports was checked present on the `grb_research` public API (`hasattr`) under `PYTHONPATH=src` alone — no longer needs `PYTHONPATH=.:src`. Not run end-to-end, since these scripts take minutes each; `pip install -e .` was not done (out of scope for this fix, still an option).

**Original diagnosis:**

`codes-for-paper/amati_relationship/amati_helpers.py:8` imports `from src.grb_research import ...`, while `amati_relationship.py:16` uses `from grb_research import ...`. The two work under different working directories / `PYTHONPATH` settings. Note that `grb_research` is **not installed** into `.venv` — scripts run only with `PYTHONPATH=src` (or from an IDE with `src` marked as a source root). Worth settling on one convention, and ideally `pip install -e .` into the venv, but not touched without your say-so.

---

### BUG-18 — LAT photon energies are MeV but $\Gamma_\text{min}$ treated them as keV — **FIXED**

Found while checking whether `LAT_analysis/` could replace the hand-transcribed `LAT_PHOTONS` (OBS-08). It supersedes BUG-16, which "fixed" the same confusion in the wrong direction.

`compute_gamma_min` computes

```python
gamma_min = tau_hat**e1 * (E_max_keV / 511.0) ** e2 * (1 + z) ** e3
```

where 511 is $m_ec^2$ in **keV**, so the argument must be keV. The values passed are in **MeV**, unconverted — a factor of 1000 low.

**The units are MeV.** Three independent confirmations:

- The source field is literally named `P > 0.9 Max (E) MeV` in every `LAT_analysis/*/Ep*/*_analysis_result_*.txt`.
- The smallest values in the whole dataset — 117.7, 122.9, 194.5 — sit just above the LAT selection floor (`e_min: 100`). Read as keV they would be ~0.12 MeV, three orders of magnitude below the LAT band. The user confirmed the LAT analysis is thresholded at >100 MeV.
- `fit_results` reports `Energy Flux (0.1 - 100.0) GeV`, and §5 of the paper speaks of "the detection of GeV photons". 27428.80 keV would be 27.4 MeV, not a GeV photon.

**Impact.** $\Gamma_\text{min} \propto E_\text{max}^{(\alpha-1)/(2\alpha+2)}$, so the error propagates as $1000^{(\alpha-1)/(2\alpha+2)}$ = 3.4–4.5:

| GRB080916C | published | corrected |
|---|---|---|
| T90 | 134 | 507 |
| EX0 | 78 | 350 |
| TR1 | 86 | 380 |
| TR2 | 129 | 447 |
| TR3 | 137 | 504 |
| TR4 | 80 | 295 |
| EX1 | 83 | 287 |
| TR5 | 71 | 260 |

**Fix.** `lorentz_factor.py` is now consistently MeV — matching `tau_hat`, whose `0.511` is already MeV, and `f_1` in ph/cm²/s/MeV. `compute_gamma_min` takes `E_max_MeV` and divides by 0.511; the CSV column is `E_max_MeV`. The LaTeX table divides by $10^3$ to print GeV, which was correct before and is unchanged.

**Verified:** regenerated $\Gamma_\text{min}$ matches the analytic prediction $\Gamma_\text{old}\times1000^{(\alpha-1)/(2\alpha+2)}$ to $2\times10^{-16}$, and every other CSV column ($f_1$, $\beta$, $t_v$, $t_\text{arr}$, $z$) is bit-identical.

**Knock-on: the "conservative by construction" narrative shrinks.** The thermal $\Gamma$ (752–861) exceeded the opacity floor by "factors of five to ten"; it is now **1.5–2.4**. The two independent methods therefore agree far more closely than the paper claimed — a stronger consistency check, not a weaker one. Updated in `abstract.tex`, `section-5-data-analysis.tex` (×2) and `section-7-conclusion.tex` (×2).

The gap to \citet{Abdo2009FermiObservations080916C} also narrows: our TR3 was 137 against their 887, attributed in HANDOFF §5 to Lithwick & Sari Limit A being a simplified analytic form. At 504 the gap is a factor of 1.8, not 7 — most of it was this bug, not method.

---

### OBS-09 — `kt_bb_keV` means different things in the two Phase 1/2 CSVs — **NOTED, deliberate but undocumented until now**

Found during the Phase 4 spot-check of the `bb_fraction`/`photospheric_radius` tables against `results.json`, run 2026-08-21 via four parallel per-GRB subagents (one per burst in the sample), each independently pulling raw fit parameters from `results.json`, checking transcription into `bb_flux_fraction.csv`/`pe_er_photosphere.csv` and the two generated LaTeX tables, and hand-recomputing $f_\text{BB}$ and the Pe'er (2007) eqs. (1)/(4)/(5) outputs. All four agents flagged the same thing independently, and it was confirmed by direct inspection (not just relayed): `bb_flux_fraction.csv`'s `kt_bb_keV` column is the **raw best-fit** value from `results.json` (bit-identical, e.g. GRB080916C T90: `44.452754974365234`), while `pe_er_photosphere.csv`'s same-named column is the **Monte Carlo median** of resampled `kt_bb` draws (`44.46361948397558` for the same interval — a ~0.02% difference, consistent with ordinary MC noise, not an error). By contrast `f_bb` genuinely is numerically identical between the two CSVs, which **is** already documented (`photospheric_radius.md` §2.4: same seed, same draws).

**Verdict: not a bug.** Both values trace back correctly to the same fit; nothing downstream is computed from the wrong one. But the same column name carrying two different statistics across sibling CSVs is exactly the kind of thing that misleads a future reader diffing the files, so it did not deserve to stay silent. Documented in `bb_fraction.md` and `photospheric_radius.md`.

**No table or CSV data changed.** All other spot-check findings across all four bursts (13 BB-inclusive intervals total) were clean: exact transcription of `model_name`/`kt_bb`/errors from `results.json`, correct value+uncertainty rendering in both generated LaTeX tables (including the fiducial-$z$ dagger for the three bursts without a measured redshift), and independent hand-recomputation of $f_\text{BB}$ and eqs. (1)/(4)/(5) agreeing with the CSVs to within ≤3% (GRB080916C, GRB140206B) or ≤0.6% (GRB131014A, GRB231129C) — well inside the slop expected between an MC median and a point-estimate/central-value evaluation. Closes Phase 4 item 3 in `HANDOFF.md` §6.

---

### OBS-10 — `bb_flux_fraction.py`'s RNG bypassed `get_rng` despite `bb_fraction.md` claiming otherwise — **FIXED**

Found 2026-08-31 while extending `bb_flux_fraction.py` for the rest-frame $f_\text{BB}$ work (weakness-review #2, see `HANDOFF.md` §12). `bb_fraction.md` §2.4 stated "RNG is obtained via `get_rng(seed=..., rng=...)` and threaded through, never a global `np.random.seed()`" — true of the *intent*, but `compute_fraction`'s actual code (line 143, pre-fix) was `rng if rng is not None else np.random.default_rng(seed)`, a hand-rolled equivalent that bypassed the shared helper entirely.

**Verdict: not a correctness bug.** The hand-rolled fallback is functionally identical to `get_rng`'s behavior for the cases actually exercised (a `seed` with no `rng`, or an explicit `rng`), so no published number is affected. It is a documentation/implementation mismatch — the kind of thing that misleads a future reader who greps for `get_rng` usage to confirm the RNG convention is followed everywhere.

**Fix.** Removed the local RNG construction entirely; `seed`/`rng` are now passed straight through to `draw_model_samples`, which already calls `get_rng` internally (`grb_calculations.py:921`) — one correct call site instead of a redundant second one. Documented in `bb_fraction.md` §2.4.

---

### BUG-19 — `convert_sbpl_to_band` silently overwrote a caller-supplied `rng` whenever `seed` was also given — **FIXED**

Found 2026-09-03 while auditing every function that accepts both `seed` and `rng` parameters, as part of the RNG-seeding overhaul (`SEED_PLAN.md`). `convert_sbpl_to_band`'s old body opened with `if seed is not None: rng = np.random.default_rng(seed)` — unconditional, with no check for whether `rng` had already been supplied non-`None`. Any caller passing both `seed=` and `rng=` together would silently get the `seed`-derived generator instead of the one they explicitly passed in — the opposite of `get_rng`'s own precedence, where `rng` wins.

**Impact: none realized.** No caller in the codebase ever passed both arguments together — confirmed by checking every call site before the fix (`model_parameters/utils.py::extract_kt_epeak_from_models`, the only caller, always passed `rng=` alone). A latent precedence bug, not a wrong number in the paper.

**Fix.** Routed through `get_rng(seed=seed, rng=rng)` like every other RNG-consuming function in the project; later in the same pass `seed` was removed from the signature entirely and `rng` made required, once every caller was confirmed to always supply it. Full seeding-scheme writeup in `src/grb_research/SEEDING.md`.

---

### BUG-20 — `FluxFluenceCalculator.__init__`'s falsy `if seed:` check silently dropped `seed=0` — **FIXED**

Found 2026-09-03, same audit as BUG-19. The old body was `if seed: rng = np.random.default_rng(seed)` — a truthiness check, not an `is not None` check — so a caller passing `seed=0` (a legitimate integer) would skip the branch entirely and leave `self.rng` at its default of `None`, crashing on the first subsequent `self.rng.multivariate_normal(...)` call with an `AttributeError`, not silently producing a wrong number.

**Impact: none realized.** No caller in the codebase has ever passed `seed=0` — caught by inspection, not by a crash in the wild.

**Fix.** Routed through `get_rng(seed=seed, rng=rng)`, whose own `if seed is None and rng is None` check is correct against `None` rather than falsiness. `seed` was removed from the signature entirely in the same pass as BUG-19.

### BUG-21 — `codes-for-paper/variability_analysis/fitter*.py` import a `variability_timescale` module that no longer exists — **FIXED**

**Fix.** `codes-for-paper/variability_analysis/` now carries its own local `norris_fit.py` and `light_curves.py` (the former a proper rename of the old `norris..py`, not just a copy), and every script that previously imported from the dead `variability_timescale` package now imports these local modules instead: `fitter.py`, `fitter_GRB140206275.py`, `fitter_GRB140206275_simple.py`, and `fitter_GRB231129779.py` were repointed in commit `ce18e1c` (`[main-minor-89]`, 2026-09-16); `fitter_CLAUDE_GRB131014215.py` and `experiments/window_sensitivity_GRB231129779/window_sensitivity.py` (which also got its own local `norris_fit.py` copy) were fixed in the working tree shortly after. **Confirmed, not assumed:** all six scripts' `from light_curves import lightcurve_data` / `from norris_fit import NorrisFitter, ...` imports resolve cleanly (`PYTHONPATH=src` from `GRBResearchWork/`, `.venv` active) — no `ModuleNotFoundError`, all expected symbols present.

**Original diagnosis:** Found 2026-09-16 while auditing project `.md` files for staleness. `GRBResearchWork/codes-for-paper/variability_analysis/fitter.py`, `fitter_CLAUDE_GRB131014215.py`, `fitter_GRB140206275.py`, `fitter_GRB140206275_simple.py`, `fitter_GRB231129779.py`, and `experiments/window_sensitivity_GRB231129779/window_sensitivity.py` all did `from variability_timescale.norris_fit import NorrisFitter` and/or `from variability_timescale.light_curves import lightcurve_data`. But `GRBResearchWork/variability_timescale/norris_fit.py` and `light_curves.py` no longer existed on disk — only stale `__pycache__/*.pyc` remained, and the whole `variability_timescale/` directory was untracked (`??` in `git status`), so there was no commit to recover the missing source from. Confirmed via `ModuleNotFoundError: No module named 'variability_timescale.norris_fit'`.

**Root cause.** The automated Norris-fitting pipeline these scripts depended on (MEPSA/`scipy.signal.find_peaks` peak detection + a per-window local fit, living in `variability_timescale/`) was abandoned after stalling on GRB080916C's TR2 (zero MEPSA detections — TR2 is one broad, smoothly-declining pulse with no local excess for a spike detector to find), in favor of a manual, GRB-by-GRB joint-Norris-fit track built directly in `codes-for-paper/variability_analysis/` (see that folder's `variability_analysis.md`). The older scripts' imports were never repointed when the automated pipeline's source files were deleted. Cross-referenced in `GRBResearchWork/PHASE5_TV_PLAN.md`'s 2026-09-16 status update.

---

### BUG-22 — `light_curves.py`'s `lightcurve_data()` combined per-channel errors linearly instead of in quadrature — **FIXED**

Found 2026-09-17 while cross-checking GRB080916C's 8-pulse Norris fit in CERN ROOT (`codes-for-paper/variability_analysis/experiments/root_fit_GRB080916009/GRB080916009_norris_fit.C`), exported from `lightcurve_data(..., errors=True)`. An errors-weighted `TGraphErrors` fit came back `chi2/ndf = 33.4/1077 ≈ 0.031` — far below the ~1 a correctly-weighted fit should give, meaning the error bars fed to ROOT were too large relative to the data's real point-to-point scatter. An unweighted fit (`TGraph`, ignoring the error column) tracked the data just as well, with an implied RMS residual of `sqrt(chi2/n) ≈ 252` cts/s — far smaller than what the error column itself claimed (~830–1000 cts/s per point).

**Root cause, confirmed by direct inspection.** `lightcurve_data()`'s column-indexing is correct — verified against `light_curves/GRB080916009/n3.dat`'s own header line (`; Tstart(s), Tend(s), RateCh1(ct/s), ErrRateCh1(ct/s), BgdCh1(ct/s), ErrBgdCh1(ct/s), RateCh2(ct/s), ...`), which exactly matches the function's `4*x + 2/3/4/5` offset formula (2 leading time columns, then repeating 4-tuples per energy channel). But summing the *rate* linearly across channels is correct (physical count rates from independent channels add), while summing the *rate error*/`background error` the same way (`data_rate_error += detector_data[:, 4*x+3]`) is not — independent statistical uncertainties combine in quadrature. For a 10–400 keV band this sums **89 energy channels** (`start_idx=7` to `end_idx=96`, counted directly against the file's own channel-boundary table), so the linear sum overstates the combined error by close to $\sqrt{89}\approx9.4\times$ relative to the correct combination — the right order of magnitude to explain both the ~830–1000 cts/s vs. ~252 cts/s discrepancy above and the `chi2/ndf ≈ 0.031 ≈ (1/9.4)^2` anomaly.

**Impact: none realized on anything published.** Checked every call site of `lightcurve_data()` in the repo before fixing: `light_curves/make_lightcurve.py` (all 3 calls) and every `fitter*.py`/`window_sensitivity.py` script in `codes-for-paper/variability_analysis/` (12 call sites) all use `errors=False`, the default. `errors=True` was never actually invoked anywhere in the codebase except the ad-hoc terminal commands that built the ROOT CSV export this session — never saved to a script, never part of any pipeline.

**Fix.** Accumulate `detector_data[:, 4*x+3] ** 2` / `detector_data[:, 4*x+5] ** 2` inside the per-channel loop instead of the bare values, then `np.sqrt(...)` once after the loop, before returning. Applied identically to all six copies of this function in the repo (`light_curves/light_curves.py`, `codes-for-paper/variability_analysis/light_curves.py`, and its four `experiments/*/` copies), confirmed still byte-identical to each other (function body) and all six compile cleanly. Re-ran the ROOT weighted fit after the fix: `chi2/ndf = 1565.5/1077 ≈ 1.454` — a physically credible value, consistent with the unweighted fit's own RMS-residual estimate. Documented in `codes-for-paper/variability_analysis/variability_analysis.md`'s GRB080916C section.

---

### OBS-11 — GRB231129C's zenith-cut boundary (~100.2–100.3°) exceeds the stated 100° cut in §4.2.4, exposure-loss correction unverified — **DEFERRED, user's call**

**Where:** `tex_files/section-3-data-preparation-and-analysis.tex` §4.2.4 (GRB231129C LAT data reduction), `GRBResearchPaper`.

§4.2.4 states the 12° ROI's outer boundary reaches ~100.2–100.3° for GRB231129C, exceeding the paper's stated 100° zenith cut. The text acknowledges the overshoot but does not say what was done about it — whether the resulting exposure loss was corrected in the `gtmktime`/`gtselect` reduction.

**Checked, not assumed.** The actual `gtmktime`/`gtselect` commands for this burst's LAT reduction were run outside this repo and never checked in, so there is no source to verify what was actually done. **Confirmed directly with the user: it was not corrected.**

**Decision (user's, 2026-09-22).** Leave the paper text and data as-is. Properly addressing it means re-running the LAT data reduction with a tightened zenith cut — a data-regeneration task, not a text fix, out of proportion to a quick-fixes pass. Revisit if/when LAT data regeneration is otherwise in scope.

**Originally logged:** `quick-fixes-mid-priority.md`, 2026-09-04 (moved here 2026-09-22 once that file's other items were confirmed resolved and the file removed from `GRBResearchWork/md_file_refs/`).

---

### BUG-23 — GRB080916C's `fitter.py`/`window_sensitivity.py` silently picked `n4` instead of the documented `n3` NaI detector after a data resync — **FIXED for GRB080916C; confirmed present, unfixed, for GRB131014A and GRB140206B**

Found 2026-09-22, mid-session, while rerunning `codes-for-paper/variability_analysis/experiments/window_sensitivity_GRB080916009/window_sensitivity.py` to re-verify the already-committed 6-vs-7-pulse window-sensitivity check (`window_sensitivity_results.csv`, committed at `8e98c98`) after `light_curves/GRB080916009/`'s NaI `.dat` files were freshly synced onto this machine. The rerun's converged pulse parameters didn't match the committed CSV (e.g. pulse 1 `t_s` off by ~0.3s, well outside MC/optimizer noise) despite `n3.dat`/`n4.dat`/`window_sensitivity.py` all being byte-identical to `HEAD`.

**Root cause.** Both `fitter.py` and `window_sensitivity.py` build the single-detector NaI file list via `os.listdir(LC_DIR)` and then unconditionally take element `[0]` — `variability_analysis.md`'s own "Open cross-burst issue" note already assumed this always lands on the alphabetically-first detector ("`n3` for GRB080916C"), but `os.listdir()` order is filesystem directory-entry order, not alphabetical, and isn't guaranteed stable across a resync. Confirmed directly: `os.listdir('light_curves/GRB080916009')` returned `['n4.dat', 'n3.dat', ...]` post-resync, so `dat_NaI[0]` silently became `n4` — a genuinely different detector (peak count rate 1614.97 vs `n3`'s 1879.66 cts/s, `np.array_equal` confirms the two light curves are not the same data).

**Impact.** Every `fitter.py`/`window_sensitivity.py` output generated on this machine between the resync and this fix used `n4`, not `n3` — this includes a `norris_fit_results_GRB080916009.csv` regeneration and a `.norris_fitted_GRB080916009.png/.pdf` regeneration (with a new LAT-photon overlay) done earlier in the same session, both silently on the wrong detector. Caught before being reported as final, not after.

**Fix for GRB080916C.** `dat_NaI = sorted(...)` in both `fitter.py` and this experiment's `window_sensitivity.py`, pinning detector selection to alphabetical order regardless of filesystem listing order. Re-running both after the fix reproduces the committed `window_sensitivity_results.csv` to ~4 significant figures (small residual differences consistent with ordinary optimizer floating-point noise, not a data difference) and gives `y_max_cts_per_s = 1879.6628` in the regenerated `norris_fit_results_GRB080916009.csv`, matching `n3`'s own peak rate exactly.

**Extended to the other three bursts, 2026-09-22 — verified, not yet fixed.** Checked whether GRB131014A/GRB140206B/GRB231129C's `fitter*.py`/`window_sensitivity.py` scripts share the same unsorted `os.listdir()[0]` pattern (they do, confirmed by inspection of all 6 remaining scripts) and whether it's currently live-broken, by comparing each burst's *current* `os.listdir()[0]`-picked detector against the detector implied by its *committed* CSV's `y_max_cts_per_s` (cross-checked against each candidate `.dat` file's own computed peak count rate):

| GRB | documented/committed detector | current `os.listdir()[0]` | live-broken now? |
|---|---|---|---|
| GRB131014A | `na` (peak 44084.39) | `nb` (peak 43019.09) | **yes** |
| GRB140206B | `n3` (peak 4271.27) | `n0` (peak 4672.74) | **yes** |
| GRB231129C | `n7` (peak 10487.75) | `n7` (peak 10487.75) | no, currently correct |

**The `sorted()` fix used for GRB080916C does not generalize.** It happened to work there only because `n3 < n4` alphabetically among that burst's 2-detector set. GRB131014A has 3 detectors and the documented one (`na`) is not alphabetically first (`sorted(['n9','na','nb'])[0] == 'n9'`); GRB140206B likewise (`sorted(['n0','n1','n3'])[0] == 'n0'`, not `n3`); applying `sorted()` to GRB231129C would actively break it (`sorted(['n3','n6','n7'])[0] == 'n3'`, not the correct `n7`). The right fix is an explicit hardcoded detector name per burst, not an alphabetical sort — GRB080916C's own `sorted()` fix should be revisited to an explicit pin too, since it is currently right by coincidence, not by design.

**Not yet fixed — 6 files affected:** `fitter_GRB140206275.py`, `fitter_GRB140206275_simple.py`, `fitter_GRB231129779.py`, and the three matching `experiments/window_sensitivity_GRB*/window_sensitivity.py` copies. (`fitter_CLAUDE_GRB131014215.py` and GRB131014A's own experiment-folder copy have the same pattern too, covered by the GRB131014A row above.) User asked to document this in the `.md` files first, before any fix is applied.
