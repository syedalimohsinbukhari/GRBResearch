# Amati relation ($E_{i,p}$–$E_\text{iso}$) — method notes

Companion to `amati_relationship.py`, `amati_helpers.py`, `csv_to_latex.py` and `csv_to_latex_unknown_z.py`.

**Written 2026-10-20 (documentation sweep), retrospectively.** This folder had no method note although the project convention requires one per folder. Everything below was reconstructed from the code, the two output CSVs and the git log, not written as the work happened, so the "who decided" attributions are limited to what the code comments and commit messages record; anything else is marked as unattributed. Treat the validation section as the checks that can be evidenced today, not as a full history.

---

## 1. What the code computes

For each interval of the paper's best-fit model set (PL models removed, `grb_best = [ModelSet([i for i in j if i.name != "PL"]) ...]`, since a pure power law has no peak), one point on the Amati plane with asymmetric Monte-Carlo errors:

- **Rest-frame peak energy** $E_{i,p} = (1+z)\,E_{p,\text{obs}}$. $E_{p,\text{obs}}$ is the fitted `e_peak` for BAND/CPL-based models, and for SBPL models it is derived from the break energy and the two indices with `break_e_to_e_peak` (`amati_helpers._compute_ep_eiso`).
- **Isotropic-equivalent energy**, `mc_e_iso_sampler(..., method=2)` in `grb_calculations.py`:

$$E_\text{iso} = \frac{4\pi\,d_L^2(z)\,S_\text{bol}}{1+z},\qquad S_\text{bol} = \Delta t \int_{E_1/(1+z)}^{E_2/(1+z)} E\,N(E)\,\mathrm{d}E$$

  with $\Delta t$ the interval duration, the rest-frame band $E_1$–$E_2$ = **1 keV – 10 MeV** (`bol_min=0.0`, `bol_max=4.0` in log10 keV, shifted to observed energies by $\log_{10}(1+z)$), $d_L$ from `FlatLambdaCDM(H0=69.6, Om0=0.286)` (Fana Dirirsa 2019, the project cosmology; defaults hard-wired in `mc_e_iso_sampler`). `method=2` takes the model's own bolometric flux directly as the observed fluence, i.e. the k-correction is not applied separately (`method=1` would multiply a detector-band flux by a bolometric/detector ratio, which cancels by construction for a model-anchored flux; see `BUGS.md` OBS-05). The $(1+z)$ in the denominator is the time-dilation factor.
- **Reference relation** drawn on every panel (`amati_relationship_dirirsa2019`), Fana Dirirsa et al. (2019):

$$\log_{10}\frac{E_\text{iso}}{10^{52}\,\text{erg}} = k + m\,\log_{10}\frac{E_{i,p}}{950\,\text{keV}},\qquad \sigma_y^2 = \sigma_k^2 + x^2\sigma_m^2 + \sigma_\text{ext}^2$$

  with $k=1.67\pm0.16$, $m=1.16\pm0.37$, $\sigma_\text{ext}=0.47$, bands at 1, 2, 3 $\sigma$. The code takes the **mean** of $\sigma_y$ over the plotted range (`use_average=True`), i.e. the bands are parallel to the line, not funnel-shaped. These constants match what the paper now states after the 2026-09-06 fix (paper commit `44fe7b8`, which replaced an earlier $E_{0,p}=300$ keV, $m=0.5$, $k=0$ in the prose that contradicted both the figure and the cited source).
- **Two redshift regimes.** GRB080916C is evaluated at its measured $z=4.35$. The other three bursts (no measured $z$) are evaluated at assumed $z = 1, 2, 3, 5, 7$ (`Z_VALUES_UNKNOWN`), drawn as a dashed locus on the plot. The z grid was unified with the other sweeps (Lorentz factor, $f_\text{BB}$, photospheric radius) on 2026-10-02.

Outputs: `amati_relationship.csv` (8 rows, GRB080916C only), `amati_relationship_unknown_z.csv` (17 rows: GRB131014A 5, GRB140206B 7, GRB231129C 5), `amati_relationship.png`/`.pdf` (the four-panel figure), and the generated tables `amati_relationship_table.tex` (`tab:eiso`) and `amati_relationship_unknown_z_table.tex` (`tab:eiso-unknown-z`).

## 2. Judgement calls

- **Bolometric rest-frame band 1 keV – 10 MeV** — conventions as in `bb_fraction.md` §2.2b; attribution not recorded in the code. The same shift formula is shared with `bb_flux_fraction.py`'s rest-frame $f_\text{BB}$ by design.
- **Mean-$\sigma$ confidence bands** (`use_average=True`) — set in `amati_relationship.py` (`amati_kw`); the reason is not recorded.
- **Unknown-$z$ table design — user review, 2026-09-30** (from `csv_to_latex_unknown_z.py`'s docstring and the code comments in `amati_relationship.py`):
  - $E_{i,p}$ is **not** tabulated as a function of $z$: it is an exact identity $(1+z)E_{p,\text{obs}}$ (one draw set per episode, rescaled, never refit), so the table carries a single $E_{p,\text{obs}}$ column and the caption states the scaling. $E_{p,\text{obs}}$ is recovered from the $z=1$ entry; the $z=7$ entry is a free consistency check that is printed to the console, not stored.
  - $E_\text{iso}$ gets one column per $z$ because it runs through $d_L^2(z)$ and a bolometric integral, so it is not a fixed multiple of anything.
  - Each $E_\text{iso}$ cell shows value and MC error. An earlier draft showed bare medians on the grounds that the MC error is 0.4–3.6% of the value against a 37–650% change between adjacent $z$ columns (checked, per the docstring); the user preferred value $\pm$ error to match the other tables. The bare-median draft is kept in `_backup_2026-09-30_bare_eiso/` (referenced in the docstring; **that folder is not present in the working tree on this machine**).
- **One shared RNG, one seed.** `SEED = seed_from_name(__file__)` = **3271127181** (`seed_registry.yaml`: `amati_relation`; confirmed equal to the `seed` column of both CSVs and to `seed_table.csv`), `n_samples = 10000` (the shared `N_SAMPLES`). Every model and every $z$ draws sequentially from the same generator, so draws are independent across rows and across $z$ values within an episode. This differs from the Lorentz-factor sweeps, which deliberately reuse one draw set across $z$.
- **Variable-name collision fixed** (comment in `amati_relationship.py`): the unknown-$z$ loop's unpacking targets carry a `uz_` prefix because top-level script code once let `model_list` collide and garbled the Model column of the GRB080916C rows in `amati_relationship.csv`. Recorded so it is not "simplified" back.
- **Episode ordering** in both CSVs follows `grb_research.episode_order` (T90, EX0, TR1…TRn, EX1), not `results.json` order.

## 3. Conventions that could plausibly have gone another way

- Observer vs rest frame for $E_p$: the CSV `E_i_peak` is rest-frame; the unknown-$z$ CSV stores observer-frame `E_p_obs` instead (see above).
- $E_\text{iso}$ is model-anchored (OBS-05): it inherits the fitted model's extrapolation outside the detector band.
- GRB080916C is the only burst plotted with a point rather than a locus; the others cannot be placed on the plane without assuming $z$.
- Column names bake units in as `E_i_peak__1e+03_keV`, `E_0_iso__1e+52_erg` (normalisation in the name: $10^3$ keV, $10^{52}$ erg), per the project CSV standard; `n_samples` and `seed` are constant columns.

## 4. Validation that can be evidenced

- **BUG-11 / OBS-03:** the $E_\text{iso}$ scale was corrected (date not recorded here; `BUGS.md` has it) (rest-frame grid previously paired with an observer-frame variable underestimated $E_\text{iso}$ by about 3.7x at $z=4.35$); the CSV was regenerated 2026-08-27 (commit `9ea9913`, "re-ran … with the corrected eiso calculations"). See `BUGS.md`.
- **Constants cross-check:** paper commit `44fe7b8` verified $E_{0,p}=950$ keV, $m=1.16\pm0.37$, $k=1.67\pm0.16$ against Fana Dirirsa (2019) Table 3 (their "F10" row) and against the values hard-coded in `amati_helpers.py`. `sigma_ext = 0.47` is in the code; the commit message quotes $0.47\pm0.12$ for the paper text.
- **Seed:** `amati_relationship.csv` and `amati_relationship_unknown_z.csv` both carry `seed = 3271127181` and `n_samples = 10000`, checked 2026-10-20 against `seed_table.csv`.
- **$E_{p,\text{obs}}$ consistency at $z=1$ vs $z=7$:** claimed in the code comments to agree to <0.3% on every row (MC noise only); the console output is not saved and was not re-run in this sweep, so this is the script's own statement, not re-verified today.
- **Not independently checked:** no literature comparison of the individual $E_\text{iso}$ values for GRB080916C (e.g. against the published $E_\text{iso}\sim 10^{55}$ erg ballpark of Abdo et al. 2009) is recorded here; the T90 row is $3.6\times10^{54}$ erg in the current CSV.

## 5. Known limitations and open questions

- No MC for the reference relation's own parameters: the bands are the analytic $\pm n\sigma$ of the published fit, not resampled.
- The unknown-$z$ points only show how the plane position moves with assumed $z$; they are not a test of the relation.
- The `_backup_2026-09-30_bare_eiso/` folder cited in `csv_to_latex_unknown_z.py` is absent here (possibly on the other machine).
- Cross-references: `BUGS.md` BUG-11, OBS-03, OBS-05; `bb_fraction.md` §2.2b for the shared rest-frame shift; `lorentz_factor.md` §13 for the parallel redshift sweep.

## 6. Files

| file | role |
|---|---|
| `amati_relationship.py` | computation, both CSVs, four-panel figure |
| `amati_helpers.py` | reference relation, $E_p$/$E_\text{iso}$ sampling, plotting helpers |
| `amati_relationship.csv` | GRB080916C, one row per episode |
| `amati_relationship_unknown_z.csv` | the three no-redshift bursts, $E_{p,\text{obs}}$ and $E_\text{iso}$ at $z=1,2,3,5,7$ |
| `csv_to_latex.py` → `amati_relationship_table.tex` | `tab:eiso`; also provides `decimals_for()` etc. used by the other generators |
| `csv_to_latex_unknown_z.py` → `amati_relationship_unknown_z_table.tex` | `tab:eiso-unknown-z` |
| `amati_relationship.png` / `.pdf` | figure, synced to `images/section4/` |
| `amati_relationship.md` | this file |
