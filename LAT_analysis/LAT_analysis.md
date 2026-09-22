# LAT analysis — method notes

Companion to `txt_to_csv.py` and `csv_to_latex.py`. Written 2026-08-21 while closing `BUGS.md` OBS-08.

Unlike the `codes-for-paper/<topic>/` folders, this one holds no physics. It is a **data-plumbing layer**: it turns the per-episode gtburst/gtlike output into one machine-readable table that both the paper's LAT appendix and the $\gamma\gamma$-opacity calculation read from.

## 1. What the code computes

Nothing derived — it *transcribes*. For each episode directory it reads two files gtburst/gtlike leave behind:

| File                      | Fields taken                                                                     |
|---------------------------|----------------------------------------------------------------------------------|
| `*_analysis_result_*.txt` | `# of Events`, `# of P > 0.9`, `P > 0.9 Max (E) MeV`, `Arrival Time (s)`, `TS`   |
| `*_fit_results_*.txt`     | `Index`, `Index Error`, `Flux (0.1 - 100.0) GeV`, `Flux Error (0.1 - 100.0) GeV` |

and writes `lat_photons.csv`, one row per (GRB, episode). `csv_to_latex.py` renders that CSV as the appendix table `tab:burst_table`.

## 2. The problem this solves

`lorentz_factor.py` previously carried a hand-transcribed `LAT_PHOTONS` dictionary, copied from `appendices/appendix_LAT_info.tex`, because the LaTeX table was the only place the photon energies existed. That is the exact failure mode of BUG-10, where the hand-typed $E_\text{iso}$ table drifted until it matched no code state at all — and here the risk ran in reverse: edit the appendix and the code would silently keep the stale numbers.

Both consumers now read `lat_photons.csv`, so the table and the code cannot disagree.

**Verified before switching over.** All 26 episodes were cross-checked, twice: once by matching each hardcoded `(E_max, t_arr)` pair back to a source directory, and once by mapping every directory to its episode label and comparing against the code. Zero mismatches — the transcription had been faithful. Switching to the CSV changed $\Gamma_\text{min}$ by at most $2.9\times10^{-6}$ relative, and no rounded value moved, because the CSV carries full source precision (`301.204`) where the hand copy carried the appendix's rounded value (`301.20`).

## 3. Judgement calls

- **Units are MeV, not keV** *(settled jointly; user confirmed the LAT threshold)*. The source field is named `P > 0.9 Max (E) MeV`, the smallest values in the sample (117.7, 122.9, 194.5) sit just above the 100 MeV LAT selection floor, and the user confirmed their LAT analysis starts at >100 MeV. Read as keV they would fall three orders of magnitude below the LAT band. This overturned BUG-16, which had "corrected" the appendix header from MeV to keV on the strength of an arithmetic slip. See BUG-18 — the consequence was a factor 3.4–4.5 error in $\Gamma_\text{min}$.
- **Episode labels come from `results.json`, never from directory names** *(proposed by Claude)*. `txt_to_csv.py` matches each directory's interval bounds against the interval strings in `results.json` and takes the label from `grb_research`'s own parser. Deriving labels from directory names is wrong: `Ep5A` is `EX1` for GRB080916C, not `EX4`. Bounds are compared with a $10^{-6}$ s tolerance because directory names strip trailing zeros (`1.28` vs `1.280`). See `CLAUDE.md`, "Episode naming".
- **The `TS < 25` footnote rule is derived, not listed** *(proposed by Claude)*. Both scripts apply the same `TS_SECURE_DETECTION = 25.0` threshold rather than carrying a hardcoded list of weak episodes. `lorentz_factor.py`'s former hardcoded `LOW_SIGNIFICANCE` set is gone; it matched the derived set exactly.
- **Flux upper limits stay hand-supplied** *(settled reluctantly; see §5)*.
- **Uniform decimal places** *(proposed by Claude)*. The hand-maintained table mixed 2 and 3 decimals on photon energy (`753.113` in the time-integrated block, `753.11` in the time-resolved one — the same photon) and mixed `1.28` with `1.280` on interval bounds. The generator uses 2 dp for energy and 3 dp for times and TS throughout. Note this is *not* the `decimals_for()` rule from the other `csv_to_latex.py` scripts, which sizes decimals to the smallest error in the column — these columns have no error bars, so there is nothing to size against.

## 4. Conventions that could have gone another way

- **Energy band** is the LAT selection band of the analysis itself, 100 MeV – 100 GeV (`e_min: 100`, `e_max: 100000` in the per-episode `LAT_config.yaml`), and fluxes are quoted over 0.1–100 GeV. Observer frame throughout; no redshift correction is applied anywhere in this layer.
- **"Photon energy" is the highest-energy event with association probability $P > 0.9$**, not the highest-energy event outright. This is what `\Gamma_{\min}$ needs — an insecure association would inflate the bound.
- **No Monte Carlo, no cosmology, hence no `n_samples`/`seed`/`H0`/`Om0` columns.** `CLAUDE.md`'s CSV standard requires those for MC-derived quantities; every number here is read straight from a file, so there is nothing to seed and no cosmology assumed. Units are still baked into the column names as the standard requires.

## 5. Known limitations

- **The `\pFlux` upper limits in the appendix footnotes are not derivable from this directory.** They are 95% profile-likelihood limits from gtlike's `UpperLimits`, whose output is not among the files here. They sit 1–8% above (`Flux` + 2$\sigma$), so they cannot be reconstructed from the fit results either. They live in `FLUX_UPPER_LIMITS` in `csv_to_latex.py`, keyed on `(GRB, episode)`, and the generator **raises** if that dict and the derived `TS < 25` set ever disagree — so the failure mode is a crash, not a wrong number. This is the one remaining hand-maintained input; adding the `UpperLimits` output to each episode directory would close it.
- **`GRB_DIRECTORIES` is an explicit map** rather than a pattern, because three directories are named for the Fermi trigger id (`018__GRB080916009`) and one for the paper burst (`GRB231129C`). Renaming that directory would let the map be derived.
- **The copy into `GRBResearchPaper/tex_files/generated/` is manual**, as for the Amati and photospheric tables — the two repos are independent (BUG-10).
- **The four `lat_analysis_<start>_<stop>.tex` files** this directory used to emit are superseded by `lat_info_table.tex`. They were also stale: each was missing its `Ep0` (time-integrated) row, because those directories were added after the last run.
