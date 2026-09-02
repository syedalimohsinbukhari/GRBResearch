# GBM-only refit robustness check (GRB131014A) — method notes

## 1. What the code computes

Weakness #6 in `review-resolution.md` (project root) asks whether GRB131014A's blackbody (BB) detections could be
an artefact of the joint GBM+LAT fit, since this burst sits at a large (~70 deg) LAT off-axis angle for its entire
T90. A WebSearch confirmed Fermi-LAT's effective area is ~0 below 100 MeV regardless of off-axis angle, which rules
out the specific mechanism the weakness doc proposed (LAT-band curvature masquerading as GBM-band BB curvature,
since the BB lives at kT ~ 27-59 keV, entirely inside GBM). That leaves a softer, indirect concern: LAT data could
still bias the *other* continuum parameters through the shared joint likelihood in a way that happens to favour a
spurious BB. The only way to close that off is to drop LAT and refit GBM alone.

The user ran that refit in RMFIT (GBM/NaI+BGO only, LAT excluded from the likelihood) for all five of GRB131014A's
episodes and added it to `results.json` under a new top-level key, `GRB131014215GBM`, alongside the existing joint
fit `GRB131014215`. This folder's script, `gbm_only_refit.py`, loads both and compares them episode by episode:

- Does the same model win in both fits (`model_name_joint` vs `model_name_gbm`)?
- Is the fitted kT_BB consistent within 1sigma between the two fits?
- Does BB clear the paper's own selection threshold (Delta-C-stat >= 28.74, Delta k = 2) in the GBM-only fit, using
  the same rule stated in `section-1-introduction.tex`?
- f_BB (observer-frame, bolometric, 1 keV - 10 MeV), recomputed for the GBM-only fit with the identical MC machinery
  `bb_fraction/bb_flux_fraction.py` uses, so the two columns are on the same footing.
- How much does dropping LAT actually shift kT_BB and f_BB, as a fractional (%) difference, not just whether the two
  are consistent within error — added 2026-09-01 per user request ("a cool addition" to see the size of the effect,
  not only whether it's statistically indistinguishable from zero).

Output: `gbm_only_refit_comparison.csv` (one row per shared episode) and three figures, ranked by how directly each
answers the #6 concern (see Sec 3 for the reasoning):

- **`gbm_only_refit_kt_comparison.png/.pdf` (primary)** — kT_BB, joint vs GBM-only, per episode.
- **`gbm_only_refit_delta_cstat.png/.pdf` (primary)** — Delta-C-stat (BB significance), joint vs GBM-only, against
  the 28.74 threshold. Added 2026-09-01 per user request, after discussing that this is the paper's own statistical
  currency for "is BB required," more central to #6 than kT_BB matching.
- **`gbm_only_refit_fractional_diff.png/.pdf` (secondary / optional)** — % difference (GBM-only vs Joint) in kT_BB
  and f_BB, two stacked panels. Useful supporting detail (quantifies *how much* the two fits differ, not just
  whether they're consistent), but not needed to make the #6 case on its own.

## 2. Code walkthrough (for reference, not just the paper)

This section walks through `gbm_only_refit.py` function by function, in plain terms, since the file mixes several
distinct ideas (results.json parsing, C-stat model selection, Monte Carlo flux decomposition) that are each simple
on their own but easy to lose track of together.

**Loading the two fits (`main()`, top half).** `results.json` is a nested dict: GRB name -> episode string ->
model name -> fit numbers. `GRB.from_dictionary()` (from `grb_research.grb_core`, not written in this folder) turns
one GRB's slice of that dict into a `GRB` object whose `.intervals` is a list of `TimeInterval`s (T90, EX0, TR1, ...),
each carrying a `.models` `ModelSet` (BAND, BAND_BB, SBPL_BB, ...) with each model's fitted parameters, its
covariance matrix, its C-stat/dof, and its `_status` flag (`BEST`/`SAFE`/`UNSAFE`, as recorded by RMFIT). This script
calls `GRB.from_dictionary()` twice — once for `"GRB131014215"` (the joint GBM+LAT fit) and once for
`"GRB131014215GBM"` (the GBM-only refit) — giving two independent `GRB` objects to compare episode by episode.
`episode_label()` just turns a `TimeInterval` into the short string used throughout the paper (`"T90"`, `"TR1"`, ...).

**Comparing model selection (`base_bb_delta_cstat`, `any_bb_step_clears_threshold`).** The paper's rule for adding a
blackbody is: fit BASE and BASE+BB, and keep the BB version only if the C-stat improves by at least 28.74 (this is
the "Delta k = 2" threshold from `section-1-introduction.tex`, since BB adds exactly 2 parameters, `amp_bb` and
`kt_bb`). `base_bb_delta_cstat()` looks up one base/BB pair (e.g. `"BAND"` and `"BAND_BB"`) in one episode's model
set and returns their C-stat difference (or `None` if either wasn't fit). `any_bb_step_clears_threshold()` calls
that for all four base models (BAND, CPL, SBPL, PL) and returns `True` if *any* of them clears 28.74 — this is
independent of whichever model `results.json` actually marks as `BEST`, which is what let it catch the T90 labelling
mistake (Sec 4): the recorded winner had no BB in it, but a different BB model in the same episode cleared the
threshold easily.

**The Monte Carlo pieces (`analytic_bb_bolometric_flux`, `split_energy_flux`, `compute_f_bb`,
`paired_ratio_pct`).** These exist because f_BB (what fraction of the total flux is thermal) is *not* a fit
parameter — it has to be computed from the fitted spectrum, and its uncertainty has to come from resampling the fit,
not from a single formula. The mechanism (copied from `bb_fraction/bb_flux_fraction.py`, per this project's
cross-folder-import convention):
1. `draw_model_samples()` (defined in `grb_research`, not here) draws N=10,000 parameter vectors from a model's fit
   covariance matrix (a multivariate Gaussian approximation to the fit's uncertainty), then rejects and re-draws any
   vector that violates a physical constraint (e.g. a negative amplitude) until all N are valid.
2. `split_energy_flux()` evaluates the spectral model at those N parameter vectors over an energy grid and integrates
   it, returning both the blackbody component's flux and the total flux, for all N draws at once.
3. `compute_f_bb()` divides the two to get N samples of f_BB, discards any that landed outside (0, 1) (numerical
   edge cases from the Gaussian tails), and reports the 50th percentile as the value with the 16th/84th as the
   asymmetric 1-sigma interval — the same summarise-by-percentiles pattern used everywhere MC is used in this
   project.
4. `paired_ratio_pct()` does the analogous thing for the *ratio* between the two fits: draw N samples from the joint
   fit's covariance and N samples from the GBM-only fit's covariance (independently — see Sec 4 for why this
   matters), compute kT_BB and f_BB for each, and take percentiles of `(gbm_sample / joint_sample - 1) * 100`
   directly, rather than computing each side's error separately and combining them with a formula. This is what
   produces the fractional-difference plot's (correctly asymmetric) error bars.

**Assembling one row (`build_row`)**: for one episode, this pulls each fit's `_status: BEST` model, calls the pieces
above, checks whether the two fits agree on the winning model and on kT_BB (within combined error), and packages
everything into one `ComparisonRow` — one row of the output CSV.

**The three figures (`make_plot`, `make_delta_cstat_plot`, `make_fractional_diff_plot`)** are matplotlib boilerplate
on top of the CSV: point plots with error bars, coloured by fit type and shaped by episode
(`EpisodeMarkerResolver`, an existing project convention — not invented in this folder), saved to `.png` and `.pdf`
via `update_style()`'s shared rcParams. `make_delta_cstat_plot()` plots `delta_cstat_joint`/`delta_cstat_gbm` on a
log y-axis (the values span roughly an order of magnitude) against a horizontal line at `DELTA_CSTAT_THRESHOLD`.

## 3. Result

**Which figure to use where (decided 2026-09-01, jointly):** `gbm_only_refit_kt_comparison.pdf` and
`gbm_only_refit_delta_cstat.pdf` are the **primary** pair for the paper -- between them they show that the same
model wins with matching kT_BB (physical consistency) *and* that BB clears the selection threshold by a wide margin
in both fits (statistical robustness), which is the whole #6 argument in two panels. `gbm_only_refit_fractional_diff.pdf`
is **secondary / optional** -- it quantifies the size of the (already-established-as-consistent) difference, useful
supporting detail if a reviewer or the text wants a number, but not load-bearing for the argument itself.

**All 5 episodes: clean, positive confirmation.** The same BB-augmented model wins in both the joint and GBM-only
fits in every episode (BAND_BB for T90/EX0/TR1/EX1, SBPL_BB for TR2), kT_BB agrees within 1sigma throughout
(e.g. TR1: 44.05 +/- 1.62 keV joint vs 44.49 +/- 1.54 keV GBM-only; T90: 34.92 +/- 0.80 keV joint vs
35.36 +/- 0.76 keV GBM-only), and the BASE -> BASE+BB step clears the 28.74 threshold by a wide margin in the
GBM-only fit in every episode (Delta-C-stat = 46-382). `status_consistent_with_threshold_gbm` is `True` for all
5 rows. This directly answers the #6 concern: the BB detections do not depend on LAT being in the fit, closing even
the indirect joint-likelihood pathway, not just the effective-area mechanism the weakness doc originally named.

**T90 anomaly (resolved 2026-09-01):** an earlier version of `results.json`'s `GRB131014215GBM` entry had no
`BAND_BB` key for T90 at all -- `any_bb_step_clears_threshold()` (Sec 2) flagged `status_consistent_with_threshold_gbm
= False` for that one row, since `SBPL_BB` was present and cleared the threshold by Delta-C-stat = 491 while the
recorded winner was a non-BB `BAND`. Traced to a labelling mistake on the user's side: the T90 GBM-only BB fit had
been saved as `BAND_PL` instead of `BAND_BB`. Renamed and `results.json` regenerated; BAND_BB is now present
(C-stat = 672.02/448, correctly beating SBPL_BB's 690.20/448) and recorded as `_status: BEST`, matching the joint
fit's own T90 winner. Re-running the script confirms `status_consistent_with_threshold_gbm = True` for all 5
episodes with no remaining exceptions. This is exactly the kind of disagreement the independent
`any_bb_step_clears_threshold()` check exists to catch -- it flagged a real data problem, not the script's.

**BB significance (delta-C-stat plot, added 2026-09-01):** in every episode, GBM-only's Delta-C-stat is *higher*
than the joint fit's, not lower -- T90: 359.7 (joint) vs 381.7 (GBM-only); EX0: 179.6 vs 197.6; TR1: 164.2 vs 180.7;
TR2: 74.3 vs 91.0; EX1 (the weakest case in both fits): 30.4 vs 46.4, still 1.6x the 28.74 threshold. Dropping LAT
does not weaken the case for a blackbody in this burst; if anything the GBM-only likelihood prefers it slightly more
strongly in every single episode, which is the opposite of what the #6 concern (LAT systematics manufacturing a
spurious BB) would predict.

**Size of the effect (fractional-difference plot, paired-MC version):** kT_BB shifts by roughly +0.7% to +3.1%
(GBM-only relative to Joint, all consistent with zero within error) across the 5 episodes; f_BB shifts by roughly
+6.6% to +17.3% (also all consistent with zero within their larger MC errors, but systematically positive in every
episode -- GBM-only always reads *somewhat* higher, never lower). The consistent sign, even though no single episode
is individually significant, is noted as a mild pattern worth keeping in mind rather than a detected effect -- 5
small positive shifts is not itself evidence of a systematic (no combined/stacked significance test was run across
episodes), but it is the kind of pattern a fuller write-up of #6 should mention rather than only quoting the largest
per-episode error bar.

## 4. Judgement calls

- **Log y-axis for the delta-C-stat plot.** Values span roughly one order of magnitude (30-380) across episodes;
  a linear axis would compress EX1 (the case closest to threshold, and therefore the most important one to read
  precisely) into a few pixels while T90 dominates the plot. Log scale keeps the 28.74 threshold line and every
  point legible at once. Decided by Claude.
- **Fractional-difference errors use paired Monte Carlo, not linearized (delta-method) error propagation --
  upgraded 2026-09-01 per user request.** The first version combined each side's already-computed marginal error
  analytically (`sigma(a/b-1) = sqrt((sigma_a/b)^2 + (a*sigma_b/b^2)^2)`), collapsing f_BB's asymmetric MC interval
  to a symmetric proxy to do so. `paired_ratio_pct()` instead draws fresh parameter samples for *each* fit and takes
  percentiles of the ratio distribution directly, which captures whatever asymmetry each side's own posterior has.
  - **A same-seed (common-random-numbers) version was tried first and rejected.** Using the same seed for both
    fits' draws seemed like a reasonable way to "pair" them, but validation showed it was wrong: because both fits
    use the same model form with a similar covariance structure, sharing a seed made the two draws' kt_bb nearly
    perfectly correlated (measured r = 0.9996 for T90's BAND_BB), which collapsed the ratio's spread to an
    artificially tight ~0.15%-wide interval — about 20x tighter than either an independent-seed MC run or the
    original linearized estimate, which agree with each other at the few-percent level. That tightness was a
    seed-sharing artifact (both draws dominated by the same underlying random vector, not by anything physical about
    correlated fits), not a real reduction in uncertainty. Caught by directly checking the correlation coefficient
    and comparing against the earlier delta-method numbers before accepting the result -- the ~20x jump was the
    signal that something was wrong, not a discovery. Fixed by drawing the two fits' samples with different seeds
    (`seed` and `seed + 1`); the two fits are separate RMFIT optimizations with no actual joint covariance available
    to justify treating them as paired in the first place.
- **Legend split into fit-type (colour) and episode (marker), per user request (2026-09-01).** Colour: yellow
  (`JOINT_COLOR = "gold"`, not the project's per-GRB `GRBPlotStyle` colour, since both series here are the *same*
  burst -- colour encodes fit type in this folder, not GRB identity) for Joint, black for GBM-only. Marker: taken
  from `EpisodeMarkerResolver` (`grb_research.grb_utils`), the project's existing per-episode marker convention,
  rather than invented locally -- T90 uses `"s"` to match `bb_flux_fraction.py`'s `T90_MARKERS[1]` for this burst.
  The model name (needed per CLAUDE.md's "identify burst/episode/model" legend rule) is carried on each x-tick
  instead of in either legend, since burst is fixed for the whole figure and every episode uses one model in both
  fits (`models_agree` is True throughout) -- a third legend dimension for a value that never varies would be noise.
- **f_BB and MC settings copied from `bb_fraction/bb_flux_fraction.py`**, not imported, per the cross-folder-import
  convention (CLAUDE.md) — `episode_label`, `analytic_bb_bolometric_flux`, `split_energy_flux`, and the percentile
  MC pattern in `compute_f_bb` are verbatim copies. Same energy band (1 keV - 10 MeV observer-frame), same
  `N_SAMPLES=10000`, `SEED=12345`, so the joint-fit `f_bb_joint` column here should reproduce
  `bb_flux_fraction.csv`'s values for GRB131014A up to MC noise (not bit-for-bit, since the seed produces different
  draws in a differently-shaped call, but the same distribution) — not cross-checked against that file line-by-line,
  since the two are computed independently from the same model objects for a different purpose (a comparison table,
  not the paper's fraction figure). Decided by Claude.
- **No z / cosmology columns.** f_BB in the 1 keV - 10 MeV observer-frame band is redshift- and
  cosmology-independent per-episode (established in `bb_fraction.md`), and this check never touches a rest-frame
  quantity, so those columns would be dead weight per CLAUDE.md's CSV convention. Decided by Claude.
- **`DELTA_CSTAT_THRESHOLD = 28.74` used uniformly.** Every BASE -> BASE+BB step adds exactly 2 free parameters
  (`amp_bb`, `kt_bb`), so the single Delta k = 2 threshold from `section-1-introduction.tex` applies to all four
  base models without needing the Delta k = 4 (36.86) threshold used elsewhere for the further BASE+BB+PL step,
  which isn't part of this specific BB-survives-without-LAT question. Decided by Claude.
- **`any_bb_step_clears_threshold()` checks every base/BB pair present, not just the recorded winner's pair.** An
  earlier version only checked the *recorded* winner's own base model, which silently missed the T90 case (the
  recorded winner, BAND, has no BB counterpart at all in the fit output, so the naive check saw nothing to compare
  and reported nothing wrong). Checking every present pair independently of what's recorded is what surfaced the
  anomaly. Decided by Claude, after the first version produced a false "all consistent" read on a manual spot check.
- **Plot scope kept to the 4 confirmed episodes in the kT comparison plot.** T90 has no `kt_bb_keV_gbm` to plot at
  the time this decision was made (no BB model fit for it yet), so it appeared only as the lone joint-fit point;
  adding a placeholder or annotation for the missing GBM-only T90 point was judged more likely to be misread as "T90
  has no BB without LAT" (not established at the time) than to add clarity. Superseded by the T90 fix (Sec 3) --
  all 5 episodes now plot normally -- but the reasoning is kept here since the same judgement would apply again if a
  future episode is ever missing a point.

## 5. Validation

- Re-ran the `any_bb_step_clears_threshold` check by hand for T90 GBM-only from the raw `results.json` values
  (SBPL: 1181.228515625/450, SBPL_BB: 690.1964111328125/448) — Delta-C-stat = 491.03, matches the script's output.
- Cross-checked that the joint-fit `model_name_joint` column matches the independently-established canonical
  per-episode model list for GRB131014A from `bb_fraction/bb_flux_fraction.csv` (BAND_BB for T90/EX0/TR1/EX1,
  SBPL_BB for TR2) — exact agreement, confirming `GRB.from_dictionary` on the `GRB131014215` key reproduces what
  the rest of the paper already uses.
- `kt_consistent_1sigma` uses a simple combined-error check (`|Delta kT| <= sqrt(err_joint^2 + err_gbm^2)`), not a
  full covariance-aware comparison — sufficient for a sanity check across 5 episodes where the differences are all
  well under 1 combined sigma.
- The paired-MC fractional-difference numbers were cross-checked against the original linearized (delta-method)
  numbers before being accepted (Sec 4) — both approaches agree to within about 1 percentage point on the median and
  give similarly-sized error bars, which is the expected behaviour when the underlying quantities are close enough
  to Gaussian for the linear approximation to hold; the two independent calculations agreeing is itself a check that
  neither has a sign or scale error.

## 6. Known limitations and open questions

- This folder does not update `review-resolution.md` or the paper text on its own — it produces the comparison data
  the user needs to decide how (or whether) to write up #6, per the division of labor already used for every other
  RMFIT-dependent item in this project (the user runs RMFIT, Claude parses/compares).
- `GRB131014215GBM` is a non-standard results.json key (not in `grb_constants.short_to_long`), so it's loaded here
  via `GRB.from_dictionary` directly rather than through `prepare_grbs`/`GRBCatalog`, which assumes the standard
  short-name mapping. This is fine for a one-off comparison script but would need a name-mapping entry if this key
  is ever consumed by another script in this project.
- The paired-MC ratio still doesn't use the *true* joint covariance between the two fits (which doesn't exist --
  they're separate RMFIT optimizations), so it's still an approximation: independent-seed draws are a reasonable
  default in the absence of that information, not a claim that the two fits are actually statistically independent
  (they share the same underlying GBM photon counts). No sensitivity check was run on how much the paired-MC
  intervals would change under a partial, hand-specified correlation.
