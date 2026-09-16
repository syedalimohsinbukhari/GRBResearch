# Phase 5 — Norris-profile t_v measurement (working plan)

Companion to `PLAN.md` Phase 5 (project root, one level up) and `codes-for-paper/lorentz_factor/lorentz_factor.md` §6b. This is a working note, not an approved/locked spec — update it as the implementation surfaces things the plan got wrong.

## Status update (2026-09-16) — automated pipeline superseded, plan below is stale on that point

The automated pipeline this plan describes (MEPSA/`scipy.signal.find_peaks` peak detection + a per-window local `pymultifit` fit, living in `GRBResearchWork/variability_timescale/`) stalled on GRB080916C's TR2 (zero MEPSA detections even with padding, because TR2 is one broad, smoothly-declining pulse with no local excess for a spike detector to find). Its source files (`norris_fit.py`, `variability_timescale.py`, `light_curves.py` copy, `mepsa_peaks.py`, `fit_tr234_joint.py`) have since been deleted from disk — only stale `__pycache__/*.pyc` remain, and the directory was never committed (`??` in `git status`), so there is no commit to recover them from either.

Work continues instead as a **separate, manually-driven, GRB-by-GRB joint-Norris-fit track** in `codes-for-paper/variability_analysis/` (see that folder's `variability_analysis.md`), started 2026-09-07. So far this covers GRB131014A, GRB140206B, and GRB231129C; GRB080916C (where the automated approach stalled) has not yet been redone manually. The "New files" section below (§ "New files — `GRBResearchWork/variability_timescale/`") describes the abandoned automated layout and no longer reflects what's on disk — kept here for the historical methodology reasoning (energy band, eq. 10 sourcing, EX0/EX1 window-inheritance rule, etc.), which still applies to the manual track, not for the file list.

**Broken-import state, fixed 2026-09-16 (BUG-21, now closed).** `codes-for-paper/variability_analysis/fitter*.py` and `experiments/window_sensitivity_GRB231129779/window_sensitivity.py` used to `import` from the now-deleted `variability_timescale.norris_fit` / `variability_timescale.light_curves` modules. That folder now carries its own local `norris_fit.py` (a proper rename of the old `norris..py`) and a copied-in `light_curves.py`, and every script has been repointed to import those instead — confirmed resolving cleanly in the repo's `.venv`. See `BUGS.md` BUG-21.

## Context

`lorentz_factor.py`'s γγ-opacity Γ_min calculation currently uses each episode's full duration as the variability timescale t_v — a deliberately conservative upper bound (`VARIABILITY_TIMESCALE` is an empty override dict, so every episode falls through to `interval.end - interval.start`). `lorentz_factor.md` §6b already assessed replacing this with a measured t_v from fitting the Norris et al. (2005) pulse profile to the light curve, following the precedent of Bukhari et al. 2022 (the user's own prior paper), but marked it "feasible, not done." This session worked out the actual methodology through direct back-and-forth, because the naive version (fit one Norris pulse per episode window) breaks down: several episodes visibly contain multiple pulses, T90/TR3 share the same γγ-opacity-defining LAT photon at t≈40.5s where the light curve looks nearly featureless, and the original Norris (2005) paper isn't in the project's library so even the τ1/τ2→t_v formula had to be tracked down externally.

## Decisions locked in this session

- **Scope**: all four bursts get t_v measurements (light-curve data exists for all four), even though only GRB080916C's currently feeds a published Γ_min.
- **Deliverable boundary**: this produces a standalone, documented CSV + diagnostic plots only. It does **not** touch `VARIABILITY_TIMESCALE`, does not rerun Γ_min/Limit B, and does not touch any `.tex` file. Adopting the measurement into the paper is a separate, explicit follow-up once the fit quality is reviewed — matching how every other correction in this project was verified before being propagated.
- **T90 excluded entirely**: spans the whole burst, not a single/few pulses — no fit attempted, no CSV row. `t_v_source` for T90 in the existing pipeline is untouched. (User separately floated Bayesian-Block/MVT-based local t_v for T90 as *optional future work* — not built here.)
- **Multi-pulse decomposition, not one-Norris-per-episode**: the Norris profile is inherently single-peaked, and the user's own visual read of the 10–400 keV light curve found multiple pulses in several episodes (GRB080916C: TR1-window=3, TR2=1, TR3=multiple, TR4=1, TR5-window=1). Each fit window gets a sum-of-N-Norris-pulses fit; the reported t_v for a given episode comes from whichever fitted sub-pulse's peak time is closest to that episode's own Γ_min-defining LAT photon arrival time (`t_arr_s`, already in `lorentz_results.csv`).
- **EX0/EX1 inherit from the wider window, not fit separately**: per `CLAUDE.md`'s own EX semantics (EX0 shares its paired TR's *end* and starts earlier; EX1 shares its paired TR's *start* and ends later), the EX interval's boundary is a strict superset of its paired TR's in every case checked against `results.json` for all four bursts (verified below). So the EX boundary is used as the single fit window, and the paired TR just reads its own t_v off the same decomposition.
- **t_v per Bukhari et al. (2022) eq. (10)**: `t_v = (τ2/2)·√[(ln2 + 2√(τ1/τ2))² − 4τ1/τ2]`, attributed there to Norris et al. (2005) and described as "the half width of the pulse at half maximum." This is the literal, primary-cited formula — found by reading page 10 of `2022_06_06_..._Bukhari_Urooj_Murtaza.pdf` directly (the equation had been mangled into unreadable glyphs by `pdftotext`'s text extraction, which is why an earlier grep-based search of the same PDF missed it). Verified numerically: for a test pulse, eq. (10)'s value matched FWHM/2 (computed from actual half-max crossing points on both the rise and decay sides of an asymmetric pulse) to 5+ significant figures. This supersedes the earlier `τ_peak=√(τ1τ2)` candidate (a Nemiroff-2011-sourced *proxy*, not the real formula) — eq. (10) is the actual thing, not a stand-in for it. Note `t_peak = t_s + √(τ1τ2)` is a *separate* quantity, still used for sub-pulse selection (step 5 below) — it is no longer equal to t_v.
- **Energy band: 10–400 keV only.** No fitting in other bands. Bukhari et al. 2022 also tried 400–900 keV and 250 keV–5 MeV as a cross-check for R_GeV≈R_MeV, but that's explicitly not being replicated here — one band, one fit per window.
- **Null result is a valid outcome**: if a window shows no fittable pulse near the relevant photon's arrival time (expected for TR3's 27 GeV photon at t≈40.5s, per the user's own read of that region as "nearly featureless"), that gets recorded as `t_v_source="no_reliable_fit"`, not papered over with a forced fit.
- **Fitting library: `pymultifit`**, the user's own package (`github.com/syedalimohsinbukhari/pyMultiFit`), used for the multi-peak fitting instead of a hand-rolled `scipy.optimize.curve_fit` sum — install via `uv add pymultifit` into `GRBResearchWork/.venv` (should pull v1.0.9, the latest published release).
- **New files live in a separate top-level folder** — `GRBResearchWork/variability_timescale/`, not under `codes-for-paper/` — since this is measurement/exploration work, not a paper-ready figure/table script.

## Episode grouping (verified against `results.json` for all four bursts)

15 fit windows total, producing 22 episode-level t_v rows (26 LAT-covered episodes across the sample, minus 4 T90 rows).

**GRB080916C** (`GRB080916009`):
| Window (boundary) | Episodes read | Defining photon(s): t_arr [E_max] |
|---|---|---|
| EX0 (−0.128–4.864) | EX0, TR1 | both t=3.03 s [301.2 MeV] — same photon |
| TR2 (4.864–15.040) | TR2 | t=6.86 s [2110 MeV] |
| TR3 (15.040–55.296) | TR3 | t=40.50 s [27.4 GeV] — also T90's, but T90 excluded |
| TR4 (55.296–59.520) | TR4 | t=55.58 s [464 MeV] |
| EX1 (59.520–67.904) | TR5, EX1 | TR5: t=60.75 s [340 MeV]; EX1: t=65.52 s [989 MeV] — distinct sub-pulses |

**GRB131014A** (`GRB131014215`):
| Window | Episodes | Defining photon(s) |
|---|---|---|
| EX0 (−0.192–2.432) | EX0, TR1 | both t=1.98 s [1233 MeV] |
| EX1 (2.432–6.976) | TR2, EX1 | TR2: t=3.18 s [1021 MeV]; EX1: t=4.31 s [1193 MeV] |

**GRB140206B** (`GRB140206275`, no trailing excess — matches `CLAUDE.md`):
| Window | Episodes | Defining photon(s) |
|---|---|---|
| EX0 (4.288–11.072) | EX0, TR1 | both t=7.56 s [225 MeV] |
| TR2 (11.072–20.032) | TR2 | t=17.94 s [118 MeV] |
| TR3 (20.032–26.752) | TR3 | t=24.00 s [753 MeV] — also T90's |
| TR4 (26.752–59.776) | TR4 | t=27.50 s [576 MeV] |
| TR5 (59.776–100.032) | TR5 | t=81.57 s [194 MeV] |
| TR6 (100.032–154.240) | TR6 | t=105.49 s [526 MeV] |

**GRB231129C** (`GRB231129779`):
| Window | Episodes | Defining photon(s) |
|---|---|---|
| EX0 (−0.192–3.136) | EX0, TR1 | both t=0.68 s [123 MeV] |
| EX1 (3.136–10.048) | TR2, EX1 | both t=3.84 s [723 MeV] — same photon |

## Method, per fit window

1. **Light curve**: sum 10–400 keV background-subtracted rate across that burst's NaI `.dat` files, reusing exactly what `light_curves/make_lightcurve.py::make_lightcurves` already does — call `light_curves/light_curves.py::lightcurve_data(dat_file, energy_low=10, energy_high=400, errors=True)` per detector (default band, same as the existing driver scripts), sum across detectors, then `y = rate - background`, `y_err = sqrt(rate_err**2 + background_err**2)`.
2. **Slice** to the window boundary padded by 20% of its raw duration on each side (minimum 0.5 s) — a Claude-proposed default, recorded in the method note for review, not silently baked in.
3. **Peak detection**: `scipy.signal.find_peaks` on `y`, with height/prominence thresholds set from the pre-burst (t<0) baseline noise level — reported per window as `n_pulses_detected`, cross-checked against the user's own visual counts where given (GRB080916C only) and any disagreement flagged, not silently resolved either way.
4. **Multi-pulse fit**: via `pymultifit` — sum of N Norris pulses, `I(t) = Σ_i A_i·exp(2√(τ1ᵢ/τ2ᵢ))·exp(−τ1ᵢ/(t−tsᵢ) − (t−tsᵢ)/τ2ᵢ)` for `t>tsᵢ` (else 0 for that term) plus a constant offset, seeded from the detected peaks. Confirm `pymultifit`'s actual custom-fitter API by reading its installed source before writing this — it wasn't pinned down from docs alone this session.
5. **Sub-pulse selection**: for each episode reading this window, take the fitted pulse whose `t_peak = ts + √(τ1τ2)` is closest to that episode's `t_arr_s` (from `lorentz_results.csv`); record the offset `|t_peak − t_arr|` so a bad association is visible.
6. **t_v + uncertainty**: `t_v = (τ2/2)·√[(ln2+2√(τ1/τ2))²−4τ1/τ2]` for the chosen pulse (Bukhari et al. 2022 eq. 10 — see sourcing above). Propagate via MC: draw `N_SAMPLES=10_000` parameter vectors from the fit's covariance using `get_rng(seed=12345)` (same constants as `lorentz_factor.py`), recompute t_v per draw, take 16/50/84 percentiles. Note this is a separate computation from step 5's `t_peak = t_s + √(τ1τ2)` — the two are no longer the same quantity (t_v is a width, t_peak a location), so both need to be carried through.
7. **Null result**: if the fit doesn't converge, the chosen pulse's τ1 or τ2 has >100% relative error, or no detected peak falls within ~1 s of `t_arr`, record `t_v_source="no_reliable_fit"` and leave t_v blank — a legitimate outcome, expected for TR3's window given the user's own read of it as near-featureless around t≈40.5s.

## New files — `GRBResearchWork/variability_timescale/`

- `norris_fit.py` — the Norris pulse function wired into `pymultifit`'s custom-fitter mechanism, plus the `t_v = √(τ1τ2)` derivation.
- `variability_timescale.py` — driver: builds the 15 windows from `results.json` (via existing `grb_research`/`TimeInterval` parsing — reuse, don't hand-derive boundaries) + `lorentz_results.csv` (defining photons), loads light curves via `light_curves/light_curves.py::lightcurve_data` (add `light_curves/` to `sys.path` since it's a standalone module outside `src/grb_research`), fits via `pymultifit`, and writes the CSV.
- `variability_timescale.csv` — 22 rows. Columns: `grb_name`, `episode`, `window_label`, `t_start_s`/`t_stop_s` (actual padded fit window), `energy_low_keV=10`, `energy_high_keV=400`, `n_pulses_detected`, `n_pulses_visual` (blank unless manually cross-checked), `pulse_index_used`, `t_arr_s`, `t_peak_offset_s`, `A`/`t_s`/`tau1`/`tau2` + MC errors for the chosen pulse, `t_v_s`, `t_v_err_lower_s`, `t_v_err_upper_s`, `t_v_source` (`norris_fit`/`no_reliable_fit`), `t_v_definition` (constant string: `"t_v=(tau2/2)*sqrt((ln2+2*sqrt(tau1/tau2))^2-4*tau1/tau2), Bukhari et al. 2022 Adv.Space Res. eq.10, attributed to Norris et al. 2005; = FWHM/2"`), `n_samples`, `seed`.
- Diagnostic plots (one per burst, multi-panel) — light curve + fitted multi-pulse overlay + `t_arr` markers, `update_style()`/`GRBPlotStyle`, saved `.png`+`.pdf`.
- `variability_timescale.md` — required per-folder method note: defining equation (Norris et al. 2005 pulse form, as reformulated and cited in `nermioff_mnras.pdf` §2.1/2.2, itself in the repo's Literature Review and read directly — Norris 2005 itself still isn't held locally, so this is one citation removed rather than fully primary), every judgement call above attributed to this session's decisions, and a limitations section covering the τ_peak-as-proxy-for-w caveat and the optional future Bayesian-Block/MVT approach for T90.

## Reuse

- `light_curves/light_curves.py::lightcurve_data()` — reader, called exactly as `make_lightcurve.py` already does.
- `src/grb_research/grb_calculations.py::get_rng()` — RNG convention.
- `codes-for-paper/lorentz_factor/lorentz_results.csv` — defining-photon lookup by `grb_name`+`episode`.
- `results.json` via existing `grb_research` interval parsing — boundaries.
- `src/grb_research/__init__.py::update_style()`, `grb_styles.py::GRBPlotStyle` — plot conventions.

## Explicitly out of scope

- No changes to `lorentz_factor.py`/`lorentz_factor_limit_b.py`, no table regeneration, no `.tex` edits.
- No T90 fitting, any burst.
- No Bayesian-Block/MVT implementation (noted as future work only).
- No multi-band fitting (no 400–900 keV or 250 keV–5 MeV cross-check, unlike Bukhari et al.'s R_GeV≈R_MeV exercise) — 10–400 keV only.
- No attempt to track down the original Norris (2005) paper itself — no longer needed for the t_v formula (eq. 10 is cited directly from Bukhari et al. 2022, read from the source PDF), but still absent from the project's library if anyone wants the primary derivation.

## Verification

- Unit-checked (done, 2026-08-22): `norris_pulse()`'s peak occurs at `t_peak = t_s + √(τ1τ2)` with value exactly equal to the amplitude parameter `A` (confirmed for a test case: A=60, τ1=4, τ2=25 → peak at t=10, value=60.0 exactly) — this is the expected behavior of Bukhari eq. (9)'s parameterization, where the leading `exp(2√(τ1/τ2))` prefactor is specifically there to make `A` the peak value (unlike Nemiroff's bare-exponential parameterization, where the peak is suppressed below `A`). `tv_value()` (eq. 10) verified against a direct numerical half-max-crossing search on the same test pulse: matched FWHM/2 to 5+ significant figures (15.759365 vs 15.759357).
- Run `variability_timescale.py` end-to-end (`PYTHONPATH=src`, from `GRBResearchWork/`): confirm 15 windows attempted, 22 CSV rows, no unhandled exceptions.
- For GRB080916C, compare `n_pulses_detected` against the user's visual counts and report any mismatch.
- Explicitly check the TR3 window (27 GeV photon, t≈40.5s) and report whether it comes back `no_reliable_fit` as expected, or finds a pulse — either way, report it, don't assume.
- Not a paper build — no `pdflatex`/`pdftotext`, so the `VERIFY.md` convention doesn't apply; the CSV + plots are the review artifact.