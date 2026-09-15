# Variability timescale (t_v) via manually-seeded joint Norris fit — GRB131014A

Companion to `PHASE5_TV_PLAN.md` (project root), which documented the earlier automated
(MEPSA/scipy peak-detection + per-window local fit) pipeline for this same Phase 5 goal. That pipeline's
source files (`variability_timescale/norris_fit.py`, `variability_timescale.py`, etc.) have since been
deleted — see `PHASE5_TV_PLAN.md`'s 2026-09-16 status update — so this folder is now the only live track.
This folder is a **separate, manually-driven track**, worked one GRB at a time starting with GRB131014A on
2026-09-07 after the automated pipeline stalled on GRB080916C's TR2 (zero MEPSA detections even with
padding, because TR2's true shape is one broad, smoothly-declining pulse with no local excess for a spike
detector to find — see the MEPSA-fallback discussion in this session's history). GRB140206B
(`fitter_GRB140206275.py` + a `_simple` variant) and GRB231129C (`fitter_GRB231129779.py`) were fit in
this same manual style in a later session (2026-09-15/16) but are **not yet written up below** — the
"What the code computes" / "Every judgement call" / "Results" sections that follow describe the GRB131014A
fit only; see each script's own docstring/comments for the other two bursts until this note is extended.
GRB080916C has no manual fit yet. Not yet wired into `lorentz_factor.py`'s `Gamma_min` pipeline — same
deliverable-boundary stance as the abandoned automated pipeline.

## What the code computes

Same Norris (2005) pulse and $t_v$ definition as `variability_timescale/norris_fit.py` (Bukhari et al. 2022,
Adv. Space Res., eq. 9/10) — this folder's `fitter.py`/`fitter_CLAUDE_GRB131014215.py` import that module directly rather
than redefining it. $t_v = (\tau_2/2)\sqrt{(\ln2+2\sqrt{\tau_1/\tau_2})^2-4\tau_1/\tau_2}$, the FWHM/2 of the
fitted pulse.

The fit itself: a single joint `NorrisFitter` (pymultifit `BaseFitter`, 5 pulses, 20 free parameters) over
one fixed window, `[-1, 10]` s since trigger, of the 10–400 keV NaI background-subtracted light curve
(GRB131014215, detectors `na`/`nb`/`n9`), rather than the automated pipeline's per-episode local windows.

## Every judgement call, and who made it

- **Manual joint fit over automated per-window detection.** After the automated MEPSA-based pipeline
  (`variability_timescale/`) hit real, GRB080916C-specific problems (TR2 needing a direct-max fallback,
  TR4 landing on a degenerate solution), the user instead hand-fit GRB131014A directly with carefully
  chosen $p_0$ seeds spanning the whole EX0/TR1/TR2/EX1 region at once. **User decision**, motivated by
  wanting a clean, verifiable reference fit for this burst before deciding whether/how to generalize the
  automated pipeline further.
- **Amplitude/shape decoupling.** The fit runs on the light curve normalized to peak=1 (for numerical
  conditioning — feeding raw counts/s directly into $p_0$ made the optimizer behave badly), and $A$ is
  rescaled back to physical counts/s *after* fitting via the recorded `Y_MAX_CTS_PER_S`, exact since $A$ is
  purely multiplicative in the Norris formula. **Claude decision**, prompted directly by the user's own
  observation that pre-scaling $p_0$'s amplitude by $y_\text{max}$ "screws the fitter over" — confirmed the
  root cause is optimizer conditioning, not a flaw in the rescale itself.
- **Guess-sensitivity used as the actual reliability criterion, not `mc_kept_fraction` alone.** Refitting
  with every $\tau_1=\tau_2=1$ (a neutral, uninformative shape seed, only $A$/$t_s$ roughly hinted) reproduced
  pulses 3, 4, 5 to within 0.01 s in $t_\text{peak}$ and comparable `mc_kept_fraction` (0.77–1.00 either way),
  but landed pulses 1 and 2 on genuinely different $t_\text{peak}$/$t_v$ each time, with `mc_kept_fraction`
  staying poor (0.10–0.21) under both seedings. **Jointly established**: the user proposed the neutral-seed
  test as an argument that careful guessing wasn't load-bearing; the result confirmed it only for pulses
  3–5, and additionally showed pulses 1/2 are unreliable *because* they sit in the earliest, heavily
  overlapping part of the light curve (Hakkila 2026's "overlapping pulses are not unique" failure mode), not
  because of any fitter/seeding problem that better guessing would fix.
- **$(\tau,\xi)$ reparameterization tried and found not to help, for this fit specifically.** Converting the
  same $p_0$ via `norris_reparam.tau_xi()` and refitting with `NorrisReparamFitter` landed pulses 1–3 on
  different local minima entirely (pulse 2 converged to a *worse* degenerate solution, $t_v=252$ s in an
  11 s window), and where it did land near the same point (pulses 4, 5), the $t_s$–$\tau$/$t_s$–$\xi$
  correlations stayed at 0.94–0.99 — no better than the original $t_s$–$\tau_1$ correlation it was meant to
  fix. **Negative result, recorded so it is not re-tried blindly**: Recipe C's premise
  (`variability_timescale/norris_shenanigans.md`) doesn't hold for this joint 5-pulse fit.
- **Photon-to-pulse assignment: nearest preceding pulse onset.** Each LAT photon is assigned to the pulse
  with the largest $t_s \le$ the photon's arrival time — i.e. whichever pulse has "turned on" and not yet
  been superseded by the next pulse's onset. **Claude-proposed, verified against the user's own manual
  identification** of two specific photons (445.28 MeV at $t=2.533$ s → pulse 4; 552.79 MeV at $t=4.111$ s →
  pulse 5) before being generalized to all 10 photons in the window. Not the only possible rule (nearest by
  absolute time, or nearest to $t_\text{peak}$, would be alternatives) but this one is physically motivated
  and reproduced the verified cases exactly.
- **Episode boundaries from `results.json` directly**, not eyeballed or taken from the automated pipeline's
  photon-centered heuristic windows: `EX0 -0.192_2.432`, `TR1 0.960_2.432`, `TR2 2.432_4.160`,
  `EX1 2.432_6.976`. A pulse can match more than one episode (EX0 fully contains TR1's window, EX1 fully
  contains TR2's) — matching is "does this pulse's $t_\text{peak}$ fall inside this episode's boundary", not
  a 1:1 assignment.

## Conventions that could plausibly have gone another way

- The fixed `[-1, 10]` s window and 5-pulse count were chosen by inspection for this burst specifically —
  not derived from a general rule, and not expected to transfer unchanged to another GRB.
- $p_0$ itself (in `fitter.py`) is the user's own hand-tuned values, arrived at by trial and error against
  the normalized light curve; left untouched throughout this session's work per explicit instruction.
- The twin-axis photon overlay's y-limits are set to `(0, data.max())` on both axes (no autoscale padding)
  specifically so 0 and each series' own peak align visually — a display choice, not a data transform.

## Validation performed

- **Neutral-seed reproducibility test** (above): the actual mechanism used to separate "real, data-determined"
  pulses (3, 4, 5) from "seed-dependent, non-unique" ones (1, 2) — not just asserted, run twice with visibly
  different $\tau_1,\tau_2$ starting points and compared.
- **Photon table cross-checked against `LAT_analysis/014__GRB131014215/` and `lorentz_results.csv`**: the
  1232.92 MeV (TR1/EX0), 1021.24 MeV (TR2), 1193.12 MeV (EX1) "defining" photons already used elsewhere in
  the project appear at the same arrival times in this window's full photon list, confirming
  `GRB131014215_lat.fits` (copied here from `light_curves/GRB131014215/lat.fits`, not re-derived) and the
  hardcoded `T0_MET_S = 403420143.2` (read from
  `LAT_analysis/014__GRB131014215/Ep1__0.960_2.432/GRB131014A_fit_results_0.96_2.432.txt`'s own `T_0` line)
  are correct: the FITS file's own GTI (403420144.16–403420147.36) equals `T0_MET_S+0.96` to `T0_MET_S+4.16`,
  exactly the T90 window.
- $A_\text{cts/s} = A_\text{norm}\times Y_\text{MAX\_CTS\_PER\_S}$ is an exact linear identity, not an
  approximation — verified by re-plotting both data and every fitted curve at physical scale and confirming
  the total fit still tracks the data (a scale error would show as a visible mismatch; none appeared).

## Window-widening as a pulse-count diagnostic (found via GRB231129C, applies generally)

Found while investigating a *different* burst (GRB231129C, `experiments/window_sensitivity_GRB231129779/`)
but the conclusion applies to every joint Norris fit in this project, including the GRB131014A one
documented above — recorded here since this is the method note for that fit family.

**The check:** refit the same pulses, from the same $p_0$, over two different window widths (e.g.
$[-1,10]$s vs $[-10,20]$s), and compare not just each pulse's $t_v$ but its `mc_kept_fraction`. This is
possible because of a specific, verified mechanism: `NorrisFitter.fit_boundaries()`
(`variability_timescale/norris_fit.py:49-53`) sets $t_s$'s lower bound to `x_values.min()` — the fit
window's own left edge is $t_s$'s box constraint — so widening the window directly widens how far the
optimizer can push $t_s$ along the already-established near-total $t_s$↔$\tau_1$ degeneracy (this session's
correlation matrices, TR4, the $(\tau,\xi)$ reparam test above).

**What it found for GRB231129C, concretely:** a 4-pulse fit's pulse 3 kept its $t_v$ point estimate stable
under the window change (0.505→0.490, ~3%) but its `mc_kept_fraction` collapsed (0.451→0.148) — the median
survived, the reliability didn't, and checking $t_v$ alone would have missed this entirely. Adding a 5th
pulse (there was real structure near $t\approx4.6$–$5.5$s the 4-pulse model had nowhere to put) made the
instability disappear outright: every pulse's $t_v$ shift dropped to <10%, and `mc_kept_fraction` *improved*
for all five pulses going from narrow to wide, instead of collapsing for one.

**The methodological conclusion:** window-sensitivity is a **necessary, not sufficient**, diagnostic for an
under-specified pulse count. If parameters (especially `mc_kept_fraction`, not just $t_v$) are stable
under a window change, that's real evidence the count is right. If they're not, the confirming step is
adding a pulse and checking whether the instability actually goes away — the instability alone is a red
flag, not proof; a genuinely isolated pulse could in principle be unstable for its own unrelated reasons.

**Flagged follow-up: GRB131014A's pulses 1 and 2** (this fit, above) were called "unresolved" based only on
a seed-sensitivity test (different $\tau_1,\tau_2$ seeds landing on different answers), **not** this
window-widening check. Given what this diagnostic just caught for GRB231129C, it's an open question whether
pulses 1/2's instability is a genuine, irreducible degeneracy in that heavily-overlapping early structure,
or whether it's actually a symptom of the same under-fitting problem (a 6th pulse missing from the
$[-1,10]$s window). Not yet tested — the next thing to check before trusting the "genuinely non-unique"
conclusion in the Known Limitations section below.

## Results (this session, 2026-09-07)

| pulse | episode(s) | $t_\text{peak}$ [s] | $t_v$ [s] | kept | anchoring photon |
|---|---|---|---|---|---|
| 3 | TR1 (+EX0) | 1.766–1.774 | $0.375$–$0.381$ | 0.91–0.999 | 1232.92 MeV @ 1.984 s (+6 more, all in decay) |
| 4 | TR2 (+EX1) | 2.697–2.699 | $0.552$–$0.561$ | 0.999–1.00 | 445.28 MeV @ 2.533 s |
| 5 | TR2 (+EX1) | 3.589–3.590 | $0.266$–$0.269$ | 0.77–0.85 | **1021.24 MeV @ 3.184 s** (TR2's existing $\Gamma_\text{min}$ photon) + 552.79 MeV @ 4.111 s |

## Known limitations and open questions

- **Pulses 1 and 2 are unresolved**, not merely unmeasured: two different, reasonable seedings give two
  different answers and both stay poorly MC-constrained. No LAT photon falls in either pulse's active
  window (first photon arrives at 1.852 s), so this doesn't block feeding TR1/TR2's $\Gamma_\text{min}$ —
  but it does mean the earliest part of this burst's variability structure is not characterized here.
  **Not yet re-checked against the window-widening pulse-count diagnostic** (see that section above,
  found via GRB231129C) — this "unresolved" conclusion currently rests only on seed-sensitivity, and it's
  an open question whether pulses 1/2 are a genuine irreducible degeneracy or a symptom of a missing 6th
  pulse in the $[-1,10]$s window, the same way GRB231129C's pulse 3 turned out to be.
- **TR2 now has two candidate anchor pulses**, not one: pulse 4 (445 MeV photon, extremely well-constrained,
  $t_v\approx0.56$ s) and pulse 5 (the *already-published* 1021 MeV defining photon, less tightly
  constrained, $t_v\approx0.27$ s). Which one should replace TR2's current duration-based $t_v$ in
  `lorentz_factor.py`, if either, is an open decision — not resolved here, and not yet fed back into that
  pipeline (same as the rest of Phase 5).
- **This is one GRB, hand-tuned.** Whether the neutral-seed reliability test, the photon-assignment rule, or
  the 5-pulse/[-1,10]s window choice generalize to the other three bursts is untested; the user's plan is to
  work through them manually one at a time.

## Files here

GRB131014A (documented above):
- `fitter.py` — the user's own working file, reused/repurposed across bursts as this track progressed
  (currently set up for GRB140206B, not GRB131014A — see below); not touched by Claude past the two
  changes explicitly requested for the GRB131014A stage (episode-boundary `axvline`s, fixing a broken
  `savefig` path).
- `fitter_CLAUDE_GRB131014215.py` — Claude's copy, used for the amplitude-rescale, physical-unit plot, LAT-photon overlay,
  and photon-to-pulse assignment work documented above. Diverges from `fitter.py` from that point on.
- `GRB131014215_lat.fits` — copied from `light_curves/GRB131014215/lat.fits`; the FITS-reading convention
  (`fits.open(...)[1].data`, per-source probability column named after the source) is copied from
  `light_curves/make_lightcurve.py`, not re-derived.
- `norris_fit_results_GRB131014215.csv` — one row per (pulse, episode) match, produced by `fitter_CLAUDE_GRB131014215.py`.
- `norris_fitted_GRB131014215.png/.pdf` — the decorated plot (physical units, episode boundaries, LAT
  photon overlay on twin axis).

GRB140206B and GRB231129C (files exist, **not yet written up** in this note — see "not yet written up" above):
- `fitter_GRB140206275.py`, `fitter_GRB140206275_simple.py`, `GRB140206275_lat.fits`,
  `norris_fit_results_GRB140206275.csv`, `norris_fit_results_GRB140206275_simple.csv`,
  `norris_fit_diagnostics_GRB140206275.csv`, `norris_fit_GRB140206275.png`,
  `norris_fitted_GRB140206275.png/.pdf`, `norris_fitted_GRB140206275_simple.png/.pdf`.
- `fitter_GRB231129779.py`, `GRB231129779_lat.fits`, `norris_fit_results_GRB231129779.csv`,
  `norris_fit_GRB231129779.png` (+ `__bkp` variant), `norris_fitted_GRB231129779.png/.pdf`.
- `norris..py` — a full, working `NorrisFitter`/`norris_pulse`/`tv_value`/`tv_mc_summary` implementation
  already committed in this folder (pre-Phase-5, `main-minor-75`), matching the exact interface the
  `fitter*.py` scripts need. It is almost certainly what those scripts *should* be importing locally
  (per `CLAUDE.md`'s "copy rather than fight `sys.path`" convention) instead of the broken
  `from variability_timescale.norris_fit import NorrisFitter` (see `PHASE5_TV_PLAN.md`'s bug note) — but
  the filename has a stray extra dot (`norris..py`, not `norris_fit.py`) so no current script actually
  imports from it under that name. Flagged, not fixed, in this pass: renaming it and repointing the
  `fitter*.py` imports is a code change, out of scope for a docs-staleness pass.
- `dry_run.png` — a dry-run diagnostic plot; which script/burst produced it is not documented here.
- `experiments/window_sensitivity_GRB231129779/` — the window-widening pulse-count diagnostic (see
  section above), `window_sensitivity.py` + its `window_sensitivity_results.csv`.
