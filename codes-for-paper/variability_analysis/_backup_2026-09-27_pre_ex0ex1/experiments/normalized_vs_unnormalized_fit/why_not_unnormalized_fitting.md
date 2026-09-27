# Why the Norris-pulse fits use peak-normalized light curves, not raw counts/s

Standalone reviewer-facing note, 2026-09-23. Written to stand on its own if a referee asks directly why
this project's variability-timescale pipeline (`../../fitter*.py`, all four bursts) fits
$y(t)/\max(y)$ rather than the raw background-subtracted count-rate curve $y(t)$. Full investigation,
method, and every intermediate result are in `comparison.md` in this same folder — this file is the
short answer with the evidence, not a replacement for it.

## The question

A referee could reasonably ask: absolute count rate is the physically meaningful quantity, and the fitted
amplitude $A$ has direct physical content (peak flux). Dividing by $Y_\text{max}$ before fitting looks like
an arbitrary transform of the data. Does it bias the fit, or hide something a raw-scale fit would show?

## The short answer

**No. Tested directly, not assumed.** Peak-normalization is a numerical solver convenience for
`scipy.optimize.curve_fit`'s Levenberg-Marquardt/trust-region step sizing, not a change to what is being
fit. An equally valid unnormalized method reproduces the same fit quality and the same well-constrained
pulse parameters, for all four bursts in this project's sample, without touching the data. Where the two
methods *did* disagree, it flagged two genuinely weak or degenerate pulses that needed attention on their
own merits — a diagnostic benefit of running both, not evidence against normalization.

## What was actually compared

Two ways of fitting the same sum of Norris (2005) pulses,
$I(t)=A\exp[2\sqrt{\tau_1/\tau_2}]\exp[-\tau_1/(t-t_s)-(t-t_s)/\tau_2]$, to the same background-subtracted,
summed-NaI-detector light curve:

1. **Normalized (production method).** Fit $y(t)/\max(y)$ with `NorrisFitter` (this project's wrapper
   around `pymultifit.BaseFitter`, which itself calls plain `scipy.optimize.curve_fit` with bounds but no
   scaling option).
2. **Unnormalized, explicit `x_scale`.** Fit the raw $y(t)$ directly via `scipy.optimize.least_squares`,
   telling the optimizer each parameter's natural scale explicitly (amplitude $\sim Y_\text{max}$,
   $t_s,\tau_1\sim5$s, $\tau_2\sim1$s) since raw counts/s and timing parameters otherwise differ by 3–5
   orders of magnitude, which is badly conditioned for a generic optimizer.

`scipy`'s *automatic* per-parameter scaling (`x_scale='jac'`) was tried first and **rejected** — it
converges, but to a measurably worse optimum (SSE $\approx1.03\times10^8$ vs. the normalized method's
$5.71\times10^7$ for the burst it was first tested on) — no better than no scaling at all. Only an
explicit, hand-specified scale array works, and that array cannot be passed through `pymultifit`'s
`BaseFitter.fit()` at all (it hardcodes the `curve_fit` call), so using it in production would mean
maintaining a second, parallel fitting implementation outside the one every other script in this project
already shares.

## The evidence

Overall fit quality, SSE ratio (unnormalized / normalized), computed in the same physical counts/s$^2$
units both ways:

| GRB | SSE ratio |
|---|---|
| GRB080916C (6-pulse, post BUG-23 pulse-5 drop) | 0.9999 |
| GRB131014A | 0.9991 |
| GRB140206B (production window) | 0.9986 |
| GRB231129C | 1.0000 |

Every ratio within 0.15% of 1. Most individual pulse parameters ($t_\text{peak}$, $t_v$) agree between
methods to better than 1%, across the full 23-pulse sample.

**Where the two methods disagreed, the disagreement was diagnostic, not a normalization artifact:**

- **GRB080916C's original 7th pulse** diverged sharply between methods ($t_v=0.046$ vs. $0.153$s).
  Independently, in production, that same pulse converged to a degenerate near-delta-function fit
  (`tau2`$\to0.009$, `mc_kept_fraction` $0.52\to0.05$) — confirmed against the raw data to be overfitting
  two noisy bins, not resolving a real feature. It was dropped from the model. The
  normalized-vs-unnormalized comparison caught this independently of, and before, the production
  diagnosis was cross-checked — two unrelated methods agreeing that a specific pulse was unreliable is
  stronger evidence than either alone.
- **GRB140206B's pulses 6 and 7** similarly diverged at the production fit window. Follow-up (window
  widening, then a targeted reseed — see `comparison.md`'s "GRB140206B follow-up" section) traced this to
  genuine parameter degeneracies in the light curve itself ($t_s\leftrightarrow\tau_1$ trade-offs, a
  window-edge pinning that only released at the full data range), not to which fitting method was used.
  Both methods, once given a sensible window and seed, converge to the same physically consistent answer.

In both cases, **the two fitting methods served as an independent cross-check that caught real problems
the single-method production pipeline would eventually have needed to explain anyway** — an argument for
keeping this comparison available as a diagnostic, not for switching the production default.

## Why normalization, specifically, and not some other fix

- **It is the simplest correct implementation.** The unnormalized method needs a hand-tuned `x_scale`
  array and bypasses this project's shared `NorrisFitter` class entirely. Peak-normalization needs one
  line (`y /= Y_MAX_CTS_PER_S`) and works within the existing, already-validated fitting infrastructure
  every other script in this project uses.
- **It is not an arbitrary choice of scale.** $Y_\text{max}$ is the data's own peak — the normalized
  amplitude $A\in[0,1]$ is $A_\text{physical}/Y_\text{max}$, a dimensionless fraction of peak flux, and is
  rescaled back to physical counts/s (`A_cts_per_s = A_norm * Y_MAX_CTS_PER_S`) in every results CSV this
  project produces. No physical information is discarded; it is carried alongside the normalized value.
- **Automatic alternatives don't work.** `x_scale='jac'` was tested and rejected (above) — this rules out
  "let the solver figure out the scale" as a lower-effort substitute.

## Bottom line

Peak-normalization was tested against a from-scratch alternative and found to change nothing about the
physics: same fit quality, same well-constrained parameters, for every burst in this project's sample.
The two cases where the methods disagreed were correctly attributed to weakly-constrained or degenerate
pulses in the data, not to the normalization step — confirmed by tracing each one to its root cause
independently (a raw-data check for GRB080916C, a window/seed scan for GRB140206B), not by assuming the
normalized result was right. This file and `comparison.md` are the full audit trail.
