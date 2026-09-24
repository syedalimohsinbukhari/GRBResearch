# Idea: refit the Norris pulse directly in (t_peak, t_v, r) instead of (t_s, tau1, tau2)

**Status: speculative, not implemented, not scheduled.** Written up 2026-09-24 to get outside opinions
(other AI sessions / models) before deciding whether it's worth a real trial. If nobody finds a hole in
the math or a reason it's a bad idea, it might be worth a small A/B test on one already-characterized
pulse; if not, no harm done — nothing about the production pipeline depends on this.

## 1. Context: what's being fit, and what we already observed

This project fits gamma-ray burst light curves as a sum of Norris (2005) pulses:

```
I(t) = A * exp(2*sqrt(tau1/tau2)) * exp(-tau1/(t-t_s) - (t-t_s)/tau2),   t > t_s
     = 0,                                                                 t <= t_s
```

Four parameters per pulse: amplitude `A`, onset time `t_s`, rise parameter `tau1`, decay parameter
`tau2`. Two derived quantities matter downstream, computed *after* the fit from the raw parameters:

```
t_peak = t_s + sqrt(tau1 * tau2)

t_v = (tau2/2) * sqrt[ (ln(2) + 2*sqrt(tau1/tau2))^2 - 4*tau1/tau2 ]
```

(`t_v` is attributed to Norris et al. 2005 via Bukhari et al. 2022, Adv. Space Res., eq. 10 — it's
FWHM/2 of the fitted pulse.) `t_peak` and `t_v` are the two quantities that actually feed downstream
science in this project (episode matching and the variability-timescale input to a gamma-gamma-opacity
Lorentz-factor bound, respectively). `t_s`, `tau1`, `tau2` individually feed nothing downstream except
through those two derived quantities.

**The observed phenomenon.** Refitting the same pulse from two different starting guesses (or two
different detector-summing choices, or two different fit-window widths) routinely produces wildly
different raw `t_s`/`tau1`/`tau2` while `t_peak` and `t_v` barely move. Concrete example, one pulse from
this project's GRB080916C fit, same data, two different P0 seeds landing on two different local
solutions:

| quantity | seed A | seed B | % change |
|---|---|---|---|
| `tau1` | 85.6 | 130.4 | **52%** |
| `tau2` | 0.453 | 0.393 | 13% |
| `t_peak` | 57.706 | 57.735 | **0.05%** |
| `t_v` | 1.197 | 1.185 | **1.05%** |

This isn't a one-off. The same pattern (large swings in raw `t_s`/`tau1`/`tau2`, small swings in
`t_peak`/`t_v`) has shown up independently under three different triggers in this project's fitting work:
changing the P0 seed (above), widening the fit window (a `t_s`<->`tau1` degenerate ridge that a wider
window lets the optimizer slide further along), and changing which/how many NaI detectors are summed
into the light curve (different noise realization, same underlying degeneracy).

## 2. Why this happens (the math, not just the observation)

`t_peak - t_s = sqrt(tau1 * tau2)` — the **geometric mean** of `tau1` and `tau2`, not either one
individually. Geometric mean has a property individual factors don't: it's invariant under any
rebalancing that holds the *product* fixed. If `tau1 -> 2*tau1` and `tau2 -> tau2/2`, the product
`tau1*tau2` is unchanged, so `t_peak` doesn't move *at all* — only the `tau1/tau2` ratio changed, and
that ratio governs pulse *asymmetry* (how much of the width is rise vs. decay), not *where the peak
sits*.

`t_v` isn't literally invariant the same way, but it's built from the same scale, just with a different
(still product-dominated) dependence on the ratio — hence it moves more than `t_peak` but far less than
the raw `tau1`/`tau2` do individually.

This is a textbook "sloppy parameter" situation in nonlinear model fitting: the light curve data
constrains the pulse's *timing and rough width* (the `tau1*tau2`-scale combination) far more tightly than
it constrains the *individual* rise/decay split. Two very different `(tau1, tau2)` pairs can trace nearly
the same curve near the peak if the data doesn't have the leverage to say which side is "really" rising
slowly vs. decaying slowly — which is exactly the case for a noisy tail or an overlapping neighboring
pulse.

## 3. The proposed reparametrization

If `t_peak` and `t_v` are the well-determined combinations, why not fit *directly* in terms of them,
instead of discovering after the fact that they happened to be stable? Three shape parameters
(`t_s, tau1, tau2`) describe three real degrees of freedom (location, scale, asymmetry) — `t_peak` and
`t_v` alone are only two numbers, so a third is needed to close the system. The natural third choice is
the asymmetry ratio:

```
r = tau1 / tau2      (r > 0)
```

**Derivation of the closed-form inverse map**, `(t_peak, t_v, r) -> (t_s, tau1, tau2)`:

Start from the existing formulas, substitute `tau1 = r*tau2`:

```
sqrt(tau1/tau2) = sqrt(r)
sqrt(tau1*tau2) = tau2 * sqrt(r)
```

So `t_peak - t_s = tau2 * sqrt(r)`, and the `t_v` formula becomes `t_v = tau2 * f(r)` where

```
f(r) := (1/2) * sqrt[ (ln(2) + 2*sqrt(r))^2 - 4r ]
```

**This simplifies.** Expand the bracket:

```
(ln2 + 2*sqrt(r))^2 - 4r = ln2^2 + 4*ln2*sqrt(r) + 4r - 4r = ln2^2 + 4*ln2*sqrt(r)
```

so

```
f(r) = (1/2) * sqrt( ln2^2 + 4*ln2*sqrt(r) ) = (sqrt(ln2)/2) * sqrt( ln2 + 4*sqrt(r) )
```

Nice properties, worth someone double-checking independently:
- `f(0) = (1/2)*ln2 ~ 0.3466` (finite, no singularity as `r -> 0`, i.e. as the pulse becomes a pure
  instantaneous-rise/pure-decay shape).
- `f` is smooth and strictly increasing in `r` for all `r >= 0` (the argument of the inner sqrt grows
  monotonically with `sqrt(r)`), so it's invertible on the whole physical domain — no branch/sign
  ambiguity to worry about.
- As `r -> infinity`, `f(r) ~ sqrt(ln2) * r^{1/4} -> infinity` — unbounded, so `t_v` can encode arbitrarily
  asymmetric pulses without `r` needing to blow up disproportionately.

Given `f(r)`, the full inverse map:

```
tau2 = t_v / f(r)
tau1 = r * tau2
t_s  = t_peak - sqrt(r) * tau2
```

**Implementation sketch**: a `norris_pulse_peaked(t, params)` wrapper where
`params = (A, t_peak, t_v, r)` computes `(t_s, tau1, tau2)` via the three lines above, then calls the
*existing*, unmodified `norris_pulse()` (or the project's own `norris_fit.py::norris_pulse`). No new
physics, no new pulse shape — purely a change of fitting coordinates. The forward map (needed to seed the
new parametrization from an existing fit, or to report results) is just the project's existing
`t_peak()` and `tv_value()` functions plus `r = tau1/tau2`.

## 4. Why this differs from an earlier attempt in this project that *didn't* help

This project already tried something in this spirit and it didn't work — worth being upfront about, since
it's the obvious objection. From this folder's `variability_analysis.md` (GRB131014A section): a
`(tau, xi)` reparametrization (`norris_reparam.tau_xi()`, fit via a `NorrisReparamFitter`) was tried and
"landed pulses 1-3 on different local minima entirely... and where it did land near the same point
(pulses 4, 5), the `t_s`-`tau`/`t_s`-`xi` correlations stayed at 0.94-0.99 — no better than the original
`t_s`-`tau1` correlation it was meant to fix." Recorded as a negative result, not retried since.

The key detail: that attempt reparametrized **only `tau1, tau2 -> tau, xi`** (some scale/asymmetry pair)
while **keeping `t_s` as a directly-fit parameter**. The correlations it reports afterward are explicitly
`t_s`-`tau` and `t_s`-`xi` — meaning `t_s` was never decoupled from anything; the fix targeted a
*secondary* degeneracy (how `tau1`/`tau2` individually split) while leaving the *dominant* one
(`t_s`<->`tau1`) completely untouched. That it didn't help is not evidence against the idea in this
document — it's evidence the earlier fix aimed at the wrong axis. The proposal here removes `t_s` from
the direct optimization entirely, replacing it with `t_peak` — a materially different, more aggressive
change that the earlier attempt never actually tested.

## 5. Open questions for outside review

These are the things worth another model (or a human) poking at before anyone spends real time
implementing this:

1. **Is the algebra above actually correct?** Recheck the `f(r)` simplification and the three-line
   inverse map independently, ideally by symbolic substitution back into the original `t_peak`/`t_v`
   formulas.
2. **Does isolating the sloppy direction into `r` actually help the optimizer**, or does curve_fit just
   relocate the same numerical difficulty onto a different coordinate without reducing it? The geometric
   argument (Section 2) suggests `r` should absorb most of the uncertainty cleanly, but that's a
   hypothesis, not something verified by running anything yet.
3. **Bounds behavior**: `NorrisFitter.fit_boundaries()` currently ties `t_s`'s lower bound to the fit
   window's own left edge (`x_values.min()`), which is the literal mechanism behind the
   window-widening-degeneracy findings elsewhere in this project. Under the new parametrization, what
   should `t_peak`'s bound be? A pulse's peak has to lie somewhere the data actually shows a peak, which
   is a physically tighter constraint than "anywhere the window happens to start" — plausible this avoids
   the box-constraint-pinning failure mode entirely, but untested.
4. **Does the literature already do this?** A search turned up the Kocevski-Ryde-Liang (2003) pulse
   function as a genuinely different (not a reparametrization of Norris 2005 — a separate functional form)
   pulse model built directly around a peak time, but nothing found so far that reparametrizes *this
   specific* Norris (2005) four-parameter function into `(t_peak, t_v, r)` or an equivalent. If this is
   already published somewhere, worth knowing before treating it as novel.
5. **Multi-pulse fits**: this project almost never fits a single isolated pulse — most fits are joint,
   5-7 pulses simultaneously, with neighboring pulses' tails overlapping. Does the reparametrization still
   help (or even still make sense) when the "which side is well-constrained" story is complicated by a
   neighbor's tail contributing flux under the pulse being fit? The single-pulse argument in Section 2
   doesn't obviously survive that complication unchanged.
6. **Practical test, if this clears review**: refit one already-characterized fragile pulse (e.g. this
   project's GRB080916C pulse 5, or GRB140206B's pulse 6 — both already have a precise "before" picture
   of how badly the raw parameters wander) under the new parametrization, and check (a) a neutral-seed
   reproducibility test (does `r` converge to the same point from an uninformative seed the way `t_peak`
   and `t_v` already do?), and (b) the same window-widening diagnostic already used elsewhere in this
   project, watching whether `r` now rides to its own box edge the way `t_s` used to.

## 6. Bottom line

Mathematically well-defined, cheap to implement if it clears review (the wrapper is a handful of lines
around existing functions), and targets a mechanism (Section 2) that's more specific than "reparametrize
and hope" — it explicitly removes `t_s` from direct optimization, which the project's one prior attempt
at something adjacent never actually did. Whether it *helps in practice*, on this project's actual noisy,
multi-pulse, overlapping data, is genuinely unknown and would need the small test in item 6 above to find
out. Not on any current to-do list; this document exists to get a second (and third, and fourth) opinion
before deciding whether it's worth that test at all.
