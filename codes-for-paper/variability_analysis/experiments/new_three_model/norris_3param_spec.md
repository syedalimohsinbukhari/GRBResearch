# Implementation spec: 3-parameter Norris pulse fit (fixed-asymmetry reduction)

Purpose: fit GRB light-curve pulses in coordinates (A, t_peak, t_v) with the rise/decay
asymmetry ratio r = tau1/tau2 fixed at a constant r0. This removes the degenerate fitting
direction (the product-preserving tau1<->tau2 rebalance at fixed t_peak/t_v) while keeping
the two quantities used by downstream science (episode matching via t_peak; variability
timescale via t_v).

Everything below is self-contained. No step requires the original 4-parameter fitter except
as an optional comparison.

--------------------------------------------------------------------------------
1. Definitions
--------------------------------------------------------------------------------

Norris (2005) pulse, 4 raw parameters (A, t_s, tau1, tau2):

    x    = t - t_s
    I(t) = A * exp( 2*mu - tau1/x - x/tau2 )   for x > 0
         = 0                                    for x <= 0
    mu   = sqrt(tau1/tau2)

In this normalization A is exactly the peak flux: at x = sqrt(tau1*tau2) the exponent is 0.

Derived quantities (computed from raw parameters, never fitted directly in 4-param mode):

    t_peak = t_s + sqrt(tau1*tau2)                       (pulse peak time)
    r      = tau1 / tau2                                  (asymmetry ratio, r > 0)
    t_v    = (tau2/2) * sqrt( (ln2 + 2*sqrt(r))^2 - 4r )  (= FWHM/2 of the pulse)
           = tau2 * f(r)

Width function f(r) — closed form after expanding the bracket:

    f(r) = (sqrt(ln2)/2) * sqrt( ln2 + 4*sqrt(r) )

Properties (verified): f(0) = ln2/2 = 0.34657; strictly increasing on r >= 0;
f(r) ~ sqrt(ln2) * r^(1/4) as r -> inf; closed-form inverse
r(f) = [ (4 f^2 - ln2^2) / (4 ln2) ]^2.
ln2 = 0.6931471805599453.

--------------------------------------------------------------------------------
2. The 3-parameter model
--------------------------------------------------------------------------------

Fix r = r0 (constant per fit; see section 5 for how to choose it). Precompute:

    f0 = (sqrt(ln2)/2) * sqrt( ln2 + 4*sqrt(r0) )      # scalar constant

The map from fitted parameters (A, t_peak, t_v) to raw Norris parameters is LINEAR
in t_v — three numeric constants, no root-finding:

    c_tau1 = r0 / f0        # tau1 = c_tau1 * t_v
    c_tau2 = 1  / f0        # tau2 = c_tau2 * t_v
    c_ts   = sqrt(r0) / f0  # t_s  = t_peak - c_ts * t_v

    pulse3(t; A, t_peak, t_v, r0) = norris(t; A, t_peak - c_ts*t_v, c_tau1*t_v, c_tau2*t_v)

Example: r0 = 200 -> f0 = 3.138, tau1 = 63.73*t_v, tau2 = 0.3187*t_v, t_s = t_peak - 4.503*t_v.

Numerically safe evaluation of norris(): compute the exponent argument
arg = 2*mu - tau1/x - x/tau2 with x = t - t_s, return 0 for x <= 0, and evaluate
A*exp(arg) with arg clipped to a maximum of ~+50 (exp(50) ~ 5e21; any real light-curve
rate is far below this, and clipping prevents overflow when the optimizer probes
extreme tau values). Do NOT use the equivalent form
A*exp(2*mu)*exp(-tau1/x)*exp(-x/tau2): the leading exp(2*mu) overflows for mu > ~355.

Optional 4-param wrapper (for comparison only): parameters (A, t_peak, t_v, r) with
tau2 = t_v / f(r), tau1 = r * tau2, t_s = t_peak - sqrt(r)*tau2.

--------------------------------------------------------------------------------
3. Bounds and seeding
--------------------------------------------------------------------------------

Per-pulse bounds for the 3-param fit (fit-window [xmin, xmax], time bin width dt):

    A      in [0, inf)
    t_peak in [xmin, xmax]                 (a real peak must be visible in the window)
    t_v    in [2*dt, (xmax - xmin)]        (below 2 bins = unresolved; above = window scale)

Neutral seed (works for essentially any pulse):

    A0      = max(0, max(counts_per_bin)/dt - background_rate)   (or just max flux in window)
    t_peak0 = t[argmax(residual light curve)]
    t_v0    = ~2-3 dt, or half the FWHM eyeballed from the residuals

With these bounds and seeds the 3-param fit has no bad directions left: multi-start
reproducibility should be exact (5 neutral seeds -> identical solution) and is itself
a useful acceptance test.

--------------------------------------------------------------------------------
4. Forward conversion for reporting
--------------------------------------------------------------------------------

After the fit, report raw Norris parameters via the same linear map:

    tau2 = t_v / f0
    tau1 = r0 * tau2
    t_s  = t_peak - sqrt(r0) * tau2
    r    = r0 (fixed; record the flat-zone range from section 5 as its uncertainty)

Uncertainties: use the fit covariance directly for (A, t_peak, t_v). For t_peak/t_v also
propagate the r0-choice systematic: spread of best-fit t_peak/t_v across the flat zone
(section 5); add in quadrature with the statistical term.

--------------------------------------------------------------------------------
5. Choosing r0: the chi^2 profile scan (decision rule)
--------------------------------------------------------------------------------

Do this once per pulse before committing to the 3-param fit.

    1. Grid: r0 over logspace(-1, 4, 25)  (0.1 ... 10000; extend if the profile is still
       falling at the edges).
    2. For each r0: fit (A, t_peak, t_v) from the single neutral seed of section 3;
       record chi2_min(r0), t_peak(r0), t_v(r0).
    3. Plot chi2 vs r0 (log x-axis). Two shapes:

    FLAT BOTTOM  — chi2 within +10 of its minimum over a contiguous r0 range spanning
    factors of several or more. The data has no asymmetry leverage; the free r in a
    4-param fit would only slide along its degenerate valley. Action: fix r0 at the
    population median of archived 4-param r values for this project (or geometric
    center of the flat zone). Record the flat zone.

    SHARP MINIMUM — a clear chi2 minimum with curvature. The data does constrain r.
    Action: either keep that pulse 4-param (wrapper in section 2), or fix r0 at the
    profile minimum and note the smaller flat zone around it.

    4. Cost report: always compare chi2_3param(best r0) against the best multi-seed
    4-param chi2 for the same pulse. Delta-chi2 ~ a few or less = the 4th parameter was
    dead weight; Delta-chi2 >> 10 = keep r free for this pulse. This Delta-chi2 is the
    quantitative justification for the parameter reduction, pulse by pulse.

Joint multi-pulse fits: default is independent 3-param pulses (3n parameters). Optional
hierarchical variant — one r shared per episode (3n+1 parameters): use only if archived
per-pulse r values (where individually well-determined) cluster; if two well-measured
pulses in the same episode disagree on r, keep r per-pulse via the profile rule above.

--------------------------------------------------------------------------------
6. Consistency requirement (READ BEFORE IMPLEMENTING)
--------------------------------------------------------------------------------

The width function f(r) is the inverse of whatever function the pipeline uses downstream
as t_v. Before building the wrapper, VALIDATE the pipeline's existing tv_value() against

    tv_value(tau1, tau2) == (tau2/2) * sqrt( (ln2 + 2*sqrt(tau1/tau2))^2 - 4*tau1/tau2 )

on a sample of archived fits (e.g. tau1=85.6, tau2=0.453 must give t_v = 1.407). If the
existing implementation differs, the f(r) used here must be rebuilt as the exact inverse
of the implementation actually in use, or every 3-param result will be subtly biased in
(t_peak, t_v)-space.

--------------------------------------------------------------------------------
7. Implementation validation tests
--------------------------------------------------------------------------------

Run these before touching real data:

    T1  Round-trip: pick (tau1, tau2) = (85.6, 0.453); compute t_peak, t_v, r; apply
        section-2 map with r0 = r; recover tau1, tau2 to machine precision.
    T2  f(r): check f(0) = 0.34657, monotone increasing on logspace(-6, 8),
        f(1e8) ~ 83.26; check r(f) inverse round-trips.
    T3  FWHM/2: for (tau1, tau2) = (85.6, 0.453), numerically locate the two half-max
        points of I(t) (brentq on I(t) - I(t_peak)/2); (x_right - x_left)/2 must equal
        1.40710035884524.
    T4  Degeneracy check (synthetic): simulate a Norris pulse (r ~ 190) with Poisson
        noise, window starting partway up the rise. Multi-seed 4-param fits should show
        large r scatter (possibly a runaway pinned to a box edge) while t_peak is stable;
        3-param fits from 5 neutral seeds must all converge to the identical solution
        for every r0 in the grid.
    T5  Profile sanity: the chi2 profile from section 5 on the T4 simulation must be
        flat-bottomed (best 3-param chi2 within ~1 of the best multi-seed 4-param chi2).

--------------------------------------------------------------------------------
8. Post-fit audit checklist
--------------------------------------------------------------------------------

    [ ] No parameter sitting on a box edge (esp. t_v at 2*dt or t_peak at window edge)
    [ ] Multi-start reproducibility exact (3-param) — if not, something is wrong upstream
    [ ] Delta-chi2 vs 4-param recorded (section 5.4)
    [ ] r0 flat-zone range recorded; t_peak/t_v r0-systematic added in quadrature
    [ ] Window-widening diagnostic rerun: result (t_peak, t_v) should be invariant under
        window changes; if t_peak rides to a window edge, the window itself is wrong,
        not the fit
