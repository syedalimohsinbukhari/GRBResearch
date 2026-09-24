"""Section-3 implementation (norris_3param_spec.md): per-pulse bounds and neutral seeding for
the 3-parameter Norris fit (A, t_peak, t_v) at fixed r0.

Self-contained, matching pulse3.py -- does not import norris_fit.py. Bounds/seed are returned in
the (lb, ub) / p0 shape scipy.optimize.curve_fit expects, so a fit is just:

    model = make_pulse3(r0)                       # pulse3.py
    lb, ub = pulse3_bounds(t_window, dt)
    p0 = pulse3_seed(t_window, y_window, dt)
    popt, pcov = curve_fit(model, t_window, y_window, p0=p0, bounds=(lb, ub))
"""
import numpy as np


def median_dt(t) -> float:
    """Median spacing of a (possibly unsorted) time array -- the `dt` the spec's bounds/seed
    formulas need when the caller doesn't already track the light curve's bin width."""
    t = np.asarray(t, dtype=float)
    if t.size < 2:
        raise ValueError(f"need at least 2 samples to infer dt, got {t.size}")
    return float(np.median(np.diff(np.sort(t))))


def pulse3_bounds(t_window, dt: float | None = None):
    """Per-pulse bounds for the 3-param fit (section 3), for a fit window t_window = [xmin, xmax].

        A      in [0, inf)
        t_peak in [xmin, xmax]           -- a real peak must be visible in the window
        t_v    in [2*dt, xmax - xmin]    -- below 2 bins = unresolved; above = window scale

    Returns (lb, ub) as 3-tuples, ready for scipy.optimize.curve_fit(..., bounds=(lb, ub)).
    """
    t_window = np.asarray(t_window, dtype=float)
    x_min, x_max = float(t_window.min()), float(t_window.max())
    window = x_max - x_min
    if dt is None:
        dt = median_dt(t_window)

    if not (2 * dt < window):
        raise ValueError(
            f"fit window ({window!r}) is not wider than 2*dt ({2 * dt!r}) -- the t_v bound "
            f"[2*dt, window] would be empty; widen the fit window or check dt"
        )

    lb = (0.0, x_min, 2 * dt)
    ub = (np.inf, x_max, window)
    return lb, ub


def pulse3_seed(t_window, y_window, dt: float | None = None, background_rate: float = 0.0):
    """Neutral seed for the 3-param fit (section 3) -- works for essentially any pulse.

        A0      = max(0, max(y_window) - background_rate)
        t_peak0 = t[argmax(y_window)]
        t_v0    = 2.5 * dt   (spec: "~2-3 dt, or half the FWHM eyeballed from the residuals";
                               2.5*dt is the midpoint of that range)

    y_window is the light curve to seed from -- the raw data for a single-pulse fit, or the
    residual after subtracting already-fitted pulses for pulse k>1 of an N-pulse decomposition
    (spec: "t[argmax(residual light curve)]"). This pipeline's light curves are already
    background-subtracted count *rates* (see norris_fit.py's NorrisFitter docstring), so
    y_window is taken directly as a rate and background_rate defaults to 0 -- pass a nonzero
    background_rate only if y_window is still counts/dt with residual background in it.

    t_v0 is clamped into the section-3 bounds [2*dt, window] so the seed is always fit-bounds
    -valid even for an unusually narrow window; t_peak0 is architecturally always in-bounds
    since it's argmax over t_window itself.
    """
    t_window = np.asarray(t_window, dtype=float)
    y_window = np.asarray(y_window, dtype=float)
    if t_window.shape != y_window.shape:
        raise ValueError(f"t_window and y_window shape mismatch: {t_window.shape} vs {y_window.shape}")
    if dt is None:
        dt = median_dt(t_window)

    window = float(t_window.max() - t_window.min())
    idx_peak = int(np.argmax(y_window))
    t_peak0 = float(t_window[idx_peak])
    a0 = max(0.0, float(y_window[idx_peak]) - background_rate)
    t_v0 = min(max(2.5 * dt, 2 * dt), window)

    return a0, t_peak0, t_v0


def jittered_seeds(
    t_window,
    y_window,
    rng: np.random.Generator,
    n: int = 5,
    dt: float | None = None,
    background_rate: float = 0.0,
    amplitude_jitter_frac: float = 0.3,
    t_peak_jitter_dt: float = 5.0,
    t_v_jitter_frac: float = 0.5,
):
    """n seeds around the neutral seed of pulse3_seed, for multi-start reproducibility testing
    (spec section 7 T4: "3-param fits from 5 neutral seeds must all converge to the identical
    solution"). The first seed returned is always the exact, unperturbed neutral seed; the rest
    perturb each parameter relative to *its own* neutral-seed scale (A by +/- amplitude_jitter_frac,
    t_peak by +/- t_peak_jitter_dt bins, t_v by +/- t_v_jitter_frac), then clamp into bounds.

    Deliberately NOT a fraction of each parameter's full bound range: for a fit window with a lot
    of near-empty baseline around a short pulse (bound range set by the whole window, not by
    where the signal actually is), that would scatter t_peak/t_v seeds into uninformative regions
    unrelated to "a slightly different but still reasonable guess at this pulse", the thing this
    function is meant to probe.

    rng is required, not defaulted -- this project's RNG convention (src/grb_research/SEEDING.md)
    is that every draw traces back to a `seed_from_name(__file__)`-derived, caller-supplied
    generator; a silent `np.random.default_rng()` fallback here would make jittered_seeds()
    non-reproducible by default.
    """
    if dt is None:
        dt = median_dt(t_window)
    lb, ub, seed = pulse3_bounds_and_seed(t_window, y_window, dt=dt, background_rate=background_rate)
    a0, t_peak0, t_v0 = seed
    seeds = [seed]
    for _ in range(n - 1):
        a = a0 * (1 + rng.uniform(-amplitude_jitter_frac, amplitude_jitter_frac))
        t_peak = t_peak0 + rng.uniform(-t_peak_jitter_dt, t_peak_jitter_dt) * dt
        t_v = t_v0 * (1 + rng.uniform(-t_v_jitter_frac, t_v_jitter_frac))
        clamped = (float(np.clip(a, lb[0], ub[0])), float(np.clip(t_peak, lb[1], ub[1])), float(np.clip(t_v, lb[2], ub[2])))
        seeds.append(clamped)
    return seeds


def pulse3_bounds_and_seed(t_window, y_window, dt: float | None = None, background_rate: float = 0.0):
    """Convenience wrapper: bounds and seed together, both from the same inferred dt, ready to
    unpack straight into curve_fit(..., p0=p0, bounds=(lb, ub))."""
    if dt is None:
        dt = median_dt(t_window)
    lb, ub = pulse3_bounds(t_window, dt=dt)
    p0 = pulse3_seed(t_window, y_window, dt=dt, background_rate=background_rate)
    for value, lo, hi, name in zip(p0, lb, ub, ("A0", "t_peak0", "t_v0")):
        if not (lo <= value <= hi):
            raise AssertionError(f"seed {name}={value!r} falls outside bounds [{lo!r}, {hi!r}]")
    return lb, ub, p0
