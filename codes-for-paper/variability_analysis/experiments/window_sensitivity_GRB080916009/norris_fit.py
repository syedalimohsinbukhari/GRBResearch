"""Norris (2005) pulse fitter, built on pymultifit's BaseFitter N-fit mechanism.

t_v is defined per Bukhari et al. (2022, Adv. Space Res.) eq. (10):

    t_v = (tau2/2) * sqrt[(ln(2) + 2*sqrt(tau1/tau2))^2 - 4*tau1/tau2]

attributed there to Norris et al. (2005), and described as "the half width of the pulse at
half maximum" -- verified numerically (see variability_timescale.md) to equal FWHM/2 of the
fitted pulse, which for an asymmetric pulse (tau1 != tau2) is neither the rise-side nor the
decay-side half-width alone, but their symmetrized combination.
"""

import numpy as np
from pymultifit.fitters.backend import BaseFitter

EXP_ARG_CLIP = 50.0  # exp(50) ~ 5e21, far above any real light-curve rate; prevents
# overflow when the optimizer probes extreme tau values during a fit (BUG-26). Same
# value and same fix as experiments/new_three_model/pulse3.py's norris_raw().

from grb_research.grb_calculations import get_rng
from grb_research.grb_constants import N_SAMPLES


def norris_pulse(x: np.ndarray, params) -> np.ndarray:
    """Single Norris pulse: A * exp(2*sqrt(tau1/tau2) - tau1/(t-ts) - (t-ts)/tau2), t > ts, else 0.

    Combines the exponent into one argument (clipped to EXP_ARG_CLIP before exp()), not the
    factored A*exp(2*mu)*exp(-tau1/x)*exp(-x/tau2) -- that form's leading exp(2*mu) overflows for
    mu=sqrt(tau1/tau2) > ~355, which a fit's optimizer can reach probing extreme tau (BUG-26).
    Bit-identical to the old factored form for any mu that doesn't already overflow it (verified
    to floating-point noise, ~1e-14, on ordinary parameter values) -- this only changes behavior
    in the overflow regime, from NaN/inf to a large-but-finite, correctly-decaying value.
    """
    amplitude, t_s, tau1, tau2 = params
    dt = x - t_s
    safe_dt = np.where(dt > 0, dt, 1.0)  # placeholder to avoid division by zero; masked out below
    mu = np.sqrt(tau1 / tau2)
    arg = 2 * mu - tau1 / safe_dt - safe_dt / tau2
    arg = np.minimum(arg, EXP_ARG_CLIP)
    value = amplitude * np.exp(arg)
    return np.where(dt > 0, value, 0.0)


def t_peak(t_s: float, tau1: float, tau2: float) -> float:
    """Pulse peak time (absolute), t_s + sqrt(tau1*tau2)."""
    return t_s + np.sqrt(tau1 * tau2)


def tv_value(tau1, tau2):
    """t_v per Bukhari et al. (2022) eq. (10) -- see module docstring for sourcing."""
    ratio = tau1 / tau2
    return (tau2 / 2) * np.sqrt((np.log(2) + 2 * np.sqrt(ratio)) ** 2 - 4 * ratio)


class NorrisFitter(BaseFitter):
    """Fits a sum of N Norris (2005) pulses to background-subtracted count-rate data."""

    def __init__(self, x_values, y_values, max_iterations: int = 5000):
        super().__init__(x_values=x_values, y_values=y_values, max_iterations=max_iterations)
        self.n_par = 4  # amplitude, t_s, tau1, tau2

    def fit_boundaries(self):
        x_min, x_max = self.x_values.min(), self.x_values.max()
        lb = (0.0, x_min, 1e-4, 1e-4)
        ub = (np.inf, x_max, np.inf, np.inf)
        return lb, ub

    @staticmethod
    def fitter(x, params) -> np.ndarray:
        return norris_pulse(x, params)


def tv_mc_summary(fitter: NorrisFitter, pulse_index: int, seed: int, n_samples: int = N_SAMPLES):
    """MC-propagate t_v and its 16/50/84 percentiles for one fitted pulse (1-indexed, matching pymultifit's own convention).

    Draws from the full fit covariance's 4x4 sub-block for this pulse, discards draws with
    non-positive tau1/tau2 (unphysical), and reports the fraction kept.

    `seed` is required, not defaulted, per SEEDING.md's "never a bare literal" convention -- every
    caller must supply its own seed_from_name(__file__)-derived value rather than silently falling
    back to a shared constant (this module used to default to a hardcoded SEED=12345 here; BUG-25).
    """
    n_par = fitter.n_par
    start = (pulse_index - 1) * n_par
    mean = fitter.params[start : start + n_par]
    cov = fitter.covariance[start : start + n_par, start : start + n_par]

    rng = get_rng(seed=seed)
    draws = rng.multivariate_normal(mean=mean, cov=cov, size=n_samples)

    tau1_draws, tau2_draws = draws[:, 2], draws[:, 3]
    valid = (tau1_draws > 0) & (tau2_draws > 0)
    kept_fraction = valid.mean()

    tv_draws = tv_value(tau1_draws[valid], tau2_draws[valid])
    p16, p50, p84 = np.percentile(tv_draws, [16.0, 50.0, 84.0])

    return {
        "t_v_s": p50,
        "t_v_err_lower_s": p50 - p16,
        "t_v_err_upper_s": p84 - p50,
        "n_samples": n_samples,
        "seed": seed,
        "kept_fraction": kept_fraction,
    }
