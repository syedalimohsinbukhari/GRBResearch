"""Created on Aug 22 14:47:20 2026"""

import numpy as np
from pymultifit.fitters.backend import BaseFitter

from grb_research.grb_calculations import get_rng

N_SAMPLES = 10_000
SEED = 12345


def norris_pulse(x: np.ndarray, params) -> np.ndarray:
    """Single Norris pulse: A * exp(2*sqrt(tau1/tau2)) * exp(-tau1/(t-ts) - (t-ts)/tau2), t > ts, else 0."""
    amplitude, t_s, tau1, tau2 = params
    dt = x - t_s
    safe_dt = np.where(dt > 0, dt, 1.0)  # placeholder to avoid division by zero; masked out below
    value = amplitude * np.exp(2 * np.sqrt(tau1 / tau2)) * np.exp(-tau1 / safe_dt - safe_dt / tau2)
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


def tv_mc_summary(fitter: NorrisFitter, pulse_index: int, seed: int = SEED, n_samples: int = N_SAMPLES):
    """MC-propagate t_v and its 16/50/84 percentiles for one fitted pulse (1-indexed, matching pymultifit's own convention).

    Draws from the full fit covariance's 4x4 sub-block for this pulse, discards draws with
    non-positive tau1/tau2 (unphysical), and reports the fraction kept.
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
