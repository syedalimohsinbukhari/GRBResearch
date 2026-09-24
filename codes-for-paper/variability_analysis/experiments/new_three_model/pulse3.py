"""Section-2 implementation (norris_3param_spec.md): 3-parameter Norris pulse (A, t_peak, t_v)
at fixed asymmetry ratio r0 = tau1/tau2.

Self-contained -- does not import norris_fit.py or anything outside this file, per the spec's
"Everything below is self-contained" instruction. Section 6 (../check_tv_consistency.py) already
confirmed the pipeline's tv_value() matches width_function() below, so no rebuild was needed there.
"""
import numpy as np

LN2 = np.log(2)
EXP_ARG_CLIP = 50.0  # exp(50) ~ 5e21, far above any real light-curve rate; prevents overflow
# when the optimizer probes extreme tau values (section 2).


def width_function(r):
    """f(r) = (sqrt(ln2)/2) * sqrt(ln2 + 4*sqrt(r)) -- t_v = tau2 * f(r). Section 1.

    f(0) = ln2/2 = 0.34657, strictly increasing on r >= 0, f(r) ~ sqrt(ln2)*r^(1/4) as r -> inf.
    """
    r = np.asarray(r, dtype=float)
    return (np.sqrt(LN2) / 2) * np.sqrt(LN2 + 4 * np.sqrt(r))


def width_function_inverse(f):
    """r(f) = [(4f^2 - ln2^2) / (4 ln2)]^2 -- closed-form inverse of width_function. Section 1."""
    f = np.asarray(f, dtype=float)
    return ((4 * f**2 - LN2**2) / (4 * LN2)) ** 2


def norris_raw(t, amplitude, t_s, tau1, tau2):
    """Norris (2005) pulse from raw parameters, numerically safe (section 2).

    arg = 2*mu - tau1/x - x/tau2 (x = t - t_s), clipped to <= EXP_ARG_CLIP before exp().
    Deliberately NOT the factored form A*exp(2*mu)*exp(-tau1/x)*exp(-x/tau2): that form's
    leading exp(2*mu) overflows for mu = sqrt(tau1/tau2) > ~355, well within the range the
    fixed-r0 wrapper below can reach (r0 up to 1e4 in the spec's profile-scan grid).
    """
    t = np.asarray(t, dtype=float)
    x = t - t_s
    safe_x = np.where(x > 0, x, 1.0)  # placeholder to avoid division by zero; masked out below
    mu = np.sqrt(tau1 / tau2)
    arg = 2 * mu - tau1 / safe_x - safe_x / tau2
    arg = np.minimum(arg, EXP_ARG_CLIP)
    value = amplitude * np.exp(arg)
    return np.where(x > 0, value, 0.0)


class Pulse3Mapping:
    """Precomputed linear-map constants (c_tau1, c_tau2, c_ts) for one fixed r0. Section 2.

    The map from (A, t_peak, t_v) to raw (A, t_s, tau1, tau2) is linear in t_v -- three
    numeric constants, no root-finding -- once r0 is fixed:

        tau1 = c_tau1 * t_v
        tau2 = c_tau2 * t_v
        t_s  = t_peak - c_ts * t_v
    """

    __slots__ = ("r0", "f0", "c_tau1", "c_tau2", "c_ts")

    def __init__(self, r0: float):
        if r0 <= 0:
            raise ValueError(f"r0 must be > 0, got {r0}")
        self.r0 = float(r0)
        self.f0 = float(width_function(r0))
        self.c_tau1 = r0 / self.f0
        self.c_tau2 = 1.0 / self.f0
        self.c_ts = np.sqrt(r0) / self.f0

    def to_raw(self, amplitude, t_peak, t_v):
        """(A, t_peak, t_v) -> raw Norris (A, t_s, tau1, tau2) for this r0. Sections 2 and 4."""
        tau1 = self.c_tau1 * t_v
        tau2 = self.c_tau2 * t_v
        t_s = t_peak - self.c_ts * t_v
        return amplitude, t_s, tau1, tau2


def make_pulse3(r0: float):
    """Build a model function f(t, amplitude, t_peak, t_v) -> flux for a single fixed r0.

    Precomputes the Pulse3Mapping once (not per call), so this is what should be handed to
    scipy.optimize.curve_fit / least_squares as the model for a given r0 -- the returned
    callable's `.mapping` attribute exposes c_tau1/c_tau2/c_ts for the section-4 forward
    conversion of the fitted (A, t_peak, t_v) back to raw Norris parameters after the fit.
    """
    mapping = Pulse3Mapping(r0)

    def _pulse3(t, amplitude, t_peak, t_v):
        _, t_s, tau1, tau2 = mapping.to_raw(amplitude, t_peak, t_v)
        return norris_raw(t, amplitude, t_s, tau1, tau2)

    _pulse3.mapping = mapping
    return _pulse3


def pulse3(t, amplitude, t_peak, t_v, r0):
    """Convenience one-off form of make_pulse3, for when r0 varies call to call (e.g. the
    section-5 chi^2 profile scan, round-trip tests) rather than being fixed across a fit."""
    return make_pulse3(r0)(t, amplitude, t_peak, t_v)


def reduced_from_raw(t_s, tau1, tau2):
    """Inverse of Pulse3Mapping.to_raw: raw Norris (t_s, tau1, tau2) -> (t_peak, t_v, r).

    Section 1 forward definitions, used for round-trip tests (section 7, T1) and for reading
    r off archived 4-param fits when choosing r0 (section 5).
    """
    t_peak = t_s + np.sqrt(tau1 * tau2)
    r = tau1 / tau2
    t_v = tau2 * width_function(r)
    return t_peak, t_v, r


def pulse4(t, amplitude, t_peak, t_v, r):
    """Optional 4-parameter wrapper (A, t_peak, t_v, r) -- comparison only (section 2, last
    paragraph). NOT the primary 3-param model; r is a free fit parameter here, not fixed."""
    tau2 = t_v / width_function(r)
    tau1 = r * tau2
    t_s = t_peak - np.sqrt(r) * tau2
    return norris_raw(t, amplitude, t_s, tau1, tau2)
