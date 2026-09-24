"""Section-4 implementation (norris_3param_spec.md): forward conversion of the fitted
(A, t_peak, t_v) back to raw Norris parameters for reporting, with uncertainty propagation.

    tau2 = t_v / f0
    tau1 = r0 * tau2
    t_s  = t_peak - sqrt(r0) * tau2
    r    = r0 (fixed; the flat-zone range from section 5 is recorded as its uncertainty)

Uncertainties: fit covariance directly for (A, t_peak, t_v); for t_peak/t_v, an r0-choice
systematic is added in quadrature with the statistical term, then propagated through the
(linear) raw-parameter map.

Review finding (2026-09-25): that r0-choice systematic is NOT the coarse section-5 flat zone
(Delta-chi2<=10 over a 25-point, 5-decade grid) -- reusing it verbatim conflated "r0 has no
meaningful uncertainty" with "the grid was too coarse to resolve the basin" for sharp profiles,
where the flat zone collapses to a single grid point and reports an artificially exact-looking
sys=0.0. Instead this module always runs profile_scan.refine_r0_zone (a dedicated, finer,
Delta-chi2<=1 local scan around r0_chosen) and uses *that* zone's t_peak(r0)/t_v(r0) spread as
the systematic, for both flat and sharp profiles alike. See refine_r0_zone's own docstring for
why +10 (classification) and Delta-chi2<=1 (this systematic) are deliberately different
thresholds answering different questions.

Builds on pulse3.py (Pulse3Mapping) and profile_scan.py (R0SelectionResult, fit_single_r0,
refine_r0_zone).
"""
from dataclasses import dataclass

import numpy as np

from profile_scan import R0SelectionResult, fit_single_r0, refine_r0_zone
from pulse3 import Pulse3Mapping


@dataclass
class ReportedPulse:
    r0: float
    r0_flat_zone: tuple  # section 5's coarse (Delta-chi2<=10) flat zone -- section 4's literal
    # "record the flat-zone range from section 5 as its uncertainty" instruction for r itself.
    r0_refined_zone: tuple  # Delta-chi2<=1 local zone actually used for the t_peak/t_v systematic
    # below (see module docstring) -- tighter and more meaningful than r0_flat_zone for a sharp
    # profile, where r0_flat_zone can collapse to a single point.

    amplitude: float
    amplitude_err: float

    t_peak: float
    t_peak_err_stat: float
    t_peak_err_sys: float
    t_peak_err: float  # stat, sys combined in quadrature

    t_v: float
    t_v_err_stat: float
    t_v_err_sys: float
    t_v_err: float

    t_s: float
    t_s_err: float
    tau1: float
    tau1_err: float
    tau2: float
    tau2_err: float


def _zone_spread(scan, zone_lo: float, zone_hi: float):
    """Half the (max - min) of best-fit t_peak(r0)/t_v(r0) over [zone_lo, zone_hi] within `scan`.
    0 if fewer than 2 grid points fall in the zone (nothing to spread over)."""
    in_zone = (scan.r0_grid >= zone_lo) & (scan.r0_grid <= zone_hi)
    finite = in_zone & np.isfinite(scan.t_peak) & np.isfinite(scan.t_v)
    if np.sum(finite) < 2:
        return 0.0, 0.0
    t_peak_zone = scan.t_peak[finite]
    t_v_zone = scan.t_v[finite]
    t_peak_sys = float(np.max(t_peak_zone) - np.min(t_peak_zone)) / 2
    t_v_sys = float(np.max(t_v_zone) - np.min(t_v_zone)) / 2
    return t_peak_sys, t_v_sys


def finalize_pulse(
    t_window, y_window, selection: R0SelectionResult, sigma=None, dt=None, background_rate=0.0
) -> ReportedPulse:
    """Section 4: forward-convert + propagate uncertainty at r0 = selection.r0_chosen (the
    section-5 decision rule's output), which may or may not land on a profile-scan grid point --
    refits at r0_chosen directly if not, so the reported (A, t_peak, t_v) and its covariance are
    always for the actual r0 being used downstream, not the nearest grid point.
    """
    r0 = selection.r0_chosen
    grid_idx = np.where(np.isclose(selection.scan.r0_grid, r0))[0]
    if grid_idx.size:
        fit_result = selection.scan.results[grid_idx[0]]
    else:
        fit_result = fit_single_r0(t_window, y_window, r0, sigma=sigma, dt=dt, background_rate=background_rate)
    if not fit_result.success:
        raise RuntimeError(f"fit at chosen r0={r0!r} did not converge -- cannot finalize this pulse")

    mapping = Pulse3Mapping(r0)
    amplitude, t_peak, t_v = fit_result.amplitude, fit_result.t_peak, fit_result.t_v
    _, t_s, tau1, tau2 = mapping.to_raw(amplitude, t_peak, t_v)

    pcov = fit_result.pcov
    if pcov is not None:
        a_err, t_peak_err_stat, t_v_err_stat = np.sqrt(np.diag(pcov))
        cov_tp_tv = pcov[1, 2]
    else:
        a_err = t_peak_err_stat = t_v_err_stat = np.nan
        cov_tp_tv = 0.0

    refined_scan, refined_classification = refine_r0_zone(
        t_window, y_window, r0, sigma=sigma, dt=dt, background_rate=background_rate
    )
    r0_refined_zone = refined_classification.flat_zone
    t_peak_err_sys, t_v_err_sys = _zone_spread(refined_scan, *r0_refined_zone)
    t_peak_err = float(np.hypot(t_peak_err_stat, t_peak_err_sys))
    t_v_err = float(np.hypot(t_v_err_stat, t_v_err_sys))

    # raw-parameter propagation: tau1, tau2 are linear in t_v (exact); t_s is affine in
    # (t_peak, t_v), so its statistical variance uses the fit's t_peak-t_v covariance directly,
    # while the r0-choice systematic on t_peak/t_v is combined without a cross-term (the spec
    # does not specify a joint flat-zone systematic covariance, so this is conservative).
    tau1_err = float(mapping.c_tau1 * t_v_err)
    tau2_err = float(mapping.c_tau2 * t_v_err)
    if pcov is not None:
        var_t_s_stat = t_peak_err_stat**2 + mapping.c_ts**2 * t_v_err_stat**2 - 2 * mapping.c_ts * cov_tp_tv
        t_s_err_stat = float(np.sqrt(max(var_t_s_stat, 0.0)))
    else:
        t_s_err_stat = np.nan
    t_s_err_sys = float(np.hypot(t_peak_err_sys, mapping.c_ts * t_v_err_sys))
    t_s_err = float(np.hypot(t_s_err_stat, t_s_err_sys))

    return ReportedPulse(
        r0=r0,
        r0_flat_zone=selection.classification.flat_zone,
        r0_refined_zone=r0_refined_zone,
        amplitude=float(amplitude),
        amplitude_err=float(a_err),
        t_peak=float(t_peak),
        t_peak_err_stat=float(t_peak_err_stat),
        t_peak_err_sys=float(t_peak_err_sys),
        t_peak_err=t_peak_err,
        t_v=float(t_v),
        t_v_err_stat=float(t_v_err_stat),
        t_v_err_sys=float(t_v_err_sys),
        t_v_err=t_v_err,
        t_s=float(t_s),
        t_s_err=t_s_err,
        tau1=float(tau1),
        tau1_err=tau1_err,
        tau2=float(tau2),
        tau2_err=tau2_err,
    )
