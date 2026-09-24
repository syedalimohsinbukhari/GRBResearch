"""Section-5 implementation (norris_3param_spec.md): the chi^2 profile scan over r0 -- the
decision rule for choosing a per-pulse fixed asymmetry ratio, run once per pulse before
committing to the 3-param fit.

Builds on pulse3.py (sections 1-2: model) and bounds_seeding.py (section 3: bounds/seed).
Self-contained otherwise -- does not import norris_fit.py.
"""
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from bounds_seeding import pulse3_bounds_and_seed
from pulse3 import make_pulse3, pulse4

DEFAULT_R0_GRID = np.logspace(-1, 4, 25)  # section 5 step 1: 0.1 ... 10000
FLAT_ZONE_DELTA_CHI2 = 10.0  # "chi2 within +10 of its minimum" -- section 5 step 3
FLAT_ZONE_MIN_SPAN = 3.0  # "spanning factors of several or more"; 3x picked as the threshold
DEGENERACY_DELTA_CHI2 = 10.0  # "Delta-chi2 >> 10 = keep r free" -- section 5 step 4
REFINEMENT_DELTA_CHI2 = 1.0  # Delta-chi2<=1 zone used for the section-4 r0-systematic (not the
# section-5 flat/sharp classification, which stays at +10 -- see refine_r0_zone).
DEFAULT_R_BOUNDS = (1e-3, 1e6)  # free-r box for the section-5.4 comparison fit


def _require_sigma(sigma) -> None:
    """Review finding (2026-09-25): chi_square() silently falling back to unweighted SSE when
    sigma=None must never reach the decision logic -- the flat/sharp +10 threshold and the
    Delta-chi2=10 "dead weight" verdict are only meaningful in chi^2 units, and on unweighted SSE
    both scale with the arbitrary units/magnitude of the rate, making every downstream verdict
    meaningless. All real data has per-bin uncertainties available (light_curves.py's
    data_rate_error); every decision-driving entry point below requires sigma explicitly rather
    than defaulting to SSE. (chi_square() itself keeps the sigma=None fallback -- it's a generic
    utility also used for plain residual reporting, e.g. plot_diagnostics.py's fit overlays.)
    """
    if sigma is None:
        raise ValueError(
            "sigma is required here: unweighted SSE must not reach the section-5 decision logic "
            "(the flat/sharp classification and Delta-chi2 verdict are only meaningful in chi^2 "
            "units). Pass real per-bin uncertainties explicitly."
        )


def chi_square(y_data, y_model, sigma=None) -> float:
    """Sum of squared residuals, chi^2-weighted if per-point sigma is given (else plain SSE).
    Generic utility -- sigma is optional here even though it is mandatory at the decision-driving
    entry points below (fit_single_r0, profile_scan, fit_4param_all_seeds, select_r0)."""
    resid = np.asarray(y_data, dtype=float) - np.asarray(y_model, dtype=float)
    if sigma is not None:
        resid = resid / np.asarray(sigma, dtype=float)
    return float(np.sum(resid**2))


@dataclass
class R0FitResult:
    r0: float
    success: bool
    chi2: float
    amplitude: float
    t_peak: float
    t_v: float
    pcov: np.ndarray | None


def fit_single_r0(t_window, y_window, r0, sigma=None, dt=None, background_rate=0.0, p0=None) -> R0FitResult:
    """One curve_fit at fixed r0, using the section-3 bounds/seed (spec section 5 step 2:
    "from the single neutral seed of section 3") unless p0 is overridden."""
    _require_sigma(sigma)
    lb, ub, seed = pulse3_bounds_and_seed(t_window, y_window, dt=dt, background_rate=background_rate)
    if p0 is None:
        p0 = seed
    model = make_pulse3(r0)
    try:
        popt, pcov = curve_fit(
            model,
            t_window,
            y_window,
            p0=p0,
            bounds=(lb, ub),
            sigma=sigma,
            absolute_sigma=sigma is not None,
            maxfev=20000,
        )
    except RuntimeError:
        return R0FitResult(r0=r0, success=False, chi2=np.inf, amplitude=np.nan, t_peak=np.nan, t_v=np.nan, pcov=None)
    y_fit = model(t_window, *popt)
    chi2 = chi_square(y_window, y_fit, sigma=sigma)
    return R0FitResult(r0=r0, success=True, chi2=chi2, amplitude=popt[0], t_peak=popt[1], t_v=popt[2], pcov=pcov)


@dataclass
class ProfileScan:
    r0_grid: np.ndarray
    results: list  # list[R0FitResult], same order as r0_grid

    @property
    def chi2(self) -> np.ndarray:
        return np.array([r.chi2 for r in self.results])

    @property
    def t_peak(self) -> np.ndarray:
        return np.array([r.t_peak for r in self.results])

    @property
    def t_v(self) -> np.ndarray:
        return np.array([r.t_v for r in self.results])


def profile_scan(t_window, y_window, r0_grid=None, sigma=None, dt=None, background_rate=0.0) -> ProfileScan:
    """Section 5 steps 1-2: fit (A, t_peak, t_v) at each r0 in the grid, recording chi2_min(r0),
    t_peak(r0), t_v(r0)."""
    _require_sigma(sigma)
    if r0_grid is None:
        r0_grid = DEFAULT_R0_GRID
    r0_grid = np.asarray(r0_grid, dtype=float)
    results = [
        fit_single_r0(t_window, y_window, r0, sigma=sigma, dt=dt, background_rate=background_rate)
        for r0 in r0_grid
    ]
    return ProfileScan(r0_grid=r0_grid, results=results)


@dataclass
class ProfileClassification:
    shape: str  # "flat" or "sharp"
    r0_at_min: float
    chi2_min: float
    flat_zone: tuple  # (lo, hi) in r0 -- contiguous run around the minimum within flat_tol
    flat_zone_span: float  # flat_zone[1] / flat_zone[0]
    r0_recommended: float  # geometric center of flat_zone if flat, else r0_at_min


def classify_profile(
    scan: ProfileScan, flat_tol: float = FLAT_ZONE_DELTA_CHI2, min_span: float = FLAT_ZONE_MIN_SPAN
) -> ProfileClassification:
    """Section 5 step 3: FLAT BOTTOM (chi2 within +flat_tol of its minimum over a contiguous r0
    range spanning >= min_span) vs. SHARP MINIMUM (a clear minimum with curvature)."""
    chi2 = scan.chi2
    r0_grid = scan.r0_grid
    valid = np.isfinite(chi2)
    if not np.any(valid):
        raise RuntimeError("profile scan: no r0 in the grid produced a successful fit")

    chi2_min = float(np.min(chi2[valid]))
    idx_min = int(np.argmin(np.where(valid, chi2, np.inf)))
    r0_at_min = float(r0_grid[idx_min])

    within = valid & (chi2 <= chi2_min + flat_tol)
    lo, hi = idx_min, idx_min
    while lo - 1 >= 0 and within[lo - 1]:
        lo -= 1
    while hi + 1 < len(within) and within[hi + 1]:
        hi += 1
    flat_zone = (float(r0_grid[lo]), float(r0_grid[hi]))
    span = flat_zone[1] / flat_zone[0] if flat_zone[0] > 0 else np.inf

    shape = "flat" if span >= min_span else "sharp"
    r0_recommended = float(np.sqrt(flat_zone[0] * flat_zone[1])) if shape == "flat" else r0_at_min

    return ProfileClassification(
        shape=shape,
        r0_at_min=r0_at_min,
        chi2_min=chi2_min,
        flat_zone=flat_zone,
        flat_zone_span=span,
        r0_recommended=r0_recommended,
    )


def refine_r0_zone(
    t_window,
    y_window,
    r0_center: float,
    sigma=None,
    dt=None,
    background_rate: float = 0.0,
    n_points: int = 10,
    span_factor: float = 3.0,
    delta_chi2_tol: float = REFINEMENT_DELTA_CHI2,
):
    """Review finding (2026-09-25): the coarse 25-point section-5 grid (log-spaced over 5
    decades) is too coarse to resolve a *sharp* profile's own curvature-based r0 uncertainty --
    classify_profile's flat_tol=10 zone collapses to a single grid point there, reporting
    t_peak_err_sys = t_v_err_sys = 0.0 exactly, which conflates "r0 has no meaningful
    uncertainty" with "the grid was too coarse to resolve the basin". A sharp minimum still has
    genuine, non-zero r0 uncertainty from the chi^2 curvature near it.

    Fix: always run a dedicated, finer local scan around r0_center (n_points log-spaced over
    [r0_center/span_factor, r0_center*span_factor]) and take its own Delta-chi2<=delta_chi2_tol
    zone (1.0 by default, a genuine ~1-sigma-flavored interval -- NOT the section-5
    classification's +10, which stays as-is; that threshold is deliberately wide because it
    answers a different question, "does this pulse's data meaningfully constrain r at all", not
    "what is r0's local uncertainty"). Used by report_pulse.py for the section-4 r0-choice
    systematic, for both flat AND sharp profiles alike -- not only where section 5's coarse
    classification happened to already resolve multiple grid points.

    Returns (refined_scan, refined_classification); refined_classification.flat_zone is the
    Delta-chi2<=1 interval to use for the systematic, regardless of its own .shape label (which
    is a byproduct of reusing classify_profile, not itself meaningful at this tolerance).
    """
    r0_grid = np.geomspace(r0_center / span_factor, r0_center * span_factor, n_points)
    refined_scan = profile_scan(t_window, y_window, r0_grid=r0_grid, sigma=sigma, dt=dt, background_rate=background_rate)
    refined_classification = classify_profile(refined_scan, flat_tol=delta_chi2_tol, min_span=np.inf)
    return refined_scan, refined_classification


def fit_4param_all_seeds(
    t_window, y_window, r0_seeds=None, sigma=None, dt=None, background_rate=0.0, r_bounds=DEFAULT_R_BOUNDS
) -> list:
    """Multi-seed 4-param (A, t_peak, t_v, r free) fit via pulse4, returning every seed's
    converged result (not just the best) -- used by section 7 T4 to check the free-r fit's
    scatter across seeds ("large r scatter ... while t_peak is stable"), the contrast case for
    the 3-param fit's exact reproducibility.
    """
    _require_sigma(sigma)
    if r0_seeds is None:
        # Spans within a decade of the full r_bounds range by default (not just 1-1000 against
        # bounds that reach 1e6): delta_chi2_report's pinned-comparison flag (fix for review
        # finding #3) can only ever fire if at least one seed actually lands near the bound, so
        # the default range needs to reach there rather than requiring the caller to widen
        # r0_seeds by hand before a real box-edge runaway would even be detectable.
        lo, hi = r_bounds
        r0_seeds = np.geomspace(max(lo * 10, 1e-2), min(hi / 10, 1e5), 7)
    lb3, ub3, seed3 = pulse3_bounds_and_seed(t_window, y_window, dt=dt, background_rate=background_rate)
    lb = (*lb3, r_bounds[0])
    ub = (*ub3, r_bounds[1])

    results = []
    for r_seed in r0_seeds:
        p0 = (*seed3, r_seed)
        try:
            popt, pcov = curve_fit(
                pulse4,
                t_window,
                y_window,
                p0=p0,
                bounds=(lb, ub),
                sigma=sigma,
                absolute_sigma=sigma is not None,
                maxfev=20000,
            )
        except RuntimeError:
            continue
        y_fit = pulse4(t_window, *popt)
        chi2 = chi_square(y_window, y_fit, sigma=sigma)
        results.append(
            {
                "r_seed": r_seed,
                "chi2": chi2,
                "amplitude": popt[0],
                "t_peak": popt[1],
                "t_v": popt[2],
                "r": popt[3],
                "pcov": pcov,
            }
        )
    return results


def fit_4param_multistart(
    t_window, y_window, r0_seeds=None, sigma=None, dt=None, background_rate=0.0, r_bounds=DEFAULT_R_BOUNDS
):
    """Best-chi2 seed from fit_4param_all_seeds, for the section-5.4 cost report ("always compare
    chi2_3param(best r0) against the best multi-seed 4-param chi2"). Returns None if every seed
    fails to converge."""
    results = fit_4param_all_seeds(
        t_window, y_window, r0_seeds=r0_seeds, sigma=sigma, dt=dt, background_rate=background_rate, r_bounds=r_bounds
    )
    if not results:
        return None
    return min(results, key=lambda r: r["chi2"])


def _is_pinned_r(value: float, r_bounds, rtol: float = 1e-6) -> bool:
    lo, hi = r_bounds
    return bool(np.isclose(value, lo, rtol=rtol) or np.isclose(value, hi, rtol=rtol))


def delta_chi2_report(classification: ProfileClassification, fit4: dict | None, r_bounds=DEFAULT_R_BOUNDS) -> dict:
    """Section 5 step 4: Delta-chi2 = chi2_3param(best r0) - chi2_4param(best multi-seed).
    ~ a few or less = the 4th parameter was dead weight; >> 10 = keep r free for this pulse.

    Review finding (2026-09-25): if the best-of-seeds 4-param comparison fit itself ran away and
    pinned to an r-bound (T4b's box-edge case: all 9 seeds converged to r=1e6), the verdict is
    comparing the 3-param fit against a comparison fit that is itself untrustworthy -- "dead
    weight" would be misleadingly reassuring (the 4-param fit isn't meaningfully constraining r
    either) and "keep r free" would be misleadingly alarming (there's no reliable free-r estimate
    to keep). Both are replaced with an explicit "inconclusive" verdict when this happens, and
    r_pinned=True is recorded so callers (e.g. audit_checklist.py's Delta-chi2-recorded check)
    surface it rather than silently trusting either normal verdict.
    """
    chi2_3param_best = classification.chi2_min
    if fit4 is None:
        return {
            "chi2_3param": chi2_3param_best,
            "chi2_4param": None,
            "delta_chi2": None,
            "r_pinned": None,
            "verdict": "4-param multi-start fit did not converge from any seed",
        }
    r_pinned = _is_pinned_r(fit4["r"], r_bounds)
    delta = chi2_3param_best - fit4["chi2"]
    if r_pinned:
        verdict = (
            f"inconclusive: 4-param comparison fit ran away and pinned at r={fit4['r']:.4g} "
            f"(bounds {r_bounds}) -- Delta-chi2 does not reflect a trustworthy free-r estimate"
        )
    else:
        verdict = "4th parameter is dead weight" if delta <= DEGENERACY_DELTA_CHI2 else "keep r free for this pulse"
    return {
        "chi2_3param": chi2_3param_best,
        "chi2_4param": fit4["chi2"],
        "delta_chi2": delta,
        "r_pinned": r_pinned,
        "verdict": verdict,
    }


@dataclass
class R0SelectionResult:
    scan: ProfileScan
    classification: ProfileClassification
    fit4: dict | None
    delta_chi2: dict
    r0_chosen: float
    r0_choice_reason: str


def select_r0(
    t_window,
    y_window,
    r0_grid=None,
    sigma=None,
    dt=None,
    background_rate=0.0,
    archived_r_median=None,
    r0_seeds_4param=None,
    r_bounds=DEFAULT_R_BOUNDS,
) -> R0SelectionResult:
    """Run the full section-5 decision rule for one pulse: profile scan, flat/sharp
    classification, multi-seed 4-param cost comparison, and the r0 choice itself.

    archived_r_median: population median of archived 4-param r values for this project (spec
    section 5 step 3's preferred flat-bottom action); if None, falls back to the flat zone's
    geometric center. Review finding (2026-09-25): taken on faith previously -- if it falls
    outside this pulse's own measured flat zone (possible, since the archived value comes from
    4-param fits that may themselves have been degenerate), that would silently fix r0 at a point
    this pulse's own profile says is measurably worse. Now checked: an out-of-zone median warns
    and falls back to the flat-zone geometric center instead of being trusted blindly.
    """
    scan = profile_scan(t_window, y_window, r0_grid=r0_grid, sigma=sigma, dt=dt, background_rate=background_rate)
    classification = classify_profile(scan)
    fit4 = fit_4param_multistart(
        t_window, y_window, r0_seeds=r0_seeds_4param, sigma=sigma, dt=dt, background_rate=background_rate, r_bounds=r_bounds
    )
    delta = delta_chi2_report(classification, fit4, r_bounds=r_bounds)

    if classification.shape == "flat":
        zone_lo, zone_hi = classification.flat_zone
        if archived_r_median is not None and zone_lo <= archived_r_median <= zone_hi:
            r0_chosen = float(archived_r_median)
            reason = (
                f"flat-bottom profile (flat zone {classification.flat_zone}, span "
                f"x{classification.flat_zone_span:.1f}); fixed to archived population median r={r0_chosen!r}"
            )
        else:
            if archived_r_median is not None:
                warnings.warn(
                    f"archived_r_median={archived_r_median!r} falls outside this pulse's own flat zone "
                    f"{classification.flat_zone} -- this pulse's profile says that value is measurably "
                    f"worse than the flat-zone center; falling back to the flat-zone geometric center "
                    f"instead of trusting the archived value.",
                    stacklevel=2,
                )
            r0_chosen = classification.r0_recommended
            reason = (
                f"flat-bottom profile (flat zone {classification.flat_zone}, span "
                f"x{classification.flat_zone_span:.1f}); "
                + (
                    "no archived median given, "
                    if archived_r_median is None
                    else f"archived median r={archived_r_median!r} fell outside the flat zone, "
                )
                + f"fixed to flat-zone geometric center r0={r0_chosen!r}"
            )
    else:
        r0_chosen = classification.r0_at_min
        reason = (
            f"sharp minimum at r0={r0_chosen!r} (flat zone {classification.flat_zone}, span only "
            f"x{classification.flat_zone_span:.1f})"
        )

    return R0SelectionResult(
        scan=scan, classification=classification, fit4=fit4, delta_chi2=delta, r0_chosen=r0_chosen, r0_choice_reason=reason
    )
