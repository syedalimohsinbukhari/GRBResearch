"""Section-8 implementation (norris_3param_spec.md): the post-fit audit checklist, run once a
pulse's r0 has been chosen (section 5, profile_scan.py) and the pulse finalized (section 4,
report_pulse.py).

    [ ] No parameter sitting on a box edge (esp. t_v at 2*dt or t_peak at window edge)
    [ ] Multi-start reproducibility exact (3-param) -- if not, something is wrong upstream
    [ ] Delta-chi2 vs 4-param recorded (section 5.4)
    [ ] r0 flat-zone range recorded; t_peak/t_v r0-systematic added in quadrature
    [ ] Window-widening diagnostic rerun: result (t_peak, t_v) should be invariant under window
        changes; if t_peak rides to a window edge, the window itself is wrong, not the fit

Every check below returns an AuditItem; run_audit() assembles them into one AuditReport.
"""
from dataclasses import dataclass

import numpy as np

from bounds_seeding import jittered_seeds, pulse3_bounds
from profile_scan import R0SelectionResult, fit_single_r0
from report_pulse import ReportedPulse

EDGE_REL_TOL = 1e-3  # "sitting on a box edge" -- fraction of the bound's own range
REPRO_REL_TOL = 1e-3  # multi-start "exact" reproducibility tolerance (see test_multistart_t4.py)
WINDOW_WIDEN_REL_TOL = 0.05  # "invariant under window changes" tolerance


@dataclass
class AuditItem:
    name: str
    passed: bool
    detail: str
    data: object = None  # optional auxiliary payload (e.g. window-widening rows) for plotting


@dataclass
class AuditReport:
    items: list  # list[AuditItem]

    @property
    def all_passed(self) -> bool:
        return all(item.passed for item in self.items)

    def __str__(self) -> str:
        lines = [f"[{'x' if i.passed else ' '}] {i.name}: {i.detail}" for i in self.items]
        return "\n".join(lines)


def _is_pinned(value: float, lo: float, hi: float, rel_tol: float = EDGE_REL_TOL) -> bool:
    """Bug found and fixed (2026-09-25, GRB231129C_joint's full-x-span window-sensitivity test):
    the original version used a tolerance relative to the bound SPAN (hi - lo). That's fine for a
    normal-sized window, but breaks down once the span becomes large relative to the fitted value
    -- e.g. fitting over a much wider window makes t_v's upper bound huge, so rel_tol*span becomes
    an enormous absolute slop, enough to falsely flag an ordinary value nowhere near either edge as
    "pinned". Fixed to match profile_scan.py's `_is_pinned_r`: tolerance relative to the bound
    VALUE itself (np.isclose-style), with a small absolute floor so a bound of exactly 0 (e.g.
    amplitude's lower bound) still works.
    """
    near_lo = np.isclose(value, lo, rtol=rel_tol, atol=rel_tol * max(abs(lo), 1e-3))
    near_hi = np.isfinite(hi) and np.isclose(value, hi, rtol=rel_tol, atol=rel_tol * max(abs(hi), 1e-3))
    return bool(near_lo or near_hi)


def check_edge_pinning(t_window, reported: ReportedPulse, dt: float | None = None) -> AuditItem:
    """"No parameter sitting on a box edge (esp. t_v at 2*dt or t_peak at window edge)"."""
    lb, ub = pulse3_bounds(t_window, dt=dt)
    t_peak_pinned = _is_pinned(reported.t_peak, lb[1], ub[1])
    t_v_pinned = _is_pinned(reported.t_v, lb[2], ub[2])
    passed = not (t_peak_pinned or t_v_pinned)
    detail = (
        f"t_peak={reported.t_peak:.4g} (bounds [{lb[1]:.4g}, {ub[1]:.4g}]{' PINNED' if t_peak_pinned else ''}), "
        f"t_v={reported.t_v:.4g} (bounds [{lb[2]:.4g}, {ub[2]:.4g}]{' PINNED' if t_v_pinned else ''})"
    )
    return AuditItem(name="no parameter on a box edge", passed=passed, detail=detail)


def check_multistart_reproducibility(
    t_window, y_window, r0: float, rng: np.random.Generator, sigma=None, dt: float | None = None, n: int = 5
) -> AuditItem:
    """"Multi-start reproducibility exact (3-param) -- if not, something is wrong upstream".
    rng is required (see jittered_seeds' docstring -- src/grb_research/SEEDING.md)."""
    seeds = jittered_seeds(t_window, y_window, rng, n=n, dt=dt)
    solutions = []
    for p0 in seeds:
        result = fit_single_r0(t_window, y_window, r0, sigma=sigma, dt=dt, p0=p0)
        if not result.success:
            return AuditItem(
                name="multi-start reproducibility",
                passed=False,
                detail=f"a seed failed to converge at r0={r0!r} (p0={p0})",
            )
        solutions.append((result.amplitude, result.t_peak, result.t_v))

    arr = np.array(solutions)
    spread = arr.max(axis=0) - arr.min(axis=0)
    scale = np.abs(arr).mean(axis=0)
    passed = bool(spread[0] < REPRO_REL_TOL * scale[0] and spread[1] < REPRO_REL_TOL and spread[2] < REPRO_REL_TOL * scale[2])
    detail = f"{n} jittered seeds at r0={r0:.4g}: spread(A, t_peak, t_v)={spread}"
    return AuditItem(name="multi-start reproducibility", passed=passed, detail=detail)


def check_delta_chi2_recorded(selection: R0SelectionResult) -> AuditItem:
    """"Delta-chi2 vs 4-param recorded (section 5.4)"."""
    d = selection.delta_chi2
    passed = d["delta_chi2"] is not None
    detail = (
        f"delta_chi2={d['delta_chi2']:.4g} ({d['verdict']})"
        if passed
        else f"NOT recorded: {d['verdict']}"
    )
    return AuditItem(name="Delta-chi2 vs 4-param recorded", passed=passed, detail=detail)


def check_flat_zone_recorded(reported: ReportedPulse) -> AuditItem:
    """"r0 flat-zone range recorded; t_peak/t_v r0-systematic added in quadrature"."""
    lo, hi = reported.r0_flat_zone
    zone_valid = lo <= reported.r0 <= hi or np.isclose(reported.r0, lo) or np.isclose(reported.r0, hi) or lo <= hi
    quadrature_applied = reported.t_peak_err >= reported.t_peak_err_stat - 1e-12 and reported.t_v_err >= reported.t_v_err_stat - 1e-12
    passed = bool(zone_valid and quadrature_applied)
    detail = (
        f"flat_zone={reported.r0_flat_zone}, t_peak_err={reported.t_peak_err:.4g} "
        f"(stat={reported.t_peak_err_stat:.4g}, sys={reported.t_peak_err_sys:.4g}), "
        f"t_v_err={reported.t_v_err:.4g} (stat={reported.t_v_err_stat:.4g}, sys={reported.t_v_err_sys:.4g})"
    )
    return AuditItem(name="r0 flat-zone recorded + systematic in quadrature", passed=passed, detail=detail)


def check_window_widening(
    make_window_data,
    t_window_base,
    r0: float,
    dt: float,
    widen_factors=(1.5, 2.0),
    background_rate: float = 0.0,
    rel_tol: float = WINDOW_WIDEN_REL_TOL,
) -> AuditItem:
    """"Window-widening diagnostic rerun: result (t_peak, t_v) should be invariant under window
    changes; if t_peak rides to a window edge, the window itself is wrong, not the fit".

    make_window_data(t_min, t_max) -> (t_window, y_window, sigma) is the caller's data source
    (regenerated synthetic data for a wider window in tests; a slice of the same real light curve
    in production), kept data-source-agnostic on purpose.
    """
    t_min0, t_max0 = float(t_window_base.min()), float(t_window_base.max())
    t0, y0, sigma0 = make_window_data(t_min0, t_max0)
    base = fit_single_r0(t0, y0, r0, sigma=sigma0, dt=dt, background_rate=background_rate)
    if not base.success:
        return AuditItem(name="window-widening invariance", passed=False, detail="base-window fit did not converge", data=None)

    rows = [("base", t_min0, t_max0, base.t_peak, base.t_v)]
    all_ok = True
    for factor in widen_factors:
        half_extra = 0.5 * (factor - 1.0) * (t_max0 - t_min0)
        t_min, t_max = t_min0 - half_extra, t_max0 + half_extra
        t_w, y_w, sigma_w = make_window_data(t_min, t_max)
        result = fit_single_r0(t_w, y_w, r0, sigma=sigma_w, dt=dt, background_rate=background_rate)
        if not result.success:
            all_ok = False
            rows.append((f"x{factor:.1f}", t_min, t_max, np.nan, np.nan))
            continue
        rows.append((f"x{factor:.1f}", t_min, t_max, result.t_peak, result.t_v))

        t_peak_ok = abs(result.t_peak - base.t_peak) <= rel_tol * max(abs(base.t_peak), 1.0)
        t_v_ok = abs(result.t_v - base.t_v) <= rel_tol * base.t_v
        edge_ok = not _is_pinned(result.t_peak, t_min, t_max, rel_tol=1e-3)
        all_ok &= t_peak_ok and t_v_ok and edge_ok

    detail = "; ".join(f"{name} [{lo:.3g},{hi:.3g}]: t_peak={tp!r} t_v={tv!r}" for name, lo, hi, tp, tv in rows)
    return AuditItem(name="window-widening invariance", passed=bool(all_ok), detail=detail, data=rows)


def run_audit(
    t_window,
    y_window,
    selection: R0SelectionResult,
    reported: ReportedPulse,
    rng: np.random.Generator,
    sigma=None,
    dt: float | None = None,
    make_window_data=None,
    widen_factors=(1.5, 2.0),
    background_rate: float = 0.0,
) -> AuditReport:
    """Run all 5 section-8 checklist items for one finalized pulse. `make_window_data` is
    optional -- if omitted, the window-widening check is skipped (nothing to widen into) rather
    than reported as a failure, since it requires a data source beyond the fit window itself.
    `rng` is required (see jittered_seeds' docstring -- src/grb_research/SEEDING.md); pass
    `get_rng(seed=seed_from_name(__file__))` from the calling script.
    """
    items = [
        check_edge_pinning(t_window, reported, dt=dt),
        check_multistart_reproducibility(t_window, y_window, selection.r0_chosen, rng, sigma=sigma, dt=dt),
        check_delta_chi2_recorded(selection),
        check_flat_zone_recorded(reported),
    ]
    if make_window_data is not None:
        items.append(
            check_window_widening(
                make_window_data,
                np.asarray(t_window),
                selection.r0_chosen,
                dt=dt,
                widen_factors=widen_factors,
                background_rate=background_rate,
            )
        )
    return AuditReport(items=items)
