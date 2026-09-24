"""Joint N-pulse extension of pulse3.py: sum N pulse3 models, each pulse keeping its own FIXED
r0 (chosen beforehand, e.g. by the single-pulse residual pipeline in ../GRB231129C/), and fit
every pulse's (A, t_peak, t_v) simultaneously against the raw light curve in one curve_fit call --
no residual-subtraction approximation, unlike ../GRB231129C/run_pipeline.py's per-pulse method.

Kept in its own folder (GRB231129C_joint/, not GRB231129C/) per user request: same burst, genuinely
different fitting method (simultaneous vs. sequential-residual), not a modification of the
already-validated single-pulse run.

Self-contained apart from importing the already-validated building blocks one directory up
(pulse3.py's make_pulse3/Pulse3Mapping, bounds_seeding.py's pulse3_bounds, profile_scan.py's
chi_square/refine_r0_zone/_require_sigma) -- no new model math here, only the N-pulse orchestration.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))

from bounds_seeding import pulse3_bounds  # noqa: E402
from pulse3 import Pulse3Mapping, make_pulse3  # noqa: E402


def make_joint_pulse3(r0_list):
    """Build a joint model f(t, A_1, t_peak_1, t_v_1, ..., A_N, t_peak_N, t_v_N) -> summed flux,
    one fixed r0 per pulse (r0_list, length N). Ready for scipy.optimize.curve_fit with a flat
    3N-parameter vector. Exposes .n_pulses and .mappings (one Pulse3Mapping per pulse, for
    section-4-style forward conversion of each pulse's fitted (A, t_peak, t_v) back to raw
    Norris params after the joint fit).
    """
    models = [make_pulse3(r0) for r0 in r0_list]
    mappings = [Pulse3Mapping(r0) for r0 in r0_list]
    n = len(r0_list)

    def joint_model(t, *params):
        if len(params) != 3 * n:
            raise ValueError(f"expected {3 * n} params (3 per pulse x {n} pulses), got {len(params)}")
        t = np.asarray(t, dtype=float)
        total = np.zeros_like(t)
        for i in range(n):
            amplitude, t_peak, t_v = params[3 * i : 3 * i + 3]
            total += models[i](t, amplitude, t_peak, t_v)
        return total

    joint_model.n_pulses = n
    joint_model.mappings = mappings
    joint_model.pulse_model = lambda i, t, amplitude, t_peak, t_v: models[i](t, amplitude, t_peak, t_v)
    return joint_model


def joint_bounds(t_window, dt, n_pulses):
    """Same section-3 bounds for every pulse (all share the same fit window), tiled N times."""
    lb1, ub1 = pulse3_bounds(t_window, dt=dt)
    return lb1 * n_pulses, ub1 * n_pulses


def flatten_seed(per_pulse_seed):
    """per_pulse_seed: list of (A, t_peak, t_v) tuples, one per pulse -> flat 3N-tuple for curve_fit."""
    flat = []
    for a, tp, tv in per_pulse_seed:
        flat.extend((a, tp, tv))
    return tuple(flat)


def unflatten_params(popt, n_pulses):
    """Flat 3N-vector -> list of N (A, t_peak, t_v) tuples."""
    return [tuple(popt[3 * i : 3 * i + 3]) for i in range(n_pulses)]


def pulse_block_cov(pcov, pulse_index):
    """3x3 covariance sub-block for pulse `pulse_index` (0-indexed) out of a joint 3N x 3N pcov."""
    start = 3 * pulse_index
    return pcov[start : start + 3, start : start + 3]


def residual_excluding(joint_model, t, y, popt, exclude_index):
    """y minus every OTHER pulse's joint-best-fit model (all pulses except exclude_index) --
    isolates pulse `exclude_index` using the joint fit's own self-consistent neighbor estimates,
    for the post-fit r0-systematic refinement (profile_scan.refine_r0_zone needs a single-pulse
    (t, y_isolated) series, not the full joint sum)."""
    params = unflatten_params(popt, joint_model.n_pulses)
    total_others = np.zeros_like(np.asarray(t, dtype=float))
    for i, (a, tp, tv) in enumerate(params):
        if i == exclude_index:
            continue
        total_others += joint_model.pulse_model(i, t, a, tp, tv)
    return np.asarray(y, dtype=float) - total_others
