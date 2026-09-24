"""Shared synthetic-pulse dataset builders for this folder's test suite, seeded via this
project's RNG convention (src/grb_research/SEEDING.md): a deterministic seed derived from a name
string via `seed_from_name`, never a bare literal int.

Centralized here rather than duplicated per test file specifically because test_multistart_t4.py
(spec section 7 T4) and test_profile_sanity_t5.py (T5, "the chi2 profile ... on the T4
simulation") need to share the *exact same* noise realization by construction. Before this file
existed, that was done by independently typing the same literal seed integer (3, 7) in both
files -- exactly the fragile, easy-to-drift pattern SEEDING.md's overhaul eliminated everywhere
else in the codebase. Here, each dataset's seed is derived once from a fixed logical name
(`seed_from_name(f"{__file__}::<dataset>")`), so any caller, in any process, gets a bit-identical
dataset -- reproducible by construction, not by coincidence of two files agreeing on a number.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[3]  # new_three_model -> experiments -> variability_analysis -> codes-for-paper -> GRBResearchWork
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from grb_research import get_rng, seed_from_name  # noqa: E402

from pulse3 import make_pulse3  # noqa: E402

R_TRUE = 190.0
A_TRUE, T_PEAK_TRUE, T_V_TRUE = 20.0, 15.0, 1.2
DT = 0.03


def _seed(suffix: str) -> int:
    """Deterministic seed for one named dataset, derived from this module's own path plus a
    logical suffix -- so the dataset is identified by name, not by an incidental integer."""
    return seed_from_name(f"{__file__}::{suffix}")


def sharp_scenario():
    """Full pulse window, low noise -- well-resolved rise+decay, sharp chi2(r0) minimum
    expected. Used by section 5/7 (T4a: 3-param multi-start reproducibility) and section 8's
    "clean fit" reference case."""
    rng = get_rng(seed=_seed("sharp"))
    model = make_pulse3(R_TRUE)
    t_window = np.arange(0.0, 30.0, DT)
    y_clean = model(t_window, A_TRUE, T_PEAK_TRUE, T_V_TRUE)
    sigma = np.full_like(y_clean, 0.02 * A_TRUE)
    y_window = y_clean + rng.normal(scale=sigma)
    return t_window, y_window, sigma, DT


def flat_scenario():
    """Window starts partway up the rise (like moderate_truncation_scenario), but with constant
    8%-of-amplitude noise instead of sqrt(counts) -- the simpler flat-noise degenerate scenario
    shared between test_profile_scan.py's own flat_bottom_scenario and test_report_pulse.py's
    flat_case (identical construction, previously duplicated with two different literal seeds)."""
    rng = get_rng(seed=_seed("flat"))
    model = make_pulse3(R_TRUE)
    t_window = np.arange(14.5, 30.0, DT)
    y_clean = model(t_window, A_TRUE, T_PEAK_TRUE, T_V_TRUE)
    sigma = np.full_like(y_clean, 0.08 * A_TRUE)
    y_window = y_clean + rng.normal(scale=sigma)
    return t_window, y_window, sigma, DT


def moderate_truncation_scenario():
    """Window starts partway up the rise -- spec's own T4 setup; flat chi2(r0) bottom expected.
    Shared verbatim between test_multistart_t4.py (T4a) and test_profile_sanity_t5.py (T5)."""
    rng = get_rng(seed=_seed("moderate_truncation"))
    model = make_pulse3(R_TRUE)
    t_window = np.arange(14.5, 30.0, DT)
    y_clean = model(t_window, A_TRUE, T_PEAK_TRUE, T_V_TRUE)
    sigma = np.sqrt(np.maximum(y_clean, 1.0))
    y_window = y_clean + rng.normal(scale=sigma)
    return t_window, y_window, sigma, DT


def severe_truncation_scenario():
    """Window starts almost at the peak and stops well before a full decay -- more severe than
    moderate_truncation_scenario; the free-r 4-param fit runs away to a box edge here. Shared
    verbatim between test_multistart_t4.py (T4b) and test_profile_sanity_t5.py (T5)."""
    rng = get_rng(seed=_seed("severe_truncation"))
    model = make_pulse3(R_TRUE)
    t_window = np.arange(14.9, 20.0, DT)
    y_clean = model(t_window, A_TRUE, T_PEAK_TRUE, T_V_TRUE)
    sigma = np.sqrt(np.maximum(y_clean, 1.0)) * 1.5
    y_window = y_clean + rng.normal(scale=sigma)
    return t_window, y_window, sigma, DT


def audit_source(t_full_min: float = -20.0, t_full_max: float = 60.0, noise_frac: float = 0.03):
    """Wide-baseline noisy source for section 8's window-widening audit: one fixed noise
    realization over a wide baseline, returned as a make_window_data(t_min, t_max) closure that
    slices it -- so widening the window adds more of the *same* underlying data rather than
    fresh noise, the correct setup for a window-widening invariance check."""
    rng = get_rng(seed=_seed(f"audit_source:{t_full_min}:{t_full_max}:{noise_frac}"))
    model = make_pulse3(R_TRUE)
    t_full = np.arange(t_full_min, t_full_max, DT)
    y_clean_full = model(t_full, A_TRUE, T_PEAK_TRUE, T_V_TRUE)
    sigma_full = np.sqrt(np.maximum(y_clean_full, 1.0)) * noise_frac
    y_full = y_clean_full + rng.normal(scale=sigma_full)

    def make_window_data(t_min, t_max):
        mask = (t_full >= t_min) & (t_full <= t_max)
        return t_full[mask], y_full[mask], sigma_full[mask]

    return make_window_data
