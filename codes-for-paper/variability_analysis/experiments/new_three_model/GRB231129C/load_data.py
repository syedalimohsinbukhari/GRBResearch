"""Real-data loader for the GRB231129C end-to-end run (norris_3param_reduction.md's "Known
limitations" -- the first real multi-pulse light curve through the section 1-8 pipeline).

Light-curve loading, episode bounds and the archived 6-pulse 4-param seed/fit are copied from
../../../GRB231129C/_common.py (per CLAUDE.md's "copy rather than fight sys.path" convention --
that folder is a sibling under variability_analysis/, not an ancestor already on this folder's
import chain), NOT re-derived: same energy band (10-400 keV), same three NaI detectors (n3, n6,
n7), same episode bounds and P0. Extended here to also load per-bin uncertainties (errors=True),
which GRB231129C/_common.py's own loader does not request -- this folder's pipeline requires real
sigma at every decision-driving step (see profile_scan.py's _require_sigma, review point 2).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
VARIABILITY_DIR = NEW_THREE_MODEL_DIR.parents[1]  # codes-for-paper/variability_analysis/
PROJECT_ROOT = VARIABILITY_DIR.parents[1]  # GRBResearchWork/
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from light_curves import lightcurve_data  # noqa: E402

GRB_PAPER_NAME = "GRB231129C"
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB231129779"
ENERGY_LOW, ENERGY_HIGH = 10, 400
DAT_NAI = ("n3", "n6", "n7")  # copied from ../../../GRB231129C/_common.py

# Episode boundaries, copied verbatim from ../../../GRB231129C/_common.py (sourced there from
# results.json, not eyeballed).
EPISODE_BOUNDS = {
    "EX0": (-0.192, 3.136),
    "TR1": (0.384, 3.136),
    "TR2": (3.136, 7.296),
    "EX1": (3.136, 10.048),
}

# Archived unnormalized 4-param fit results (6 pulses, full range fit_window (-138.496, 475.84)),
# already on disk from prior production work -- read directly rather than re-copied as literals.
ARCHIVED_CSV = VARIABILITY_DIR / "GRB231129C" / "norris_fit_results_GRB231129779_unnormalized.csv"

# T05/T95 (T90 start/end) and plot padding -- copied from ../../../GRB231129C/_common.py (same
# CLAUDE.md "copy rather than fight sys.path" convention as this module's own docstring notes).
# Plot-axis bounding only, same as there -- never used to restrict what gets fitted.
T05, T95 = 0.384, 7.296
PLOT_PAD_S = 5.0

# Population of archived 4-param r = tau1/tau2 values across every archived fit in this project
# (all 4 GRBs' norris_fit_results*.csv, deduplicated per unique pulse) -- the "population median of
# archived 4-param r values for this project" the spec's section 5 step 3 asks for as the
# flat-bottom action, and review point 4's other half (look at the shape of this distribution).
ALL_ARCHIVED_CSVS_GLOB = "norris_fit_results*.csv"


def load_light_curve():
    """Background-subtracted, NaI-summed count rate and its per-bin uncertainty.

    Uncertainty treatment: each detector's own within-detector, across-energy-channel
    uncertainty is already quadrature-combined by lightcurve_data() itself (see that function's
    own comment on the ~9.4x inflation bug this avoids). Summing N *independent* detectors'
    already-quadrature-combined rate/background estimates combines their uncertainties in
    quadrature again, one level up -- summing counts linearly, summing variances (not
    uncertainties) linearly, matching how this project's data_rate_error is documented to behave
    (light_curves.py's own docstring/comment).
    """
    per_detector = [
        lightcurve_data(str(LC_DIR / f"{d}.dat"), ENERGY_LOW, ENERGY_HIGH, errors=True) for d in DAT_NAI
    ]
    t = per_detector[0][0]
    for det, (t_i, _, _, _) in zip(DAT_NAI[1:], per_detector[1:]):
        assert np.array_equal(t, t_i), f"{det}'s time grid differs from {DAT_NAI[0]}'s -- cannot sum"

    rate = np.sum([r for _, r, _, _ in per_detector], axis=0)
    background = np.sum([b for _, _, b, _ in per_detector], axis=0)
    rate_var = np.sum([errs[0] ** 2 for _, _, _, errs in per_detector], axis=0)
    background_var = np.sum([errs[1] ** 2 for _, _, _, errs in per_detector], axis=0)

    y = rate - background
    sigma = np.sqrt(rate_var + background_var)
    return t, y, sigma


def load_archived_pulses() -> pd.DataFrame:
    """One row per pulse (6 rows), deduplicated from the archived CSV's per-(pulse, episode) rows
    -- (A_cts_per_s, t_s, tau1, tau2, t_peak_s) at the full-range 4-param fit, plus r = tau1/tau2."""
    df = pd.read_csv(ARCHIVED_CSV)
    df = df.drop_duplicates(subset="pulse_index").sort_values("pulse_index").reset_index(drop=True)
    df["r"] = df["tau1"] / df["tau2"]
    return df[["pulse_index", "A_cts_per_s", "t_s", "tau1", "tau2", "t_peak_s", "r"]]


def load_project_r_population() -> np.ndarray:
    """r = tau1/tau2 across every archived 4-param fit in the whole project (all norris_fit_results*.csv
    under codes-for-paper/variability_analysis/), deduplicated per unique (grb, episode-independent)
    pulse where possible. This is the "population median of archived 4-param r values for this
    project" the spec's section 5 step 3 names as the preferred flat-bottom r0 choice."""
    csv_paths = sorted(VARIABILITY_DIR.rglob(ALL_ARCHIVED_CSVS_GLOB))
    r_values = []
    for path in csv_paths:
        df = pd.read_csv(path)
        if not {"tau1", "tau2"}.issubset(df.columns):
            continue
        dedup_cols = [c for c in ("grb_name", "pulse_index") if c in df.columns] or None
        if dedup_cols:
            df = df.drop_duplicates(subset=dedup_cols)
        r_values.append((df["tau1"] / df["tau2"]).to_numpy())
    return np.concatenate(r_values)
