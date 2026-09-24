"""Raw-data-only loader for GRB131014215's blind discovery test (PROGRESS.md: "blank slate test
for the pipeline" -- user explicitly withheld any already-fitted parameters for this burst; this
file must not, and does not, read any archived P0, fitted (A, t_s, tau1, tau2), or episode-bound
file for this burst -- only the raw RMFIT .dat products).

Detector list (n9, na, nb) and energy band (10-400 keV) are read directly from what .dat files
exist on disk / this project's universal energy-band convention (used identically for all 4
archived GRBs' light curves) -- data-selection metadata, not a fitted parameter.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NEW_THREE_MODEL_DIR = HERE.parent
VARIABILITY_DIR = NEW_THREE_MODEL_DIR.parents[1]  # codes-for-paper/variability_analysis/
PROJECT_ROOT = VARIABILITY_DIR.parents[1]  # GRBResearchWork/
sys.path.insert(0, str(NEW_THREE_MODEL_DIR))
sys.path.insert(0, str(VARIABILITY_DIR))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from light_curves import lightcurve_data  # noqa: E402

GRB_PAPER_NAME = "GRB131014A"
LC_DIR = PROJECT_ROOT / "light_curves" / "GRB131014215"
ENERGY_LOW, ENERGY_HIGH = 10, 400
DAT_NAI = tuple(sorted(f.stem for f in LC_DIR.glob("*.dat") if f.stem.startswith("n")))


def load_light_curve():
    """Background-subtracted, NaI-summed count rate and its per-bin uncertainty -- same
    quadrature-combination treatment as GRB231129C/load_data.py's load_light_curve (independent
    detectors' variances add linearly, not their uncertainties)."""
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
