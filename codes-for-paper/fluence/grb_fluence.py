"""Created on Apr 28 07:04:07 2026"""

import numpy as np
import pandas as pd
from black import Path

from grb_research import EpisodeTypes, FluxFluenceCalculator, find_project_root, get_rng, prepare_grbs, seed_from_name
from grb_research.grb_constants import N_SAMPLES

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SOURCE_ROOT = find_project_root()
RESULT_FILE = SOURCE_ROOT / "results.json"
GRB_LIST = ["080916C", "140206B", "131014A", "231129C"]
RANDOM_SEED = seed_from_name(__file__)

# GBM selection-criterion band, matching section-1-introduction.tex's stated
# "8 keV--40 MeV" range for the paper's "high GBM fluence" sample-selection
# criterion (fixed 2026-09-01: this previously silently used the
# FluxFluenceCalculator default of 10 keV-1 MeV, a narrower band than the
# criterion it was meant to substantiate -- see review-resolution.md
# Priority 3 item 8).
E_MIN_KEV = 8.0
E_MAX_KEV = 40_000.0  # 40 MeV
LOG_ENERGY_RANGE = (np.log10(E_MIN_KEV), np.log10(E_MAX_KEV))

cur_dir = Path(__file__).parent


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------


def episode_label(model) -> str:
    """Produce the legend label for a single model's episode."""
    if model.interval.kind in (EpisodeTypes.T90, EpisodeTypes.EX0, EpisodeTypes.EX1):
        return str(model.interval.kind)
    return f"{model.interval.kind}{model.interval.index}"


def compute_flux_fluence(
    model,
    rng: np.random.Generator,
    n_samples: int = N_SAMPLES,
    get_energy_flux: bool = False,
    log_energy_range: tuple[float, float] = LOG_ENERGY_RANGE,
) -> dict:
    """Calculate flux and fluence with uncertainties for a single model, over `log_energy_range`."""
    fc = FluxFluenceCalculator(model, rng=rng, n_samples=n_samples, log_energy_range=log_energy_range)

    flux_val, flux_lo, flux_hi = fc.calculate("flux")
    if get_energy_flux:
        fluence_val, fluence_lo, fluence_hi = fc.calculate("fluence", in_ergs=True, energy_flux=True)
    else:
        fluence_val, fluence_lo, fluence_hi = fc.calculate("fluence", in_ergs=True)

    return {
        "grb_name": None,  # Will be filled in loop
        "ep_type": episode_label(model),
        "model_name": model.name,
        "e_min_keV": 10.0 ** log_energy_range[0],
        "e_max_keV": 10.0 ** log_energy_range[1],
        "flux_ph_cm2_s": flux_val,
        "flux_err_lower_ph_cm2_s": flux_lo,
        "flux_err_upper_ph_cm2_s": flux_hi,
        "fluence_erg_cm2": fluence_val,
        "fluence_err_lower_erg_cm2": fluence_lo,
        "fluence_err_upper_erg_cm2": fluence_hi,
        "n_samples": n_samples,
        "seed": RANDOM_SEED,
    }


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

# Prepare GRB data
gc, grb_list_long, grb_objs, grb_best = prepare_grbs(GRB_LIST, RESULT_FILE, get_best=True)

rng = get_rng(seed=RANDOM_SEED)

s1, s2 = rng.spawn(2)

# Count total models for progress bar
total_models = sum(len(models) for models in grb_best)

# Collect results
results = []
results2 = []

# Zip against the short GRB_LIST (not grb_list_long, which holds the raw
# directory-style names, e.g. GRB080916009) so grb_name is the paper name
# (e.g. GRB080916C), matching the identity-column convention used by every
# other codes-for-paper CSV (see lorentz_factor.py's identical pattern).
for short_name, models in zip(GRB_LIST, grb_best):
    grb_name = f"GRB{short_name}"
    print(f"{grb_name}")

    for model in models:
        print(f"{model.name}")

        row = compute_flux_fluence(model, s1, n_samples=N_SAMPLES)
        row2 = compute_flux_fluence(model, s2, n_samples=N_SAMPLES, get_energy_flux=True)
        row["grb_name"] = grb_name
        row2["grb_name"] = grb_name
        results.append(row)
        results2.append(row2)

# Create DataFrame
flux_fluence_dataframe = pd.DataFrame(results)
flux_fluence_dataframe2 = pd.DataFrame(results2)

# Save to CSV
output_path = cur_dir / "flux_fluence.csv"
output_path2 = cur_dir / "flux_energy_flux.csv"
flux_fluence_dataframe.to_csv(output_path, index=False)
flux_fluence_dataframe2.to_csv(output_path2, index=False)
