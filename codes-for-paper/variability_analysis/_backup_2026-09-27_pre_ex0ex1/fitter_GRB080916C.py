"""Created on Sep 26 03∶55∶51 2026."""

from pathlib import Path

import numpy as np

from grb_research import update_style, seed_from_name
from norris_fit import NorrisFitter
from shared_utilities import (
    initialize_work,
    grb080916C_bounds,
    grb_name_resolution,
    get_summed_lc,
    lat_details,
    save_data,
    get_pulse_summary,
    norris_plotter,
)

SEED = seed_from_name(__file__)
update_style()

GRB_NAME = "080916C"
FULL_NAME, GRB_PAPER_NAME = grb_name_resolution(grb_name=GRB_NAME)
GRB_PAPER_NAME = f"GRB{GRB_PAPER_NAME}"

ROOT = Path(__file__).resolve()
LAT_FITS = Path(__file__).parent / f"{GRB_PAPER_NAME}_lat.fits"

PROJECT_ROOT, LC_DIR, GRB_LC = initialize_work(root_path=ROOT, grb_full_name=FULL_NAME)
_, TIME, _ = grb080916C_bounds()

T05, T95, T0_MET_S = TIME

PHOTON_E_MIN_MEV = 1000
PLOT_PAD_S = 25.0

t, (y, max_y_cps), dat_NaI = get_summed_lc(GRB_LC, norm=True)
t_min, t_max = float(np.min(t)), float(np.max(t))

fitter = NorrisFitter(t, y, max_iterations=20000)

p0 = [
    (0.4, -0.7, 2.34, 0.471),
    (0.7, 0.3, 0.88, 6),
    (0.2, 5.3, 0.6, 6),
    (0.3, 1.3, 43, 14),
    (0.3, 52, 18, 0.7),
    (0.3, 61, 0.3, 3),
]

p0 = np.array(p0)
fitter.fit(p0=p0)
parameters = np.reshape(fitter.params, p0.shape)
n_pulses = parameters.shape[0]

photon_energy_MeV, photon_t_arr_s = lat_details(
    lat_fit_path=LAT_FITS, t0_met=T0_MET_S, min_photon_energy=PHOTON_E_MIN_MEV
)

photon_summary = get_pulse_summary(
    lat_t_arr=photon_t_arr_s, parameters=parameters, n_pulses=n_pulses, photon_energy=photon_energy_MeV
)

save_data(
    fitter=fitter,
    parameters=parameters,
    grb_paper_name=GRB_PAPER_NAME,
    t_range=(t_min, t_max),
    max_y_cps=max_y_cps,
    photon_summary=photon_summary,
    n_pulses=n_pulses,
    seed_number=SEED,
)

SAVE_NAME = f"norris_fitted_{GRB_PAPER_NAME}"

norris_plotter(
    fitter=fitter,
    y_max=y.max(),
    nai_dat_files=dat_NaI,
    photon_t_arr_s=photon_t_arr_s,
    photon_energy=photon_energy_MeV,
    t05=T05,
    t95=T95,
    plot_pad=PLOT_PAD_S,
    save_name=SAVE_NAME,
)
