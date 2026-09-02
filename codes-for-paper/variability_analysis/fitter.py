"""Created on Aug 22 14:47:37 2026"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from variability_timescale.light_curves import lightcurve_data
from variability_timescale.norris_fit import NorrisFitter

ROOT = Path(__file__).resolve()
PROJECT_ROOT = ROOT.parent.parent.parent
LC_DIR = PROJECT_ROOT / "light_curves"
GRB_080916C = LC_DIR / "GRB080916009"
ENERGY_LOW, ENERGY_HIGH = 10, 400
BIN_WIDTH_S = 0.064
PAD_FRACTION = 0.20
PAD_MIN_S = 0.5
NULL_REL_ERR_THRESHOLD = 1.0  # >100% relative error on tau1 or tau2 => no_reliable_fit
T_ARR_MATCH_TOLERANCE_S = 1.0

START1 = -0.192
END1 = 15.040

dat = [f for f in os.listdir(f"{GRB_080916C}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]
dat_NaI = [i.split(".")[0] for i in dat if "n" in i]

nai_data = [lightcurve_data(f"{GRB_080916C}/{i}.dat", ENERGY_LOW, ENERGY_HIGH) for i in dat_NaI]

t1, r1, b1 = nai_data[0]
# t2, r2, b2 = nai_data[1]

mask_ = np.logical_and(t1 > START1, t1 < END1)

nf = NorrisFitter(t1[mask_], (r1 - b1)[mask_], max_iterations=10000)
# nf = NorrisFitter(t1, (r1 - b1), max_iterations=10000)
# nf.dry_run()
# plt.show()
nf.fit(
    p0=[
        (1517, -0.1, 0.5, 8),
        (500, 5.2, 0.9, 0.3),
        # (614, 8, 19, 16),
        # (250, 47, 360, 0.2), (320, 60, 8, 0.8)]
    ]
)
nf.plot_fit(show_individuals=True)

print(nf.covariance)

plt.show()
