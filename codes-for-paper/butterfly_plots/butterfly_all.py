"""Created on Jan 10 00:09:45 2026"""

import json

import numpy as np

from grb_research import GRBCatalog, find_project_root, update_style
from grb_research.grb_calculations import plot_all_models
from grb_research.grb_constants import short_to_long

update_style()

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"

grb_name = ["080916C", "131014A", "140206B", "231129C"]
# grb_name = ["150210A"]

with open(result_file, "r") as f:
    data = json.load(f)

grb = GRBCatalog.from_iterable(grb_list=grb_name, data=data, name_mapping=short_to_long)
grb_best = [i.get_all_best_models() for i in grb]

is_ex = [[i.interval.is_ex for i in j] for j in grb_best]
is_ex = [sum(i) for i in is_ex]

rng = np.random.default_rng(seed=42)

plot_all_models(grb_best, grb_name, 2, 2, save=True, rng=rng)
