"""
Stage 2 sensitivity test: how does 1 + sigma_min move when R_0 (hence t_v) is scaled?
=====================================================================================
A validation run, not a pipeline script: ``zhang_peer_suppression.py`` is imported with ``n_samples = 1000`` (so the
output is NOT the production CSV) and its ``base_radius_cm`` is scaled by 1, 0.3 and 0.1, with the same random stream
for every scale so that only R_0 differs. Written to ``tv_sensitivity_stage2.csv`` next to this script.

Why it exists (user, 2026-10-25): five episodes use their duration as t_v. For Gamma_min that is a conservative upper
bound; for sigma_min t_v is a physical input (R_0 = c t_v / (1 + z)) and this test shows the substitution has no
consistent direction (outflow_composition.md section 6.7).

Run:  python tv_sensitivity_stage2.py        (about 8 minutes)
"""

import sys

from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
import grb_research.grb_constants as constants

constants.N_SAMPLES = 1000
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
import zhang_peer_suppression as m
from grb_research import get_rng

assert m.N_SAMPLES == 1000
original = m.base_radius_cm
frames = []
for scale in (1.0, 0.3, 0.1):
    m.base_radius_cm = lambda t, z, s=scale: original(t, z) * s
    m.rng = get_rng(seed=m.SEED)  # same stream for every scale, so only R_0 differs
    rows = m.collect_results()
    frame = pd.DataFrame(
        [
            {
                "grb": r["grb_name"], "episode": r["episode"], "z": r["z"], "l_band": r["l_band"], "band": r["band"],
                "duration": r["t_v_is_duration_fallback"], "scale": scale,
                "edge_fraction": r["binding_at_lower_edge_fraction"],
                "one_plus_sigma_min": r["summary"]["one_plus_sigma_min"][0],
                "kT_max": r["summary"]["kt_max_keV"][0],
            }
            for r in rows
        ]
    )
    frames.append(frame)
out = pd.concat(frames)
out.to_csv(HERE / "tv_sensitivity_stage2.csv", index=False)
print(f"Saved: {HERE / 'tv_sensitivity_stage2.csv'} ({len(out)} rows)")
