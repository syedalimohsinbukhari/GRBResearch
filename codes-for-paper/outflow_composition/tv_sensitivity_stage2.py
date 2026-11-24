"""
Stage 2 sensitivity test: how does 1 + sigma_min move when R_0 (hence t_v) is scaled?
=====================================================================================
``zhang_peer_suppression.py`` is run at the project's shared ``N_SAMPLES`` (``grb_constants``) with its
``base_radius_cm`` scaled by 1, 0.3 and 0.1, i.e. R_0 = c * (scale * t_v) / (1 + z).  The random generator is reset to
the production script's own seed (``zhang_peer_suppression.SEED``) before every scale, so the three runs share one
random stream and only R_0 differs; the scale-1 run is therefore the production run again, and is checked against
``zhang_peer_suppression.csv`` at the end (the script raises if they disagree beyond float noise).

Why it exists (user, 2026-10-25): five episodes use their duration as t_v.  For Gamma_min that is a conservative upper
bound; for sigma_min t_v is a physical input (R_0 = c t_v / (1 + z)) and this test shows the substitution has no
consistent direction (outflow_composition.md section 6.7).

Output, ``tv_sensitivity_stage2.csv`` next to this script: the production CSV's schema (one self-describing row per
episode, redshift, L_gamma band, test band; identity, assumptions, medians with asymmetric 16-84 % errors, ``n_samples``
and ``seed`` as constant columns) plus ``t_v_scale``.  ``t_v_s`` is the unscaled timescale (the Norris median or the
episode duration); the timescale actually used is ``t_v_scale * t_v_s``.  Because the seed is the production seed, the
``t_v_scale == 1`` rows reproduce the production rows.

Run:  python tv_sensitivity_stage2.py     (the former 1000-sample version took about 8 minutes; at 10 000 samples
expect roughly ten times that, not yet timed)
"""

import sys

from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
import numpy as np
import pandas as pd
import zhang_peer_suppression as m
from grb_research import get_rng
from grb_research.grb_constants import N_SAMPLES

SCALES = (1.0, 0.3, 0.1)
OUTPUT = HERE / "tv_sensitivity_stage2.csv"
PRODUCTION = HERE / "zhang_peer_suppression.csv"
KEY = ["grb_name", "episode", "z", "l_band", "band"]
RTOL = 1e-6  # the scale-1 run must reproduce the production CSV; this allows only float-library noise

assert m.N_SAMPLES == N_SAMPLES, "zhang_peer_suppression.py must run at the shared N_SAMPLES"


def to_records(rows, scale):
    """Flatten ``collect_results()`` rows into the production CSV schema, with ``t_v_scale`` after the t_v columns."""
    records = []
    for row in rows:
        record = {}
        for key, value in row.items():
            if key == "summary":
                continue
            record[key] = value
            if key == "t_v_is_duration_fallback":
                record["t_v_scale"] = scale
        for key, column, unit in m.QUANTITIES:
            median, lower, upper = row["summary"][key]
            record[f"{column}{unit}"] = median
            record[f"{column}_err_lower{unit}"] = lower
            record[f"{column}_err_upper{unit}"] = upper
        record["n_samples"] = m.N_SAMPLES
        record["seed"] = m.SEED
        records.append(record)
    return records


def check_against_production(frame):
    """Raise if the t_v_scale == 1 rows differ from the production CSV beyond float noise."""
    if not PRODUCTION.exists():
        print(f"WARNING: {PRODUCTION.name} not found, so the scale-1 reproduction check was skipped.")
        return
    production = pd.read_csv(PRODUCTION)
    base = frame[frame.t_v_scale == 1.0]
    assert len(base) == len(production), f"row count {len(base)} against production {len(production)}"
    left = base.assign(z=base.z.round(6)).set_index(KEY).sort_index()
    right = production.assign(z=production.z.round(6)).set_index(KEY).sort_index()
    assert left.index.equals(right.index), "scale-1 rows and production rows do not cover the same episodes"
    worst = 0.0
    for _, column, unit in m.QUANTITIES:
        for suffix in ("", "_err_lower", "_err_upper"):
            name = f"{column}{suffix}{unit}"
            a, b = left[name].to_numpy(float), right[name].to_numpy(float)
            both = np.isfinite(a) & np.isfinite(b)
            if not np.array_equal(np.isfinite(a), np.isfinite(b)):
                raise SystemExit(f"scale-1 run and production disagree on which {name} entries are finite")
            worst = max(worst, float(np.max(np.abs(a[both] - b[both]) / np.maximum(np.abs(b[both]), 1e-300))))
    if worst > RTOL:
        raise SystemExit(f"scale-1 run does not reproduce {PRODUCTION.name}: worst relative difference {worst:.3g}")
    print(f"Scale-1 rows reproduce {PRODUCTION.name} (worst relative difference {worst:.2g}).")


def main():
    original = m.base_radius_cm
    frames = []
    try:
        for scale in SCALES:
            print(f"\nt_v scale = {scale}")
            m.base_radius_cm = lambda t, z, s=scale: original(t, z) * s
            m.rng = get_rng(seed=m.SEED)  # same stream for every scale, so only R_0 differs
            frames.append(pd.DataFrame(to_records(m.collect_results(), scale)))
    finally:
        m.base_radius_cm = original
    frame = pd.concat(frames, ignore_index=True)
    frame.to_csv(OUTPUT, index=False)
    print(f"Saved: {OUTPUT} ({len(frame)} rows)")
    check_against_production(frame)


if __name__ == "__main__":
    main()
