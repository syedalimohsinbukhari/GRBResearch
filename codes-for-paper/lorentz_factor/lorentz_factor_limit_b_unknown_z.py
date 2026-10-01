"""Minimum Lorentz factor (Limit B), redshift sweep for the three bursts without a spectroscopic redshift
=========================================================================================================
Limit B counterpart of ``lorentz_factor_unknown_z.py`` (see its docstring for the design).  The sweep itself is the
shared ``run_sweep()``; this script only supplies Limit B and its own seed, which differs from the Limit A script's
because ``seed_from_name(__file__)`` hashes the script name -- the same decorrelation ``lorentz_factor_limit_b.py``
has from ``lorentz_factor.py``.

Outputs:
    - lorentz_results_limit_b_unknown_z.csv  -- all computed values, at z = 1, 3, 5, 7
    - lorentz_curves_limit_b_unknown_z.csv   -- the same on a dense z grid, for the comparison figure
    - lorentz_table_limit_b_unknown_z.tex    -- rendered by generate_lorentz_table_unknown_z.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from grb_research import get_rng, seed_from_name

from lorentz_factor_unknown_z import run_sweep

SEED = seed_from_name(__file__)


def main():
    """Limit B sweep: write the CSV, then render its table."""
    out_dir = Path(__file__).parent
    table, curve = run_sweep("B", get_rng(seed=SEED), SEED)
    table.to_csv(out_dir / "lorentz_results_limit_b_unknown_z.csv", index=False)
    curve.to_csv(out_dir / "lorentz_curves_limit_b_unknown_z.csv", index=False)
    print("\nSaved: lorentz_results_limit_b_unknown_z.csv, lorentz_curves_limit_b_unknown_z.csv")

    subprocess.run([sys.executable, str(out_dir / "generate_lorentz_table_unknown_z.py"), "B"], check=True)


if __name__ == "__main__":
    main()
