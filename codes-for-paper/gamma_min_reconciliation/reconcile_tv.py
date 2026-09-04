"""Reconcile this paper's conservative Gamma_min against a literature value from finer time binning.

This paper adopts each episode's full duration as t_v -- a deliberately conservative upper bound on
the true variability timescale (see lorentz_factor.py, lorentz_factor.md), which makes Gamma_min a
conservative (weaker) lower limit. Where a literature source reports a higher Gamma_min for the same
burst from finer time bins, the natural question is whether the two are consistent once t_v is
accounted for. This script answers that with the Lithwick & Sari (2001) scaling relation used
throughout this paper,

    Gamma_min ~ t_v^{-1/(2*alpha+2)}                          (eq:gamma_min, this paper; sec-5)

by computing the t_v the literature source's finer bin would need, for its stated alpha and
Gamma_min, given our own (t_v, alpha, Gamma_min) at the episode-duration t_v.

This is a scaling *check*, not a full reconciliation: it holds f_1 (the photon flux entering
tau_hat, eq:tau_hat) fixed at our own episode-duration value, when in reality f_1 also varies
between time bins and is not itself reconciled here. See the caveat this produced in the paper text
(GRBResearchPaper/tex_files/section-6-discussion.tex, GRB080916C paragraph).

Our own (t_v, alpha, Gamma_min) are read directly from codes-for-paper/lorentz_factor/lorentz_results.csv
-- the same values already in the paper's own Gamma_min table -- rather than re-derived or hardcoded here,
so this script can never silently drift out of sync with that table. Only the literature comparison values
(Gamma_min, and which paper/time-bin they come from) are not derivable from this project's own data and
are therefore the one thing hardcoded below, in RECONCILIATIONS.

Usage: add a case to RECONCILIATIONS below for any future literature Gamma_min comparison, then run
    .venv/bin/python codes-for-paper/gamma_min_reconciliation/reconcile_tv.py
"""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

LORENTZ_RESULTS_CSV = Path(__file__).parent.parent / "lorentz_factor" / "lorentz_results.csv"


@dataclass
class Reconciliation:
    grb: str
    episode: str
    gamma_min_literature: float
    literature_citation: str
    literature_note: str


# Literature Gamma_min values this project's own conservative estimate is being checked against.
# Not derivable from this repo's data -- see module docstring.
RECONCILIATIONS = [
    Reconciliation(
        grb="GRB080916C",
        episode="T90",
        gamma_min_literature=887.0,
        literature_citation="Abdo2009FermiObservations080916C",
        literature_note="finer time bin (this paper's cited higher of their two limits)",
    ),
    Reconciliation(
        grb="GRB080916C",
        episode="T90",
        gamma_min_literature=608.0,
        literature_citation="Abdo2009FermiObservations080916C",
        literature_note="finer time bin (this paper's cited lower of their two limits)",
    ),
]


def required_t_v(gamma_min_ours, t_v_ours_s, alpha_ls, gamma_min_literature):
    """t_v the literature source's Gamma_min would require, via Gamma_min ~ t_v^{-1/(2*alpha+2)}."""
    exponent = -1.0 / (2.0 * alpha_ls + 2.0)
    ratio = gamma_min_literature / gamma_min_ours
    return t_v_ours_s * ratio ** (1.0 / exponent)


def main():
    results = pd.read_csv(LORENTZ_RESULTS_CSV)

    for r in RECONCILIATIONS:
        row = results[(results.GRB == r.grb) & (results.episode == r.episode)]
        if row.empty:
            raise KeyError(f"{r.grb} {r.episode} not found in {LORENTZ_RESULTS_CSV} -- has it been renamed?")
        row = row.iloc[0]

        t_v_needed = required_t_v(row.Gamma_min, row.t_v_s, row.alpha_LS, r.gamma_min_literature)

        print(f"{r.grb} {r.episode} vs {r.literature_citation} ({r.literature_note})")
        print(f"  ours:       Gamma_min = {row.Gamma_min:.1f} at t_v = {row.t_v_s:.3f} s, alpha = {row.alpha_LS:.3f}")
        print(f"  literature: Gamma_min = {r.gamma_min_literature:.1f}")
        print(f"  requires t_v ~= {t_v_needed:.3f} s")
        print()


if __name__ == "__main__":
    main()
