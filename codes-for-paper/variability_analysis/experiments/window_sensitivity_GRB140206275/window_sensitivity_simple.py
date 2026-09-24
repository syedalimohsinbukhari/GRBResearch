"""GRB140206B SIMPLE-model (5-pulse) window-sensitivity check -- split out of the combined
window_sensitivity.py (2026-09-24, backup of the original at _backup_2026-09-24_pre_split/) so SIMPLE
and COMPLEX can be run independently. See window_sensitivity_common.py for the shared light-curve
loading, window definitions, and the full diagnostic rationale (unedited, carried over verbatim).

SIMPLE_P0 = that model's own narrow-window converged parameters, from
norris_fit_results_GRB140206275_simple.csv (fitter_GRB140206275_simple.py) -- copied, not re-derived.
"""
from window_sensitivity_common import run_model

SIMPLE_P0 = [
    (0.13821359231007244, -0.31068299950172557, 0.0786997853378577, 1.7179861641051446),
    (0.35672897543605114, 4.137546935144213, 5.976466850654499, 13.932789445290558),
    (0.5506659874868615, 11.076896556433432, 5.533197632966581, 1.4756390342157466),
    (0.2323517995016442, 26.73173459922519, 2.535166747797389, 3.786359998026647),
    (0.050309612777835425, -0.9518746258166284, 3439.7640412860924, 4.462751768530366),
]

if __name__ == "__main__":
    run_model("SIMPLE", SIMPLE_P0)
