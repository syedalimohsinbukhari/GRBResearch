"""GRB140206B COMPLEX-model (7-pulse) window-sensitivity check -- split out of the combined
window_sensitivity.py (2026-09-24, backup of the original at _backup_2026-09-24_pre_split/) so SIMPLE
and COMPLEX can be run independently. See window_sensitivity_common.py for the shared light-curve
loading, window definitions, and the full diagnostic rationale (unedited, carried over verbatim).

COMPLEX_P0 = that model's own narrow-window converged parameters, from
norris_fit_results_GRB140206275.csv (fitter_GRB140206275.py) -- copied, not re-derived.

Pulse 7 is tracked alongside pulse 5 (also_track=[7]): pulse 5 here is one of the narrow 23-28s-region
pulses (t_s~23.08, tau1~1.1), a different physical object from SIMPLE's pulse 5; COMPLEX's actual broad
pedestal is pulse 7, reported for context per the original combined script's convention.
"""
from window_sensitivity_common import run_model

COMPLEX_P0 = [
    (0.13715836098096673, -0.3082778877479421, 0.07280291527795804, 1.7779333409054536),
    (0.3485522723217019, 4.358824486438564, 5.423853787049376, 13.077114856867508),
    (0.5572846305208764, 11.202577021209974, 4.484415169647049, 1.6621123907230426),
    (0.2771673583293726, 23.06415137350121, 46.192192621651564, 1.0927353385142597),
    (0.0635676507074934, 23.081058799329575, 1.0587226502022993, 0.9800958462581607),
    (0.07553834749655068, 26.10103412529695, 68.9065407772954, 2.016129355318782),
    (0.05008650875532962, -0.6470589696141678, 3280.0479468729336, 4.642381720487216),
]

if __name__ == "__main__":
    run_model("COMPLEX", COMPLEX_P0, also_track=[7])
