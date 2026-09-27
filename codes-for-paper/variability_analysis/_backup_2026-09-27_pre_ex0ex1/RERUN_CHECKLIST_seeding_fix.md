# Rerun checklist — BUG-25 seeding fix (temp, delete when done)

Run each from `GRBResearchWork/`. Top-level canonical scripts were renamed/restructured
(`fitter_GRBNAME.py` + `shared_utilities.py`) after this list was first written — updated to match.

## Already rerun since the seeding fix (confirmed by CSV/script mtime, can skip)

| file | command |
|---|---|
| `fitter_GRB080916C.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/fitter_GRB080916C.py` |
| `fitter_GRB131014A.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/fitter_GRB131014A.py` |
| `fitter_GRB140206B.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/fitter_GRB140206B.py` |
| `fitter_GRB231129C.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/fitter_GRB231129C.py` |
| `GRB080916009/fitter_normalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB080916009/fitter_normalized.py` |
| `GRB080916009/fitter_unnormalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB080916009/fitter_unnormalized.py` |

## Still pending (CSV predates the seeding fix — not yet rerun)

| file | command |
|---|---|
| `GRB131014215/fitter_normalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB131014215/fitter_normalized.py` |
| `GRB131014215/fitter_unnormalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB131014215/fitter_unnormalized.py` |
| `GRB140206275/fitter_normalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB140206275/fitter_normalized.py` |
| `GRB140206275/fitter_unnormalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB140206275/fitter_unnormalized.py` |
| `GRB140206275/fitter_normalized_p0.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB140206275/fitter_normalized_p0.py` |
| `GRB140206275/fitter_unnormalized_p0.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB140206275/fitter_unnormalized_p0.py` |
| `GRB231129779/fitter_normalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB231129779/fitter_normalized.py` |
| `GRB231129779/fitter_unnormalized.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/GRB231129779/fitter_unnormalized.py` |
| `experiments/window_sensitivity_GRB080916009/window_sensitivity.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/experiments/window_sensitivity_GRB080916009/window_sensitivity.py` |
| `experiments/window_sensitivity_GRB131014215/window_sensitivity.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/experiments/window_sensitivity_GRB131014215/window_sensitivity.py` |
| `experiments/window_sensitivity_GRB140206275/window_sensitivity_simple.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/experiments/window_sensitivity_GRB140206275/window_sensitivity_simple.py` |
| `experiments/window_sensitivity_GRB140206275/window_sensitivity_complex.py` | `PYTHONPATH=src .venv/bin/python codes-for-paper/variability_analysis/experiments/window_sensitivity_GRB140206275/window_sensitivity_complex.py` |

## Deferred — later work (2026-09-27, user decision), not part of this checklist's remaining scope

| file | status |
|---|---|
| `experiments/window_sensitivity_GRB231129779/window_sensitivity.py` | **Not run.** BUG-26 (`norris_pulse()` exponent overflow) is fixed and verified, but this script's `widest_full_range` window will still fail to converge — pulse 3 is genuinely degenerate there (`t_s`<->`tau1` runaway), a separate, real finding, not a seeding or numerics bug. Deferred because resolving it (bound `tau1`, drop the pulse, or accept non-convergence as the finding) needs more compute than is worth spending outside the Xeon machine, which isn't reachable remotely right now. Full writeup: `variability_analysis.md`'s "Session update, 2026-09-27" section. The `narrow_-1_10`/`wide_-10_20` windows in this same script are unaffected and would converge cleanly if run — only `widest_full_range` is blocked. |
