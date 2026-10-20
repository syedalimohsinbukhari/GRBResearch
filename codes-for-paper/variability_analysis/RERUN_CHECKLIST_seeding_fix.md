# Rerun checklist — BUG-25 seeding fix (temp, delete when done)

Run each from `GRBResearchWork/`. Top-level canonical scripts were renamed/restructured
(`fitter_GRBNAME.py` + `shared_utilities.py`) after this list was first written — updated to match.

> **Status check, 2026-10-20 (documentation sweep): the "Still pending" table below is stale.** Every output CSV it lists now exists with a `seed` column equal to the script's `seed_from_name(__file__)` (and to `seed_registry.yaml` / `seed_table.csv`), i.e. all were rerun after the fix: the `GRB131014215/`, `GRB140206275/` (including both `_p0` variants) and `GRB231129779/` normalized/unnormalized CSVs (mtimes 2026-09-27 05:24-07:01, seeds 3447222311 / 1602551691 / 2091438651 / 3792274645), and the `experiments/window_sensitivity_GRB{080916009,131014215,231129779}/window_sensitivity_results.csv` (seed 493797354) and `window_sensitivity_GRB140206275/window_sensitivity_results_{simple,complex}.csv` (seed 2437793349). `experiments/window_sensitivity_GRB140206275/window_sensitivity_results.csv` (seed 12345, 2026-09-16) is the pre-split, pre-fix file and is a leftover. Only the "Deferred" row (GRB231129779 `widest_full_range`) is a genuine open item as far as this check could tell; its CSV nevertheless carries the post-fix seed, so whether that window was rerun is worth a glance. The file header says "temp, delete when done"; deleting it is the user's call.

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
