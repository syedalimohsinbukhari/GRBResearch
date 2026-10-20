# GRBResearchWork

Python analysis code and data for a paper on photospheric (blackbody) emission in four *Fermi* GRBs: GRB080916C, GRB131014A, GRB140206B and GRB231129C. The LaTeX manuscript lives in the sibling repo `../GRBResearchPaper`; this repo produces the numbers, tables and figures it includes.

Four more burst directories are present (`GRB110721200`, `GRB110731465`, `GRB150210935`, `GRB190114873`). They are **not** part of this paper's sample; see `md_file_refs/CLAUDE.md` before using them.

## Running scripts

Run everything from the repo root with the project virtualenv. The `grb_research` package under `src/` is not installed into `.venv`, so put it on the path with `PYTHONPATH=src`:

```bash
PYTHONPATH=src .venv/bin/python -c "import grb_research"      # sanity check
PYTHONPATH=src .venv/bin/python codes-for-paper/runner_all.py --dry-run   # list the pipeline, run nothing
PYTHONPATH=src .venv/bin/python codes-for-paper/runner_all.py             # rebuild every registered script
PYTHONPATH=src .venv/bin/python codes-for-paper/sync_paper_assets.py      # copy registered figures/tables into ../GRBResearchPaper
```

`runner_all.py` sets `PYTHONPATH` itself and runs each script with its own folder as the working directory; do the same when running a single script by hand (`cd codes-for-paper/<topic>` and use `PYTHONPATH=../../src ../../.venv/bin/python <script>.py`), because several scripts write outputs relative to the current directory. Dependencies are listed in `pyproject.toml` / `uv.lock` (Python 3.10).

## Folder map

| Path | Contents |
|---|---|
| `src/grb_research/` | the `grb_research` package: fit-result parsing, time intervals, flux/fluence Monte Carlo, styling and constants (`SEEDING.md` and `CRITERIA.md` document the RNG-seeding and model-selection conventions) |
| `codes-for-paper/` | one folder per analysis topic (`fluence/`, `bb_fraction/`, `photospheric_radius/`, `lorentz_factor/`, `amati_relationship/`, `variability_analysis/`, ...), each with its scripts, CSVs, figures, generated `.tex` tables and a method note (`<topic>.md`) |
| `codes-for-paper/*_registry.yaml` | `runner_registry.yaml` (scripts `runner_all.py` runs), `figure_registry.yaml` / `table_registry.yaml` (what `sync_paper_assets.py` copies into the paper repo), `seed_registry.yaml` (scripts in the paper's RNG-seed table) |
| `results.json` | fitted-model database (all bursts, all intervals and models) that most scripts read |
| `GRB080916009/`, `GRB131014215/`, `GRB140206275/`, `GRB231129779/` (and the other four) | raw per-interval fit directories; `GRB131014215GBM/` and `GRB131014215GBM40keV/` are GBM-only refits of GRB131014A |
| `LAT_analysis/` | LAT photon extraction (`lat_photons.csv`, the LAT-info table); see `LAT_analysis.md` |
| `light_curves/` | light-curve scripts and plots per burst |
| `080916C.tex`, `131014A.tex`, `140206B.tex`, `231129C.tex` | per-burst appendix tables generated from the fit logs by `generate_table_from_log.py` (the other four `.tex` files at the root are for bursts outside the sample) |
| `working.py`, `generate_table_from_log.py` | the fitting driver (writes `results.json` and a `cstat_run_*.log`) and the log-to-LaTeX table generator; both use `src.grb_research` imports and are meant to be run from the repo root |
| `md_file_refs/` | project working documents (see below) |
| `ToDo.md`, `PHASE5_TV_PLAN.md`, `NORRIS_TV_LORENTZ_INTEGRATION_PLAN.md` | science checklist and two Norris-`t_v` working plans (historical; each carries a dated status header); the Phase 6 plan is in `codes-for-paper/outflow_composition/` |

## Where to read next

- `md_file_refs/CLAUDE.md`: project rules (sample, episode naming, plotting and CSV conventions, per-folder method notes). Read before changing anything.
- `md_file_refs/HANDOFF.md`: the current state of the project and what is pending. Read first when resuming work.
- Also in `md_file_refs/`: `PLAN.md` (phased plan), `BUGS.md` (bug and observation log), `VERIFY.md` (pending PDF-build checks).
