"""Build a machine-readable table of the per-episode Fermi-LAT analysis products.

Source of truth is the pair of files gtburst/gtlike leaves in each episode
directory:

    LAT_analysis/<grb>/Ep*/<name>_analysis_result_<start>_<stop>.txt
    LAT_analysis/<grb>/Ep*/<name>_fit_results_<start>_<stop>.txt

Everything the paper's LAT appendix table reports lives in those two files, so
this script collects them into `lat_photons.csv` and both consumers -- the
appendix table (`csv_to_latex.py`) and the gamma-gamma opacity calculation
(`codes-for-paper/lorentz_factor/lorentz_factor.py`) -- read that one CSV.
Before this existed the photon energies were transcribed by hand into
`lorentz_factor.py`, which is the drift trap recorded as BUGS.md OBS-08.

Units: photon energies are **MeV**. The field is named "P > 0.9 Max (E) MeV" in
the source, and the LAT selection floor here is 100 MeV -- the smallest values
in the sample (117.7, 122.9, 194.5) sit just above it. See BUGS.md BUG-18.

Created on Jun 13 19:21:38 2026
"""

import re
from pathlib import Path
from typing import Any

import pandas as pd

from grb_research import find_project_root
from grb_research.grb_time import TimeInterval

# Trigger-id directory -> the burst's name in the paper. The 231129C directory
# is already named for the paper burst; the other three are named for the
# Fermi trigger id, hence an explicit map rather than a pattern.
GRB_DIRECTORIES: dict[str, tuple[str, str]] = {
    "018__GRB080916009": ("GRB080916C", "GRB080916009"),
    "014__GRB131014215": ("GRB131014A", "GRB131014215"),
    "007__GRB140206275": ("GRB140206B", "GRB140206275"),
    "GRB231129C": ("GRB231129C", "GRB231129779"),
}

# Interval bounds are written to 3 dp in results.json but with trailing zeros
# stripped in the directory names (1.28 vs 1.280), so they are matched on value.
BOUND_TOLERANCE_S = 1e-6


def parse_key_value_file(file_path: Path) -> dict[str, str]:
    """Parse comma-separated key/value records from a file."""
    data = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if ', ' in line:
                key, value = line.split(', ', 1)
                data[key.strip()] = value.strip()
    return data


def parse_interval_from_directory_name(dir_name: str) -> tuple[str, str]:
    """
    Parse interval bounds from directory name like Ep*__{start}_{end}.
    Convert leading 'm' to '-'.
    Returns (start, end) as strings preserving original precision.
    """
    match = re.match(r'^Ep\w*__([m\-]?[\d.]+)_([\d.]+)$', dir_name)
    if not match:
        raise ValueError(f"Malformed directory name: {dir_name}")

    start_str = match.group(1)
    end_str = match.group(2)

    # Convert m prefix to minus sign
    if start_str.startswith('m'):
        start_str = '-' + start_str[1:]

    return start_str, end_str


def find_companion_files(epoch_dir: Path) -> tuple[Path, Path]:
    """
    Find the analysis and fit result files in the epoch directory.
    Raises error if missing or duplicate files are found.
    """
    analysis_files = list(epoch_dir.glob('*_analysis_result_*.txt'))
    fit_files = list(epoch_dir.glob('*_fit_results_*.txt'))

    if len(analysis_files) == 0:
        raise FileNotFoundError(f"Missing analysis file in {epoch_dir}")
    if len(analysis_files) > 1:
        raise ValueError(f"Duplicate analysis files in {epoch_dir}: {analysis_files}")

    if len(fit_files) == 0:
        raise FileNotFoundError(f"Missing fit file in {epoch_dir}")
    if len(fit_files) > 1:
        raise ValueError(f"Duplicate fit files in {epoch_dir}: {fit_files}")

    return analysis_files[0], fit_files[0]


def process_epoch(epoch_dir: Path) -> dict[str, Any]:
    """Process a single epoch directory and extract required data."""
    start_str, end_str = parse_interval_from_directory_name(epoch_dir.name)

    analysis_file, fit_file = find_companion_files(epoch_dir)

    analysis_data = parse_key_value_file(analysis_file)
    fit_data = parse_key_value_file(fit_file)

    required_analysis = ['# of Events', '# of P > 0.9', 'P > 0.9 Max (E) MeV', 'Arrival Time (s)', 'TS']
    for field in required_analysis:
        if field not in analysis_data:
            raise ValueError(f"Missing required field '{field}' in {analysis_file}")

    required_fit = ['Index', 'Index Error', 'Flux (0.1 - 100.0) GeV', 'Flux Error (0.1 - 100.0) GeV']
    for field in required_fit:
        if field not in fit_data:
            raise ValueError(f"Missing required field '{field}' in {fit_file}")

    return {
        'start': start_str,
        'end': end_str,
        't_start_s': float(start_str),
        't_stop_s': float(end_str),
        'n_events': int(analysis_data['# of Events']),
        'n_events_high_prob': int(analysis_data['# of P > 0.9']),
        'e_max_MeV': float(analysis_data['P > 0.9 Max (E) MeV']),
        't_arr_s': float(analysis_data['Arrival Time (s)']),
        'ts': float(analysis_data['TS']),
        'photon_index': float(fit_data['Index']),
        'photon_index_err': float(fit_data['Index Error']),
        'flux_0p1_100_GeV_ph_cm2_s': float(fit_data['Flux (0.1 - 100.0) GeV']),
        'flux_err_0p1_100_GeV_ph_cm2_s': float(fit_data['Flux Error (0.1 - 100.0) GeV']),
    }


def episode_labels_from_results(results_path: Path) -> dict[str, dict[tuple[float, float], str]]:
    """Map each GRB's (start, stop) interval to its canonical episode label.

    The labels (T90, TR1, EX0, ...) are defined by the interval strings in
    results.json and parsed by grb_research, so this script never invents a
    naming rule of its own. See CLAUDE.md, "Episode naming".
    """
    import json

    with open(results_path) as handle:
        results = json.load(handle)

    labels: dict[str, dict[tuple[float, float], str]] = {}
    for trigger_id, episodes in results.items():
        by_bounds = {}
        for key in episodes:
            interval = TimeInterval.from_string(key)
            if interval.start is None:
                continue
            kind = interval.kind.name
            label = f"{kind}{interval.index}" if kind in ("TR", "SP") else kind
            by_bounds[(interval.start, interval.end)] = label
        labels[trigger_id] = by_bounds
    return labels


def lookup_label(by_bounds: dict[tuple[float, float], str], start: float, stop: float, where: Path) -> str:
    """Find the episode label whose interval matches these bounds."""
    for (a, b), label in by_bounds.items():
        if abs(a - start) < BOUND_TOLERANCE_S and abs(b - stop) < BOUND_TOLERANCE_S:
            return label
    raise ValueError(f"No results.json interval matches {start}_{stop} for {where}")


def main():
    """Collect every episode of every GRB into a single CSV."""
    script_dir = Path(__file__).parent
    root = find_project_root()
    labels = episode_labels_from_results(root / "results.json")

    rows = []
    for dir_name, (grb_name, trigger_id) in sorted(GRB_DIRECTORIES.items()):
        grb_dir = script_dir / dir_name
        if not grb_dir.is_dir():
            raise FileNotFoundError(f"Missing GRB directory: {grb_dir}")

        epoch_dirs = [d for d in grb_dir.iterdir() if d.is_dir() and d.name.startswith('Ep')]
        if not epoch_dirs:
            raise ValueError(f"No epoch directories found in {grb_dir}")

        for epoch_dir in epoch_dirs:
            data = process_epoch(epoch_dir)
            data['grb_name'] = grb_name
            data['grb_trigger_id'] = trigger_id
            data['episode'] = lookup_label(labels[trigger_id], data['t_start_s'], data['t_stop_s'], epoch_dir)
            rows.append(data)
        print(f"{dir_name}: {len(epoch_dirs)} episodes")

    columns = ['grb_name', 'grb_trigger_id', 'episode', 't_start_s', 't_stop_s',
               'n_events', 'n_events_high_prob', 'e_max_MeV', 't_arr_s', 'ts',
               'photon_index', 'photon_index_err',
               'flux_0p1_100_GeV_ph_cm2_s', 'flux_err_0p1_100_GeV_ph_cm2_s']
    frame = pd.DataFrame(rows)[columns]
    frame = frame.sort_values(['grb_name', 't_start_s', 't_stop_s']).reset_index(drop=True)

    out_path = script_dir / "lat_photons.csv"
    frame.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}  ({len(frame)} rows)")


if __name__ == '__main__':
    main()
