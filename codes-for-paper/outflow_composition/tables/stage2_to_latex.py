"""
Phase 6 stage 2 LaTeX table (Zhang & Pe'er 2009 minimum-magnetisation bound)
===========================================================================
Pure CSV-in, tex-out: reads ``../zhang_peer_suppression.csv`` and writes ``outflow_stage2_table.tex`` next to this
script: 1 + sigma_min of the 13 episodes whose BEST model has no blackbody, over the assumed-redshift grid (one column per
z) and at the spectroscopic redshift of GRB080916C. Only the primary L_gamma band and test band are tabulated; the
sensitivity rows stay in the CSV and in ``zhang_peer_tables.generated.md``. No randomness, hence no seed.

Per-cell and per-episode markers: a dagger on the episode when the episode duration stands in for t_v, a double dagger on
a value whose binding blackbody/spectrum ratio sits at the test band's lower edge in at least half of the draws.

Run:  python stage2_to_latex.py        (after zhang_peer_suppression.py; not with ``-I``, it imports ``latex_common``)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from latex_common import (
    BIB_ZHANG_PEER,
    band_si,
    decimals_for,
    format_grb_name,
    format_model_name,
    format_value,
    wrap_table,
)

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "zhang_peer_suppression.csv"
OUTPUT = HERE / "outflow_stage2_table.tex"
SCRIPT = "codes-for-paper/outflow_composition/tables/stage2_to_latex.py"

PRIMARY_L_BAND, PRIMARY_TEST_BAND = "zp09", "fit"
# Same threshold as `stage2_tables.py`: a value is band-limited when at least this fraction of the draws bind at the
# lower edge of the test band.
EDGE_LIMITED_THRESHOLD = 0.5
CONSTANT_COLUMNS = ("branch", "r0_convention", "l_w_over_l_gamma", "n_samples", "seed")


def load():
    """Primary-band rows of the CSV; asserts the columns the caption states once are constant."""
    frame = pd.read_csv(CSV)
    for column in CONSTANT_COLUMNS:
        assert frame[column].nunique() == 1, f"{CSV.name}: column {column} is not constant"
    assert frame.l_w_over_l_gamma.iloc[0] == 1.0, "caption states L_w = L_gamma"
    assert frame.branch.iloc[0] == "eta>=eta_star" and frame.r0_convention.iloc[0] == "c*t_v/(1+z)"
    return frame, frame[(frame.l_band == PRIMARY_L_BAND) & (frame.band == PRIMARY_TEST_BAND)]


def band_of(frame, column_prefix, name_column, name):
    """(lo, hi) keV of a named band."""
    row = frame[frame[name_column] == name].iloc[0]
    return row[f"{column_prefix}_lo_keV"], row[f"{column_prefix}_hi_keV"]


def build(frame, primary):
    """The table text."""
    spectroscopic = sorted(primary[primary.z_source == "spectroscopic"].z.unique())
    assumed = sorted(primary[primary.z_source != "spectroscopic"].z.unique())
    z_columns = spectroscopic + assumed
    n_cols = 2 + len(z_columns)
    # Decimals per z column, from the smallest error displayed in that column.
    decimals = {}
    for z in z_columns:
        column = primary[np.isclose(primary.z, z)]
        decimals[z] = decimals_for(list(column.one_plus_sigma_min_err_lower) + list(column.one_plus_sigma_min_err_upper))

    rows, current, duration, edge = "", None, False, False
    for (grb, episode, model), group in primary.groupby(["grb_name", "episode", "model_name"], sort=False):
        if grb != current:
            if current is not None:
                rows += "\\midrule\n"
            rows += f"\\multicolumn{{{n_cols}}}{{l}}{{\\textbf{{{format_grb_name(grb)}}}}} \\\\\n"
            current = grb
        uses_duration = bool(group.t_v_is_duration_fallback.any())
        duration |= uses_duration
        label = f"{episode}$^{{\\dagger}}$" if uses_duration else episode
        cells = []
        for z in z_columns:
            match = group[np.isclose(group.z, z)]
            if match.empty:
                cells.append("")
                continue
            row = match.iloc[0]
            cell = format_value(
                row.one_plus_sigma_min, row.one_plus_sigma_min_err_lower, row.one_plus_sigma_min_err_upper, decimals[z]
            )
            if row.binding_at_lower_edge_fraction >= EDGE_LIMITED_THRESHOLD:
                cell += "$^{\\ddagger}$"
                edge = True
            cells.append(cell)
        rows += f"{label} & {format_model_name(model)} & " + " & ".join(cells) + " \\\\\n"

    notes = []
    if duration:
        notes.append(
            (r"$\dagger$", "No pulse-based variability timescale is available and the episode duration stands in for "
                           "$t_\\mathrm{v}$. Unlike $\\Gamma_{\\min}$, $1+\\sigma_\\mathrm{min}$ is not monotonic in "
                           "$t_\\mathrm{v}$, so the entry may lie above or below its value for the true timescale.")
        )
    if edge:
        low, _ = band_of(frame, "band", "band", PRIMARY_TEST_BAND)
        notes.append(
            (r"$\ddagger$", f"The binding blackbody-to-spectrum ratio lies at the lower edge of the test band "
                            f"(\\SI{{{low:g}}}{{\\kilo\\electronvolt}}) in at least "
                            f"\\SI{{{100 * EDGE_LIMITED_THRESHOLD:g}}}{{\\percent}} of the Monte Carlo draws.")
        )

    row0 = primary.iloc[0]
    l_lo, l_hi = band_of(frame, "l_band", "l_band", PRIMARY_L_BAND)
    t_lo, t_hi = band_of(frame, "band", "band", PRIMARY_TEST_BAND)
    caption = (
        f"Minimum magnetisation following \\citet{{{BIB_ZHANG_PEER}}}: the smallest $1+\\sigma$ for which the "
        "photospheric blackbody of the fireball stays below the best-fit spectrum across the test band, for the "
        "episodes whose best-fit model has no blackbody.\n"
        f"The test band is the observer-frame {band_si(t_lo, t_hi)}, the gamma-ray luminosity $L_\\gamma$ is integrated "
        f"over the observer-frame {band_si(l_lo, l_hi)}, and $L_\\mathrm{{w}}=L_\\gamma$.\n"
        "The base radius is $R_0=c\\,t_\\mathrm{v}/(1+z)$ and the branch is $\\eta\\geq\\eta_*$.\n"
        f"Values are quoted at the spectroscopic redshift $z={spectroscopic[0]:g}$ for "
        "\\grbzeroeightzeroninesixteenC\\ and at the assumed redshifts "
        f"$z={','.join(f'{z:g}' for z in assumed)}$ for the other bursts.\n"
        f"Errors are asymmetric $1\\sigma$ intervals from ${int(row0.n_samples)}$ Monte Carlo samples "
        f"(seed ${int(row0.seed)}$)."
    )
    head = "Episode & Model & " + " & ".join(f"$z={z:g}$" for z in z_columns) + " \\\\\n"
    return wrap_table(
        header_comment=SCRIPT,
        caption=caption,
        label="tab:outflow_stage2",
        colspec="ll" + "c" * len(z_columns),
        head=head,
        rows=rows,
        notes=notes,
        star=True,
    )


def main():
    frame, primary = load()
    OUTPUT.write_text(build(frame, primary))
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
