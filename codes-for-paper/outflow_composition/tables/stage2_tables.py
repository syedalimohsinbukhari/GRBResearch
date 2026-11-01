"""
Phase 6 stage 2 (Zhang & Pe'er 2009 suppression bound) tables, generated from the CSV
======================================================================================
Reads ``../zhang_peer_suppression.csv`` and writes ``zhang_peer_tables.generated.md`` next to this script: the stage 2
results table and every number ``outflow_composition.md`` quotes from it (edge-limited rows, the extrapolated-band
spread, the exact L_gamma scaling check). The method note points at the generated file instead of carrying typed copies.

Kept separate from ``zhang_peer_suppression.py`` and from ``snapshot_tables.py`` (stage 1) on purpose (user,
2026-10-25). No randomness, hence no seed.

Run from anywhere:  python stage2_tables.py        (after zhang_peer_suppression.py)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CSV = HERE.parent / "zhang_peer_suppression.csv"
OUTPUT = HERE / "zhang_peer_tables.generated.md"

PRIMARY_L_BAND, SENSITIVITY_L_BAND = "zp09", "rest_1_1e4"
PRIMARY_TEST_BAND, SENSITIVITY_TEST_BAND = "fit", "extrapolated"
# A row is edge-limited when the binding blackbody/spectrum ratio sits at the test band's lower edge in at least this
# fraction of the MC draws (`binding_at_lower_edge_fraction`).
EDGE_LIMITED_THRESHOLD = 0.5
# Columns that must be constant over the CSV for the header to state them once.
CONSTANT_COLUMNS = ("H0", "Om0", "branch", "r0_convention", "l_w_over_l_gamma", "n_samples", "seed")


def number(value):
    """Three significant figures, no trailing zeros."""
    return f"{value:.3g}"


def span(values):
    """'a – b' for a range, 'a' when all values agree to three figures."""
    low, high = number(min(values)), number(max(values))
    return low if low == high else f"{low} – {high}"


def energy_label(kev):
    """keV, MeV or GeV label for an energy given in keV."""
    for scale, unit in ((1.0e6, "GeV"), (1.0e3, "MeV"), (1.0, "keV")):
        if kev >= scale:
            return f"{number(kev / scale)} {unit}"


def load():
    """CSV, asserting the columns the header states once are constant."""
    frame = pd.read_csv(CSV)
    for column in CONSTANT_COLUMNS:
        assert frame[column].nunique() == 1, f"{CSV.name}: column {column} is not constant"
    return frame


def band_description(frame, prefix):
    """'observer-frame 10 keV – 10 GeV' style label of the (frame, lo, hi) columns for each value of the band column."""
    out = []
    for name, group in frame.groupby(prefix, sort=False):
        first = group.iloc[0]
        lo, hi = first[f"{prefix}_lo_keV" if prefix == "band" else "l_band_lo_keV"], first[
            f"{prefix}_hi_keV" if prefix == "band" else "l_band_hi_keV"
        ]
        frame_name = "observer" if prefix == "band" else first.l_band_frame
        out.append(f"`{name}` = {frame_name}-frame {energy_label(lo)} – {energy_label(hi)}")
    return "; ".join(out)


def redshift_label(group):
    """'z = 4.35' for a spectroscopic burst, 'z = 1 … 7' for an assumed grid."""
    redshifts = sorted(group.z.unique())
    return f"z = {number(redshifts[0])}" if len(redshifts) == 1 else f"z = {number(redshifts[0])} … {number(redshifts[-1])}"


def select(frame, l_band, band):
    return frame[(frame.l_band == l_band) & (frame.band == band)]


def results_table(frame):
    """Per (burst, episode): 1 + sigma_min in the primary and the extrapolated test band, over the z grid."""
    primary = select(frame, PRIMARY_L_BAND, PRIMARY_TEST_BAND.replace("extrapolated", "fit"))
    extrapolated = select(frame, PRIMARY_L_BAND, SENSITIVITY_TEST_BAND)
    lines = [
        "| burst (z) | episode (model) | t_v [s] (source) | draws needing suppression | "
        "edge-limited draws | 1 + sigma_min, fit band | 1 + sigma_min, extrapolated band | kT_max [keV] | eta* |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for (burst, episode, model), group in primary.groupby(["grb_name", "episode", "model_name"], sort=False):
        other = extrapolated[(extrapolated.grb_name == burst) & (extrapolated.episode == episode)]
        sources = "/".join(sorted(group.t_v_source.unique()))
        lines.append(
            f"| {burst} ({redshift_label(group)}) | {episode} ({model}) | {span(group.t_v_s)} ({sources}) | "
            f"{span(group.needs_suppression_fraction)} | {span(group.binding_at_lower_edge_fraction)} | "
            f"{span(group.one_plus_sigma_min)} | {span(other.one_plus_sigma_min)} | {span(group.kT_max_keV)} | "
            f"{span(group.eta_star)} |"
        )
    return lines


def edge_limited(frame):
    """Lines on the rows whose binding point sits at the lower edge of the primary test band."""
    primary = select(frame, PRIMARY_L_BAND, PRIMARY_TEST_BAND)
    extrapolated = select(frame, PRIMARY_L_BAND, SENSITIVITY_TEST_BAND).set_index(["grb_name", "episode", "z"])
    flagged = primary[primary.binding_at_lower_edge_fraction >= EDGE_LIMITED_THRESHOLD]
    lines = [
        f"{len(flagged)} of {len(primary)} fit-band rows are edge-limited (`binding_at_lower_edge_fraction` >= "
        f"{number(EDGE_LIMITED_THRESHOLD)}).",
        "",
        "| burst | episode | rows (z values) |",
        "|---|---|---|",
    ]
    for (burst, episode), group in flagged.groupby(["grb_name", "episode"], sort=False):
        lines.append(f"| {burst} | {episode} | {len(group)} ({', '.join(number(z) for z in sorted(group.z))}) |")

    keys = list(zip(primary.grb_name, primary.episode, primary.z))
    ratio = pd.Series(
        [extrapolated.loc[k].one_plus_sigma_min / row for k, row in zip(keys, primary.one_plus_sigma_min)],
        index=primary.index,
    )
    worst = ratio.idxmax()
    lines += [
        "",
        f"Largest extrapolated-band / fit-band ratio of 1 + sigma_min over all {len(primary)} rows: {number(ratio.max())} "
        f"({primary.loc[worst].grb_name} {primary.loc[worst].episode}, z = {number(primary.loc[worst].z)}: "
        f"{number(primary.loc[worst].one_plus_sigma_min)} fit vs "
        f"{number(extrapolated.loc[keys[list(primary.index).index(worst)]].one_plus_sigma_min)} extrapolated).",
        f"Among the {len(primary) - len(flagged)} rows that are not edge-limited the largest ratio is "
        f"{number(ratio[~primary.index.isin(flagged.index)].max())}.",
    ]
    return lines


def l_gamma_scaling(frame):
    """The exact proportionality 1 + sigma_min ~ L_gamma between the two L_gamma bands, and the factor between them."""
    keys = ["grb_name", "episode", "z", "band"]
    primary = select(frame, PRIMARY_L_BAND, PRIMARY_TEST_BAND).merge(
        select(frame, SENSITIVITY_L_BAND, PRIMARY_TEST_BAND), on=keys, suffixes=("_zp09", "_rest")
    )
    both = pd.concat(
        [
            select(frame, PRIMARY_L_BAND, b).merge(select(frame, SENSITIVITY_L_BAND, b), on=keys, suffixes=("_zp09", "_rest"))
            for b in (PRIMARY_TEST_BAND, SENSITIVITY_TEST_BAND)
        ]
    )
    sigma_ratio = both.one_plus_sigma_min_zp09 / both.one_plus_sigma_min_rest
    l_ratio = both.L_gamma_erg_s_zp09 / both.L_gamma_erg_s_rest
    exact = sigma_ratio / l_ratio
    primary["factor"] = primary.one_plus_sigma_min_zp09 / primary.one_plus_sigma_min_rest
    lines = [
        f"Over all {len(both)} (episode, z, test band) pairs, [1 + sigma_min({PRIMARY_L_BAND}) / "
        f"1 + sigma_min({SENSITIVITY_L_BAND})] divided by [L_gamma({PRIMARY_L_BAND}) / L_gamma({SENSITIVITY_L_BAND})] is "
        f"{number(exact.median())} (range {number(exact.min())} – {number(exact.max())}).",
        "",
        f"Factor 1 + sigma_min({PRIMARY_L_BAND}) / 1 + sigma_min({SENSITIVITY_L_BAND}) in the fit band, per episode:",
        "",
        "| burst | episode | factor |",
        "|---|---|---|",
    ]
    for (burst, episode), group in primary.groupby(["grb_name", "episode"], sort=False):
        lines.append(f"| {burst} | {episode} | {span(group.factor)} |")
    return lines


def build():
    frame = load()
    row = frame.iloc[0]
    lines = [
        f"<!-- AUTO-GENERATED by tables/stage2_tables.py from {CSV.name}; do not edit by hand, rerun the script. -->",
        "",
        "# Phase 6 stage 2 — Zhang & Pe'er suppression bound (generated)",
        "",
        f"`{CSV.name}`: {len(frame)} rows; branch `{row.branch}`, R_0 = `{row.r0_convention}`, "
        f"L_w / L_gamma = {number(row.l_w_over_l_gamma)}; H0 = {number(row.H0)}, Om0 = {number(row.Om0)}; "
        f"n_samples = {int(row.n_samples)}, seed = {int(row.seed)}.",
        "",
        f"L_gamma bands: {band_description(frame, 'l_band')}. Test bands: {band_description(frame, 'band')}.",
        f"Primary: `{PRIMARY_L_BAND}` L_gamma, `{PRIMARY_TEST_BAND}` test band; every other row is a sensitivity row.",
        "Cells that span the redshift grid read `min – max`.",
        "",
        "## Results (primary L_gamma band; 1 + sigma_min, median)",
        "",
        *results_table(frame),
        "",
        "## Edge-limited rows and the test-band lower edge",
        "",
        *edge_limited(frame),
        "",
        "## Exact scaling with L_gamma",
        "",
        *l_gamma_scaling(frame),
        "",
    ]
    return "\n".join(lines)


def main():
    OUTPUT.write_text(build())
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
