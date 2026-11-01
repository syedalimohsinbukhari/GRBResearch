"""
Phase 6 stage 1 snapshot tables, generated from the CSVs
========================================================
Reads ``../gao_zhang_inversion.csv`` and ``../gao_zhang_inversion_partner.csv`` and writes
``outflow_composition_tables.generated.md`` next to this script: the results-snapshot tables and every number that
``outflow_composition.md`` quotes from them. The method note points at the generated file instead of carrying typed
copies, so the numbers cannot drift from the CSVs (CLAUDE.md: generated content is never hand-retyped).

Kept separate from ``gao_zhang_inversion.py`` on purpose (user, 2026-10-25): the computation and its presentation are
different scripts, as with the project's ``csv_to_latex.py`` helpers. No randomness, hence no seed.

Run from anywhere:  python snapshot_tables.py        (after gao_zhang_inversion.py)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUTFLOW_DIR = HERE.parent
MAIN_CSV = OUTFLOW_DIR / "gao_zhang_inversion.csv"
PARTNER_CSV = OUTFLOW_DIR / "gao_zhang_inversion_partner.csv"
OUTPUT = HERE / "outflow_composition_tables.generated.md"

# Regimes the Gao & Zhang inversion can solve; anything else ('none', 'ambiguous') has no 1 + sigma_0.
SOLVED_REGIMES = ("II", "III", "V", "VI")
# Redshift shown for the bursts without a measured z (the project's fiducial, CLAUDE.md "Output conventions").
Z_FIDUCIAL = 2.0
# regime_fraction below this marks a row whose MC draws straddle a regime boundary.
REGIME_FRACTION_FLAG = 0.9
R0_GZ2015_CM = (1.0e8, 1.0e9, 1.0e10)
# Columns that must be constant over a CSV for the header to state them once.
CONSTANT_COLUMNS = ("H0", "Om0", "e_min_keV", "e_max_keV", "f_gamma", "delta", "n_samples", "seed")


def number(value):
    """Three significant figures, no trailing zeros."""
    return f"{value:.3g}"


def span(values):
    """'a–b' for a range, 'a' when all values agree to three figures, 'none' for no values."""
    if len(values) == 0:
        return "none"
    low, high = number(min(values)), number(max(values))
    return low if low == high else f"{low}–{high}"


def energy_label(kev):
    """keV, MeV or GeV label for an energy given in keV."""
    for scale, unit in ((1.0e6, "GeV"), (1.0e3, "MeV"), (1.0, "keV")):
        if kev >= scale:
            return f"{number(kev / scale)} {unit}"


def exponent(value):
    """10^k label for a power-of-ten r_0."""
    return f"$10^{{{int(round(np.log10(value)))}}}$ cm"


def load(path):
    """CSV with a 'solved' flag; asserts the columns the header states once are in fact constant."""
    frame = pd.read_csv(path)
    for column in CONSTANT_COLUMNS:
        assert frame[column].nunique() == 1, f"{path.name}: column {column} is not constant"
    frame["solved"] = frame.regime.isin(SOLVED_REGIMES)
    return frame


def provenance(frame, name):
    """One line of the assumptions and MC settings that produced a CSV."""
    row = frame.iloc[0]
    return (
        f"`{name}`: {len(frame)} rows; observer-frame {energy_label(row.e_min_keV)} – {energy_label(row.e_max_keV)}; "
        f"f_gamma = {number(row.f_gamma)}, delta = {number(row.delta)}; H0 = {number(row.H0)}, Om0 = {number(row.Om0)}; "
        f"n_samples = {int(row.n_samples)}, seed = {int(row.seed)}"
    )


def reference_redshift(burst_frame):
    """The measured z for a spectroscopic burst, otherwise the fiducial one."""
    sources = burst_frame.z_source.unique()
    assert len(sources) == 1, f"mixed z_source for {burst_frame.grb_name.iloc[0]}"
    if sources[0] == "spectroscopic":
        return float(burst_frame.z.iloc[0])
    assert Z_FIDUCIAL in set(burst_frame.z), f"fiducial z = {Z_FIDUCIAL} missing"
    return Z_FIDUCIAL


def regime_label(frame):
    """'II', 'III' or 'II/III' for the regimes present among solved rows."""
    return "/".join(sorted(frame[frame.solved].regime.unique(), key=len))


def cell(frame):
    """Table cell for one (burst, r_0): range of 1 + sigma_0 with the regime(s), and the episodes with no regime."""
    solved, unsolved = frame[frame.solved], frame[~frame.solved]
    parts = []
    if len(solved):
        parts.append(f"{span(solved.one_plus_sigma0)} ({regime_label(frame)})")
    if len(unsolved):
        parts.append("none" if not len(solved) else f"none: {', '.join(unsolved.episode)}")
    return "; ".join(parts)


def counts(frame):
    """'N regime X, ...' summary of the regime column, solved regimes first."""
    solved = frame.regime.value_counts()
    order = [r for r in SOLVED_REGIMES if r in solved] + [r for r in solved.index if r not in SOLVED_REGIMES]
    return ", ".join(f"{solved[r]} {r}" for r in order)


def straddling(frame):
    """(number of solved rows below the regime_fraction flag, lowest regime_fraction among solved rows)."""
    solved = frame[frame.solved]
    return int((solved.regime_fraction < REGIME_FRACTION_FLAG).sum()), float(solved.regime_fraction.min())


def main_tables(frame):
    """Markdown for the main (BEST-BB) results snapshot."""
    gz = frame[frame.r0_set == "GZ2015"]
    lines = [
        "## Main table (BEST-BB episodes)",
        "",
        f"{int(frame.solved.sum())} of {len(frame)} trials have a regime ({counts(frame[frame.solved])}); "
        f"{int((~frame.solved).sum())} have none.",
        "{} of the solved trials have `regime_fraction` < {} (lowest {}).".format(
            *[number(v) for v in (straddling(frame)[0], REGIME_FRACTION_FLAG, straddling(frame)[1])]
        ),
        "",
        f"### 1 + sigma_0, `GZ2015` set, at the measured redshift (spectroscopic bursts) or z = {number(Z_FIDUCIAL)}",
        "",
        "| burst (episodes) | z | " + " | ".join(f"r_0 = {exponent(r)}" for r in R0_GZ2015_CM) + " |",
        "|---|---|" + "---|" * len(R0_GZ2015_CM),
    ]
    for burst, burst_frame in gz.groupby("grb_name", sort=False):
        z_ref = reference_redshift(burst_frame)
        at_z = burst_frame[np.isclose(burst_frame.z, z_ref)]
        cells = [cell(at_z[np.isclose(at_z.r0_cm, r)]) for r in R0_GZ2015_CM]
        lines.append(f"| {burst} ({at_z.episode.nunique()}) | {number(z_ref)} | " + " | ".join(cells) + " |")

    anchor = frame[(frame.r0_set == "project_anchor") & np.isclose(frame.r0_multiplier, 1.0) & frame.solved]
    low, high = anchor.loc[anchor.one_plus_sigma0.idxmin()], anchor.loc[anchor.one_plus_sigma0.idxmax()]
    lines += [
        "",
        "### Project anchor x 1 (the pure-fireball r_0 of each burst)",
        "",
        f"Per-episode 1 + sigma_0 over all {len(anchor)} solved rows: {number(low.one_plus_sigma0)} "
        f"({low.grb_name} {low.episode}, z = {number(low.z)}) to {number(high.one_plus_sigma0)} "
        f"({high.grb_name} {high.episode}, z = {number(high.z)}).",
        "",
        "### Redshift dependence, `GZ2015` set, r_0 = $10^{9}$ cm (bursts with an assumed z)",
        "",
    ]
    swept = gz[(gz.z_source != "spectroscopic") & np.isclose(gz.r0_cm, 1.0e9) & gz.solved]
    redshifts = sorted(swept.z.unique())
    for quantity, header in (("one_plus_sigma0", "1 + sigma_0"), ("eta", "eta")):
        lines += [
            f"{header}, range over the burst's solved episodes:",
            "",
            "| burst | " + " | ".join(f"z = {number(z)}" for z in redshifts) + " |",
            "|---|" + "---|" * len(redshifts),
        ]
        for burst, burst_frame in swept.groupby("grb_name", sort=False):
            cells = [span(burst_frame[np.isclose(burst_frame.z, z)][quantity]) for z in redshifts]
            lines.append(f"| {burst} | " + " | ".join(cells) + " |")
        lines.append("")
    return lines


def partner_tables(frame):
    """Markdown for the partner (non-BEST blackbody) models."""
    gz = frame[frame.r0_set == "GZ2015"]
    lines = [
        "## Partner models (X_BB of episodes whose BEST model has no blackbody)",
        "",
        f"{int(frame.solved.sum())} of {len(frame)} trials have a regime ({counts(frame[frame.solved])}); "
        f"{int((~frame.solved).sum())} have none; {int(frame.physical.sum())} are `physical`.",
        "",
        f"### 1 + sigma_0, `GZ2015` set, at the measured redshift (spectroscopic bursts) or z = {number(Z_FIDUCIAL)}",
        "",
        "| episode (partner) | kT_BB [keV] | delta C-stat vs BEST | "
        + " | ".join(f"r_0 = {exponent(r)}" for r in R0_GZ2015_CM)
        + " |",
        "|---|---|---|" + "---|" * len(R0_GZ2015_CM),
    ]
    for (burst, episode, model), group in gz.groupby(["grb_name", "episode", "model_name"], sort=False):
        z_ref = reference_redshift(group)
        at_z = group[np.isclose(group.z, z_ref)]
        first = at_z.iloc[0]
        status = "" if first.bb_status == "SAFE" else f", {first.bb_status}"
        cells = []
        for r0 in R0_GZ2015_CM:
            row = at_z[np.isclose(at_z.r0_cm, r0)].iloc[0]
            cells.append(f"{number(row.one_plus_sigma0)} ({row.regime})" if row.solved else "none")
        lines.append(
            f"| {burst} {episode} ({model}{status}) | {number(first.kt_bb_keV)} | "
            f"{number(first.delta_cstat_vs_best)} | " + " | ".join(cells) + " |"
        )
    lines.append("")
    return lines


def build():
    """Full markdown text of the generated file."""
    main, partner = load(MAIN_CSV), load(PARTNER_CSV)
    header = [
        f"<!-- AUTO-GENERATED by tables/snapshot_tables.py from {MAIN_CSV.name} and {PARTNER_CSV.name}; "
        "do not edit by hand, rerun the script. -->",
        "",
        "# Phase 6 stage 1 — snapshot tables (generated)",
        "",
        "Sources and the assumptions that produced them:",
        "",
        f"- {provenance(main, MAIN_CSV.name)}",
        f"- {provenance(partner, PARTNER_CSV.name)}",
        "",
        "A cell reads `range of 1 + sigma_0 (regime)`; `none` means no Gao & Zhang regime applies at that r_0, "
        "and a list after `none:` names the episodes without one.",
        "",
    ]
    return "\n".join(header + main_tables(main) + partner_tables(partner)) + "\n"


def main():
    OUTPUT.write_text(build())
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
