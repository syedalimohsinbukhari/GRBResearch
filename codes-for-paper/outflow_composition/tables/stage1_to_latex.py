"""
Phase 6 stage 1 LaTeX tables (Gao & Zhang 2015 hybrid-outflow inversion)
===========================================================================
Pure CSV-in, tex-out: reads ``../gao_zhang_inversion.csv`` and ``../gao_zhang_inversion_partner.csv`` and writes three
tables next to this script, in the layout of ``photospheric_radius/csv_to_latex.py``:

    outflow_stage1_table.tex            BEST-BB episodes, 1 + sigma_0 at three assumed r_0 (the paper's main table)
    outflow_stage1_unknown_z_table.tex  the same episodes of the three bursts without a redshift, over z = 1, 2, 3, 5, 7
    outflow_partner_table.tex           the SAFE/MARGINAL X_BB partner of every non-BB episode, reported separately

The model-set rules (user, 2026-10-24): the main table holds only episodes whose BEST model has a blackbody; non-BEST
blackbody models are a separate table and MARGINAL ones are flagged in red, as the appendix tables do; no UNSAFE
blackbody model is used. No randomness, hence no seed. Never hand-edit the generated .tex.

Run:  python stage1_to_latex.py        (after gao_zhang_inversion.py; not with ``-I``, it imports ``latex_common``)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from latex_common import (
    BIB_GAO,
    MISSING,
    band_si,
    decimals_for,
    exponent,
    format_grb_name,
    format_model_name,
    format_value,
    red,
    wrap_table,
)

HERE = Path(__file__).resolve().parent
MAIN_CSV = HERE.parent / "gao_zhang_inversion.csv"
PARTNER_CSV = HERE.parent / "gao_zhang_inversion_partner.csv"
SCRIPT = "codes-for-paper/outflow_composition/tables/stage1_to_latex.py"

SOLVED_REGIMES = ("II", "III", "V", "VI")
Z_FIDUCIAL = 2.0
R0_GZ2015_CM = (1.0e8, 1.0e9, 1.0e10)
# r_0 at which the redshift sweep is tabulated (the middle trial; the sweep over r_0 is the main table's job).
R0_SWEEP_CM = 1.0e9
CONSTANT_COLUMNS = ("f_gamma", "delta", "e_min_keV", "e_max_keV", "n_samples", "seed")


def load(path):
    """CSV with a ``solved`` flag; asserts the columns the captions state once are constant."""
    frame = pd.read_csv(path)
    for column in CONSTANT_COLUMNS:
        assert frame[column].nunique() == 1, f"{path.name}: column {column} is not constant"
    assert abs(frame.delta.iloc[0] - 1 / 3) < 1e-9, "captions state delta = 1/3"
    frame["solved"] = frame.regime.isin(SOLVED_REGIMES)
    return frame


def reference_redshift(group):
    """The measured z for a spectroscopic burst, otherwise the fiducial one."""
    if (group.z_source == "spectroscopic").all():
        return float(group.z.iloc[0])
    assert Z_FIDUCIAL in set(group.z), f"fiducial z = {Z_FIDUCIAL} missing"
    return Z_FIDUCIAL


def at(frame, *, z, r0):
    """The single row of an episode's group at one (z, r_0) of the GZ2015 set."""
    rows = frame[(frame.r0_set == "GZ2015") & np.isclose(frame.z, z) & np.isclose(frame.r0_cm, r0)]
    assert len(rows) == 1, f"expected one row, got {len(rows)}"
    return rows.iloc[0]


def sigma_cell(row, decimals, flag=lambda text: text):
    """1 + sigma_0 with its asymmetric errors, or the missing marker when no regime applies."""
    if not row.solved:
        return flag(MISSING)
    return flag(
        format_value(
            row.one_plus_sigma0, row.one_plus_sigma0_err_lower, row.one_plus_sigma0_err_upper, decimals
        )
    )


def sigma_decimals(frame, r0_values, reference):
    """Per-column decimals from the smallest error among the rows actually displayed."""
    out = {}
    for r0 in r0_values:
        errors = []
        for _, group in frame[frame.r0_set == "GZ2015"].groupby(["grb_name", "episode", "model_name"], sort=False):
            z = reference(group)
            rows = group[np.isclose(group.z, z) & np.isclose(group.r0_cm, r0) & group.solved]
            errors += list(rows.one_plus_sigma0_err_lower) + list(rows.one_plus_sigma0_err_upper)
        out[r0] = decimals_for(errors)
    return out


def provenance(frame):
    """(band, f_gamma, n_samples, seed) for the captions."""
    row = frame.iloc[0]
    return band_si(row.e_min_keV, row.e_max_keV), row.f_gamma, int(row.n_samples), int(row.seed)


def r0_header(r0):
    return f"$r_0=10^{{{exponent(r0)}}}$ cm"


def sigma_head(lead_labels):
    """Two-row header: the leading columns span both rows (multirow), 1 + sigma_0 spans the r_0 columns, r_0 below it."""
    n_lead, n_r0 = len(lead_labels), len(R0_GZ2015_CM)
    first = " & ".join(f"\\multirow{{2}}{{*}}{{{label}}}" for label in lead_labels)
    first += f" & \\multicolumn{{{n_r0}}}{{c}}{{$1+\\sigma_0$}} \\\\\n"
    rule = f"\\cmidrule(lr){{{n_lead + 1}-{n_lead + n_r0}}}\n"
    second = " & " * n_lead + " & ".join(r0_header(r) for r in R0_GZ2015_CM) + " \\\\\n"
    return first + rule + second


def z_label(z, measured):
    """Redshift for a burst header row; the dagger marks an assumed fiducial value."""
    return f"$z={z:.2f}$" if measured else f"$z={z:.2f}^{{\\dagger}}$"


def grb_header(grb_name, n_cols, z=None, measured=True):
    """Burst header row spanning the table, with the redshift the burst's values are quoted at."""
    suffix = "" if z is None else f" ({z_label(z, measured)})"
    return f"\\multicolumn{{{n_cols}}}{{l}}{{\\textbf{{{format_grb_name(grb_name)}}}{suffix}}} \\\\\n"


def caption_tail(frame):
    """Sentences shared by every stage 1 caption: band, constants, MC settings."""
    band, f_gamma, n_samples, seed = provenance(frame)
    return (
        f"The blackbody temperature and the thermal and total fluxes are integrated over the observer-frame {band}.\n"
        f"The radiative efficiency is $f_\\gamma={f_gamma:g}$ and the slow-acceleration index is $\\delta=1/3$.\n"
        f"Errors are asymmetric $1\\sigma$ intervals from ${n_samples}$ Monte Carlo samples (seed ${seed}$)."
    )


def main_table(frame):
    """BEST-BB episodes at the reference redshift, 1 + sigma_0 for the three assumed r_0."""
    decimals = sigma_decimals(frame, R0_GZ2015_CM, reference_redshift)
    rows, current, assumed = "", None, False
    for (grb, episode, model), group in frame.groupby(["grb_name", "episode", "model_name"], sort=False):
        if grb != current:
            if current is not None:
                rows += "\\midrule\n"
            z = reference_redshift(group)
            measured = (group.z_source == "spectroscopic").all()
            rows += grb_header(grb, 2 + len(R0_GZ2015_CM), z, measured)
            current = grb
        z = reference_redshift(group)
        measured = (group.z_source == "spectroscopic").all()
        assumed |= not measured
        cells = " & ".join(sigma_cell(at(group, z=z, r0=r0), decimals[r0]) for r0 in R0_GZ2015_CM)
        rows += f"{format_model_name(model)} & {episode} & {cells} \\\\\n"

    notes = []
    if assumed:
        notes.append((r"$\dagger$", f"Assumed fiducial redshift; no spectroscopic measurement exists for this burst."))
    if MISSING in rows.replace("\\\\", ""):
        notes.append(
            (
                MISSING,
                "No photospheric regime of \\citet{" + BIB_GAO + "} applies at this $r_0$ "
                "($r_\\mathrm{ph}$ below the radius at which acceleration ends, or no regime criterion is met).",
            )
        )
    caption = (
        f"Hybrid-outflow inversion of the thermal component following \\citet{{{BIB_GAO}}}: the magnetisation "
        "parameter $1+\\sigma_0$ of the episodes whose best-fit model contains a blackbody, for an assumed "
        "jet-base radius $r_0$.\n"
        "Only \\grbzeroeightzeroninesixteenC\\ has a spectroscopic redshift; for the remaining bursts values are "
        f"quoted at a fiducial $z={Z_FIDUCIAL:g}$.\n" + caption_tail(frame)
    )
    head = sigma_head(["Model", "Episode"])
    return wrap_table(
        header_comment=SCRIPT,
        caption=caption,
        label="tab:outflow_stage1",
        colspec="ll" + "c" * len(R0_GZ2015_CM),
        head=head,
        rows=rows,
        notes=notes,
    )


def sweep_table(frame):
    """The assumed-redshift bursts over the z grid, at the single r_0 = R0_SWEEP_CM."""
    swept = frame[frame.z_source != "spectroscopic"]
    z_values = sorted(swept.z.unique())
    errors = swept[
        (swept.r0_set == "GZ2015") & np.isclose(swept.r0_cm, R0_SWEEP_CM) & swept.solved
    ]
    decimals = decimals_for(list(errors.one_plus_sigma0_err_lower) + list(errors.one_plus_sigma0_err_upper))
    rows, current = "", None
    for (grb, episode, model), group in swept.groupby(["grb_name", "episode", "model_name"], sort=False):
        if grb != current:
            if current is not None:
                rows += "\\midrule\n"
            rows += grb_header(grb, 2 + len(z_values))
            current = grb
        cells = " & ".join(sigma_cell(at(group, z=z, r0=R0_SWEEP_CM), decimals) for z in z_values)
        rows += f"{format_model_name(model)} & {episode} & {cells} \\\\\n"

    notes = []
    if MISSING in rows.replace("\\\\", ""):
        notes.append((MISSING, f"No photospheric regime of \\citet{{{BIB_GAO}}} applies at this redshift."))
    caption = (
        f"Hybrid-outflow inversion following \\citet{{{BIB_GAO}}} for the three bursts without a measured redshift: "
        f"$1+\\sigma_0$ at an assumed $r_0=10^{{{exponent(R0_SWEEP_CM)}}}$ cm and the assumed redshifts "
        f"$z={','.join(f'{z:g}' for z in z_values)}$.\n" + caption_tail(frame)
    )
    head = "Model & Episode & " + " & ".join(f"$z={z:g}$" for z in z_values) + " \\\\\n"
    return wrap_table(
        header_comment=SCRIPT,
        caption=caption,
        label="tab:outflow_stage1_unknown_z",
        colspec="ll" + "c" * len(z_values),
        head=head,
        rows=rows,
        notes=notes,
    )


def partner_table(frame):
    """The X_BB partner of every non-BB episode: separate from the main table, MARGINAL rows in red."""
    decimals = sigma_decimals(frame, R0_GZ2015_CM, reference_redshift)
    rows, current, assumed, marginal = "", None, False, False
    for (grb, episode, model), group in frame.groupby(["grb_name", "episode", "model_name"], sort=False):
        if grb != current:
            if current is not None:
                rows += "\\midrule\n"
            z = reference_redshift(group)
            measured = (group.z_source == "spectroscopic").all()
            rows += grb_header(grb, 4 + len(R0_GZ2015_CM), z, measured)
            current = grb
        z = reference_redshift(group)
        measured = (group.z_source == "spectroscopic").all()
        assumed |= not measured
        first = group.iloc[0]
        flag = red if first.bb_status == "MARGINAL" else (lambda text: text)
        marginal |= first.bb_status == "MARGINAL"
        cells = " & ".join(sigma_cell(at(group, z=z, r0=r0), decimals[r0], flag) for r0 in R0_GZ2015_CM)
        lead = " & ".join(
            flag(text)
            for text in (
                episode,
                format_model_name(first.best_model),
                format_model_name(model),
                f"${first.delta_cstat_vs_best:.1f}$",
            )
        )
        rows += f"{lead} & {cells} \\\\\n"

    notes = []
    if assumed:
        notes.append((r"$\dagger$", "Assumed fiducial redshift; no spectroscopic measurement exists for this burst."))
    if marginal:
        notes.append(("red", "Blackbody model classed MARGINAL; all other partner models are SAFE."))
    if MISSING in rows.replace("\\\\", ""):
        notes.append((MISSING, f"No photospheric regime of \\citet{{{BIB_GAO}}} applies at this $r_0$."))
    threshold = frame.bb_threshold_delta_cstat.iloc[0]
    caption = (
        f"Hybrid-outflow inversion following \\citet{{{BIB_GAO}}} applied to the blackbody-bearing partner of the "
        "episodes whose best-fit model has no blackbody: the best-fit continuum plus a blackbody.\n"
        f"Every partner improves the fit by less than the blackbody-selection threshold of $\\Delta C={threshold:g}$ "
        "(column $\\Delta C$ gives the improvement over the best-fit model).\n" + caption_tail(frame)
    )
    head = sigma_head(["Episode", "Best model", "Partner", "$\\Delta C$"])
    assert (frame.delta_cstat_vs_best < threshold).all(), "a partner above the selection threshold would be BEST"
    return wrap_table(
        header_comment=SCRIPT,
        caption=caption,
        label="tab:outflow_partner",
        colspec="lllc" + "c" * len(R0_GZ2015_CM),
        head=head,
        rows=rows,
        notes=notes,
        star=True,
    )


def main():
    main_frame, partner_frame = load(MAIN_CSV), load(PARTNER_CSV)
    outputs = {
        "outflow_stage1_table.tex": main_table(main_frame),
        "outflow_stage1_unknown_z_table.tex": sweep_table(main_frame),
        "outflow_partner_table.tex": partner_table(partner_frame),
    }
    for name, text in outputs.items():
        (HERE / name).write_text(text)
        print(f"Saved: {HERE / name}")


if __name__ == "__main__":
    main()
