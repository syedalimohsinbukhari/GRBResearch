"""
Lorentz factor LaTeX table generator, redshift sweep (Limits A and B)
=====================================================================
Pure CSV-in, tex-out.  Renders one table per limit for the three bursts without a spectroscopic redshift, with one
Gamma_min column per swept z:

    python generate_lorentz_table_unknown_z.py A   lorentz_results_unknown_z.csv         -> lorentz_table_unknown_z.tex
    python generate_lorentz_table_unknown_z.py B   lorentz_results_limit_b_unknown_z.csv -> lorentz_table_limit_b_unknown_z.tex

With no argument both are rendered (each CSV must then exist).  The CSVs come from ``lorentz_factor_unknown_z.py`` and
``lorentz_factor_limit_b_unknown_z.py`` respectively, which shell out to this script by subprocess.  Same layout idea as
``amati_relationship/csv_to_latex_unknown_z.py``.  No dependency on ``grb_research`` or ``results.json``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from generate_lorentz_table import TEX_NAMES, episode_order

HERE = Path(__file__).parent
INPUTS = {
    "A": (HERE / "lorentz_results_unknown_z.csv", HERE / "lorentz_table_unknown_z.tex"),
    "B": (HERE / "lorentz_results_limit_b_unknown_z.csv", HERE / "lorentz_table_limit_b_unknown_z.tex"),
}

Z_VALUES = (1, 2, 3, 5, 7)
MISSING = "-"


def load_results(csv_path):
    """Read the CSV into list-of-dicts, converting pandas' NaN back to None (as the other generators do)."""
    df = pd.read_csv(csv_path)
    return [{k: (None if pd.isna(v) else v) for k, v in row.items()} for row in df.to_dict("records")]


def fmt_plain(val, fmt_str):
    """Math-mode number, or the missing marker."""
    return MISSING if val is None else f"${format(val, fmt_str)}$"


def fmt_t_v(r):
    """The t_v cell, value only (``v`` or ``v^{+hi}_{-lo}`` for a Norris measurement) -- no source marker.

    The duration-as-t_v marker lives on the episode name instead (see ``episode_label``), not on this value.
    """
    if r["t_v_source"] == "norris":
        return f"${r['t_v_s']:.3f}^{{+{r['t_v_err_upper_s']:.3f}}}_{{-{r['t_v_err_lower_s']:.3f}}}$"
    return f"${r['t_v_s']:.3f}$"


def fmt_gamma(r, limit, z):
    """One Gamma_min cell: ``value_{-lo}^{+hi}``, or the missing marker."""
    key = f"Gamma_min_{limit}_z{z}"
    if r[key] is None:
        return MISSING
    return f"${round(r[key])}_{{-{r[key + '_err_lower']:.0f}}}^{{+{r[key + '_err_upper']:.0f}}}$"


def episode_label(r, limit):
    """Episode name carrying its markers: $\\dagger$ if t_v is the episode duration, $\\ddagger$ if the LAT detection has TS < 25.

    Markers are on the episode name, not on the values, since each describes the episode as a whole.  The TS < 25 flag
    is Limit A only: it marks a weak photon association, which Limit B does not use.
    """
    marks = []
    if r["t_v_source"] == "duration":
        marks.append(r"\dagger")
    if limit == "A" and r["low_significance"]:
        marks.append(r"\ddagger")
    return r["episode"] + (f"$^{{{','.join(marks)}}}$" if marks else "")


def build_rows(results, limit):
    """Rows grouped by burst (own header row each), episodes in temporal order."""
    results = sorted(results, key=lambda r: (r["GRB"], episode_order(r["episode"])))
    n_cols = 8 if limit == "B" else 9

    rows, current = "", None
    for r in results:
        if r["GRB"] != current:
            if current is not None:
                rows += "    \\midrule\n"
            tex_name = TEX_NAMES[r["GRB"].removeprefix("GRB")]
            rows += f"    \\multicolumn{{{n_cols}}}{{l}}{{\\textbf{{{tex_name}}}}} \\\\\n"
            current = r["GRB"]

        gammas = " & ".join(fmt_gamma(r, limit, z) for z in Z_VALUES)
        lead = f"{fmt_plain(r['E_max_MeV'] / 1e3, '.2f')} & " if limit == "A" else ""
        rows += f"    {episode_label(r, limit)} & {lead}{fmt_t_v(r)} & {fmt_plain(r['beta'], '.3f')} & {gammas} \\\\\n"
    return rows


def tablenotes(rows):
    """Only document markers actually printed in ``rows``."""
    notes = ""
    if r"\dagger" in rows:
        notes += (r"\item[$\dagger$] Episode duration adopted as an upper bound on the variability timescale, "
                  r"giving a conservative $\Gamma_{\min}$." "\n")
    if r"\ddagger" in rows:
        notes += (r"\item[$\ddagger$] \ac{LAT} detection with $\mathrm{TS} < 25$; the highest-energy photon "
                  r"association is not secure." "\n")
    if f" {MISSING} " in rows:
        notes += (r"\item[$-$] Not computed: the episode's best-fit model has no usable high-energy photon index "
                  r"($\beta<-1$ required)." "\n")
    return notes


def build_latex_table(results, limit):
    """Render the Limit-A or Limit-B redshift-sweep table."""
    rows = build_rows(results, limit)
    seed = int(results[0]["seed"])
    if limit == "A":
        label, colspec, n_lead = "tab:lorentz_unknown_z", "lcccccccc", 4
        what = (r"Minimum bulk Lorentz factor $\Gamma_{\min}$ from the gamma-gamma opacity condition (Limit A), "
                r"for the three bursts without a measured redshift, at the assumed redshifts $z=1,2,3,5,7$.")
        cols = r"Episode & $E_{\rm GeV}$ & $t_{\rm v}$ [s] & $\beta$ &"
    else:
        label, colspec, n_lead = "tab:lorentz_limit_b_unknown_z", "lccccccc", 3
        what = (r"Minimum bulk Lorentz factor $\Gamma_{\min}$ from Compton scattering off pair-produced $e^{\pm}$ "
                r"(Limit B), for the three bursts without a measured redshift, at the assumed redshifts $z=1,2,3,5,7$.")
        cols = r"Episode & $t_{\rm v}$ [s] & $\beta$ &"

    blanks = " &" * n_lead
    return (
        "% AUTO-GENERATED by codes-for-paper/lorentz_factor/generate_lorentz_table_unknown_z.py\n"
        "% Do not edit by hand — regenerate from lorentz_results_unknown_z.csv instead.\n"
        r"\begin{table}[!ht]" "\n"
        r"\centering" "\n"
        rf"\caption{{{what}" "\n"
        rf"Errors are the statistical $1\sigma$ interval from $10^{{4}}$ Monte Carlo draws (seed {seed}) "
        r"of the spectral parameters~\citep{Lithwick2001}.}" "\n"
        rf"\label{{{label}}}" "\n"
        r"\resizebox{\columnwidth}{!}{" "\n"
        r"\begin{threeparttable}" "\n"
        r"\renewcommand{\arraystretch}{1.25}" "\n"
        rf"\begin{{tabular}}{{{colspec}}}" "\n"
        r"\toprule" "\n"
        rf"{cols} \multicolumn{{5}}{{c}}{{$\Gamma_{{\min}}$}} \\" "\n"
        rf"{blanks} $z=1$ & $z=2$ & $z=3$ & $z=5$ & $z=7$ \\" "\n"
        r"\midrule" "\n"
        + rows
        + r"\bottomrule" "\n"
        r"\end{tabular}" "\n"
        r"\begin{tablenotes}" "\n"
        r"\footnotesize" "\n"
        + tablenotes(rows)
        + r"\end{tablenotes}" "\n"
        r"\end{threeparttable}" "\n"
        r"}" "\n"
        r"\end{table}" "\n"
    )


def main(limits=("A", "B")):
    for limit in limits:
        csv_path, tex_path = INPUTS[limit]
        tex_path.write_text(build_latex_table(load_results(csv_path), limit))
        print(f"Saved: {tex_path}")


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("A", "B"))
