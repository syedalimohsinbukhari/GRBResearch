"""
Lorentz factor LaTeX table generator (Limit A)
===============================================
Pure CSV-in, tex-out: reads ``lorentz_results.csv`` (written by ``lorentz_factor.py``) and renders ``lorentz_table.tex``.
Has no dependency on ``grb_research``/``results.json`` and does not recompute anything, so it can be rerun on its own
without repeating the Gamma_min Monte Carlo.

Run standalone with ``python generate_lorentz_table.py`` from this directory, or via ``lorentz_factor.py``'s ``main()``,
which shells out to this script by subprocess after writing ``lorentz_results.csv`` so the two stay in sync without
importing each other's internals directly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
RESULTS_CSV = HERE / "lorentz_results.csv"
TABLE_TEX = HERE / "lorentz_table.tex"

TEX_NAMES = {
    "080916C": r"\grbzeroeightzeroninesixteenC",
    "131014A": r"\grbthirteentenfourteenA",
    "140206B": r"\grbfourteenzerotwozerosixB",
    "231129C": r"\grbtwentythreeeleventwentynineC",
}


def episode_order(label: str) -> int:
    """Sort key putting episode labels in temporal reading order: T90, EX0, TR1..TRn, EX1.

    Duplicated from grb_research.episode_order (the canonical definition, also used by
    gamma_comparison_plot.py, pe_er_photosphere.py, amati_relationship.py) rather than imported,
    same as TEX_NAMES above -- this generator is deliberately grb_research-free so it stays
    rerunnable from just the CSV. Keep in sync if the canonical version ever changes.
    """
    fixed = {"T90": 0, "EX0": 1, "EX1": 90}
    if label in fixed:
        return fixed[label]
    if label.startswith("TR"):
        return 10 + int(label[2:])
    return 99

# Norris-measured t_v is the default method now, not a marked exception -- it's explained in the
# main text instead. Only duration (the fallback, still worth flagging inline) gets a marker.
T_V_SOURCE_MARKER = {"duration": r"\dagger", "norris": "", "literature": ""}


def load_results(csv_path=RESULTS_CSV):
    """Read lorentz_results.csv into the list-of-dicts shape build_latex_table() expects.

    `pandas` reads an empty/None cell back as NaN, not None -- fmt()/fmt_t_v() check `is None`,
    so NaN has to be converted back explicitly or every such cell would render as a stray
    number instead of the intended ``\\ldots``/bare-duration formatting.
    """
    df = pd.read_csv(csv_path)
    df["tex_name"] = df["GRB"].map(lambda g: TEX_NAMES[g.removeprefix("GRB")])
    return [{k: (None if pd.isna(v) else v) for k, v in row.items()} for row in df.to_dict("records")]


def fmt(val, fmt_str):
    """Format a value for a math-mode table cell, or an unset marker if it is None."""
    return r"\ldots" if val is None else f"${format(val, fmt_str)}$"


def fmt_t_v(r):
    """Format the t_v cell: a bare value with its source marker, or value +- error for a measured one.

    Per-cell, not per-column -- t_v_source varies row to row (duration vs. a measured Norris
    value), so a single header dagger can no longer describe every row.
    """
    if r["t_v_s"] is None:
        return r"\ldots"

    marker = T_V_SOURCE_MARKER[r["t_v_source"]]
    # `\,^{marker}` (empty-nucleus trick), not a bare trailing `^{marker}` -- the norris branch's
    # base already carries its own ^{+err}_{-err}, and a second `^{...}` glued straight onto that
    # same atom is a LaTeX "Double superscript" error. The thin space starts a fresh atom for the
    # marker to attach to instead. Same pattern build_latex_table() already uses for `\,^{\ddagger}`.
    marker_str = f"\\,^{{{marker}}}" if marker else ""
    if r["t_v_source"] == "norris":
        return f"${r['t_v_s']:.3f}^{{+{r['t_v_err_upper_s']:.3f}}}_{{-{r['t_v_err_lower_s']:.3f}}}{marker_str}$"
    return f"${r['t_v_s']:.3f}{marker_str}$"


def build_latex_table(results):
    """Render the per-episode Gamma_min table."""
    results = [r for r in results if r["z"] is not None]
    # Sort within each GRB by temporal episode order -- results.json's own row order isn't
    # guaranteed chronological (this is what let EX1 print before TR5 in an earlier draft of
    # amati_relationship.py's table; not currently visible here since GRB080916C's Limit A/B
    # tables happen to already list T90, EX0, TR1, TR2, TR3 in order, but not something to rely on).
    results = sorted(results, key=lambda r: (r["GRB"], episode_order(r["episode"])))
    seed = int(results[0]["seed"])

    rows = ""
    current = None
    for r in results:
        if r["GRB"] != current:
            if current is not None:
                rows += "    \\midrule\n"
            rows += f"    \\multicolumn{{6}}{{l}}{{\\textbf{{{r['tex_name']}}}}} \\\\\n"
            current = r["GRB"]

        gamma = None if r["Gamma_min"] is None else round(r["Gamma_min"])
        # Markers go inside the math group, not appended as a second one.
        if gamma is None:
            gamma_str = r"\ldots"
        else:
            gamma_str = f"${gamma}_{{-{r['Gamma_min_err_lower']:.0f}}}^{{+{r['Gamma_min_err_upper']:.0f}}}"
            if r["low_significance"]:
                gamma_str += r"\,^{\ddagger}"
            gamma_str += "$"

        rows += (
            f"    {r['episode']} & {fmt(r['E_max_MeV'] / 1e3, '.2f')} & "
            f"{fmt(r['t_arr_s'], '.2f')} & {fmt_t_v(r)} & {fmt(r['beta'], '.3f')} & {gamma_str} \\\\\n"
        )

    # Only document a marker actually printed in `rows` -- e.g. \ddagger (low-significance
    # detection) currently never fires for GRB080916C's surviving 5-episode table (the
    # >1 GeV photon floor, lorentz_factor.md Sec. 12, already dropped its one low-significance
    # row, TR5), and a stale tablenote for a symbol that isn't in the table is worse than no
    # tablenote at all.
    tablenote_items = ""
    if r"\dagger" in rows:
        tablenote_items += r"\item[$\dagger$] Episode duration adopted as an upper bound on the variability timescale, giving a conservative $\Gamma_{\min}$." "\n"
    if r"\ast" in rows:
        tablenote_items += r"\item[$\ast$] Measured variability timescale from a joint Norris-pulse fit to the light curve, given with its $1\sigma$ Monte Carlo uncertainty." "\n"
    if r"\ddagger" in rows:
        tablenote_items += r"\item[$\ddagger$] \ac{LAT} detection with $\mathrm{TS} < 25$; the highest-energy photon association is not secure." "\n"

    return (
        "% AUTO-GENERATED by codes-for-paper/lorentz_factor/generate_lorentz_table.py\n"
        "% Do not edit by hand — regenerate from lorentz_results.csv instead.\n"
        r"""\begin{table}[!ht]
\centering
\caption{Minimum bulk Lorentz factor $\Gamma_{\min}$ from the gamma-gamma opacity condition, evaluated for every episode with \ac{LAT} coverage.
$E_{\rm GeV}$ is the highest-energy \ac{LAT} photon of that episode, $t_{\rm arr}$ its arrival time relative to $T_0$, $t_{\rm v}$ the variability timescale, $\beta$ the high-energy photon index of the episode's best-fit model, and $\Gamma_{\min}$ the derived lower limit~\citep{Lithwick2001}.
% Errors on $\Gamma_{\min}$ are the statistical $1\sigma$ interval from $10^{4}$ Monte Carlo draws (seed __SEED__) of the spectral parameters (and, for a $\dagger$-free $t_{\rm v}$, of $t_{\rm v}$'s own measured uncertainty as well); they are far smaller than the systematic uncertainty from the analytic approximation adopted, and should not be read as the total uncertainty.
Errors in the table are the statistical $1\sigma$ interval from $10^{4}$ Monte Carlo draws (seed __SEED__) of the spectral parameters;
the redshift of \grbzeroeightzeroninesixteenC\ is $z = 4.35$.
% Only \grbzeroeightzeroninesixteenC\ has a confirmed spectroscopic redshift ($z = 4.35$); the other three bursts lack a measured redshift, so $\Gamma_{\min}$ is undetermined for them and they are omitted from this table.
}
\label{tab:lorentz}
\resizebox{\columnwidth}{!}{
\begin{threeparttable}
\renewcommand{\arraystretch}{1.25}
\begin{tabular}{lccccc}
\toprule
Episode & $E_{\rm GeV}$ [GeV] & $t_{\rm arr}$ [s] &
    $t_{\rm v}$ [s] & $\beta$ & $\Gamma_{\min}$ \\
\midrule
"""
        + rows
        + r"""\bottomrule
\end{tabular}
\begin{tablenotes}
\footnotesize
"""
        + tablenote_items
        + r"""\end{tablenotes}
\end{threeparttable}
}
\end{table}
"""
    ).replace("__SEED__", str(seed))


def main():
    results = load_results()
    TABLE_TEX.write_text(build_latex_table(results))
    print(f"Saved: {TABLE_TEX}")


if __name__ == "__main__":
    main()
