"""
Lorentz-factor LaTeX table generator (Limit A)
===============================================
Pure CSV-in, tex-out: reads ``lorentz_results.csv`` (written by ``lorentz_factor.py``)
and renders ``lorentz_table.tex``. Has no dependency on ``grb_research``/``results.json``
and does not recompute anything, so it can be rerun on its own -- e.g. to pick up a
table-formatting change -- without repeating the Gamma_min Monte Carlo.

Run standalone with ``python generate_lorentz_table.py`` from this directory, or via
``lorentz_factor.py``'s ``main()``, which shells out to this script by subprocess after
writing ``lorentz_results.csv`` so the two stay in sync without importing each other's
internals directly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
RESULTS_CSV = HERE / "lorentz_results.csv"
TABLE_TEX = HERE / "lorentz_table.tex"

# Duplicated from lorentz_factor.py rather than imported -- this project's "copy rather
# than fight sys.path" convention (see e.g. this same file's own load_norris_tv() comment).
TEX_NAMES = {
    "080916C": r"\grbzeroeightzeroninesixteenC",
    "131014A": r"\grbthirteentenfourteenA",
    "140206B": r"\grbfourteenzerotwozerosixB",
    "231129C": r"\grbtwentythreeeleventwentynineC",
}

# t_v source -> table marker. "literature" carries none today (VARIABILITY_TIMESCALE in
# lorentz_factor.py is empty), but is included so a future entry there doesn't need this
# map touched again.
T_V_SOURCE_MARKER = {"duration": r"\dagger", "norris": r"\ast", "literature": ""}


def load_results(csv_path=RESULTS_CSV):
    """Read lorentz_results.csv into the list-of-dicts shape build_latex_table() expects.

    pandas reads an empty/None cell back as NaN, not None -- fmt()/fmt_t_v() check `is None`,
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
    if r["t_v_source"] == "norris":
        return f"${r['t_v_s']:.3f}^{{+{r['t_v_err_upper_s']:.3f}}}_{{-{r['t_v_err_lower_s']:.3f}}}\\,{marker}$"
    return f"${r['t_v_s']:.3f}" + (f"\\,{marker}" if marker else "") + "$"


def build_latex_table(results):
    """Render the per-episode Gamma_min table.

    Bursts without a redshift never yield a Gamma_min (z enters d_L in tau_hat), so
    they are dropped from the table entirely rather than shown as all-ellipsis rows
    -- the CSV keeps every row regardless, this is a table-only presentation choice.
    """
    results = [r for r in results if r["z"] is not None]
    seed = int(results[0]["seed"])

    rows = ""
    current = None
    for r in results:
        if r["GRB"] != current:
            if current is not None:
                rows += "    \\midrule\n"
            rows += f"    \\multicolumn{{7}}{{l}}{{\\textbf{{{r['tex_name']}}}}} \\\\\n"
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
            f"    {r['episode']} & {fmt(r['z'], '.2f')} & {fmt(r['E_max_MeV'] / 1e3, '.2f')} & "
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
\caption{Minimum bulk Lorentz factor $\Gamma_{\min}$ from the gamma-gamma opacity
condition, evaluated for every episode with \ac{LAT} coverage.
$E_{\rm GeV}$ is the highest-energy \ac{LAT} photon of that episode,
$t_{\rm arr}$ its arrival time relative to $T_0$,
$t_{\rm v}$ the variability timescale,
$\beta$ the high-energy photon index of the episode's best-fit model,
and $\Gamma_{\min}$ the derived lower limit~\citep{Lithwick2001}.
Errors on $\Gamma_{\min}$ are the statistical $1\sigma$ interval from $10^{4}$ Monte Carlo
draws (seed __SEED__) of the spectral parameters (and, for a $\dagger$-free $t_{\rm v}$, of
$t_{\rm v}$'s own measured uncertainty as well); they are far smaller than the systematic
uncertainty from the analytic approximation adopted, and should not be read as the total
uncertainty.
Only \grbzeroeightzeroninesixteenC\ has a confirmed spectroscopic redshift ($z = 4.35$);
the other three bursts lack a measured redshift, so $\Gamma_{\min}$ is undetermined for
them and they are omitted from this table.}
\label{tab:lorentz}
\resizebox{\columnwidth}{!}{
\begin{threeparttable}
\renewcommand{\arraystretch}{1.25}
\begin{tabular}{lcccccc}
\toprule
Episode & $z$ & $E_{\rm GeV}$ [GeV] & $t_{\rm arr}$ [s] &
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
