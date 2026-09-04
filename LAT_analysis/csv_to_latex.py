"""Generate the paper's LAT appendix table from `lat_photons.csv`.

Reference implementation for this convention is
`codes-for-paper/amati_relationship/csv_to_latex.py`: read the CSV, emit the
paper's `\\grb...` macros, never retype numbers into the `.tex`.

Output goes to `lat_info_table.tex` here and is copied into
`GRBResearchPaper/tex_files/generated/`; the two repos are independent, so that
copy is manual (same as the Amati and photospheric tables -- BUGS.md BUG-10).
"""

from pathlib import Path

import pandas as pd

CSV_PATH = Path(__file__).parent / "lat_photons.csv"
OUT_PATH = Path(__file__).parent / "lat_info_table.tex"

# Order the bursts as the paper does, not alphabetically.
GRB_ORDER = ["GRB080916C", "GRB131014A", "GRB140206B", "GRB231129C"]

TEX_NAMES = {
    "GRB080916C": r"\grbzeroeightzeroninesixteenC",
    "GRB131014A": r"\grbthirteentenfourteenA",
    "GRB140206B": r"\grbfourteenzerotwozerosixB",
    "GRB231129C": r"\grbtwentythreeeleventwentynineC",
}

# A LAT detection below this TS is not a secure association; those rows carry a
# footnote giving a flux upper limit. Same threshold as lorentz_factor.py.
TS_SECURE_DETECTION = 25.0

# 95% flux upper limits [ph/cm^2/s] for the TS < 25 episodes.
#
# NOT DERIVABLE from LAT_analysis/: these come from a profile-likelihood upper
# limit computed separately by gtlike's UpperLimits, and that output is not among
# the files this directory holds. They sit 1-8% above (Flux + 2 sigma), so they
# cannot be reconstructed from the fit results either. Keyed on (GRB, episode) so
# a stale entry surfaces as a KeyError rather than a silently wrong number.
FLUX_UPPER_LIMITS = {
    ("GRB080916C", "TR5"): "4.11e-04",
    ("GRB140206B", "TR5"): "4.75e-05",
    ("GRB140206B", "TR6"): "3.75e-05",
    ("GRB231129C", "EX0"): "4.73e-04",
    ("GRB231129C", "TR1"): "5.81e-04",
}


# The table caption, one LaTeX sentence per entry. Emitted one-per-line so the
# generated .tex follows the project's one-sentence-per-line LaTeX rule while
# keeping every line of this file inside the 120-character Python limit.
CAPTION_SENTENCES = [
    r"Time-resolved spectral analysis episodes for four Fermi-GBM gamma-ray bursts.",
    r"For each GRB, we show the episode type, time interval relative to trigger, "
    r"total number of events and high-probability events, characteristic photon "
    r"energy, arrival time, test statistic (TS), and LAT photon index.",
    r"Time-integrated episodes (\tnty) represent the standard burst duration, "
    r"while time-resolved episodes (TR, and EX) enable spectral evolution studies "
    r"within individual bursts.",
]

CAPTION = "".join(f"        {sentence}\n" for sentence in CAPTION_SENTENCES)


def episode_tex(episode: str) -> str:
    """Render an episode label for the appendix table.

    Unified with the EX0/EX1 convention used everywhere else in this project (results.json,
    grb_time.py, log_to_latex_parser.py) -- this table previously wrote the excess episodes as
    EX--A/EX--B, a second naming scheme that only existed here; unified 2026-09-04 at the user's
    direction so the whole paper uses one scheme.
    """
    if episode == "T90":
        return r"\tnty"
    if episode == "EX0":
        return "EX0"
    if episode == "EX1":
        return "EX1"
    if episode.startswith("TR"):
        return f"TR {episode[2:]}"
    raise ValueError(f"Unhandled episode label: {episode}")


def data_cells(row, note_letter=None) -> str:
    """The seven data columns shared by both blocks of the table."""
    ts = f"{row.ts:.3f}"
    if note_letter is not None:
        ts += rf"\tnote{{{note_letter}}}"
    return (
        rf"\sirangeDuration{{{row.t_start_s:.3f}}}{{{row.t_stop_s:.3f}}} & "
        f"{row.n_events} & {row.n_events_high_prob} & "
        f"{row.e_max_MeV:.2f} & {row.t_arr_s:.3f} & {ts} & "
        f"{row.photon_index:.2f}({row.photon_index_err:.2f})"
    )


def build_table(frame: pd.DataFrame) -> str:
    """Render the whole appendix table."""
    notes = {}
    for letter, key in zip("abcdefghij", sorted(FLUX_UPPER_LIMITS, key=lambda k: (GRB_ORDER.index(k[0]), k[1]))):
        notes[key] = letter

    # --- time-integrated block: one T90 row per burst -------------------------
    integrated = ""
    for grb in GRB_ORDER:
        row = next(frame[(frame.grb_name == grb) & (frame.episode == "T90")].itertuples(index=False))
        integrated += f"                {TEX_NAMES[grb]} & {episode_tex('T90')} & {data_cells(row)} \\\\\n"

    # --- time-resolved block: everything else, grouped by burst ---------------
    resolved = ""
    for position, grb in enumerate(GRB_ORDER):
        rows = frame[(frame.grb_name == grb) & (frame.episode != "T90")]
        rows = rows.sort_values(["t_start_s", "t_stop_s"])
        if position:
            resolved += "                \\cmidrule{2-9}\n"
        for offset, row in enumerate(rows.itertuples(index=False)):
            letter = notes.get((grb, row.episode)) if row.ts < TS_SECURE_DETECTION else None
            lead = rf"\multirow{{{len(rows)}}}{{*}}{{{TEX_NAMES[grb]}}} " if offset == 0 else ""
            resolved += f"                {lead}& {episode_tex(row.episode)} & {data_cells(row, letter)} \\\\\n"

    tablenotes = "".join(
        f"                \\item[\\textit{{{letter}}}] \\pFlux{{{FLUX_UPPER_LIMITS[key]}}}\n"
        for key, letter in sorted(notes.items(), key=lambda kv: kv[1])
    )

    return (
        "% AUTO-GENERATED by LAT_analysis/csv_to_latex.py from LAT_analysis/lat_photons.csv\n"
        "% Do not hand-edit: regenerate that script and copy the result here.\n"
        r"""\begin{table*}[!ht]
    \centering
    \caption{
"""
        + CAPTION
        + r"""    }
    \label{tab:burst_table}
    \resizebox{\textwidth}{!}{
        \begin{threeparttable}
            \begin{tabular}{
                    @{}lll
                    ccccccc@{}
                }
                \toprule
                \multirow{2.5}{*}{GRB name} &
                \multirow{2.5}{*}{Ep. Type} &
                \multirow{2.5}{*}{Time interval (s)} &
                \multicolumn{2}{c}{Events} &
                \multirow{2.5}{*}{\twoRs{Photon energy}{(MeV)}} &
                \multirow{2.5}{*}{\twoRs{Arrival time}{(s)}} &
                \multirow{2.5}{*}{TS} &
                \multirow{2.5}{*}{Index} \\
                \cmidrule(lr){4-5}
                &&& {Total} & {High Prob.} & & & & \\
                \midrule
                \multicolumn{9}{l}{\textbf{Time-integrated}} \\
"""
        + integrated
        + "                \\midrule\n"
        + "                \\multicolumn{9}{l}{\\textbf{Time-resolved}}\\\\\n"
        + resolved
        + r"""                \bottomrule
            \end{tabular}
            \begin{tablenotes}[para, flushleft]
"""
        + tablenotes
        + r"""            \end{tablenotes}
        \end{threeparttable}
    }
\end{table*}
"""
    )


def main():
    """Read the CSV and write the appendix table."""
    frame = pd.read_csv(CSV_PATH)
    missing = set(FLUX_UPPER_LIMITS) - set(zip(frame.grb_name, frame.episode))
    if missing:
        raise KeyError(f"FLUX_UPPER_LIMITS refers to episodes not in the CSV: {sorted(missing)}")

    weak = {(r.grb_name, r.episode) for r in frame.itertuples(index=False) if r.ts < TS_SECURE_DETECTION}
    if weak != set(FLUX_UPPER_LIMITS):
        raise KeyError(
            "TS < 25 episodes and FLUX_UPPER_LIMITS disagree.\n"
            f"  missing a limit: {sorted(weak - set(FLUX_UPPER_LIMITS))}\n"
            f"  limit no longer needed: {sorted(set(FLUX_UPPER_LIMITS) - weak)}"
        )

    OUT_PATH.write_text(build_table(frame))
    print(f"Saved: {OUT_PATH}  ({len(frame)} episodes, {len(weak)} with TS < {TS_SECURE_DETECTION:.0f})")


if __name__ == "__main__":
    main()
