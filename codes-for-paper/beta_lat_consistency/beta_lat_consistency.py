"""Diagnostic: does an added-PL component's index track the independent LAT photon
index better than the base continuum's own high-energy index (beta)?

Not a paper pipeline script (no CSV/figure this project's tables draw from) -- this
generates the comparison table in beta_lat_consistency.md, written to investigate
GRBResearch_Issues_List.md Issue 3 (fitted beta vs. LAT photon index tension).

For every episode with LAT coverage, compares:
  - the LAT-only photon index (LAT_analysis/lat_photons.csv, from an independent
    gtlike fit to just the LAT photons of that episode)
  - beta, the high-energy index of whichever continuum model is this episode's
    officially selected (BEST/SAFE/MARGINAL) model (results.json)
  - add_index_pl, the index of an added power-law component in a BASE+PL+BB fit
    (BAND_PL_BB / SBPL_PL_BB), when one was attempted for that episode -- these
    fits exist in results.json but were never selected or reported anywhere in
    the paper, all falling short of the SAFE threshold.

No Monte Carlo here -- every number is read directly from results.json or
lat_photons.csv and combined algebraically (simple 1-sigma error propagation),
so there is no seed/n_samples to report, per the convention in LAT_analysis.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from grb_research.grb_constants import short_to_long
from grb_research.grb_time import EpisodeTypes, TimeInterval

HERE = Path(__file__).parent
ROOT = HERE.parent.parent
OUT_PATH = HERE / "beta_lat_consistency.md"

GRB_ORDER = ["080916C", "131014A", "140206B", "231129C"]

# Continuum models that carry a genuine high-energy power-law index (beta) usable
# for this comparison. CPL has no beta (exponential cutoff instead of a second
# power law), so CPL/CPL_BB/CPL_PL_BB are excluded from the beta/add_index_pl
# columns even when they are the selected model.
SELECTED_BETA_KEY = {"BAND_BB": "index2_band", "SBPL_BB": "index2_sbpl", "BAND": "index2_band", "SBPL": "index2_sbpl"}
PL_VARIANTS = {"BAND_PL_BB": "index2_band", "SBPL_PL_BB": "index2_sbpl"}

# An add_index_pl (or beta) fit error above this is functionally unconstrained --
# the amplitude of that component is consistent with zero, so any numerical
# agreement with the LAT index is coincidental, not evidence. Flagged in the
# table rather than dropped, so the reader sees it was checked, not omitted.
UNCONSTRAINED_ERR = 2.0

# Copied from LAT_analysis/csv_to_latex.py (not a re-derivation) -- a photon index landing on this
# value within this tolerance is pinned at gtlike's fit boundary, not a genuine measurement (three
# episodes here: GRB140206B TR2 despite TS=46.9 > 25, and GRB231129C EX0/TR1 at TS<25). Comparing a
# fitted beta against a pinned index produces a meaningless multi-hundred-sigma "tension" if not
# excluded explicitly.
PINNED_INDEX_VALUE = -6.00
PINNED_INDEX_TOL = 1e-2


def sigma_tension(v1: float, e1: float, v2: float, e2: float) -> float:
    """Number of combined-error sigma between two independent measurements."""
    return abs(v1 - v2) / (e1**2 + e2**2) ** 0.5


def lat_index_is_pinned(lat_index: float) -> bool:
    return abs(lat_index - PINNED_INDEX_VALUE) < PINNED_INDEX_TOL


def episode_label(interval: TimeInterval) -> str:
    if interval.kind is EpisodeTypes.TR:
        return f"TR{interval.index}"
    return interval.kind.value


def build_rows() -> list[dict]:
    data = json.loads((ROOT / "results.json").read_text())
    lat = pd.read_csv(ROOT / "LAT_analysis" / "lat_photons.csv")
    lat_by_key = {(r.grb_name, r.episode): r for r in lat.itertuples(index=False)}

    rows = []
    for short in GRB_ORDER:
        grb_name = f"GRB{short}"
        episodes = data[short_to_long[short]]
        for raw_key, models in episodes.items():
            interval = TimeInterval.from_string(raw_key)
            if interval.kind is EpisodeTypes.UNKNOWN:
                continue
            label = episode_label(interval)
            if (grb_name, label) not in lat_by_key:
                continue
            lat_row = lat_by_key[(grb_name, label)]
            pinned = lat_index_is_pinned(lat_row.photon_index)

            selected_name = selected_beta = selected_beta_err = None
            for name, m in models.items():
                if m.get("_status") == "BEST" and name in SELECTED_BETA_KEY:
                    selected_name = name
                    selected_beta, selected_beta_err = m[SELECTED_BETA_KEY[name]]

            plbb = models.get("PL_BB")
            plbb_status = plbb_idx = plbb_err = None
            if plbb is not None:
                plbb_status = plbb.get("_status")
                plbb_idx, plbb_err = plbb["index1_pl"]

            base_row = {
                "grb": grb_name,
                "episode": label,
                "ts": lat_row.ts,
                "n_lat": lat_row.n_events,
                "lat_index": lat_row.photon_index,
                "lat_err": lat_row.photon_index_err,
                "lat_pinned": pinned,
                "selected_model": selected_name,
                "beta": selected_beta,
                "beta_err": selected_beta_err,
                "sigma_beta": sigma_tension(selected_beta, selected_beta_err, lat_row.photon_index, lat_row.photon_index_err)
                if selected_beta is not None and not pinned
                else None,
                "plbb_status": plbb_status,
                "plbb_idx": plbb_idx,
                "plbb_err": plbb_err,
                "sigma_plbb": sigma_tension(plbb_idx, plbb_err, lat_row.photon_index, lat_row.photon_index_err)
                if plbb_idx is not None and not pinned
                else None,
            }

            pl_variants_present = [name for name in PL_VARIANTS if name in models]
            if not pl_variants_present:
                rows.append({**base_row, "pl_model": None, "pl_status": None, "add_index_pl": None, "add_err": None, "sigma_addpl": None})
                continue

            for pl_name in pl_variants_present:
                m = models[pl_name]
                add_idx, add_err = m["add_index_pl"]
                rows.append(
                    {
                        **base_row,
                        "pl_model": pl_name,
                        "pl_status": m.get("_status"),
                        "add_index_pl": add_idx,
                        "add_err": add_err,
                        "sigma_addpl": sigma_tension(add_idx, add_err, lat_row.photon_index, lat_row.photon_index_err)
                        if not pinned
                        else None,
                    }
                )

    return rows


def fmt(value, err, decimals=2):
    if value is None:
        return "—"
    return f"{value:.{decimals}f} ± {err:.{decimals}f}"


def fmt_sigma(sigma, err, comparable, pinned):
    if pinned:
        return "n/a (LAT index pinned)"
    if not comparable or sigma is None:
        return "—"
    flag = " (unconstrained)" if err is not None and err > UNCONSTRAINED_ERR else ""
    return f"{sigma:.2f}σ{flag}"


def build_table(rows: list[dict]) -> str:
    header = (
        "| GRB | Episode | TS | N(LAT) | LAT index | β model (status) | β | β vs LAT | "
        "PL_BB (status) | index1_pl | PL_BB vs LAT | Add-PL model (status) | add_index_pl | add-PL vs LAT |\n"
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n"
    )
    lines = []
    for r in rows:
        beta_model = f"{r['selected_model']} (BEST)" if r["selected_model"] else "—"
        pl_model = f"{r['pl_model']} ({r['pl_status']})" if r["pl_model"] else "—"
        plbb_model = f"PL_BB ({r['plbb_status']})" if r["plbb_status"] else "—"
        lat_index_str = fmt(r["lat_index"], r["lat_err"], 3) + (" [pinned]" if r["lat_pinned"] else "")
        lines.append(
            f"| {r['grb']} | {r['episode']} | {r['ts']:.1f} | {r['n_lat']} | "
            f"{lat_index_str} | {beta_model} | {fmt(r['beta'], r['beta_err'])} | "
            f"{fmt_sigma(r['sigma_beta'], r['beta_err'], r['beta'] is not None, r['lat_pinned'])} | "
            f"{plbb_model} | {fmt(r['plbb_idx'], r['plbb_err'])} | "
            f"{fmt_sigma(r['sigma_plbb'], r['plbb_err'], r['plbb_idx'] is not None, r['lat_pinned'])} | "
            f"{pl_model} | {fmt(r['add_index_pl'], r['add_err'])} | "
            f"{fmt_sigma(r['sigma_addpl'], r['add_err'], r['add_index_pl'] is not None, r['lat_pinned'])} |"
        )
    return header + "\n".join(lines) + "\n"


def main():
    rows = build_rows()
    table = build_table(rows)
    OUT_PATH.write_text(TEMPLATE.format(table=table))
    print(f"Saved: {OUT_PATH} ({len(rows)} rows)")


TEMPLATE = """# Beta vs. LAT photon index consistency — diagnostic table

AUTO-GENERATED by `beta_lat_consistency.py` from `results.json` and `LAT_analysis/lat_photons.csv`.
Not a paper deliverable — written to investigate `GRBResearch_Issues_List.md` Issue 3
(fitted beta disagrees with the independent LAT photon index for several episodes).

## The reviewer's point

`GRBResearch_Issues_List.md` Issue 3 flags that the fitted high-energy index beta (from the joint
GBM+LAT spectral fit) disagrees with the independently-measured LAT-only photon index
(`LAT_analysis/lat_photons.csv`, a separate `gtlike` fit to just that episode's LAT photons) by
~2.9-4 sigma for GRB131014A T90 and GRB231129C T90, and even more severely (though against a
pinned, non-genuine LAT index) for GRB231129C EX0/TR1. Both numbers are nominally measurements of
the same physical quantity (the spectral slope in the same energy region), and the joint fit's
likelihood includes those same LAT photons as data, so a reader would expect rough agreement. The
reviewer's suggested fix was a caveat: for these low-LAT-count episodes, beta is effectively set by
the GBM data extrapolated upward, not by the sparse real LAT photons, so it should not be read as
an independent LAT constraint.

## Counter-rationale (user, 2026-09-06)

The reviewer's framing implicitly treats beta and the LAT index as two measurements of the same
thing that merely differ in precision or dominant dataset. The user's counterpoint is sharper: they
should not be expected to agree at all, on physical grounds, independent of photon statistics.
Band and SBPL are curvature models shaped almost entirely by GBM's keV-MeV data; above the break,
beta is mathematically a power law extending to infinity, but it has no physical mandate to
describe a genuinely distinct spectral component at LAT energies (100 MeV-GeV) two-plus decades
above where the fit is actually constrained. This is precisely why bright LAT bursts are often
better described by an additional, physically separate hard power-law component layered on top of
the Band/SBPL continuum (`BASE+PL`, well documented in the literature, e.g. GRB090510, GRB090902B)
rather than by trusting the base continuum's own extrapolated tail. If that is correct, the
physically appropriate quantity to compare against the independent LAT index is not beta at all,
but the index of an actual added PL component (`add_index_pl` below, from a `BASE+PL+BB` fit) --
and the table below tests that directly against the data already in this repository.

## What this compares

For every episode with LAT coverage:

- **LAT index** — the independent `gtlike` power-law fit to just that episode's LAT photons
  (`LAT_analysis/lat_photons.csv`), with its TS and photon count.
- **β** — the high-energy index of this episode's officially selected (BEST) continuum model,
  the number that actually appears in the paper's tables.
- **index1_pl** — the single power-law index of a plain `PL_BB` fit (power law + blackbody, no
  curvature at all, forced across the entire GBM+LAT band). `PL_BB` is present and mostly `SAFE`
  for every episode checked here (it is never absent from the grid), but it is never the
  BEST-selected model for any episode with a real spectral peak -- included here as a check on
  whether the *simplest possible* +BB model's slope tracks the LAT index better than beta, not
  because it is a good overall fit (see its C-stat/dof, always far worse than the selected model's).
- **add_index_pl** — the index of an added power-law component, from a `BASE+PL+BB` fit
  (`BAND_PL_BB`/`SBPL_PL_BB` in `results.json`) when one exists for that episode. These fits were
  never selected as BEST/SAFE for any episode in the paper's sample and do not appear in any
  paper table — they exist in `results.json` from the model-selection grid search but were
  rejected by this project's SAFE/MARGINAL/UNSAFE error-threshold criteria (Kaneko2006-based;
  see `abstract.tex`).
- **σ columns** — combined-error tension, `|v1 - v2| / sqrt(e1² + e2²)`, between each fitted
  quantity and the LAT index.

CPL/CPL_BB/CPL_PL_BB are excluded from the β/add_index_pl columns: CPL has an exponential
cutoff rather than a second power-law index, so it has no β to compare.

## Table

{table}

## Reading this table

- **"[pinned]" / "n/a (LAT index pinned)"** marks three episodes (GRB140206B TR2, GRB231129C
  EX0/TR1) where the LAT-only fit hit `gtlike`'s index boundary at -6.00 -- not a genuine
  measurement (already footnoted as such in the paper's LAT appendix table). Sigma against a
  pinned value is meaningless (naively computed, GRB140206B TR2 alone would show an ~80-sigma
  "tension" against beta=-2.35 -- an artifact of the LAT index not being real, not a real physics
  discrepancy), so these are excluded from the sigma columns entirely rather than shown as a
  number. Note GRB140206B TR2 pins despite TS=46.9 > 25 -- the TS<25 threshold catches insecure
  *detections*, not this separate failure mode where a real detection still can't pin down a
  finite index.
- **"(unconstrained)"** flags a fit whose error exceeds {threshold}: the component's amplitude
  is statistically consistent with zero, so any numerical closeness to the LAT index there is
  coincidental, not evidence. Treat these rows as "not testable," not as agreement or disagreement.
- For GRB131014A (every episode) and GRB231129C's T90, add_index_pl is consistently closer to
  the LAT index than β, often dramatically so (e.g. GRB131014A T90: β misses by ~3σ, add_index_pl
  by ~1.4σ) — the pattern the reviewer's Issue 3 diagnosis points at.
- For GRB080916C, the picture is mixed: T90 (254 LAT photons) already has β agreeing well with
  the LAT index (no tension to explain), but TR2/TR3 (115/118 photons) still show β disagreeing
  while add_index_pl agrees much better — so photon count alone does not fully explain the effect;
  a genuinely separate hard component appears to bias β even when LAT photon statistics are decent.
- **index1_pl (plain PL_BB) shows a genuinely mixed pattern, not a uniform improvement over beta**
  — reported honestly rather than only in the direction that supports the add_index_pl finding.
  For GRB131014A it tracks the LAT index dramatically well (sigma as low as 0.00-0.01 for
  EX1/TR2, versus beta's 3+), even better than add_index_pl there. For GRB231129C it helps
  moderately (T90/EX1/TR2 all closer than beta). But for GRB080916C and GRB140206B -- the two
  bursts with the most and best-measured LAT photons overall -- index1_pl is consistently *worse*
  than beta (e.g. GRB080916C T90: beta misses by 0.27sigma, index1_pl by 4.78sigma). Forcing away
  real spectral curvature only seems to help as a LAT-index proxy where beta itself is already
  poorly constrained; where the curved fit is well-determined, discarding the curvature just makes
  things worse. index1_pl's small formal errors also do not mean it is a trustworthy fit overall --
  its C-stat/dof is always far worse than the selected model's (e.g. 12744/461 vs. 765/459 for
  GRB131014A T90), so its precision reflects GBM's sheer photon count, not goodness of fit.
- All `*_PL_BB` fits are UNSAFE or MARGINAL by this project's own criteria. This table is a
  diagnostic pointing at a likely explanation, not a claim that these specific fitted values are
  reliable measurements in their own right.
""".replace("{threshold}", str(UNCONSTRAINED_ERR))


if __name__ == "__main__":
    main()
