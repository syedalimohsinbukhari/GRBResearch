# GBM fluence and flux — method notes

This folder previously had two scripts (`grb_fluence.py`, `flux_energyFlux.py`) computing photon
flux/fluence and energy flux respectively, with no method note — a pre-existing gap, not one
introduced here. This file covers `grb_fluence.py` and its two outputs, `flux_fluence.csv` and
`flux_energy_flux.csv`; `flux_energyFlux.py`/`flux_energy_flux.py` (the flux-vs-energy-flux
scatter plot) is unaffected by the 2026-09-01 fix below and not otherwise documented here.

## 1. What is computed

For each burst's BEST-fit models, `grb_fluence.py` uses `FluxFluenceCalculator`
(`src/grb_research/grb_calculations.py`) to Monte Carlo sample the fitted spectral parameters and
compute:

- **Photon flux** $F = \int_{E_\text{min}}^{E_\text{max}} N(E)\,dE$ [ph cm$^{-2}$ s$^{-1}$], the
  fitted photon spectrum integrated over energy.
- **Fluence** $S = \int_{E_\text{min}}^{E_\text{max}} N(E)\,E\,dE \times \Delta t$ [erg cm$^{-2}$],
  the energy flux integrated over both energy and the episode duration $\Delta t$
  (`calculate("fluence", in_ergs=True)`; the `energy_flux=True` path used for
  `flux_energy_flux.csv` instead sets $\Delta t = 1$, giving an energy flux rather than a
  time-integrated fluence despite the shared column name — a pre-existing naming quirk, not
  changed here since only `flux_fluence.csv` is used by the paper's fluence table).

Both are purely observer-frame quantities from the best-fit model; no redshift or cosmology
enters either calculation.

## 2. Decisions and who made them

### 2.1 Energy band: 8 keV–40 MeV, not the calculator's 10 keV–1 MeV default — *Claude, fixing a bug found 2026-09-01*

`FluxFluenceCalculator.__init__` defaults `log_energy_range=(1, 3)`, i.e. 10 keV–1 MeV. Before
2026-09-01, `grb_fluence.py` never overrode this, so `flux_fluence.csv` silently reported fluence
over a narrower band than the paper's own stated sample-selection criterion — "High GBM fluence,
ensuring well-constrained spectral fits across the full 8 keV–40 MeV energy range"
(`section-1-introduction.tex`). Found while building a fluence table to substantiate that
criterion (`grb_paper_weaknesses_and_fixes.md` Priority 3, "bright framing with no fluence
table"): the table would have shown numbers for the wrong band. Fixed by passing
`log_energy_range=(log10(8), log10(40000))` explicitly, matching the criterion's own stated
range. See `review-resolution.md` Priority 3 item 8.

### 2.2 No redshift/cosmology columns — *Claude, consequence of the definition*

Both flux and fluence here are integrated in the observed frame from the fitted (observer-frame)
model — the same reasoning `bb_fraction.md` §2.2 gives for its observer-frame $f_\text{BB}^\text{obs}$
column. CLAUDE.md's CSV convention calls for recording "the swept value, when a z-range is used
instead of a fixed z" — not applicable here, since no z enters the calculation at all, so the
column is omitted rather than filled with a placeholder.

### 2.3 T90-only in the paper table — *Claude*

`flux_fluence.csv` retains one row per BEST-fit model per episode (all 26 across the four
bursts), consistent with every other CSV in this project. The generated LaTeX table
(`csv_to_latex.py`) filters to `ep_type == "T90"` — one row per burst — since the table exists
only to substantiate the "high GBM fluence" selection criterion for the sample as a whole, not to
report a per-episode fluence breakdown (nothing in the paper currently needs the latter).

### 2.4 Uncertainties: Monte Carlo, 10 000 samples — *existing project convention*

Sample count unchanged from before this fix; matches every other CSV in `codes-for-paper/`. **Updated in the RNG/seeding overhaul (Phase B, 2026-09-02):** the literal seed `12345` was replaced with `seed_from_name(__file__)`, a per-file value deterministically derived from the project master seed (`src/grb_research/SEEDING.md`), decorrelating this script's draws from every other script's. This file's own `rng.spawn(2)` pattern (independent child streams for the two flux/fluence columns) was already correct and needed no other change — it is the reference pattern the rest of Phase B was ported from.

## 3. Validation

- **Directional sanity check**: widening the band from 10 keV–1 MeV to 8 keV–40 MeV increased
  every burst's T90 fluence, as expected for a strictly wider integration range over a
  monotonically-supported spectrum.
- **Literature rank-order cross-check — user, 2026-09-01.** Checked the four bursts against the
  Fermi LAT GRB catalog (the catalog this sample was originally drawn from). The catalog's GBM
  fluences are not numerically identical to `tab:fluence`'s values (expected — different band,
  and the catalog is a standard on-board/ground pipeline value rather than this project's own
  RMFIT best-fit integration) but the **rank order matches exactly** for the three bursts present
  in that catalog: GRB231129C > GRB131014A > GRB140206B > GRB080916C in `tab:fluence`, and
  GRB131014A > GRB140206B > GRB080916C for the same three in the LAT catalog (GRB231129C isn't in
  that catalog at all — it's the more recent of the four bursts). GRB231129C also being the
  fluence-brightest burst here is independently consistent with it being the burst that isn't yet
  catalogued, i.e. plausibly a more recent/less-processed detection rather than a computational
  artifact. This is a real independent check, not just an internal consistency check.

## 4. Limitations and open questions

- **No catalogue cross-check for GRB231129C.** GRB231129C is absent from the Fermi LAT GRB
  catalog used for the rank-order check above, so its fluence has no independent cross-check at
  all — only the other three bursts got one, and only at the rank-order level, not a numeric
  comparison (see §3).
- **`flux_energy_flux.csv`'s "fluence" columns are actually an energy flux** (duration fixed to
  1 s via `energy_flux=True`), despite sharing column names with `flux_fluence.csv`. Not
  currently used by the paper table; flagged here so a future reader doesn't conflate the two
  files.
