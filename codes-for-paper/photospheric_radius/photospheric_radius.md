# Photospheric radius and bulk Lorentz factor — method notes

Companion to `pe_er_photosphere.py`. Phase 2 of `PLAN.md`. Produced 2026-08-21.

---

## 1. What is computed

From the observed temperature and thermal flux of the blackbody component, following Pe'er, Ryde, Wijers, Mészáros & Rees (2007), ApJ **664**, L1:

| symbol | meaning |
|---|---|
| $\mathcal{R}$ | $(F^\text{ob}_\text{BB}/\sigma T_\text{ob}^4)^{1/2}$, their eq. (1) — dimensionless |
| $r_0$ | size at the base of the flow |
| $\Gamma = \eta$ | bulk Lorentz factor in the coasting phase; $\eta$ is also the dimensionless entropy (baryon loading) |
| $r_\text{ph}$ | photospheric radius |
| $r_s = \eta\,r_0$ | saturation radius |

The equations, in their numbering:

$$\mathcal{R} = \left(\frac{F^\text{ob}_\text{BB}}{\sigma T_\text{ob}^4}\right)^{1/2} = 1.06\,\frac{(1+z)^2 r_\text{ph}}{d_L\,\Gamma} \tag{1}$$

$$\eta = \left[1.06\,(1+z)^2 d_L \frac{Y F^\text{ob}\sigma_T}{2 m_p c^3 \mathcal{R}}\right]^{1/4} \tag{4}$$

$$r_0\,(r_\text{ph} > r_s) = \frac{4^{3/2}}{1.48^6\,1.06^4}\,\frac{d_L}{(1+z)^2}\left(\frac{F^\text{ob}_\text{BB}}{Y F^\text{ob}}\right)^{3/2}\mathcal{R} \tag{5}$$

Note that eq. (5) contains $F^\text{ob}_\text{BB}/(Y F^\text{ob}) = f_\text{BB}/Y$ — **the Phase 1 quantity feeds directly into it**, which is why $f_\text{BB}$ was defined in Pe'er's convention there.

---

## 2. Decisions and who made them

### 2.1 Redshift handling — *user, at the planning stage*

Only GRB080916C has a spectroscopic redshift ($z = 4.35$). Rather than assume a value for the other three, $z$ is **swept** over a log grid and the sensitivity reported.

Grid: $z \in [0.5, 5]$, 25 points, fiducial $z = 2$ (close to the long-GRB median). Chosen by the user from three options; the wider $[0.1, 8]$ and tighter $[1, 4]$ alternatives were rejected as, respectively, dominated by physically unlikely tails and an under-exploration of the uncertainty.

The CSV carries one row per interval **per redshift** (263 rows) so the full curve is reproducible. The LaTeX table collapses to one row per interval at the measured or fiducial $z$, flagged with a dagger where assumed.

### 2.2 $Y$ fixed at 1, with scalings stated — *Claude, following Pe'er et al.*

$Y \equiv \epsilon/\epsilon_\gamma \ge 1$ is the ratio of total fireball energy to energy radiated in $\gamma$-rays. It is **not measurable from the prompt spectrum alone**, so it cannot be derived here.

Pe'er et al. handle this by quoting results at $Y = 1$ with the scaling explicit, e.g. $\Gamma = (305 \pm 28)Y_0^{1/4}$. The same is done here:

$$\Gamma \propto Y^{1/4},\qquad r_0 \propto Y^{-3/2},\qquad r_\text{ph} \propto Y^{1/4},\qquad r_s \propto Y^{-5/4}$$

These four scalings were derived from eqs. (1), (4) and (5) and then **checked against the exponents Pe'er et al. print** in their GRB 970828 and GRB 990510 results — all four match. That is an independent confirmation that the equations were transcribed correctly.

*Rejected alternative:* adopting a literature value such as $Y \sim 2$–$5$ from afterglow efficiency arguments. That would bury an assumption inside a quoted number; the scaling form lets a reader apply their own $Y$.

### 2.3 Which branch of $r_0$ — *determined from the data, not assumed*

Pe'er et al. give two cases: eq. (2) for $r_\text{ph} < r_s$ (saturated) and eq. (5) for $r_\text{ph} > r_s$. Only in the second case can $r_0$ be determined; in the first, $\Gamma(r_\text{ph})$ is not recoverable.

The code computes both $r_\text{ph}$ and $r_s$ per Monte Carlo draw and records `saturated_fraction`, the fraction of draws falling in the saturated regime. **It is 0% for every interval in the sample**, so eq. (5) is the correct branch throughout and no interval needed to be dropped. This is checked at runtime rather than assumed.

### 2.4 Correlated errors from shared draws — *Claude*

$\mathcal{R}$, $f_\text{BB}$, $F^\text{ob}$ and $kT$ all derive from the **same** parameter draws, so their uncertainties stay correlated when combined into $\Gamma$ and $r_0$. Recomputing Phase 1 here rather than reading its CSV is what makes this possible — combining two independent MC runs would have treated correlated quantities as independent and mis-stated the error bars.

Same sample count (10 000), energy band (observer-frame 1 keV – 10 MeV) and grid as Phase 1, so `f_bb` in this CSV is the same quantity as in `bb_flux_fraction.csv` up to MC noise. **Updated in the RNG/seeding overhaul (Phase B, 2026-09-02):** the two files no longer share a literal seed (`12345`) — each derives its own via `seed_from_name(__file__)`, so a rerun of this file reproduces itself exactly but is no longer bit-for-bit against `bb_flux_fraction.csv`'s draws; only the underlying distribution is unchanged.

### 2.5 Shared decomposition code — *Claude*

`component_energy_fluxes` and `draw_model_samples` were promoted into `grb_research` (`grb_calculations.py`) so Phase 1 and Phase 2 use one implementation of the component split rather than two copies that could drift. Phase 1 was refactored onto it and verified unchanged.

### 2.6 Third panel for $r_\text{ph}(z)$ — *user request, 2026-09-04*

`make_plot()` originally plotted only $r_0(z)$ and $\Gamma(z)$ (two panels); $r_\text{ph}$ was in the CSV and the `PhotosphereResult` dataclass (`r_ph` field) but had no figure of its own. The gap surfaced while addressing a referee-style weakness (`quick-fixes-high-priority.md` item 4): the paper claims GRB131014A/GRB231129C's $r_\text{ph}$ offset from GRB080916C is "largely a redshift artifact," backed only by a table and prose numbers. A quick, unstyled one-off plot (`_quicklook_rph.png`, not committed, deleted after this change) made the effect immediately legible in a way the numbers alone hadn't — the user's own words were "there's no way I'd have understood the r_ph correction without this plot." That plot is now promoted into the real, styled figure: `make_plot()` takes a third `(axes[2], "r_ph", ...)` panel identical in structure to the other two, `figsize` widened from `(12.5, 5.0)` to `(18.0, 5.0)`, and the shared-legend `tight_layout` rect fraction scaled by `12.5/18.0` to keep the legend column's absolute width constant across the figure-width change. No new computation — `r_ph` was already being computed and written to the CSV every run, only the plot was missing it.

Paper side: `fig:photospheric`'s caption updated from "Base radius $r_0$ (left) and bulk Lorentz factor $\Gamma$ (right)" to name all three panels; no other section needed a change since the surrounding prose already discussed $r_\text{ph}$ at length (`subsec:photospheric`, including the item-4 rewrite) without previously being able to point at a figure for it.

---

## 3. Validation

### 3.1 Reproducing Pe'er et al.'s own worked example

`validate_against_peer2007()` runs at every invocation. Using only their published inputs for GRB 970828 ($z = 0.9578$, $T_\text{ob} = 78.5$ keV, $\mathcal{R} = 1.88\times10^{-19}$, $f_\text{BB} = 0.64$, $d_L = 1.94\times10^{28}$ cm):

| quantity | computed | published |
|---|---|---|
| $\Gamma$ | 286 | $305 \pm 28$ |
| $r_0$ | $2.94\times10^{8}$ cm | $(2.9 \pm 1.8)\times10^{8}$ cm |
| $r_\text{ph}$ | $2.57\times10^{11}$ cm | $2.7\times10^{11}$ cm |
| $r_s$ | $8.40\times10^{10}$ cm | $9.0\times10^{10}$ cm |

$r_0$ matches to three significant figures; $r_\text{ph}$ and $r_s$ to ~5%; $\Gamma$ sits inside their quoted $1\sigma$. The residuals are consistent with their inputs being printed to 3 significant figures. **This is the primary evidence that the implementation is correct**, and it runs as part of `main()` so it cannot silently rot.

### 3.2 Comparison with the $\gamma\gamma$-opacity bound — and what it constrains

**These are not two estimates of the same quantity.** $\Gamma_\text{min}$ is a one-sided bound. We observed a 27.4 GeV photon from GRB080916C; had $\Gamma$ been below $\Gamma_\text{min}$, that photon would have pair-produced on the burst's own softer photons and never reached us. The observation therefore establishes $\Gamma \geq \Gamma_\text{min}$ and nothing more — the comparison can only *fail*, by returning a thermal $\Gamma$ below the floor. It cannot "agree".

Our own bounds have since been recomputed **per episode**, on exactly the interval boundaries used here (see `lorentz_factor.md`). That makes the comparison internal and self-consistent:

| episode | $\Gamma_\text{min}$ (opacity) | $\Gamma$ (thermal, $Y{=}1$) | clears the floor? |
|---|---|---|---|
| T90 | $134 \pm 1$ | $752^{+13}_{-12}$ | yes |
| EX0 | $78^{+5}_{-4}$ | $852^{+27}_{-25}$ | yes |
| TR1 | $86^{+5}_{-5}$ | $861^{+34}_{-29}$ | yes |

The thermal values clear the opacity floor everywhere, by factors of five to ten. **But the floors are far too weak to constrain $Y$**: from EX0, $\Gamma \geq \Gamma_\text{min}$ requires only $Y \geq (78/852)^4 \approx 7\times10^{-5}$, and $Y \geq 1$ by definition. The check therefore confirms consistency and yields no information about $Y$.

**An earlier version of this note claimed $Y \gtrsim 1.9$. That claim is withdrawn.** It came from comparing against \citet{Abdo2009FermiObservations080916C}'s much stricter published bound ($887 \pm 21$), which is computed on *different time bins* — their bin b spans 3.6–7.7 s, overlapping our BB-bearing episodes by only ~1.3 s, and the interval that does cover most of bin b (our TR2, 4.864–15.040 s) has no blackbody component at all. Since $\Gamma$ demonstrably varies between episodes, that comparison mixed different quantities. Recomputing our own bound on matching intervals removed the apparent tension entirely.

Note also that our bounds are systematically weaker than published treatments — our TR2 gives 129 where Abdo et al. obtain 887 for an overlapping interval containing the same ~2 GeV photon. The gap is too large to be explained by the variability timescale and is methodological: Lithwick & Sari's Limit A is a simplified analytic form. Our $\Gamma_\text{min}$ should be read as a conservative limit, not as competitive with the literature.

**The fix, and it is worth doing:** `codes-for-paper/lorentz_factor/lorentz_factor.py` currently computes $\Gamma_\text{min}$ only for T90. Extending it to every episode would put the opacity bound and the thermal estimate on *identical* interval boundaries, making the comparison internal, self-consistent, and free of Abdo et al.'s binning. That would turn this from an indicative remark into a defensible constraint on $Y$. Tracked in `BUGS.md`.

---

## 4. Results

At the measured or fiducial redshift ($Y = 1$):

| GRB | $\Gamma$ | $r_0$ [cm] | $r_\text{ph}$ [cm] |
|---|---|---|---|
| GRB080916C ($z=4.35$) | 752 – 861 | $0.9$–$1.4\times10^{7}$ | $6.1\times10^{11}$ – $1.4\times10^{12}$ |
| GRB131014A ($z=2^\dagger$) | 458 – 714 | $0.7$–$4.7\times10^{8}$ | $3.4$–$7.6\times10^{12}$ |
| GRB140206B ($z=2^\dagger$) | 597 – 600 | $0.8$–$1.3\times10^{7}$ | $5.6$–$7.5\times10^{11}$ |
| GRB231129C ($z=2^\dagger$) | 385 – 426 | $6.2$–$8.9\times10^{8}$ | $3.5$–$4.8\times10^{12}$ |

Redshift dependence is mild and analytically understandable: $\Gamma$ rises monotonically with $z$, while $r_0 \propto d_L/(1+z)^2$ is non-monotonic and peaks near $z \sim 1.7$, because $d_L/(1+z)^2$ is the angular-diameter-like distance.

**Bearing on the existing text (`BUGS.md` OBS-04).** §6 and §7 both claim photospheric radii "of order $10^{11}$–$10^{12}$ cm", cited from Ryde 2010 / Guiriec 2011. Our own values only partly support this: GRB080916C ($6\times10^{11}$–$1.4\times10^{12}$) and GRB140206B ($\sim7\times10^{11}$) sit in that range, but GRB131014A and GRB231129C come out at $3$–$8\times10^{12}$ at $z = 2$, i.e. an order of magnitude higher. The sentences must be rewritten against these numbers in Phase 3, not left as a borrowed range.

---

## 4.5 `kt_bb_keV` here is an MC median, not the raw fit value

This CSV's `kt_bb_keV`/error columns are the median (and 16th/84th-percentile spread) of the resampled `kt_bb` draws used to propagate uncertainty into $\Gamma$, $r_0$, $r_\text{ph}$ — consistent with §2.4's point that everything here shares one draw set with Phase 1. `bb_flux_fraction.csv` has a same-named `kt_bb_keV` column that is instead the raw point-estimate fit value, so the two disagree at the ~0.01–0.02% level despite the shared name. Found and confirmed during the Phase 4 spot-check, `BUGS.md` OBS-09 — not a bug, but worth knowing before diffing the two files.

**Surfaced to the reader, 2026-09-01:** `grb_paper_weaknesses_and_fixes.md`'s Priority 3 editorial pass flagged the resulting kT mismatch between `tab:photospheric` (this table) and `tab:bbfraction` (44.46 vs 44.45 keV for GRB080916C T90) as a possible inconsistency. `csv_to_latex.py`'s caption now cross-references `tab:bbfraction` explicitly, stating the two are different statistics of the same fit and are expected to agree only within their stated uncertainties — see `review-resolution.md` Priority 3 item 5.

## 5. Limitations and open questions

- **$Y$ is unconstrained.** Everything scales with it as above. If $Y \sim 3$, $\Gamma$ rises by only $\sim 32\%$ but $r_0$ falls by a factor $\sim 5$ — $r_0$ is by far the more $Y$-sensitive quantity.
- **Planck, not Wien.** Pe'er et al. §4 note that Comptonisation-dominated emission produces a Wien spectrum, giving $T \to 3k_BT/2.7k_BT$, a $\sim 10\%$ temperature shift and $\sim 5\%$ systematic in $\eta$. The fits assume a Planck function throughout.
- **The high-latitude interpretation** underlying the $1.06$ and $1.48$ prefactors assumes emission is dominated by photons from $\theta \simeq 0$ at $t^\text{ob} < t_\text{break}$. Our intervals are not resolved finely enough to locate a temporal break, so this is assumed rather than verified.
- **Errors are statistical only.** They come from the fit covariance and do not include the systematic from $Y$, the Planck/Wien choice, or the assumed redshift for three of four bursts.
- **Magnetisation $\sigma_0$ not yet computed.** The remaining `ToDo.md` stretch item needs the hybrid-outflow framework of Gao & Zhang (2015), now in the literature folder; Pe'er's pure-fireball equations cannot give it.

---

---

## Table formatting: normalisation and precision

*Decided by the user, 2026-08-21, after reading the rendered tables.*

Columns are divided by a common power of ten carried in the header rather than repeating `\times 10^{n}` on every row:

| table | column | norm |
|---|---|---|
| Table 4 | $F_\text{BB}$, $F_\text{total}$ | $10^{-6}$ erg cm$^{-2}$ s$^{-1}$ |
| Table 4, Table 5 | $f_\text{BB}$ | $10^{-3}$ |
| Table 5 | $r_0$ | $10^{7}$ cm |
| Table 5 | $r_\text{ph}$ | $10^{11}$ cm |

$kT$ and $\Gamma$ are left unnormalised — both are already of order $10$–$10^{3}$, and normalising them made the numbers *less* readable ($\Gamma = 3.854 \times 10^2$ rather than $385$).

**Decimal places are derived numerically inside `csv_to_latex.py`, not chosen by eye.** The rule, implemented as `decimals_for()`:

$$\text{decimals} = \left\lceil \text{sig} - 1 - \lfloor \log_{10}(\min(\text{error})/\text{norm}) \rfloor \right\rceil, \quad \text{sig} = 2$$

i.e. enough decimals that the *smallest* error in the column still shows two significant figures. It is computed from the dataframe at generation time, so the formatting adapts if the data changes instead of going stale.

This was worth doing numerically: eyeballing had produced 3 decimals where 1 was right ($f_\text{BB}$, whose errors are $\sim5$ in units of $10^{-3}$) and 4 decimals on the Amati table where the errors are of order 7. The same helper now runs in all three `csv_to_latex.py` scripts.

## 6. Files

| file | role |
|---|---|
| `pe_er_photosphere.py` | computation, CSV and figure, plus the Pe'er validation |
| `pe_er_photosphere.csv` | one row per interval per redshift (263 rows) |
| `pe_er_photosphere.png` / `.pdf` | $r_0(z)$ and $\Gamma(z)$ |
| `csv_to_latex.py` → `photospheric_table.tex` | paper table at measured/fiducial $z$ |
| `photospheric_radius.md` | this file |
