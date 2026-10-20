# Draft Plan for GRB Time-Resolved Analysis

### Scatter Plots

- [ ] Produce scatter plots of key spectral parameters across the time-resolved intervals of each GRB to visualize overall temporal trends. 
- [ ] Generate comparison scatter plots showing how the inclusion of additional spectral components (e.g., a power-law [PL] or blackbody [BB]) affects the base model parameters. In particular, compare $E_\text{peak}$, $\alpha$, and $\beta$ before and after adding the extra PL/BB component.
- [ ] Create scatter plots of energy fluence versus the same parameters to investigate possible correlations between spectral evolution and total emitted energy.
- [ ] Because this work focuses on the additional BB component, examine the relationship between $\text{BB}_\text{AMP}$ and $\text{BB}_\text{kT}$ over the time-resolved intervals to look for systematic trends. Similar checks can be performed for other parameter pairs if warranted by the data.

### Additional Correlation Studies

- [ ] **Hardness–intensity/flux correlations:** Track $E_\text{peak}$ or $\text{kT}$ versus flux to search for classic “hard-to-soft” or “tracking” behaviors.
- [ ] **Isotropic energy vs. BB temperature:** Test theoretical scaling relations predicted by photospheric models.
- [x] **Redshift dependence:** For GRBs with unknown redshift, explore how inferred physical parameters (e.g., $E_\text{iso}$, photospheric radius) vary across a plausible redshift range.
  <br>Done (2026-10-01/02). $E_\text{iso}$ at $z = 1, 2, 3, 5, 7$ in `codes-for-paper/amati_relationship/amati_relationship_unknown_z.csv`; photospheric radius / $\Gamma$ / $r_0$ swept over $z \in [1, 7]$ in `codes-for-paper/photospheric_radius/` (`pe_er_photosphere.py`, `Z_MIN, Z_MAX = 1.0, 7.0`); $f_\text{BB}$ vs $z$ in `codes-for-paper/bb_fraction/bb_flux_fraction_rest_vs_z.*`; the $\gamma\gamma$ $\Gamma_\text{min}$ limits swept over $z$ in `codes-for-paper/lorentz_factor/*_unknown_z*` (`lorentz_factor.md` §13).

### Butterfly (Spectral-Model) Plots

- [x] Produce model spectra with full error propagation (**butterfly plots**) for both **SAFE** and **UNSAFE** fits to illustrate the statistical range of each model.
- [x] For joint time-integrated (TI) and time-resolved (TR) analyses, overplot butterflies from different model quality categories (e.g., SAFE, GOOD, BEST) on the same figure to enable direct visual comparison.

### Blackbody Component Diagnostics

- [x] **Photospheric radius & Lorentz factor:** Use BB temperature and flux to estimate the emitting radius and bulk Lorentz factor (e.g., Pe'er 2007 method).
  <br>Done in `codes-for-paper/photospheric_radius/`. Implementation validated by reproducing Pe'er et al.'s own GRB 970828 worked example ($r_0$ to 3 significant figures). See `photospheric_radius.md`.
- [x] **Multi-component decomposition:** Quantify the fractional BB contribution to total energy flux and track its evolution throughout the burst.
  <br>Done in `codes-for-paper/bb_fraction/`. $f_\text{BB}$ computed for all 13 BB-inclusive episodes, in Pe'er's bolometric convention so it feeds the photospheric calculation directly. See `bb_fraction.md`.
- [ ] **Fireball parameter estimation:** From $kT$ and flux, estimate additional fireball parameters such as baryon loading, initial radius, and possible jet magnetization.
  <br>*Partially done.* Baryon loading $\eta$ and initial radius $r_0$ both fall out of the Pe'er method ($\eta = \Gamma$ in the coasting phase) and are reported. **Magnetization $\sigma_0$ is still outstanding** — Pe'er's pure-fireball equations cannot give it; it needs the hybrid-outflow framework of Gao \& Zhang (2015), now in the literature folder. *Status 2026-10-20:* scoped as Phase 6 (`PLAN.md`); `codes-for-paper/outflow_composition/PHASE6_PLAN.md` exists but is a plan only — no code, CSV or table has been written, and the folder is untracked. Left unticked.
- [x] **Lorentz factor lower limits:** Use opacity arguments (e.g., $\gamma\gamma$ pair production) to constrain the minimum bulk Lorentz factor.
  <br>Done in `codes-for-paper/lorentz_factor/`, now evaluated per episode rather than for T90 alone, with MC errors. Note the limits are conservative and weaker than published treatments of the same burst — see `lorentz_factor.md` §5 for the uncertainty budget.

### Model Robustness and Alternatives

- [ ] **Alternative models:** Fit synchrotron, Comptonized, or two-BB models to test whether the BB signature is genuine or an artifact of model choice.
  <br>*Not done (checked 2026-10-20).* The fitted-model database (`results.json`) contains only PL, CPL, BAND and SBPL, each with and without a BB; no synchrotron, physical-Comptonization or two-BB fits exist. The related robustness checks that *have* been done — the GBM-only vs joint refit (`codes-for-paper/gbm_only_refit/`) and the LAT high-energy-index consistency check (`codes-for-paper/beta_lat_consistency/`) — test the data choice, not the model choice, so this item stays open.
- [x] Compare fit statistics (likelihood ratio, AIC/BIC) across all tested models to establish the significance of the BB component.

### GRB Characterization

- [x] Test each burst against well-known empirical relations such as the **Amati** and ~~**Yonetoku**~~ relationships.
- [x] Derive physical parameters such as isotropic energy ($E_\text{iso}$), peak luminosity, and, where feasible, constraints on jet overflow or jet-break signatures.

### Population Context and Comparative Analysis

- [ ] **Comparative plots:** Place the four BB-dominated GRBs within the context of a larger GRB population (e.g., the Fermi-GBM catalog) to determine whether they occupy distinct regions in $E_\text{peak}$–$E_\text{iso}$ space or in $kT$ distributions.
- [ ] Examine whether BB-dominated bursts differ systematically in duration, fluence, or spectral hardness compared with the general GRB sample.

### Literature Context

- [ ] Survey recent GRB spectral and temporal studies to identify comparable analyses and benchmark results.
- [ ] Highlight key findings from the literature that can guide the physical interpretation of trends observed in the current dataset.

---

### Status note (2026-10-20)

Unticked items that remain: the scatter/correlation studies (partial plots exist in `codes-for-paper/model_parameters/` and `epeak_vs_kt.py`, but no item above has been confirmed as fully covered, so none is ticked), $\sigma_0$ (Phase 6, plan only), alternative spectral models, the population comparison, and the literature survey. For priorities, ordering and decisions pending from the user, see `md_file_refs/PLAN.md` and `md_file_refs/HANDOFF.md`, not this file. Work not listed above but done since: Norris-profile $t_v$ fits (`codes-for-paper/variability_analysis/`) feeding the $\gamma\gamma$ opacity limit in `codes-for-paper/lorentz_factor/` (`lorentz_factor.md` §11, §15).
