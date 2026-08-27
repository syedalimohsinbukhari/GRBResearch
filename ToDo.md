# Draft Plan for GRB Time-Resolved Analysis

### Scatter Plots

- [ ] Produce scatter plots of key spectral parameters across the time-resolved intervals of each GRB to visualize overall temporal trends. 
- [ ] Generate comparison scatter plots showing how the inclusion of additional spectral components (e.g., a power-law [PL] or blackbody [BB]) affects the base model parameters. In particular, compare $E_\text{peak}$, $\alpha$, and $\beta$ before and after adding the extra PL/BB component.
- [ ] Create scatter plots of energy fluence versus the same parameters to investigate possible correlations between spectral evolution and total emitted energy.
- [ ] Because this work focuses on the additional BB component, examine the relationship between $\text{BB}_\text{AMP}$ and $\text{BB}_\text{kT}$ over the time-resolved intervals to look for systematic trends. Similar checks can be performed for other parameter pairs if warranted by the data.

### Additional Correlation Studies

- [ ] **Hardness–intensity/flux correlations:** Track $E_\text{peak}$ or $\text{kT}$ versus flux to search for classic “hard-to-soft” or “tracking” behaviors.
- [ ] **Isotropic energy vs. BB temperature:** Test theoretical scaling relations predicted by photospheric models.
- [ ] **Redshift dependence:** For GRBs with unknown redshift, explore how inferred physical parameters (e.g., $E_\text{iso}$, photospheric radius) vary across a plausible redshift range.

### Butterfly (Spectral-Model) Plots

- [x] Produce model spectra with full error propagation (**butterfly plots**) for both **SAFE** and **UNSAFE** fits to illustrate the statistical range of each model.
- [x] For joint time-integrated (TI) and time-resolved (TR) analyses, overplot butterflies from different model quality categories (e.g., SAFE, GOOD, BEST) on the same figure to enable direct visual comparison.

### Blackbody Component Diagnostics

- [x] **Photospheric radius & Lorentz factor:** Use BB temperature and flux to estimate the emitting radius and bulk Lorentz factor (e.g., Pe'er 2007 method).
  <br>Done in `codes-for-paper/photospheric_radius/`. Implementation validated by reproducing Pe'er et al.'s own GRB 970828 worked example ($r_0$ to 3 significant figures). See `photospheric_radius.md`.
- [x] **Multi-component decomposition:** Quantify the fractional BB contribution to total energy flux and track its evolution throughout the burst.
  <br>Done in `codes-for-paper/bb_fraction/`. $f_\text{BB}$ computed for all 13 BB-inclusive episodes, in Pe'er's bolometric convention so it feeds the photospheric calculation directly. See `bb_fraction.md`.
- [ ] **Fireball parameter estimation:** From $kT$ and flux, estimate additional fireball parameters such as baryon loading, initial radius, and possible jet magnetization.
  <br>*Partially done.* Baryon loading $\eta$ and initial radius $r_0$ both fall out of the Pe'er method ($\eta = \Gamma$ in the coasting phase) and are reported. **Magnetization $\sigma_0$ is still outstanding** — Pe'er's pure-fireball equations cannot give it; it needs the hybrid-outflow framework of Gao \& Zhang (2015), now in the literature folder.
- [x] **Lorentz factor lower limits:** Use opacity arguments (e.g., $\gamma\gamma$ pair production) to constrain the minimum bulk Lorentz factor.
  <br>Done in `codes-for-paper/lorentz_factor/`, now evaluated per episode rather than for T90 alone, with MC errors. Note the limits are conservative and weaker than published treatments of the same burst — see `lorentz_factor.md` §5 for the uncertainty budget.

### Model Robustness and Alternatives

- [ ] **Alternative models:** Fit synchrotron, Comptonized, or two-BB models to test whether the BB signature is genuine or an artifact of model choice.
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

