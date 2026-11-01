# Outflow composition (η, σ₀) — method note

Phase 6 of `md_file_refs/PLAN.md`. This is the per-folder method note required by `CLAUDE.md`; it records *how and why*, as the
work happens. `PHASE6_PLAN.md` (2026-09-27) is the plan this started from and `CI.md` its review; where this note and the
plan disagree, **this note wins** (the disagreements are listed in §2).

**Status (2026-10-24): stage 1 and stage 2 built and validated.** Stage 1 (`gao_zhang_inversion.py`): the Gao & Zhang inversion for the 13
BEST-BB episodes (main table, `gao_zhang_inversion.csv`) and, reported separately, for the SAFE/MARGINAL `X_BB` partner of 7 non-BB episodes
(`gao_zhang_inversion_partner.csv`). Stage 2 (`zhang_peer_suppression.py`, `zhang_peer_suppression.csv`): the Zhang & Pe'er suppression bound for the 13
episodes whose BEST model has no blackbody. Not built: the LaTeX table, the figure, the paper text. Nothing in the paper references this folder yet.
Three band decisions are still open (§8); the defaults in use are stated there. Two validation figures (§6.6) are drawn by `zhang_peer_suppression_plot.py`.

## 1. What stage 1 computes

For each BB-detected BEST model (13 episodes: GRB080916C T90/EX0/TR1; GRB131014A T90/EX0/TR1/EX1/TR2; GRB140206B EX0/TR1;
GRB231129C T90/EX0/TR1), the observed blackbody temperature $kT_\text{ob}$, blackbody flux $F_\text{BB}$ and total flux
$F_\text{ob}$ (integrated over observer-frame 10 keV – 10 GeV since 2026-10-25, a wider band than `pe_er_photosphere.py`'s 1 keV – 10 MeV; $f_\text{th}=F_\text{BB}/F_\text{ob}$)
are inverted to the central-engine pair $(\eta, 1+\sigma_0)$, the photospheric radius $r_\text{ph}$ and Lorentz factor
$\Gamma_\text{ph}$, for an *assumed* jet-base radius $r_0$ and $f_\gamma=L_\gamma/L_w$.

Source: Gao & Zhang (2015), ApJ 801, 103 (arXiv:1409.3584), §4. Table 2 selects the regime from the observables (II, III, V or VI); eqs.
(38), (39), (40) give the inversion in regimes II, III/VI and V. Each printed row is one `PowerLaw` entry in the script, in the
paper's own form $c\,(1+z)^a (kT/kT_0)^b (F_\text{BB}/F_0)^c r_{0,9}^d (f_\text{th}/0.1)^e f_\gamma^f d_{L,28}^g$, so the transcription
can be audited line by line against the PDF. Regimes I and IV ($r_\text{ph}<r_\text{ra}$) are degenerate and cannot be inverted
(the paper's own statement); a trial where no Table 2 regime holds is reported as `none`, not defaulted.

Regime of a row: that of the *median* inputs. The summaries (median, 16/84 percentiles) use only the MC draws that fall in the same
regime, and `regime_fraction` records how many do, so a row straddling a boundary is visible rather than averaged across two formulae.
A draw satisfying two regimes' criteria at once (possible within ~1% of a boundary, because the paper's coefficients are rounded to
three figures) is `ambiguous` and excluded from the summary. `physical` is `1+σ₀ ≥ 0.97`, the tolerance coming from the closure test (§4).

## 2. Corrections to `PHASE6_PLAN.md`, found 2026-10-24 by reading the PDFs

Verified against the papers themselves (text extracted and the table/equation pages rendered and read), not from the earlier summary:

1. **Regime selection uses Table 2, not Table 1.** The plan said to classify with Table 1 ("$\eta^{12/5}(1+\sigma_0)^{3/5}$ vs
   $7.22\times10^5L_{w,52}^{3/5}r_{0,9}^{-3/5}$"). Table 1 is written in $\eta$, $\sigma_0$ and $L_w$, which are the unknowns;
   Table 2 is the observational version and is what the inversion needs. Table 1 is used only to *generate* test points (§4).
2. **The top-down formulae exist only for δ = 1/3.** The paper: "For the top-down approach, we only derived the formulae for the δ = 1/3
   case" (§6). The plan's Verification #5 ("check the regime classification is stable under δ = 1/3 ± ...") cannot be done with the
   paper's equations; Table 3 gives general-δ criteria in $(\eta,\sigma_0,L_w)$, not in observables. Resolved in §8 (option a). The paper's
   own statement (§6) is that $T_\text{ob}$ and $F_\text{BB}$ at δ = 1/3 are upper limits for δ < 1/3 over the range relevant to GRBs.
3. **Li et al. (2024)'s Table 1 is not the spot-check the plan hoped for.** Its GRB080916C row is a featureless-Band slice with an
   *assumed* kT of 5 keV (a lower limit on $1+\sigma_0$ from a simulated thermal component), not a measured blackbody. The GRB110721A
   rows list $kT$ and $f_\text{th}$ but not $F_\text{ob}$ or $z$, so $1+\sigma_0$ cannot be reproduced. Only the $r_0$ scaling can be
   checked: $22.206/4.784 = 4.6417$ against $10^{2/3}=4.6416$ ✓ (`check_li_r0_scaling()`). Also, Li et al.'s prefactors are the same group's
   restatement of Gao & Zhang's, so they are a weak independent check; the strong ones are in §4.
4. **The "quiet assumption" of the plan does not arise.** The plan worried that the project's pure-fireball $r_0$ (derived at σ₀ = 0)
   was being reused as an input to a calculation about σ₀ > 0. §4 shows the inversion at σ₀ = 0 *is* Pe'er et al. (2007), so the
   pure-fireball $r_0$ is not an independent input to be justified: it is the value of $r_0$ at which the inversion returns σ₀ = 0.
   In regime III $1+\sigma_0 \propto r_0^{2/3}/f_\gamma$ (all regimes: $r_0^{2/3}$), so *any assumed $r_0$ above that threshold requires σ₀ > 0*.
   That is a statement about the assumption, and is how the result should be worded in the paper: composition conditional on $r_0$.

## 3. Judgement calls

| Choice | Value | Decided by | Rejected / note |
|---|---|---|---|
| Scope | all four bursts; assumed-$z$ grid $z=1,2,3,5,7$ for the three without a redshift, $z=4.35$ for GRB080916C | **user, 2026-10-24** | GRB080916C only (the plan). Extends the "composition" result to bursts whose $d_L$ is assumed: every number for them is conditional on $z$ (the paper's section-5 convention) |
| $r_0$ trials | two sets: `GZ2015` = 1e8, 1e9, 1e10 cm; `project_anchor` = threshold × 0.1, 1, 10 | **user, 2026-10-24** ("both sets") | the other two options were one set each |
| Anchor definition | Pe'er+07 $r_0$ (Y = 1, median over the burst's BB episodes at the same $z$) × $f_\gamma^{3/2}$ | Claude (refinement of the plan's anchor) | Pe'er's $r_0\propto Y^{-3/2}$ and $f_\gamma=0.5$ is $Y=2$; without the rescaling the anchor sits at $1+\sigma_0\approx1/f_\gamma=2$ by construction (seen in the closure run) rather than at 1. A per-burst constant, because the method assumes a constant $r_0$ through a burst |
| $f_\gamma$ | 0.5 | paper's and Li+24's fiducial; **user confirmed, 2026-10-24** | $1+\sigma_0\propto1/f_\gamma$ in every regime, so any other value is an exact rescale |
| $\delta$ | 1/3 | the paper's only available value; **user decided 2026-10-24: option (a)**, cite the paper's §6 statement and drop the numerical δ check | option (b), deriving general-δ top-down formulae, is recorded in §8 as a possible later experimental setup, not part of Phase 6 |
| Flux band | observer-frame 10 keV – 10 GeV for $F_\text{BB}$ and $F_\text{ob}$ (stage 1 only; Phase 1/2 stay at 1 keV – 10 MeV) | **user, 2026-10-25** ("stage 1 only, anchor recomputed") | 1 keV – 10 MeV inherited from `pe_er_photosphere.py` (so the inputs were the Pe'er+07 ones); widening Phase 1/2 too was rejected (it changes published $f_\text{BB}$, $r_\text{ph}$, $\Gamma_\text{ph}$). Consequence: the `project_anchor` $r_0$ is no longer read from `pe_er_photosphere.csv` but recomputed by Pe'er+07 eq. (5) on stage 1's own draws in the new band (`compute_pure_fireball_anchors`) |
| Non-BB episodes | Zhang & Pe'er root-find (stage 2); duration-fallback episodes kept, labelled `t_v_source = duration` | **user, 2026-10-24** | forced-BB refit + Li et al. closed form; BB episodes only |
| Which models feed which table | see §3a | **user, 2026-10-24** | the partner could have been any BB-bearing model of the episode (e.g. BAND_BB when the BEST is SBPL); rejected as a cross-continuum comparison |
| Stage 2 branch / $R_0$ | $\eta\ge\eta_*$ branch; $R_0=c\,t_v/(1+z)$ | **user, 2026-10-24** ("go with those defaults") | the $\eta<\eta_*$ curves need an extra parameter; $R_0=c\,t_v$ is excluded by `validate_r0_convention()` (§6.2) |
| Stage 2 $L_\gamma$ band | observer-frame 10 keV – 10 GeV (primary), rest-frame 1 – 10⁴ keV (sensitivity) | Claude, on the evidence of §6.2; **primary band confirmed by the user, 2026-10-25**; sensitivity row not addressed, left as is | the Amati rest-frame band is excluded by the Abdo epoch-b reproduction (2× low in $L_\gamma$); 8 keV – 40 MeV is also consistent with it (§6.2) |
| Stage 2 test band | observer-frame 10 keV – 10 MeV (primary, "fit"), 1 keV – 10 MeV (sensitivity, "extrapolated") | Claude; **lower edge 10 keV confirmed by the user, 2026-10-25** (8 keV rejected) | no basis in Zhang & Pe'er; see §6.4 |

### 3a. Model sets (user, 2026-10-24)

1. **Main table (paper): only episodes whose BEST model contains a blackbody** — 13 episodes, stage 1, `gao_zhang_inversion.csv`.
2. **Blackbody models that are not BEST are used separately and reported separately**, flagged the way MARGINAL models are marked in the appendix tables.
3. **The partner of a non-BB BEST model X is X_BB and nothing else** (SBPL → SBPL_BB, not BAND_BB): the comparison is within the BEST model's own continuum.
   **UNSAFE partners are excluded outright**; SAFE and MARGINAL are admitted (user confirmed SAFE is allowed).

Applied to the 13 non-BB episodes (status and $\Delta C$-stat from `results.json`; all $\Delta C<28.74$, which is why none is BEST):

| episode | BEST | partner | status | ΔC vs BEST | kT_BB | stage 1 on partner? |
|---|---|---|---|---|---|---|
| GRB080916C TR2 | BAND | BAND_BB | SAFE | 11.4 | 46.0 keV | yes |
| GRB080916C TR3 | SBPL | SBPL_BB | SAFE | 24.1 | 49.1 keV | yes |
| GRB080916C TR4, EX1, TR5 | SBPL | SBPL_BB | UNSAFE | 0.7, 5.4, 1.4 | 26, 36, 44 keV | no |
| GRB140206B T90 | BAND | BAND_BB | SAFE | 24.9 | 27.1 keV | yes |
| GRB140206B TR2 | BAND | BAND_BB | SAFE | 22.0 | 31.5 keV | yes |
| GRB140206B TR3 | BAND | BAND_BB | MARGINAL | 6.7 | 16.9 keV | yes |
| GRB140206B TR4 | BAND | BAND_BB | SAFE | 23.6 | 5.0 keV | yes |
| GRB140206B TR5 | PL | PL_BB | UNSAFE | 28.7 | 44.0 keV | no |
| GRB140206B TR6 | CPL | CPL_BB | not in `results.json` | — | — | no |
| GRB231129C EX1 | SBPL | SBPL_BB | SAFE | 21.2 | 53.1 keV | yes |
| GRB231129C TR2 | SBPL | SBPL_BB | UNSAFE | 5.2 | 8.6 keV | no |

Stage 2 runs on all 13 (it uses the BEST spectrum, so no blackbody model is involved); stage 1 runs on the 7 usable partners only. The partner rows carry
`best_model`, `bb_status`, `delta_cstat_vs_best`, `delta_dof_vs_best`, `bb_threshold_delta_cstat` and `bb_required` (False for all 162 rows).

## 4. Validation (run by `main()` before any output; all numbers from the 2026-10-24 run)

1. **Round trip, bottom-up → top-down.** The paper prints two independent sets of equations: the bottom-up eqs. (24)–(30) giving
   $kT_\text{ob}$, $F_\text{BB}$, $r_\text{ph}$, $\Gamma_\text{ph}$ from $(L_w,r_0,\eta,\sigma_0)$ with Table 1's regime criteria, and the top-down eqs. (38)–(40) with
   Table 2. 40 000 random central-engine parameter sets were pushed through the first and recovered with the second:

   | regime | n | max \|Δ(1+σ₀)\| | max \|Δη\| | max \|Δr_ph\| | max \|ΔΓ\| |
   |---|---|---|---|---|---|
   | II | 13 801 | 1.0e-3 | 2.3e-4 | 3.1e-4 | 1.1e-3 |
   | III | 8 956 | 1.8e-4 | 2.6e-3 | 1.0e-3 | 2.8e-4 |
   | V | 4 186 | 1.8e-3 | 2.2e-3 | 4.2e-3 | 2.2e-3 |
   | VI | 1 187 | 1.8e-4 | 2.6e-3 | 1.0e-3 | 2.8e-4 |

   11 853 points fell in the degenerate regimes I/IV (correctly not inverted); Table 1 and Table 2 disagreed on the regime for 17 of
   the rest, all within ~1% of a regime boundary (rounded coefficients). Agreement ≤ 0.5% is the floor set by 3-figure prefactors;
   it tests the transcription of both sets of equations and the consistency of Table 2 with Table 1.
2. **Regime-boundary continuity** (bottom-up observables of II and III on $r_\text{ph}=r_\text{c}$, V and VI on theirs): max relative difference
   6.3e-4 over 1 676 points (II/III) and 2.0e-3 over 990 points (V/VI). The plan asked for this; it tests the bottom-up transcription.
3. **Closure against Pe'er et al. (2007)** (`check_pure_fireball_closure()`): with $f_\gamma=1$ ($Y=1$) and $r_0$ set to each episode's own Pe'er $r_0$
   from `pe_er_photosphere.csv`, the inversion returns **$1+\sigma_0=0.971$–$0.980$, $\eta/\Gamma_\text{Pe'er}=1.013$–$1.031$, $r_\text{ph}$ ratio $1.021$–$1.033$** for all 13 episodes
   (regime III throughout). This is the strongest check available: two published methods and this project's own independent
   Pe'er implementation agree to 3% at σ₀ = 0, as they must. The 3% is the rounding of the paper's coefficients, and is why
   `physical` allows 1+σ₀ ≥ 0.97. Gao & Zhang's paper gives no worked numeric example to reproduce (their case study is figure-only), so this closure takes its place.
4. **$r_0$ scaling against Li+24 Table 1**: 4.6417 vs 10^(2/3) = 4.6416 ✓ (§2 item 3; the only part of that table that can be checked).

**Not validated:** an absolute reproduction of any published (η, σ₀) for a real burst. No such tuple is available in either paper (§2).

## 5. Results snapshot (CSV `gao_zhang_inversion.csv`)

The tables and every number read off the two CSVs live in `tables/outflow_composition_tables.generated.md`, written by `tables/snapshot_tables.py` (rerun it after `gao_zhang_inversion.py`; it is never edited or retyped by hand).
The sample counts, regime counts, `GZ2015` $1+\sigma_0$ per burst and $r_0$, the project-anchor range and the redshift dependence are all there, with the band, $f_\gamma$, $\delta$, cosmology, `n_samples` and seed that produced them.
The reading below is qualitative and does not quote values.

Values below 1 mean the assumed $r_0$ is *below* the pure-fireball threshold for that episode: no magnetisation is needed and the data prefer a larger $r_0$ (consistent with the Phase 2 Pe'er $r_0$ of those bursts, `pe_er_photosphere.csv`).
The inferred σ₀ rises with the assumed $r_0$ (∝ $r_0^{2/3}$) until the regime degenerates (`none`, i.e. $r_\text{ph}$ falls below $r_\text{ra}$ or no criterion holds).
The $z$ dependence of $1+\sigma_0$ at fixed $r_0$ is modest while η and $\Gamma_\text{ph}$ rise strongly with $z$.

At the project anchor × 1 the per-episode $1+\sigma_0$ is not exactly 1: the anchor is one value per burst while Pe'er's $r_0$ varies across a burst (e.g. GRB131014A), which is precisely the behaviour Gao & Zhang cite against a varying $r_0$.

### 5a. Partner models (CSV `gao_zhang_inversion_partner.csv`)

Reported separately because the blackbody is not required by the data (§3a); a result here is *conditional on a component the fit does not need*.
For stage 1 a non-significant blackbody amplitude should be read as an upper limit on $F_\text{BB}$: in all of regimes II, III and V, $1+\sigma_0\propto f_\text{th}^{-1}F_\text{BB}^{-1/3}\propto F_\text{BB}^{-4/3}$ at fixed $F_\text{ob}$, so these are, if anything, *lower* limits on $1+\sigma_0$ (the best-fit amplitude is used as given; a profile-likelihood upper limit has not been computed).
The partner table (kT, ΔC-stat against the BEST model, $1+\sigma_0$ per $r_0$) and the regime counts are in `tables/outflow_composition_tables.generated.md`; MARGINAL partners are flagged there.
One-off comparison, 2026-10-25: moving from the 1 keV – 10 MeV band to 10 keV – 10 GeV changed the regime of three `GZ2015` partner trials, all III to II, all in GRB140206B.

### 5b. Figures (`gao_zhang_inversion_plot.py`; reads both CSVs)

`gao_zhang_inversion.png/.pdf` (main table, one panel per burst) and `gao_zhang_inversion_partner.png/.pdf` (partners, separately): $1+\sigma_0$ against the assumed $r_0$, one curve per episode through its six trials
(both trial sets share the axis) with the 16–84% band, at the measured $z$ (GRB080916C) or $z=2$; the line $1+\sigma_0=1$ and the burst's pure-fireball $r_0$ (median over its BEST-BB episodes) are marked, and grey crosses mark trial $r_0$ values at which an episode has no valid regime.
Every curve has slope 2/3 in log–log (§5: $1+\sigma_0\propto r_0^{2/3}$ in every regime), so the figure shows the conditionality on $r_0$ directly: each episode crosses $1+\sigma_0=1$ near the burst's threshold and no regime exists at the largest $r_0$.
Both figures are drawn for checking; they are not yet in the paper.

## 6. Stage 2 — Zhang & Pe'er (2009) suppression bound (`zhang_peer_suppression.py`, CSV `zhang_peer_suppression.csv`)

### 6.1 What it computes

For each of the 13 episodes whose BEST model has no blackbody, the smallest $1+\sigma$ for which a pure baryonic fireball's photospheric blackbody lies at or below the *best-fit* spectrum
(Zhang & Pe'er 2009, §3): $L_\text{th}=L_w/(1+\sigma)$ and $kT_\text{ob}=kT_\text{max}(1+\sigma)^{-1/4}$ with $kT_\text{max}=(L_w/4\pi R_0^2 c\,a)^{1/4}/(1+z)$, $L_w=L_\gamma$ (100% efficiency, as the paper), $\eta\ge\eta_*$,
$R_0=c\,t_v/(1+z)$. Per MC draw (10 000, the episode's covariance; $t_v$ from `lorentz_results.csv`, the Norris value resampled split-normally where measured), $L_\gamma$, $kT_\text{max}$ and $\eta_*$ are computed
and $\sigma_\text{min}$ is found by bisection on $\ln(1+\sigma)$ against $\max_E\,[\text{blackbody}/\text{observed }\nu F_\nu]\le1$. A draw that already satisfies this at σ = 0 returns exactly 1; a draw that cannot be bracketed raises.
The CSV has $2\times2$ rows per episode and redshift: two $L_\gamma$ bands × two test bands (§6.3, §6.4). 180 rows; 13 episodes; $z$ = 4.35 for GRB080916C and 1, 2, 3, 5, 7 for the other two.
`t_v_is_duration_fallback` marks the five episodes that use their duration (no Norris $t_v$): GRB080916C TR4, EX1, TR5; GRB140206B T90, TR5. A duration overstates $R_0$ and so lowers $kT_\text{max}$; §6.7 shows that the effect on $1+\sigma_\text{min}$ has no consistent direction.

**Correction to an earlier statement of this note:** the blackbody's Rayleigh–Jeans $\nu F_\nu$ is *not* σ-independent. At fixed energy it scales as $(1+\sigma)^{-1/4}$ (the coefficient $F/T^4$ is
σ-independent but the Rayleigh–Jeans tail goes as $T$). The ratio blackbody/observed is strictly decreasing in σ at every energy — $d\ln\nu F_\nu/d\ln(1+\sigma)=-(x/4)\,e^x/(e^x-1)<0$ — so the root is unique
and a bisection is exact (`check_monotonic()` confirms it numerically).

### 6.2 Validation (`main()` runs these before any output)

1. Blackbody $\nu F_\nu$ integrates to its flux: max error 5×10⁻¹¹.
2. Monotonicity of the max ratio in ln(1+σ) over 3 temperatures × 3 spectra: largest relative step −1.8×10⁻² (must be < 0).
3. Solver: the ratio at the root is 1 (to 10⁻⁶); a spectrum far above the blackbody gives 1+σ = 1 exactly.
4. **$R_0$ convention** (`validate_r0_convention()`): the paper quotes $kT_\text{max}=50$ keV and $\eta_*=850$ for $\delta t_\text{ob}=0.5$ s at $z=4.35$. $R_0=c\,\delta t/(1+z)$ gives $\eta_*=812$ (−5%); $R_0=c\,\delta t$ gives 1235 (+45%). The $(1+z)$ convention is the paper's.
5. **Abdo et al. (2009) Table 1**, from the literature-review folder. The five epochs' Band parameters evaluated with the project's `band_function` reproduce Abdo's 50–300 keV photon fluxes to ≤1.1% (epoch b: 5.57 vs 5.63) and the total 10 keV–10 GeV fluence to 1% (2.43×10⁻⁴ vs 2.4×10⁻⁴ erg cm⁻²).
   This checks `band_function`'s normalisation and $E_\text{peak}$ convention, **not** the project's fits.
6. **Epoch b reproduction** (`validate_epoch_b()`; Abdo's own spectrum for 3.58–7.68 s, $\delta t_\text{ob}=0.5$ s, $z=4.35$). The paper's two quoted numbers imply $L_w=2.1\times10^{54}$ (from $kT_\text{max}$) and $2.5\times10^{54}$ erg/s (from $\eta_*$). Results by $L_\gamma$ band:

   | $L_\gamma$ band | $L_\gamma$ [erg/s] | kT_max (paper 50) | η* (paper 850) | σ_min (paper ≈20; 15–20 elsewhere) | $E_\text{iso}$, 5 epochs (paper ≈9×10⁵⁴) |
   |---|---|---|---|---|---|
   | observer 10 keV – 10 GeV | 2.68e54 | 53.3 | 866 | 22.2 | 8.87e54 |
   | observer 8 keV – 40 MeV | 1.96e54 | 49.3 | 800 | 15.9 | 6.90e54 |
   | observer 1 keV – 10 MeV | 1.61e54 | 46.9 | 762 | 12.9 | 5.99e54 |
   | rest-frame 1 – 10⁴ keV | 1.01e54 | 41.7 | 678 | 7.7 | 4.35e54 |

   The rest-frame (Amati) band is excluded. The observer-frame 10 keV – 10 GeV band also matches $E_\text{iso}$ and σ; 8 keV – 40 MeV matches $kT_\text{max}$ and $\eta_*$ just as well (they scale as $L^{1/4}$) and gives σ in the paper's 15–20 range.
   Zhang & Pe'er do not state their $L_\gamma$ band; 10 keV – 10 GeV is what Abdo et al. quote for the fluence, and it is inferred, not documented.
   `validate_epoch_b()` asserts the 10 keV – 10 GeV result against $kT_\text{max}$ (10%), $\eta_*$ (10%), σ (±25% of 20) and the Table 1 flux (3%).
   Before Table 1 was available, GRB080916C TR2 (4.9–15 s, the nearest project interval) with the same $\delta t$ gave σ_min = 6.5 on the rest-frame band, a factor 2.3–3 below the paper; a scan of Band shapes at fixed $L_\gamma=2.1\times10^{54}$ erg/s gave
   $1+\sigma$ = 3–11, so the interval and the spectrum together explain the gap.

### 6.3 $L_\gamma$ band

Primary (`l_band = zp09`): observer-frame 10 keV – 10 GeV. Sensitivity (`rest_1_1e4`): rest-frame 1 – 10⁴ keV (the Amati $E_\text{iso}$ band). Because the project's fits are joint GBM + LAT, the 10 MeV – 10 GeV part of
$L_\gamma$ is partly constrained for the LAT-covered bursts, but between ~40 and 100 MeV neither instrument constrains it, and the high-energy index is loosely pinned. Fraction of the 1 keV – 10 GeV flux below 10 MeV:
0.54–0.83 for GRB080916C, 0.78–0.92 for GRB140206B (PL TR5 0.84, the high-energy index −2.16 extrapolated), 0.97–1.0 for GRB131014A and GRB231129C. Because σ_min is a *lower limit that grows with $L_\gamma$*,
an over-estimated extrapolation errs toward a stronger claim; the 8 keV – 40 MeV band is the data-bounded alternative. On the rest-frame band the GRB080916C values are lower by the per-episode factors in `tables/zhang_peer_tables.generated.md`.

**Exact scaling with $L_\gamma$.** At fixed $R_0$ (hence fixed $t_v$ and $z$) and spectrum, $F_\text{th}/T^4\propto R_0^2$ is independent of $L_\gamma$, so the whole σ-family of
blackbodies lies on one curve $F=aT^4$ and the tangent (maximal allowed) blackbody is the *same* for any $L_\gamma$; only the starting point on that curve moves. Hence
$1+\sigma_\text{min}\propto L_\gamma$ exactly; the check over every (episode, $z$, test band) pair, and the factor between the two bands per episode, are in `tables/zhang_peer_tables.generated.md`.
The choice of $L_\gamma$ band is therefore a pure multiplicative factor on $1+\sigma_\text{min}$, and the sensitivity row can be replaced by one sentence giving that factor.
`zhang_peer_epoch_b.png` shows it: the two suppressed curves coincide.

### 6.4 Test band

Primary (`band = fit`): observer-frame 10 keV – 10 MeV; sensitivity (`extrapolated`): 1 keV – 10 MeV. Neither is in Zhang & Pe'er (they say only that the blackbody must be "unobservable"; their Fig. 2 spans 10⁰–10⁷ keV).
The upper edge is irrelevant (a blackbody is exponentially suppressed there); the lower edge matters. 8 keV is the GBM edge used for the fluence, so 10 keV is not quite consistent with it.

### 6.5 Results (zp09 band, fit band; $1+\sigma_\text{min}$, median)

The results table (per burst and episode: $t_v$ and its source, the fraction of draws needing suppression and the fraction that are edge-limited, $1+\sigma_\text{min}$ in the fit and extrapolated test bands over the redshift grid, $kT_\text{max}$, $\eta_*$), the edge-limited rows and the $L_\gamma$ scaling check are in `tables/zhang_peer_tables.generated.md`, written by `tables/stage2_tables.py` (rerun it after `zhang_peer_suppression.py`; never edited or retyped by hand).

**Edge-limited rows:** the maximum blackbody/spectrum ratio sits exactly at the 10 keV lower edge because $kT_\text{max}$ is only 1–7 keV.
For those rows the answer is set by the Wien tail at the band edge and the extrapolated band can differ by a large factor, so they should be reported as band-limited, not as measurements.
The rows that are not edge-limited are bound by the spectrum inside the band and do not depend on the edge.
`binding_energy_keV` and `binding_at_lower_edge_fraction` are in the CSV; the generated file lists the flagged rows and the largest fit-to-extrapolated ratio.

**Blackbody detection status in those episodes** (`results.json`; ΔC vs BEST; UNSAFE models excluded per §3a): no blackbody is detected in any of them. GRB080916C TR3 has SAFE SBPL_BB (ΔC 24.1) and BAND_BB (21.6) below the 28.74 threshold,
GRB140206B T90 a SAFE BAND_BB (24.9); GRB080916C EX1, TR5 and GRB140206B TR5 have only UNSAFE candidates; GRB140206B TR6 has no BB-bearing partner in the stored fits. A note for the paper: the project's test is
statistical detection, while σ_min is a deterministic comparison with the best-fit curve and is therefore a *floor* — a blackbody equal to the best-fit curve at its peak would be detected immediately, so a statistical-exclusion version would give a larger σ (§8).

### 6.6 Diagnostic figures (`zhang_peer_suppression_plot.py`, reads the CSV; validation figures, not yet in the paper)

- `zhang_peer_epoch_b.png/.pdf`: Abdo et al. epoch-b Band spectrum with the fireball blackbody at σ = 0 and σ_min for both $L_\gamma$ bands; the open circle is the binding energy (70 keV). The two suppressed curves (kT ≈ 24 keV) coincide (§6.3).
- `zhang_peer_suppression_spectra.png/.pdf`: one panel per non-BB episode ($z=4.35$ for GRB080916C, 2 for the others), best-fit spectrum, σ = 0 and σ_min blackbodies (primary configuration), edge-limited rows marked.
  In the edge-limited panels (GRB140206B T90, TR5, TR6; GRB080916C TR3, EX1) the "suppressed" blackbody peaks at or below the 10 keV edge and lies *above* the best-fit spectrum below it (GRB140206B T90: by roughly an order of magnitude at 4 keV), because the criterion is only tested from 10 keV.
  The best-fit curves use the stored best-fit parameters and the blackbodies the CSV medians, so they can differ from per-draw solutions by Monte Carlo spread.

### 6.7 Sensitivity to $t_v$ (user request, 2026-10-25)

Five episodes use their duration as $t_v$ (`t_v_is_duration_fallback`). For $\Gamma_\text{min}$ that is a valid upper-bound substitute ($\Gamma_\text{min}\propto t_v^{-1/(2\alpha+2)}$); for $1+\sigma_\text{min}$, where $t_v$ is a physical input through $R_0=c\,t_v/(1+z)$, nothing makes it a floor.
`tv_sensitivity_stage2.py` (validation, not a pipeline script; 1000 MC samples, the same random stream at every scale) rescales $R_0$ by 1, 0.3 and 0.1 and writes `tv_sensitivity_stage2.csv`. Primary $L_\gamma$ band, fit test band, $1+\sigma_\text{min}$ at scale 0.1 and 0.3 divided by scale 1; $kT_\text{max}$ rises by $\sqrt{10}$ at 0.1 by construction:

| episode (duration-based) | edge-limited at scale 1 | 0.3 / 1 | 0.1 / 1 |
|---|---|---|---|
| GRB080916C TR4, EX1, TR5 | partly (TR4: 26%), yes (EX1), no (TR5) | 0.41–0.60 | 0.19–0.36 |
| GRB140206B T90 | yes, all $z$ | 2.8–3.8 | 3.8–8.9 |
| GRB140206B TR5 | yes, all $z$ | 1.6–2.7 | 1.7–3.8 |

For the episodes with a measured $t_v$ the ratio at 0.1 ranges over 0.21–1.7 and also changes sign with $z$ (GRB231129C EX1 and TR2: 0.43–0.48 at $z=1$, 1.2–1.4 at $z=7$). The direction is therefore episode- and redshift-dependent: shortening $t_v$ lowers $1+\sigma_\text{min}$ for every GRB080916C episode and raises it for GRB140206B T90, TR5 and (at low $z$) TR6. Being edge-limited does not explain the sign, since GRB080916C TR3 and EX1 are edge-limited and fall; $kT_\text{max}$ is $\approx$6–10 keV there and $\approx$1–2 keV for GRB140206B T90 and TR5. The mechanism has not been worked out.
Consequence for the paper: the duration-based entries are neither floors nor ceilings, only an indication of the scale, and the $t_v$ footnote and section 5 text say so. The 1000-sample scale-1 values agree with the production CSV to Monte Carlo noise (e.g. GRB080916C TR2 19.2 against 19.4, TR3 63.3 against 63.4).

## 7. What the paper has to say (wording constraints from the above)

- Composition results are **conditional on $r_0$** (stage 1; $1+\sigma_0\propto r_0^{2/3}/f_\gamma$) and, for three bursts, on the **assumed redshift**; the $\sigma_0=0$ pure-fireball $r_0$ is the reference.
- The main table contains BEST-BB episodes only; blackbody models that are not BEST are reported separately and flagged (§3a); UNSAFE blackbody models are never used.
- Stage 2 is a **floor** ("the minimum Poynting-flux fraction for the fireball blackbody to fall below the best-fit spectrum"), conditional on the $\eta\ge\eta_*$ branch, $L_w=L_\gamma$, $R_0=c\,t_v/(1+z)$, and on the stated $L_\gamma$ and test bands;
  edge-limited rows (§6.5) and duration-fallback $t_v$ episodes are labelled.
- Stage 1 takes $f_\gamma=0.5$ ($L_w=2L_\gamma$); stage 2 takes $L_w=L_\gamma$ as Zhang & Pe'er do. The two stages' luminosity assumptions differ and the text must say so (for stage 2, $L_w=2L_\gamma$ would raise $kT_\text{max}$ by $2^{1/4}$ and so σ_min).
- Stage 1's $F_\text{BB}$ and $F_\text{ob}$ are observer-frame 10 keV – 10 GeV, not the Phase 1/2 band (1 keV – 10 MeV), and the text must say so. For GRB080916C, where 17–46% of the 1 keV – 10 GeV flux lies above 10 MeV, the wider band raised $1+\sigma_0$ by 20–50% at $r_0=10^8$ cm
  (T90 9.6 → 14.2, EX0 7.3 → 8.8, TR1 8.2 → 10.2); GRB140206B rose by about 16% and the other two bursts by 3% or less. The $L_\gamma$ band of stage 2 is a different quantity from the Amati $E_\text{iso}$ (rest-frame 1 – 10⁴ keV).

## 8. Open items and decisions

**Open (user's call; the defaults in use are stated):**
1. *(Decided 2026-10-25: 10 keV, no change to the CSVs.)* *Lower edge of the stage 2 test band*: 10 keV in use; 8 keV is the GBM edge used for the fluence. It only changes the 17 edge-limited rows.
2. *(Primary band decided 2026-10-25: observer-frame 10 keV – 10 GeV, no change to the CSVs. Sensitivity row decided the same day: replace the rows by one sentence in the paper giving the proportionality and the factor; they stay in the CSV.)* *Stage 2 sensitivity row for $L_\gamma$*: the rest-frame 1 – 10⁴ keV row is in use; the alternative is observer-frame 8 keV – 40 MeV (data-bounded, conservative, also consistent with the reproduction). The rest-frame row needs extrapolation to 0.125 keV at $z=7$. Since $1+\sigma_\text{min}\propto L_\gamma$ exactly (§6.3), either row is a pure rescaling and a one-sentence factor can replace the extra CSV rows.
3. *(Decided 2026-10-25: stage 1 uses 10 keV – 10 GeV, Phase 1/2 unchanged, anchor recomputed. Code changed and `gao_zhang_inversion.py` and `gao_zhang_inversion_plot.py` rerun at the full 10 000 samples on 2026-10-25: both CSVs and all four figures are the 10 keV – 10 GeV run. GRB080916C $1+\sigma_0$ at $r_0=10^8$ cm is now 14.2 (T90), 8.8 (EX0), 10.2 (TR1) against 9.6, 7.3, 8.2 before; no regime change among the 159 main-table `GZ2015` rows.)* *Stage 1 $F_\text{ob}$ band*: was observer-frame 1 keV – 10 MeV for consistency with Phase 1/2 (published in the paper).
4. *(LaTeX tables generated 2026-10-25 by `tables/stage1_to_latex.py` and `tables/stage2_to_latex.py`: main, redshift-sweep and partner tables for stage 1 and one table for stage 2, written next to the scripts and registered as `unused` in `table_registry.yaml`. Their citation keys `GaoZhang2015` and `ZhangPeer2009` were added to `ref.bib` the same day (title, authors, journal and DOI checked against the arXiv pages). The earlier source line in this note and in `gao_zhang_inversion.py` read "Gao, He & Zhang (2015), ApJ 810, 160 (arXiv:1507.01604)"; arXiv:1507.01604 is an unrelated paper and the reference is Gao & Zhang (2015), ApJ 801, 103 (arXiv:1409.3584), now corrected. The section, table and equation numbers cited from it were not re-checked against the PDF. The tables follow the model-set rules of §3a; the choice of columns, of the three `GZ2015` $r_0$ values and of $r_0=10^9$ cm for the sweep is Claude's and not yet the user's decision. Paper content decided 2026-10-25: main, sweep and stage 2 tables in the main text, the partner table in an appendix, the main figure only, no kT or $t_v$ columns (captions cite the photospheric and $\Gamma_\text{min}$ tables).)* LaTeX table(s) and figure (design: one marker per episode, $1+\sigma_0$ vs $r_0$ with the pure-fireball threshold marked, in the style of `pe_er_photosphere.py`; the partner results and stage 2 as separate, flagged tables), and the paper text.
5. A statistical-exclusion version of stage 2: replace "below the best-fit curve" by "below the largest blackbody the data allow at that temperature" (a profile-likelihood upper limit from frozen-blackbody refits in counts space). It needs the detector responses and RMFIT, which are not in the repository; the same upper limit would turn the stage 1 partner results into proper lower limits on $1+\sigma_0$. Recorded as a future experimental setup.
6. `seed_table_to_latex.py` not rerun (both scripts are registered `unused` in `seed_registry.yaml`, so the paper's seed table is unchanged but `seed_table.csv` is stale).

**Decided (user, 2026-10-24):**
- *δ sensitivity: option (a).* The paper cites Gao & Zhang's §6 statement that $T_\text{ob}$ and $F_\text{BB}$ at δ = 1/3 are upper limits for δ < 1/3 over the range relevant to GRBs; Phase 6 carries no numerical δ check. **Option (b)**, deriving general-δ top-down formulae
  (new, unvalidated algebra; Gao & Zhang's Table 3 gives the general-δ regime criteria in $(\eta,\sigma_0,L_w)$ to check against), is a possible later experimental setup, an `experiments/`-style study outside the paper's pipeline; nothing is built for it.
- *$f_\gamma=0.5$.* Confirmed.
- *Stage 2 defaults:* $\eta\ge\eta_*$ branch, $R_0=c\,t_v/(1+z)$, fit band from 10 keV with the 1 keV extension as sensitivity.
- *Model sets:* §3a (all three rules and the SAFE/MARGINAL admission confirmed).

## 9. Known limitations

- $f_\gamma$, $\delta$ and $r_0$ are assumed, not measured; $1+\sigma_0\propto1/f_\gamma$ exactly.
- A 2–3% systematic floor on $1+\sigma_0$ near 1 from the paper's rounded coefficients (§4).
- Regime selection on the median inputs; 21 rows have a minority of draws in another regime (reported via `regime_fraction`).
- Observer-frame 10 keV – 10 GeV band for $F_\text{BB}$ and $F_\text{ob}$ (from 2026-10-25), no longer the Phase 2 band, so stage 1's $F_\text{BB}$, $f_\text{th}$ and $kT$ are not the quantities in the published Pe'er+07 table; the paper does not specify a band. For the LAT-covered bursts the 10 MeV – 10 GeV part is an extrapolation of the joint-fit model where the data do not constrain it (§6.3).
- Three bursts use assumed redshifts; none of their composition numbers are measurements.
- MC draws are independent of `pe_er_photosphere.py`'s (own seed), so $kT$ and $F_\text{BB}$ differ from its medians by ordinary MC noise.
- Stage 2: $\sigma_\text{min}$ compares the suppressed blackbody with the best-fit curve, not with the data's statistical sensitivity (§6.5, §8 item 5); the $\eta<\eta_*$ branch is not computed.
- Stage 2: $L_\gamma$ over 10 keV – 10 GeV extrapolates the joint-fit model beyond the constrained energies for episodes with shallow high-energy indices (§6.3).
- Stage 1 partner models: the blackbody amplitude is the best-fit value of a component that is not required; treat the results as conditional (§5a).
