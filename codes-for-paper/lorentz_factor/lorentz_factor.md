# Minimum bulk Lorentz factor from $\gamma\gamma$ opacity — method notes

Companion to `lorentz_factor.py` and `lorentz_factor_limit_b.py`. Produced 2026-08-21 (Phase 0 fixes, extended to per-episode in Phase 2; Limit B added and integrated into the paper during Phase 4 QA — §8).

---

## 1. What is computed

The lower limit on the bulk Lorentz factor set by pair-production transparency, following Lithwick & Sari (2001), ApJ **555**, 540, "Limit A", with redshift corrections:

$$\hat\tau = 2.1\times10^{11}\,\frac{(d_L/7\,\text{Gpc})^2\,(0.511)^{-\alpha+1} f_1}{(\Delta T/0.1\,\text{s})(\alpha-1)}$$

$$\Gamma_\text{min} = \hat\tau^{\frac{1}{2\alpha+2}}\left(\frac{E_\text{max}}{511\,\text{keV}}\right)^{\frac{\alpha-1}{2\alpha+2}}(1+z)^{\frac{\alpha-1}{\alpha+1}}$$

where $\alpha = -\beta$ is the (positive) high-energy photon index, $f_1$ the photon flux at 1 MeV, $E_\text{max}$ the highest-energy LAT photon, and $\Delta T$ the variability timescale.

**This is a one-sided bound, not an estimate.** We observed a high-energy LAT photon; had $\Gamma$ been below $\Gamma_\text{min}$, that photon would have pair-produced on the burst's own softer photons and never reached us. So the measurement establishes $\Gamma \geq \Gamma_\text{min}$ and nothing more.

---

## 2. Decisions and who made them

### 2.1 Per-episode evaluation — *user, Phase 2*

Originally computed for the T90 interval only. Because $\Gamma_\text{min}$ grows with photon flux, a burst-averaged spectrum gives a systematically **weak** bound. It is now evaluated for every episode with LAT coverage, so it can be compared with the Phase 2 thermal $\Gamma$ on *identical* interval boundaries rather than against an external result with different binning (`BUGS.md` OBS-07).

### 2.2 Inputs read from the fitted-model database — *user, Phase 0*

$f_1$ and $\beta$ come from `results.json` through the `grb_research` class API, not from hardcoded numbers. Only quantities absent from that database are tabulated in the script: the per-episode LAT photon energy and arrival time (from the paper's LAT appendix), the redshifts, and any literature variability timescale.

### 2.3 $f_1$ from the actual spectrum, non-thermal only — *Claude, Phase 0*

$f_1$ is the continuum photon flux evaluated **at 1 MeV on the fitted model**, with any blackbody component excluded — the opacity is set by the power-law photons that pair-produce with the LAT photon, not by the thermal component.

The previous implementation extrapolated from the 100 keV pivot using the steep `index2`, which is wrong for SBPL (whose amplitude is defined at the pivot, where the *shallow* index applies). That underestimated $f_1$ by 5.9× and $\Gamma_\text{min}$ by 1.3×; see `BUGS.md` BUG-02.

### 2.4 Variability timescale — *user; see §4 for the decision and its consequences*

### 2.5 Where the LAT photon data comes from — *and a unit bug found on the way*

The per-episode $E_\text{max}$ and arrival times in `LAT_PHOTONS` were originally transcribed by hand from the paper's own LAT appendix, `GRBResearchPaper/appendices/appendix_LAT_info.tex` (`tab:burst_table`).

**Bug found while transcribing (historical).** That table's photon-energy column was headed "(MeV)", but its values include `27428.80` — which in MeV would be 27 TeV. The paper text and this script both treat them as keV (27428.80 keV = 27.43 GeV), so the header was wrong. Fixed to (keV); logged as `BUGS.md` BUG-16.

**OBS-08 fixed (status current as of this session, 2026-09-16).** The hand-transcription risk flagged here has since been resolved: `lorentz_factor.py` no longer hardcodes a `LAT_PHOTONS` dict. It now loads via `load_lat_photons()` from `LAT_analysis/lat_photons.csv`, which is itself built from the `gtburst`/`gtlike` output files by `LAT_analysis/txt_to_csv.py` — a machine-derived CSV, not a copy of the paper's `.tex` table. Photon energies in that CSV are MeV (per its own `e_max_MeV` column, with the LAT selection floor at 100 MeV; see `BUGS.md` BUG-18 for the earlier keV/MeV mixup this project hit elsewhere). Check `LAT_analysis/LAT_analysis.md` for that pipeline's own provenance/validation notes rather than duplicating them here.

### 2.6 Low-significance detections flagged, not dropped — *Claude*

Five episodes have LAT detections with TS < 25 (footnotes a–e of the appendix table): GRB080916C TR5, GRB140206B TR5/TR6, GRB231129C EX0/TR1. Their "highest-energy photon" is not a secure association, so the derived $\Gamma_\text{min}$ carries a $\ddagger$ in the table rather than being silently tabulated as equivalent.

---

## 3. Results

**STALE — kept for historical narrative only, do not read these numbers as current.** Confirmed 2026-09-16: the $\Gamma_\text{min}$ values below (134/78/86/129/137/80/83/71) are the pre-BUG-18 numbers. The current `lorentz_results.csv` (and §8.3's Limit A column) instead give 507/350/380/447/504/295/287/260 for the same episodes in the same order — see §8.3's own note ("Limit A values here are post-BUG-18 ... not the pre-fix numbers in §3 above"), which this session's check confirms still describes the live discrepancy correctly. The thermal-$\Gamma$ comparison and $Y$-constraint below were computed from the *stale* column and have not been redone with the current values; the qualitative conclusion (bounds too weak to constrain $Y$) almost certainly still holds since $(350/852)^4\approx0.03$ is still $\ll1$, but the quoted figure ($7\times10^{-5}$) is off by roughly two orders of magnitude and should be recomputed before being cited again.

Only GRB080916C has a spectroscopic redshift, so only it yields values:

All rows use the episode duration as $t_v$ (§4). Errors are statistical only — see §5.

| Episode | $t_v$ [s] | $E_\text{max}$ [keV] | $\Gamma_\text{min}$ |
|---|---|---|---|
| T90 | 62.976 | 27428.8 | $134 \pm 1$ |
| EX0 | 4.992 | 301.2 | $78^{+5}_{-4}$ |
| TR1 | 3.584 | 301.2 | $86^{+5}_{-5}$ |
| TR2 | 10.176 | 2110.1 | $129^{+2}_{-2}$ |
| TR3 | 40.256 | 27428.8 | $137^{+1}_{-1}$ |
| TR4 | 4.224 | 464.5 | $80^{+6}_{-6}$ |
| EX1 | 8.384 | 989.4 | $83^{+3}_{-3}$ |
| TR5 | 4.736 | 340.0 | $71^{+6}_{-5}\,^\ddagger$ |

**Comparison with the Phase 2 thermal $\Gamma$**, now on identical boundaries:

| Episode | $\Gamma_\text{min}$ (opacity) | $\Gamma$ (thermal, $Y{=}1$) | consistent? |
|---|---|---|---|
| T90 | 134 | 752 | yes |
| EX0 | 78 | 852 | yes |
| TR1 | 86 | 861 | yes |

The thermal $\Gamma$ clears the opacity floor everywhere. **But the bounds are far too weak to constrain $Y$**: from EX0, $Y \geq (78/852)^4 \approx 7\times10^{-5}$, which is no constraint at all given $Y \geq 1$ by definition.

So the earlier suggestion that the two methods together imply $Y \gtrsim 1.9$ **does not survive** once the comparison is made internally. That figure came entirely from Abdo et al.'s much stricter published bound (887 ± 21), computed on different time bins and evidently with a stricter treatment than Lithwick & Sari Limit A.

---

## 4. The variability timescale — decision and rationale

**Decided by the user, 2026-08-21: use each episode's own duration, uniformly. This must be stated explicitly in the paper's methods text.**

### What the problem was

Three conventions were in play and gave different answers for the same burst:

| convention | $t_v$ for T90 | $\Gamma_\text{min}$ |
|---|---|---|
| the code as inherited (a literature value) | 0.9 s | 258 |
| what `section-5-data-analysis.tex` *said* — "the duration of the shortest well-resolved emission episode" | 3.584 s (TR1) | 209 |
| **each interval's own duration (adopted)** | 62.976 s | **134** |

Two separate defects. First, the code and the text disagreed: 0.9 s corresponds to no episode of GRB080916C, so §5 did not describe what the code actually did. Second, once the calculation went per-episode, T90 was using a literature value while every other row used its duration — **two conventions in one column**.

### Why the duration was chosen

- **It is honest about direction.** The episode duration is an *upper bound* on the true variability timescale. Since $\Gamma_\text{min} \propto \Delta T^{-1/(2\alpha+2)}$, an over-estimated $\Delta T$ yields an under-estimated $\Gamma_\text{min}$ — i.e. a conservative lower limit that remains valid. Given our bounds are already weaker than published treatments (§5), erring toward conservatism is the defensible posture.
- **It is uniform.** Every row is computed the same way, so episodes can be compared with each other and with the Phase 2 thermal $\Gamma$ without a footnote explaining which rows differ.
- **It is self-contained.** The duration is read from the interval itself, so nothing has to be sourced externally for the three bursts that have no published variability timescale.
- **The sensitivity is weak.** The exponent is $\approx -0.15$, so even the factor-70 spread between 0.9 s and 62.976 s moves $\Gamma_\text{min}$ by less than a factor of two (258 → 134).

*Rejected:* keeping 0.9 s for T90 preserved the previously quoted number but left the table internally inconsistent and the text describing neither convention.

### What this obliges the paper to say

$\Gamma_\text{min}$ for GRB080916C's T90 interval changes from **258 to 134**, and §5's description of $t_v$ must be rewritten — the "shortest well-resolved emission episode" wording is no longer what is done. The text should state that $t_v$ is taken as the episode duration and that this makes the limit conservative. `VARIABILITY_TIMESCALE` in the script remains as an override map should a published per-episode value later be adopted.

---

## 5. Uncertainties — statistical vs systematic

**Added 2026-08-21 at the user's prompting: $\Gamma_\text{min}$ originally carried no error bars at all.** It was the one calculation in the project still evaluated at best-fit values only, while everything else propagated the fit covariance. That was an oversight, and it violated the project rule requiring error bars wherever an MC uncertainty exists.

### Statistical error (now computed)

The fit covariance is propagated through **both** $f_1$ and $\alpha$ from the same 10 000 draws, so their correlation is preserved — they are not independent, since both derive from the same spectral parameters. Draws with $\alpha \leq 1$ are rejected, as the Lithwick & Sari expression is singular there.

**Same §3 staleness applies here**: the point estimates below are the pre-BUG-18 values (see §3's note); the *percentage* errors are still representative of the current fit's relative precision (the fix was a systematic shift in $f_1$, not a change in how well-constrained it is), but the absolute $\Gamma_\text{min}$ values themselves should be read from §8.3/`lorentz_results.csv`, not from this table.

| Episode | $\Gamma_\text{min}$ | stat. error |
|---|---|---|
| T90 | 134 | $\pm 1$ (0.7%) |
| EX0 | 78 | $\pm 4.5$ (5.8%) |
| TR1 | 86 | $\pm 5.1$ (6%) |
| TR2 | 129 | $\pm 1.9$ (1.4%) |
| TR3 | 137 | $\pm 1.1$ (0.8%) |
| TR4 | 80 | $\pm 6$ (7.5%) |

### Why they are so small — and why that is misleading

Every input enters through a fractional power:

$$\Gamma_\text{min} \propto f_1^{\frac{1}{2\alpha+2}} \, E_\text{max}^{\frac{\alpha-1}{2\alpha+2}} \, \Delta T^{-\frac{1}{2\alpha+2}}$$

With $\alpha \approx 2.2$ the exponent on $f_1$ and $\Delta T$ is only $\approx 0.15$, and on $E_\text{max}$ only $\approx 0.19$. So a 10% error in the photon flux moves $\Gamma_\text{min}$ by 1.5%. The calculation is intrinsically insensitive, which is exactly why it makes a robust *bound* — and exactly why tight error bars on it mean very little.

**The systematic uncertainty is one to two orders of magnitude larger:**

| source | effect on $\Gamma_\text{min}$ |
|---|---|
| statistical (fit covariance) | 1–8% |
| choice of $t_v$ convention (§4: 0.9 s / 3.584 s / 62.976 s) | factor $\approx 2$ (134 ↔ 258) |
| analytic approximation — Limit A vs. a fuller treatment | factor $\approx 7$ (our 129 vs. Abdo et al.'s 887 on an overlapping interval) |
| LAT energy resolution on $E_\text{max}$ (~10%, not propagated — no per-photon errors available) | ~2% |

**Consequence for the paper: quoting $134^{+1}_{-1}$ without qualification would be actively misleading.** The table caption now states explicitly that these are statistical only and are far smaller than the systematic uncertainty. The honest reading is that $\Gamma_\text{min}$ is known to within a factor of a few, not to 1%.

---

## 6. Limitations

- **Our bounds are much weaker than published ones for the same burst.** Abdo et al. (2009) obtain $\Gamma_\text{min} = 887 \pm 21$ (their bin b, 3.6–7.7 s) and $608 \pm 15$ (bin d); our TR2, which contains the same ~2 GeV photon, gives 129 pre-BUG-18 / **447 in the current, post-BUG-18 `lorentz_results.csv`** (see §3's staleness note). Even at the corrected value the gap to Abdo et al. remains large and is still methodological, not explained by $\Delta T$ alone — Lithwick & Sari's Limit A is a simplified analytic form and Abdo et al. use a fuller treatment. **Our values should be presented as conservative lower limits, not as competitive with published ones.**
- **Only GRB080916C yields values**, since the other three lack redshifts.
- **T90 rows are not independent** of the TR rows that tile them.
- Episodes whose best-fit model has $\beta \geq -1$ are skipped: $\alpha \leq 1$ makes the expression singular.

---

## 6b. A better $t_v$: the Norris-profile fit — **assessed 2026-08-21, feasible; measurement now in progress, see status note**

**Status update (2026-09-16):** this section's "not yet done" is now only true of *this file* (`VARIABILITY_TIMESCALE` below is still empty and $\Gamma_\text{min}$ is still computed at episode duration). The actual measurement work described here started under `PHASE5_TV_PLAN.md` (project root): an automated MEPSA-based pipeline was tried and abandoned (stalled on GRB080916C's TR2), and a manual, GRB-by-GRB joint-Norris-fit track is now underway in `codes-for-paper/variability_analysis/` (see that folder's `variability_analysis.md`), covering GRB131014A, GRB140206B, and GRB231129C so far — GRB080916C (the burst that actually feeds a published $\Gamma_\text{min}$ here) does not yet have a manual fit. None of this has been fed back into `VARIABILITY_TIMESCALE` or this script; that remains a separate, explicit follow-up per `PHASE5_TV_PLAN.md`'s own deliverable-boundary decision.

Raised by the user after re-reading \citet{Bukhari2022} (`GRBResearchPaper/Literature Review/2022_06_06_Dr_Saeeda_Sajjad_Syed_Ali_Mohsin_Bukhari_Urooj_Murtaza.pdf`), §6.5 of which derives $t_v$ by fitting the light curve rather than adopting an interval duration.

### The method

They fit a single pulse of the 10–400 keV light curve with the Norris et al. (2005) profile (their eq. 9):

$$I(t) = F \exp\!\left[2\left(\frac{\tau_1}{\tau_2}\right)^{1/2}\right] \exp\!\left[-\frac{\tau_1}{t-t_s} - \frac{t-t_s}{\tau_2}\right], \qquad t > t_s$$

with free parameters $F$, $t_s$, $\tau_1$, $\tau_2$. For GRB110721A they fit the 1.472–2.56 s window and obtained $t_v = 0.7 \pm 0.02$ s. They repeated the fit in 400–900 keV and 250 keV–5 MeV, getting $0.8 \pm 1.0$ and $0.8 \pm 0.05$ s, and used the agreement to argue $R_\text{GeV} \simeq R_\text{MeV}$ — a second, independent use of the same fit.

### Is it feasible for us? **Yes — the data is in hand and is good.**

`GRBResearchWork/light_curves/` already holds background-subtracted, energy-resolved RMFIT `.dat` files for all four bursts, with a reader (`lightcurve_data()`) and per-GRB driver scripts:

| GRB | detectors | bins | binning | peak rate [s$^{-1}$] |
|---|---|---|---|---|
| GRB080916C | b0, n3, n4 | 5103 | 64 ms | 2 915 |
| GRB131014A | b1, n9, na, nb | 9600 | 64 ms | 44 341 |
| GRB140206B | b0, n0, n1, n3 | 9600 | 64 ms | 5 585 |
| GRB231129C | b0, n3, n6, n7 | 9600 | 64 ms | 12 785 |

64 ms binning resolves a $\sim$0.7 s pulse with $\sim$10 bins across the rise alone, so the fit would be well constrained. GRB131014A in particular is very bright.

### Why it is worth doing

**$t_v$ is currently the single largest systematic in $\Gamma_\text{min}$** — larger than everything in §5 except the analytic approximation itself. The three candidate conventions span a factor of two (134 / 209 / 258 for T90). We deliberately chose the most conservative, the episode duration, which is an *upper bound* on the true variability timescale and therefore *depresses* $\Gamma_\text{min}$.

A Norris fit would replace that upper bound with a measurement. Since $\Gamma_\text{min} \propto t_v^{-1/(2\alpha+2)}$, the effect is predictable: for TR3, moving $t_v$ from its 40.256 s duration to a Norris-like $\sim$0.7 s would raise $\Gamma_\text{min}$ by a factor $(40.256/0.7)^{0.154} \approx 1.87$ — from the pre-BUG-18 137 to roughly 256 as originally estimated here, or, applied to the current post-BUG-18 value (504, §3/§8.3), to roughly $504\times1.87\approx943$. That closes part — not all, in the pre-BUG-18 estimate, or plausibly all, at the current value — of the gap to \citet{Abdo2009FermiObservations080916C}'s 887; this illustrative number has not been recomputed properly with the current $f_1$/$\alpha$, only rescaled by the same factor, so treat it as indicative only. It would let us state that our limits are conservative *for a quantified reason* rather than by construction.

### Caveats to settle first

- **The Norris profile describes one pulse.** Several of our episodes are Bayesian-block intervals containing multiple pulses, and T90 spans the whole burst (63 s for GRB080916C) — a single-pulse fit there is meaningless. Pulse identification per episode is a prerequisite, and it is a judgement call, not a mechanical step.
- Which pulse defines $t_v$ for a multi-pulse episode needs a stated rule (brightest? narrowest? first?).
- \citet{Bukhari2022} used $H_0 = 67.4$, $\Omega_M = 0.315$ (Planck 2020), **not** this project's 69.6/0.286. Fine for borrowing the method; do not import their numbers.
- This changes a published number again ($\Gamma_\text{min}$ moved 197 → 258 → 134 already). Worth doing once, deliberately.

**Not required for the paper as it stands** — the current treatment is internally consistent and honestly labelled conservative. This is an improvement, not a correction.

---

## 7. Files

| file | role |
|---|---|
| `lorentz_factor.py` | Limit A computation, CSV and LaTeX table |
| `lorentz_results.csv` | one row per episode with LAT coverage (Limit A) |
| `lorentz_table.tex` | generated paper table (Limit A) |
| `lorentz_factor_limit_b.py` | Limit B computation, CSV and LaTeX table (§8) |
| `lorentz_results_limit_b.csv` | one row per episode with LAT coverage (Limit B) |
| `lorentz_table_limit_b.tex` | generated paper table (Limit B) |
| `gamma_comparison_plot.py` | comparison figure: Limit A, Limit B, thermal $\Gamma$ (§8.6) |
| `gamma_comparison.png` / `.pdf` | the figure itself, also copied to `GRBResearchPaper/images/section5/` |
| `lorentz_factor.md` | this file |

---

## 8. Limit B — Compton scattering off pair-produced $e^{\pm}$ — **added 2026-08-21**

Raised by the user after the Limit A cross-check against the paper PDF (§1). Implemented in a separate module, `lorentz_factor_limit_b.py`, with its own CSV and LaTeX table — deliberately not folded into `lorentz_factor.py`'s output, so a future edit to one cannot silently corrupt the other.

### 8.1 What is computed

Lithwick & Sari (2001), Table 2, eq. (8): the burst is optically thin if the $e^\pm$ pairs created by photon annihilation are themselves too few to Compton-scatter the burst's photons appreciably:

$$\Gamma_\text{min} = \hat\tau^{\frac{1}{\alpha+3}}(1+z)^{\frac{\alpha-1}{\alpha+3}}$$

using the *same* $\hat\tau$ as Limit A (eq. 4/9) — the two limits share the optical-depth bookkeeping quantity and differ only in how it's combined with $E_\text{max}$/$z$ afterward. **Limit B does not depend on $E_\text{max}$ at all** — no highest-energy photon needs to be identified or trusted, only the fitted continuum's normalisation and slope ($f_1$, $\alpha$) through $\hat\tau$.

### 8.2 Judgement calls

- **$\hat\tau$ factored out and shared, not re-derived — *Claude*.** `compute_tau_hat(alpha_LS, f_1, delta_T_s, z)` was extracted from `lorentz_factor.py::compute_gamma_min` and imported by both modules, so the formula exists in exactly one place. **Validated**: joining the two output CSVs on `GRB`+`episode` gives `tau_hat` bit-identical (max abs difference `0.0`) across all 26 rows, as it must for two functions computing the same closed-form expression from the same best-fit point values ($f_1$, $\alpha$) — the CSV's `tau_hat` column is deterministic, not MC-drawn, so this identity holds regardless of each script's seed (each now derives its own via `seed_from_name(__file__)`, deliberately decorrelated between Limit A and Limit B — see RNG/seeding overhaul, Phase B). Only `Gamma_min_err_lower/upper` depend on the seed.
- **Same episode set as Limit A — *Claude*, deliberate, revisit if wanted.** The loop still restricts to `episode in photons` (i.e. episodes with a LAT-catalogued photon at all), mirroring Limit A and the paper's own Table 3, which lists both limits side by side for the same bursts. This was **not** widened to "every BEST-model episode" even though Limit B's formula would technically allow it (no $E_\text{max}$ needed) — doing so would be a scope decision (what does "LAT coverage" mean when no photon energy is used?) that wasn't asked for. If broader coverage is wanted later, this is the line to revisit.
- **`low_significance` (TS < 25) flag dropped from the Limit B CSV/table — *Claude*.** That flag marks an insecure *photon-energy* association, which Limit B never uses. Carrying it into Limit B's output would misleadingly suggest a caveat that doesn't apply here — this is precisely the practical benefit of Limit B raised when it was proposed: it gives a usable bound for episodes where Limit A's photon association is shaky.
- **`E_max_MeV` and `t_arr_s` dropped from the Limit B CSV — *Claude***, per the CSV minimalism rule (`CLAUDE.md`): they are not inputs to this formula, so including them would be dead weight. A reader who wants to cross-reference against the observed photon can join on `GRB`+`episode` with `lorentz_results.csv`.
- **Only GRB080916C yields values — same redshift limitation as Limit A**, since $\hat\tau$ needs $d_L(z)$ regardless of which limit is taken afterward. Limit B does **not** unlock the other three bursts.

### 8.3 Results, GRB080916C

| Episode | Limit A $\Gamma_\text{min}$ | Limit B $\Gamma_\text{min}$ | larger |
|---|---|---|---|
| T90 | $507 \pm 1$ | $112 \pm 1$ | A |
| EX0 | $350^{+5}_{-4}$ | $190 \pm 8$ | A |
| TR1 | $380^{+5}_{-5}$ | $213^{+10}_{-9}$ | A |
| TR2 | $447^{+2}_{-2}$ | $190 \pm 3$ | A |
| TR3 | $504^{+1}_{-1}$ | $115 \pm 1$ | A |
| TR4 | $295^{+6}_{-6}$ | $154^{+10}_{-9}$ | A |
| EX1 | $287^{+3}_{-3}$ | $131 \pm 5$ | A |
| TR5 | $260^{+6}_{-5}\,^\ddagger$ | $141 \pm 9$ | A |

(Limit A values here are post-BUG-18, i.e. the corrected MeV-consistent table, not the pre-fix numbers in §3 above.) For this sample, **Limit A is the larger (binding) bound in every episode**, so $\max(A,B) = A$ throughout and this does not change the paper's reported $\Gamma_\text{min}$. This is a real result, not a wasted computation: it confirms Limit A's dominance rather than assuming it, and Limit B remains available as the fallback for the low-significance-flagged episodes should Limit A's photon association ever be reconsidered.

### 8.4 Limitations

- Same $t_v$-convention and analytic-approximation caveats as Limit A (§4, §5, §6) apply identically, since both share every input except the final combination step.

### 8.5 Paper integration — **added 2026-08-21**

Integrated into `GRBResearchPaper/tex_files/section-5-data-analysis.tex`, immediately after the existing Limit A paragraph (`sec:lorentz`), as new prose rather than a silent replacement of the Limit A table:

- $\Gamma_{\min,B}$ given as a numbered, labelled equation (`eq:gamma_min_limitB`), reusing the already-defined $\hat\tau$ (`eq:tau_hat`) via `\cref`, matching this file's own established citation/cross-reference style.
- `lorentz_table_limit_b.tex` copied manually into `GRBResearchPaper/tex_files/generated/` (the two repos are independent, per `CLAUDE.md` — this copy step does not happen automatically) and `\input`, with the same `% Generated by ... — do not hand-edit` header convention as the Limit A table.
- Prose states the actual result (§8.3): Limit A dominates in every episode of this sample, so the paper's adopted $\Gamma_\text{min}$ does not change; Limit B is presented as an independent cross-check, plus a forward-looking note that it would be the fallback for a low-significance-flagged episode, should one ever need it.
- **Verified, not assumed**: full `latexmk` rebuild (0 errors, 0 undefined refs/citations, 23→24 pages), followed by `pdftotext` grep confirming the equation renders as `Γmin,B = τ̂ 1/(α+3) (1 + z)(α−1)/(α+3)` and the new Table 4 renders with the exact same $\Gamma_{\min,B}$ values as `lorentz_results_limit_b.csv` (112, 190, 213, 190, 115, 154, 131, 141) — not just that the build succeeded, per the project's standing rule that a clean build proves nothing about content on its own.

### 8.6 Comparison figure — **added 2026-08-21**

Raised by the user: `lorentz_factor/` was the only Phase 1/2 topic folder without a plot. `gamma_comparison_plot.py` puts Limit A, Limit B, and the Phase 2 thermal $\Gamma$ (`photospheric_radius/pe_er_photosphere.csv`, $Y=1$) on one log-scaled axis, per episode of GRB080916C, full $1\sigma$ MC error bars on every point.

**Reads the three existing CSVs rather than recomputing — *Claude*, deliberate.** Every other plotting script in this project computes and plots from the same in-memory data (`amati_relationship.py`, `pe_er_photosphere.py`), but this figure spans three sibling scripts' outputs; recomputing any of the three formulas here would be a fourth independent implementation risking exactly the drift `CLAUDE.md`'s "Verify, don't assume" section and this file's own §8.2 warn about. Each of the three source CSVs is already self-contained and defensible per `CLAUDE.md`'s CSV standard, so reading them is not a shortcut.

**Episode markers/order come from live `TimeInterval` objects, not re-parsed CSV strings — *Claude*.** `collect_intervals()` calls `prepare_grbs` the same way `lorentz_factor.py::main()` does and resolves markers via the real `EpisodeMarkerResolver.resolve(interval)`, so marker shapes can never drift from what the resolver itself assigns. `episode_order()` (temporal x-axis ordering) is a small intentional duplicate of `photospheric_radius/pe_er_photosphere.py`'s function of the same name — a display-ordering convenience, not a formula, so the drift risk that justifies sharing `compute_tau_hat` doesn't apply here.

**Two real bugs, caught by the user in the rendered image, not by the script exiting cleanly:**
1. **In-axes legend covering data.** The first version placed the `Method` legend at `loc="upper left"` inside the axes; since the thermal-$\Gamma$ points are the highest values on the log axis, the legend box sat directly on top of them — the same class of defect as `BUGS.md` BUG-15 (`pe_er_photosphere.py`'s earlier legend/marker bug). Fixed by moving both legends outside the axes via `bbox_to_anchor`.
2. **Legend text clipped at the canvas edge.** Even after moving the legends outside, their text was cut off at the figure's right border across three separate widening attempts (figsize 8→11→13). Root cause: `update_style()` already sets `rcParams["savefig.bbox"]="tight"`, but a legend added via `axis.add_artist()` (as opposed to the axes' own current legend) is not reliably included in matplotlib's automatic tight-bbox artist search. Fixed with an explicit `bbox_extra_artists=(legend1, legend2)` passed to `savefig`. **Verified by cropping and re-inspecting the actual saved PNG pixels** (`PIL.Image.crop` on the right 1000–1200 px) after each attempt, not by re-running the script and assuming the fix worked — the first `bbox_inches="tight"` attempt looked identical to the broken version until checked this way.

**Final layout, per the user's follow-up request** ("previous size was correct; shrink the axes; use newline in legend"): `figsize=(8, 5.5)` (back down from the 13-wide version used mid-debugging), `figure.subplots_adjust(right=0.62)` to reserve room for the legends within that fixed canvas, and `METHOD_LABELS` given embedded `\n` line breaks so the long method descriptions don't force the legend box wide in the first place.

**Result, verified numerically:** thermal $\Gamma$ exceeds Limit A by $1.5$–$2.4\times$ and Limit B by $4.0$–$6.7\times$ across T90/EX0/TR1 (the only BB-inclusive episodes); Limit A exceeds Limit B in all 8 episodes with LAT coverage. The $4.0$–$6.7\times$ figure was computed from the CSVs (`Gamma / Gamma_min_B`) before being quoted in the paper, not estimated from the plot.

**Paper integration**: copied to `GRBResearchPaper/images/section5/gamma_comparison.png`, inserted as Figure 9 (`fig:gamma_comparison`) in `section-5-data-analysis.tex` right after the Limit A/B prose, with a cross-reference added from §6's pre-existing thermal-vs-opacity sentence (which now also states the Limit B ratio). **Verified without `draft` mode** — a draft build stubs figures and proves nothing about whether the image path resolves — using a throwaway `out_nodraft/` output directory (removed after): 0 missing-figure warnings, `pdftotext`-confirmed caption text, and the actual rendered page 12 visually inspected as a PNG to confirm both legends display fully and don't overlap data. `main.tex`'s `draft` flag and the normal `out/` build were restored/rebuilt afterward, unchanged from the documented workflow.

---

## 9. Table formatting fixes — **added 2026-08-21, at the user's direction**

The user hand-edited `GRBResearchPaper/tex_files/generated/lorentz_table_limit_b.tex` directly to fix a width problem, then asked for the same two fixes to be made properly in both generator scripts (`lorentz_factor.py` and `lorentz_factor_limit_b.py`), so the tables stay generated rather than diverging from their scripts on the next run:

1. **Wrap the table in `threeparttable`** (already loaded via `preamble.tex`), with footnotes in a `tablenotes` block rather than as `\multicolumn` rows inside the `tabular`. Applied identically to both `build_latex_table()` functions.
2. **Drop bursts without a redshift from the table entirely**, rather than showing them as all-`\ldots` rows. `results = [r for r in results if r["z"] is not None]` at the top of each `build_latex_table()` — the CSV is untouched and still carries every row (all four bursts), this is a presentation-only filter. Since $z$ is set per-burst (not per-episode) and only GRB080916C has one, this cleanly drops exactly the three bursts with no values to show.

**One consequence of (2) that needed a deliberate call, not just a mechanical filter — *Claude*.** The user's manual edit moved the "$^\dagger$ episode duration adopted as $t_v$" marker from a per-row superscript to a single marker in the column header, since every row now left in the table happens to be duration-sourced. Replicating that verbatim would hardcode an assumption that stops being true the moment `VARIABILITY_TIMESCALE` (currently empty) gets an entry for one of these episodes — silently mislabeling a literature-sourced $t_v$ as duration-derived. Rather than silently trust it, both scripts now **assert** `{r["t_v_source"] for r in results} == {"duration"}` before building the header this way, and raise with a clear message telling the reader to move the dagger back to per-cell if it ever fires. This is the same pattern as `LAT_analysis/csv_to_latex.py`'s `FLUX_UPPER_LIMITS` assertion: an input that could plausibly go stale gets a loud failure mode, not a silent wrong table.

**Verified, not assumed:** regenerated both tables from the fixed scripts and diffed the Limit B output against the user's manual edit — identical in substance (all values, `threeparttable` structure, dropped-GRB behavior), differing only in indentation style and one caption sentence reworded to state explicitly that the other three bursts are *omitted* (rather than the old wording, which described them as present-but-undetermined — no longer accurate once they're actually dropped from the table). Both scripts' assertions passed on this run. Full rebuild via the new `pdflatex`+`bibtex` sequence (§10): 0 errors, 0 undefined refs, both tables confirmed via `pdftotext` to show exactly 8 GRB080916C rows each with the `threeparttable` footnote marker in place.

## 10. Build tooling — **switched from `latexmk` to `pdflatex`+`bibtex`, 2026-08-21**

At the user's request, after past friction with `latexmk` elsewhere. Full command sequence and the `BIBINPUTS`/`BSTINPUTS` gotcha this surfaced (plain `bibtex out/main` cannot find `ref.bib`/`bibtex/aa.bst`, both of which live at the repo root, not `out/`) are recorded in `HANDOFF.md` §2, not duplicated here — that file is the build-instructions source of truth for the whole paper, not just this topic's tables.
