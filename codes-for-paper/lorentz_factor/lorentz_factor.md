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

**STALE — kept for historical narrative only, do not read these numbers as current.** Confirmed 2026-09-16: the $\Gamma_\text{min}$ values below (134/78/86/129/137/80/83/71) are the pre-BUG-18 numbers. `lorentz_results.csv` moved on twice since: to 507/350/380/447/504/295/287/260 post-BUG-18 (§8.3's original note), then, as of 2026-09-27 (§11), to **507/361/374/741/585/360/375/308** once measured Norris-fit $t_v$ replaced duration for every episode but T90. See §11 for the full before/after table and the decisions behind it. The thermal-$\Gamma$ comparison and $Y$-constraint below were computed from the doubly-stale column and have not been redone with the current values; the qualitative conclusion (bounds too weak to constrain $Y$) almost certainly still holds, but the quoted figures should be recomputed from §11's numbers before being cited again.

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

- **Our bounds are much weaker than published ones for the same burst.** Abdo et al. (2009) obtain $\Gamma_\text{min} = 887 \pm 21$ (their bin b, 3.6–7.7 s) and $608 \pm 15$ (bin d); our TR2, which contains the same ~2 GeV photon, gives 129 pre-BUG-18 / 447 post-BUG-18 / **741 in the current `lorentz_results.csv`, since §11's Norris-measured $t_v$** (see §3's staleness note). The gap to Abdo et al. has narrowed substantially with each correction but is still methodological where it remains, not fully explained by $\Delta T$ alone — Lithwick & Sari's Limit A is a simplified analytic form and Abdo et al. use a fuller treatment. **Our values should be presented as conservative lower limits, not as competitive with published ones.**
- **Only GRB080916C yields values**, since the other three lack redshifts.
- **T90 rows are not independent** of the TR rows that tile them.
- Episodes whose best-fit model has $\beta \geq -1$ are skipped: $\alpha \leq 1$ makes the expression singular.

---

## 6b. A better $t_v$: the Norris-profile fit — **assessed 2026-08-21, feasible; wired into Gamma_min 2026-09-27, see §11**

**Status update (2026-09-27): done — see §11.** Everything below this line describes the assessment and the measurement work as it stood before wiring; §11 records what was actually fed into `Gamma_min` and the resulting numbers. Kept here for the historical reasoning (why it was worth doing, what the caveats were before a real fit existed), not as the current state.

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

A Norris fit would replace that upper bound with a measurement. Since $\Gamma_\text{min} \propto t_v^{-1/(2\alpha+2)}$, the effect is predictable: for TR3, moving $t_v$ from its 40.256 s duration to a Norris-like $\sim$0.7 s would raise $\Gamma_\text{min}$ by a factor $(40.256/0.7)^{0.154} \approx 1.87$ — from the pre-BUG-18 137 to roughly 256 as originally estimated here, or, applied to the then-current post-BUG-18 value (504, §3/§8.3), to roughly $504\times1.87\approx943$. **This illustrative projection is now superseded by an actual measurement, §11**: TR3's real fitted $t_v$ is 15.385 s (not the ~0.7 s guessed here from Bukhari et al.'s GRB110721A), giving $\Gamma_\text{min}=585$ — real, but a smaller jump than this section's illustrative estimate assumed, since GRB080916C's TR3 pulse turned out substantially wider than the guessed proxy value. It would let us state that our limits are conservative *for a quantified reason* rather than by construction.

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
| `lorentz_factor_unknown_z.py` | Limit A redshift sweep for the three no-redshift bursts, plus the shared `run_sweep()` (§13) |
| `lorentz_factor_limit_b_unknown_z.py` | Limit B redshift sweep, own seed; imports `run_sweep()` (§13) |
| `generate_lorentz_table_unknown_z.py` | renders either sweep table from its CSV (`A` / `B` argument) |
| `lorentz_results_unknown_z.csv`, `lorentz_results_limit_b_unknown_z.csv` | one row per no-redshift-burst episode, Gamma columns per swept $z$ |
| `lorentz_table_unknown_z.tex`, `lorentz_table_limit_b_unknown_z.tex` | generated paper tables, Limits A / B at $z=1,3,5,7$ |
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

(Limit A/B values here are post-BUG-18 but pre-§11, i.e. duration-sourced $t_v$ — **superseded by §11's Norris-measured values**, not the pre-fix numbers in §3 above.) For this sample, **Limit A is the larger (binding) bound in every episode**, so $\max(A,B) = A$ throughout and this does not change the paper's reported $\Gamma_\text{min}$ — re-checked with §11's numbers, same conclusion holds. This is a real result, not a wasted computation: it confirms Limit A's dominance rather than assuming it, and Limit B remains available as the fallback for the low-significance-flagged episodes should Limit A's photon association ever be reconsidered.

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

## 11. Norris-measured t_v wired into Gamma_min — **added 2026-09-27**

Full plan in `NORRIS_TV_LORENTZ_INTEGRATION_PLAN.md` (project root) — this section records what
was actually implemented and the resulting numbers; the plan file records the reasoning behind
each decision before implementation and isn't duplicated here.

**What changed.** §4/§6b above describe using each episode's duration as $t_v$ because no
measurement existed. Phase 5 (`PLAN.md`, `codes-for-paper/variability_analysis/`) has since
produced one: a manual, per-burst joint Norris-pulse fit to the 10–400 keV summed light curve,
with per-pulse $t_v$ and its MC-propagated uncertainty in `norris_fit_results_GRB<name>.csv`.
`lorentz_factor.py`/`lorentz_factor_limit_b.py` now use this measurement wherever one exists and
passes a quality gate, falling back to duration otherwise. **This is the fourth time
$\Gamma_\text{min}$ has changed (197→258→134→507→see below)** — every prior change is recorded in
`BUGS.md`/`HANDOFF.md` §4; this one is a systematic improvement (a measurement replacing a
deliberately conservative upper bound), not a bug fix.

### 11.1 Decisions (user, 2026-09-27)

1. **Multi-pulse → episode selection.** Several episodes have more than one Norris-fitted pulse
   whose `t_peak_s` falls inside their window (e.g. `EX0` is a strict superset of `TR1`'s window,
   so `TR1`'s own pulse is also a candidate for `EX0`). The candidate whose `t_peak_s` is closest
   to that episode's own $\Gamma_\text{min}$-defining LAT photon arrival time (`t_arr_s`) is
   selected — the exact rule already locked in `PHASE5_TV_PLAN.md` for the (abandoned) automated
   pipeline, applied here to the manual joint-fit CSVs instead. Implemented in
   `load_norris_tv()`.
2. **GRB080916C's missing `EX0`/`EX1` bounds.** `variability_analysis/shared_utilities.py`'s
   `grb080916C_bounds()` was the only one of the four bursts' bounds functions without `EX0`/`EX1`
   entries (confirmed against `results.json`: `EX0 -0.128_4.864`, `EX1 59.520_67.904`) — so before
   this fix, only `TR1`–`TR5` were eligible for a measured $t_v$ for this burst, the one that
   actually matters. Added, and `fitter_GRB080916C.py` rerun to regenerate
   `norris_fit_results_GRB080916C.csv` (8 rows, was 6) — diffed against the pre-fix CSV first:
   exactly two new rows (`EX0`/pulse 1 relabelled from blank, `EX0`/pulse 2 and `EX1`/pulse 6
   newly added), every pre-existing row's fit parameters unchanged to ~1e-7 relative (ordinary
   optimizer floating-point noise, same class as documented for BUG-23).
3. **Quality gate.** A Norris measurement is only used if the selected pulse's `mc_kept_fraction`
   is at least `MC_KEPT_FRACTION_MIN = 0.5`; below it, the episode falls back to duration with the
   rejection recorded (`load_norris_tv()`'s `rejected_reason`), never silently dropped. Verified
   the gate actually triggers (not a no-op) by temporarily raising the threshold to 0.99 in a
   throwaway check: it correctly flagged exactly the lower-confidence candidates
   (`GRB080916C` EX1/TR2/TR4/TR5, `GRB140206B` TR3, `GRB231129C` EX1/TR2 — all in the
   0.53–0.86 `mc_kept_fraction` range) and none of the high-confidence ones. At the adopted
   threshold (0.5), every one of GRB080916C's episodes passes.
4. **Error propagation.** A measured $t_v$ has its own MC uncertainty
   (`t_v_err_lower_s`/`upper_s`, 16/50/84 percentiles from the Norris fit's own MC), which should
   show up in $\Gamma_\text{min}$'s error bars (`CLAUDE.md`'s "every MC-derived quantity carries an
   error"). The underlying per-draw $t_v$ samples aren't persisted (only the percentiles are), and
   re-deriving them exactly would mean importing the fitted `NorrisFitter` object across folders —
   against `CLAUDE.md`'s "copy rather than fight `sys.path`" convention. **Approximation, stated
   explicitly here rather than left implicit:** `sample_split_normal()` draws $t_v$ per MC
   iteration from a two-piece (split) normal built from the median and its asymmetric
   lower/upper errors, then clips to stay positive (t_v is a divisor in `compute_gamma_min`). This
   is not a reproduction of the Norris fit's own MC draws, and is acceptable only because
   $\Gamma_\text{min}$ is already known to depend weakly on $t_v$ (exponent $\approx0.15$, §5
   above). Duration- and literature-sourced rows are unaffected (a bound/adopted value is treated
   as exact, no resampling).

### 11.2 Results, GRB080916C — supersedes §3/§8.3

Every episode except `T90` (Norris deliberately never fits it — a multi-pulse span, §6b) moves
from duration to a measured $t_v$; every candidate passed the quality gate at the adopted 0.5
threshold.

| Episode | $t_v$ old (duration) [s] | $t_v$ new (Norris) [s] | source | Limit A $\Gamma_\text{min}$: old → new | Limit B $\Gamma_{\text{min},B}$: old → new |
|---|---|---|---|---|---|
| T90 | 62.976 | 62.976 | duration (unchanged) | 507 → 507 | 112 → 112 |
| EX0 | 4.992 | 4.030 | norris | 350 → 361 | 190 → 197 |
| TR1 | 3.584 | 4.030 | norris | 380 → 374 | 213 → 209 |
| TR2 | 10.176 | 0.432 | norris | 447 → **741** | 190 → **352** |
| TR3 | 40.256 | 15.385 | norris | 504 → 585 | 115 → 138 |
| TR4 | 4.224 | 1.192 | norris | 295 → 360 | 154 → 197 |
| EX1 | 8.384 | 1.599 | norris | 287 → 375 | 131 → 182 |
| TR5 | 4.736 | 1.599 | norris | 260 → 308 | 141 → 174 |

`TR2` moves the most (a much shorter measured pulse than its 10.176 s duration), narrowing the gap
to Abdo et al.'s 887 (overlapping interval) further than any prior correction alone. **Limit A
still dominates Limit B in every episode** — re-verified after this change, same conclusion as
§8.3. The shared `tau_hat` bit-identity between the two limits (§8.2) was re-checked and still
holds exactly (max abs diff `0.0`) with the new Norris-derived `delta_T`, as it must for two
functions computing the same closed form from the same inputs.

`gamma_comparison.png`/`.pdf` regenerated from the new CSVs; visually confirmed clean (no gaps,
legend intact, both limits and the thermal $\Gamma$ correctly ordered).

### 11.3 Table changes

`build_latex_table()` in both scripts previously carried a single `\dagger` in the column header,
asserting every shown row was duration-sourced. That assertion would now fail (rows are a mix of
`duration` and `norris`), so the marker moved per-cell (`fmt_t_v()`): `\dagger` for
duration-sourced (unchanged meaning), a new `\ast` for Norris-sourced, printed with its
$1\sigma$ MC uncertainty (`$4.030^{+0.168}_{-0.181}\,\ast$`, not a bare median) per `CLAUDE.md`'s
uncertainty rule. Both tablenotes blocks gained a `$\ast$` entry.

### 11.4 Scope boundary

This work is `GRBResearchWork`-only: `lorentz_results.csv`/`lorentz_results_limit_b.csv`/both
`.tex` tables/`gamma_comparison.png`/`.pdf` are regenerated, and `norris_fit_results_GRB080916C.csv`
now has 8 rows. **`GRBResearchPaper` is untouched** (confirmed via `git status`) — paper
integration (new prose, a `VERIFY.md` entry, copying the regenerated tables) is a deliberate
follow-up once these numbers are reviewed, matching the precedent set by BUG-18 and the RNG
overhaul (§13 of `HANDOFF.md`).

### 11.5 Files

- `codes-for-paper/lorentz_factor/_backup_2026-09-27_pre_norris_tv/` — every file in this folder
  as it stood before this change.
- `codes-for-paper/variability_analysis/_backup_2026-09-27_pre_ex0ex1/` — `shared_utilities.py`
  and `norris_fit_results_GRB080916C.csv` as they stood before the `EX0`/`EX1` bounds fix.
- New in `lorentz_factor.py`: `load_norris_tv()`, `sample_split_normal()`, `fmt_t_v()`,
  `T_V_SOURCE_MARKER`, `VARIABILITY_ANALYSIS_DIR`, `MC_KEPT_FRACTION_MIN`, `NORRIS_TV`.
  `variability_timescale()` now returns a 4-tuple (value, source, err_lower, err_upper).
- `lorentz_factor_limit_b.py` imports the new shared helpers from `lorentz_factor.py` (same
  pattern as the existing `compute_tau_hat` sharing) rather than duplicating them.

## 12. A `>1 GeV` LAT-photon floor on GRB080916C's episode set — **added 2026-09-27**

Follow-up to §11, same session. Raised by the user asking which LAT photons actually anchor each
episode's Norris-measured $t_v$ selection.

### 12.1 The photon inventory that motivated this

Pulling GRB080916C's full LAT photon list (`shared_utilities.py::lat_details()`, `>1 GeV` floor —
matching `fitter_GRB080916C.py`'s own `PHOTON_E_MIN_MEV=1000` convention for photon-to-pulse
assignment) found only **14 photons above 1 GeV in the whole burst**, clustering into exactly two
windows:

| Window (pulse) | Photons `>1 GeV` | Energies [MeV] |
|---|---|---|
| TR2 (pulse 3, $t_\text{peak}=5.85$ s) | 3 | 1693.6, 2110.1, 1500.6 |
| TR3 (pulse 4, $t_\text{peak}=25.91$ s) | 11 | 1106.1, 12421.5, 1230.4, 1382.1, 2568.1, 1649.1, 2499.7, 1742.3, 6721.3, 27428.8, 5707.0 |

No `>1 GeV` photon falls anywhere near `TR1`/`EX0`, `TR4`, or `TR5`/`EX1` — their current
`LAT_analysis/lat_photons.csv` "defining photon" is sub-GeV (301.2/464.5/340.0/989.4 MeV
respectively). One photon (1106.1 MeV, $t=10.215$s) sits inside `TR2`'s *time window* but is
`assign_pulse()`-assigned to pulse 4 (TR3's pulse) rather than pulse 3 — the two classification
schemes (simple time-bounds vs. which pulse's profile actually dominates at that instant) disagree
for this one photon. Doesn't change which photon is TR2's highest-energy one either way (2110.1 >
1106.1), so it doesn't affect anything below.

### 12.2 Decision (user, 2026-09-27)

A sub-GeV highest-energy photon gives too weak a $\gamma\gamma$ pair-production constraint to be
worth reporting as a Limit A/B bound. Episodes without any `>1 GeV` photon are dropped from
`lorentz_factor.py`/`lorentz_factor_limit_b.py`'s output entirely — no CSV row, no table row —
**for GRB080916C only** (the other three bursts' rows are informational, no redshift, no real
$\Gamma_\text{min}$, and are left untouched).

**Conflict found and resolved before implementing:** `gamma_comparison_plot.py` (paper Figure 9)
derives its entire x-axis from `lorentz_results.csv`'s episode set, and GRB080916C's thermal
$\Gamma$ (Pe'er 2007, §8.6) only exists for **T90, EX0, TR1** — the three BB-inclusive episodes.
Dropping `EX0`/`TR1` under this cut would have silently invalidated the paper's existing
comparison ("thermal $\Gamma$ exceeds Limit A by 1.5–2.4$\times$ and Limit B by 4.0–6.7$\times$
across T90/EX0/TR1", already in `section-5-data-analysis.tex`/Figure 9). **Resolved: `EX0`/`TR1`
are exempted from the cut** (`PHOTON_ENERGY_CUT_EXEMPT_EPISODES`) — they keep their existing
sub-GeV-photon Limit A/B bound so that comparison is unaffected. `TR4`/`TR5`/`EX1` have no
BB/thermal-$\Gamma$ counterpart at all (non-BB best-fit models), so nothing downstream depends on
keeping them; they're dropped cleanly.

### 12.3 Result

GRB080916C's Limit A/B tables go from 8 rows to **5**: `T90`, `EX0`, `TR1`, `TR2`, `TR3` survive;
`TR4`, `TR5`, `EX1` are dropped. **Every surviving row's values are bit-identical to before this
change** (max abs diff `0.0` across `t_v_s`/`Gamma_min`/`Gamma_min_B`/both error columns) — this is
a pure row filter, not a recomputation, verified by diffing against the pre-cut CSVs rather than
assumed. The three dropped episodes print a `skipped: E_max=... MeV < 1000 MeV cut` console line
each run, so their absence is never silent.

`gamma_comparison_plot.py` needed **no code changes** — thermal $\Gamma$'s episode set
(T90/EX0/TR1) remains a full subset of the new, smaller Limit A/B episode set, so there's no
`KeyError` risk. Regenerated and visually confirmed: 5 columns instead of 8, T90/EX0/TR1 still
show all three series (thermal, Limit A, Limit B), no gaps or crashes.

### 12.4 Scope

Same boundary as §11: `GRBResearchWork`-only. `GRBResearchPaper` is untouched (confirmed via
`git status`) — the paper's Figure 9/Table 3/prose still reflect the pre-cut 8-episode data until
this is reviewed and integrated deliberately.

---

## 13. Redshift sweep for the three no-redshift bursts — **added 2026-10-01**

Limits A and B were pinned to `REDSHIFTS`, so GRB131014A, GRB140206B and GRB231129C never appeared in any Lorentz table.
Mirroring the Amati unknown-$z$ table, both limits are now evaluated at assumed $z=1,3,5,7$ and given their own tables
(`tab:lorentz_unknown_z`, `tab:lorentz_limit_b_unknown_z`).

### 13.1 Decisions (user, 2026-10-01)

- **Limit B gets its own second table**, as for GRB080916C, rather than extra columns in one wide table.
- **Every LAT episode is kept**, with no photon-energy floor and no TS cut (the §12 floor is GRB080916C-only). Episodes with
  TS < 25 carry a $\ddagger$ in the Limit A table (Limit B does not use the photon, so no flag there). Cutting would have left
  GRB140206B nearly empty.
- **Missing values are shown as "-"**, not an ellipsis. Only GRB140206B TR6 has none: it is a CPL fit, so it has no
  high-energy index and the Lithwick & Sari expression cannot be evaluated.
- **Separate seeds for the two limits**, for consistency with GRB080916C (§8). A first version drew both limits from one
  shared draw set under one seed; the user asked for that to be split.

### 13.2 Implementation

- One draw set per episode, reused across the four redshifts, so the $z$ columns of a row are correlated (same as the Amati
  sweep). $z$ enters only through $d_L$ in $\hat\tau$ and the $(1+z)$ factors; nothing is refit. $t_v$ follows the usual
  precedence (Norris if it passes the quality gate, else duration), and its uncertainty is resampled as in §11.
- `run_sweep(limit, rng, seed)` in `lorentz_factor_unknown_z.py` is shared; each script supplies its own `rng`/`seed`, so the
  two limits use independent draws. Seeds: Limit A 1051076921, Limit B 3962326526.
- Output rows: 18 (5 + 8 + 5), of which 17 have limits.

### 13.3 Results

From $z=1$ to $z=7$ Limit A grows by $\approx3.5$–$3.7\times$ in every episode ($\Gamma_{\min}\approx46$–$341$ to $162$–$1225$);
Limit B by $\approx3.4$–$3.6\times$. Limit A is the larger bound in every episode at every swept $z$. These are illustrations of
the $z$ dependence, not measured limits.

### 13.4 Registries and the Norris fitters

`seed_registry.yaml`, `table_registry.yaml` and `runner_registry.yaml` were updated. The three no-redshift Norris fitters
(`fitter_GRB131014A/140206B/231129C.py`) moved from `unused` to `active` in the seed registry, since their $t_v$ now feeds a
paper table.

### 13.5 Scope

Paper prose and the two `\input` lines are in `section-5-data-analysis.tex`; a LaTeX build was left to the user and has not
been run for these tables. `gamma_comparison` is unchanged (still GRB080916C only); comparing these limits against the
thermal $\Gamma$ at the fiducial $z=2$ is a possible follow-up.
