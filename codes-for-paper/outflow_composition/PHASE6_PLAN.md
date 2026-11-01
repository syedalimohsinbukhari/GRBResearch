# Phase 6 — outflow composition (η, σ0): plan, post-literature-review

> ## Update 2026-10-25 — decisions settled, tables generated, paper subsection drafted
>
> Band decisions are settled (stage 2 test band from 10 keV, $L_\gamma$ observer-frame 10 keV–10 GeV; stage 1 $F_\text{BB}$/$F_\text{ob}$ 10 keV–10 GeV with Phase 1/2 unchanged) and stage 1 was rerun at 10 000 samples. LaTeX and markdown tables are generated in `tables/`; the paper has a new section 5.7, an appendix partner table, and matching abstract, §6 and §7 text (`VERIFY.md` entry 2026-10-25, PDF checks pending). This plan is historical; `outflow_composition.md` is the current record.
>
> ## Update 2026-10-24 (later) — stage 2 built too
>
> `zhang_peer_suppression.py` (Zhang & Pe'er minimum-σ bound for the non-BB episodes) is built and validated; `gao_zhang_inversion.py` also runs on the SAFE/MARGINAL `X_BB`
> partners of non-BB episodes, reported separately. User rules: the paper's main table holds BEST-BB episodes only; non-BEST blackbody models are separate and flagged; UNSAFE blackbody models are never used.
> Stage 2 uses the η ≥ η* branch with R₀ = c·t_v/(1+z), as the paper's own numbers imply. Details: `outflow_composition.md` §3a, §6, §8.
>
> ## Update 2026-10-24 — stage 1 built; read `outflow_composition.md` first
>
> The Gao & Zhang inversion (`gao_zhang_inversion.py`) is implemented and validated; its method note is `outflow_composition.md`, which
> supersedes this plan where they differ (Table 2 not Table 1 for regime selection; top-down formulae exist only for δ = 1/3; the Li et al. Table 1
> spot-check is not possible as planned; the pure-fireball $r_0$ is the σ₀ = 0 threshold, not an independent input). User decisions of that date:
> scope extended to all four bursts on the $z=1,2,3,5,7$ grid, both $r_0$ sets, Zhang & Pe'er root-find primary for the non-BB episodes, fallback episodes
> kept and labelled. Stage 2 (Zhang & Pe'er) is not built. The status block below is the 2026-10-20 text.
>
> ## Status as of 2026-10-20
>
> - **Not implemented.** `gao_zhang_inversion.py` and `zhang_peer_suppression.py` do not exist; there is no CSV, LaTeX table or figure for this phase. The folder holds only this plan, `PHASE6_PLAN.pdf` (a render of an earlier revision: it predates the last edit of this file), `CI.md` and `literature/`. The whole folder is **untracked** in git (`??`), and nothing in the paper or registries references it.
> - **The `CI.md` review's requests, and what this plan already answers** (checked against the text below):
>   - Choose between Li et al. (2024)'s closed form and the Zhang & Pe'er root-find *before coding*: answered, in the "Resolved (not adopted as primary)" note under the Li et al. bullets. The Zhang & Pe'er root-find stays primary because no forced-BB upper-limit fits exist for TR2-TR5/EX1.
>   - State the $r_0$ independence assumption (a pure-fireball $r_0$ derived at $\sigma_0=0$ reused as an anchor): answered, in the "Quiet assumption" paragraph; still to be cited or argued, and to carry into the paper and a future "Known limitations" section.
>   - Justify the anchor sitting an order of magnitude below Gao & Zhang's $10^8$ cm: not resolved; recorded as an open modelling choice with no literature support found.
>   - Label the $t_v$ vs duration-fallback non-uniformity: answered, via the required `r0_source`/`variability_source` CSV column. Still to be carried through to the LaTeX table.
>   - Monotonicity check of $L_\text{th}(\sigma)$ and a loud "no solution" failure: answered, in the two "gaps" bullets and Verification #6.
>   - $\delta$ sensitivity: answered, as Verification #5 (a stability check, not a full general-$\delta$ implementation).
>   - State the GRB080916C-only scope in the paper's methods: noted as a "paper-writing reminder" here; must still be done in the `.tex`.
>   - Li et al. Table 1 spot-check: **not done**; this plan marks it `STATUS: not done`.
> - **Decisions still needed from the user before coding** (this plan does not settle them):
>   1. Whether to start Phase 6 at all, given it is optional and the paper's other open items. `PLAN.md` and `HANDOFF.md` both say it is scoped only and ask for confirmation.
>   2. The $r_0$ trial set: the proposed project-specific anchor ($\sim0.9$-$1.4\times10^7$ cm, times 0.1 to 10) versus Gao & Zhang's $10^8$/$10^9$/$10^{10}$ cm, or both. Also whether the independence argument is acceptable.
>   3. Whether the Zhang & Pe'er non-detection criterion (the algebraic root-find, "proposed, open to revision") is acceptable as the primary method, or whether a forced-BB upper-limit refit pass for TR2-TR5/EX1 (which would enable the Li et al. closed form) is worth running.
>   4. The search bracket for $\sigma$ in the root-find, and how to treat episodes on the duration fallback (keep with a label, or exclude).
>   5. Whether to extend beyond GRB080916C using the fiducial-$z$/z-sweep convention. The plan restricts to GRB080916C.
>   6. Whether `PHASE6_PLAN.pdf` should be regenerated, and whether the folder should be committed (it is untracked; `git add` is the user's call).
>   7. The $f_\gamma$ (0.5) and $\delta$ (1/3) fiducials: adopted from the paper; confirm.
>
> Nothing below this header has been edited; it is the 2026-09-27 plan text.

Companion to `PLAN.md` Phase 6 (project root) and `codes-for-paper/lorentz_factor/lorentz_factor.md` §6b/§8. 
Written 2026-09-27, not yet implemented — this is the plan for the outflow-composition work, refined after a literature-review pass found a real problem with Phase 6's originally-stated premise. 
Nothing in this folder exists yet except this plan file.

## Context

`PLAN.md` Phase 6 proposed using Gao & Zhang (2015)'s "top-down" inversion for BB-detected episodes and Zhang & Pe'er (2009)'s suppression-bound for non-BB episodes, to determine GRB080916C's outflow composition (η, σ0). 
A literature-review subagent read both papers (plus Yassine et al. 2017) in full and extracted the actual equations, regime table, and worked
examples — replacing what had previously only been recalled from a prior session's summary.

**This surfaced a real problem with Phase 6's stated premise, not just missing detail.**

## Critical finding: the Γ_min-for-r0 substitution doesn't hold up

`PLAN.md`'s original reasoning was: *"Their own worked example (GRB110721A) supplied the missing
constraint as an assumed r0 — ours doesn't have to, because Limit A/B's Γ_min is already sitting in
lorentz_results.csv/lorentz_results_limit_b.csv and can serve as that second constraint instead."*

The subagent's careful reading of Gao & Zhang's actual equations (38)– (40) shows this doesn't work
as stated, for two concrete reasons:

1. **r0 is a length (jet-base radius), not a Lorentz factor.** It enters every regime formula as a
   fixed multiplicative power (e.g. η ∝ r0,9^{5/24} in Regime II, r0,9^{-2/3} in Regime III/VI).
   $\Gamma_\text{ph}$ is a *derived output* of the equations (also given explicitly in eqs. 38-40),
   not a free input slot — there's no equation where a Lorentz factor gets plugged in directly in
   r0's place.
2. **Equality vs. inequality mismatch.** Even inverting $\Gamma_\text{ph} (r_0) = \Gamma_\text{min}$
   to back out an implied r0 would mean treating Lithwick & Sari's *lower bound*
   ($\Gamma \geq \Gamma_\text{min}$) as if it were the *exact* $\Gamma_\text{ph}$ — an assumption
   Gao & Zhang's own math doesn't sanction, with no guarantee of which direction it biases (η, σ0), and that direction
   would need working out per regime.

**Gao & Zhang's own closure choice is instead: assume a fixed r0 outright** (their GRB110721A
application tries three trial values — $10^8$, $10^9$, $10^{10}$ cm — and reports η (t)/σ0 (t) for
each, purely graphically, no single quoted number). This is much closer to what Pe'er (2007)
itself already does (assume r0, or rather solve for it under the pure-fireball assumption) than to
"use Γ_min instead of guessing."

**Also found**: no numeric worked example exists in Gao & Zhang (2015) analogous to Pe'er (2007)'s
appendix (`pe_er_photosphere.py`'s `validate_against_peer2007()` has no equivalent anchor to check
against) — their own case study is presented only as digitized-figure curves, no quoted (kT_BB, F_BB, z) → (η, σ0)
tuple. And Zhang & Pe'er (2009)'s suppression bound is a graphical
per-epoch comparison (predicted suppressed spectrum vs. actual fitted spectrum), not a single
closed-form σ ≥ f (...) inequality — implementing it means defining our own concrete algebraic
non-detection criterion, not porting one equation.

## Verified equations and tables (full detail - see the literature-review subagent's report)

reproduced in full further below in this file's "Appendix" section, or re-run the same review)

- **Regime table is Table 1** (p.34, δ=1/3 fixed), not "Table 3" as `PLAN.md` guessed — Table 3 is the general-δ version; Table 2 is the observational-quantities version.
- Six regimes, keyed on $\eta \gtrless (1+\sigma_0)^{1/2}$ and where $r_\text{ph}$ falls relative to $r_\text{ra}$/$r_\text{c}$. Only Regimes II, III, V, VI close (I and IV are degenerate — same degeneracy as Pe'er 2007's own $r_\text{ph} \leq r_\text{c}$ case); Regime VI shares Regime III's scalings.
- Full closed-form (η, σ0, r_ph, Γ_ph, σ_ph, σ_r15) formulas for Regimes II (eq. 38), III/VI (eq. 39), V (eq. 40) — all in terms of measurables ($kT_\text{ob}$, $F_\text{BB}$, $z$, $d_L$) plus assumed constants ($r_0$, $f_\gamma$=0.5 fiducial).
- δ (slow-acceleration index, \Gamma\propto r^\delta for r_\text{ra} < r < r_\text{c}): fixed at **1/3** in every main-text result (the paper's own fiducial/canonical choice, "achievable for an impulsive, non-dissipative magnetic shell") — adopting 1/3 is following the paper's own convention, not inventing a new assumption.
- Zhang & Pe'er (2009): pure-fireball $R_\text{ph}$/$T_\text{ph}^\text{ob}$/$L_\text{th}$
  predictions (unnumbered but cited from Mészáros & Rees 2000 / Pe'er 2008), eq. (3) for
  $L_\text{th}$, and the suppression rule ($L_w \to L_w/ (1+\sigma)$ substituted throughout, then
  increase $\sigma$ until the suppressed prediction stops exceeding the observed Band spectrum) —
  needs $L_w=L_\gamma$ (100% radiative efficiency assumption) and $R_0=c\,\delta t$ (assumed
  variability timescale) as closure inputs, same structural role as Gao & Zhang's $r_0$.
- Yassine et al. (2017) confirmed as direct precedent for per-episode (not burst-averaged)
  presentation: a table + a plot with one color/marker per time interval, exactly matching this
  project's existing `GRBPlotStyle` per-episode convention.

## Decision (confirmed with user, 2026-09-27)

**Assume a fixed r0 outright, following Gao & Zhang's own precedent** — try trial values and
report η (t)/σ0 (t) for each, the same way their own GRB110721A case study does. Abandons the
original "skip guessing r0 because we have Γ_min" motivation, since that didn't hold up anyway;
this is the literature-faithful choice.

## Follow-up literature check (added 2026-09-27, post-decision)

A second research pass, done specifically to find independent numeric anchors for validation (see
"Validation" below), turned up two usable cross-checks. Both papers were read full-text, not
abstract-only.

- **Zhang & Pe'er (2009) itself, on GRB080916C, quotes a number in the text** (not just Fig. 2):
  assuming $R_0=c\,\delta t_\text{ob}$ with $\delta t_\text{ob}=0.5$ s and $L_w=L_\gamma$, their
  epoch "b" needs **$\sigma_\text{min}\approx20$** to suppress the thermal prediction below the
  observed spectrum ($\approx15$–20 under two more conservative trial temperatures). Episode
  binning differs from ours (their epochs a–e vs. our TR2/TR3/TR4/TR5/EX1) and $R_0$ differs from
  this project's Norris-$t_v$ choice, so this is an **order-of-magnitude sanity target**, not an
  exact reproduction target, for `zhang_peer_suppression.py`'s $\sigma_\text{min}$ root-find.
- **Li, Gao, Lan & Zhang (2024, arXiv:2408.01161, ApJ, accepted)** — a systematic Fermi-sample
  application of Gao & Zhang's own method:
    - Their eq. (8) + Table 2 independently restate the Regime II/III/V/VI closed forms with explicit
      numeric prefactors (14.8, 0.24, $1.43\times10^{-5}$ for II; 8.28, $9.42\times10^{-2}$ for III;
      41.4, 5.28, $1.16\times10^{-5}$ for V; 8.28, $5.63\times10^{-3}$ for VI) — **a genuine
      second-source check** on our transcription of eqs. 38–40, to diff against once
      `gao_zhang_inversion.py` is implemented.
    - Their Table 1 tabulates $(kT_\text{ob}, f_\text{th})\to (1+\sigma_0)$ at
      $r_0=10^8/10^9/10^{10}$ cm for many time slices, including **GRB110721A** (Gao & Zhang's own
      case-study burst — e.g. $kT=34.419$ keV, $f_\text{th}=0.0854 \to 1+\sigma_0=4.784$ at
      $r_0=10^8$, $22.206$ at $r_0=10^9$, no solution at $r_0=10^{10}$) and **GRB080916C** itself. Not
      a blind plug-in — $F_\text{ob}$ and $z$ per row aren't given in the table, only $f_\text{th}$ —
      but the regime-boundary values and prefactors are directly checkable now.
    - They apply Gao & Zhang's own closed-form to *non-thermal* slices too, via a simulated
      upper-limit thermal component, rather than inventing a new root-find — a **published,
      closed-form precedent that could replace this plan's "proposed, open to revision" suppression
      criterion** for TR2/TR3/TR4/TR5/EX1.
      **Resolved (not adopted as primary — checked, not assumed)**: Li et al.'s method needs a
      per-episode upper limit on a BB normalization/flux, from a spectral fit that *includes* a BB
      component even when it isn't preferred. Checked `bb_fraction/bb_flux_fraction.csv` — it only
      has rows for T90/EX0/TR1 (the 3 BB-detected episodes); no such fit, upper-limit or otherwise,
      exists anywhere in this project's pipeline for TR2/TR3/TR4/TR5/EX1. Adopting Li et al.'s
      approach would mean running a brand-new spectral-fitting pass (force-fit a BB component,
      extract its normalization upper bound) for those 5 episodes — a new, independently-risked
      derivation, not a reuse of already-computed values, which is exactly what `CLAUDE.md`'s
      "verify, don't assume"/reuse convention (and Verification item 4 below) says to avoid this
      late in the project. **Zhang & Pe'er's root-find stays the primary method for the non-BB
      episodes**; Li et al.'s closed-form remains a validation cross-check (via their in-text
      GRB080916C number, see above) rather than the implementation path. Revisit if a future phase
      adds forced-BB upper-limit fits for these episodes.
    - Leads flagged but not verified: GRB100724B (Guiriec et al. 2011) and GRB130323A (Guiriec et al.
        2013) as further independent applications of one or both methods.
- PDF saved: `literature/li_gao_lan_zhang_2024_arxiv2408.01161.pdf` (this folder). Zhang & Pe'er (2009) and Gao & Zhang
  (2015) are already present in
  `../../../GRBResearchPaper/Literature Review/` (`Zhang2009 - ...pdf`, `GaoHeZhang2015.pdf`) — no
  need to re-save those.

## Scope

**GRB080916C only**, matching every other decision this session (only burst with a redshift):

- **Gao & Zhang inversion** → the 3 BB-inclusive episodes: `T90`, `EX0`, `TR1` (their
  $kT_\text{BB}$/$F_\text{BB}$ already computed in `photospheric_radius/pe_er_photosphere.csv`;
  their best-fit models are `SBPL_BB`/`BAND_BB`/`BAND_BB` respectively).
- **Zhang & Pe'er suppression bound** → the non-BB episodes: `TR2`, `TR3`, `TR4`, `TR5`, `EX1`.
  This method needs no LAT photon coverage at all (unlike Limit A/B), so it is **not** subject to
  the `>1 GeV` floor added earlier this session — that cut was specific to the opacity pipeline,
  irrelevant here. All 5 non-BB episodes are in scope regardless of their LAT photon energy.
- **Paper-writing reminder**: "GRB080916C only because it's the only burst with a redshift" needs
  to be stated plainly in the paper's methods text, not left implicit — otherwise a reader will
  wonder why outflow composition doesn't appear for the other three GRBs in the sample.

## Implementation plan

### New module: `codes-for-paper/outflow_composition/`

Following the established template (`bb_fraction/`, `photospheric_radius/`, `lorentz_factor/`):
one script per method (kept separate, like Limit A/Limit B, since they apply to disjoint episode
sets and share no formula), a CSV + LaTeX table + plot each, and this file playing the role of the
required per-folder method note (`CLAUDE.md`'s "Per-folder method notes" convention) once it's
extended with the "what was actually built and validated" sections as implementation proceeds.

#### `gao_zhang_inversion.py` (BB episodes: T90, EX0, TR1)

- Inputs per episode, all already computed: $kT_\text{BB}$, $F_\text{BB}$ (from
  `pe_er_photosphere.csv`), $F_\text{ob}$/$f_\text{th}=F_\text{BB}/F_\text{ob}$ (from the fitted
  model), $z=4.35$, $d_L (z)$ (same cosmology as everywhere else in this project — $H_0=69.6$,
  $\Omega_m=0.286$).
- Fixed constants, adopted as the paper's own fiducial choices (not new assumptions this project
  is introducing): $\delta=1/3$, $f_\gamma=0.5$.
- **r0 trial values — proposed, open to revision at review**: reuse this project's own
  already-computed pure-fireball $r_0$ range for these exact episodes (`pe_er_photosphere.csv`,
  GRB080916C: $\sim0.9$–$1.4\times10^7$ cm) as the anchor, plus a
  bracketing sweep (e.g. $\times0.1$ to $\times10$ around that range) to show sensitivity — more
  physically grounded for *this* burst than blindly reusing Gao & Zhang's own $10^{8-10}$ cm trial
  values, which were chosen for a different burst (GRB110721A) at a different physical scale.
  **Quiet assumption, now stated explicitly**: that pure-fireball $r_0$ was derived *under
  $\sigma_0=0$*, and is now being reused to anchor a calculation whose point is to test whether
  $\sigma_0>0$. This isn't circular only if $r_0$ is treated as a central-engine/launch-radius
  scale set by the progenitor geometry — roughly independent of how the launched energy later
  splits between kinetic (baryonic) and Poynting-flux (magnetic) components — rather than as a
  quantity the fireball assumption itself determined. That independence assumption is physically
  standard but **not independently cited here** — flag as a known limitation, not a settled point,
  until a citation or argument for it is added. Separately: this anchor sits roughly an order of
  magnitude below Gao & Zhang's smallest GRB110721A trial value ($10^8$ cm); no independent
  citation was found to justify that departure (checked while researching the "Follow-up literature
  check" cross-checks above, without success) — state as an open modeling choice, not a
  literature-supported one, in this file's "Known limitations" section once implementation starts.
- Implement regime selection via **Table 1** (p.34, $\delta=1/3$ fixed): classify each (episode, trial-$r_0$) pair into
  Regime II, III/VI, or V using the boundary inequalities on
  $\eta^{a} (1+\sigma_0)^{b}$ vs. $7.22\times10^5 L_{w,52}^{3/5} r_{0,9}^{-3/5}$; **Regimes I/IV are
  degenerate** (per the paper itself) — detect and report as "cannot be inferred" rather than
  silently defaulting to another regime or crashing.
- Implement eqs. (38)/ (39)/ (40) for $(1+\sigma_0, \eta, r_\text{ph}, \Gamma_\text{ph})$
  (+ $\sigma_\text{ph}$, $\sigma_{r15}$ where given) per regime, exactly as extracted.
- MC uncertainty: propagate the same $kT_\text{BB}$/$F_\text{BB}$ draws already used in
  `pe_er_photosphere.py` (reuse its resampling, don't re-derive), holding $r_0$ fixed per trial
  value (the trial sweep is the $r_0$ sensitivity study; MC is for the *measured* inputs).

#### `zhang_peer_suppression.py` (non-BB episodes: TR2, TR3, TR4, TR5, EX1)

- Inputs: each episode's fitted total flux/model (for the "actual observed spectrum" comparison),
  $z=4.35$, $d_L (z)$.
- Assumed constants (same structural role as Gao & Zhang's $r_0$): $L_w=L_\gamma$ (100% radiative
  efficiency — the paper's own assumption, carried over, itself a lower-limit-type premise) and
  $R_0=c\,\delta t$ using each episode's own variability timescale — **reuse the Norris-measured
  $t_v$ from this session's Limit A/B work** (`lorentz_factor.py`'s `load_norris_tv()` output) as
  $\delta t$ where available, falling back to duration otherwise, for exactly the same reason it
  was adopted there. **This fallback makes episodes non-uniform in a way that must be labeled, not
  just implemented**: $\sigma_\text{min}$ depends on $R_0$, so a $t_v$-branch episode and a
  duration-fallback episode are answering slightly different physical questions, not two directly
  comparable numbers. Add an explicit `r0_source` (or `variability_source`) column to the output
  CSV (`norris_tv` / `duration_fallback`) so the eventual LaTeX table doesn't present all 5
  $\sigma_\text{min}$ values as methodologically uniform when some aren't.
- **Concrete non-detection criterion — proposed, open to revision at review, since the paper gives
  none as a closed formula**: solve for the minimum $\sigma$ such that the suppressed predicted
  $L_\text{th} (\sigma)$ (via $L_w\to L_w/ (1+\sigma)$ substituted through the pure-fireball
  $R_\text{ph}$/$T_\text{ph}^\text{ob}$ formulas) no longer exceeds that episode's actual fitted
  total flux at the predicted $T_\text{ph}^\text{ob} (\sigma)$ — an algebraic translation of Zhang &
  Pe'er's own graphical criterion (their Fig. 2: suppressed curve vs. observed Band spectrum),
  solved numerically (root-find on $\sigma$) rather than read off a plot. See "Method choice
  resolved" above for why this stays the primary path over Li et al.'s closed-form alternative.
- **Two gaps in the root-find spec itself, to close before/during implementation, not after**:
  1. No stated check that $L_\text{th}(\sigma)$ is monotonic in $\sigma$ over the search range
     before root-finding on it. Verify monotonicity numerically first (e.g. evaluate on a coarse
     grid); if it isn't monotonic somewhere physical, a naive root-finder can converge to the wrong
     root or fail silently — don't assume monotonicity without checking.
  2. The $\sigma_\text{min}\to0$ edge case is already covered (Verification #3), but not the
     opposite failure: suppression never brings the prediction below the observed flux even as
     $\sigma\to\infty$ within the search bracket. That must fail loudly and report "no solution",
     analogous to how Regime I/IV degeneracy is already handled for the Gao & Zhang side — not
     silently return the bracket's edge value as if it were a real $\sigma_\text{min}$.
- Report $\sigma_\text{min}$ per episode as a **lower limit**, explicitly labeled as such (matches
  this project's existing convention for one-sided bounds, e.g. $\Gamma_\text{min}$'s own framing).

### Validation

- **Gao & Zhang (2015) itself gives no single quoted (kT_BB, F_BB, z) → (η, σ0) tuple** (their
  GRB110721A case study is figure-only) — still true, still flagged honestly.
- **But two independent anchors now exist (see "Follow-up literature check" above), so validation
  is not purely self-consistency**:
    1. Regime-formula prefactors from Li, Gao, Lan & Zhang (2024) eq. (8)/Table 2 — diff our
       transcription of eqs. 38–40 against theirs directly.
    2. Zhang & Pe'er (2009)'s own in-text $\sigma_\text{min}\approx15$–20 for GRB080916C — an
       order-of-magnitude sanity target for `zhang_peer_suppression.py`'s output (not exact, since
       episode binning and $R_0$ assumption differ).
    3. Li et al.'s Table 1 GRB110721A/GRB080916C rows, once $F_\text{ob}/z$ are pulled from their
       methods section — a closer numeric spot-check, still to be completed before implementation
       is considered validated. **STATUS: not done.** This is the one validation step that is
       currently just a prose plan, not a completed check — call it out again explicitly once
       implementation starts so it doesn't quietly slip off the list.
- **Still required regardless of the above**: confirm Regime II and III/VI formulas agree at the
  $r_\text{ph}=r_\text{c}$ boundary (they must, by construction — a real bit-identity check, same
  spirit as this project's existing `tau_hat` cross-check between Limit A/B). This is a check on
  our own algebra, not a substitute for the external anchors above.
- State whatever residual gap remains (e.g. if the Table 1 spot-check can't be completed for lack
  of $F_\text{ob}/z$) explicitly in this file's future "Known limitations" section, the same way
  `lorentz_factor.md` already documents $\Gamma_\text{min}$'s own systematic-uncertainty budget —
  a number nobody independently checked against a published value should be labelled as such, per
  `CLAUDE.md`.

### Output conventions

Same as every other module: `update_style()`, `GRBPlotStyle`, CSV columns carrying
`grb_name`/`episode`/`z`/`H0`/`Om0`/`n_samples`/`seed` plus the derived quantities with asymmetric
MC errors, a generated LaTeX table via the established `% AUTO-GENERATED` convention, and a plot
following Yassine et al. (2017)'s precedent (confirmed by the literature review): one color/marker
per episode, quantity vs. episode (matching this project's existing per-episode plotting style,
e.g. `lorentz_factor/gamma_comparison_plot.py`), not a burst-averaged single point.

## Verification

1. Regime-boundary continuity check (Regime II ↔ III/VI at $r_\text{ph}=r_\text{c}$) — must agree
   to numerical precision, not approximately.
2. Sanity-check the $r_0$ sensitivity sweep: confirm $\eta$/$\sigma_0$ move in the directions the
   eq. (38)- (40) exponents predict as $r_0$ is varied (a numeric check, not a read-through).
3. Zhang & Pe'er suppression bound: confirm the root-find converges and that $\sigma_\text{min}\to0$
   in the limit where the unsuppressed prediction already sits below the observed flux (a sanity
   edge case). Also check the resulting $\sigma_\text{min}$ values land in the same
   order-of-magnitude as Zhang & Pe'er's own in-text $\approx15$–20 for this burst (see "Follow-up
   literature check") — a large discrepancy would indicate an error, not just a binning difference.
4. Confirm every number traces back to already-computed CSVs (`pe_er_photosphere.csv`, Norris
   $t_v$) rather than re-deriving $kT_\text{BB}$/$F_\text{BB}$/$t_v$ a second, independently-risked
   way — per `CLAUDE.md`'s "Verify, don't assume" and the drift risk it warns about.
5. $\delta$ sensitivity: $\delta=1/3$ is Gao & Zhang's own fiducial choice (not new), but unlike
   $r_0$ it currently gets no sensitivity sweep, even though Table 3 (general-$\delta$) exists to
   check against. At minimum, confirm numerically that regime classification is stable under
   $\delta=1/3\pm$ a reasonable spread before treating $\delta=1/3$ as settled — a one-paragraph
   "checked, doesn't change the classification" result, not a full general-$\delta$
   re-implementation.
6. `zhang_peer_suppression.py` root-find: confirm $L_\text{th}(\sigma)$ is monotonic over the
   search range before relying on the root-finder, and confirm the "no solution within the search
   bracket" failure mode reports explicitly (not silently) — see the two root-find gaps under
   "Implementation plan" above.
7. No commits — leave the working tree for review, same as every other change this session.