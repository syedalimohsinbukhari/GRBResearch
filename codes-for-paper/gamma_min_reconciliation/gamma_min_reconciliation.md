# Γ_min literature reconciliation

## What it computes

`reconcile_tv.py` answers: "given our own conservative Γ_min (computed at our episode-duration t_v),
what t_v would a literature source's finer time bin need, to produce *its* higher reported Γ_min?"

Defining equation, from Lithwick & Sari (2001) (also this paper's `eq:gamma_min`, section-5-data-analysis.tex):

$$\Gamma_{\min} \propto t_\text{v}^{-1/(2\alpha+2)}$$

Rearranged and solved for the literature source's implied $t_\text{v}$, holding $\alpha$ fixed at our
own fitted value:

$$t_{\text{v,lit}} = t_{\text{v,ours}} \left(\frac{\Gamma_{\min,\text{lit}}}{\Gamma_{\min,\text{ours}}}\right)^{-(2\alpha+2)}$$

## Why this exists

This project adopts each episode's full duration as $t_\text{v}$ — a deliberately conservative upper
bound, not a measurement of the true variability timescale (see `codes-for-paper/lorentz_factor/lorentz_factor.md`).
A literature source using finer time bins will report a shorter $t_\text{v}$ and therefore a higher
$\Gamma_{\min}$ for the same physical burst — the two are not directly comparable numbers, and the
question a reader (or reviewer) will ask is whether the resulting *difference* is even plausible, or a
sign of an actual inconsistency. First raised for GRB080916C, where our $\Gamma_{\min} = 507$ (Limit A,
T90) sits well below Abdo et al. (2009)'s $\Gamma_{\min} = 887\pm21$ and $608\pm15$ from their own
finer time-resolved analysis of the same burst — `section-6-discussion.tex` originally noted this
without a quantitative check ("we do not attempt to reconcile the two $t_\text{v}$ conventions
quantitatively here"). The user asked for the actual number (2026-09-04), since it's a short
calculation and materially strengthens the argument if it checks out.

## Judgement calls

- **Decided by the user**: do the calculation and add it to the paper as a footnote/aside, replacing
  the earlier disclaimer — but keep a caveat alongside it (see below), not a bare number. Order:
  reconciliation first, then the (reworded) disclaimer.
- **Decided by Claude, and this script's actual design choice**: read $(\Gamma_{\min}, t_\text{v},
  \alpha)$ for "our" side directly from `lorentz_factor/lorentz_results.csv` (the same values already
  in the paper's own Γ_min table) rather than hardcoding them here, so this script can never silently
  drift out of sync with that table if the fits are ever rerun. Only the literature comparison values
  themselves (which aren't derivable from this project's data) are hardcoded, in `RECONCILIATIONS`.

## Validation performed

Cross-checked by hand (a plain Python one-liner, not this script) before the script was written:
$\alpha = 2.251$ (this episode's fitted $\beta = -2.251$), exponent $-1/(2\alpha+2) = -0.1538$,
giving $t_\text{v} \approx 1.66$–$1.67$~s for the $\Gamma_{\min}=887$ case — the script's more precise
$\alpha_\text{LS}=2.2507238\ldots$ (read from the CSV rather than the paper's rounded 2.251) gives
$1.669$~s, matching to the precision quoted in the paper ("$\approx 1.7$~s"). Both agree.

## Known limitations

- **$f_1$ is not reconciled.** The scaling relation used here holds $f_1$ (the 1~MeV photon flux
  entering $\hat\tau$, `eq:tau_hat`) fixed at our own value; in reality $f_1$ also differs between our
  full-duration bin and the literature's finer bin, and that difference is not accounted for. This is
  exactly why the paper text calls this a "scaling check rather than a rigorous reconciliation" — see
  the caveat sentence immediately following the reconciliation in `section-6-discussion.tex`.
- **The second Abdo et al. (2009) value (608) is not discussed in the paper text.** The script computes
  a reconciliation for it too (requires $t_\text{v}\approx19.5$~s, notably different from the 887 case)
  — included for completeness since it was already in the source paper, but the main text only
  reconciles against the *higher* Abdo value (887), since that's the one that sits above this paper's
  thermal $\Gamma$ range and is the actual point of tension being addressed.
