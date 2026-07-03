"""
Created on Mar 17 20:57:08 2026

GRB Spectral Model Comparison
------------------------------
Statistical tests for comparing nested spectral models using your existing
Model dataclass. Requires only cstat, dof, parameters, and covariance_matrix.

Tests implemented
-----------------
  - Likelihood Ratio Test (LRT) via Δcstat ~ χ²(Δk)
  - AIC, AICc, BIC and their deltas
  - AIC model weights
  - BIC Bayes factor approximation
  - BB normalization significance (z-score)
  - Parameter stability (shift in units of combined sigma)
  - Covariance/correlation matrix analysis (condition number, BB correlations)

Usage
-----
    result = ModelComparison(simple_model, complex_model, bb_param_names=["kT", "norm_bb"])
    result.report()
    result.report(save_path="model_comparison_report.txt")
"""

import sys
from contextlib import redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import stats

# -----------------------------------------------------------------------------
# Thresholds (edit to match your pipeline conventions)
# -----------------------------------------------------------------------------

DELTA_CSTAT_THRESHOLD = 28.74  # pipeline detection threshold
BB_SIGMA_THRESHOLD = 3.0  # minimum sigma for BB normalization
HIGH_CORR_THRESHOLD = 0.9  # |ρ| above this flags degeneracy


class _ReportWriter:
    """Write report output to multiple text streams."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


# -----------------------------------------------------------------------------
# Result containers
# -----------------------------------------------------------------------------


@dataclass
class LRTResult:
    delta_cstat: float
    delta_k: int
    p_value: float
    sigma: float
    detected: bool  # True if delta_cstat > threshold


@dataclass
class ICResult:
    aic_simple: float
    aic_complex: float
    aic_c_simple: float
    aic_c_complex: float
    bic_simple: float
    bic_complex: float
    delta_aic: float  # positive → complex is better
    delta_aic_c: float
    delta_bic: float
    w_simple: float  # AIC model weight
    w_complex: float
    bayes_factor: float  # exp(ΔBIC/2) approximation


@dataclass
class EvidenceSummary:
    lrt_supports_complex: bool
    aic_supports_complex: bool
    bic_supports_complex: bool
    aic_evidence_ratio: float
    bic_evidence_ratio: float
    selection_conflict: bool
    recommended_model: str
    recommendation_reason: str


@dataclass
class BBSignificanceResult:
    param_name: str
    value: float
    error: float
    z_score: float
    significant: bool  # z > BB_SIGMA_THRESHOLD


@dataclass
class ParamShift:
    name: str
    v_simple: float
    v_complex: float
    sigma: float  # shift in units of combined error
    flag: bool  # True if shift > 3σ


@dataclass
class CovarianceAnalysis:
    condition_number: float  # raw — sensitive to parameter scaling
    condition_number_scaled: float  # from correlation matrix — scaling-independent
    is_ill_conditioned: bool  # now based on the scaled version
    correlation_matrix: np.ndarray
    param_names: list
    flagged_pairs: list  # pairs with |ρ| > HIGH_CORR_THRESHOLD
    bb_correlations: Optional[dict]  # {(bb_param, continuum_param): ρ}


# -----------------------------------------------------------------------------
# Main class
# -----------------------------------------------------------------------------


def scaled_condition_number(cov):
    """
    Normalize by parameter scales before computing condition number.
    Equivalent to computing the condition number of the correlation matrix.
    """
    std = np.sqrt(np.diag(cov))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(std, std)
    corr = np.nan_to_num(corr, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    try:
        cond = np.linalg.cond(corr)
    except np.linalg.LinAlgError:
        cond = np.inf
    return corr, cond


class ModelComparison:
    """
    Compare a simple (nested) model against its complex extension.

    Parameters
    ----------
    simple : Model
        The simpler model (fewer parameters), e.g., Band.
    complex_ : Model
        The model with extra components, e.g., Band+BB.
    bb_param_names : list[str], optional
        Parameter names belonging to the extra BB component in complex_.
        Matching is case-insensitive and partial, e.g. "norm_bb" will match
        a parameter named "norm_BB_1".
    """

    def __init__(self, simple, complex_, bb_param_names: Optional[list] = None, report_unsafe: bool = False):
        self.simple = simple
        self.complex = complex_
        self.bb_names = [s.lower() for s in (bb_param_names or [])]
        self._warnings = [
            "LRT assumes the compared models are nested.",
            "BB normalization near a physical boundary can make the chi-square LRT approximation optimistic.",
            "AIC/BIC measure model-selection evidence, not direct physical detection.",
        ]

        self._k_simple = len(simple.parameters)
        self._k_complex = len(complex_.parameters)
        self._delta_k = self._k_complex - self._k_simple

        if self._delta_k <= 0:
            raise ValueError(
                f"complex_ must have more parameters than simple. "
                f"Got Δk = {self._delta_k}."
            )
        if self.complex.is_unsafe and not report_unsafe:
            raise ValueError(
                f"Complex model {complex_.name} is not good. " f"Cannot compare to it."
            )

        self._delta_cstat = simple.cstat - complex_.cstat
        if self._delta_cstat < 0:
            self._warnings.append(
                "Complex model has higher cstat than the simple model despite having more parameters."
            )
        if not self.bb_names:
            self._warnings.append(
                "No BB parameter name was provided; BB significance diagnostics were skipped."
            )

        # run all tests on construction
        self._lrt = self._run_lrt()
        self._ic = self._run_ic()
        self._bb_sig = self._run_bb_significance()
        self._shifts = self._run_param_stability()
        self._cov = self._run_covariance_analysis()
        self._cov2 = self._run_covariance_analysis(False)
        self._evidence = self._run_evidence_summary()

    # -- public properties ----------------------------------------------------

    @property
    def lrt(self) -> LRTResult:
        return self._lrt

    @property
    def information_criteria(self) -> ICResult:
        return self._ic

    @property
    def bb_significance(self) -> list:
        return self._bb_sig

    @property
    def parameter_stability(self) -> list:
        return self._shifts

    @property
    def covariance_analysis(self) -> Optional[CovarianceAnalysis]:
        return self._cov

    @property
    def simple_covariance_analysis(self) -> Optional[CovarianceAnalysis]:
        return self._cov2

    @property
    def evidence_summary(self) -> EvidenceSummary:
        return self._evidence

    @property
    def assumption_warnings(self) -> list[str]:
        return self._warnings.copy()

    # -- internal helpers -----------------------------------------------------

    def _n(self, model) -> int:
        """Effective sample size for a model."""
        return model.dof + len(model.parameters)

    def _aic(self, model) -> float:
        return model.cstat + 2 * len(model.parameters)

    def _aicc(self, model, model_label: str) -> float:
        k = len(model.parameters)
        n = self._n(model)
        denominator = n - k - 1
        if denominator <= 0:
            self._warnings.append(
                f"AICc is undefined for {model_label} model because n - k - 1 <= 0; reported as NaN."
            )
            return np.nan
        return self._aic(model) + (2 * k * (k + 1)) / (n - k - 1)

    def _bic(self, model) -> float:
        k = len(model.parameters)
        n = self._n(model)
        return model.cstat + k * np.log(n)

    def _is_bb_param(self, name: str) -> bool:
        if not self.bb_names:
            return False
        return any(b in name.lower() for b in self.bb_names)

    # -- test implementations -------------------------------------------------

    def _run_lrt(self) -> LRTResult:
        dc = self._delta_cstat
        dk = self._delta_k
        p = 1.0 - stats.chi2.cdf(dc, df=dk)
        p = max(p, 1e-300)
        sigma = abs(stats.norm.ppf(p / 2.0))
        return LRTResult(
            delta_cstat=dc,
            delta_k=dk,
            p_value=p,
            sigma=sigma,
            detected=dc > DELTA_CSTAT_THRESHOLD,
        )

    def _run_ic(self) -> ICResult:
        aic_s = self._aic(self.simple)
        aic_c = self._aic(self.complex)
        aic_c_s = self._aicc(self.simple, "simple")
        aicc_c = self._aicc(self.complex, "complex")
        bic_s = self._bic(self.simple)
        bic_c = self._bic(self.complex)

        d_aic = aic_s - aic_c  # positive → complex wins
        d_aic_c = aic_c_s - aicc_c
        d_bic = bic_s - bic_c

        # AIC model weights
        min_aic = min(aic_s, aic_c)
        e_s = np.exp(-0.5 * (aic_s - min_aic))
        e_c = np.exp(-0.5 * (aic_c - min_aic))
        w_s = e_s / (e_s + e_c)
        w_c = e_c / (e_s + e_c)

        # BIC evidence ratio / Bayes-factor approximation (Kass & Raftery)
        bf = self._safe_exp_half(d_bic)  # > 1 favours complex

        return ICResult(
            aic_simple=aic_s,
            aic_complex=aic_c,
            aic_c_simple=aic_c_s,
            aic_c_complex=aicc_c,
            bic_simple=bic_s,
            bic_complex=bic_c,
            delta_aic=d_aic,
            delta_aic_c=d_aic_c,
            delta_bic=d_bic,
            w_simple=w_s,
            w_complex=w_c,
            bayes_factor=bf,
        )

    @staticmethod
    def _safe_exp_half(delta: float) -> float:
        with np.errstate(over="ignore", invalid="ignore"):
            return float(np.exp(delta / 2.0))

    def _run_evidence_summary(self) -> EvidenceSummary:
        lrt_supports_complex = self._lrt.detected
        aic_supports_complex = self._ic.delta_aic > 0
        bic_supports_complex = self._ic.delta_bic > 0
        support_flags = [
            lrt_supports_complex,
            aic_supports_complex,
            bic_supports_complex,
        ]
        selection_conflict = any(support_flags) and not all(support_flags)

        complex_recommended = lrt_supports_complex and (
            aic_supports_complex or bic_supports_complex
        )
        recommended_model = (
            self.complex.name if complex_recommended else self.simple.name
        )
        if complex_recommended:
            reason = "LRT passes the detection threshold and AIC or BIC supports the complex model."
        elif not lrt_supports_complex:
            reason = "LRT does not pass the detection threshold."
        else:
            reason = "LRT passes, but neither AIC nor BIC supports the complex model."

        return EvidenceSummary(
            lrt_supports_complex=lrt_supports_complex,
            aic_supports_complex=aic_supports_complex,
            bic_supports_complex=bic_supports_complex,
            aic_evidence_ratio=self._safe_exp_half(self._ic.delta_aic),
            bic_evidence_ratio=self._safe_exp_half(self._ic.delta_bic),
            selection_conflict=selection_conflict,
            recommended_model=recommended_model,
            recommendation_reason=reason,
        )

    def _run_bb_significance(self) -> list:
        results = []
        for p in self.complex.parameters:
            if self._is_bb_param(p.name):
                z = (p.value / p.error) if p.error else np.inf
                results.append(
                    BBSignificanceResult(
                        param_name=p.name,
                        value=p.value,
                        error=p.error,
                        z_score=z,
                        significant=abs(z) >= BB_SIGMA_THRESHOLD,
                    )
                )
        if self.bb_names and not results:
            self._warnings.append(
                "No complex-model parameter matched the provided BB parameter name."
            )
        if len(results) > 1:
            self._warnings.append(
                "More than one BB-like parameter matched; BB diagnostics assume a single BB component."
            )
        return results

    def _run_param_stability(self) -> list:
        """
        For every parameter shared between both models (matched by name,
        case-insensitive), compute the shift in units of combined sigma.
        """
        simple_map = {p.name.lower(): p for p in self.simple.parameters}
        complex_map = {p.name.lower(): p for p in self.complex.parameters}

        shifts = []
        for name_l, ps in simple_map.items():
            if name_l in complex_map:
                pc = complex_map[name_l]
                denom = np.sqrt(ps.error**2 + pc.error**2)

                # both fixed
                if denom == 0:
                    if ps.value == pc.value:
                        continue  # same fixed value — no shift, skip
                    else:
                        shift = np.inf  # fixed to different values — flag it
                        shifts.append(
                            ParamShift(
                                name=ps.name,
                                v_simple=ps.value,
                                v_complex=pc.value,
                                sigma=shift,
                                flag=True,  # always flag — physically suspicious
                            )
                        )
                        continue

                shift = abs(ps.value - pc.value) / denom
                shifts.append(
                    ParamShift(
                        name=ps.name,
                        v_simple=ps.value,
                        v_complex=pc.value,
                        sigma=shift,
                        flag=shift > 3.0,
                    )
                )

        return shifts

    def _run_covariance_analysis(
        self, analyze_complex=True
    ) -> Optional[CovarianceAnalysis]:
        """Analyze the covariance matrix of the complex model."""
        cov_ = self.complex if analyze_complex else self.simple
        model_label = "complex" if analyze_complex else "simple"
        if cov_.covariance_matrix is None:
            self._warnings.append(
                f"Covariance diagnostics skipped for {model_label} model: covariance is missing."
            )
            return None

        cov = cov_.covariance_matrix_value  # symmetric ndarray
        names = [p.name for p in cov_.parameters]
        errors = [p.error for p in cov_.parameters]

        # ── exclude fixed parameters (σ = 0) ────────────────────────────
        free_idx = [i for i, e in enumerate(errors) if e > 0]

        if len(free_idx) < 2:
            self._warnings.append(
                f"Covariance diagnostics skipped for {model_label} model: fewer than two free parameters."
            )
            return None  # nothing meaningful to analyze

        cov = cov[np.ix_(free_idx, free_idx)]
        names = [names[i] for i in free_idx]

        n = len(names)

        if cov.shape != (n, n):
            self._warnings.append(
                f"Covariance diagnostics skipped for {model_label} model: covariance shape does not match parameters."
            )
            return None

        # condition number — high value signals near-degeneracy
        try:
            cond_raw = np.linalg.cond(cov)
        except np.linalg.LinAlgError:
            cond_raw = np.inf
            self._warnings.append(
                f"Raw covariance condition number failed for {model_label} model; reported as inf."
            )
        corr, cond_scaled = scaled_condition_number(cov)
        # flag highly correlated off-diagonal pairs

        flagged = []
        for i in range(n):
            for j in range(i + 1, n):
                if abs(corr[i, j]) >= HIGH_CORR_THRESHOLD:
                    flagged.append((names[i], names[j], float(corr[i, j])))

        # BB ↔ continuum correlations
        if analyze_complex:
            bb_corr = {}
            if self.bb_names:
                bb_idx = [i for i, nm in enumerate(names) if self._is_bb_param(nm)]
                cont_idx = [
                    i for i, nm in enumerate(names) if not self._is_bb_param(nm)
                ]
                for bi in bb_idx:
                    for ci in cont_idx:
                        bb_corr[(names[bi], names[ci])] = float(corr[bi, ci])

            return CovarianceAnalysis(
                condition_number=cond_raw,
                condition_number_scaled=cond_scaled,
                is_ill_conditioned=cond_scaled > 1e4,
                correlation_matrix=corr,
                param_names=names,
                flagged_pairs=flagged,
                bb_correlations=bb_corr if bb_corr else None,
            )
        else:
            return CovarianceAnalysis(
                condition_number=cond_raw,
                condition_number_scaled=cond_scaled,
                is_ill_conditioned=cond_scaled > 1e4,
                correlation_matrix=corr,
                param_names=names,
                flagged_pairs=flagged,
                bb_correlations=None,
            )

    # -- verdict helpers ------------------------------------------------------

    @staticmethod
    def _aic_label(delta: float) -> str:
        """Burnham & Anderson ΔAIC evidence label (magnitude)."""
        if not np.isfinite(delta):
            return "undefined"
        a = abs(delta)
        if a < 2:
            return "negligible"
        if a < 6:
            return "considerable"
        if a < 10:
            return "strong"
        return "decisive"

    @staticmethod
    def _bic_label(delta: float) -> str:
        """Kass & Raftery ΔBIC evidence label (magnitude)."""
        if not np.isfinite(delta):
            return "undefined"
        a = abs(delta)
        if a < 2:
            return "not worth mentioning"
        if a < 6:
            return "positive"
        if a < 10:
            return "strong"
        return "very strong"

    def _overall_verdict(self) -> str:
        evidence = self._evidence
        conflict = "; selection conflict" if evidence.selection_conflict else ""
        return f"{evidence.recommended_model}  ({evidence.recommendation_reason}{conflict})"

    @staticmethod
    def _covariance_to_dict(cov_a: Optional[CovarianceAnalysis]) -> Optional[dict]:
        if cov_a is None:
            return None
        return {
            "condition_number": cov_a.condition_number,
            "condition_number_scaled": cov_a.condition_number_scaled,
            "is_ill_conditioned": cov_a.is_ill_conditioned,
            "correlation_matrix": cov_a.correlation_matrix.tolist(),
            "param_names": list(cov_a.param_names),
            "flagged_pairs": [
                {"param_1": p1, "param_2": p2, "rho": rho}
                for p1, p2, rho in cov_a.flagged_pairs
            ],
            "bb_correlations": [
                {"bb_param": b, "continuum_param": c, "rho": rho}
                for (b, c), rho in (cov_a.bb_correlations or {}).items()
            ],
        }

    @staticmethod
    def _json_safe(value):
        if isinstance(value, dict):
            return {key: ModelComparison._json_safe(val) for key, val in value.items()}
        if isinstance(value, list):
            return [ModelComparison._json_safe(item) for item in value]
        if isinstance(value, tuple):
            return tuple(ModelComparison._json_safe(item) for item in value)
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        return value

    def summary_dict(self) -> dict:
        """Return model-comparison results in a machine-readable form."""
        summary = {
            "simple_model": {
                "name": self.simple.name,
                "k": self._k_simple,
                "cstat": self.simple.cstat,
                "dof": self.simple.dof,
            },
            "complex_model": {
                "name": self.complex.name,
                "k": self._k_complex,
                "cstat": self.complex.cstat,
                "dof": self.complex.dof,
            },
            "lrt": {
                "delta_cstat": self._lrt.delta_cstat,
                "delta_k": self._lrt.delta_k,
                "p_value": self._lrt.p_value,
                "sigma": self._lrt.sigma,
                "detected": self._lrt.detected,
            },
            "information_criteria": {
                "aic_simple": self._ic.aic_simple,
                "aic_complex": self._ic.aic_complex,
                "aic_c_simple": self._ic.aic_c_simple,
                "aic_c_complex": self._ic.aic_c_complex,
                "bic_simple": self._ic.bic_simple,
                "bic_complex": self._ic.bic_complex,
                "delta_aic": self._ic.delta_aic,
                "delta_aic_c": self._ic.delta_aic_c,
                "delta_bic": self._ic.delta_bic,
                "w_simple": self._ic.w_simple,
                "w_complex": self._ic.w_complex,
                "bic_evidence_ratio": self._ic.bayes_factor,
            },
            "bb_significance": [
                {
                    "param_name": bb.param_name,
                    "value": bb.value,
                    "error": bb.error,
                    "z_score": bb.z_score,
                    "significant": bb.significant,
                }
                for bb in self._bb_sig
            ],
            "parameter_stability": [
                {
                    "name": shift.name,
                    "v_simple": shift.v_simple,
                    "v_complex": shift.v_complex,
                    "sigma": shift.sigma,
                    "flag": shift.flag,
                }
                for shift in self._shifts
            ],
            "covariance_analysis": {
                "complex": self._covariance_to_dict(self._cov),
                "simple": self._covariance_to_dict(self._cov2),
            },
            "evidence_summary": {
                "lrt_supports_complex": self._evidence.lrt_supports_complex,
                "aic_supports_complex": self._evidence.aic_supports_complex,
                "bic_supports_complex": self._evidence.bic_supports_complex,
                "aic_evidence_ratio": self._evidence.aic_evidence_ratio,
                "bic_evidence_ratio": self._evidence.bic_evidence_ratio,
                "selection_conflict": self._evidence.selection_conflict,
                "recommended_model": self._evidence.recommended_model,
                "recommendation_reason": self._evidence.recommendation_reason,
            },
            "assumption_warnings": self.assumption_warnings,
        }
        return self._json_safe(summary)

    # -- public report --------------------------------------------------------

    def report(self, decimals: int = 3, save_path: Optional[str | Path] = None) -> None:
        """Print a formatted summary of all statistical tests, optionally saving it to disk."""
        if save_path is not None:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            with save_path.open("w", encoding="utf-8") as report_file:
                with redirect_stdout(_ReportWriter(sys.stdout, report_file)):
                    self.report(decimals=decimals)
            return

        sep = "-" * 68
        sep2 = "═" * 68

        def row(label, value, note=""):
            note_str = f"  [{note}]" if note else ""
            print(f"  {label:<34} {value}{note_str}")

        print(sep2)
        print(f"  MODEL COMPARISON")
        print(
            f"  Simple  :  {self.simple.name}   "
            f"(k={self._k_simple}, cstat={self.simple.cstat:.{decimals}f}, dof={self.simple.dof})"
        )
        print(
            f"  Complex :  {self.complex.name}  "
            f"(k={self._k_complex}, cstat={self.complex.cstat:.{decimals}f}, dof={self.complex.dof})"
        )
        print(sep2)

        # -- LRT --------------------------------------------------------------
        lrt = self._lrt
        print(f"\n  LIKELIHOOD RATIO TEST")
        print(sep)
        row(
            "Δcstat",
            f"{lrt.delta_cstat:.{decimals}f}",
            f"threshold = {DELTA_CSTAT_THRESHOLD}",
        )
        row("Δk", str(lrt.delta_k))
        row("p-value", f"{lrt.p_value:.3e}")
        row("Significance", f"{lrt.sigma:.2f}σ")
        row("Pipeline detect?", "YES ✓" if lrt.detected else "NO  ✗")

        # -- IC ---------------------------------------------------------------
        ic = self._ic
        print(f"\n  INFORMATION CRITERIA")
        print(sep)
        print(f"  {'':34} {'Simple':>12}   {'Complex':>12}")
        for label, vs, vc in [
            ("AIC", ic.aic_simple, ic.aic_complex),
            ("AICc", ic.aic_c_simple, ic.aic_c_complex),
            ("BIC", ic.bic_simple, ic.bic_complex),
        ]:
            print(f"  {label:<34} {vs:>12.{decimals}f}   {vc:>12.{decimals}f}")
        print(sep)
        row(
            "ΔAIC  (+ favours complex)",
            f"{ic.delta_aic:+.{decimals}f}",
            self._aic_label(ic.delta_aic),
        )
        row(
            "ΔAICc (+ favours complex)",
            f"{ic.delta_aic_c:+.{decimals}f}",
            self._aic_label(ic.delta_aic_c),
        )
        row(
            "ΔBIC  (+ favours complex)",
            f"{ic.delta_bic:+.{decimals}f}",
            self._bic_label(ic.delta_bic),
        )
        row("AIC weight — simple", f"{ic.w_simple * 100:.1f}%")
        row("AIC weight — complex", f"{ic.w_complex * 100:.1f}%")
        row("BIC evidence ratio (≈)", f"{ic.bayes_factor:.2f}", "> 1 favours complex")

        # -- evidence summary --------------------------------------------------
        evidence = self._evidence
        print(f"\n  EVIDENCE SUMMARY")
        print(sep)
        row("LRT supports complex?", "YES" if evidence.lrt_supports_complex else "NO")
        row("AIC supports complex?", "YES" if evidence.aic_supports_complex else "NO")
        row("BIC supports complex?", "YES" if evidence.bic_supports_complex else "NO")
        row(
            "AIC evidence ratio (≈)",
            f"{evidence.aic_evidence_ratio:.2f}",
            "> 1 favours complex",
        )
        row(
            "BIC evidence ratio (≈)",
            f"{evidence.bic_evidence_ratio:.2f}",
            "> 1 favours complex",
        )
        row("Selection conflict?", "YES" if evidence.selection_conflict else "NO")
        row(
            "Recommended model",
            evidence.recommended_model,
            evidence.recommendation_reason,
        )

        # -- BB significance ---------------------------------------------------
        # if self._bb_sig:
        #     print(f"\n  BB COMPONENT SIGNIFICANCE  (threshold = {BB_SIGMA_THRESHOLD}σ)")
        #     print(sep)
        #     for bb in self._bb_sig:
        #         flag = "✓" if bb.significant else "✗  below threshold"
        #         row(
        #             bb.param_name,
        #             f"{bb.value:.3e} ± {bb.error:.3e}",
        #             f"{bb.z_score:.2f}σ  {flag}",
        #         )

        # -- parameter stability -----------------------------------------------
        if self._shifts:
            print(f"\n  PARAMETER STABILITY")
            print(sep)
            print(f"  {'Parameter':<18} {'Simple':>14}  {'Complex':>14}  {'Shift':>8}")
            for s in self._shifts:
                flag = "  ← flag" if s.flag else ""
                print(
                    f"  {s.name:<18} {s.v_simple:>14.4g}  {s.v_complex:>14.4g}  "
                    f"{s.sigma:>6.2f}σ{flag}"
                )

        # -- covariance analysis -----------------------------------------------
        cov_a = self._cov
        if cov_a is not None:
            print(f"\n  COVARIANCE ANALYSIS  (complex model)")
            print(sep)
            row(
                "Condition number (scaled)",
                f"{cov_a.condition_number_scaled:.2e}",
                (
                    "ILL-CONDITIONED ✗"
                    if cov_a.is_ill_conditioned
                    else "well-conditioned ✓"
                ),
            )

            if cov_a.flagged_pairs:
                print(f"\n  High-correlation pairs  (|ρ| ≥ {HIGH_CORR_THRESHOLD}):")
                for p1, p2, rho in cov_a.flagged_pairs:
                    print(f"    {p1}  ↔  {p2}:  ρ = {rho:+.3f}  ← degeneracy flag")
            else:
                print("  No high-correlation pairs found.")

            if cov_a.bb_correlations:
                print(f"\n  BB ↔ continuum correlations:")
                for (b, c), rho in cov_a.bb_correlations.items():
                    warn = "  ← strong" if abs(rho) >= HIGH_CORR_THRESHOLD else ""
                    print(f"    {b}  ↔  {c}:  ρ = {rho:+.3f}{warn}")

        # -- covariance analysis -----------------------------------------------
        cov_a = self._cov2
        if cov_a is not None:
            print(f"\n  COVARIANCE ANALYSIS  (simple model)")
            print(sep)
            row(
                "Condition number (scaled)",
                f"{cov_a.condition_number_scaled:.2e}",
                (
                    "ILL-CONDITIONED ✗"
                    if cov_a.is_ill_conditioned
                    else "well-conditioned ✓"
                ),
            )

            if cov_a.flagged_pairs:
                print(f"\n  High-correlation pairs  (|ρ| ≥ {HIGH_CORR_THRESHOLD}):")
                for p1, p2, rho in cov_a.flagged_pairs:
                    print(f"    {p1}  ↔  {p2}:  ρ = {rho:+.3f}  ← degeneracy flag")
            else:
                print("  No high-correlation pairs found.")

            if cov_a.bb_correlations:
                print(f"\n  BB ↔ continuum correlations:")
                for (b, c), rho in cov_a.bb_correlations.items():
                    warn = "  ← strong" if abs(rho) >= HIGH_CORR_THRESHOLD else ""
                    print(f"    {b}  ↔  {c}:  ρ = {rho:+.3f}{warn}")

        # if self._warnings:
        #     print(f"\n  ASSUMPTIONS / WARNINGS")
        #     print(sep)
        #     for warning in self._warnings:
        #         print(f"  - {warning}")

        # -- overall verdict ---------------------------------------------------
        print(f"\n{sep2}")
        print(f"  VERDICT:  {self._overall_verdict()}")
        print(f"{sep2}\n")
