"""
Minimum Lorentz Factor Calculator — Limit B
=============================================
Computes Gamma_min using the Compton-scattering-off-pair-produced-e+/- method
(Lithwick & Sari 2001, Limit B), as a separate, independent bound alongside
Limit A (``lorentz_factor.py``).

Reference: Lithwick & Sari (2001), ApJ, 555, 540

Formula (Lithwick & Sari 2001, Limit B, Table 2 eq. 8, with redshift corrections):

    tau_hat = 2.1e11 * [(d_L/7Gpc)^2 * (0.511)^(-alpha+1) * f_1]
              / [(delta_T/0.1s) * (alpha-1)]                          -- identical to Limit A

    Gamma_min = tau_hat^(1/(a+3)) * (1+z)^((a-1)/(a+3))

    where:
        alpha    = Lithwick & Sari photon index (positive) = -beta (high-energy index)
        f_1      = photon flux at 1 MeV [ph/cm^2/s/MeV]
        delta_T  = variability timescale [s]
        z        = redshift
        d_L      = luminosity distance [cm]

Unlike Limit A, Limit B does **not** depend on E_max: it bounds Gamma by requiring
the e+/- pairs created by photon annihilation to be Compton-thin, not by requiring
a *specific* observed photon to escape annihilation. The paper's own convention
(Table 3) is to report max(Limit A, Limit B) per burst -- this script keeps Limit B
in a separate CSV/table rather than folding it into lorentz_factor.py's output, so
neither can silently overwrite the other; combining them is left as a later, explicit
step (join on GRB + episode).

tau_hat is imported from ``lorentz_factor.py`` (``compute_tau_hat``) rather than
re-derived here, since it is algebraically identical between the two limits --
duplicating it would risk exactly the kind of drift BUGS.md already logs for other
formulas in this project.

Outputs:
    lorentz_results_limit_b.csv   — all computed values
    lorentz_table_limit_b.tex     — LaTeX table for paper
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from grb_research import draw_model_samples, find_project_root, prepare_grbs

from lorentz_factor import (
    GRB_LIST,
    LAT_PHOTONS,
    N_SAMPLES,
    PERCENTILES,
    REDSHIFTS,
    SEED,
    TEX_NAMES,
    H0,
    OM0,
    compute_tau_hat,
    episode_label,
    f1_from_values,
    high_energy_index,
    high_energy_index_param_name,
    variability_timescale,
)


def compute_gamma_min_limit_b(alpha_LS, f_1, delta_T_s, z):
    """
    Compute a minimum Lorentz factor (Lithwick & Sari 2001, Limit B).

    Parameters
    ----------
    alpha_LS :
        L&S photon index (positive) = -beta (high-energy index)
    f_1 :
        photon flux at 1 MeV [ph/cm^2/s/MeV]
    delta_T_s :
        variability timescale [s]
    z :
        redshift

    Returns
    -------
    gamma_min : float
    tau_hat   : float
    """
    tau_hat = compute_tau_hat(alpha_LS, f_1, delta_T_s, z)

    e1 = 1.0 / (alpha_LS + 3)
    e2 = (alpha_LS - 1) / (alpha_LS + 3)

    gamma_min = tau_hat**e1 * (1 + z) ** e2
    return gamma_min, tau_hat


def main():
    """Compute Limit-B Gamma_min for every episode with LAT coverage, then tabulate."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    results = []
    header = f"{'GRB':<12}{'Ep.':<6}{'model':<10}{'alpha':>7}{'f_1':>12}{'dT [s]':>9}{'Gamma_min_B':>21}"
    print(f"\n{header}\n" + "-" * len(header))

    for short_name, grb in zip(GRB_LIST, grb_objects):
        redshift = REDSHIFTS[short_name]
        photons = LAT_PHOTONS[short_name]

        for model in grb.get_all_best_models():
            episode = episode_label(model.interval)
            if episode not in photons:
                # Same episode set as Limit A, for direct per-episode comparability
                # (mirrors the paper's own Table 3, which lists both limits side by side).
                continue

            delta_t, delta_t_source = variability_timescale(short_name, episode, model.interval)
            beta = high_energy_index(model)

            if redshift is None or beta is None or beta >= -1.0:
                f_1 = alpha_ls = tau_hat = gamma = gamma_lo = gamma_hi = None
                print(f"GRB{short_name:<9}{episode:<6}{model.name:<10}{'—':>7}{'—':>12}{delta_t:>9.3f}{'—':>11}")
            else:
                f_1 = f1_from_values(model.name, [p.value for p in model.parameters])[0]
                alpha_ls = -beta
                gamma, tau_hat = compute_gamma_min_limit_b(alpha_ls, f_1, delta_t, redshift)

                samples = draw_model_samples(model, n_samples=N_SAMPLES, seed=SEED)
                index_position = [q.name for q in model.parameters].index(high_energy_index_param_name(model))
                alpha_draws = -samples[:, index_position]
                f1_draws = f1_from_values(model.name, samples)

                usable = np.isfinite(f1_draws) & (f1_draws > 0) & (alpha_draws > 1.0)
                gamma_draws, _ = compute_gamma_min_limit_b(alpha_draws[usable], f1_draws[usable], delta_t, redshift)
                gamma_draws = gamma_draws[np.isfinite(gamma_draws)]
                lo, med, hi = np.percentile(gamma_draws, PERCENTILES)
                gamma_lo, gamma_hi = med - lo, hi - med

                print(
                    f"GRB{short_name:<9}{episode:<6}{model.name:<10}{alpha_ls:>7.3f}"
                    f"{f_1:>12.4e}{delta_t:>9.3f}{gamma:>8.0f} -{gamma_lo:<5.1f}+{gamma_hi:<5.1f}"
                )

            results.append(
                {
                    "GRB": f"GRB{short_name}",
                    "tex_name": TEX_NAMES[short_name],
                    "episode": episode,
                    "H0": H0,
                    "Om0": OM0,
                    "z": redshift,
                    "t_v_s": delta_t,
                    "t_v_source": delta_t_source,
                    "model": model.name,
                    "beta": beta,
                    "alpha_LS": alpha_ls,
                    "f_1": f_1,
                    "tau_hat": tau_hat,
                    "Gamma_min_B": gamma,
                    "Gamma_min_B_err_lower": gamma_lo,
                    "Gamma_min_B_err_upper": gamma_hi,
                    "n_samples": N_SAMPLES,
                    "seed": SEED,
                }
            )

    df = pd.DataFrame(results)
    df.drop(columns=["tex_name"]).to_csv("lorentz_results_limit_b.csv", index=False)
    print("\nSaved: lorentz_results_limit_b.csv")

    with open("lorentz_table_limit_b.tex", "w") as handle:
        handle.write(build_latex_table(results))
    print("Saved: lorentz_table_limit_b.tex")


def fmt(val, fmt_str):
    """Format a value for a math-mode table cell, or an unset marker if it is None."""
    return r"\ldots" if val is None else f"${format(val, fmt_str)}$"


def build_latex_table(results):
    """Render the per-episode Limit-B Gamma_min table.

    Bursts without a redshift never yield a Gamma_min (z enters d_L in tau_hat), so
    they are dropped from the table entirely rather than shown as all-ellipsis rows
    -- the CSV keeps every row regardless, this is a table-only presentation choice.
    """
    results = [r for r in results if r["z"] is not None]

    # See lorentz_factor.py::build_latex_table for why this is asserted rather than
    # assumed: the header carries one dagger for all rows, valid only while every
    # remaining row's t_v is duration-sourced.
    t_v_sources = {r["t_v_source"] for r in results}
    assert t_v_sources == {"duration"}, (
        f"build_latex_table assumes every shown row's t_v is duration-sourced (header carries a single "
        f"dagger accordingly); found {t_v_sources}. Move the dagger back to per-cell if this is no longer true."
    )

    rows = ""
    current = None
    for r in results:
        if r["GRB"] != current:
            if current is not None:
                rows += "    \\midrule\n"
            rows += f"    \\multicolumn{{5}}{{l}}{{\\textbf{{{r['tex_name']}}}}} \\\\\n"
            current = r["GRB"]

        gamma = None if r["Gamma_min_B"] is None else round(r["Gamma_min_B"])
        if gamma is None:
            gamma_str = r"\ldots"
        else:
            gamma_str = f"${gamma}_{{-{r['Gamma_min_B_err_lower']:.0f}}}^{{+{r['Gamma_min_B_err_upper']:.0f}}}$"

        rows += (
            f"    {r['episode']} & {fmt(r['z'], '.2f')} & "
            f"{fmt(r['t_v_s'], '.3f')} & {fmt(r['beta'], '.3f')} & {gamma_str} \\\\\n"
        )

    return (
        r"""\begin{table}[!ht]
\centering
\caption{
Minimum bulk Lorentz factor $\Gamma_{\min}$ from Compton scattering off pair-produced $e^{\pm}$ (Lithwick \& Sari 2001, Limit B), evaluated for every episode with \ac{LAT} coverage -- the same episode set as the Limit A table (\ref{tab:lorentz}).
Unlike Limit A, this bound does not depend on the LAT photon energy, so it applies equally to episodes whose highest-energy photon is a low-significance association.
$t_{\rm v}$ the variability timescale, $\beta$ the high-energy photon index of the episode's best-fit model, and $\Gamma_{\min}$ the derived lower limit~\citep{Lithwick2001}.
Errors are the statistical $1\sigma$ interval from $10^{4}$ Monte Carlo draws of the spectral parameters; as with Limit A they are far smaller than the systematic uncertainty from the choice of $t_{\rm v}$ and the analytic approximation adopted.
Only \grbzeroeightzeroninesixteenC\ has a confirmed spectroscopic redshift ($z = 4.35$); the other three bursts lack a measured redshift, so $\Gamma_{\min}$ is undetermined for them and they are omitted from this table.
}
\label{tab:lorentz_limit_b}
\begin{threeparttable}
\renewcommand{\arraystretch}{1.25}
\resizebox{\columnwidth}{!}{
\begin{tabular}{lcccc}
\toprule
Episode & $z$ & $t_{\rm v}$ [s]$^{\dagger}$ & $\beta$ & $\Gamma_{\min}$ \\
\midrule
"""
        + rows
        + r"""\bottomrule
\end{tabular}
}
\begin{tablenotes}
\footnotesize
\item[$\dagger$] Episode duration adopted as an upper bound on the variability timescale, giving a conservative $\Gamma_{\min}$.
\end{tablenotes}
\end{threeparttable}
\end{table}
"""
    )


if __name__ == "__main__":
    main()
