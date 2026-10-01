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
    lorentz_table_limit_b.tex     — LaTeX table for paper, rendered by
                                     generate_lorentz_table_limit_b.py (invoked here via
                                     subprocess; rerun that script directly to regenerate
                                     the table alone, without repeating the MC)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from grb_research import draw_model_samples, find_project_root, get_rng, prepare_grbs, seed_from_name

from lorentz_factor import (
    GRB_LIST,
    LAT_PHOTONS,
    N_SAMPLES,
    PERCENTILES,
    PHOTON_E_MIN_MEV,
    PHOTON_ENERGY_CUT_EXEMPT_EPISODES,
    REDSHIFTS,
    TEX_NAMES,
    H0,
    OM0,
    compute_tau_hat,
    episode_label,
    f1_from_values,
    high_energy_index,
    high_energy_index_param_name,
    sample_split_normal,
    variability_timescale,
)

# Independent of lorentz_factor.py's SEED -- __file__ differs naturally, decorrelating Limit A vs Limit B draws.
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)


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

            e_max_mev, _ = photons[episode]
            if (
                short_name == "080916C"
                and e_max_mev < PHOTON_E_MIN_MEV
                and episode not in PHOTON_ENERGY_CUT_EXEMPT_EPISODES
            ):
                # Same energy floor as Limit A, so both tables keep the same episode set -- see
                # lorentz_factor.py's identical check and lorentz_factor.md sec 12.
                print(f"GRB{short_name:<9}{episode:<6}  skipped: E_max={e_max_mev:.1f} MeV < {PHOTON_E_MIN_MEV:.0f} MeV cut")
                continue

            delta_t, delta_t_source, delta_t_err_lo, delta_t_err_hi = variability_timescale(
                short_name, episode, model.interval
            )
            beta = high_energy_index(model)

            if redshift is None or beta is None or beta >= -1.0:
                f_1 = alpha_ls = tau_hat = gamma = gamma_lo = gamma_hi = None
                print(f"GRB{short_name:<9}{episode:<6}{model.name:<10}{'—':>7}{'—':>12}{delta_t:>9.3f}{'—':>11}")
            else:
                f_1 = f1_from_values(model.name, [p.value for p in model.parameters])[0]
                alpha_ls = -beta
                gamma, tau_hat = compute_gamma_min_limit_b(alpha_ls, f_1, delta_t, redshift)

                samples = draw_model_samples(model, n_samples=N_SAMPLES, rng=rng)
                index_position = [q.name for q in model.parameters].index(high_energy_index_param_name(model))
                alpha_draws = -samples[:, index_position]
                f1_draws = f1_from_values(model.name, samples)

                # Same t_v error propagation as Limit A -- see lorentz_factor.py::main() and
                # lorentz_factor.md; both limits share every t_v caveat (lorentz_factor.md sec 8.4).
                if delta_t_source == "norris":
                    delta_t_for_draws = sample_split_normal(delta_t, delta_t_err_lo, delta_t_err_hi, N_SAMPLES, rng)
                else:
                    delta_t_for_draws = delta_t

                usable = np.isfinite(f1_draws) & (f1_draws > 0) & (alpha_draws > 1.0)
                delta_t_usable = delta_t_for_draws[usable] if delta_t_source == "norris" else delta_t
                gamma_draws, _ = compute_gamma_min_limit_b(alpha_draws[usable], f1_draws[usable], delta_t_usable, redshift)
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
                    "t_v_err_lower_s": delta_t_err_lo,
                    "t_v_err_upper_s": delta_t_err_hi,
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

    subprocess.run(
        [sys.executable, str(Path(__file__).parent / "generate_lorentz_table_limit_b.py")], check=True
    )


if __name__ == "__main__":
    main()
