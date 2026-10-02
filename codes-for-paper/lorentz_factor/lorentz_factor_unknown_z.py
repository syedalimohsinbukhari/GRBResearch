"""Minimum Lorentz factor (Limit A), redshift sweep for the three bursts without a spectroscopic redshift
=========================================================================================================
Companion to ``lorentz_factor.py`` (Limit A) and ``lorentz_factor_limit_b.py`` (Limit B), which are
pinned to ``REDSHIFTS`` and so only ever cover GRB080916C.  This script evaluates Limit A for
GRB131014A, GRB140206B and GRB231129C at the assumed redshifts z = 1, 2, 3, 5, 7 -- the same sweep as
``amati_relationship.py``'s unknown-z table -- the way that table sweeps E_iso.

``lorentz_factor_limit_b_unknown_z.py`` is the Limit B counterpart.  It reuses ``run_sweep()`` from here but has
its own script name and therefore its own, independent seed, exactly as Limit A and B are decorrelated for GRB080916C.

Design, mirroring ``amati_relationship.py``:
  - One Monte Carlo draw set per episode (spectral parameters and, for a Norris t_v, t_v itself),
    reused across every z, so the z columns of a row are correlated rather than independently
    resampled.  z only enters through d_L in tau_hat and the (1+z) factors; nothing is refit.
  - Every episode with LAT coverage is kept: no photon-energy floor and no TS cut (unlike the
    GRB080916C tables), since cutting would leave GRB140206B nearly empty.  Episodes with
    TS < 25 are flagged (``low_significance``) and carry a double dagger in the Limit A table.
  - Episodes with no usable high-energy index (e.g. a CPL best fit) get empty Gamma columns.

Outputs:
    - lorentz_results_unknown_z.csv   -- all computed values, at z = 1, 2, 3, 5, 7
    - lorentz_curves_unknown_z.csv    -- the same on a dense z grid, for the comparison figure
    - lorentz_table_unknown_z.tex     -- rendered by generate_lorentz_table_unknown_z.py
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
    H0,
    LAT_PHOTONS,
    LOW_SIGNIFICANCE,
    N_SAMPLES,
    OM0,
    PERCENTILES,
    REDSHIFTS,
    compute_gamma_min,
    episode_label,
    f1_from_values,
    high_energy_index,
    high_energy_index_param_name,
    sample_split_normal,
    variability_timescale,
)
from lorentz_factor_limit_b import compute_gamma_min_limit_b

Z_VALUES = (1, 2, 3, 5, 7)
# Dense grid for the redshift curves in gamma_comparison_unknown_z_plot.py: log-spaced over the photospheric sweep's range
# (1-7), plus the five tabulated z so every table value is also an exact curve point.  The curves reuse each episode's
# existing draw set, so they consume no extra random numbers and the table values are unaffected.
Z_CURVE = np.unique(np.append(np.logspace(np.log10(1.0), np.log10(7.0), 60), Z_VALUES))

# Independent of lorentz_factor_limit_b_unknown_z.py's seed: that script's __file__ differs.
SEED = seed_from_name(__file__)


def _gamma(limit, alpha, f_1, e_max_mev, delta_t, z):
    """Gamma_min for the given limit ('A' or 'B'); returns only the Lorentz factor."""
    if limit == "A":
        return compute_gamma_min(alpha, f_1, e_max_mev, delta_t, z)[0]
    return compute_gamma_min_limit_b(alpha, f_1, delta_t, z)[0]


def run_sweep(limit, rng, seed):
    """Gamma_min for ``limit`` ('A' or 'B') at every swept z, for each no-redshift burst's LAT episodes.

    Shared by this script (Limit A) and ``lorentz_factor_limit_b_unknown_z.py`` (Limit B), each passing
    its own ``rng``/``seed`` so the two limits use independent draws.

    Returns ``(table, curve)``: ``table`` is one row per episode with Gamma at each of ``Z_VALUES``; ``curve`` is the same
    quantity (best fit, with the same percentile-based errors) on the dense ``Z_CURVE`` grid, one row per episode and z.
    """
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    results = []
    curve_rows = []
    header = f"{'GRB':<12}{'Ep.':<6}{'model':<10}{'alpha':>7}{'dT [s]':>9}" + "".join(
        f"{'z=' + str(z):>12}" for z in Z_VALUES
    )
    print(f"\nLimit {limit}\n{header}\n" + "-" * len(header))

    for short_name, grb in zip(GRB_LIST, grb_objects):
        if REDSHIFTS[short_name] is not None:
            continue  # known-z bursts are covered by lorentz_factor.py / lorentz_factor_limit_b.py

        photons = LAT_PHOTONS[short_name]

        for model in grb.get_all_best_models():
            episode = episode_label(model.interval)
            if episode not in photons:
                continue

            e_max_mev, t_arr = photons[episode]
            delta_t, delta_t_source, delta_t_err_lo, delta_t_err_hi = variability_timescale(
                short_name, episode, model.interval
            )
            beta = high_energy_index(model)

            row = {
                "GRB": f"GRB{short_name}",
                "episode": episode,
                "H0": H0,
                "Om0": OM0,
                "E_max_MeV": e_max_mev,
                "t_arr_s": t_arr,
                "t_v_s": delta_t,
                "t_v_source": delta_t_source,
                "t_v_err_lower_s": delta_t_err_lo,
                "t_v_err_upper_s": delta_t_err_hi,
                "model": model.name,
                "beta": beta,
                "alpha_LS": None,
                "f_1": None,
                "n_samples": N_SAMPLES,
                "seed": seed,
                "low_significance": (short_name, episode) in LOW_SIGNIFICANCE,
            }
            for z in Z_VALUES:
                row[f"Gamma_min_{limit}_z{z}"] = None
                row[f"Gamma_min_{limit}_z{z}_err_lower"] = None
                row[f"Gamma_min_{limit}_z{z}_err_upper"] = None

            if beta is None or beta >= -1.0:
                # alpha_LS <= 1 makes the Lithwick & Sari expression singular, as in lorentz_factor.py.
                print(f"GRB{short_name:<9}{episode:<6}{model.name:<10}{'-':>7}{delta_t:>9.3f}")
                results.append(row)
                continue

            f_1 = f1_from_values(model.name, [p.value for p in model.parameters])[0]
            alpha_ls = -beta
            row["alpha_LS"], row["f_1"] = alpha_ls, f_1

            # One draw set per episode, reused for every z.
            samples = draw_model_samples(model, n_samples=N_SAMPLES, rng=rng)
            index_position = [q.name for q in model.parameters].index(high_energy_index_param_name(model))
            alpha_draws = -samples[:, index_position]
            f1_draws = f1_from_values(model.name, samples)

            if delta_t_source == "norris":
                delta_t_for_draws = sample_split_normal(delta_t, delta_t_err_lo, delta_t_err_hi, N_SAMPLES, rng)
            else:
                delta_t_for_draws = np.full(N_SAMPLES, delta_t)

            usable = np.isfinite(f1_draws) & (f1_draws > 0) & (alpha_draws > 1.0)
            line = f"GRB{short_name:<9}{episode:<6}{model.name:<10}{alpha_ls:>7.3f}{delta_t:>9.3f}"

            for z in Z_VALUES:
                gamma = _gamma(limit, alpha_ls, f_1, e_max_mev, delta_t, z)
                draws = _gamma(limit, alpha_draws[usable], f1_draws[usable], e_max_mev, delta_t_for_draws[usable], z)
                draws = draws[np.isfinite(draws)]
                lo, med, hi = np.percentile(draws, PERCENTILES)

                key = f"Gamma_min_{limit}_z{z}"
                row[key] = gamma
                row[f"{key}_err_lower"] = med - lo
                row[f"{key}_err_upper"] = hi - med
                line += f"{gamma:>12.0f}"

            print(line + ("  *" if limit == "A" and row["low_significance"] else ""))
            results.append(row)

            for z in Z_CURVE:
                gamma = _gamma(limit, alpha_ls, f_1, e_max_mev, delta_t, z)
                draws = _gamma(limit, alpha_draws[usable], f1_draws[usable], e_max_mev, delta_t_for_draws[usable], z)
                lo, med, hi = np.percentile(draws[np.isfinite(draws)], PERCENTILES)
                curve_rows.append(
                    {
                        "GRB": f"GRB{short_name}",
                        "episode": episode,
                        "model": model.name,
                        "z": z,
                        f"Gamma_min_{limit}": gamma,
                        f"Gamma_min_{limit}_err_lower": med - lo,
                        f"Gamma_min_{limit}_err_upper": hi - med,
                        "seed": seed,
                    }
                )

    return pd.DataFrame(results), pd.DataFrame(curve_rows)


def main():
    """Limit A sweep: write the CSV, then render its table."""
    out_dir = Path(__file__).parent
    table, curve = run_sweep("A", get_rng(seed=SEED), SEED)
    table.to_csv(out_dir / "lorentz_results_unknown_z.csv", index=False)
    curve.to_csv(out_dir / "lorentz_curves_unknown_z.csv", index=False)
    print("\nSaved: lorentz_results_unknown_z.csv, lorentz_curves_unknown_z.csv")

    subprocess.run([sys.executable, str(out_dir / "generate_lorentz_table_unknown_z.py"), "A"], check=True)


if __name__ == "__main__":
    main()
