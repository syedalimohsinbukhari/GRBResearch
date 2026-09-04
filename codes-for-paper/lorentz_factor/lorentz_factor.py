"""
Minimum Lorentz Factor Calculator
===================================
Computes Gamma_min using the gamma-gamma opacity method.

Reference: Lithwick & Sari (2001), ApJ, 555, 540
           Abdo et al. (2009), Science, 323, 1688

Formula (Lithwick & Sari 2001, Limit A, Table 1 with redshift corrections):

    tau_hat = 2.1e11 * [(d_L/7Gpc)^2 * (0.511)^(-alpha+1) * f_1]
              / [(delta_T/0.1s) * (alpha-1)]

    Gamma_min = tau_hat^(1/(2a+2))
                * (E_max/0.511)^((a-1)/(2a+2))
                * (1+z)^((a-1)/(a+1))

    where:
        alpha    = Lithwick & Sari photon index (positive) = -beta (high-energy index)
        f_1      = photon flux at 1 MeV [ph/cm^2/s/MeV]
        E_max    = highest LAT photon energy [MeV]; 0.511 MeV = m_e c^2
        delta_T  = variability timescale [s]
        z        = redshift
        d_L      = luminosity distance [cm]

Spectral parameters are read from ``results.json`` through the ``grb_research``
class API (``prepare_grbs`` -> ``GRB`` -> ``Model``) rather than being restated
here, so this script cannot drift out of sync with the fitted-model database.
Only quantities that do not live in ``results.json`` — the LAT photon
properties and the redshifts — are tabulated below.

Outputs:
    lorentz_results.csv   — all computed values
    lorentz_table.tex     — LaTeX table for paper
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from grb_research import draw_model_samples, find_project_root, get_rng, prepare_grbs, seed_from_name
from grb_research.grb_constants import model_n_pars, N_SAMPLES
from grb_research.grb_enums import GRBModelsCombinations as gmC
from grb_research.grb_sed import MODEL_MAP
from grb_research.grb_seds import band_function, cutoff_powerlaw, powerlaw, smoothly_broken_power_law

# Photon-spectrum function for each single component, keyed as in SpectralModels.
SED_FUNCTIONS = {
    gmC.PL: powerlaw,
    gmC.CPL: cutoff_powerlaw,
    gmC.BAND: band_function,
    gmC.SBPL: smoothly_broken_power_law,
}

# ─── Cosmology ────────────────────────────────────────────────────────────────
# Fana Dirirsa et al. (2019) — the single cosmology used across this project.
# Must stay reconciled with tex_files/section-5-data-analysis.tex and
# grb_calculations.py::mc_e_iso_sampler; do not introduce a third value.
H0 = 69.6
OM0 = 0.286
cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)

# ─── Sample ───────────────────────────────────────────────────────────────────
# The four GRBs of this paper. Everything else is keyed off this list.
GRB_LIST = ["080916C", "131014A", "140206B", "231129C"]

TEX_NAMES = {
    "080916C": r"\grbzeroeightzeroninesixteenC",
    "131014A": r"\grbthirteentenfourteenA",
    "140206B": r"\grbfourteenzerotwozerosixB",
    "231129C": r"\grbtwentythreeeleventwentynineC",
}

E_1MEV_KEV = 1000.0

# Parameter holding the high-energy photon index, per continuum model.
# The Lithwick & Sari index is its negative.
HIGH_ENERGY_INDEX_PARAM = {gmC.BAND: "index2_band", gmC.SBPL: "index2_sbpl", gmC.PL: "index1_pl"}

# Per-episode LAT products, read from LAT_analysis/lat_photons.csv, which is built from the gtburst/gtlike output files by LAT_analysis/txt_to_csv.py.
# These values used to be transcribed by hand from the paper's appendix table, so an edit to the table left this script silently stale -- BUGS.md OBS-08.
#
# Photon energies are MeV (the source field is "P > 0.9 Max (E) MeV"; the LAT selection floor here is 100 MeV).
# See BUGS.md BUG-18.
LAT_PHOTONS_CSV = find_project_root() / "LAT_analysis" / "lat_photons.csv"

# LAT detection below this TS is not a secure association, so the derived Gamma_min is flagged rather than tabulated alongside the secure ones.
# This is the threshold quoted in the table caption and in the appendix footnotes.
TS_SECURE_DETECTION = 25.0


def load_lat_photons(csv_path=LAT_PHOTONS_CSV):
    """Read the LAT table, returning per-episode photon data and the weak detections.

    Returns
    -------
    photons          : {short_name: {episode: (E_max [MeV], arrival time [s])}}
    low_significance : {(short_name, episode)} for every episode with TS < 25
    """
    table = pd.read_csv(csv_path)
    photons: dict[str, dict[str, tuple[float, float]]] = {}
    low_significance: set[tuple[str, str]] = set()

    for row in table.itertuples(index=False):
        short_name = row.grb_name.removeprefix("GRB")
        photons.setdefault(short_name, {})[row.episode] = (row.e_max_MeV, row.t_arr_s)
        if row.ts < TS_SECURE_DETECTION:
            low_significance.add((short_name, row.episode))

    return photons, low_significance


LAT_PHOTONS, LOW_SIGNIFICANCE = load_lat_photons()

# Spectroscopic redshifts; None means Gamma_min cannot be computed.
REDSHIFTS = {"080916C": 4.35, "131014A": None, "140206B": None, "231129C": None}

# Variability timescale: the episode's own duration is used for every episode, uniformly.
# It is an upper bound on the true variability timescale, so the resulting Gamma_min is a conservative (weaker) lower limit -- the right posture given these bounds is already weaker than published treatments.
# Populate this map to override an episode with a published value; mixing sources within one table was previously an inconsistency (see `lorentz_factor.md` section 4).
VARIABILITY_TIMESCALE: dict = {}

# Monte-Carlo settings, matching the rest of the project.
SEED = seed_from_name(__file__)
rng = get_rng(seed=SEED)
PERCENTILES = (16.0, 50.0, 84.0)


# ─── Spectral quantities from the fitted model ───────────────────────────────


def continuum_photon_flux(model_name, values, energy_kev):
    """Non-thermal continuum photon flux [ph/cm^2/s/keV] for one or many draws.

    Any blackbody component is dropped: the gamma-gamma opacity is set by the
    power-law photons that pair-produce with the LAT photon, not by the thermal
    component. Parameters are sliced in declaration order exactly as
    ``SpectralModels._evaluate_components`` does.

    ``values`` may be ``(n_pars,)`` or ``(n, n_pars)``; the result is ``(n,)``
    for a scalar energy.
    """
    values = np.atleast_2d(np.asarray(values, dtype=float))
    energy = np.atleast_1d(np.asarray(energy_kev, dtype=float))[None, :]

    key = gmC(model_name.lower())
    components = MODEL_MAP.get(key, (key,))

    total = np.zeros((values.shape[0], energy.shape[1]))
    idx = 0
    for component in components:
        n_pars = model_n_pars[component]
        pars = [values[:, idx + offset, None] for offset in range(n_pars)]
        idx += n_pars
        if component is gmC.BB:
            continue
        total = total + SED_FUNCTIONS[component](energy, *pars)

    return total[:, 0]


def f1_from_values(model_name, values):
    """Continuum photon flux at 1 MeV [ph/cm^2/s/MeV] for one or many draws."""
    return continuum_photon_flux(model_name, values, E_1MEV_KEV) * 1000.0


def high_energy_index_param_name(model):
    """Name of the parameter holding the continuum's high-energy index."""
    key = gmC(model.name.lower())
    continuum = key
    if key in MODEL_MAP:
        curved = [c for c in MODEL_MAP[key] if c in (gmC.BAND, gmC.SBPL)]
        continuum = curved[0] if curved else next(c for c in MODEL_MAP[key] if c is not gmC.BB)
    return HIGH_ENERGY_INDEX_PARAM.get(continuum)


def high_energy_index(model):
    """High-energy photon index beta of the continuum, read from the fitted model."""
    param = high_energy_index_param_name(model)
    return model.get_parameter_value(param) if param else None


# ─── Gamma_min formula ───────────────────────────────────────────────────────


def compute_tau_hat(alpha_LS, f_1, delta_T_s, z):
    """
    Dimensionless optical-depth quantity tau-hat (Lithwick & Sari 2001, eq. 4 / eq. 9).

    Shared by both Limit A (this module) and Limit B (``lorentz_factor_limit_b.py``) --
    the two limits differ only in how tau-hat is combined with E_max and z afterwards,
    not in tau-hat itself, so it is computed once here to avoid the formula drifting
    between the two files.

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
    tau_hat : float
    """
    d_L_cm = cosmo.luminosity_distance(z).cgs.value
    d_7Gpc = d_L_cm / (7.0 * 3.0857e27)

    return 2.1e11 * d_7Gpc ** 2 * (0.511) ** (-alpha_LS + 1) * f_1 / ((delta_T_s / 0.1) * (alpha_LS - 1))


def compute_gamma_min(alpha_LS, f_1, E_max_MeV, delta_T_s, z):
    """
    Compute a minimum Lorentz factor (Lithwick & Sari 2001, Limit A).

    Parameters
    ----------
    alpha_LS :
        L&S photon index (positive) = -beta (high-energy index)
    f_1 :
        photon flux at 1 MeV [ph/cm^2/s/MeV]
    E_max_MeV :
        highest LAT photon energy [MeV]
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

    e1 = 1.0 / (2 * alpha_LS + 2)
    e2 = (alpha_LS - 1) / (2 * alpha_LS + 2)
    e3 = (alpha_LS - 1) / (alpha_LS + 1)

    # 0.511 MeV = m_e c^2, matching the MeV convention used for f_1 and tau_hat above.
    gamma_min = tau_hat ** e1 * (E_max_MeV / 0.511) ** e2 * (1 + z) ** e3
    return gamma_min, tau_hat


# ─── COMPUTE ─────────────────────────────────────────────────────────────────


def episode_label(interval):
    """Short episode label, e.g. T90, EX0, TR1."""
    kind = interval.kind.name
    if kind in ("TR", "SP"):
        return f"{kind}{interval.index}"
    return kind


def variability_timescale(short_name, episode, interval):
    """Variability timescale for an episode, and how it was obtained.

    A published value is used where one exists. Otherwise the episode duration
    is adopted: it is an upper bound on the true variability timescale, so the
    resulting Gamma_min is a conservative lower limit rather than an optimistic
    one. Gamma_min depends only weakly on it, as delta_T^(-1/(2*alpha+2)).
    """
    key = (short_name, episode)
    if key in VARIABILITY_TIMESCALE:
        return VARIABILITY_TIMESCALE[key], "literature"
    return interval.end - interval.start, "duration"


def main():
    """Compute Gamma_min for every episode with LAT coverage, then tabulate."""
    root = find_project_root()
    _, _, grb_objects, _ = prepare_grbs(grb_list=GRB_LIST, result_file=root / "results.json", get_best=True)

    results = []
    header = f"{'GRB':<12}{'Ep.':<6}{'model':<10}{'alpha':>7}{'f_1':>12}{'dT [s]':>9}{'Gamma_min':>21}"
    print(f"\n{header}\n" + "-" * len(header))

    for short_name, grb in zip(GRB_LIST, grb_objects):
        redshift = REDSHIFTS[short_name]
        photons = LAT_PHOTONS[short_name]

        for model in grb.get_all_best_models():
            episode = episode_label(model.interval)
            if episode not in photons:
                continue

            e_max_mev, t_arr = photons[episode]
            delta_t, delta_t_source = variability_timescale(short_name, episode, model.interval)
            beta = high_energy_index(model)

            if redshift is None or beta is None or beta >= -1.0:
                # No redshift or no usable high-energy index: alpha_LS <= 1 makes the Lithwick & Sari expression singular.
                f_1 = alpha_ls = tau_hat = gamma = gamma_lo = gamma_hi = None
                print(f"GRB{short_name:<9}{episode:<6}{model.name:<10}{'—':>7}{'—':>12}{delta_t:>9.3f}{'—':>11}")
            else:
                f_1 = f1_from_values(model.name, [p.value for p in model.parameters])[0]
                alpha_ls = -beta
                gamma, tau_hat = compute_gamma_min(alpha_ls, f_1, e_max_mev, delta_t, redshift)

                # Statistical uncertainty: propagate the fit covariance through both f_1 and alpha, which are correlated because they come from the same spectral parameters.
                samples = draw_model_samples(model, n_samples=N_SAMPLES, rng=rng)
                index_position = [q.name for q in model.parameters].index(high_energy_index_param_name(model))
                alpha_draws = -samples[:, index_position]
                f1_draws = f1_from_values(model.name, samples)

                usable = np.isfinite(f1_draws) & (f1_draws > 0) & (alpha_draws > 1.0)
                gamma_draws, _ = compute_gamma_min(alpha_draws[usable], f1_draws[usable], e_max_mev, delta_t, redshift)
                gamma_draws = gamma_draws[np.isfinite(gamma_draws)]
                lo, med, hi = np.percentile(gamma_draws, PERCENTILES)
                gamma_lo, gamma_hi = med - lo, hi - med

                flag = "*" if (short_name, episode) in LOW_SIGNIFICANCE else ""
                print(
                    f"GRB{short_name:<9}{episode:<6}{model.name:<10}{alpha_ls:>7.3f}"
                    f"{f_1:>12.4e}{delta_t:>9.3f}{gamma:>8.0f} -{gamma_lo:<5.1f}+{gamma_hi:<5.1f}{flag}"
                )

            results.append(
                {
                    "GRB": f"GRB{short_name}",
                    "tex_name": TEX_NAMES[short_name],
                    "episode": episode,
                    "H0": H0,
                    "Om0": OM0,
                    "z": redshift,
                    "E_max_MeV": e_max_mev,
                    "t_arr_s": t_arr,
                    "t_v_s": delta_t,
                    "t_v_source": delta_t_source,
                    "model": model.name,
                    "beta": beta,
                    "alpha_LS": alpha_ls,
                    "f_1": f_1,
                    "tau_hat": tau_hat,
                    "Gamma_min": gamma,
                    "Gamma_min_err_lower": gamma_lo,
                    "Gamma_min_err_upper": gamma_hi,
                    "n_samples": N_SAMPLES,
                    "seed": SEED,
                    "low_significance": (short_name, episode) in LOW_SIGNIFICANCE,
                }
            )

    df = pd.DataFrame(results)
    df.drop(columns=["tex_name"]).to_csv("lorentz_results.csv", index=False)
    print("\nSaved: lorentz_results.csv")

    with open("lorentz_table.tex", "w") as handle:
        handle.write(build_latex_table(results))
    print("Saved: lorentz_table.tex")


def fmt(val, fmt_str):
    """Format a value for a math-mode table cell, or an unset marker if it is None."""
    return r"\ldots" if val is None else f"${format(val, fmt_str)}$"


def build_latex_table(results):
    """Render the per-episode Gamma_min table.

    Bursts without a redshift never yield a Gamma_min (z enters d_L in tau_hat), so
    they are dropped from the table entirely rather than shown as all-ellipsis rows
    -- the CSV keeps every row regardless, this is a table-only presentation choice.
    """
    results = [r for r in results if r["z"] is not None]

    # The dagger marking "duration adopted as t_v" is carried once, in the column
    # header, rather than repeated per cell -- true today because every remaining
    # row (all GRB080916C, VARIABILITY_TIMESCALE is empty) shares the same source.
    # Asserted rather than assumed: if a literature t_v is ever added for one of
    # these episodes, this must fail loudly instead of silently mislabeling it.
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
            rows += f"    \\multicolumn{{7}}{{l}}{{\\textbf{{{r['tex_name']}}}}} \\\\\n"
            current = r["GRB"]

        gamma = None if r["Gamma_min"] is None else round(r["Gamma_min"])
        # Markers go inside the math group, not appended as a second one.
        if gamma is None:
            gamma_str = r"\ldots"
        else:
            gamma_str = f"${gamma}_{{-{r['Gamma_min_err_lower']:.0f}}}^{{+{r['Gamma_min_err_upper']:.0f}}}"
            if r["low_significance"]:
                gamma_str += r"\,^{\ddagger}"
            gamma_str += "$"

        rows += (
            f"    {r['episode']} & {fmt(r['z'], '.2f')} & {fmt(r['E_max_MeV'] / 1e3, '.2f')} & "
            f"{fmt(r['t_arr_s'], '.2f')} & {fmt(r['t_v_s'], '.3f')} & {fmt(r['beta'], '.3f')} & {gamma_str} \\\\\n"
        )

    return (
        r"""\begin{table}[!ht]
\centering
\caption{Minimum bulk Lorentz factor $\Gamma_{\min}$ from the gamma-gamma opacity
condition, evaluated for every episode with \ac{LAT} coverage.
$E_{\rm GeV}$ is the highest-energy \ac{LAT} photon of that episode,
$t_{\rm arr}$ its arrival time relative to $T_0$,
$t_{\rm v}$ the variability timescale,
$\beta$ the high-energy photon index of the episode's best-fit model,
and $\Gamma_{\min}$ the derived lower limit~\citep{Lithwick2001}.
Errors on $\Gamma_{\min}$ are the statistical $1\sigma$ interval from $10^{4}$ Monte Carlo
draws of the spectral parameters; they are far smaller than the systematic
uncertainty from the choice of $t_{\rm v}$ and from the analytic approximation
adopted, and should not be read as the total uncertainty.
Only \grbzeroeightzeroninesixteenC\ has a confirmed spectroscopic redshift ($z = 4.35$);
the other three bursts lack a measured redshift, so $\Gamma_{\min}$ is undetermined for
them and they are omitted from this table.}
\label{tab:lorentz}
\resizebox{\columnwidth}{!}{
\begin{threeparttable}
\renewcommand{\arraystretch}{1.25}
\begin{tabular}{lcccccc}
\toprule
Episode & $z$ & $E_{\rm GeV}$ [GeV] & $t_{\rm arr}$ [s] &
    $t_{\rm v}$ [s]$^{\dagger}$ & $\beta$ & $\Gamma_{\min}$ \\
\midrule
"""
        + rows
        + r"""\bottomrule
\end{tabular}
\begin{tablenotes}
\footnotesize
\item[$\dagger$] Episode duration adopted as an upper bound on the variability timescale, giving a conservative $\Gamma_{\min}$.
\item[$\ddagger$] \ac{LAT} detection with $\mathrm{TS} < 25$; the highest-energy photon association is not secure.
\end{tablenotes}
\end{threeparttable}
}
\end{table}
"""
    )


if __name__ == "__main__":
    main()
