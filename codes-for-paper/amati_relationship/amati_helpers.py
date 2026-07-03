"""Created on Jan 24 02:32:18 2026 — refactored Mar 27 2026"""

from typing import List, Sequence

import numpy as np
from numpy.typing import ArrayLike

from src.grb_research.grb_calculations import mc_e_iso_sampler, ModelResampler
from src.grb_research.grb_utils import break_e_to_e_peak, EpisodeMarkerResolver

# ---------------------------------------------------------------------------
# Normalization — single source of truth for axis units
# ---------------------------------------------------------------------------

EP_NORM = 1e3  # E_{i,peak} plotted in units of keV x 1e3 → x-axis in 1e3 keV
EI_NORM = 1e52  # E_iso plotted in units of erg x 1e52 → y-axis in 1e52 erg


def _episode_label(m) -> str:
    """Produce the legend label for a single model's episode."""
    kind_name = m.interval.kind.name
    if kind_name in ("T90", "EX0", "EX1"):
        return str(m.interval.kind)
    return f"{m.interval.kind}{m.interval.index}"


# Define base models that have an e_peak parameter
_EP_MODELS = {"band", "cpl"}
_NO_EP_MODELS = {"pl"}
_SBPL = "sbpl"


def _get_base_model_name(m_name: str) -> str:
    """
    Extract the base spectral model name from a compound model string.

    Supported compounds: BAND, BAND_BB, BAND_PL_BB,
                         CPL,  CPL_BB,  CPL_PL_BB

    Raises
    ------
    ValueError
        If the model is PL-only or PL+BB (no e_peak defined),
        or if no recognised base model is found.
    """
    components = [c.lower() for c in m_name.split("_")]

    # Check for PL-only or PL+BB — no e_peak exists
    base_candidates = [c for c in components if c not in ("bb", "pl")]

    if not base_candidates:
        raise ValueError(
            f"Model '{m_name}' contains no base spectral component with a "
            f"defined e_peak (PL-only or PL+BB models do not have e_peak). "
            f"Cannot compute Ep-Eiso for this model."
        )

    base = base_candidates[0]  # first non-PL, non-BB component

    if base not in _EP_MODELS:
        raise ValueError(f"Model '{m_name}': unrecognised base component '{base}'. " f"Expected one of {_EP_MODELS}.")

    return base


def _compute_ep_eiso(
    m, redshift: float, n_sample: int, n_grid: int, seed_number: int, rng
) -> tuple[ArrayLike, ArrayLike]:
    m_name = m.name
    pc = m.get_parameter_set
    pc_names = [p.name for p in pc]
    cov_ = 0.5 * (m.covariance_matrix_value + m.covariance_matrix_value.T)

    print(f"{m_name}")

    raw = pc.get_populated_values(cov_, size=n_sample, rng=rng)
    m_res = ModelResampler(model=m, samples=raw, rng=rng, destroy=True)
    raw = m_res.run_resampler()
    mvd = {v: raw[:, i] for i, v in enumerate(pc_names)}

    if "sbpl" in m_name.lower():
        ep_samples = break_e_to_e_peak(
            index1_sbpl=mvd["index1_sbpl"], break_energy_sbpl=mvd["e_break_sbpl"], index2_sbpl=mvd["index2_sbpl"]
        )
    else:
        base = _get_base_model_name(m_name)
        ep_samples = mvd[f"e_peak_{base.lower()}"]

    samples_arr = np.array(list(mvd.values())).T

    ############################################

    eiso_samples = mc_e_iso_sampler(
        m, redshift, n_samples=n_sample, n_grid=n_grid, method=2, samples=samples_arr, seed_number=seed_number, rng=rng
    )

    ep_intrinsic = ep_samples * (1 + redshift)
    return ep_intrinsic, eiso_samples


def _plot_model_point(
    m,
    redshift: int | float,
    marker: str,
    color: str,
    n_grid: int,
    n_sample: int,
    seed_number: int,
    rng,
    alpha: int | float,
    label: str,
    axis=None,
) -> tuple[int | float, int | float, np.ndarray, np.ndarray]:
    """
    Compute and draw a single (E_peak, E_iso) point with error bars.

    Returns
    -------
    p50_ep, p50_eiso :
        Median values, useful for building redshift tracks.
    x_err, y_err :
        Asymmetric 1-sigma errors shaped (2, 1) for E_peak and E_iso.
        First row is lower error, second row is upper error.
    """
    ep_s, ei_s = _compute_ep_eiso(m, redshift, n_sample, n_grid, seed_number, rng)

    p50_ei, p50_ep, x_err, y_err = percentile_calculator(ei_s, ep_s)

    # Normalise to plotting units
    p50_ep /= EP_NORM
    p50_ei /= EI_NORM
    x_err /= EP_NORM
    y_err /= EI_NORM

    if axis is not None:
        axis.scatter(p50_ep, p50_ei, marker=marker, s=25, color=color, alpha=alpha, label=label, zorder=3)
        axis.errorbar(p50_ep, p50_ei, xerr=x_err, ms=0, color=color, alpha=alpha, zorder=2, capsize=5)

    return p50_ep, p50_ei, x_err, y_err


def percentile_calculator(ei_s: ArrayLike, ep_s: ArrayLike, return_percentiles_only=False):
    p16_ep, p50_ep, p84_ep = np.percentile(ep_s, [16, 50, 84])
    p16_ei, p50_ei, p84_ei = np.percentile(ei_s, [16, 50, 84])

    x_err = np.array([[p50_ep - p16_ep], [p84_ep - p50_ep]])
    y_err = np.array([[p50_ei - p16_ei], [p84_ei - p50_ei]])

    if return_percentiles_only:
        return p16_ep, p50_ep, p84_ep, p16_ei, p50_ei, p84_ei
    else:
        return p50_ei, p50_ep, x_err, y_err


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def amati_relationship_dirirsa2019(
    e_iso_norm: int | float = 1e52,
    e_i_peak_norm: int | float = 950.0,
    k: int | float = 1.67,
    sigma_k: int | float = 0.16,
    m: int | float = 1.16,
    sigma_m: int | float = 0.37,
    sigma_ext: int | float = 0.47,
    sigmas: Sequence[int] = (1, 2, 3),
    x_lim: tuple = (10, 1e5),
    y_lim: tuple = (1e50, 1e55),
    num_points: int = 1_000,
    use_average: bool = False,
    axis=None,
) -> None:
    """Plot the Amati relation with confidence bands (FanaDirirsa et al. 2019).

    All values are plotted in normalised units (EP_NORM, EI_NORM) to match
    the data points produced by _plot_model_point.
    """
    if axis is None:
        raise ValueError("An axis must be provided.")

    e_i_peak = np.logspace(np.log10(x_lim[0]), np.log10(x_lim[1]), num=num_points)
    x = np.log10(e_i_peak / e_i_peak_norm)
    y = k + m * x
    sigma_y = np.sqrt(sigma_k ** 2 + x ** 2 * sigma_m ** 2 + sigma_ext ** 2)

    if use_average:
        sigma_y = np.mean(sigma_y)
    e_isotropic = (10 ** y) * e_iso_norm

    # Normalize to plotting units
    e_i_peak_plot = e_i_peak / EP_NORM
    e_isotropic_plot = e_isotropic / EI_NORM

    axis.loglog(e_i_peak_plot, e_isotropic_plot, lw=1, alpha=0.45, color="k")

    band_colors = ["#FFD166", "#FF9F43", "#FF6B6B"]
    for i, n_sigma in enumerate(sigmas):
        c = band_colors[i % len(band_colors)]
        e_upper = (10 ** (y + n_sigma * sigma_y)) * e_iso_norm / EI_NORM
        e_lower = (10 ** (y - n_sigma * sigma_y)) * e_iso_norm / EI_NORM
        axis.fill_between(e_i_peak_plot, e_lower, e_upper, color=c, alpha=0.1)
        axis.plot(e_i_peak_plot, e_lower, color=c, ls="--")
        axis.plot(e_i_peak_plot, e_upper, color=c, ls="--")

    axis.set_xlim(np.array(x_lim) / EP_NORM)
    axis.set_ylim(np.array(y_lim) / EI_NORM)


def plot_grbs_over_amati_relationship(
    best_model_list,
    redshift_list: List[float],
    t90_marker_list: List[str],
    n_grid: int = 10_000,
    n_sample: int = 10_000,
    seed_number: int = 0,
    alpha: float = 1.0,
    axis=None,
) -> tuple[list, list, list, list, list, list]:
    """
    Plot one or more GRBs on the Amati plane.

    Parameters
    ----------
    best_model_list :
        One ModelSet per GRB; each ModelSet is an iterable of Model objects.
    redshift_list :
        Redshift for each GRB, aligned with best_model_list.
    t90_marker_list :
        Matplotlib marker string for T90 of each GRB (the only per-GRB marker dimension).
        All other episode types use fixed global markers defined in EpisodeMarkerResolver.
    n_grid :
        Grid resolution for E_iso integration.
    n_sample :
        Posterior sample count.
    seed_number :
        Base RNG seed; incremented per model to keep runs reproducible.
    alpha :
        Scatter / errorbar opacity.
    axis :
        Target axes object. Required.

    Returns
    -------
    ep_total :
        E_peak median values.
    ei_total :
        E_iso median values.
    ep_labels :
        Episode labels.
    _models :
        Model names.
    ep_err_total :
        E_peak error arrays (2, 1) shaped.
    ei_err_total :
        E_iso error arrays (2, 1) shaped.
    """
    if axis is None:
        raise ValueError("An axis must be provided.")

    rng = np.random.default_rng(seed_number)

    ep_total, ei_total = [], []
    ep_err_total, ei_err_total = [], []
    ep_labels = []
    _models = []

    for index, (models, redshift, t90_marker) in enumerate(zip(best_model_list, redshift_list, t90_marker_list)):
        resolver = EpisodeMarkerResolver(t90_marker=t90_marker)
        for index2, m in enumerate(models):
            ep, ei, ep_err, ei_err = _plot_model_point(
                m=m,
                redshift=redshift,
                marker=resolver.resolve(m.interval),
                color=resolver.get_color(m.interval),
                n_grid=n_grid,
                n_sample=n_sample,
                seed_number=seed_number + index2,
                rng=rng,
                alpha=alpha,
                label=_episode_label(m),
                axis=axis,
            )
            ep_total.append(ep)
            ei_total.append(ei)
            ep_err_total.append(ep_err)
            ei_err_total.append(ei_err)
            ep_labels.append(_episode_label(m))
            _models.append(m.name)

    return ep_total, ei_total, ep_labels, _models, ep_err_total, ei_err_total


def plot_unknown_redshift_grb(
    models,
    t90_marker: str,
    z_values: Sequence[float] = (1, 3, 5, 7),
    n_grid: int = 10_000,
    n_sample: int = 10_000,
    seed_number: int = 0,
    axis=None,
) -> tuple[list, list, list, list]:
    """
    Plot a GRB with unknown redshift across several assumed z values.

    For each episode in *models*, four points are drawn (one per entry in
    z_values) and connected by a dashed line to show the locus of the burst
    on the Amati plane as a function of assumed redshift. Each episode uses
    the same marker scheme as in plot_grbs_over_amati_relationship so the
    reader can match episodes across panels.

    Parameters
    ----------
    models : ModelSet
        All episode models for the unknown-redshift GRB.
    t90_marker : str
        Matplotlib marker for the T90 episode of this GRB.
    z_values : sequence of float
        Redshift values to evaluate. Defaults to (1, 3, 5, 7).
    n_grid, n_sample, seed_number : int
        Passed through to the sampler.
    axis : matplotlib Axes
        Target axes object. Required.
    """
    if axis is None:
        raise ValueError("An axis must be provided.")

    rng = np.random.default_rng(seed_number)
    resolver = EpisodeMarkerResolver(t90_marker=t90_marker)

    ep_all, ei_all = [], []
    ep_name = []
    model_list = []

    for ep_idx, m in enumerate(models):
        marker = resolver.resolve(m.interval)
        color = resolver.get_color(m.interval)

        ep_track: List[int | float] = []
        ei_track: List[int | float] = []

        for z_idx, z in enumerate(z_values):
            # Label only on the first redshift so the legend has one entry per episode, not one per (episode × z).
            label = _episode_label(m) if z_idx == 0 else ""

            ep, ei, _, _ = _plot_model_point(
                m=m,
                redshift=z,
                marker=marker,
                color=color,
                n_grid=n_grid,
                n_sample=n_sample,
                seed_number=seed_number + ep_idx * len(z_values) + z_idx,
                rng=rng,
                alpha=0.75,
                label=label,
                axis=axis,
            )
            ep_track.append(ep)
            ei_track.append(ei)

            ep_all.append(ep)
            ei_all.append(ei)
            model_list.append(m.name)
            ep_name.append(_episode_label(m))

        # Connect the z-track for this episode on the correct axis
        axis.plot(ep_track, ei_track, ls="--", color=color, alpha=0.5, zorder=1)

    return ep_all, ei_all, ep_name, model_list
