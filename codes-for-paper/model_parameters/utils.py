"""Created on Apr 02 14:28:01 2026"""

from typing import Optional, Tuple

import numpy as np
from scipy.odr import ODR, Model as ODRModel, RealData

from src.grb_research import GRBCatalog, find_project_root  # noqa: F401
from src.grb_research.grb_constants import (  # noqa: F401
    short_to_long,
    LEGEND_FONT_SIZE,
    LABEL_FONT_SIZE,
    LEGEND_TITLE_FONT_SIZE,
    TICK_FONT_SIZE,
)
from src.grb_research.grb_core import prepare_grbs  # noqa: F401
from src.grb_research.grb_model import ModelSet  # noqa: F401
from src.grb_research.grb_time import EpisodeTypes  # noqa: F401
from src.grb_research.grb_utils import break_e_to_e_peak, EpisodeMarkerResolver
from src.grb_research.grb_utils import plot_per_episode, save_value_error_as_parquet  # noqa: F401

# -- Re-exports (commonly needed across scripts) -----------------------------


# -- General parameter extractor ----------------------------------------------


def extract_parameter(model, param_pattern: str, *, return_asymmetric: bool = False):
    """Extract the first parameter whose name contains *param_pattern*.

    Parameters
    ----------
    model : Model
        GRB spectral model.
    param_pattern : str
        Substring to match in parameter names.
    return_asymmetric : bool
        If ``True``  → ``(error_lo, value, error_hi)`` (for errorbar plots).
        If ``False`` → ``(value, error)``.

    Returns
    -------
    tuple or None
        Extracted values, or ``None`` if no matching parameter is found.
    """
    for p in model.parameters:
        if param_pattern in p.name:
            if return_asymmetric:
                return p.error, p.value, p.error
            return p.value, p.error
    return None


# -- SBPL → Band conversion --------------------------------------------------


def sbpl_mask(index1_sbpl, index2_sbpl, e_break_sbpl):
    """Physical-validity mask for SBPL Monte-Carlo samples."""
    return np.logical_and(np.abs((index1_sbpl + index2_sbpl + 4) / (index1_sbpl - index2_sbpl)) < 1, e_break_sbpl > 0)


def convert_sbpl_to_band(model, n_sample: int = 10_000, seed=None, rng=None):
    """Convert SBPL break energy to Band E_peak via Monte-Carlo.

    Draws ``1.5 × n_sample`` multivariate-normal samples from the parameter
    covariance, applies a physical-validity filter, then computes percentiles
    of the derived E_peak distribution.

    Returns
    -------
    tuple of (error_lo, median, error_hi)
        From the 16 / 50 / 84 percentiles.
    """
    if seed is not None:
        rng = np.random.default_rng(seed)

    parameters = model.parameters
    cov_matrix = model.covariance_matrix_value
    raw = model.get_parameter_set.get_populated_values(cov_matrix, size=int(1.5 * n_sample), rng=rng)

    mvd = {p.name: raw[:, i] for i, p in enumerate(parameters)}

    mask = sbpl_mask(mvd["index1_sbpl"], mvd["index2_sbpl"], mvd["e_break_sbpl"])
    mvd_f = {k: v[mask] for k, v in mvd.items()}
    if mvd_f["index1_sbpl"].shape[0] < n_sample:
        raise ValueError("Not enough valid SBPL samples after physical filter.")

    idx = rng.choice(mvd_f["index1_sbpl"].shape[0], size=n_sample, replace=False)
    mvd_s = {k: v[idx] for k, v in mvd_f.items()}

    ep_samples = break_e_to_e_peak(
        index1_sbpl=mvd_s["index1_sbpl"], break_energy_sbpl=mvd_s["e_break_sbpl"], index2_sbpl=mvd_s["index2_sbpl"]
    )
    p = np.percentile(ep_samples, [16, 50, 84])
    return p[1] - p[0], p[1], p[2] - p[1]


# -- ODR linear-fit utilities ------------------------------------------------


def linear(params, x):
    """Linear model ``y = params[0] * x + params[1]`` for ODR."""
    return params[0] * x + params[1]


def fit_and_plot_odr(
    x_data,
    y_data,
    ax,
    *,
    mask=None,
    color: str = "#8B0000",
    linestyle: str = "--",
    annotation_xy: Tuple[float, float] = (0.05, 0.92),
    fontsize: float = LEGEND_FONT_SIZE,
    y_min_clip: Optional[float] = None,
    x_symbol: str = "kT",
    y_symbol: str = r"E_{\rm peak}",
):
    """Perform an ODR linear fit and draw the result with an uncertainty band.

    Uses full covariance-based uncertainty propagation for the confidence band.

    Parameters
    ----------
    x_data, y_data : np.ndarray, shape ``(N, 3)``
        Columns are ``(error_lo, value, error_hi)``.
    ax : matplotlib Axes
        Target axes for the plot.
    mask : array-like of bool, optional
        If given, only the selected rows are used for the fit.
    color : str
        Colour for fit line, fill, and annotation.
    linestyle : str
        Line style for the fit line.
    annotation_xy : tuple
        ``(x, y)`` in *axes fraction* for the equation annotation.
    fontsize : float
        Font size for the annotation.
    y_min_clip : float, optional
        If given, clips the lower uncertainty band to this value.
    x_symbol, y_symbol : str
        LaTeX math symbols for the independent/dependent variable in the
        annotated fit equation.

    Returns
    -------
    scipy.odr.Output
        The ODR result object.
    """
    x, y = x_data.copy(), y_data.copy()
    if mask is not None:
        x, y = x[mask], y[mask]

    x_centers = x[:, 1]
    y_centers = y[:, 1]
    x_errors = 0.5 * (x[:, 0] + x[:, 2])
    y_errors = 0.5 * (y[:, 0] + y[:, 2])

    data = RealData(x_centers, y_centers, sx=x_errors, sy=y_errors)
    odr = ODR(data, ODRModel(linear), beta0=[1, 1])
    result = odr.run()

    x_fine = np.linspace(x_centers.min(), x_centers.max(), 200)
    cov = result.cov_beta
    y_fit = linear(result.beta, x_fine)
    y_var = x_fine**2 * cov[0, 0] + cov[1, 1] + 2 * x_fine * cov[0, 1]
    y_err = np.sqrt(np.maximum(y_var, 0))

    ax.plot(x_fine, y_fit, color=color, ls=linestyle)

    lower = y_fit - y_err
    if y_min_clip is not None:
        lower = np.maximum(lower, y_min_clip)
    ax.fill_between(x_fine, lower, y_fit + y_err, color=color, alpha=0.1)

    ax.annotate(
        f"${y_symbol} = {result.beta[0]:+.1f}({result.sd_beta[0]:.1f})"
        f"\\cdot {x_symbol} {result.beta[1]:+.1f}({result.sd_beta[1]:.1f})$",
        xy=annotation_xy,
        xycoords="axes fraction",
        fontsize=fontsize,
        color=color,
    )

    return result


# -- Batch kT / E_peak extraction --------------------------------------------


def extract_kt_epeak_from_models(models, t90_marker="o", seed=1234):
    """Extract kT and E_peak arrays, markers, colours, and labels from models.

    For SBPL models the E_peak is derived via :func:`convert_sbpl_to_band`;
    for BAND / CPL models it is read directly from the ``e_peak`` parameter.

    Parameters
    ----------
    models : list[Model]
        Spectral models to process.
    t90_marker : str
        Marker to use for the T90 episode.
    seed : int
        Random seed for the SBPL Monte-Carlo conversion.

    Returns
    -------
    kt_values : np.ndarray, shape ``(N, 3)``
    ep_values : np.ndarray, shape ``(N, 3)``
    markers   : list[str]
    colors    : list[str]
    labels    : list[str]
    """
    resolver = EpisodeMarkerResolver(t90_marker=t90_marker)
    kt_values, ep_values = [], []
    markers, colors, labels = [], [], []

    for model_ in models:
        kt_values.append(extract_parameter(model_, "kt", return_asymmetric=True))

        if "SBPL" in model_.name:
            ep_values.append(convert_sbpl_to_band(model_, seed=seed))
        elif "BAND" in model_.name or "CPL" in model_.name:
            ep_values.append(extract_parameter(model_, "e_peak", return_asymmetric=True))
        else:
            raise ValueError(f">{model_.name}< is not expected for this GRB.")

        markers.append(resolver.resolve(model_.interval))
        colors.append(resolver.get_color(model_.interval))

        idx = "" if model_.interval.index is None else f"{model_.interval.index}"
        labels.append(f"{model_.interval.kind.value}{idx}" + r"$_\text{" + model_.name.replace("_", "+") + r"}$")

    return np.array(kt_values), np.array(ep_values), markers, colors, labels
