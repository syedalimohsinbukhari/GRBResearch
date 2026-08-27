"""Created on Apr 02 14:28:01 2026"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from scipy.odr import ODR, Model as ODRModel, RealData

from grb_research import EpisodeMarkerResolver, LEGEND_FONT_SIZE, break_e_to_e_peak

# -- Named parameter extraction ------------------------------------------------

# Parameter names in results.json carry the spectral-shape suffix directly
# (e.g. "index1_band", "index1_sbpl", "e_peak_cpl") -- Model.from_dictionary
# keeps the raw JSON key as the parameter name verbatim (grb_model.py:104), no
# normalisation happens. A couple of names have no family variant at all
# ("kt_bb", "add_index_pl") and are used exactly as given. Confirmed against
# every parameter name actually present in results.json (2026-08-27).
_FAMILY_SUFFIXES = ("band", "sbpl", "cpl", "pl")  # order matters: "sbpl"/"cpl" both contain "pl"
_NO_SUFFIX = {"add_index_pl"}
_ALIASES = {"kt": "kt_bb"}


def _resolve_parameter_name(model_name: str, param_pattern: str) -> str:
    """Map a base parameter pattern (e.g. ``"index2"``) to the exact results.json
    parameter name for *model_name*'s spectral family (e.g. ``"index2_band"``)."""
    if param_pattern in _NO_SUFFIX:
        return param_pattern
    if param_pattern in _ALIASES:
        return _ALIASES[param_pattern]

    lowered = model_name.lower()
    for family in _FAMILY_SUFFIXES:
        if family in lowered:
            return f"{param_pattern}_{family}"
    raise ValueError(f"Cannot resolve parameter {param_pattern!r} for model {model_name!r}")


def extract_parameter(model, param_pattern: str, *, return_asymmetric: bool = False):
    """Extract a named parameter's value and error via the model's ``get_parameter_set``.

    Parameters
    ----------
    model : Model
        GRB spectral model.
    param_pattern : str
        Base parameter name, resolved to the exact results.json name for
        *model*'s spectral family via :func:`_resolve_parameter_name`.
    return_asymmetric : bool
        If ``True``  → ``(error_lo, value, error_hi)`` (for errorbar plots).
        If ``False`` → ``(value, error)``.

    Returns
    -------
    tuple or None
        Extracted values, or ``None`` if no matching parameter exists for this model.
    """
    exact_name = _resolve_parameter_name(model.name, param_pattern)
    try:
        p = model.get_parameter_set[exact_name]
    except KeyError:
        return None

    if return_asymmetric:
        return p.error, p.value, p.error
    return p.value, p.error


# -- Per-GRB panel-plot boilerplate -------------------------------------------


# Per-GRB T90 marker, matching pe_er_photosphere.py's/gamma_comparison_plot.py's convention
# (EX0/EX1/TRn markers are already GRB-independent via EpisodeMarkerResolver itself).
T90_MARKERS = ["o", "s", "X", "D"]


@dataclass
class PanelData:
    """Per-GRB arrays needed by the value-vs-time panel plots.

    Replaces the interval-array / has_BB / episode-label / marker boilerplate that
    was previously duplicated near-verbatim in high_index.py, low_index.py and
    peak_energy.py.
    """

    start: np.ndarray
    end: np.ndarray
    diff: np.ndarray
    midpoint: np.ndarray
    has_bb: List[bool]
    episode_labels: List[str]
    model_names: List[str]
    markers: List[str]


def prepare_panel_data(grb_objs) -> List[PanelData]:
    """One :class:`PanelData` per GRB in *grb_objs*.

    Relies on ``grb.intervals`` iteration order matching both
    ``extract_interval_arrays``'s array order and ``get_all_best_models()``'s
    model order -- already an implicit assumption of the pre-rewrite code
    (has_BB and the value/error arrays were already zipped against the
    interval arrays positionally), just centralised here.
    """
    panels = []
    for i, grb in enumerate(grb_objs):
        start, end, diff, midpoint = grb.intervals.extract_interval_arrays(return_include=("diff", "midpoint"))
        best_models = grb.get_all_best_models()
        episode_labels = [
            interval.kind.value if interval.index is None else f"{interval.kind.value}{interval.index}"
            for interval in grb.intervals
        ]
        resolver = EpisodeMarkerResolver(t90_marker=T90_MARKERS[i % len(T90_MARKERS)])
        panels.append(
            PanelData(
                start=start,
                end=end,
                diff=diff,
                midpoint=midpoint,
                has_bb=["BB" in m.name for m in best_models],
                episode_labels=episode_labels,
                model_names=[m.name for m in best_models],
                markers=[resolver.resolve(interval) for interval in grb.intervals],
            )
        )
    return panels


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
    is_mc     : np.ndarray of bool
        True where ``ep_values`` came from the SBPL -> Band Monte-Carlo conversion
        (so its `seed`/n_samples provenance is meaningful); False for BAND/CPL rows
        read directly from a fitted ``e_peak`` parameter. ``kt_values`` is never MC.
    """
    resolver = EpisodeMarkerResolver(t90_marker=t90_marker)
    kt_values, ep_values, is_mc = [], [], []
    markers, colors, labels = [], [], []

    for model_ in models:
        kt_values.append(extract_parameter(model_, "kt", return_asymmetric=True))

        if "SBPL" in model_.name:
            ep_values.append(convert_sbpl_to_band(model_, seed=seed))
            is_mc.append(True)
        elif "BAND" in model_.name or "CPL" in model_.name:
            ep_values.append(extract_parameter(model_, "e_peak", return_asymmetric=True))
            is_mc.append(False)
        else:
            raise ValueError(f">{model_.name}< is not expected for this GRB.")

        markers.append(resolver.resolve(model_.interval))
        colors.append(resolver.get_color(model_.interval))

        idx = "" if model_.interval.index is None else f"{model_.interval.index}"
        labels.append(f"{model_.interval.kind.value}{idx}" + r"$_\text{" + model_.name.replace("_", "+") + r"}$")

    return np.array(kt_values), np.array(ep_values), markers, colors, labels, np.array(is_mc)
