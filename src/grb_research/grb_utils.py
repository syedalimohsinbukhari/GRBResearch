"""Created on Sep 20 12:48:54 2025"""

import json
import os
from itertools import chain
from pathlib import Path
from typing import Dict, Optional, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
from matplotlib.patches import Ellipse
from scipy import stats

from .grb_constants import MODEL_PARAMETERS, NOK_THRESHOLD, OK_THRESHOLD, model_n_pars
from .grb_enums import GRBModelsCombinations, ModelStatus
from .grb_model import Model

m_style = seaborn.color_palette("deep6")


def get_directories_in_current_folder(cur_dir=None):
    """Get a sorted list of relevant directories in the current folder."""
    current_directory = os.getcwd() if cur_dir is None else cur_dir
    all_entries = os.listdir(current_directory)
    directories = []
    for entry in all_entries:
        full_path = os.path.join(current_directory, entry)
        if np.logical_and(os.path.isdir(full_path), np.logical_or(entry.startswith("GRB"), entry.startswith("Ep"))):
            if "Research" not in entry:
                directories.append(entry)
    directories.sort()
    return directories


def make_json_safe(obj):
    """Recursively convert an object to a JSON-serializable format."""
    if isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [make_json_safe(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(make_json_safe(v) for v in obj)
    elif hasattr(obj, "tolist"):  # catches numpy arrays/scalars
        return obj.tolist()
    else:
        return obj


def deep_merge(d, u):
    """Deep merge dictionary u into dictionary d."""
    for k_, v_ in u.items():
        if isinstance(v_, dict):
            # ensure destination has a dict to merge into
            node = d.setdefault(k_, {})
            if isinstance(node, dict):
                deep_merge(d=node, u=v_)
            else:
                d[k_] = v_
        else:
            d[k_] = v_
    return d


def covariance_to_correlation(cov):
    """Convert a covariance matrix to a correlation matrix."""
    cov = np.asarray(a=cov, dtype=float)
    d = np.sqrt(np.diag(cov))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(a=d, b=d)
        corr[~np.isfinite(corr)] = 0.0  # handle 0/0 or inf cases
    return corr


def flatten_results(res_total, include_covariance=True):
    """Flatten nested results dictionary into a pandas DataFrame."""
    rows = []
    for grb, epochs in res_total.items():
        for ep, models in epochs.items():
            for model, params in models.items():
                status = params.get("_status", "NA")
                for p_name, val in params.items():
                    if p_name == "_status":
                        continue
                    if p_name in ["covariance_matrix", "correlation_matrix"] and not include_covariance:
                        continue  # skip it
                    if isinstance(val, (list, tuple)) and len(val) == 2 and all(np.isscalar(x) for x in val):
                        v, e = val
                    else:
                        v, e = val, None
                    rows.append(
                        {
                            "GRB": grb,
                            "epoch": ep,
                            "model": model,
                            "status": status,
                            "param": p_name,
                            "value": v,
                            "error": e,
                        }
                    )
    return pd.DataFrame(rows)


def filter_covariance(cov_matrix, param_names):
    """Remove parameters with ~0 variance values from covariance matrix and names."""
    cov_matrix = np.asarray(cov_matrix)
    diagonals = np.diag(cov_matrix)
    keep_idx = [i for i, v in enumerate(diagonals.tolist()) if not float(v) == 0.0]
    filtered_cov = cov_matrix[np.ix_(keep_idx, keep_idx)]
    filtered_names = [param_names[i] for i in keep_idx]
    return filtered_cov, filtered_names, keep_idx


def plot_covariance_corner(
    model: Model, seed: Optional[int] = None, rng: Optional[np.random.Generator] = None
):
    """
    Corner-style plot with histograms on the diagonal and covariance ellipses off-diagonal.

    Parameters
    ----------
    means : array-like, shape (N,)
        Mean values of the parameters.
    cov_matrix : array-like, shape (N, N)
        Covariance matrix.
    param_names : list of str
        Names of parameters (length N).
    seed : int, optional
        Random seed for reproducibility. Ignored if rng is provided.
    rng : np.random.Generator, optional
        Random number generator instance for reproducibility.
    """
    means = np.asarray([i.value for i in model.parameters], dtype=float)
    cov_matrix = np.asarray(model.covariance_matrix_value, dtype=float)
    param_names = [i.name for i in model.parameters]

    if cov_matrix.ndim != 2 or cov_matrix.shape[0] != cov_matrix.shape[1]:
        raise ValueError("cov_matrix must be a square 2D array")
    if len(means) != cov_matrix.shape[0]:
        raise ValueError("means length must match cov_matrix dimensions")
    if len(param_names) != cov_matrix.shape[0]:
        raise ValueError("param_names length must match cov_matrix dimensions")

    cov_matrix, param_names, keep_idx = filter_covariance(cov_matrix, param_names)
    means = means[keep_idx]
    n_params = len(param_names)

    if n_params == 0:
        raise ValueError("cov_matrix has no non-zero variance parameters to plot")

    # Get or create RNG
    if rng is not None:
        rng_instance = rng
    else:
        rng_instance = np.random.default_rng(seed)

    stds = np.sqrt(np.diag(cov_matrix))
    n_stds = range(1, 4)
    max_std = max(n_stds) + 1

    fig, axes = plt.subplots(n_params, n_params, figsize=(3 * n_params, 2.5 * n_params), constrained_layout=False)
    axes = np.atleast_2d(axes)

    for i in range(n_params):
        for j in range(n_params):
            ax = axes[i, j]

            if i == j:
                # 1D Gaussian histogram centered on mean
                vals = rng_instance.normal(loc=means[i], scale=stds[i], size=10_000)
                ax.hist(vals, bins=24, fc="w", ec="k", histtype="step")
                ax.set_xlim(means[i] - max_std * stds[i], means[i] + max_std * stds[i])
                m, s = np.mean(vals), np.std(vals)
                spec = "g" if abs(s) < 1 else "f"
                ax.set_title(rf"{param_names[i]}" + "\n" + rf"${m:.3{spec}}\pm{s:.3{spec}}$", fontsize=14)
            elif j < i:
                # covariance ellipse centered at (mean[j], mean[i])
                cov_2d = cov_matrix[np.ix_([j, i], [j, i])]
                cov_2d = (cov_2d + cov_2d.T) / 2

                vals_, vectors = np.linalg.eigh(cov_2d)
                order = vals_.argsort()[::-1]
                vals_, vectors = vals_[order], vectors[:, order]

                theta = np.degrees(np.arctan2(*vectors[:, 0][::-1]))

                samples = rng_instance.multivariate_normal(mean=[means[j], means[i]], cov=cov_2d, size=5_000)

                ax.scatter(samples[:, 0], samples[:, 1], s=1, color="k", alpha=0.3, zorder=1)
                for k, n_std in enumerate(n_stds[::-1]):
                    if n_std < max(n_stds):
                        width, height = 2 * n_std * np.sqrt(vals_)
                        ellipse = Ellipse(
                            xy=(float(means[j]), float(means[i])),
                            width=width,
                            height=height,
                            angle=theta,
                            ec="r" if k == 1 else "g",
                            fc="w",
                            lw=1.5,
                            zorder=2 + k,
                            ls="--" if k == 1 else "-.",
                        )
                        ax.add_patch(ellipse)

                ax.set_xlim(means[j] - max_std * stds[j], means[j] + max_std * stds[j])
                ax.set_ylim(means[i] - max_std * stds[i], means[i] + max_std * stds[i])
            else:
                ax.axis("off")

            if i == 0 and j == 0:
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")
            else:
                if i != n_params - 1:
                    ax.set_xticks([])
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel(param_names[j])

                if j != 0:
                    ax.set_yticks([])
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel(param_names[i] if i != 0 else "")

            ax.tick_params(axis="both", labelrotation=45)

    top, right = 0.96 if n_params > 4 else 0.94, 0.97

    fig.subplots_adjust(left=right * 0.1, bottom=top * 0.1, right=right, top=top, wspace=0.03, hspace=0.03)
    return fig, axes


def flattened_json(file_path="results.json"):
    """Load and flatten results from a JSON file."""
    with open(f"{file_path}", "r") as f:
        data = json.load(f)
    return flatten_results(res_total=data, include_covariance=True)


def query_data(data, grb_name: str, m_name: str, status: str = "both", epoch: str = None):
    """
    Query the dataset for GRB, model, and optional filters.

    Parameters
    ----------
    data : pd.DataFrame
        The dataset to query.
    grb_name : str
        Required GRB identifier.
    m_name : str
        Required model name.
    status : str, optional
        Status filter: 'SAFE', 'UNSAFE', or 'both' (default).
    epoch : str, optional
        Epoch filter (default is None).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    query_str = f"GRB == '{grb_name}' and model == '{m_name}'"

    if status.upper() in ("SAFE", "UNSAFE"):
        query_str += f" and status == '{status.upper()}'"
    elif status.lower() != "both":
        raise ValueError("status must be 'SAFE', 'UNSAFE', or 'both'")

    if epoch is not None:
        query_str += f" and epoch == '{epoch}'"

    result = data.query(query_str).reset_index(drop=True)
    return result


def two_scatter(
    start_list,
    end_list,
    val1,
    val2,
    err1,
    err2,
    ti_fmt="s",
    ti_color="r",
    tr_fmt="o",
    x_time=False,
    plot_axis=None,
    remove_extra=False,
):
    """Create a scatter plot comparing two sets of values with error bars."""
    for index, (v1, v2, e1, e2) in enumerate(zip(val1, val2, err1, err2)):
        err_crit2 = e2 / abs(v2)
        fmt_, color_ = tr_fmt, m_style[index % len(m_style)]

        if index == 0:
            fmt_, color_ = ti_fmt, ti_color

        if not x_time:
            err_crit1 = e1 / abs(v1)
            err_bad = np.logical_or(err_crit1 > NOK_THRESHOLD, err_crit2 > NOK_THRESHOLD)
            err_warn = np.logical_or(err_crit1 > OK_THRESHOLD, err_crit2 > OK_THRESHOLD)
        else:
            err_bad = err_crit2 > NOK_THRESHOLD
            err_warn = err_crit2 > OK_THRESHOLD

        alpha_mask = np.logical_or(start_list[index] < start_list[0], end_list[index] > end_list[0])

        if index == 0:
            alpha_ = 1.0
        else:
            alpha_ = 0.75 if alpha_mask else 1.0

        if err_bad:
            l_lim, u_lim = 1, 1
        elif err_warn:
            l_lim, u_lim = 0, 1
        else:
            l_lim, u_lim = 0, 0

        if remove_extra:
            if alpha_mask:
                continue

        plot_axis.errorbar(
            x=v1,
            y=v2,
            xerr=e1,
            yerr=e2,
            color=color_,
            fmt=fmt_,
            capsize=3,
            lolims=l_lim,
            uplims=u_lim,
            alpha=alpha_,
            label=f"{start_list[index]}_{end_list[index]}",
        )


def epoch_to_time(epochs, differences=False):
    """Convert epoch strings to numerical start and end times."""
    label, start, end = [], [], []
    for time_ in epochs:
        label.append(time_.split(" ")[0])
        ts_, te_ = map(float, time_.split(" ")[1].split("_"))
        start.append(ts_)
        end.append(te_)
    start, end = np.array(start), np.array(end)
    if differences:
        return {
            "episode_label": label,
            "start": start,
            "end": end,
            "difference": 0.5 * np.diff(a=[start, end], axis=0)[0],
            "midpoint": 0.5 * (start + end),
        }
    else:
        return label, start, end


def grb_characteristics(grb_df, model_name, epoch_difference=False):
    """Get characteristics of a GRB for a specific model."""
    name_enum = GRBModelsCombinations(model_name.lower())
    unique_epochs = grb_df["epoch"].unique()
    model_n_par = model_n_pars[name_enum]
    model_labels = MODEL_PARAMETERS[name_enum]
    epoch = epoch_to_time(epochs=unique_epochs, differences=epoch_difference)

    return unique_epochs, model_n_par, model_labels, epoch


def break_e_to_e_peak(index1_sbpl, break_energy_sbpl, index2_sbpl):
    f1 = index1_sbpl + 2
    f2 = -(index2_sbpl + 2)
    f3 = (0.3 / 2) * np.log(f1 / f2)
    return break_energy_sbpl * 10 ** f3


def plot_per_episode(values, errors, m_name, start, end, difference, midpoints, axes, has_BB=None,
                      episode_labels=None, model_names=None, markers=None):
    """Plot one GRB's values-vs-time panel: index 0 is the T90 reference band, the rest
    are per-episode points.

    ``episode_labels``/``model_names``, when given, must be parallel arrays to
    ``start``/``end``/``values`` (same index -> same episode) -- if provided, each
    point's legend entry identifies its episode and fitted model, per this project's
    legend convention (a bare "GRB<name>" entry alone doesn't say which episode/model a
    point is). Omit both to keep the old behaviour (single GRB-name legend entry only).

    ``markers``, when given, is a parallel array of per-episode marker shapes (e.g. from
    ``EpisodeMarkerResolver``), so EX0/EX1/TRn are visually distinguishable rather than
    all sharing the same dot. Omit to keep the old fixed "." marker for every point.
    """
    errors = np.asarray(errors)

    def _point_label(i):
        if episode_labels is None or model_names is None:
            return None
        return f"{episode_labels[i]}" + r"$_\text{" + model_names[i].replace("_", "+") + r"}$"

    axes.plot([start[0], end[0]], [values[0], values[0]], c="k", ls="--", lw=2, label=_point_label(0))

    if errors.ndim == 1:
        y_low = values[0] - errors[0]
        y_high = values[0] + errors[0]
    else:
        y_low = values[0] - errors[0, 0]
        y_high = values[0] + errors[1, 0]

    axes.fill_between(x=[start[0], end[0]], y1=y_low, y2=y_high, color="k", alpha=0.15)

    # --- Episode points ---
    for i, x in enumerate(midpoints[1:], start=1):
        if errors.ndim == 1:
            y_err = errors[i]
        else:
            y_err = errors[:, i: i + 1]  # (2, 1), symmetric or asymmetric

        col = "b" if (start[i] < start[0] or end[i] > end[0] + 0.064) else "g"

        axes.errorbar(
            x,
            values[i],
            xerr=difference[i],
            yerr=y_err,
            color=col,
            marker=markers[i] if markers is not None else ".",
            mfc='w' if has_BB[i] else col,
            # ms=10,
            capsize=5,
            label=_point_label(i),
        )


def _categorize_models_by_extension(base_name, base_containing_models):
    """Categorize models by number of extension components.

    Parameters
    ----------
    base_name : str
        Name of the base model.
    base_containing_models : dict
        Dictionary of models that contain the base model.

    Returns
    -------
    tuple
        (single_extension, double_extension) dictionaries.
    """
    single_extension = {}
    double_extension = {}

    for model_name, value in base_containing_models.items():
        additional_components = model_name.count("_") - base_name.count("_")

        if additional_components == 1:
            single_extension[model_name] = value
        elif additional_components == 2:
            double_extension[model_name] = value

    return single_extension, double_extension


def _evaluate_single_extension_models(base_value, single_extensions):
    """Evaluate single extension models against base model.

    Parameters
    ----------
    base_value : float
        C-stat value of base model.
    single_extensions : dict
        Dictionary of single extension models and their c-stat values.

    Returns
    -------
    dict
        Dictionary mapping model names to status values.
    """
    results = {}
    for model_name, value in single_extensions.items():
        diff = base_value - value
        if diff > 25:
            results[model_name] = ModelStatus.ACCEPTED.value
        else:
            results[model_name] = ModelStatus.REJECTED.value
    return results


def _evaluate_double_extension_models(base_name, base_value, double_extensions, single_results, base_containing_models):
    """Evaluate double extension models using hierarchy rules.

    Parameters
    ----------
    base_name : str
        Name of the base model.
    base_value : float
        C-stat value of base model.
    double_extensions : dict
        Dictionary of double extension models and their c-stat values.
    single_results : dict
        Results from single extension evaluation.
    base_containing_models : dict
        All models containing the base model.

    Returns
    -------
    dict
        Dictionary mapping model names to status values.
    """
    results = {}

    for model_name, value in double_extensions.items():
        # Extract the two extensions
        suffix = model_name.replace(base_name + "_", "")
        extensions = suffix.split("_")

        # Find corresponding single extension models
        relevant_singles = []
        for ext in extensions:
            single_name = f"{base_name}_{ext}"
            if single_name in base_containing_models:
                relevant_singles.append(single_name)

        if len(relevant_singles) == 2:
            # Check if any single extension was accepted
            single1_status = single_results.get(relevant_singles[0], ModelStatus.REJECTED.value)
            single2_status = single_results.get(relevant_singles[1], ModelStatus.REJECTED.value)

            if single1_status == ModelStatus.ACCEPTED.value or single2_status == ModelStatus.ACCEPTED.value:
                # Compare against best accepted single (threshold: 25)
                best_single_value = float("inf")
                for single_name in relevant_singles:
                    if single_results.get(single_name) == ModelStatus.ACCEPTED.value:
                        best_single_value = min(best_single_value, base_containing_models[single_name])

                if best_single_value - value > 25:
                    results[model_name] = ModelStatus.ACCEPTED.value
                else:
                    results[model_name] = ModelStatus.REJECTED.value
            else:
                # Neither single accepted, compare against BASE (threshold: 50)
                if base_value - value > 50:
                    results[model_name] = ModelStatus.ACCEPTED.value
                else:
                    results[model_name] = ModelStatus.REJECTED.value
        else:
            # Not a clear combination, use default comparison
            if base_value - value > 50:
                results[model_name] = ModelStatus.ACCEPTED.value
            else:
                results[model_name] = ModelStatus.REJECTED.value

    return results


def _determine_base_status(base_name, all_results, single_extensions):
    """Determine the status of the base model.

    Parameters
    ----------
    base_name : str
        Name of the base model.
    all_results : dict
        All model evaluation results.
    single_extensions : dict
        Single extension models.

    Returns
    -------
    int
        Status value for base model.
    """
    # Check if any model was accepted
    any_accepted = any(status == ModelStatus.ACCEPTED.value for status in all_results.values())

    if any_accepted:
        return ModelStatus.UNNECESSARY.value
    else:
        return ModelStatus.ACCEPTED.value


def cstat_threshold(delta_k: int, sigma: float = 5.0) -> float:
    """
    Compute the Δcstat detection threshold for a likelihood ratio test
    at a given significance level.

    Uses the fact that under H₀, Δcstat ~ χ²(Δk), so the threshold
    is the inverse χ² CDF at the desired significance level.

    Parameters
    ----------
    delta_k : int
        Difference in number of free parameters between complex and simple model.
    sigma : float
        Desired significance level in Gaussian sigma. Default is 5.0.

    Returns
    -------
    float
        Δcstat threshold value.

    Examples
    --------
    >>> cstat_threshold(delta_k=2, sigma=5.0)
    28.74
    >>> cstat_threshold(delta_k=2, sigma=4.56)
    ~25.0
    >>> cstat_threshold(delta_k=4, sigma=5.0)
    36.86
    """
    p_value = 2 * stats.norm.sf(sigma)  # two-tailed p-value for a given sigma
    threshold = stats.chi2.ppf(1 - p_value, df=delta_k)
    return threshold


def analyze_model_hierarchy(is_good: Dict) -> Dict[str, ModelStatus]:
    """
    Analyze model hierarchy with custom flags based on comparison rules.

    Rules:
    1. BASE must be present in the model name for comparison
    2. BASE -> BASE_XX requires > 25 improvement
    3. BASE -> BASE_XX_YY requires:
       - >50 if both BASE_XX and BASE_YY are REJECTED
       - >25 from best BASE_XX/BASE_YY if either is ACCEPTED
    """
    # Extract BASE information
    base_name, base_value = is_good["BASE"]

    # Initialize result dictionary
    results = {}

    # Separate models that contain BASE from those that don't
    base_containing_models = {}
    other_models = {}

    for key, (model_name, value) in is_good.items():
        if key == "BASE":
            continue

        if base_name in model_name:
            base_containing_models[model_name] = value
        else:
            other_models[model_name] = value

    # Mark non-BASE containing models as INVALID
    for model_name in other_models:
        results[model_name] = ModelStatus.INVALID.value

    # If no models contain BASE, mark BASE as ACCEPTED and return
    if not base_containing_models:
        results[base_name] = ModelStatus.ACCEPTED.value
        return results

    # Categorize models by extension count
    single_extension, double_extension = _categorize_models_by_extension(base_name, base_containing_models)

    # Evaluate single extension models (BASE_XX)
    single_results = _evaluate_single_extension_models(base_value, single_extension)
    results.update(single_results)

    # Evaluate double extension models (BASE_XX_YY)
    double_results = _evaluate_double_extension_models(
        base_name, base_value, double_extension, single_results, base_containing_models
    )
    results.update(double_results)

    # Determine BASE status
    results[base_name] = _determine_base_status(base_name, results, single_extension)

    # Handle any remaining models aren't yet processed
    for model_name, value in base_containing_models.items():
        if model_name not in results:
            additional_components = model_name.count("_") - base_name.count("_")
            # STATISTICAL THRESHOLD FOR delta_k = 2 and delta_k = 4
            threshold = cstat_threshold(delta_k=2 * additional_components, sigma=5.0)

            if base_value - value > threshold:
                results[model_name] = ModelStatus.ACCEPTED.value
            else:
                results[model_name] = ModelStatus.REJECTED.value

    return results


def save_value_error_as_parquet(list_of_ep, grb_names, list_of_values, list_of_errors, list_of_names, filename,
                                asym_errs=False):
    if asym_errs:
        if not isinstance(list_of_errors, tuple) or len(list_of_errors) != 2:
            raise ValueError("list_of_errors must be a tuple of two lists when asym_errs is True")
        if not all(isinstance(err, list) for err in list_of_errors):
            raise ValueError("list_of_errors must contain only lists")
        temp_ = [np.column_stack((i, j, k, l)) for i, j, k, l in
                 zip(list_of_ep, list_of_values, list_of_errors[0], list_of_errors[1])]
        df = pd.DataFrame(np.vstack([*temp_]), columns=["ep", "value", "error_low", "error_high"])
    else:
        temp_ = [np.column_stack((i, j, k)) for i, j, k in zip(list_of_ep, list_of_values, list_of_errors)]
        df = pd.DataFrame(np.vstack([*temp_]), columns=["ep", "value", "error"])

    df["grb_name"] = np.repeat(grb_names, [len(i) for i in list_of_values])
    df["best_model_name"] = list(chain.from_iterable(list_of_names))

    if asym_errs:
        df = df[["ep", "grb_name", "best_model_name", "value", "error_low", "error_high"]]
    else:
        df = df[["ep", "grb_name", "best_model_name", "value", "error"]]

    df.to_parquet(Path.cwd() / filename, index=False)


class EpisodeMarkerResolver:
    """Maps an episode interval to a matplotlib marker and colour.

    T90 — GRB-specific marker passed at construction (the only dimension
          that differs across GRBs).
    EX0 — star (``"*"``)
    EX1 — pentagon (``"p"``)
    TR *n* — shape from ``TR_MARKERS``, indexed by ``interval.index``
    SP *n* — shape from ``SP_MARKERS``, indexed by ``interval.index``

    Parameters
    ----------
    t90_marker : str
        Marker to use for T90 episodes of this GRB.
    """

    TR_MARKERS: List = ["v", "<", ">", "^", "d", "8", "h"]
    EX_MARKERS: List[str] = ["*", "p"]
    SP_MARKERS: List[str] = ["s", "D", "P"]

    def __init__(self, t90_marker: str) -> None:
        self.t90_marker = t90_marker

    def resolve(self, interval) -> str:
        """Return the marker for *interval* based on its kind and index."""
        kind_name = interval.kind.name
        if kind_name == "T90":
            return self.t90_marker
        if kind_name == "EX0":
            return self.EX_MARKERS[0]
        if kind_name == "EX1":
            return self.EX_MARKERS[1]
        if kind_name == "TR":
            return self.TR_MARKERS[interval.index % len(self.TR_MARKERS)]
        if kind_name == "SP":
            return self.SP_MARKERS[interval.index % len(self.SP_MARKERS)]
        raise ValueError(f"Unrecognised EpisodeType: {interval.kind}")

    @staticmethod
    def get_color(interval) -> str:
        """Return the colour string for *interval*.

        T90 → ``"r"``, TR → ``"g"``, SP → ``"k"``, anything else → ``"b"``.
        """
        kind_name = interval.kind.name
        if kind_name == "T90":
            return "r"
        if kind_name == "TR":
            return "g"
        if kind_name == "SP":
            return "k"
        return "b"
