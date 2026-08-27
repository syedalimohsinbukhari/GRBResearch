"""Created on Aug 30 15:14:03 2025"""

from typing import Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt
from uncertainties import unumpy as unp

from .grb_constants import MODEL_COLORS, kev_to_erg
from .grb_enums import GRBModelsCombinations as gmC


def powerlaw(energy, amp, e_piv, index1):
    return amp * (energy / e_piv) ** index1


def smoothly_broken_power_law(energy, amp, e_piv, index1, break_energy, delta, index2):
    m = (index2 - index1) / 2
    b = (index1 + index2) / 2

    # break_energy may be an array of Monte-Carlo draws, so the invalid case is
    # masked rather than short-circuited: a scalar negative break energy still
    # yields an all-NaN result, exactly as before. The placeholder keeps the
    # log of a non-positive number from being evaluated at all.
    break_energy = np.asarray(break_energy, dtype=float)
    invalid = break_energy <= 0
    safe_break = np.where(invalid, 1.0, break_energy)

    a = np.log10(energy / safe_break) / delta
    beta = m * delta * np.log(0.5 * (np.exp(a) + np.exp(-a)))

    a_piv = np.log10(e_piv / safe_break) / delta
    beta_piv = m * delta * np.log(0.5 * (np.exp(a_piv) + np.exp(-a_piv)))

    spectrum = amp * (energy / e_piv) ** b * 10.0 ** (beta - beta_piv)
    return np.where(invalid, np.nan, spectrum)


def band_function(energy, amp, e_peak, index1, index2):
    e_piv = 100.0
    i1_minus_i2 = index1 - index2
    break_ = e_peak / (index1 + 2.0)
    transition_condition = i1_minus_i2 * break_

    f1 = _cpl_one(energy=energy, e_peak=e_peak, index1=index1, e_piv=e_piv)

    f21 = (i1_minus_i2 * break_) / e_piv
    f2 = f21**i1_minus_i2 * np.exp(-i1_minus_i2) * _pl_one(energy=energy, e_piv=e_piv, index1=index2)

    return amp * np.where(energy < transition_condition, f1, f2)


def cutoff_powerlaw(energy, amp, e_peak, index1, e_piv):
    exp_ = -energy * (2 + index1)
    return amp * _pl_one(energy=energy, e_piv=e_piv, index1=index1) * np.exp(exp_ / e_peak)


def black_body(energy, amp, temperature):
    kt_clip = np.clip(a=energy / temperature, a_min=0, a_max=325)
    return amp * energy**2 / (np.exp(kt_clip) - 1)


def plot_model(
    x,
    model_values,
    model_strings,
    styles,
    plot_labels=None,
    x_lims=None,
    x_label=None,
    y_lims=None,
    y_label=None,
    axis=None,
    use_ergs=False,
):
    kev_to_ergs = kev_to_erg if use_ergs else 1
    if axis is None:
        f, ax = plt.subplots(figsize=(8, 6))
    else:
        ax = axis

    # -- convert keV to ergs if requested
    max_y = 1
    model_values = _swap(model_values)
    model_strings = _swap(model_strings)
    styles = _swap(styles)

    plot_labels = model_strings if plot_labels is None else plot_labels

    for index, (mv, label, style) in enumerate(zip(model_values, plot_labels, styles)):
        y_nom = unp.nominal_values(mv)
        y_err = unp.std_devs(mv)

        y_plot = y_nom * kev_to_ergs * x**2
        y_upper = (y_nom + y_err) * kev_to_ergs * x**2
        y_lower = (y_nom - y_err) * kev_to_ergs * x**2

        max_y = max(max_y, y_upper.max())

        try:
            # Use Enum lookup for colors
            clean_name = model_strings[index].replace("+", "_").lower()
            color = MODEL_COLORS[gmC(clean_name)]
        except (KeyError, ValueError):
            color = "k"
        ax.loglog(x, y_plot, style, color=color, label=label)
        ax.fill_between(x=x, y1=y_lower, y2=y_upper, color=color, alpha=0.2)

    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(alpha=0.25, ls="--")
    ax.legend(loc="best")


def _swap(list_of_values):
    """
    Rotate the list to the right by one position.
    Used to move the 'total' model (usually the last in the list) to the first position.
    """
    if isinstance(list_of_values, tuple):
        list_of_values = list(list_of_values)
    list_of_values.insert(0, list_of_values[-1])
    list_of_values = list_of_values[:-1]
    return list_of_values


def plot_single_model(
    x,
    model_values,
    model_string: str,
    plot_labels: str = None,
    x_lims: Optional[Tuple] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    use_ergs: bool = False,
    axis=None,
):
    plot_model(
        x=x,
        model_values=[model_values],
        model_strings=[model_string.replace("_", "+").upper()],
        styles=["-"],
        plot_labels=plot_labels,
        x_lims=x_lims,
        x_label=x_label,
        y_label=y_label,
        axis=axis,
        use_ergs=use_ergs,
    )


def plot_double_model(
    x, model_values, model_string, plot_labels=None, x_lims=None, x_label=None, y_label=None, axis=None, use_ergs=False
):
    m1, m2 = model_string.split("_")
    if m2 == "pl":
        m1, m2 = m2, m1

    m_string = [m1.upper(), m2.upper(), model_string.replace("_", "+").upper()]
    styles = ["--", "--", "-"]

    plot_model(
        x=x,
        model_values=model_values,
        model_strings=m_string,
        styles=styles,
        plot_labels=plot_labels,
        x_lims=x_lims,
        x_label=x_label,
        y_label=y_label,
        axis=axis,
        use_ergs=use_ergs,
    )


def plot_triple_model(
    x, model_values, model_string, plot_labels=None, x_lims=None, x_label=None, y_label=None, use_ergs: bool = False
):
    m1, m2, m3 = model_string.split("_")

    m_string = [m2.upper(), m1.upper(), m3.upper(), model_string.replace("_", "+").upper()]
    styles = ["--", "--", "--", "-"]

    plot_model(
        x=x,
        model_values=model_values,
        model_strings=m_string,
        styles=styles,
        plot_labels=plot_labels,
        x_lims=x_lims,
        x_label=x_label,
        y_label=y_label,
        use_ergs=use_ergs,
    )


###############################################################################
# AUXILIARY FUNCTIONS
###############################################################################


def _pl_one(energy, e_piv, index1):
    return (energy / e_piv) ** index1


def _cpl_one(energy, e_peak, index1, e_piv):
    return cutoff_powerlaw(energy=energy, amp=1, e_peak=e_peak, index1=index1, e_piv=e_piv)
