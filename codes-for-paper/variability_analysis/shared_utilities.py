"""Created on Sep 26 04:15:52 2026."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.axes import Axes
from matplotlib.ticker import AutoMinorLocator
from plotez.typing import NDArray

from grb_research import short_to_long, long_to_short, MARKER_SIZE
from grb_research.grb_utils import save_fig
from light_curves import lightcurve_data
from norris_fit import t_peak, tv_mc_summary, NorrisFitter, norris_pulse

GRB_PROPERTIES = tuple[dict[str, tuple[float, float]], tuple[float, float, float], tuple[float, float]]

ACTIVE_THRESHOLD_FRAC = 0.01
N_Y_TICKS = 6
Y_BUFFER_FRAC = 0.25
MINOR_TICKS_PER_MAJOR = N_Y_TICKS - 1
E_LOW, E_HIGH = 10.0, 400.0


def grb080916C_bounds(low: float = E_LOW, high: float = E_HIGH) -> GRB_PROPERTIES:
    """Return GRB080916C's episode boundaries, T90 plot bounds, and energy band.

    Parameters
    ----------
    low :
        Lower bound of the energy band in keV, by default `E_LOW`.
    high :
        Upper bound of the energy band in keV, by default `E_HIGH`.

    Returns
    -------
    tuple :
        - episode_bounds
           Maps each episode label (e.g. "TR1") to its `(start, end)` time bounds in seconds since trigger, read
           directly from `results.json`.
        - t05_t95_t0
            `(T05, T95, T0_MET_S)` -- T90 start/end (plot-axis bounding only) and the trigger time in Fermi MET seconds.
        - energy_band
            `(low, high)`, echoed back unchanged.
    """
    episode_bounds = {
        "EX0": (-0.128, 4.864),
        "TR1": (1.280, 4.864),
        "TR2": (4.864, 15.040),
        "TR3": (15.040, 55.296),
        "TR4": (55.296, 59.52),
        "TR5": (59.52, 64.256),
        "EX1": (59.520, 67.904),
    }

    t05, t95 = 1.280, 64.256
    t0_met = 243216766.62

    return episode_bounds, (t05, t95, t0_met), (low, high)


def grb131014A_bounds(low: float = E_LOW, high: float = E_HIGH) -> GRB_PROPERTIES:
    """Return GRB131014A's episode boundaries, T90 plot bounds, and energy band.

    Parameters
    ----------
    low :
        Lower bound of the energy band in keV, by default `E_LOW`.
    high :
        Upper bound of the energy band in keV, by default `E_HIGH`.

    Returns
    -------
    episode_bounds :
        Maps each episode label (e.g. "TR1") to its `(start, end)` time bounds in seconds since trigger,
        read directly from `results.json`.
    t05_t95_t0 :
        `(T05, T95, T0_MET_S)` -- T90 start/end (plot-axis bounding only) and the trigger time in Fermi MET seconds.
    energy_band :
        `(low, high)`, echoed back unchanged.
    """
    episode_bounds = {"EX0": (-0.192, 2.432), "TR1": (0.960, 2.432), "TR2": (2.432, 4.160), "EX1": (2.432, 6.976)}

    t05, t95 = 0.960, 4.160
    t0_met = 403420143.2

    return episode_bounds, (t05, t95, t0_met), (low, high)


def grb140206B_bounds(low: float = E_LOW, high: float = E_HIGH) -> GRB_PROPERTIES:
    """Return GRB140206B's episode boundaries, T90 plot bounds, and energy band.

    Parameters
    ----------
    low :
        Lower bound of the energy band in keV, by default `E_LOW`.
    high :
        Upper bound of the energy band in keV, by default `E_HIGH`.

    Returns
    -------
    episode_bounds :
        Maps each episode label (e.g. "TR1") to its `(start, end)` time bounds in seconds since trigger, read directly
        from `results.json`.
    t05_t95_t0 :
        `(T05, T95, T0_MET_S)` -- T90 start/end (plot-axis bounding only) and the trigger time in Fermi MET seconds.
    energy_band :
        `(low, high)`, echoed back unchanged.
    """
    episode_bounds = {
        "EX0": (4.288, 11.072),
        "TR1": (7.488, 11.072),
        "TR2": (11.072, 20.032),
        "TR3": (20.032, 26.752),
        "TR4": (26.752, 59.776),
        "TR5": (59.776, 100.032),
        "TR6": (100.032, 154.240),
    }

    t05, t95 = 7.488, 154.240
    t0_met = 413361375.84

    return episode_bounds, (t05, t95, t0_met), (low, high)


def grb231129C_bounds(low: float = E_LOW, high: float = E_HIGH) -> GRB_PROPERTIES:
    """Return GRB231129C's episode boundaries, T90 plot bounds, and energy band.

    Parameters
    ----------
    low :
        Lower bound of the energy band in keV, by default `E_LOW`.
    high :
        Upper bound of the energy band in keV, by default `E_HIGH`.

    Returns
    -------
    episode_bounds :
        Maps each episode label (e.g. "TR1") to its `(start, end)` time bounds in seconds since trigger, read directly
        from `results.json`.
    t05_t95_t0 :
        `(T05, T95, T0_MET_S)` -- T90 start/end (plot-axis bounding only) and the trigger time in Fermi MET seconds.
    energy_band :
        `(low, high)`, echoed back unchanged.
    """
    episode_bounds = {"EX0": (-0.192, 3.136), "TR1": (0.384, 3.136), "TR2": (3.136, 7.296), "EX1": (3.136, 10.048)}

    t05, t95 = 0.384, 7.296
    t0_met = 722977823.114

    return episode_bounds, (t05, t95, t0_met), (low, high)


def initialize_work(root_path: Path, grb_full_name: str) -> tuple[Path, Path, Path]:
    """Resolve the project root, light-curve directory, and this burst's data folder.

    Parameters
    ----------
    root_path :
        The calling script's own resolved path (`Path(__file__).resolve()`).
    grb_full_name :
        The light-curve directory name (e.g. "GRB080916009").

    Returns
    -------
    project_root :
        `GRBResearchWork/`.
    lc_directory :
        `GRBResearchWork/light_curves/`.
    grb_folder :
        `GRBResearchWork/light_curves/<grb_full_name>/`.
    """
    project_root = root_path.parent.parent.parent
    lc_directory = project_root / "light_curves"
    grb_folder = lc_directory / grb_full_name

    return project_root, lc_directory, grb_folder


def grb_name_resolution(grb_name: str) -> tuple[str, str]:
    """Resolve a burst's short paper-name suffix to its full and short names.

    Parameters
    ----------
    grb_name :
        Short paper-name suffix (e.g. "080916C"), a key in `grb_research.short_to_long`.

    Returns
    -------
    full :
        The light-curve directory name (e.g. "GRB080916009").
    paper_name :
        `grb_name` round-tripped through `long_to_short[full]` (e.g."080916C", without the "GRB" prefix).
    """
    full = short_to_long[grb_name]
    return full, long_to_short[full]


def get_summed_lc(grb_lc, low: float = E_LOW, high: float = E_HIGH, norm: bool = False):
    """Load and sum every NaI detector's background-subtracted light curve.

    Detector time grids are asserted identical before summing (BUG-23 convention: raw sum, no per-detector normalization).

    Parameters
    ----------
    grb_lc :
        The burst's light-curve directory, containing one `.dat` file per detector.
    low :
        Lower energy bound in keV, by default `E_LOW`.
    high :
        Upper energy bound in keV, by default `E_HIGH`.
    norm :
        If True, additionally peak-normalize the summed rate -- see `Returns`.
        Default False.

    Returns
    -------
    t :
        Time grid, seconds since trigger.
    y :
        - If `norm` is False, the summed background-subtracted count rate (counts/s).
        - If `norm` is True, `(y / y.max(), y.max())` -- the peak-normalized rate and the peak value, kept for rescaling
          fitted amplitudes back to physical units afterward.
    nai_dat_name :
        The detector names summed, e.g. `["n3", "n4"]`.
    """
    dat = [f for f in os.listdir(f"{grb_lc}") if f.endswith(".dat")]
    dat = [i.split(".")[0] for i in dat]
    nai_dat_name = sorted(i for i in dat if "n" in i)

    nai_data = [lightcurve_data(f"{grb_lc}/{i}.dat", low, high) for i in nai_dat_name]
    t1 = nai_data[0][0]

    for _det, (_t, _, _) in zip(nai_dat_name[1:], nai_data[1:]):
        assert np.array_equal(t1, _t), f"{_det}'s time grid differs from {nai_dat_name[0]}'s -- cannot sum"

    r1 = np.sum([r for _, r, _ in nai_data], axis=0)
    b1 = np.sum([b for _, _, b in nai_data], axis=0)
    y = r1 - b1

    if not norm:
        return t1, y, nai_dat_name
    else:
        max_y_cps = np.max(y)
        y_norm = y / max_y_cps
        return t1, (y_norm, max_y_cps), nai_dat_name


def lat_details(lat_fit_path, t0_met, min_photon_energy):
    r"""Load a burst's LAT photon list, keeping only photons above an energy floor.

    Parameters
    ----------
    lat_fit_path :
        Path to the burst's LAT FITS file.
    t0_met :
        Trigger time in Fermi Mission Elapsed Time (seconds), subtracted from each photon's arrival time so
        `photon_t_arr_s` is seconds since trigger.
    min_photon_energy :
        Energy floor in MeV; photons at or below this are dropped.

    Returns
    -------
    photon_energy :
        Photon energies in MeV, above `min_photon_energy`, sorted by arrival time.
    photon_t_arr_s :
        Matching photon arrival times, seconds since trigger.
    """
    lat = fits.open(lat_fit_path)[1].data
    photon_energy = np.asarray(lat["ENERGY"])
    photon_t_arr_s = np.asarray(lat["TIME"]) - t0_met
    _e_mask = photon_energy > min_photon_energy
    photon_energy, photon_t_arr_s = photon_energy[_e_mask], photon_t_arr_s[_e_mask]
    _order = np.argsort(photon_t_arr_s)
    photon_energy, photon_t_arr_s = photon_energy[_order], photon_t_arr_s[_order]

    return photon_energy, photon_t_arr_s


grb_bounds = {
    "GRB080916C": grb080916C_bounds,
    "GRB131014A": grb131014A_bounds,
    "GRB140206B": grb140206B_bounds,
    "GRB231129C": grb231129C_bounds,
}


def save_data(
    fitter: NorrisFitter,
    parameters,
    grb_paper_name,
    t_range: tuple,
    max_y_cps,
    photon_summary,
    n_pulses,
    seed_number: int,
):
    """Build the per-(pulse, episode) results table and write it to CSV.

    Parameters
    ----------
    fitter :
        The converged joint fit, supplying `.params` and `.covariance`.
    parameters :
        Fitted parameters reshaped to `(n_pulses, 4)`, one `(A, t_s, tau1, tau2)` row per pulse.
    grb_paper_name :
        The burst's paper name (e.g. "GRB080916C") -- used both as the CSV's `grb_name` column value and as the lookup
        key into `grb_bounds` for this burst's episode boundaries and energy band.
    t_range :
        `(fit_window_start_s, fit_window_end_s)` actually fit.
    max_y_cps :
        Peak count rate used to rescale normalized amplitudes back to counts/s.
    photon_summary :
        Per-pulse photon summary from `get_pulse_summary()`, keyed by 1-indexed pulse number.
    n_pulses :
        Number of pulses in the joint fit.
    seed_number :
        Per-script deterministic seed (`grb_research.seed_from_name`) for the MC propagation of `t_v`.
    """
    episode_bounds, _, _energy = grb_bounds[grb_paper_name]()
    e_low, e_high = _energy

    rows = []
    for i in range(n_pulses):
        A, ts, tau1, tau2 = parameters[i]
        tp = t_peak(t_s=ts, tau1=tau1, tau2=tau2)
        mc = tv_mc_summary(fitter=fitter, pulse_index=i + 1, seed=seed_number)
        matched = [name for name, (start, end) in episode_bounds.items() if start <= tp <= end]
        for episode in matched or [None]:
            rows.append(
                {
                    "grb_name": grb_paper_name,
                    "episode": episode,
                    "pulse_index": i + 1,
                    "n_pulses_total": n_pulses,
                    "fit_window_start_s": t_range[0],
                    "fit_window_end_s": t_range[1],
                    "energy_low_keV": e_low,
                    "energy_high_keV": e_high,
                    "A_norm": A,
                    "A_cts_per_s": A * max_y_cps,
                    "y_max_cts_per_s": max_y_cps,
                    "t_s": ts,
                    "tau1": tau1,
                    "tau2": tau2,
                    "t_peak_s": tp,
                    "t_v_s": mc["t_v_s"],
                    "t_v_err_lower_s": mc["t_v_err_lower_s"],
                    "t_v_err_upper_s": mc["t_v_err_upper_s"],
                    "mc_kept_fraction": mc["kept_fraction"],
                    "n_samples": mc["n_samples"],
                    "seed": mc["seed"],
                    "n_lat_photons_assigned": photon_summary[i + 1]["n_photons"],
                    "photon_e_max_MeV": photon_summary[i + 1]["photon_e_max_MeV"],
                    "photon_t_arr_at_e_max_s": photon_summary[i + 1]["photon_t_arr_at_e_max_s"],
                }
            )

    results_df = pd.DataFrame(rows)
    csv_path = Path(__file__).parent / f"norris_fit_results_{grb_paper_name}.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"wrote {csv_path} ({len(results_df)} rows)")


def assign_pulse(t_arr: float, parameters: NDArray, n_pulses: int):
    """Assign one photon arrival time to a pulse by temporal proximity: the pulse with the latest onset `t_s <= t_arr`
    among those still active at `t_arr` (own value >= `ACTIVE_THRESHOLD_FRAC` of own peak).

    (Not "dominant flux": that earlier rule was replaced -- see variability_analysis.md. `lorentz_factor.py` reimplements
    this rule in `assigned_pulse()` and must stay in sync with it.)

    Parameters
    ----------
    t_arr :
        Photon arrival time, seconds since trigger.
    parameters :
        Fitted parameters reshaped to `(n_pulses, 4)`.
    n_pulses :
        Number of pulses.

    Returns
    -------
    int or None :
        1-indexed pulse number the photon is assigned to, or None if no pulse is active at `t_arr`.
    """
    pulse_params = [tuple(parameters[i]) for i in range(n_pulses)]
    pulse_peak_values = [
        norris_pulse(np.array([t_peak(ts, tau1, tau2)]), p)[0]
        for p, (A, ts, tau1, tau2) in zip(pulse_params, pulse_params)
    ]

    candidates = []
    for i, p in enumerate(pulse_params):
        ts = p[1]
        if ts > t_arr:
            continue
        val = norris_pulse(np.array([t_arr]), p)[0]
        if pulse_peak_values[i] > 0 and (val / pulse_peak_values[i]) >= ACTIVE_THRESHOLD_FRAC:
            candidates.append(i)
    return (max(candidates, key=lambda j: pulse_params[j][1]) + 1) if candidates else None


def get_pulse_summary(lat_t_arr: NDArray, parameters: NDArray, n_pulses: int, photon_energy: NDArray):
    """Summarize LAT photons per pulse: count and the highest-energy photon.

    Parameters
    ----------
    lat_t_arr :
        Photon arrival times, seconds since trigger.
    parameters :
        Fitted parameters reshaped to `(n_pulses, 4)`.
    n_pulses :
        Number of pulses.
    photon_energy :
        Photon energies in MeV, same order as `lat_t_arr`.

    Returns
    -------
    dict :
        Each value has keys `"n_photons"`, `"photon_e_max_MeV"`, and  `"photon_t_arr_at_e_max_s"`
        (the latter two are `None` when no photon is assigned to that pulse).
    """
    photon_pulse = [assign_pulse(t_arr=t, parameters=parameters, n_pulses=n_pulses) for t in lat_t_arr]

    photon_summary = {}
    for i in range(1, n_pulses + 1):
        idx = [k for k, pi in enumerate(photon_pulse) if pi == i]
        if not idx:
            photon_summary[i] = {"n_photons": 0, "photon_e_max_MeV": None, "photon_t_arr_at_e_max_s": None}
            continue
        best = max(idx, key=lambda k: photon_energy[k])
        photon_summary[i] = {
            "n_photons": len(idx),
            "photon_e_max_MeV": photon_energy[best],
            "photon_t_arr_at_e_max_s": lat_t_arr[best],
        }

    return photon_summary


def norris_plotter(
    fitter: NorrisFitter, y_max, nai_dat_files: str, photon_t_arr_s, photon_energy, t05, t95, plot_pad, save_name
):
    """Build and save the fitted light curve with a twin-axis LAT photon overlay.

    Parameters
    ----------
    fitter :
        The converged joint fit; `fitter.plot_fit()` draws the data, total fit, and individual pulses.
    y_max :
        Peak of the (normalized) count-rate data, for the left axis's tick alignment.
    nai_dat_files :
        Detector names summed into the light curve, for the data-series legend label.
    photon_t_arr_s :
        LAT photon arrival times, seconds since trigger.
    photon_energy :
        LAT photon energies in MeV, same order as `photon_t_arr_s`.
    t05 :
        T90 start, seconds since trigger -- plot-axis bounding only.
    t95 :
        T90 end, seconds since trigger -- plot-axis bounding only.
    plot_pad :
        Seconds of padding added on each side of `(t05, t95)` for the
        x-axis limits.
    save_name :
        Output file stem (no extension); saved as `<save_name>.png` and `.pdf` next to this module.
    """
    data_label = f"10-400 keV NaI\n({'+'.join(nai_dat_files)})"
    fig, ax = plt.subplots(figsize=(12, 7.5))
    fitter.plot_fit(
        show_individuals=True,
        x_label="Time since trigger [s]",
        y_label="Count rate [counts/s]",
        data_label=data_label,
        title=" ",
        axis=ax,
    )

    ax_photon = ax.twinx()
    ax_photon.scatter(
        photon_t_arr_s,
        photon_energy,
        marker="o",
        facecolors="none",
        edgecolors="red",
        linewidths=1.2,
        s=MARKER_SIZE**2,
        label="LAT photons (E > 1 GeV)",
    )
    ax_photon.set_ylabel("Photon energy [MeV]")
    photon_top = photon_energy.max()

    _axes_shenanigans(
        ax=ax, ax_photon=ax_photon, y_max=y_max, photon_top=photon_top, t05=t05, t95=t95, plot_pad=plot_pad
    )

    fig_path = Path(__file__).parent / save_name
    save_fig(fig_=fig, path=fig_path)


def _axes_shenanigans(
    ax: Axes, ax_photon: Axes, y_max: float, photon_top: float, t05: float, t95: float, plot_pad: float
):
    """Align the count-rate and photon-energy twin axes' zero, top, and gridlines.

    Parameters
    ----------
    ax :
        Left axis (count rate).
    ax_photon :
        Right, twinned axis (photon energy).
    y_max :
        Left axis's own data maximum.
    photon_top :
        Right axis's own data maximum.
    t05 :
        T90 start, seconds since trigger.
    t95 :
        T90 end, seconds since trigger.
    plot_pad :
        Seconds of x-axis padding on each side of `(t05, t95)`.
    """
    ax.set_ylim(-Y_BUFFER_FRAC * y_max, y_max * (1 + Y_BUFFER_FRAC))
    ax_photon.set_ylim(-Y_BUFFER_FRAC * photon_top, photon_top * (1 + Y_BUFFER_FRAC))

    ax.set_yticks(np.linspace(0, y_max, N_Y_TICKS))
    ax_photon.set_yticks(np.linspace(0, photon_top, N_Y_TICKS))

    ax.yaxis.set_minor_locator(AutoMinorLocator(MINOR_TICKS_PER_MAJOR))
    ax_photon.yaxis.set_minor_locator(AutoMinorLocator(MINOR_TICKS_PER_MAJOR))

    ax.set_xlim(t05 - plot_pad, t95 + plot_pad)
