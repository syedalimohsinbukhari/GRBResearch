"""Created on Jan 07 15:37:00 2026"""

from __future__ import annotations

import hashlib
import os
import warnings
from multiprocessing import Pool, cpu_count
from typing import Optional, Tuple, Literal, Callable

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from matplotlib import pyplot as plt
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import simpson
from tqdm import tqdm

from .grb_constants import kev_to_erg, MASTER_SEED, model_n_pars, N_SAMPLES, N_GRID, FIGURE_SIZE_4x4
from .grb_enums import GRBModelsCombinations as gmC
from .grb_fits_io import build_composite_schema
from .grb_model import Model
from .grb_sed import MODEL_MAP, SpectralModels
from .grb_seds import band_function, black_body, cutoff_powerlaw, powerlaw, smoothly_broken_power_law
from .grb_time import EpisodeTypes
from .grb_utils import save_fig


def get_rng(seed: int | None = None, rng: np.random.Generator | None = None) -> np.random.Generator:
    """
    Get or create a NumPy random number generator.

    Parameters
    ----------
    seed : int, optional
        Seed for creating a new RNG. Ignored if rng is provided.
    rng : np.random.Generator, optional
        Existing RNG instance to use.

    Returns
    -------
    np.random.Generator
        RNG instance to use for random sampling.
    """
    if seed is None and rng is None:
        raise ValueError("Either seed or rng must be provided.")

    if rng is not None:
        return rng
    return np.random.default_rng(seed)


def seed_from_name(name: str, master_seed: int = MASTER_SEED) -> int:
    """Derive a deterministic 32-bit seed from a file name and a master seed.

    Notes
    -----
    Hashes the basename of `name` together with `master_seed` so different scripts get different, reproducible seeds
    without any script hardcoding a literal seed value.
    Uses `os.path.basename` rather than full path so the derived seed is stable across machines and checkout locations.

    Parameters
    ----------
    name :
        Path or file name to derive the seed from.
        Typically, a script's `__file__`.
    master_seed :
        Project-wide master seed mixed into the hash (default: `MASTER_SEED`).

    Returns
    -------
    int :
        A deterministic seed in [0, 2**32).
    """
    digest = hashlib.sha256(f"{os.path.basename(name)}-{master_seed}".encode()).hexdigest()
    return int(digest, 16) % (2 ** 32)


def legacy_build_mp(pars):
    """
    Multiprocessing worker wrapper for `SpectralModels.legacy_build`.

    Parameters
    ----------
    pars :
        - m_name : Model name.
        - interval : Model interval object (opaque to this function).
        - m_keys : List of parameter names.
        - sample : Parameter values for this sample.
        - covar : Covariance matrix.
        - model_type : String passed through to `legacy_build`.
        - e_range : Energy range for the model.
        - n_sample : Number of samples.
        - n_grid : Number of grid points.

    Returns
    -------
    NDArray
        The evaluated model values.
    """
    m_name, interval, m_keys, sample, covar, model_type, e_range, n_sample, n_grid = pars
    built = SpectralModels.legacy_build(
        m_name,
        interval,
        m_keys,
        sample,
        covar,
        n_samples=n_sample,
        n_grid=n_grid,
        model_type=model_type,
        e_range=e_range,
    ).get_values()

    # the original behavior returned element 1 when the name contains an underscore
    if "_" in m_name:
        return built[1]
    return built


def mc_spectra_sampler(
    model: Model,
    model_type="counts",
    e_range=(1, 7),
    n_samples: int = N_SAMPLES,
    n_grid: int = N_GRID,
    n_workers: int | None = None,
    samples: ArrayLike | None = None,
    *,
    rng: np.random.Generator,
):
    """
    Parallel Monte-Carlo sampler for spectral model parameter estimation.

    Parameters
    ----------
    model :
        The model to sample from.
    model_type :
        Type of model to generate (default: 'counts').
    e_range :
        Energy range for the model (default: (1, 7)).
    n_samples :
        Number of MC samples to draw (default: 10,000).
    n_grid :
        Number of grid points for numerical integration (default: 10,000).
    n_workers :
        Number of parallel workers (default: CPU count).
    samples :
        Pre-generated samples. If None, samples will be generated.
    rng :
        Random number generator instance for reproducibility.

    Returns
    -------
    list
        A list of evaluated model values for each sample.
    """
    m_keys = [i.name for i in model.parameters]
    m_vals = [i.value for i in model.parameters]

    covar_ = model.covariance_matrix_value
    covar_ = 0.5 * (covar_ + covar_.T)

    if samples is None:
        rng_instance = rng
        samples = rng_instance.multivariate_normal(mean=m_vals, cov=covar_, size=n_samples)
        m_res = ModelResampler(model=model, samples=samples, rng=rng_instance, destroy=True)
        samples = m_res.run_resampler()

    if n_workers is None:
        n_workers = cpu_count()

    args_list = [
        (model.name, model.interval, m_keys, v, covar_, model_type, e_range, n_samples, n_grid) for v in samples
    ]

    with Pool(n_workers) as pool:
        results = list(tqdm(iterable=pool.imap(func=legacy_build_mp, iterable=args_list), total=n_samples))

    return results


class ModelResampler:

    def __init__(self, model: Model, samples: NDArray, *, rng: np.random.Generator, destroy: bool = True):
        self.model = model
        self._samples = samples if destroy else samples.copy()
        self.rng: np.random.Generator = rng

        self.m_val = [i.value for i in model.parameters]
        self.errs = np.sqrt(np.diag(model.covariance_matrix_value))
        self.err_ratio = [(j / abs(i)) * 100 for i, j in zip(self.m_val, self.errs)]

    @staticmethod
    def _cond_check(schema=None) -> Tuple[bool, NDArray, NDArray]:
        pos_mask = np.array([p[-1] for p in schema], dtype=bool)
        neg_mask = ~pos_mask
        return True, pos_mask, neg_mask

    def _resampler(
        self,
        samples: NDArray,
        pos_mask: NDArray,
        neg_mask: NDArray,
        extra_mask_fn: Optional[Callable[[NDArray], NDArray]] = None,
    ) -> NDArray:
        """Resample invalid rows until all masked constraints are satisfied.

        Parameters
        ----------
        samples :
            Working copy of MC parameter draws, shape (n_samples, n_params).
        pos_mask :
            Boolean column mask; flagged parameters must be strictly positive.
        neg_mask :
            Boolean column mask; flagged parameters must be strictly negative.
        extra_mask_fn :
            Optional callable `(samples) -> bool array of shape (n_samples,)`
            returning `True` for rows that are *valid* under the model-specific
            constraint (e.g. Kaneko condition).  Recomputed from the current
            `samples` on every iteration so stale row assignments do not
            prevent convergence.  Pass `None` when there is no extra constraint.
        """
        max_rounds = 100
        for _ in range(max_rounds):
            mask = np.any(samples[:, pos_mask] < 0, axis=1) | np.any(samples[:, neg_mask] > 0, axis=1)
            if extra_mask_fn is not None:
                # recompute from current samples each round — avoids stale mask bug
                mask |= ~extra_mask_fn(samples)

            n_invalid = int(np.sum(mask))
            if n_invalid == 0:
                return samples

            print(f"The number of resampled parameters: {n_invalid}")
            samples[mask] = self.rng.multivariate_normal(self.m_val, self.model.covariance_matrix_value, n_invalid)

        warnings.warn(
            f"Reached max resampling rounds ({max_rounds}) for {self.model.name}; " "some invalid samples may remain."
        )
        return samples

    def __runner(
        self, samples: NDArray, extra_mask_fn: Callable[[NDArray], NDArray] | None = None
    ) -> NDArray:
        schema = build_composite_schema(self.model.name)
        check, pos_mask, neg_mask = self._cond_check(schema)
        if check:
            samples = self._resampler(samples, pos_mask, neg_mask, extra_mask_fn=extra_mask_fn)
        return samples

    def _pl_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _cpl_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _band_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _sbpl_resampler(self, samples: NDArray) -> NDArray:
        def _sbpl_valid(s: NDArray) -> NDArray:
            l1, l2 = s[:, 2], s[:, 5]
            return np.logical_and(l1 > -2.0, l2 < -2.05)

        return self.__runner(samples, extra_mask_fn=_sbpl_valid)

    def _pl_bb_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _cpl_bb_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _band_bb_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _sbpl_bb_resampler(self, samples: NDArray) -> NDArray:
        return self._sbpl_resampler(samples)

    def _cpl_pl_bb_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _band_pl_bb_resampler(self, samples: NDArray) -> NDArray:
        return self.__runner(samples)

    def _sbpl_pl_bb_resampler(self, samples: NDArray) -> NDArray:
        # amp_pl, index1_pl, e_piv_pl
        # amp_sbpl, e_piv_sbpl, index1_sbpl, e_break_sbpl, break_scale_sbpl, index2_sbpl
        # amp_bb, kT_bb
        def _sbpl_pl_bb_valid(s: NDArray) -> NDArray:
            # rows where the SBPL physical condition is NOT satisfied are invalid
            return ~np.logical_and(s[:, 5] > -2, s[:, 8] < -2.05)

        return self.__runner(samples, extra_mask_fn=_sbpl_pl_bb_valid)

    def run_resampler(self) -> NDArray:
        """Return a corrected copy of sampled parameters after model-specific resampling."""
        dispatcher: dict[str, Callable[[NDArray], NDArray]] = {
            gmC.PL.name_upper: self._pl_resampler,
            gmC.PL_BB.name_upper: self._pl_bb_resampler,
            gmC.CPL.name_upper: self._cpl_resampler,
            gmC.CPL_BB.name_upper: self._cpl_bb_resampler,
            gmC.CPL_PL_BB.name_upper: self._cpl_pl_bb_resampler,
            gmC.BAND.name_upper: self._band_resampler,
            gmC.BAND_BB.name_upper: self._band_bb_resampler,
            gmC.BAND_PL_BB.name_upper: self._band_pl_bb_resampler,
            gmC.SBPL.name_upper: self._sbpl_resampler,
            gmC.SBPL_BB.name_upper: self._sbpl_bb_resampler,
            gmC.SBPL_PL_BB.name_upper: self._sbpl_pl_bb_resampler,
        }

        resampler = dispatcher.get(self.model.name, None)
        if resampler is not None:
            return resampler(self._samples.copy())
        else:
            print(f"Warning: No resampler for {self.model.name} model.")
            return self._samples.copy()


def credible_interval_partition(samples: NDArray) -> Tuple[NDArray, NDArray, NDArray]:
    """
    Compute 16th, 50th (median), and 84th percentiles per parameter from MC samples.

    Parameters
    ----------
    samples :
        2D input array with shape (n_samples, n_parameters).
        Each row is an independent sample of the parameter vector.

    Returns
    -------
    tuple
        A tuple (median, lower, upper), each a 1D array of shape (n_parameters,) containing the 50th, 16th, and 84th
        percentiles respectively.
    """
    s = samples.T
    part = np.nanpercentile(s, [16, 50, 84], axis=1)

    return np.asarray(part[1], dtype=float).T, np.asarray(part[0], dtype=float).T, np.asarray(part[2], dtype=float).T


def mc_e_iso_sampler(
    model: Model,
    z: float = 1.0,
    n_samples: int = N_SAMPLES,
    n_grid: int = N_GRID,
    det_min: float = 1.0,
    det_max: float = 7.0,
    bol_min: float = 0.0,
    bol_max: float = 4.0,
    h0: float = 69.6,
    omega_m: float = 0.286,
    method: int = 1,
    samples=None,
    *,
    rng: np.random.Generator,
) -> float:
    """
    Draw MC samples and compute isotropic-equivalent energy (E_iso).

    Parameters
    ----------
    model :
        Spectral model container providing sampling and interval duration.
    z :
        Redshift used for K-correction and luminosity distance (default: 1.0).
    n_samples :
        Number of MC samples to draw (default: 100).
    n_grid :
        Number of energy grid points for numerical integration (default: 100).
    det_min :
        Log10 lower bound for detector energy grid (keV) (default: 1.0).
    det_max :
        Log10 upper bound for detector energy grid (keV) (default: 7.0).
    bol_min :
        Log10 lower bound for bolometric energy grid (keV) (default: -1.0).
    bol_max :
        Log10 upper bound for bolometric energy grid (keV) (default: 4.0).
    h0 :
        Hubble constant (default: 70.0).
    omega_m :
        Matter density parameter (default: 0.315).
    method :
        Method to use for calculation (default: 1).
    samples :
        Pre-generated samples. If None, samples will be generated.
    rng :
        Random number generator instance for reproducibility.

    Returns
    -------
    NDArray
        Array of E_iso samples in erg with shape (1, n_samples).
    """
    rng_instance = rng
    bolometric_fluence = 0

    redshift_shift = np.log10(1 + z)
    bol_range_observed = (bol_min - redshift_shift, bol_max - redshift_shift)
    e_observed = np.logspace(start=bol_range_observed[0], stop=bol_range_observed[1], num=n_grid)

    bolometric_samples = np.asarray(
        mc_spectra_sampler(
            model=model,
            model_type="energy",
            e_range=bol_range_observed,
            n_samples=n_samples,
            n_grid=n_grid,
            samples=samples,
            rng=rng_instance,
        )
    )
    bolometric_flux = simpson(y=bolometric_samples, x=e_observed, axis=1)

    if method == 1:
        energy_detector = np.logspace(start=det_min, stop=det_max, num=n_grid)
        detector_samples = np.asarray(
            mc_spectra_sampler(
                model=model,
                model_type="energy",
                e_range=(det_min, det_max),
                n_samples=n_samples,
                n_grid=n_grid,
                rng=rng_instance,
            )
        )
        detector_flux = simpson(y=detector_samples, x=energy_detector, axis=1)

        s_obs = detector_flux * model.interval.duration
        k_correction = bolometric_flux / detector_flux
        bolometric_fluence = np.asarray(s_obs * k_correction, dtype=float) * kev_to_erg
    elif method == 2:
        bolometric_fluence = bolometric_flux * kev_to_erg * model.interval.duration

    lum_distance = FlatLambdaCDM(h0, omega_m).luminosity_distance(z).cgs.value

    return 4 * np.pi * lum_distance ** 2 * np.asarray(bolometric_fluence).reshape(1, -1) / (1 + z)


def plot_all_models(
    best_models,
    grb_name,
    n_rows: int = 2,
    n_cols: int | None = None,
    n_grid: int = N_GRID,
    n_samples: int = N_SAMPLES,
    fig_size: tuple[float, float] = FIGURE_SIZE_4x4,
    save: bool = False,
    *,
    rng: np.random.Generator,
):
    x = np.logspace(1, 7, n_grid)

    legend_col = {0: 3, 1: 1, 2: 3, 3: 1}

    figure, ax = plt.subplots(n_rows, n_cols, figsize=fig_size, sharey=True, sharex=True)
    ax = ax.flatten()

    for i, v in enumerate(best_models):
        print(f"processing {grb_name[i]}")
        is_ex = sum([ep.interval.is_ex for ep in v])
        if is_ex == 2:
            v[-1], v[-2] = v[-2], v[-1]

        for j, w in enumerate(v):
            print(f"processing {w.name}: {w.interval}")
            samples = mc_spectra_sampler(w, n_samples=n_samples, n_grid=n_grid, rng=rng)
            samples = np.array(samples)

            med, low, high = credible_interval_partition(samples)
            med, low, high = med * kev_to_erg, low * kev_to_erg, high * kev_to_erg

            if j == 0:
                ax[i].loglog(
                    x,
                    med * x ** 2,
                    "k-",
                    zorder=1000,
                    label=f"{w.interval.kind}" + r"$_\text{" + f'{w.name.replace("_", "+")}' + r"}$",
                )
                ax[i].fill_between(x, low * x ** 2, high * x ** 2, zorder=1000, color="k", alpha=0.2)
            else:
                sub = (
                    f"{w.interval.kind}{w.interval.index}"
                    if w.interval.kind in [EpisodeTypes.TR, EpisodeTypes.SP]
                    else w.interval.kind
                )
                ax[i].loglog(x, med * x ** 2, "--",
                             label=f"{sub}" + r"$_\text{" + f'{w.name.replace("_", "+")}' + r"}$")
                ax[i].fill_between(x, low * x ** 2, high * x ** 2, alpha=0.2)

            ax[i].set_ylim(bottom=3.2e-10, top=1.65e-4)

        ax[i].legend(ncols=legend_col.get(i, 3), title=f"GRB{grb_name[i]}", loc="upper right")

        if i % n_cols == 0:  # fixed: was hardcoded % 2, now uses n_cols
            ax[i].set_ylabel("Energy Flux\n" + r"[erg/cm$^2$/s]")

    [
        ax[i].set_xlabel("Energy [keV]") for i in range(len(best_models) - n_cols, len(best_models))
    ]

    if save:
        save_fig(figure, "butterfly_all")
    else:
        plt.show()


def relative_error(
    value_true: float,
    value_approx: float,
    absolute: bool = True,
    as_percent: bool = False,
    zero_handling: Literal["ignore", "inf", "raise"] = "ignore",
):
    """
    Calculate the relative error between a true/reference value and an approximation.

    Parameters
    ----------
    value_true :
        The true or reference value.
    value_approx :
        The approximate or measured value.
    absolute :
        If True, returns the absolute relative error (always non-negative).
        If False, returns the signed error (positive if approximation > true).
        Default is True.
    as_percent :
        If True, multiply the result by 100 to return a percentage.
        Default is False.
    zero_handling :
        How to handle the case when value_true is zero.
        Options:
            - 'ignore' : return NaN (default)
            - 'inf' : return inf (if absolute) or with sign depending on value_approx
            - 'raise' : raise ZeroDivisionError

    Returns
    -------
    float
        The relative error. May be NaN, inf, or finite depending on inputs and options.
    """
    if value_true == 0:
        if zero_handling == "ignore":
            return float("nan")
        elif zero_handling == "inf":
            # Signed infinite: sign(value_approx) * infinity; if absolute, just inf
            if absolute:
                return float("inf")
            else:
                return float("inf") if value_approx > 0 else -float("inf")
        elif zero_handling == "raise":
            raise ZeroDivisionError("Cannot compute relative error with value_true = 0.")
        else:
            raise ValueError("zero_handling must be 'ignore', 'inf', or 'raise'.")

    if absolute:
        err = abs(value_approx - value_true) / abs(value_true)
    else:
        err = (value_approx - value_true) / value_true

    if as_percent:
        err *= 100.0

    return err


class FluxFluenceCalculator:
    """Calculates flux and fluence based on a spectral model using Monte Carlo sampling.

    Attributes
    ----------
    spectral_model :
        The spectral model used for flux and fluence calculations.
    log_energy_range :
        The logarithmic bounds of the energy range over which calculations are performed.
    n_samples :
        The number of Monte Carlo samples to generate.
    n_grid :
        The number of bins in the logarithmic energy grid.
    rng :
        The random number generator used for Monte Carlo sampling.
    """

    def __init__(
        self,
        spectral_model: Model,
        log_energy_range: tuple[float, float] = (1, 3),
        n_samples: int = N_SAMPLES,
        n_grid: int = N_GRID,
        *,
        rng: np.random.Generator,
    ):
        self.spectral_model = spectral_model
        self.log_energy_range = log_energy_range
        self.n_samples = n_samples
        self.n_grid = n_grid

        self.rng = rng

    def _flux(self) -> NDArray:
        """Generate flux values based on a given spectral model within a specified energy range.

        Returns
        -------
        NDArray
            An array containing the integrated flux values corresponding to the specified energy grid.
            Each element represents the calculated flux for the associated energy range.
        """
        x = np.logspace(*self.log_energy_range, self.n_grid)
        n_of_e = mc_spectra_sampler(
            self.spectral_model,
            "counts",
            e_range=self.log_energy_range,
            n_samples=self.n_samples,
            n_grid=self.n_grid,
            rng=self.rng,
        )
        return np.asarray(simpson(np.array(n_of_e), x))

    def _fluence(self, in_ergs: bool = False, energy_flux: bool = False) -> NDArray:
        """Calculates the fluence over a specified energy range using Monte Carlo sampling and numerical integration.

        Parameters
        ----------
        in_ergs :
            If True, convert the fluence results to ergs using an energy conversion factor.
            If False, results will use the default unit (keV).
            Default is False.

        Returns
        -------
        NDArray
            The computed fluence over the specified energy range.
        """
        converter = kev_to_erg if in_ergs else 1
        duration = 1 if energy_flux else self.spectral_model.interval.duration
        x = np.logspace(*self.log_energy_range, self.n_grid)
        n_of_e = mc_spectra_sampler(
            self.spectral_model,
            "energy",
            e_range=self.log_energy_range,
            n_samples=self.n_samples,
            n_grid=self.n_grid,
            rng=self.rng,
        )
        return np.asarray(simpson(np.array(n_of_e), x) * converter) * duration

    def calculate(
        self,
        calculation_type: Literal["flux", "fluence"] = "flux",
        get_percentiles: bool = False,
        in_ergs: bool = True,
        energy_flux: bool = False,
        get_errors: bool = True,
    ) -> NDArray | tuple[NDArray, NDArray, NDArray]:
        """Performs a calculation based on the specified type.

        Parameters
        ----------
        calculation_type :
            The type of calculation to perform. Defaults to 'flux'.
        get_percentiles :
            If True, returns the 16th, 50th, and 84th percentiles of the calculated output.
            Defaults to False.
        in_ergs :
            If True, the fluence will be returned in ergs instead of keV.
            This parameter is only relevant when `calculation_type` is set to 'fluence'.
            Defaults to True.
        get_errors :
            If True, returns the median value along with the upper and lower error margins.
            If provided, it overrides `get_percentiles`.
            Defaults to True.

        Returns
        -------
        NDArray :
            The returned value depends on the parameters:

            - If `get_percentiles`: Returns a numpy array containing the 16th, 50th, and 84th percentiles of the output.
            - If `get_errors`: Returns a tuple consisting of the median value, upper margin, and lower margin.
            - Otherwise, returns a numpy array of the calculated results for 'flux' or 'fluence'.
        """

        if get_percentiles and get_errors:
            warnings.warn("Cannot return both percentiles and errors. Using `get_errors`")
            get_percentiles = False

        if calculation_type == "flux":
            output = self._flux()
        elif calculation_type == "fluence":
            output = self._fluence(in_ergs, energy_flux=energy_flux)
        else:
            raise ValueError("Invalid calculation type. Must be 'flux' or 'fluence'.")

        if get_percentiles:
            return np.percentile(output, [16, 50, 84])
        if get_errors:
            percentiles = np.percentile(output, [16, 50, 84])
            return percentiles[1], percentiles[2] - percentiles[1], percentiles[1] - percentiles[0]

        return output


_SED_FUNCTIONS = {
    gmC.PL: powerlaw,
    gmC.CPL: cutoff_powerlaw,
    gmC.BAND: band_function,
    gmC.SBPL: smoothly_broken_power_law,
    gmC.BB: black_body,
}


def _component_energy_flux_block(model_name, values, energy):
    """Vectorized per-component energy-flux integration for a block of draws."""
    components = MODEL_MAP.get(gmC(model_name.lower()), (gmC(model_name.lower()),))
    energy = np.asarray(energy)
    batched = energy.ndim == 2

    if batched:
        energy_row = energy[None, :, :]
        x_for_simpson = np.broadcast_to(energy, (values.shape[0],) + energy.shape)
        integration_axis = -1
        total = np.zeros((values.shape[0],) + energy.shape)
    else:
        energy_row = energy[None, :]
        x_for_simpson = energy
        integration_axis = 1
        total = np.zeros((values.shape[0], energy.size))

    fluxes = {}
    idx = 0
    for component in components:
        n_pars = model_n_pars[component]
        if batched:
            pars = [values[:, idx + offset, None, None] for offset in range(n_pars)]
        else:
            pars = [values[:, idx + offset, None] for offset in range(n_pars)]
        idx += n_pars

        contribution = energy_row * _SED_FUNCTIONS[component](energy_row, *pars)
        total += contribution
        fluxes[component] = simpson(y=contribution, x=x_for_simpson, axis=integration_axis)

    return fluxes, simpson(y=total, x=x_for_simpson, axis=integration_axis)


def component_energy_fluxes(model_name, values, energy, chunk: int = 2_000):
    """Integrate each spectral component's energy flux, and their total.

    Notes
    -----
    Components are sliced in declaration order, matching `SpectralModels._evaluate_components`.
    All draws are evaluated at once, with parameters of shape `(n, 1)` broadcast against energy of shape `(1, n_grid)`,
    then integrated via a single `simpson(..., axis=1)` (same pattern as :func:`mc_e_iso_sampler`).

    Draws are processed in blocks to bound peak memory at roughly `chunk * n_grid * 8` bytes per component.
    For 2D `energy` (`n_bands` batched grids), `chunk` is divided by `n_bands` internally, since a block's working
    array is `(block_size, n_bands, n_grid)` -- without this, peak memory would scale with `n_bands`.

    Parameters
    ----------
    model_name :
        Composite model name, e.g. `"SBPL_BB"`.
    values :
        Parameter values, shape `(n_pars,)` or `(n, n_pars)`.
    energy :
        Energy grid in keV.
        Either 1D `(n_grid,)` or 2D `(n_bands, n_grid)` to integrate several bands in one batched pass;
        see :func:`_component_energy_flux_block`.
    chunk :
        Number of draws evaluated per block, at one band.
        Automatically reduced (never below 1) when `energy` is 2D, so peak memory per block is independent of `n_bands`.

    Returns
    -------
    fluxes :
        Maps each component enum to an array of shape `(n, )` or `(n, n_bands)` if `energy` is 2D.
    total :
        The total combined flux of the components.
        Has the same shape as `fluxes`.
    """
    values = np.atleast_2d(np.asarray(values, dtype=float))
    energy = np.asarray(energy)
    if energy.ndim == 2:
        chunk = max(1, chunk // energy.shape[0])

    parts, totals = [], []
    for start in range(0, values.shape[0], chunk):
        block_fluxes, block_total = _component_energy_flux_block(model_name, values[start: start + chunk], energy)
        parts.append(block_fluxes)
        totals.append(block_total)

    merged = {component: np.concatenate([p[component] for p in parts]) for component in parts[0]}
    return merged, np.concatenate(totals)


def draw_model_samples(model, n_samples: int = 10_000, *, rng: np.random.Generator):
    """Draw fit-covariance parameter samples, filtered by physical constraints.

    Notes
    -----
    Mirrors the sampling half of :func:`mc_spectra_sampler` but returns the parameter draws instead of evaluated
    spectra, so callers can compute their own derived quantities from correlated samples.
    """
    rng_instance = rng

    covariance = model.covariance_matrix_value
    covariance = 0.5 * (covariance + covariance.T)
    means = [p.value for p in model.parameters]

    samples = rng_instance.multivariate_normal(mean=means, cov=covariance, size=n_samples)
    return ModelResampler(model=model, samples=samples, rng=rng_instance, destroy=True).run_resampler()
