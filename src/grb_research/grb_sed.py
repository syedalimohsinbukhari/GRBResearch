"""Created on Dec 31 16:27:26 2025"""

from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np

from .grb_atomic import Parameter
from .grb_constants import kev_to_erg, model_n_pars
from .grb_enums import GRBModelsCombinations as gmC
from .grb_model import Model
from .grb_seds import band_function, black_body, cutoff_powerlaw, powerlaw, smoothly_broken_power_law

MODEL_MAP = {
    gmC.PL_BB: (gmC.PL, gmC.BB),
    # gmC.CPL_PL: (gmC.PL, gmC.CPL),
    gmC.CPL_BB: (gmC.CPL, gmC.BB),
    gmC.CPL_PL_BB: (gmC.PL, gmC.CPL, gmC.BB),
    # gmC.BAND_PL: (gmC.PL, gmC.BAND),
    gmC.BAND_BB: (gmC.BAND, gmC.BB),
    gmC.BAND_PL_BB: (gmC.PL, gmC.BAND, gmC.BB),
    # gmC.SBPL_PL: (gmC.PL, gmC.SBPL),
    gmC.SBPL_BB: (gmC.SBPL, gmC.BB),
    gmC.SBPL_PL_BB: (gmC.PL, gmC.SBPL, gmC.BB),
}


class GRBExceptions(Exception):
    """Base exception class for GRB-related errors."""


class ModelNotFoundInDataError(GRBExceptions):
    """Custom exception for handling missing models in SpectralModels."""

    def __init__(self):
        super().__init__(f"Model not found in the data.")


@dataclass
class SpectralModels:
    """Class to handle spectral models for GRB analysis."""

    model: Model
    model_type: str = "counts"

    n_sample: int = 10_000
    energy_range: Tuple[int, int] = (1, 7)
    n_grid: int = 10_000

    redshift: float = 0.0

    def __post_init__(self):
        self.interval = self.model.interval
        self.t_start = self.interval.start
        self.t_stop = self.interval.end

        # self.energy = np.logspace(*list(self.energy_range), self.n_points)

        self.SINGLE_COMPONENTS = {
            gmC.PL: (self._powerlaw, model_n_pars[gmC.PL]),
            gmC.BB: (self._black_body, model_n_pars[gmC.BB]),
            gmC.CPL: (self._cutoff_powerlaw, model_n_pars[gmC.CPL]),
            gmC.BAND: (self._band_grb_function, model_n_pars[gmC.BAND]),
            gmC.SBPL: (self._smoothly_broken_power_law, model_n_pars[gmC.SBPL]),
        }

    def __repr__(self):
        return (
            f"SpectralModels[\n"
            f"\tx=({self.energy_range[0]:g}, {self.energy_range[-1]}) Log[keV],\n"
            f"\t{self.model.name}, {self.interval}, '{self.model_type}', z={self.redshift:g}\n"
            f"\n]"
        )

    @classmethod
    def legacy_build(
        cls,
        m_name,
        interval_instance,
        p_name,
        p_vals,
        cov_,
        n_samples=10_000,
        n_grid=10_000,
        model_type="counts",
        e_range=(1, 7),
        redshift=0,
    ):
        """Build a SpectralModel from legacy data."""
        errors = np.sqrt(np.diag(cov_))
        p = [Parameter(i, j, k) for i, j, k in zip(p_name, p_vals, errors)]
        model = Model(m_name, p, interval_instance)

        return cls(model, model_type, n_sample=n_samples, n_grid=n_grid, energy_range=e_range, redshift=redshift)

    def _evaluate_components(self, components):
        values = [p.value for p in self.model.parameters]
        spectra = []

        idx = 0
        for name in components:
            func_name, n_pars = self.SINGLE_COMPONENTS[name]
            pars = values[idx: idx + n_pars]
            idx += n_pars

            spectra.append(func_name(pars))

        total = np.sum(spectra, axis=0)

        return spectra, total

    @property
    def _model_params(self):
        """Cache model parameters - only changes when the model changes."""
        return self.model.get_parameter_values()

    def _compute_model(self, spectral_func: Callable, joint_pars=None):
        st, ed = self.energy_range
        e_ = np.logspace(st, ed, self.n_grid)
        pars = joint_pars if joint_pars is not None else self._model_params

        spectrum = spectral_func(e_, *pars)

        # Dispatch based on model_type
        if self.model_type == "counts":
            return spectrum

        elif self.model_type == "energy":
            return e_ * spectrum

        elif self.model_type in ["nuFnu", "nfn"]:
            return e_ ** 2 * spectrum

        else:
            raise LookupError(
                f"Invalid model_type: {self.model_type}. Only 'counts', 'energy', 'nuFnu', and 'nfn' are supported."
            )

    def _powerlaw(self, joint_pars=None):
        return self._compute_model(spectral_func=powerlaw, joint_pars=joint_pars)

    def _black_body(self, joint_pars=None):
        return self._compute_model(spectral_func=black_body, joint_pars=joint_pars)

    def _cutoff_powerlaw(self, joint_pars=None):
        return self._compute_model(spectral_func=cutoff_powerlaw, joint_pars=joint_pars)

    def _smoothly_broken_power_law(self, joint_pars=None):
        return self._compute_model(spectral_func=smoothly_broken_power_law, joint_pars=joint_pars)

    def _band_grb_function(self, joint_pars=None):
        return self._compute_model(spectral_func=band_function, joint_pars=joint_pars)

    def _evaluate_model(self, model_key):
        """Evaluate a composite model based on the model key."""
        return self._evaluate_components(components=MODEL_MAP[model_key])

    def get_values(self, in_ergs=False):
        """Get the total spectral model values based on the model type."""
        convert = 1 if not in_ergs else kev_to_erg
        if self.model is None:
            raise ModelNotFoundInDataError()

        m_name = gmC(self.model.name.lower())

        singular_models = [gmC.PL, gmC.BB, gmC.CPL, gmC.BAND, gmC.SBPL]

        if m_name in singular_models:
            return self.SINGLE_COMPONENTS[m_name][0]() * convert
        else:
            seq = self._evaluate_model(model_key=m_name)
            return [i * convert if isinstance(i[0], float) else [j * convert for j in i] for i in seq]
