"""Created on Dec 26 14:20:28 2025"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .grb_atomic import CovarianceMatrix, Parameter, ParameterSet
from .grb_enums import GoodnessOfFit

if TYPE_CHECKING:
    from .grb_time import TimeInterval


@dataclass
class Model:
    """Represents a GRB model with its parameters and fit statistics."""

    name: str
    parameters: list[Parameter]
    interval: TimeInterval | None = None

    status: GoodnessOfFit = GoodnessOfFit.UNKNOWN
    cstat: float = 0.0
    dof: int = 0
    covariance_matrix: CovarianceMatrix | None = None

    def get_parameter_value(self, par_name):
        """Get parameter value by name."""
        for p in self.parameters:
            if p.name == par_name:
                return p.value
        return None

    def get_parameter_values(self, get_errors=False, get_both=False):
        """Get parameter values as a numpy array."""
        if get_both:
            return np.array([[v.value, v.error] for v in self.parameters])
        if get_errors:
            return np.array([v.error for v in self.parameters])
        return np.array([v.value for v in self.parameters])

    @property
    def is_best(self):
        """Check if the model is the best fit."""
        return self.status is GoodnessOfFit.BEST

    @property
    def is_good(self):
        """Check if the model is a good fit."""
        return self.status is GoodnessOfFit.GOOD

    @property
    def is_marginal(self):
        """Check if the model is a marginal fit."""
        return self.status is GoodnessOfFit.MARGINAL

    @property
    def is_unsafe(self):
        """Check if the model is unsafe."""
        return self.status is GoodnessOfFit.UNSAFE

    @property
    def get_reduced_cstat(self):
        """Get the reduced c-statistic (cstat/dof)."""
        if self.dof == 0:
            return np.inf
        return self.cstat / self.dof

    @property
    def get_parameter_set(self):
        """Retrieve the parameter set associated with current parameters.

        This property provides access to a ParameterSet object representing the current parameters of the instance.

        Returns:
            ParameterSet: An object encapsulating the current parameters.
        """
        return ParameterSet(self.parameters)

    @property
    def covariance_matrix_value(self):
        """Get the covariance matrix."""
        return 0.5 * (self.covariance_matrix.matrix + self.covariance_matrix.matrix.T)

    @classmethod
    def from_dictionary(cls, name: str, data: dict, interval: TimeInterval) -> "Model":
        """Create a Model from its dictionary representation."""
        internal_dict = copy.deepcopy(data)

        status = GoodnessOfFit(internal_dict["_status"])
        cstat = internal_dict["c-stat/dof"][0]
        dof = int(internal_dict["c-stat/dof"][1])
        cov_matrix = np.array(internal_dict["covariance_matrix"])
        aux_keys1 = ["_status", "c-stat/dof", "covariance_matrix"]
        for aux_ in aux_keys1:
            internal_dict.pop(aux_)

        return cls(
            name,
            [Parameter(k, v, e) for k, (v, e) in internal_dict.items()],
            interval,
            status,
            cstat,
            dof,
            CovarianceMatrix(cov_matrix),
        )

    def __str__(self) -> str:
        if self.is_best:
            safety = GoodnessOfFit.BEST
        elif self.is_good:
            safety = GoodnessOfFit.GOOD
        elif self.is_marginal:
            safety = GoodnessOfFit.MARGINAL
        else:
            safety = GoodnessOfFit.UNSAFE

        return (
            f"model: {self.name},\n"
            f"_status: {safety}\n"
            f"n_parameters: {len(self.parameters)},\n"
            f"cstat/dof: {self.cstat:.4f}/{self.dof},\n"
            f"covariance_matrix: {self.covariance_matrix}"
        )

    def __repr__(self) -> str:
        if self.is_best:
            safety = GoodnessOfFit.BEST
        elif self.is_good:
            safety = GoodnessOfFit.GOOD
        elif self.is_marginal:
            safety = GoodnessOfFit.MARGINAL
        else:
            safety = GoodnessOfFit.UNSAFE
        params = ParameterSet(self.parameters)
        # params = ",\n        ".join(repr(p) for p in self.parameters)

        return (
            "Model[\n"
            f"    name={self.name!r},\n"
            f"    _status: {safety}\n"
            f"    parameters=(\n"
            f"        {params}\n"
            f"    ),\n"
            f"    cstat={self.cstat},\n"
            f"    dof={self.dof},\n"
            f"    covariance_matrix={self.covariance_matrix}\n"
            "]"
        )


@dataclass
class ModelSet:
    """A container for GRB spectral models."""

    _models: list[Model]

    def __post_init__(self):
        self._by_name: dict[str, Model] = {m.name: m for m in self._models}

    def __repr__(self) -> str:
        if not self._models:
            return "ModelSet(empty)"

        lines = ["ModelSet("]
        for m in self._models:
            int_ = m.interval.to_string().split(" ")[0]
            lines.append(
                f"\tModel({m.name} ({int_}), status={m.status.value}, " f"cstat/dof={m.cstat:.3f}/{m.dof}),"
            )
        lines.append(")")
        return "\n".join(lines)

    def __getitem__(self, key: int | str | slice) -> Model | ModelSet:
        if isinstance(key, int):
            return self._models[key]
        elif isinstance(key, slice):
            out_ = self._models[key]
            return ModelSet(out_)
        return self._by_name[key]

    def __setitem__(self, key, value):
        self._models[key] = value

    def __iter__(self):
        return iter(self._models)

    def __len__(self):
        return len(self._models)

    @property
    def best(self) -> Model:
        """Return the BEST model."""
        for m in self._models:
            if m.status is GoodnessOfFit.BEST:
                return m
        raise LookupError("No model with status BEST found.")

    @property
    def safe(self) -> "ModelSet":
        """Return all SAFE models."""
        safe_models: list[Model] = [m for m in self._models if m.status is GoodnessOfFit.SAFE]
        return ModelSet(safe_models)

    @property
    def good(self) -> "ModelSet":
        """Return all GOOD models."""
        good_models: list[Model] = [
            m for m in self._models if m.status not in [GoodnessOfFit.UNSAFE, GoodnessOfFit.UNSAFE]
        ]
        return ModelSet(good_models)

    @property
    def unsafe(self) -> "ModelSet":
        """Return all UNSAFE models."""
        unsafe_models: list[Model] = [m for m in self._models if m.status is GoodnessOfFit.UNSAFE]
        return ModelSet(unsafe_models)

    @property
    def names(self) -> tuple[str, ...]:
        """Get the names of all models in the set."""
        return tuple(self._by_name.keys())

    def get(self, name: str):
        """Get a model by name, or None if not found."""
        return self._by_name.get(name.upper())
