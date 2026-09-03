"""Created on Dec 26 00:30:38 2025"""

from dataclasses import dataclass

import numpy as np
from numpy.linalg import LinAlgError, cholesky, inv

from .grb_constants import N_SAMPLES

FOUR_TABS = "\t\t\t\t"


# make the class as dataclass
@dataclass
class ParameterSet:
    """Container for a set of parameters."""

    parameters: list

    def __repr__(self):
        l = "ParameterSet[\n"
        for param in self.parameters:
            l += f"{FOUR_TABS}{param.name}," + "\n"
        l += f"{FOUR_TABS}]"
        return l
        # return f"ParameterSet[{', '.join(param.name for param in self.parameters)}\n]"

    def __getitem__(self, key):
        if isinstance(key, str):
            for param in self.parameters:
                if param.name == key:
                    return param
            raise KeyError(f"Parameter '{key}' not found in ParameterSet")
        else:
            return self.parameters[key]

    def get_populated_values(
        self,
        cov_matrix,
        parameter_name=None,
        size: int = N_SAMPLES,
        *,
        rng: np.random.Generator,
    ):
        """
        Returns a multivariate normal sample from the parameter set with a given covariance matrix.

        Parameters
        ----------
        cov_matrix :
            Covariance matrix for the parameters.
        parameter_name :
            Specific parameter name(s) to return.
            If None, returns all parameters.
        size :
            Number of samples to generate (default: 10_000).
        rng :
            Random number generator instance for reproducibility.

        Returns
        -------
        np.ndarray
            Multivariate normal samples.
        """
        cov_ = 0.5 * (cov_matrix + cov_matrix.T)
        values = [i.value for i in self.parameters]
        names = [i.name for i in self.parameters]

        rng_instance = rng

        if parameter_name is not None:
            if isinstance(parameter_name, str):
                parameter_name = [parameter_name]
            index = [names.index(i) for i in parameter_name]
            return rng_instance.multivariate_normal(values, cov_, size=size)[:, index]

        return rng_instance.multivariate_normal(values, cov_, size=size)


@dataclass(frozen=True)
class Parameter:
    """Scalar parameter with uncertainty."""

    name: str
    value: float
    error: float

    def __post_init__(self):
        if self.error < 0:
            raise ValueError("Parameter error must be non-negative")

    def __repr__(self):
        return f"Parameter[\n\tname={self.name}, value={self.value:.4g}, error={self.error:.4g}, rel_error={self.relative_error:.2%}\n]"

    @property
    def relative_error(self):
        """Returns the relative error (error / |value|)."""
        return self.error / abs(self.value) if self.value != 0 else np.inf

    @property
    def is_unconstrained(self):
        """Returns True if relative error >= 0.4."""
        return self.relative_error >= 0.4

    @property
    def is_moderately_constrained(self):
        """Returns True if 0.1 <= relative error < 0.4."""
        return 0.1 <= self.relative_error <= 0.4

    @property
    def is_well_constrained(self):
        """Returns True if relative error < 0.1."""
        return self.relative_error < 0.1

    def is_safe(self, threshold=0.4):
        """Returns True if relative error < threshold."""
        return self.relative_error < threshold


@dataclass
class CovarianceMatrix:
    """Class for handling covariance matrices."""

    input_matrix: np.ndarray

    def __post_init__(self):
        self.matrix = np.asarray(self.input_matrix, dtype=float)

        if self.matrix.ndim != 2:
            raise ValueError("Covariance matrix must be 2-dimensional")

        if self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError("Covariance matrix must be square")

    @property
    def variances(self):
        """Returns the variances (diagonal elements)."""
        return np.diag(self.matrix)

    @property
    def errors(self):
        """Returns standard errors (sqrt of diagonal)."""
        return np.sqrt(np.diag(self.matrix))

    def is_symmetric(self, tol=1e-10):
        """Checks if the matrix is symmetric."""
        return np.allclose(self.matrix, self.matrix.T, atol=tol)

    def is_positive_definite(self):
        """Checks if the matrix is positive definite."""
        try:
            cholesky(self.matrix)
            return True
        except LinAlgError:
            return False

    def corr(self, i, j):
        """Returns correlation coefficient between parameters i and j."""
        denominator = self.errors[i] * self.errors[j]
        if denominator == 0:
            raise ZeroDivisionError("Zero variance encountered")
        return self.matrix[i, j] / denominator

    def inverse(self):
        """Returns the inverse covariance (precision matrix)."""
        return inv(self.matrix)

    def __repr__(self):
        n = self.matrix.shape[0]
        return f"CovarianceMatrix({n}x{n})"
