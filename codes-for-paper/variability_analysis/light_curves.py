"""Copied verbatim from GRBResearchWork/light_curves/light_curves.py for standalone access
from this folder (per user request, 2026-08-22) -- not re-derived, keep in sync by hand if
the original changes."""

import numpy as np


def lightcurve_data(dat_file: str, energy_low: float, energy_high: float, errors=False):
    """
    Get the time, rate, and background from the RMFIT .dat file within a specific energy boundary.

    Parameters
    ----------
    dat_file: str
        RMFIT .dat file.
    energy_low: float
        The lower bound on energy.
    energy_high: float
        The upper bound on energy.
    errors: bool

    Returns
    -------
    time:
        Time measurements for the burst duration.
    data_rate:
        Rate of data accumulation within the energy bounds for the burst duration
    data_background:
        Estimated background for the energy bounds specified.
    errs:
        Errors on data_rate, and data_background (if errors=True)
    """

    data_rate_error = []
    data_background_error = []
    channel_number = 128

    # Read the entire file as a numpy array
    with open(dat_file) as f:
        lines = f.readlines()

    energy_channels = np.array([line.split() for line in lines[11 : 11 + channel_number]], dtype=float)
    low_energy_channel, high_energy_channel = energy_channels[:, 0], energy_channels[:, 1]

    detector_data = np.array([np.fromstring(line, sep=" ") for line in lines[12 + channel_number :]])

    # Determine which channels fall into the requested energy range
    start_idx: int = int(np.sum(low_energy_channel < energy_low))
    end_idx = int(np.sum(high_energy_channel < energy_high)) + 1

    # Initialize outputs
    data_rate = np.zeros(detector_data.shape[0])
    data_background = np.zeros(detector_data.shape[0])
    if errors:
        data_rate_error = np.zeros(detector_data.shape[0])
        data_background_error = np.zeros(detector_data.shape[0])

    # Sum over selected channels (vectorized)
    for x in range(start_idx, end_idx):
        data_rate += detector_data[:, 4 * x + 2]
        data_background += detector_data[:, 4 * x + 4]
        if errors:
            data_rate_error += detector_data[:, 4 * x + 3]
            data_background_error += detector_data[:, 4 * x + 5]  # fixed index

    time = detector_data[:, 0]

    if errors:
        return time, data_rate, data_background, (data_rate_error, data_background_error)
    else:
        return time, data_rate, data_background
