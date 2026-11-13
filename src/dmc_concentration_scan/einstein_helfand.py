"""Collective electrical conductivity from the Einstein-Helfand charge moment."""

from __future__ import annotations

import numpy as np


ELEMENTARY_CHARGE_C = 1.602176634e-19
BOLTZMANN_J_K = 1.380649e-23


def charge_moment_msd(moment_eA: np.ndarray) -> np.ndarray:
    """Return time-origin averaged |M(t+tau)-M(t)|² in e² Å² using FFT O(N log N)."""
    moment = np.asarray(moment_eA, dtype=float)
    if moment.ndim != 2 or moment.shape[1] != 3 or moment.shape[0] < 2:
        raise ValueError("moment must have shape (n_frames>=2, 3)")
    if not np.isfinite(moment).all():
        raise ValueError("charge moment contains non-finite values")
    n_frames = len(moment)
    n_fft = 1 << (2 * n_frames - 1).bit_length()
    spectrum = np.fft.rfft(moment, n=n_fft, axis=0)
    correlation = np.fft.irfft(np.sum(spectrum * np.conjugate(spectrum), axis=1), n=n_fft)[:n_frames]
    squared = np.einsum("ij,ij->i", moment, moment)
    prefix = np.concatenate(([0.0], np.cumsum(squared)))
    lag = np.arange(n_frames)
    origins = n_frames - lag
    left = prefix[n_frames - lag]
    right = prefix[n_frames] - prefix[lag]
    msd = (left + right - 2.0 * correlation) / origins
    msd[0] = 0.0
    return np.maximum(msd, 0.0)


def conductivity_from_slope(slope_e2_A2_ps: float, volume_A3: float, temperature_K: float) -> float:
    """Convert EH-MSD slope to S/m (collective ionic conductivity)."""
    if slope_e2_A2_ps <= 0 or volume_A3 <= 0 or temperature_K <= 0:
        raise ValueError("slope, volume, and temperature must be positive")
    numerator = slope_e2_A2_ps * ELEMENTARY_CHARGE_C**2 * 1.0e-20 / 1.0e-12
    denominator = 6.0 * volume_A3 * 1.0e-30 * BOLTZMANN_J_K * temperature_K
    return float(numerator / denominator)


def linear_fit(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
        raise ValueError("linear fit requires aligned 1D arrays with at least three points")
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    denom = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(residual**2)) / denom if denom > 0 else np.nan
    return float(slope), float(intercept), float(r2)
