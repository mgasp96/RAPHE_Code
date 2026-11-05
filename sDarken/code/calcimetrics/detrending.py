"""Detrending functions for dF/F traces."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import signal
from scipy.optimize import curve_fit


def detrend_linear(traces: pd.DataFrame) -> pd.DataFrame:
    """Remove linear trend from each ROI."""
    if traces.empty:
        return traces.copy()
    out = pd.DataFrame(index=traces.index)
    for col in traces.columns:
        out[col] = signal.detrend(traces[col].astype(float).to_numpy())
    return out


def detrend_polynomial(traces: pd.DataFrame, order: int = 2) -> pd.DataFrame:
    """Remove polynomial trend from each ROI."""
    if traces.empty:
        return traces.copy()
    x = np.arange(len(traces), dtype=float)
    out = pd.DataFrame(index=traces.index)
    for col in traces.columns:
        y = traces[col].astype(float).to_numpy()
        coeffs = np.polyfit(x, y, deg=order)
        fit = np.polyval(coeffs, x)
        out[col] = y - fit
    return out


def _exp_decay(x: np.ndarray, a: float, tau: float, c: float) -> np.ndarray:
    return a * np.exp(-x / tau) + c


def detrend_exponential(
    traces: pd.DataFrame,
    fit_window: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """Fit and subtract mono-exponential trend from each ROI."""
    if traces.empty:
        return traces.copy()

    x_full = np.arange(len(traces), dtype=float)
    if fit_window is None:
        fit_idx = slice(0, len(traces))
    else:
        fit_idx = slice(max(0, fit_window[0]), min(len(traces), fit_window[1]))
    x_fit = x_full[fit_idx]

    out = pd.DataFrame(index=traces.index)
    for col in traces.columns:
        y_full = traces[col].astype(float).to_numpy()
        y_fit = y_full[fit_idx]
        p0 = [float(y_fit[0] - y_fit[-1]), max(len(y_fit) / 2.0, 1.0), float(y_fit[-1])]
        popt, _ = curve_fit(_exp_decay, x_fit, y_fit, p0=p0, maxfev=10000)
        fit_curve = _exp_decay(x_full, *popt)
        out[col] = y_full - fit_curve
    return out


def detrend_differencing(traces: pd.DataFrame) -> pd.DataFrame:
    """First-order differencing."""
    if traces.empty:
        return traces.copy()
    return traces.diff().fillna(0.0)


def detrend_highpass(
    traces: pd.DataFrame,
    cutoff_hz: float = 0.01,
    fps: float = 1.0,
    order: int = 3,
) -> pd.DataFrame:
    """Butterworth high-pass filtering."""
    if traces.empty:
        return traces.copy()
    nyquist = fps / 2.0
    if cutoff_hz <= 0 or cutoff_hz >= nyquist:
        raise ValueError("cutoff_hz must be between 0 and Nyquist frequency")

    b, a = signal.butter(order, cutoff_hz / nyquist, btype="highpass")
    out = pd.DataFrame(index=traces.index)
    for col in traces.columns:
        out[col] = signal.filtfilt(b, a, traces[col].astype(float).to_numpy())
    return out
