"""Smoothing functions for fluorescence traces."""

from __future__ import annotations

import pandas as pd
from scipy.signal import savgol_filter


def gaussian_smooth(
    traces: pd.DataFrame,
    window_size: int = 5,
    sigma: float = 2.5,
) -> pd.DataFrame:
    """Apply centered Gaussian rolling smoothing."""
    if traces.empty:
        return traces.copy()
    if window_size < 2:
        raise ValueError("window_size must be >= 2")
    return traces.rolling(window=window_size, win_type="gaussian", center=True, min_periods=1).mean(std=sigma)


def savgol_smooth(
    traces: pd.DataFrame,
    window_length: int = 11,
    polyorder: int = 3,
) -> pd.DataFrame:
    """Apply Savitzky-Golay smoothing per ROI."""
    if traces.empty:
        return traces.copy()
    if window_length % 2 == 0:
        window_length += 1
    out = pd.DataFrame(index=traces.index)
    for col in traces.columns:
        out[col] = savgol_filter(
            traces[col].astype(float).to_numpy(),
            window_length=window_length,
            polyorder=polyorder,
            mode="interp",
        )
    return out


def moving_average(
    traces: pd.DataFrame,
    window_size: int = 5,
) -> pd.DataFrame:
    """Apply simple moving average."""
    if traces.empty:
        return traces.copy()
    return traces.rolling(window=window_size, center=True, min_periods=1).mean()
