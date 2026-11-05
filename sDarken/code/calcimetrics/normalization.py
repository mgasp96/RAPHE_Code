"""Normalization utilities for dF/F computations."""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_dff(
    traces: pd.DataFrame,
    method: str = "median",
    percentile: float = 0.08,
    frame_range: tuple[int, int] | None = None,
    rolling_window: int | None = None,
) -> pd.DataFrame:
    """Compute dF/F = (F - F0) / F0 for each ROI column."""
    if traces.empty:
        return traces.copy()

    valid_methods = {"median", "mean", "percentile", "frame_range", "rolling"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {sorted(valid_methods)}")

    df = traces.astype(float).copy()

    if method == "median":
        f0 = df.median(axis=0)
        return (df - f0) / f0

    if method == "mean":
        f0 = df.mean(axis=0)
        return (df - f0) / f0

    if method == "percentile":
        f0 = df.quantile(percentile, axis=0)
        return (df - f0) / f0

    if method == "frame_range":
        if frame_range is None:
            raise ValueError("frame_range must be provided when method='frame_range'")
        start, end = frame_range
        f0 = df.iloc[start:end, :].mean(axis=0)
        return (df - f0) / f0

    # rolling baseline (windowed percentile)
    if rolling_window is None or rolling_window < 2:
        raise ValueError("rolling_window must be >= 2 when method='rolling'")
    f0 = (
        df.rolling(window=rolling_window, min_periods=1, center=True)
        .quantile(percentile)
        .replace(0, np.nan)
    )
    return (df - f0) / f0


def to_percent(traces: pd.DataFrame) -> pd.DataFrame:
    """Convert dF/F to percent dF/F."""
    return traces * 100.0
