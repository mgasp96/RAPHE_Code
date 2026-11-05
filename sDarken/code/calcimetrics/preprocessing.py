"""Preprocessing operations on raw fluorescence traces."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def subtract_background(
    df: pd.DataFrame,
    background_col: int | str = -1,
    drop_background: bool = True,
) -> pd.DataFrame:
    """Subtract one background ROI column from all other columns."""
    if df.empty:
        return df.copy()

    if isinstance(background_col, int):
        bg_name = df.columns[background_col]
    else:
        bg_name = background_col
    background = df[bg_name]
    corrected = df.sub(background, axis=0)
    if drop_background:
        corrected = corrected.drop(columns=[bg_name])
    return corrected


def _exp_decay(x: np.ndarray, a: float, tau: float, c: float) -> np.ndarray:
    return a * np.exp(-x / tau) + c


def _bi_exp_decay(
    x: np.ndarray, a1: float, tau1: float, a2: float, tau2: float, c: float
) -> np.ndarray:
    return a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2) + c


def correct_photobleach(
    traces: pd.DataFrame,
    method: str = "exponential",
    fit_window: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """Correct photobleaching on raw fluorescence values."""
    if traces.empty:
        return traces.copy()

    valid_methods = {"exponential", "biexponential", "linear", "ratio"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {sorted(valid_methods)}")

    corrected = pd.DataFrame(index=traces.index)
    x_full = np.arange(len(traces), dtype=float)

    if fit_window is None:
        fit_idx = slice(0, len(traces))
    else:
        fit_idx = slice(max(fit_window[0], 0), min(fit_window[1], len(traces)))

    x_fit = x_full[fit_idx]
    if len(x_fit) < 3:
        raise ValueError("fit_window must include at least 3 frames")

    for col in traces.columns:
        y_full = traces[col].astype(float).to_numpy()
        y_fit = y_full[fit_idx]

        if method in {"exponential", "ratio"}:
            p0 = [float(y_fit[0] - y_fit[-1]), max(len(y_fit) / 2.0, 1.0), float(y_fit[-1])]
            popt, _ = curve_fit(_exp_decay, x_fit, y_fit, p0=p0, maxfev=10000)
            fit_curve = _exp_decay(x_full, *popt)
            if method == "ratio":
                baseline = np.nanmedian(fit_curve)
                corrected[col] = (y_full / fit_curve) * baseline
            else:
                corrected[col] = y_full - fit_curve

        elif method == "biexponential":
            p0 = [
                float((y_fit[0] - y_fit[-1]) * 0.6),
                max(len(y_fit) / 3.0, 1.0),
                float((y_fit[0] - y_fit[-1]) * 0.4),
                max(len(y_fit), 2.0),
                float(y_fit[-1]),
            ]
            popt, _ = curve_fit(_bi_exp_decay, x_fit, y_fit, p0=p0, maxfev=20000)
            fit_curve = _bi_exp_decay(x_full, *popt)
            corrected[col] = y_full - fit_curve

        else:  # linear
            coeffs = np.polyfit(x_fit, y_fit, deg=1)
            fit_curve = np.polyval(coeffs, x_full)
            corrected[col] = y_full - fit_curve

    return corrected
