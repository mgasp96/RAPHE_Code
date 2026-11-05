"""Matplotlib plotting helpers."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from .core import Recording


def _x_axis(rec: Recording) -> np.ndarray:
    return rec.time_axis


def plot_all_traces(
    rec: Recording,
    ncols: int = 4,
    sharey: bool = True,
    ylabel: str = "ΔF/F",
    title_per_roi: bool = True,
    figsize: tuple[float, float] | None = None,
) -> Figure:
    n_rois = rec.n_rois
    nrows = max(1, math.ceil(n_rois / ncols))
    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=sharey, figsize=figsize or (4 * ncols, 2.4 * nrows))
    axes_arr = np.atleast_1d(axes).reshape(nrows, ncols)
    x = _x_axis(rec)

    for idx, col in enumerate(rec.traces.columns):
        r, c = divmod(idx, ncols)
        ax = axes_arr[r, c]
        ax.plot(x, rec.traces[col], linewidth=0.8)
        if title_per_roi:
            ax.set_title(str(col), fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Time (s)")

    for idx in range(n_rois, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes_arr[r, c].axis("off")

    fig.tight_layout()
    return fig


def plot_overlay(
    rec: Recording,
    ylabel: str = "ΔF/F",
    alpha: float = 0.5,
) -> Figure:
    fig, ax = plt.subplots(figsize=(10, 4))
    x = _x_axis(rec)
    for col in rec.traces.columns:
        ax.plot(x, rec.traces[col], alpha=alpha, linewidth=0.8)
    ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{rec.name} overlay ({rec.n_rois} ROIs)")
    fig.tight_layout()
    return fig


def plot_mean_with_error(
    rec: Recording,
    error: str = "sd",
    stim_window: tuple[float, float] | None = None,
    ylabel: str = "ΔF/F",
    color: str = "black",
) -> Figure:
    valid = {"sd", "sem", "ci95"}
    if error not in valid:
        raise ValueError(f"error must be one of {sorted(valid)}")

    mean = rec.traces.mean(axis=1)
    if error == "sd":
        err = rec.traces.std(axis=1)
    elif error == "sem":
        err = rec.traces.sem(axis=1)
    else:
        err = rec.traces.sem(axis=1) * 1.96

    x = _x_axis(rec)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, mean, color=color, linewidth=1.2)
    ax.fill_between(x, mean - err, mean + err, alpha=0.3, color=color)
    ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
    if stim_window is not None:
        ax.axvspan(stim_window[0], stim_window[1], color="gray", alpha=0.25, lw=0)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{rec.name} mean ± {error}")
    fig.tight_layout()
    return fig


def plot_heatmap(
    rec: Recording,
    ylabel: str = "ROI",
    cbar_label: str = "ΔF/F",
    stim_lines: list[float] | None = None,
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> Figure:
    fig, ax = plt.subplots(figsize=figsize)
    data = rec.traces.T.to_numpy()
    x = _x_axis(rec)
    x_edges = np.concatenate([x, [x[-1] + (x[1] - x[0] if len(x) > 1 else 1.0)]])
    y_edges = np.arange(0, rec.n_rois + 1)

    # Default to symmetric colour scale centred on 0
    if vmin is None and vmax is None:
        absmax = float(np.nanpercentile(np.abs(data), 99))
        vmin, vmax = -absmax, absmax

    mesh = ax.pcolormesh(x_edges, y_edges, data, shading="auto",
                         cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)
    if stim_lines:
        for line in stim_lines:
            ax.axvline(line, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{rec.name} heatmap")
    fig.tight_layout()
    return fig


def plot_comparison(
    raw: pd.Series,
    smoothed: pd.Series,
    detrended: pd.Series,
    kernel: np.ndarray | None = None,
    fps: float = 1.0,
) -> Figure:
    """Diagnostic figure for one ROI trace."""
    x = np.arange(len(raw), dtype=float) / fps
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.ravel()

    axes[0].plot(x, raw, linewidth=0.9)
    axes[0].set_title("Raw")
    axes[0].set_xlabel("Time (s)")

    axes[1].plot(x, smoothed, linewidth=0.9, color="tab:orange")
    axes[1].set_title("Smoothed")
    axes[1].set_xlabel("Time (s)")

    axes[2].plot(x, detrended, linewidth=0.9, color="tab:green")
    axes[2].set_title("Detrended")
    axes[2].set_xlabel("Time (s)")

    if kernel is not None:
        axes[3].plot(kernel, linewidth=1.0, color="tab:red")
        axes[3].set_title("Kernel")
        axes[3].set_xlabel("Sample")
    else:
        axes[3].axis("off")

    fig.tight_layout()
    return fig
