"""Input/output helpers for ImageJ ROI CSV workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .core import Dataset, Recording
from .preprocessing import subtract_background


def _looks_like_frame_index(col_name: Any) -> bool:
    s = str(col_name).strip()
    if s == "":
        return True
    try:
        float(s)
        return True
    except ValueError:
        return False


def load_imagej_csv(
    path: str | Path,
    fps: float = 1.0,
    biological_replicate: str | None = None,
    condition: str | None = None,
    first_col_is_index: bool = True,
    last_col_is_background: bool = True,
    background_subtract: bool = True,
    name: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> Recording:
    """Load one ImageJ CSV into a Recording."""
    csv_path = Path(path)
    df = pd.read_csv(csv_path)

    if first_col_is_index and df.shape[1] > 0 and _looks_like_frame_index(df.columns[0]):
        df = df.drop(columns=df.columns[0])

    raw = df.copy()
    if background_subtract and last_col_is_background:
        traces = subtract_background(df, background_col=-1, drop_background=True)
    else:
        traces = df.copy()

    rec_name = name or csv_path.stem
    bio_rep = biological_replicate or rec_name
    rec = Recording(
        name=rec_name,
        raw=raw,
        traces=traces,
        fps=float(fps),
        biological_replicate=bio_rep,
        condition=condition,
        metadata=metadata or {},
        history=[],
    )
    rec.history.append(
        f"loaded: {csv_path.name} ({rec.n_frames} frames, {rec.n_rois} ROIs, fps={rec.fps})"
    )
    if background_subtract and last_col_is_background:
        rec.history.append("background_subtracted: last column")
    return rec


def load_folder(
    folder: str | Path,
    fps: float = 1.0,
    glob_pattern: str = "Results*.csv",
    replicate_map: dict[str, str] | None = None,
    condition: str | None = None,
    **kwargs: Any,
) -> Dataset:
    """Load all matching CSV files from a folder into a Dataset."""
    folder_path = Path(folder)
    dataset = Dataset()
    replicate_map = replicate_map or {}

    for csv_path in sorted(folder_path.glob(glob_pattern)):
        stem = csv_path.stem
        rec = load_imagej_csv(
            csv_path,
            fps=fps,
            biological_replicate=replicate_map.get(stem, stem),
            condition=condition,
            **kwargs,
        )
        dataset.add(rec)
    return dataset


def export_cleaned_csv(
    rec: Recording,
    path: str | Path,
    use_time_index: bool = True,
) -> None:
    """Export current traces to CSV."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = rec.traces.copy()
    if use_time_index:
        df.index = rec.time_axis
        df.index.name = "time_s"
    df.to_csv(out_path)


def load_roi_set(path: str | Path) -> dict[str, Any]:
    """Parse an ImageJ RoiSet.zip and return ROI metadata."""
    try:
        from read_roi import read_roi_zip  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "read-roi is required for load_roi_set(). Install with `pip install read-roi`."
        ) from exc

    roi_path = Path(path)
    return read_roi_zip(str(roi_path))
