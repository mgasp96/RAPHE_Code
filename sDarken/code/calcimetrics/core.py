"""Core data containers for calcimetrics."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class Recording:
    """One imaging session (one ImageJ CSV / one field of view)."""

    name: str
    raw: pd.DataFrame
    traces: pd.DataFrame
    fps: float
    excluded_rois: list[str] = field(default_factory=list)
    biological_replicate: str = ""
    condition: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    history: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.fps <= 0:
            raise ValueError("fps must be > 0")
        if not self.biological_replicate:
            self.biological_replicate = self.name
        # Avoid accidental aliasing.
        self.raw = self.raw.copy()
        self.traces = self.traces.copy()

    @property
    def n_rois(self) -> int:
        return int(self.traces.shape[1])

    @property
    def n_frames(self) -> int:
        return int(self.traces.shape[0])

    @property
    def duration_s(self) -> float:
        return float(self.n_frames / self.fps)

    @property
    def time_axis(self) -> np.ndarray:
        return np.arange(self.n_frames, dtype=float) / self.fps

    def exclude_rois(self, rois: list[int] | list[str]) -> "Recording":
        """Return a copy with selected ROI columns removed."""
        if not rois:
            return self.copy()

        columns = list(self.traces.columns)
        drop_names: list[str] = []
        for roi in rois:
            if isinstance(roi, int):
                if roi < 0 or roi >= len(columns):
                    raise IndexError(f"ROI index out of range: {roi}")
                drop_names.append(str(columns[roi]))
            else:
                if roi not in columns:
                    raise KeyError(f"ROI column not found: {roi}")
                drop_names.append(str(roi))

        kept = self.traces.drop(columns=drop_names)
        rec = self.copy()
        rec.traces = kept
        rec.excluded_rois = list(dict.fromkeys([*self.excluded_rois, *drop_names]))
        rec.history.append(
            f"excluded_rois: {drop_names} -> {rec.n_rois} ROIs remaining"
        )
        return rec

    def copy(self) -> "Recording":
        return replace(
            self,
            raw=self.raw.copy(),
            traces=self.traces.copy(),
            excluded_rois=list(self.excluded_rois),
            metadata=dict(self.metadata),
            history=list(self.history),
        )


@dataclass
class Dataset:
    """Collection of recordings with replicate-aware aggregation."""

    recordings: dict[str, Recording] = field(default_factory=dict)

    def add(self, rec: Recording) -> None:
        if self.recordings and rec.n_frames != self._first().n_frames:
            raise ValueError("All recordings in a Dataset must share the same frame count")
        if self.recordings and abs(rec.fps - self._first().fps) > 1e-12:
            raise ValueError("All recordings in a Dataset must share the same fps")
        self.recordings[rec.name] = rec

    def _first(self) -> Recording:
        return next(iter(self.recordings.values()))

    @property
    def n_recordings(self) -> int:
        return len(self.recordings)

    @property
    def replicate_ids(self) -> list[str]:
        return sorted({rec.biological_replicate for rec in self.recordings.values()})

    @property
    def n_replicates(self) -> int:
        return len(self.replicate_ids)

    def group_by_replicate(self) -> dict[str, list[Recording]]:
        grouped: dict[str, list[Recording]] = {}
        for rec in self.recordings.values():
            grouped.setdefault(rec.biological_replicate, []).append(rec)
        return grouped

    def combine_traces(self, level: str = "roi", rename_cols: bool = True) -> pd.DataFrame:
        if not self.recordings:
            return pd.DataFrame()
        if level not in {"roi", "replicate"}:
            raise ValueError("level must be either 'roi' or 'replicate'")

        if level == "roi":
            frames: list[pd.DataFrame] = []
            for rec in self.recordings.values():
                df = rec.traces.copy()
                if rename_cols:
                    df.columns = [f"{rec.name}__{c}" for c in df.columns]
                frames.append(df)
            return pd.concat(frames, axis=1)

        grouped = self.group_by_replicate()
        rep_series: dict[str, pd.Series] = {}
        for rep_id, recs in grouped.items():
            per_rec_means = [r.traces.mean(axis=1) for r in recs]
            rep_series[rep_id] = pd.concat(per_rec_means, axis=1).mean(axis=1)
        return pd.DataFrame(rep_series)

    def summary(self) -> pd.DataFrame:
        rows = []
        for rec in self.recordings.values():
            rows.append(
                {
                    "name": rec.name,
                    "biological_replicate": rec.biological_replicate,
                    "condition": rec.condition,
                    "n_rois": rec.n_rois,
                    "n_frames": rec.n_frames,
                    "duration_s": rec.duration_s,
                    "fps": rec.fps,
                    "history_steps": len(rec.history),
                }
            )

        summary_df = pd.DataFrame(rows)
        if summary_df.empty:
            return summary_df

        total = {
            "name": "TOTAL",
            "biological_replicate": "-",
            "condition": "-",
            "n_rois": int(summary_df["n_rois"].sum()),
            "n_frames": np.nan,
            "duration_s": np.nan,
            "fps": np.nan,
            "history_steps": np.nan,
            "n_recordings": self.n_recordings,
            "n_replicates": self.n_replicates,
        }
        summary_df["n_recordings"] = np.nan
        summary_df["n_replicates"] = np.nan
        return pd.concat([summary_df, pd.DataFrame([total])], ignore_index=True)
