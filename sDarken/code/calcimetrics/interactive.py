"""Interactive notebook helpers (optional dependencies)."""

from __future__ import annotations

from dataclasses import replace

import plotly.graph_objects as go
import ipywidgets as widgets
from IPython.display import display

from .core import Dataset, Recording


def roi_selector(rec: Recording) -> Recording:
    """Interactive ROI checkbox selector that returns a filtered Recording."""
    checks = [widgets.Checkbox(value=True, description=str(c)) for c in rec.traces.columns]
    box = widgets.VBox(checks)
    out = widgets.Output()
    done = widgets.Button(description="Apply selection")

    selected: Recording = rec.copy()

    def _render() -> None:
        with out:
            out.clear_output(wait=True)
            chosen = [c.description for c in checks if c.value]
            fig = go.Figure()
            for col in chosen:
                fig.add_scatter(y=rec.traces[col], mode="lines", name=col, showlegend=False)
            fig.update_layout(height=300, margin=dict(l=20, r=20, t=30, b=20))
            fig.show()

    def _apply(_):
        nonlocal selected
        chosen = [c.description for c in checks if c.value]
        selected = replace(selected, traces=rec.traces[chosen].copy())
        selected.history.append(f"roi_selector: kept {len(chosen)} of {rec.n_rois} ROIs")

    for cb in checks:
        cb.observe(lambda _: _render(), names="value")
    done.on_click(_apply)

    _render()
    display(widgets.HBox([box, widgets.VBox([done, out])]))
    return selected


def trace_browser(dataset: Dataset) -> None:
    """Interactive browser for recordings in a Dataset."""
    if not dataset.recordings:
        print("Dataset is empty.")
        return

    dropdown = widgets.Dropdown(options=list(dataset.recordings.keys()), description="Recording:")
    out = widgets.Output()

    def _render(_=None) -> None:
        with out:
            out.clear_output(wait=True)
            rec = dataset.recordings[dropdown.value]
            fig = go.Figure()
            for col in rec.traces.columns:
                fig.add_scatter(y=rec.traces[col], mode="lines", name=str(col), showlegend=False)
            fig.update_layout(
                title=f"{rec.name} ({rec.n_rois} ROIs)",
                xaxis_title="Frame",
                yaxis_title="Signal",
                height=350,
            )
            fig.show()

    dropdown.observe(_render, names="value")
    _render()
    display(widgets.VBox([dropdown, out]))
