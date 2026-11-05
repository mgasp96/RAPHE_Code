"""Public API for calcimetrics."""

from importlib.metadata import PackageNotFoundError, version

from .core import Dataset, Recording
from .io import export_cleaned_csv, load_folder, load_imagej_csv

try:
    __version__ = version("calcimetrics")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    "Dataset",
    "Recording",
    "export_cleaned_csv",
    "load_folder",
    "load_imagej_csv",
]
