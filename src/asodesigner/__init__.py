"""ASOdesigner package.

Utilities and feature calculators for designing antisense oligonucleotides.
"""

from importlib import metadata as _metadata

try:  # pragma: no cover - best effort during local development
    __version__ = _metadata.version("asodesigner")
except _metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "1.1.7"

__all__ = ["__version__"]
