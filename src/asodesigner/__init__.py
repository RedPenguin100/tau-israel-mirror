"""ASOdesigner package.

Utilities and feature calculators for designing antisense oligonucleotides.
"""

from importlib import metadata as _metadata

try:  # pragma: no cover - best effort during local development
    __version__ = _metadata.version("asodesigner")
except _metadata.PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.1.2"

from .genome import get_hg38_genome_path

__all__ = ["__version__", "get_hg38_genome_path"]
