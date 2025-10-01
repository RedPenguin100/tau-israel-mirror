"""Minimal helper for downloading and caching the human hg38 genome."""
import gzip
import os
import shutil
from pathlib import Path
from urllib import request

from .consts import HG38_CACHE_DIR

HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
HG38_URL_ENV = "ASODESIGNER_HG38_URL"
HG38_PATH_ENV = "ASODESIGNER_HG38_PATH"

_ARCHIVE_NAME = "hg38.fa.gz"
_FASTA_NAME = "hg38.fa"


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request.urlretrieve(url, destination)


def _extract(archive: Path, fasta: Path) -> None:
    with gzip.open(archive, "rb") as src, fasta.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def get_hg38_genome_path() -> Path:
    """Return a local hg38 FASTA file, downloading it once if needed."""
    override = os.environ.get(HG38_PATH_ENV)
    if override:
        path = Path(override).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"ASODESIGNER_HG38_PATH points to missing file: {path}")
        return path

    archive = HG38_CACHE_DIR / _ARCHIVE_NAME
    fasta = HG38_CACHE_DIR / _FASTA_NAME

    if fasta.is_file():
        return fasta

    url = os.environ.get(HG38_URL_ENV, HG38_URL)
    if not archive.is_file():
        _download(url, archive)

    _extract(archive, fasta)
    return fasta


__all__ = ["get_hg38_genome_path"]
