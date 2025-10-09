"""Lightweight smoke test for :func:`get_locus_to_data_dict`."""

import gzip
import os
import shutil
import sys
from pathlib import Path

import requests
from tqdm import tqdm

if __package__ in (None, ""):
    PACKAGE_ROOT = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(PACKAGE_ROOT))
    from .read_human_genome import get_locus_to_data_dict
    from .consts import HUMAN_GTF, HUMAN_GTF_GZ, HUMAN_DB_BASIC_INTRONS
else:  # pragma: no cover - exercised when imported as part of the package
    from .read_human_genome import get_locus_to_data_dict
    from .consts import HUMAN_GTF, HUMAN_GTF_GZ, HUMAN_DB_BASIC_INTRONS


GRCH38_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.basic.annotation.gtf.gz"

DEFAULT_FASTA = Path(__file__).resolve().parents[1] / "cache" / "genomes" / "hg38" / "hg38.fa.gz"


def _download_with_progress(url: str, target: Path, label: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        bar_kwargs = {
            "total": total or None,
            "unit": "B",
            "unit_scale": True,
            "unit_divisor": 1024,
            "desc": label,
        }
        with target.open("wb") as fh, tqdm(**bar_kwargs) as bar:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                fh.write(chunk)
                bar.update(len(chunk))


def download_hg38(target: Path) -> None:
    print(f"Downloading hg38 to {target}")
    _download_with_progress(GRCH38_URL, target, "hg38.fa.gz")


def ensure_hg38_path() -> Path:
    dest = Path(os.environ.get("ASODESIGNER_HG38_PATH", DEFAULT_FASTA)).expanduser()
    if not dest.exists():
        download_hg38(dest)
    else:
        print(f"Using hg38 at {dest}")
    return dest


def _decompress_gzip(src: Path, dest: Path) -> None:
    print(f"Decompressing {src} -> {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as gz, dest.open("wb") as out_fh:
        shutil.copyfileobj(gz, out_fh)


def ensure_human_gtf() -> Path:
    if HUMAN_GTF.is_file():
        print(f"Using GTF at {HUMAN_GTF}")
        return HUMAN_GTF

    if HUMAN_GTF_GZ.is_file():
        _decompress_gzip(HUMAN_GTF_GZ, HUMAN_GTF)
        return HUMAN_GTF

    print(f"Downloading GTF to {HUMAN_GTF_GZ}")
    _download_with_progress(GTF_URL, HUMAN_GTF_GZ, HUMAN_GTF_GZ.name)
    _decompress_gzip(HUMAN_GTF_GZ, HUMAN_GTF)
    return HUMAN_GTF


def run_smoke_test(gene_subset=None) -> None:
    fasta_path = ensure_hg38_path()
    gtf_path = ensure_human_gtf()
    HUMAN_DB_BASIC_INTRONS.parent.mkdir(parents=True, exist_ok=True)
    genes = list(gene_subset or ("TP53", "BRCA1", "EGFR"))
    locus_to_data = get_locus_to_data_dict(gene_subset=genes, create_db=True)
    print(f"Resolved {len(locus_to_data)} loci using FASTA {fasta_path} and GTF {gtf_path}")
    for gene in genes:
        data = locus_to_data.get(gene)
        exon_count = len(getattr(data, "exons", []) or []) if data else 0
        print(f"  {gene}: {exon_count} exons")


if __name__ == "__main__":
    run_smoke_test()
