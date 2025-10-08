#!/usr/bin/env python3
"""Download and extract canonical Mus musculus gene names.

This mirrors the yeast gene extractor but pulls the NCBI GeneInfo bundle for
mouse. By default it produces:

* A newline-delimited text file of canonical gene symbols.
* A JavaScript module suitable for the existing frontend data format.

Example:
    python scripts/extract_mouse_genes.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List
from urllib.error import URLError
from urllib.request import urlopen


DEFAULT_GENE_INFO_URL = (
    "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
)
DEFAULT_TXT_OUTPUT = Path("client/app/data/genes_lst/mouse_canonical.txt")
DEFAULT_JS_OUTPUT = Path("client/app/data/genes_lst/mouse.js")


def download_gene_info(url: str) -> Path:
    """Download the compressed NCBI gene info table to a temporary file."""
    with tempfile.NamedTemporaryFile(suffix=".gene_info.gz", delete=False) as tmp_handle:
        with urlopen(url) as response:
            shutil.copyfileobj(response, tmp_handle)
        temp_path = Path(tmp_handle.name)
    return temp_path


def parse_canonical_names(gene_info_path: Path) -> List[str]:
    """Extract canonical gene symbols, falling back to locus tags as needed."""
    canonical_names: set[str] = set()

    def add_name(candidate: str) -> None:
        cleaned = candidate.strip()
        if cleaned and cleaned != "-":
            canonical_names.add(cleaned)

    open_kwargs = {"encoding": "utf-8"}
    opener = gzip.open if gene_info_path.suffix == ".gz" else open

    with opener(gene_info_path, "rt", newline="", **open_kwargs) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            # NCBI gene_info schema:
            # 0: tax_id, 1: GeneID, 2: Symbol, 3: LocusTag, 10: Symbol_from_nomenclature_authority
            if len(row) < 11:
                continue
            symbol_authority = row[10].strip()
            symbol = row[2].strip()
            locus_tag = row[3].strip()

            if symbol_authority:
                add_name(symbol_authority)
            elif symbol:
                add_name(symbol)
            elif locus_tag:
                add_name(locus_tag)

    return sorted(canonical_names)


def write_txt(names: Iterable[str], destination: Path) -> None:
    """Write gene names to a newline-delimited text file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for name in names:
            handle.write(f"{name}\n")


def write_js_module(names: Iterable[str], destination: Path) -> None:
    """Write gene names to the frontend JavaScript list expected by the UI."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    names_list = list(names)
    with destination.open("w", encoding="utf-8") as handle:
        handle.write("export default [\n")
        for index, name in enumerate(names_list):
            suffix = "," if index < len(names_list) - 1 else ""
            handle.write(f'    "{name}"{suffix}\n')
        handle.write("];\n")


def run(
    gene_info_url: str,
    gene_info_path: Path | None,
    txt_output: Path,
    js_output: Path | None,
) -> None:
    """Coordinate download, parsing, and artifact generation."""
    cleanup_path: Path | None = None
    try:
        source_path = gene_info_path
        if source_path is None:
            cleanup_path = download_gene_info(gene_info_url)
            source_path = cleanup_path
        elif not source_path.exists():
            raise FileNotFoundError(f"Gene info table not found at {source_path}")

        gene_names = parse_canonical_names(source_path)
        write_txt(gene_names, txt_output)
        if js_output is not None:
            write_js_module(gene_names, js_output)
    finally:
        if cleanup_path and cleanup_path.exists():
            cleanup_path.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract canonical Mus musculus gene names."
    )
    parser.add_argument(
        "--gene-info-url",
        default=DEFAULT_GENE_INFO_URL,
        help="URL for the Mus musculus gene_info.gz file (default: %(default)s)",
    )
    parser.add_argument(
        "--gene-info-path",
        type=Path,
        default=None,
        help="Use a local Mus musculus gene_info.gz file instead of downloading.",
    )
    parser.add_argument(
        "--txt-output",
        type=Path,
        default=DEFAULT_TXT_OUTPUT,
        help="Destination for the newline-delimited gene list.",
    )
    parser.add_argument(
        "--js-output",
        type=Path,
        default=DEFAULT_JS_OUTPUT,
        help="Destination for the JavaScript gene list (omit to skip).",
    )
    parser.add_argument(
        "--skip-js",
        action="store_true",
        help="Only create the text file, leave the JS module untouched.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    js_output = None if args.skip_js else args.js_output

    try:
        run(
            gene_info_url=args.gene_info_url,
            gene_info_path=args.gene_info_path,
            txt_output=args.txt_output,
            js_output=js_output,
        )
    except (URLError, OSError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    sys.exit(main())
