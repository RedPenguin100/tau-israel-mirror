#!/usr/bin/env python3
"""Download and extract canonical Saccharomyces cerevisiae gene names.

The script pulls the SGD feature table, isolates canonical (standard) gene names
for ORFs, and emits two artifacts by default:

* A plain-text file with one gene name per line (useful for quick lookup).
* A JavaScript module mirroring the existing `yeast.js` format used by the UI.

Example:
    python scripts/extract_yeast_genes.py
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List
from urllib.error import URLError
from urllib.request import urlopen


DEFAULT_FEATURES_URL = (
    "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
)
DEFAULT_TXT_OUTPUT = Path("client/app/data/genes_lst/yeast_canonical.txt")
DEFAULT_JS_OUTPUT = Path("client/app/data/genes_lst/yeast.js")


def download_features(url: str) -> Path:
    """Download the SGD features table to a temporary file."""
    with tempfile.NamedTemporaryFile(suffix=".tab", delete=False) as tmp_handle:
        with urlopen(url) as response:
            shutil.copyfileobj(response, tmp_handle)
        temp_path = Path(tmp_handle.name)
    return temp_path


def parse_canonical_names(features_path: Path) -> List[str]:
    """Extract canonical gene names, falling back to the systematic name."""
    canonical_names: set[str] = set()
    fallback_names: set[str] = set()

    with features_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 5:
                continue
            if row[1] != "ORF":
                continue

            systematic_name = row[3].strip()
            standard_name = row[4].strip()

            if standard_name:
                canonical_names.add(standard_name)
            elif systematic_name:
                fallback_names.add(systematic_name)

    ordered = sorted(canonical_names) + sorted(
        name for name in fallback_names if name not in canonical_names
    )
    return ordered


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
    features_url: str,
    features_path: Path | None,
    txt_output: Path,
    js_output: Path | None,
) -> None:
    """Coordinate download, parsing, and artifact generation."""
    cleanup_path: Path | None = None
    try:
        source_path = features_path
        if source_path is None:
            cleanup_path = download_features(features_url)
            source_path = cleanup_path
        elif not source_path.exists():
            raise FileNotFoundError(f"Feature table not found at {source_path}")

        gene_names = parse_canonical_names(source_path)
        write_txt(gene_names, txt_output)
        if js_output is not None:
            write_js_module(gene_names, js_output)
    finally:
        if cleanup_path and cleanup_path.exists():
            cleanup_path.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract canonical Saccharomyces cerevisiae gene names."
    )
    parser.add_argument(
        "--features-url",
        default=DEFAULT_FEATURES_URL,
        help="URL for the SGD features table (default: %(default)s)",
    )
    parser.add_argument(
        "--features-path",
        type=Path,
        default=None,
        help="Use a local SGD features table instead of downloading.",
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
            features_url=args.features_url,
            features_path=args.features_path,
            txt_output=args.txt_output,
            js_output=js_output,
        )
    except (URLError, OSError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    sys.exit(main())
