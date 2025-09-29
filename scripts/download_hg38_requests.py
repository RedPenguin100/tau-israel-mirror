#!/usr/bin/env python3
"""Download a file using requests (kept intentionally small)."""

from __future__ import annotations

import pathlib
import sys

import requests


def download(url: str, destination: pathlib.Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True) as resp:
        resp.raise_for_status()
        with destination.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    fh.write(chunk)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("Usage: download_hg38_requests.py <url> <dest>", file=sys.stderr)
        return 1
    url, dest = argv[1], pathlib.Path(argv[2]).expanduser()
    download(url, dest)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv))
