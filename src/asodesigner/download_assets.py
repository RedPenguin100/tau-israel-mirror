#!/usr/bin/env python3
import asyncio
import gzip
import tarfile
import threading
import time
import types
import zipfile
from pathlib import Path
from typing import Union
import gdown

import requests
from prompt_toolkit import output
from tqdm import tqdm

from .consts import PROJECT_PATH

# ``mega`` depends on ``tenacity`` versions that still expect ``asyncio.coroutine``.
# Python >=3.12 removed this alias, so recreate it if needed before importing mega.
if not hasattr(asyncio, "coroutine"):
    asyncio.coroutine = types.coroutine


__all__ = ["download_google", "ensure_assets"]


def download_google(task):
    """
    Download from Google Drive with gdown + optional retries and extraction.

    Expected `task` keys:
      - url (str): Google Drive share link or file/folder id URL
      - output (str|Path): file path or directory to write into
      - folder (bool, optional): if True, use gdown.download_folder
      - retries (int, optional): number of attempts (default 3)
      - use_cookies (bool, optional): pass cookies to gdown (default False)
      - quiet (bool, optional): gdown progress (default False)
      - fuzzy (bool, optional): allow non-standard URLs (default True)
      - extract (bool, optional): auto-extract zip/tar.gz/gz (default False)
      - extract_to (str|Path, optional): destination for extracted contents
      - keep_archive (bool, optional): keep the downloaded archive (default False)
      - result (str|Path, optional): override returned path
    """
    url = task["url"]
    out = Path(task["output"])
    result_path = Path(task.get("result") or out)

    # Ensure parent directory exists (file output) or the directory itself (folder mode)
    is_folder = bool(task.get("folder", False))
    target_dir = out if (is_folder or out.suffix == "") else out.parent
    target_dir.mkdir(parents=True, exist_ok=True)

    # gdown settings
    retries = int(task.get("retries", 3))
    use_cookies = bool(task.get("use_cookies", False))
    quiet = bool(task.get("quiet", False))
    fuzzy = bool(task.get("fuzzy", True))

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            if is_folder:
                # When downloading a folder, `out` should be a directory
                saved = gdown.download_folder(
                    url=url,
                    output=str(out),
                    quiet=quiet,
                    use_cookies=use_cookies,
                    remaining_ok=True
                )
                # gdown returns a list of saved paths; normalize result
                if not saved:
                    raise RuntimeError("gdown.download_folder returned no files.")
                result_path = Path(out)
            else:
                # File download
                saved = gdown.download(
                    url=url,
                    output=str(out),
                    quiet=quiet,
                    use_cookies=use_cookies,
                    fuzzy=fuzzy
                )
                if not saved:
                    raise RuntimeError("gdown.download returned None/empty path.")
                result_path = Path(saved)
            # Success
            last_err = None
            break
        except Exception as e:
            last_err = e
            if attempt < retries:
                print(f"Download failed (attempt {attempt}/{retries}), retrying...")
                time.sleep(2)
            else:
                raise

    # Optional extraction (file mode only)
    if task.get("extract") and not is_folder:
        dest = Path(task.get("extract_to") or result_path.parent)
        dest.mkdir(parents=True, exist_ok=True)
        name = result_path.name

        if name.endswith(".zip"):
            with zipfile.ZipFile(result_path) as z:
                z.extractall(dest)
            if "result" not in task:
                result_path = dest
        elif name.endswith(".tar.gz") or name.endswith(".tgz"):
            with tarfile.open(result_path, "r:gz") as z:
                z.extractall(dest)
            if "result" not in task:
                result_path = dest
        elif name.endswith(".gz") and not name.endswith(".tar.gz"):
            target = dest / result_path.stem
            with gzip.open(result_path, "rb") as src, target.open("wb") as dst:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    dst.write(chunk)
            if "result" not in task:
                result_path = target

        if not task.get("keep_archive"):
            # best-effort cleanup
            try:
                Path(saved).unlink(missing_ok=True)  # if saved is a string path
            except Exception:
                pass

    return Path(result_path)


def ensure_assets(destination: Union[str, Path, None] = None, force: bool = False):
    """
    Ensure required human assets exist at specified locations relative to this script.
    
    Downloads 3 essential files:
    1. Human genome reference (GRCh38.p13) -> ./data/human/human_v34/
    2. Human gene annotation database -> ./data/human/human_v34/dbs/
    3. Human index structure -> ./off_target/index_structure/
    """

    print("=" * 50)
    print("Starting Asset Downloads")
    print("=" * 50)
    print(f"Cache location: {PROJECT_PATH}")
    print()

    # Define paths relative to script location
    genome_dir = PROJECT_PATH / "data" / "human" / "human_v34"
    genome_dir.expanduser().mkdir(parents=True, exist_ok=True)
    db_dir = PROJECT_PATH / "data" / "human" / "human_v34" / "dbs"
    db_dir.expanduser().mkdir(parents=True, exist_ok=True)
    index_dir = PROJECT_PATH / "off_target"
    index_dir.expanduser().mkdir(parents=True, exist_ok=True)
    project_dir = PROJECT_PATH
    project_dir.expanduser().mkdir(parents=True, exist_ok=True)

    # Create directories
    genome_dir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)
    project_dir.mkdir(parents=True, exist_ok=True)

    tasks = (
        {
            "name": "human_genome",
            "url": "https://drive.google.com/file/d/1zRxJNMtipdurpHo1EZad0YJZ8u6n3X-J/view?usp=drive_link",
            "output": str(genome_dir / "GRCh38.p13.genome.fa.gz"),
            "extract": True,
            "extract_to": str(genome_dir),
            "keep_archive": False,
            "result": str(genome_dir / "GRCh38.p13.genome.fa"),
            "retries": 3,
            "timeout": 600,
        },
        {
            "name": "human_db",
            "url": "https://drive.google.com/file/d/18tmvD9NYUpoC6LCghvVkush5fNOGI29-/view?usp=drive_link",
            "output": str(db_dir / "human_gff_basic_introns.db.gz"),
            "extract": True,
            "extract_to": str(db_dir),
            "keep_archive": False,
            "result": str(db_dir / "human_gff_basic_introns.db"),
        },
        {
            "name": "human_index_structure",
            "url": "https://drive.google.com/file/d/1nmW-IEdytnYA-_zL9nomNqtwMsI2CpB2/view?usp=drive_link",
            "output": str(index_dir / "index_structure.zip"),
            "extract": True,
            "extract_to": str(index_dir),
            "keep_archive": False,
            "result": str(index_dir / "index_structure"),
        },
        {
            "name": "chromosome_stub",
            "url": "https://drive.google.com/file/d/1jwdNFqmQB_CgUgMFHCX-i5eqDSNdaBbP/view?usp=drive_link",
            "output": str(project_dir / "chromosomes.zip"),
            "extract": True,
            "extract_to": str(project_dir),
            "keep_archive": False,
            "result": str(project_dir / "chromosomes"),
        },

    )

    results = []
    for i, task in enumerate(tasks, 1):
        result_path = Path(task.get("result") or task["output"])
        print(f"[{i}/{len(tasks)}] Processing {task['name']}...")
        
        if not force and result_path.exists():
            print(f"✓ {task['name']} already present at {result_path}")
            results.append(result_path)
            print()
            continue

        print(f"Downloading {task['name']}...")
        try:
            downloaded = download_google(task)
            print(f"✓ {task['name']} successfully stored at {downloaded}")
            results.append(downloaded)
        except Exception as e:
            print(f"✗ Failed to download {task['name']}: {e}")
            raise
        print()

    print("=" * 50)
    print("✓ All assets downloaded successfully!")
    print("=" * 50)
    print(f"Total files downloaded: {len(results)}")
    for result in results:
        print(f"  - {result}")
    print()

    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Download required human genome and annotation assets"
    )
    parser.add_argument(
        "--destination",
        "-d",
        type=str,
        default=None,
        help="Destination directory (default: package data folder)",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force re-download even if files exist",
    )
    
    args = parser.parse_args()
    
    try:
        ensure_assets(destination=args.destination, force=args.force)
    except KeyboardInterrupt:
        print("\n\nDownload interrupted by user")
        exit(1)
    except Exception as e:
        print(f"\n\nError: {e}")
        exit(1)
