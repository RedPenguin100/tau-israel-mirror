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

import requests
from tqdm import tqdm

from .consts import PROJECT_PATH

# ``mega`` depends on ``tenacity`` versions that still expect ``asyncio.coroutine``.
# Python >=3.12 removed this alias, so recreate it if needed before importing mega.
if not hasattr(asyncio, "coroutine"):
    asyncio.coroutine = types.coroutine


def download(task):
    """Download a file from a URL (supports both regular URLs and mega.nz)."""
    out = Path(task["output"])
    result_path = Path(task.get("result") or out)

    url = task["url"]
    if "mega.nz" in url:
        from mega import Mega
        mega = Mega().login()
        info = mega.get_public_url_info(url)
        name = info.get("name") or "file"
        size = info.get("size")

        out = (out / name) if out.is_dir() or out.suffix == "" else out
        out.parent.mkdir(parents=True, exist_ok=True)
        download_dir = out.parent
        tmp = download_dir / name

        bar = tqdm(total=size, unit="B", unit_scale=True, desc=f"downloading {name}")
        err = []

        def worker():
            try:
                mega.download_url(url, dest_path=str(download_dir))
            except Exception as e:
                err.append(e)

        t = threading.Thread(target=worker, daemon=True)
        t.start()
        prev = 0
        while t.is_alive():
            if tmp.exists():
                cur = tmp.stat().st_size
                bar.update(cur - prev)
                prev = cur
            time.sleep(0.5)
        t.join()
        if tmp.exists():
            cur = tmp.stat().st_size
            bar.update(cur - prev)
            if tmp != out:
                tmp.rename(out)
        bar.close()
        if err:
            raise err[0]
    else:
        out.parent.mkdir(parents=True, exist_ok=True)
        
        # Support for retry logic and timeout
        max_retries = task.get("retries", 3)
        timeout = task.get("timeout", 600)
        
        for attempt in range(max_retries):
            try:
                with requests.get(url, stream=True, timeout=timeout) as r:
                    r.raise_for_status()
                    total = int(r.headers.get("content-length", 0))
                    bar = tqdm(
                        total=total or None,
                        unit="B",
                        unit_scale=True,
                        desc=f"downloading {out.name or 'file'}",
                    )
                    with out.open("wb") as f:
                        for chunk in r.iter_content(1 << 20):
                            if chunk:
                                f.write(chunk)
                                bar.update(len(chunk))
                    bar.close()
                break  # Success, exit retry loop
            except Exception as e:
                if attempt < max_retries - 1:
                    print(f"Download failed (attempt {attempt + 1}/{max_retries}), retrying...")
                    time.sleep(2)
                else:
                    raise e
    
    if task.get("extract"):
        dest = Path(task.get("extract_to") or out.parent)
        dest.mkdir(parents=True, exist_ok=True)
        name = out.name
        if name.endswith(".zip"):
            with zipfile.ZipFile(out) as z:
                for member in tqdm(z.infolist(), unit="file", desc=f"extracting {name}"):
                    z.extract(member, dest)
            if "result" not in task:
                result_path = dest
        elif name.endswith(".tar.gz") or name.endswith(".tgz"):
            with tarfile.open(out, "r:gz") as z:
                for member in tqdm(z.getmembers(), unit="file", desc=f"extracting {name}"):
                    z.extract(member, dest)
            if "result" not in task:
                result_path = dest
        elif name.endswith(".gz"):
            target = dest / out.stem
            bar = tqdm(
                total=out.stat().st_size,
                unit="B",
                unit_scale=True,
                desc=f"extracting {name}"
            )
            with gzip.open(out, "rb") as src, target.open("wb") as dst:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    dst.write(chunk)
                    bar.update(len(chunk))
            bar.close()
            if "result" not in task:
                result_path = target
        if not task.get("keep_archive"):
            out.unlink(missing_ok=True)
    
    return result_path



__all__ = ["download", "ensure_assets"]


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

    # Create directories
    genome_dir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)
    
    tasks = (
        {
            "name": "human_genome",
            "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz",
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
            "url": "https://mega.nz/file/3MERkYxY#XPQdtz-0AMhASxGFvhNliFZEdldrfrp2kYDs5e3Jd-M",
            "output": str(db_dir / "human_gff_basic_introns.db.gz"),
            "extract": True,
            "extract_to": str(db_dir),
            "keep_archive": False,
            "result": str(db_dir / "human_gff_basic_introns.db"),
        },
        # TODO: fix
        # {
        #     "name": "human_index_structure",
        #     "url": "https://mega.nz/file/vdNwgJhD#DnUqX1l7w-yt9yn2xD3lZ7UltYDzD7y4biR_Klswu64",
        #     "output": str(index_dir / "index_structure.zip"),
        #     "extract": True,
        #     "extract_to": str(index_dir),
        #     "keep_archive": False,
        #     "result": str(index_dir / "index_structure"),
        # },
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
            downloaded = download(task)
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
