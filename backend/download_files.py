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
import sys

import requests
from tqdm import tqdm

# ``mega`` depends on ``tenacity`` versions that still expect ``asyncio.coroutine``.
# Python >=3.12 removed this alias, so recreate it if needed before importing mega.
if not hasattr(asyncio, "coroutine"):
    asyncio.coroutine = types.coroutine


def download(task):
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
        
        # Final update after thread completes
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
        with requests.get(url, stream=True, timeout=60) as r:
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
    
    if task.get("extract"):
        dest = Path(task.get("extract_to") or out.parent)
        dest.mkdir(parents=True, exist_ok=True)
        name = out.name
        if name.endswith(".zip"):
            with zipfile.ZipFile(out) as z:
                members = z.infolist()
                for member in tqdm(members, unit="file", desc=f"extracting {name}"):
                    z.extract(member, dest)
            if "result" not in task:
                result_path = dest
        elif name.endswith(".tar.gz") or name.endswith(".tgz"):
            with tarfile.open(out, "r:gz") as z:
                members = z.getmembers()
                for member in tqdm(members, unit="file", desc=f"extracting {name}"):
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


def ensure_assets(force: bool = False):
    """Download and extract required assets for the application."""
    
    # Define base paths
    app_root = Path("/app")
    human_root = app_root / "aso_gen" / "data" / "human" / "human_v34"
    human_db_dir = human_root / "dbs"
    off_target_root = app_root / "off_target"
    
    # Create directories
    human_db_dir.mkdir(parents=True, exist_ok=True)
    off_target_root.mkdir(parents=True, exist_ok=True)
    
    # Define download tasks
    tasks = [
        {
            "name": "human_db",
            "url": "https://mega.nz/file/3MERkYxY#XPQdtz-0AMhASxGFvhNliFZEdldrfrp2kYDs5e3Jd-M",
            "output": str(human_db_dir / "human_gff_basic_introns.db.gz"),
            "extract": True,
            "extract_to": str(human_db_dir),
            "keep_archive": False,
            "result": str(human_db_dir / "human_gff_basic_introns.db"),
        },
        {
            "name": "human_index_structure",
            "url": "https://mega.nz/file/vdNwgJhD#DnUqX1l7w-yt9yn2xD3lZ7UltYDzD7y4biR_Klswu64",
            "output": str(off_target_root / "index_structure.zip"),
            "extract": True,
            "extract_to": str(off_target_root),
            "keep_archive": False,
            "result": str(off_target_root / "index_structure"),
        },
    ]

    results = []
    success_count = 0
    fail_count = 0
    
    print("=" * 50)
    print("Starting Python Script Downloads")
    print("=" * 50)
    
    for idx, task in enumerate(tasks, 1):
        result_path = Path(task.get("result") or task["output"])
        
        print(f"\n[{idx}/{len(tasks)}] Processing {task['name']}...")
        
        if not force and result_path.exists():
            print(f"✓ {task['name']} already present at {result_path}")
            results.append(result_path)
            success_count += 1
            continue

        try:
            print(f"Downloading {task['name']}...")
            downloaded = download(task)
            print(f"✓ Stored at {downloaded}")
            results.append(downloaded)
            success_count += 1
        except Exception as e:
            print(f"✗ Error downloading {task['name']}: {str(e)}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            fail_count += 1

    # Summary
    print("\n" + "=" * 50)
    print("Python script completed!")
    print(f"✓ Successful: {success_count}/{len(tasks)}")
    if fail_count > 0:
        print(f"✗ Failed: {fail_count}/{len(tasks)}")
    print("=" * 50)
    
    # Exit with error if any downloads failed
    if fail_count > 0:
        sys.exit(1)
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Download required assets")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files exist"
    )
    
    args = parser.parse_args()
    
    try:
        ensure_assets(force=args.force)
        print("\n✓ All assets downloaded successfully!")
    except Exception as e:
        print(f"\n✗ Failed to download assets: {str(e)}", file=sys.stderr)
        sys.exit(1)