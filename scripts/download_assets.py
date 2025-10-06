#!/usr/bin/env python3
import argparse,json,gzip,tarfile,threading,time,zipfile
from pathlib import Path
import requests
from tqdm import tqdm

def download(task,skip):
    out=Path(task.get('output','downloaded'))
    if skip and out.exists():
        return out
    url=task['url']
    if 'mega.nz' in url:
        from mega import Mega
        mega=Mega().login()
        info=mega.get_public_url_info(url)
        name=info.get('name') or 'file'
        size=info.get('size')
        out=(out/name) if out.is_dir() or out.suffix=='' else out
        out.parent.mkdir(parents=True,exist_ok=True)
        tmp=out.parent/name
        bar=tqdm(total=size,unit='B',unit_scale=True,desc=f'downloading {name}')
        err=[]
        def worker():
            try:
                mega.download_url(url,dest_path=out.parent)
            except Exception as e:
                err.append(e)
        t=threading.Thread(target=worker,daemon=True)
        t.start()
        prev=0
        while t.is_alive():
            if tmp.exists():
                cur=tmp.stat().st_size
                bar.update(cur-prev)
                prev=cur
            time.sleep(0.5)
        t.join()
        if tmp.exists():
            cur=tmp.stat().st_size
            bar.update(cur-prev)
            if tmp!=out:
                tmp.rename(out)
        bar.close()
        if err:
            raise err[0]
    else:
        out.parent.mkdir(parents=True,exist_ok=True)
        with requests.get(url,stream=True,timeout=60) as r:
            r.raise_for_status()
            total=int(r.headers.get('content-length',0))
            bar=tqdm(total=total or None,unit='B',unit_scale=True,desc=f'downloading {out.name or "file"}')
            with out.open('wb') as f:
                for chunk in r.iter_content(1<<20):
                    if chunk:
                        f.write(chunk)
                        bar.update(len(chunk))
            bar.close()
    if task.get('extract'):
        dest=Path(task.get('extract_to') or out.parent)
        dest.mkdir(parents=True,exist_ok=True)
        name=out.name
        if name.endswith('.zip'):
            with zipfile.ZipFile(out) as z:
                for member in tqdm(z.infolist(),unit='file',desc=f'extracting {name}'):
                    z.extract(member,dest)
        elif name.endswith('.tar.gz') or name.endswith('.tgz'):
            with tarfile.open(out,'r:gz') as z:
                for member in tqdm(z.getmembers(),unit='file',desc=f'extracting {name}'):
                    z.extract(member,dest)
        elif name.endswith('.gz'):
            target=dest/out.stem
            bar=tqdm(total=out.stat().st_size,unit='B',unit_scale=True,desc=f'extracting {name}')
            with gzip.open(out,'rb') as src,target.open('wb') as dst:
                while True:
                    chunk=src.read(1<<20)
                    if not chunk:
                        break
                    dst.write(chunk)
                    bar.update(len(chunk))
            bar.close()
        if not task.get('keep_archive'):
            out.unlink(missing_ok=True)
    return out

if __name__ == "__main__":
    # Download the MEGA link using the download(task, skip) function
    url = "https://mega.nz/file/3MERkYxY#XPQdtz-0AMhASxGFvhNliFZEdldrfrp2kYDs5e3Jd-M"
    from pathlib import Path

    task = {
        "url": url,
        "out": Path("3MERkYxY_downloaded_file")
    }
    skip = False
    out_path = download(task, skip)
    print(f"Downloaded to: {out_path.resolve()}")



