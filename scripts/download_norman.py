from pathlib import Path
import concurrent.futures, subprocess, time, zipfile, json, shutil, hashlib
root=Path(__file__).resolve().parents[1]
out=root/'data/norman'
size=168758985
step=256*1024
chunks=out/'download_parts_small';chunks.mkdir(exist_ok=True)
def fetch(i):
    start=i*step;end=min(size,start+step)-1
    p=chunks/f'{i:04d}.part'
    if p.exists() and p.stat().st_size==end-start+1:return p
    for attempt in range(8):
        time.sleep(attempt * 2)
        url=next(line.removeprefix('Location: ').strip() for line in (out/'range_headers.txt').read_text().splitlines() if line.startswith('Location: '))
        r=subprocess.run(['curl.exe','-L','--fail','--silent','--show-error','--max-time','120','-r',f'{start}-{end}','-o',str(p),url],capture_output=True)
        if r.returncode==0 and p.stat().st_size==end-start+1:
            print(f'chunk {i+1}/644 complete',flush=True);return p
    raise RuntimeError(f'chunk {i} failed: {r.stderr.decode(errors="replace")}')
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    parts=list(pool.map(fetch,range((size+step-1)//step)))
archive=out/'norman.zip'
with archive.open('wb') as w:
    for p in parts:
        with p.open('rb') as r:shutil.copyfileobj(r,w)
with zipfile.ZipFile(archive) as z:
    members=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'crc':i.CRC} for i in z.infolist()]
    print(json.dumps(members,indent=2),flush=True)
    (out/'archive_manifest.json').write_text(json.dumps(members,indent=2))
    targets=[i for i in z.infolist() if i.filename.endswith('perturb_processed.h5ad')]
    if len(targets)!=1:raise RuntimeError('Expected exactly one perturb_processed.h5ad')
    with z.open(targets[0]) as r,(out/'perturb_processed.h5ad').open('wb') as w:shutil.copyfileobj(r,w)
print('EXTRACTED',flush=True)
