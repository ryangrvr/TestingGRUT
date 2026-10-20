"""Pack or unpack exact raw diagnostics and subsampled spatial volumes.

Run `python archive_io.py pack` only after the campaign finishes.
Run `python archive_io.py unpack` after a clean clone to reproduce analysis.
Each compressed file is a separate Git blob of a manageable size.
"""
from pathlib import Path
import gzip, hashlib, io, json, sys, tarfile
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def archive(paths,target):
    with target.open('wb') as raw:
      with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9) as compressed:
       with tarfile.open(fileobj=compressed,mode='w') as tar:
        for p in paths:
          info=tarfile.TarInfo('raw/'+p.name);info.size=p.stat().st_size;info.mtime=0;info.mode=0o644
          with p.open('rb') as f:tar.addfile(info,f)
def pack():
    raw=sorted((P/'raw').iterdir());hist=[p for p in raw if '.spatial.' not in p.name];spatial=[p for p in raw if '.spatial.' in p.name]
    archive(hist,P/'raw_diagnostics.tar.gz');groups=[];batch=[];size=0
    for p in spatial:
      if batch and size+p.stat().st_size>2_400_000:groups.append(batch);batch=[];size=0
      batch.append(p);size+=p.stat().st_size
    if batch:groups.append(batch)
    for i,group in enumerate(groups):archive(group,P/f'raw_spatial_{i:02}.tar.gz')
    archives=sorted(P.glob('raw_*.tar.gz'))
    manifest={'raw_files':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in raw},'archives':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in archives},'note':'Float64 scalar histories and exact full-grid component diagnostics; spatial volumes are float32 subsamples for inspection. Full evolving field/momentum histories were not retained. Inputs and exact run commands permit reproduction.'}
    (P/'RAW_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest['archives'],indent=2))
def unpack():
    m=json.loads((P/'RAW_MANIFEST.json').read_text())
    for name,entry in m['archives'].items():
      p=P/name
      if sha(p)!=entry['sha256']:raise RuntimeError('archive changed '+name)
      with tarfile.open(p,'r:gz') as tar:tar.extractall(P,filter='data')
    for name,entry in m['raw_files'].items():
      p=P/'raw'/name
      if sha(p)!=entry['sha256']:raise RuntimeError('raw changed '+name)
    print('Verified',len(m['raw_files']),'raw files')
if __name__=='__main__':
    if sys.argv[1]=='pack':pack()
    elif sys.argv[1]=='unpack':unpack()
    else:raise SystemExit('pack | unpack')
