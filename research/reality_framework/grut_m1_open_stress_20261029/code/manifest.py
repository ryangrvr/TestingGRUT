"""Write or verify hashes of all deliverable files except the manifest itself."""
from pathlib import Path
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
KEEP={'axis_spec005.npz','axis_spec005_t1000.npz','axis_spec01_t1000.npz','axis_specl6_t1000.npz'}
def payload():
    out={}
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.name=='MANIFEST.json' or any(v=='.work' or v=='__pycache__' or v.startswith('.work_') for v in p.relative_to(ROOT).parts):continue
        if p.name.startswith('axis_spec') and p.suffix=='.npz' and p.name not in KEEP:continue
        out[str(p.relative_to(ROOT))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return {'algorithm':'sha256','files':out}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();v=payload();f=ROOT/'MANIFEST.json'
if args.check:
    expected=json.loads(f.read_text());assert v==expected,'artifact manifest mismatch';print('manifest verified:',len(v['files']),'files')
else:
    f.write_text(json.dumps(v,indent=2)+'\n');print('manifest written:',len(v['files']),'files')
