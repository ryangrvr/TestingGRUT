"""Generate or check SHA-256 hashes for the committed local-plume packet."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'MANIFEST.json'

def files():
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and p!=MANIFEST
                  and not any(part in ('.work','__pycache__') for part in p.relative_to(ROOT).parts)
                  and not p.name.startswith('.work_'))

def records():
    return {str(p.relative_to(ROOT)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
            for p in files()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    now=records()
    if a.check:
        old=json.loads(MANIFEST.read_text())['files']
        different=sorted(k for k in set(now)|set(old) if now.get(k)!=old.get(k))
        print(json.dumps({'files':len(now),'mismatches':different},indent=2))
        if different:raise SystemExit(1)
    else:
        MANIFEST.write_text(json.dumps({'algorithm':'sha256','scope':'This directory excluding MANIFEST.json, .work, work logs and Python caches.','files':now},indent=2)+'\n')
        print(f'Wrote {len(now)} file records')

if __name__=='__main__':main()
