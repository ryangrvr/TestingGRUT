"""Verify restored parent inputs against Git blob IDs in the pinned tree."""
from pathlib import Path
import hashlib, json
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
pins=json.loads((ROOT/'INPUT_PINS.json').read_text());bad=[]
for name,expected in pins['files'].items():
    f=REPO/name
    if not f.exists():bad.append(name+' missing');continue
    contents=f.read_bytes();actual=hashlib.sha1(b'blob '+str(len(contents)).encode()+b'\0'+contents).hexdigest()
    if actual!=expected:bad.append(name+' changed')
if bad:raise SystemExit('\n'.join(bad))
print('verified',len(pins['files']),'parent Git blob IDs from',pins['parent_commit'])
