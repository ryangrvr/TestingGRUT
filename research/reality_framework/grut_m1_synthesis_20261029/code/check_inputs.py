"""Verify restored input bytes against Git blobs in the pinned parent tree."""
from pathlib import Path
import hashlib, json
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[2]
pins=json.loads((ROOT/'INPUT_PINS.json').read_text())
for name,want in pins['files'].items():
    data=(REPO/name).read_bytes()
    got=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if got!=want: raise SystemExit(f'FAIL input {name}: {got} != {want}')
print(f"PASS {len(pins['files'])} parent input Git blob IDs at {pins['parent_commit']}")
