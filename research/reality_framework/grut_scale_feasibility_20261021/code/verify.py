"""Packet integrity and numerical reproduction, not scientific replication."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

p=Path(__file__).resolve().parents[1]
manifest=p/'MANIFEST.json'
def entries():
    return {str(x.relative_to(p)):dict(bytes=x.stat().st_size,
            sha256=hashlib.sha256(x.read_bytes()).hexdigest())
            for x in sorted(p.rglob('*')) if x.is_file() and x!=manifest
            and '__pycache__' not in x.parts}
if '--write-manifest' in sys.argv:
    manifest.write_text(json.dumps(entries(),indent=2,sort_keys=True)+'\n')
expected=json.loads(manifest.read_text())
actual=entries()
if expected!=actual:
    raise AssertionError(dict(missing=sorted(expected.keys()-actual.keys()),
                              extra=sorted(actual.keys()-expected.keys()),
                              changed=sorted(k for k in expected.keys() & actual.keys()
                                             if expected[k]!=actual[k])))
r=subprocess.run([sys.executable,str(p/'code/analyze.py'),'--check'],capture_output=True,text=True)
if r.returncode:
    raise AssertionError(r.stderr)
summary=json.loads(r.stdout)
print(json.dumps(dict(status='PASS',manifest_files=len(actual),
                     internal_arithmetic_and_archive_checks=summary['checks'],
                     failures=summary['failures'],external_scientific_replications=0),indent=2))
