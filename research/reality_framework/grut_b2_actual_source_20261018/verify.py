#!/usr/bin/env python3
"""Check final archive, executed sources and exact-control replay."""
import hashlib
import json
from pathlib import Path

p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
manifest=json.loads((p/'MANIFEST.json').read_text())
for name,h in manifest['files'].items():
    assert sha(p/name)==h,name
sources=json.loads((p/'EXECUTED_SOURCES.json').read_text())
for name,h in sources['files'].items():
    assert sha(p/name)==h,'executed source '+name
primary=json.loads((p/'PRIMARY_WORKER.json').read_text())
replay=json.loads((p/'REPLAY.json').read_text())
assert primary['exit_code']==replay['exit_code']==0
assert primary['output_sha256']==replay['output_sha256']==sha(p/'CONTROLS.json')
assert primary['wall_seconds']+replay['wall_seconds']<60
assert max(primary['max_rss_kib'],replay['max_rss_kib'])<1024**2
assert sum(f.stat().st_size for f in p.iterdir() if f.is_file())<25*1024**2
r=json.loads((p/'CONTROLS.json').read_text())
assert r['status']=='all_frozen_controls_pass'
assert r['actual_operating_gate']=='unresolved' and r['actual_violating_point'] is None
print(json.dumps({'status':'verified','checks':len(r['checks']),
                  'cases':len(r['finite_compression_cases']),
                  'manifest_files':len(manifest['files']),
                  'output_sha256':sha(p/'CONTROLS.json')}))
