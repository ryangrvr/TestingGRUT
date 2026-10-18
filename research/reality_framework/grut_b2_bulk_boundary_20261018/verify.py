import hashlib
import json
from pathlib import Path

p=Path(__file__).resolve().parent
digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
m=json.loads((p/'MANIFEST.json').read_text())
for name,h in m['files'].items():
    assert digest(p/name)==h,name
s=json.loads((p/'EXECUTED_SOURCES.json').read_text())
for name,h in s['files'].items():
    assert digest(p/name)==h,'executed source '+name
a=json.loads((p/'PRIMARY_WORKER.json').read_text())
b=json.loads((p/'REPLAY.json').read_text())
assert a['exit_code']==b['exit_code']==0
assert a['output_sha256']==b['output_sha256']==digest(p/'CONTROLS.json')
assert a['wall_seconds']+b['wall_seconds']<60
assert max(a['max_rss_kib'],b['max_rss_kib'])<1024**2
assert sum(f.stat().st_size for f in p.iterdir() if f.is_file())<25*1024**2
r=json.loads((p/'CONTROLS.json').read_text())
assert r['status']=='all_frozen_controls_pass' and r['bounds']['target_pass']
assert r['operating_gate']=='unresolved'
print(json.dumps({'status':'verified','checks':len(r['checks']),
                  'manifest_files':len(m['files']),'output_sha256':digest(p/'CONTROLS.json')}))
