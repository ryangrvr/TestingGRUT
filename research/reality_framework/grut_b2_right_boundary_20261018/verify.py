"""Verify saved artifacts and exact controls, not the all-size theorem."""
from pathlib import Path
import hashlib
import json

p=Path(__file__).resolve().parent
digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for name,h in json.loads((p/'MANIFEST.json').read_text())['files'].items():
    assert digest(p/name)==h,name
for name,h in json.loads((p/'EXECUTED_SOURCES.json').read_text())['files'].items():
    assert digest(p/name)==h,'frozen source '+name
for x in json.loads((p/'SOURCE_PINS.json').read_text())['files']:
    assert digest(p.parents[2]/x['path'])==x['sha256'],x['path']
old=p.parent/'grut_b2_consolidation_20261018'
archive=json.loads((old/'MANIFEST.json').read_text())['files']
for name,h in archive.items(): assert digest(old/name)==h,'parent '+name
a=json.loads((p/'PRIMARY_WORKER.json').read_text())
b=json.loads((p/'REPLAY.json').read_text())
assert a['exit_code']==b['exit_code']==0
assert a['output_sha256']==b['output_sha256']==digest(p/'CONTROLS.json')
assert b['replay_identical'] and a['wall_seconds']+b['wall_seconds']<60
assert max(a['max_rss_kib'],b['max_rss_kib'])<1024**2
assert sum(f.stat().st_size for f in p.iterdir() if f.is_file())<25*1024**2
r=json.loads((p/'CONTROLS.json').read_text())
assert r['status']=='all_frozen_controls_pass' and all(c['pass'] for c in r['checks'])
assert not r['late_time_error_bound_achieved'] and r['dynamics_runs']==0
assert r['verdict']=='exact_endpoint_embedding_excluded_actual_operation_unresolved'
print(json.dumps({'status':'verified_exact_local_controls_and_saved_identity',
                  'checks':len(r['checks']),'parent_files_unchanged':len(archive),
                  'operating_gate':'unresolved','independent_review':'pending'}))
