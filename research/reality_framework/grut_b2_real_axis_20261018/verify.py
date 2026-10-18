"""Verify archived execution identity and artifact integrity, not the proof."""
import hashlib
import json
from pathlib import Path

p=Path(__file__).resolve().parent
digest=lambda f: hashlib.sha256(f.read_bytes()).hexdigest()
for name, h in json.loads((p/'MANIFEST.json').read_text())['files'].items():
    assert digest(p/name)==h, name
for name, h in json.loads((p/'EXECUTED_SOURCES.json').read_text())['files'].items():
    assert digest(p/name)==h, 'frozen source '+name
for entry in json.loads((p/'SOURCE_PINS.json').read_text())['files']:
    assert digest(p.parents[2]/entry['path'])==entry['sha256'], entry['path']
a=json.loads((p/'PRIMARY_WORKER.json').read_text())
b=json.loads((p/'REPLAY.json').read_text())
assert a['exit_code']==b['exit_code']==0
assert a['output_sha256']==b['output_sha256']==digest(p/'CONTROLS.json')
assert a['wall_seconds']+b['wall_seconds']<60
assert max(a['max_rss_kib'],b['max_rss_kib'])<1024**2
assert sum(f.stat().st_size for f in p.iterdir() if f.is_file())<25*1024**2
r=json.loads((p/'CONTROLS.json').read_text())
assert r['status']=='all_frozen_controls_pass'
assert all(c['pass'] for c in r['checks'])
assert r['operating_gate']=='author_derived_pass_pending_independent_proof_review'
print(json.dumps({'status':'verified_archival_integrity_and_execution',
                  'checks':len(r['checks']),
                  'proof_status':'author_derived_independent_review_pending',
                  'output_sha256':digest(p/'CONTROLS.json')}))
