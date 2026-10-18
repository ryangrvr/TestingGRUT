"""Archive integrity and bounded-control identity, not theorem validation."""
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
old=p.parent/'grut_b2_real_axis_20261018'
old_manifest=json.loads((old/'MANIFEST.json').read_text())['files']
for name,h in old_manifest.items():
    assert digest(old/name)==h,'B2-RA archive '+name
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
assert r['fixed_path_status']=='archived_author_derived_pass_independent_review_pending'
assert r['moving_geometry_status']=='existing_hypotheses_incompatible_operating_claim_unresolved'
assert not r['epsilon_geo_achieved'] and r['dynamics_runs']==0
print(json.dumps({'status':'verified_archival_identity_and_exact_local_controls',
                  'checks':len(r['checks']),'B2_RA_archive_files_checked':len(old_manifest),
                  'moving_geometry':'unresolved_no_error_bound_achieved',
                  'independent_proof_review':'pending'}))
