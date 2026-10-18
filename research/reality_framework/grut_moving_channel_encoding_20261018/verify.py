"""Verify saved exact results and historical files; do not rerun dynamics."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,subprocess
p=Path(__file__).resolve().parent
root=p.parents[2]
digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for name,h in json.loads((p/'MANIFEST.json').read_text())['files'].items():
    assert digest(p/name)==h,name
frozen=json.loads((p/'EXECUTED_SOURCES.json').read_text())
for name,h in frozen['files'].items():assert digest(p/name)==h,'frozen '+name
for source in json.loads((p/'SOURCE_PINS.json').read_text())['files']:
    assert digest(root/source['path'])==source['sha256'],source['path']
old=p.parent/'grut_b2_feasibility_closeout_20261018'
for name,h in json.loads((old/'MANIFEST.json').read_text())['files'].items():
    assert digest(old/name)==h,'closed '+name
a=json.loads((p/'PRIMARY_WORKER.json').read_text());b=json.loads((p/'REPLAY.json').read_text())
assert a['exit_code']==b['exit_code']==0 and b['replay_identical']
assert a['source_freeze_commit']==b['source_freeze_commit']==frozen['commit']
assert a['output_sha256']==b['output_sha256']==digest(p/'DIAGNOSTIC.json')
assert a['wall_seconds']+b['wall_seconds']<60
assert max(a['max_rss_kib'],b['max_rss_kib'])<1024**2
assert sum(f.stat().st_size for f in p.iterdir() if f.is_file())<25*1024**2
r=json.loads((p/'DIAGNOSTIC.json').read_text())
assert r['status']=='pass' and len(r['checks'])==64 and all(x['pass'] for x in r['checks'])
assert r['states_in_ball']==4059 and sum(r['inspected_directed_edges'].values())==16533
assert r['physical_receiver_modes']==1024 and r['graph_radius']==8
assert r['physical_time_trajectories']==0 and not r['actual_S9_receiver_used']
assert r['original_operating_gate']=='unresolved' and not r['all_size_formation_or_delivery_bound']
rat=lambda x:F(int(x['numerator']),int(x['denominator']))
j=F(1,3840);jp=F(1,64)
targets=r['targets']
assert rat(targets['identity_alpha']['H_power_coefficient'])==2*j*j
assert rat(targets['identity_beta']['H_power_coefficient'])==-F(3,2)*j*j
assert rat(targets['exchange_beta']['H_power_coefficient'])==-F(3,4)*j*j*jp**4
assert rat(targets['exchange_alpha']['H_power_coefficient'])==F(13,16)*j*j*jp**4
assert rat(r['strict_noise_interval_upper_t0'])>0 and rat(r['C0'])>0
print(json.dumps({'status':'saved_packet_verified','checks':64,
                  'historical_closeout_unchanged':True,'original_gate':'unresolved',
                  'external_validation':'none'}))
