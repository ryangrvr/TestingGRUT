"""Saved-output crosscheck; numerical evidence is distinct from the analytic proof."""
from pathlib import Path
import json,hashlib
p=Path(__file__).parent
read=lambda f:json.loads((p/f).read_text())
count=0
for f in sorted(p.glob('scaling_N*.json')):
 r=read(f.name)
 for x in r['rows']:
  assert x['leakage']<=x['leakage_bound']+1e-8
  for d in x['detectors']:
   assert abs(d['response_full']-d['response_block'])<=d['response_error_bound']+1e-8;count+=1
for n in [4,5]:
 a=read(f'scaling_N{n}_L3.json');b=read(f'scaling_N{n}_L4.json')
 for x,y in zip(a['rows'],b['rows']):
  diff=max(abs(u['response_full']-v['response_full']) for u,v in zip(x['detectors'],y['detectors']))
  assert diff<2e-6
  if x['T']<=1:assert diff<1e-10
for x,y in zip(read('scaling_N4_L5_long.json')['rows'],read('scaling_N4_L6_long.json')['rows']):
 assert max(abs(u['response_full']-v['response_full']) for u,v in zip(x['detectors'],y['detectors']))<1e-10
fixed=read('fixed_propagation_results.json')
for n in [6,8]:
 a=next(r for r in fixed if r['N']==n and r['L']==3);b=next(r for r in fixed if r['N']==n and r['L']==4)
 for x,y in zip(a['rows'],b['rows']):assert max(abs(u-v) for u,v in zip(x['response_by_vertex_distance'],y['response_by_vertex_distance']))<1e-7
assert read('operator_results.json')['random_block_tests']==320
assert read('operator_results.json')['max_quadratic_bound_ratio']<1
assert max(x['error'] for x in read('impulse_results.json')['controls'])<1e-10
assert abs(read('mixture_results.json')['neutral_dephasing_energy_cost']-2/15)<1e-12
assembly=read('assembly_results.json')
for x,y in zip(assembly[0]['rows'],assembly[1]['rows']):assert abs(x['connected_probability']-y['connected_probability'])<1e-9
assert assembly[1]['rows'][-1]['connected_probability']>.24
r=dict(status='CP12_OPERATIONAL_WINDOW_CONTROLS_PASS',response_bound_checks=count,random_operator_checks=320,scope='Full N=4,5 dynamics; compressed N=6,8; no inferred scaling fit or geometry formation.')
(p/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
