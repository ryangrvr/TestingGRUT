"""Validate the recorded computations, without rerunning expensive eigensolves."""
from pathlib import Path
import json,math
from fractions import Fraction
P=Path(__file__).resolve().parent
def read(name):return json.loads((P/name).read_text())
a=read('clusters_k8_L3.json');b=read('clusters_k8_L4.json');r3={r['key']:r for r in a['rows']};cutoff=[]
for d in [a,b]:
 rows={r['key']:r for r in d['rows']}
 for r in rows.values():
  assert r['residual']<1e-9
  assert abs(r['weight']-(r['relative_energy']-math.fsum(n*rows[k]['weight'] for k,n in r['subclusters'].items())))<1e-11
  assert r['square_embeddings_per_vertex']==r['thin_embeddings_per_vertex']
 for s in d['partial_sums']:
  for host in ['square','tree']:
   value=math.fsum(r['weight']*r[f'{host}_embeddings_per_vertex'] for r in rows.values() if r['edges']<=s['order'])
   assert abs(value-s[f'{host}_relative_per_vertex'])<1e-10
cutoff.append(dict(max_cluster_energy_difference=max(abs(r['relative_energy']-r3[r['key']]['relative_energy']) for r in b['rows']),max_square_partial_sum_difference=max(abs(x['square_relative_per_vertex']-y['square_relative_per_vertex']) for x,y in zip(a['partial_sums'],b['partial_sums']))))
assert cutoff[0]['max_square_partial_sum_difference']<1e-9
order_differences=[dict(order=r['order'],square_minus_tree=r['square_relative_per_vertex']-r['tree_relative_per_vertex']) for r in b['partial_sums']]
assert order_differences[-1]['square_minus_tree']<0 and order_differences[-2]['square_minus_tree']>0
boundary=read('embedding_boundary_results.json');assert len(boundary['different'])==2
assert sorted(r['density_difference'] for r in boundary['different'])==[-1.,1/9]
first=[read(f'winding_weight_l9_L{L}.json') for L in [3,4]]
assert all(r['thin_minus_balanced_first_sensitive_density']>0 for r in first)
assert abs(first[0]['thin_minus_balanced_first_sensitive_density']-first[1]['thin_minus_balanced_first_sensitive_density'])<1e-13
assert all(r['energy_difference']<1e-11 for r in read('random_start_results.json'))
expr=read('tenth_energy_expression.json');assert sum(Fraction(r['coefficient'])*r['Q'] for r in expr)==0
extra=[]
for L in [3,4]:
 f=P/f'tenth_results_L{L}.json'
 if not f.exists():continue
 d=read(f.name);value=math.fsum(float(Fraction(r['coefficient']))*r['energy'] for r in d['rows'])
 assert abs(value-d['tenth_increment'])<1e-13
 assert max(r['residual'] for r in d['rows'])<1e-10
 extra.append(dict(L=L,tenth_increment=d['tenth_increment'],cumulative_difference_through_ten=d['cumulative_difference_through_ten'],residual_indicator=d['residual_indicator']))
assert len(extra)==2
assert all(r['cumulative_difference_through_ten']<0 for r in extra)
assert abs(extra[0]['cumulative_difference_through_ten']-extra[1]['cumulative_difference_through_ten'])<1e-12
assert read('partner_verbatim_results.json')['exit_code']==0
out=dict(status='CP12_FROZEN_CLUSTER_CONTROLS_PASS',cluster_rows_L3=len(a['rows']),cluster_rows_L4=len(b['rows']),cutoff_comparisons=cutoff,square_minus_tree=order_differences,ninth_difference=first,tenth_results=extra,scope='Recorded finite-cluster numerical controls. No thermodynamic remainder or geometry-selection proof.')
(P/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
