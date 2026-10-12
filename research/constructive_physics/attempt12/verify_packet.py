"""Cross-file regression: fail if displayed claims disagree with saved executed results."""
from pathlib import Path
import json,math
import numpy as np
p=Path(__file__).parent
load=lambda name:json.loads((p/name).read_text())
groups=[]
base=load('results.json')
for r in base:
 assert r['residual']<1e-7 and abs(sum(r['graph_probability'].values())-1)<1e-9
 assert r['energies'][0]>=-r['N']/4-1e-9
for n,q in [(3,2),(3,3),(4,3),(4,4)]:
 a=next(r for r in base if (r['N'],r['Q'],r['cutoff'])==(n,q,2))
 b=next(r for r in base if (r['N'],r['Q'],r['cutoff'])==(n,q,3))
 assert abs(a['energies'][0]-b['energies'][0])<1e-8
for q in [5,6]:
 a=load(f'run_N5_Q{q}_L2.json');b=load(f'run_N5_Q{q}_L3.json')
 assert abs(a['energies'][0]-b['energies'][0])<1e-8 and b['cutoff_boundary_probability']<1e-9
 assert b['energies'][0]>=-5/4
 groups.append(f'full_N5_Q{q}_cutoff')
for r in load('reduction_validation.json'):
 assert max(r['orthogonality_error'],r['projected_matrix_error'],r['sector_leakage'])<1e-10
a=load('symmetric_N6_Q6_L2.json');b=load('symmetric_N6_Q6_L3.json');assert abs(a['energies'][0]-b['energies'][0])<1e-8
groups.append('isometry_and_N6_sector_cutoff')
a=load('quench_L4.json');b=load('quench_L5.json')
assert abs(a['controls'][-1]['cycle_probability']-b['controls'][-1]['cycle_probability'])<2e-8
assert b['controls'][-1]['electric_cutoff_boundary']<1e-6
groups.append('quench_cutoff')
a=load('thermal_N4_Q4_L4.json');b=load('thermal_N4_Q4_L5.json')
for x,y in zip(a['thermal'],b['thermal']):
 if x['beta']>=5:assert abs(x['cycle_probability']-y['cycle_probability'])<1e-10
assert b['thermal'][0]['cutoff_boundary_probability']>1e-4 # explicitly unresolved hot cutoff
groups.append('thermal_converged_cold_unresolved_hot')
for r in load('gauge_shortcut_results.json'):assert abs(abs(r['H_powers'][2])-1/60)<1e-12
assert load('distance_results.json')['total_pairs']==495
for r in load('local_reconnection_results.json'):assert r['partition_changes']==r['forbidden_partition_changes']
groups.extend(['gauge_shortcut','bare_distance','locality_component_constraint'])
partner=load('PARTNER_EXECUTION.json');assert abs(partner['full_joint'][-1]['connected_probability']-.96894089748856)<1e-10
groups.append('partner_verbatim_results')
result=dict(status='CP12_PACKET_CROSSCHECK_PASS',groups=groups,not_claimed=['thermal hot-cutoff convergence','full N6 ground state','emergent gravity','experimental validation'])
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
