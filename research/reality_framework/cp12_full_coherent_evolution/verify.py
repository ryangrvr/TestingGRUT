"""Validate recorded runs; does not rerun expensive Hamiltonian evolutions."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def read(n):return json.loads((P/n).read_text())
fields=['attached','connected','path','ring','sector_nonuniform_TV']
runs={p.name:json.loads(p.read_text()) for p in P.glob('full_N*.json')}
runs.update({p.name:json.loads(p.read_text()) for p in P.glob('partial_N*.json')})
checks=0
for name,d in runs.items():
 if 'variance' in d: assert abs(d['variance']-d['analytic_variance'])<1e-8
 for r in d['rows']:
  assert r['norm_error']<1e-8 and r['energy_error']<1e-7
  assert -1e-10<=r['connected']<=r['attached']+1e-10
  assert r['path']<=r['connected']+1e-10 and r['ring']<=r['connected']+1e-10
  checks+=1
cutoffs=[]
for n,q,lo,hi in [(4,3,4,5),(5,4,3,4),(5,4,4,5),(5,5,3,4),(5,5,4,5),(6,5,2,3),(6,5,3,4)]:
 a=runs.get(f'full_N{n}_Q{q}_L{lo}.json');b=runs.get(f'full_N{n}_Q{q}_L{hi}.json') or runs.get(f'partial_N{n}_Q{q}_L{hi}.json')
 if a is None or b is None:continue
 matched=[(r,s) for r in a['rows'] for s in b['rows'] if r['T']==s['T']]
 row=dict(N=n,Q=q,L_low=lo,L_high=hi,max_geometry_difference=max(abs(r[k]-s[k]) for r,s in matched for k in fields),
          max_electric_difference=max(abs(r['electric_square_per_vertex']-s['electric_square_per_vertex']) for r,s in matched),
          max_high_boundary=max(s['cutoff_boundary'] for r,s in matched),
          matched_times=[r['T'] for r,s in matched])
 if n==5 and hi==5:assert row['max_geometry_difference']<3e-6
 cutoffs.append(row)
mix=read('dephased_results.json');batchdiff=0
for d in mix:
 assert abs(d['mixed_energy']-d['coherent_energy']-.08)<1e-10
 ref=runs[f"full_N5_Q4_L{d['L']}.json"]
 for r in d['rows']:
  match=next(s for s in ref['rows'] if s['T']==r['T'])
  for k in ['attached','connected','path']:batchdiff=max(batchdiff,abs(r['coherent_'+k]-match[k]))
assert batchdiff<1e-10
matter=read('matter_selection_results.json')
for r in matter:
 assert abs(r['hopping_fourth_moment']-r['analytic_fourth_moment'])<1e-10
 assert abs(r['rows'][0]['inferred_quartic']-r['E4'])<.003
comp=read('compressed_competitors_results.json')
for L in [3,4]:
 r={x['label']:x for x in comp if x['L']==L}
 assert r['two_triangles']['compressed_energy']<r['ring6']['compressed_energy']
 assert max(x['residual'] for x in r.values())<1e-9
cycles=[]
for n in range(3,9):
 a=read(f'cycle_N{n}_L3.json');b=read(f'cycle_N{n}_L4.json')
 assert abs(a['energy']-b['energy'])<1e-10
 cycles.append(dict(N=n,per_vertex=b['per_vertex'],cutoff_difference=abs(a['energy']-b['energy'])))
short=read('short_time_results.json')
assert max(r['nonuniform_coefficient_L1'] for r in short[:6])<1e-10
assert short[6]['nonuniform_coefficient_L1']>8e-4
audit=read('partner_verbatim_results.json');assert all(r['exit_code']==0 for r in audit['scripts'])
interval=read('interval_composition_control.json');assert interval['max_readout_difference']<1e-10
preparation=read('preparation_geometry_results.json')
assert all(abs(r['means']['largest_fraction']-.796812)<5*r['standard_errors']['largest_fraction'] for r in preparation)
loop=read('loop_series_cutoff_control.json'); assert loop['exact_match_L2']
from fractions import Fraction
assert all(Fraction(x)==0 for x in loop['differences'][:12])
assert Fraction(loop['differences'][12])>0
assert all(r['energy_difference']>0 and r['residual']<1e-11 for r in read('loop_numeric_results.json'))
result=dict(status='CP12_FULL_COHERENT_FOLLOWUP_CONTROLS_PASS',full_timepoints_checked=checks,
            cutoff_comparisons=cutoffs,batch_single_state_max_difference=batchdiff,
            cycle_energies=cycles,scope='Author finite-system controls. No large-N selection or rigorous infinite-cutoff claim.')
(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
