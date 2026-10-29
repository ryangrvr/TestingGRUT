"""Check delivered compact artifacts without rerunning the PDE or kinetics."""
from pathlib import Path
import json, math
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
j=json.loads((ROOT/'results.json').read_text());arrays=np.load(ROOT/'field_diagnostics.npz')
checks=[]
def check(name,condition,detail=''):
    checks.append((name,bool(condition),detail))
    print(('PASS ' if condition else 'FAIL ')+name+((': '+detail) if detail else ''))
check('fine radial replay reaches t=20000',j['replay']['0.025']['end_mXt']==20000)
check('parent 0–1000 history reproduced',max(j['replay']['parent_max_abs_difference_unscaled'].values())<1e-8)
for t in (10000,15000,20000):
    d=j['replay'][f'parent_checkpoint_{t}_differences'];check(f'parent core checkpoint at t={t}',max(abs(v) for v in d.values())<.02,str(d))
for h in ('0.1','0.05','0.025'):
    q=j['replay'][h];check('radial conserved charge '+h,q['max_charge_plus_absorber_drift']<1e-4)
    check('radial energy plus absorber '+h,q['max_energy_plus_absorber_drift_mX']<2.)
e=j['early'];p=[e[h]['20']['trace_mean_mX'] for h in ('0.1','0.05','0.025')]
check('early pressure mean changes with grid',p[0]>0 and p[1]<0 and p[2]<0,str(p))
err=[abs(e[h]['20']['virial_identity_residual_mX']) for h in ('0.1','0.05','0.025')]
check('early virial residual decreases under refinement',err[0]>3*err[1] and err[1]>3*err[2],str(err))
s=[j['stationary'][h]['t0_trace60_center_mX'] for h in ('0.1','0.05','0.025')]
check('same-lattice stationary pressure approaches zero quadratically',3.9<s[0]/s[1]<4.1 and 3.9<s[1]/s[2]<4.1,str(s))
late=j['late_restart'];check('late stress restart matches the saved core state',abs(late['initial_core_E20_difference_mX'])<1e-8 and abs(late['initial_core_Q20_difference'])<1e-8)
lw=late['windows']['20'];check('late pressure mean depends on window endpoints',lw['20000_21000_dt025']['trace_mean_mX']>0 and lw['20100_21000_dt025']['trace_mean_mX']<0)
check('monopole reproduces radial solver',max(j['axis']['pure_monopole_max_differences_unscaled'].values())<1e-8)
for name,b in j['axis']['budget_controls'].items():
    check('angular budget '+name,b['max_relative_energy_plus_absorber_drift']<2e-5 and b['max_relative_charge_plus_absorber_drift']<1e-10)
fine=j['axis']['axis_spec005_t1000'];coarse=j['axis']['axis_spec01_t1000'];six=j['axis']['axis_specl6_t1000']
check('angular neutral l=2 radial refinement',abs(fine['neutral_quadrupole_energy_mX']/coarse['neutral_quadrupole_energy_mX']-1)<.001)
check('angular neutral l=2 truncation control',abs(six['neutral_quadrupole_energy_mX']/coarse['neutral_quadrupole_energy_mX']-1)<1e-5)
check('neutral quadrupole fraction changes over time',.7<fine['neutral_quadrupole_fraction']<.75 and .8<j['axis']['axis_spec005']['neutral_quadrupole_fraction']<.86)
for name in ('axis_spec005','axis_spec005_t1000','axis_spec005_t1000_w100','axis_spec005_t1000_w500'):
    check('modal Parseval '+name,j['axis'][name]['max_parseval_error']<1e-8)
for name in ('axis_direction005_t1000','axis_direction005_t1000_k02','axis_direction005_t1000_k05'):
    check('synthetic Hankel decomposition '+name,j['axis'][name]['synthetic_error']<1e-12)
a=j['axis']['axis_direction005_t1000']['gradient_proxy']['plus'];b=j['axis']['axis_direction005_t1000_k02']['gradient_proxy']['plus']
check('positive gradient proxy depends on cutoff',a/b>3,f'{a:.6g} vs {b:.6g}')
for name,v in j['feedback'].items():
    early=v['snapshots']['100'];check('feedback pre-return agreement '+name,max(abs(x) for x in early.values())<1e-6)
    check('feedback cell is much denser than assumed cosmological density '+name,v['cell_density_over_physical_ncore']>1e7)
    check('feedback energy conservation '+name,v['max_total_energy_drift_fraction']<1e-5)
ledger=j['ledger'];summed=ledger['cold_increase_mX']+ledger['gas_kinetic_increase_mX']+ledger['work_increase_mX']
check('parent ledger closes with generator defect',abs(ledger['held_drop_mX']-summed)<.001)
for filename in ('timing_sensitivity_a40.json','timing_sensitivity.json'):
    timing=json.loads((ROOT/filename).read_text());branches=timing['branches']
    check('timing injection energy matches and is covered '+filename,abs(branches[0]['held_residual_cost_mX']-branches[1]['held_residual_cost_mX'])<1e-7 and 0<branches[0]['held_residual_cost_mX']<timing['held_residual_mX'])
    check('timing ledger closes '+filename,max(abs(x['ledger_identity_residual_mX']) for x in branches)<1e-5)
    check('timing changes added work '+filename,abs(branches[0]['delta_work_mX']-branches[1]['delta_work_mX'])>0.01)
check('archived fine stress includes full t=20000',arrays['pulse0025_stress']['tau'][-1]==20000)
print(f"{sum(x[1] for x in checks)}/{len(checks)} checks passed")
if not all(x[1] for x in checks):raise SystemExit(1)
