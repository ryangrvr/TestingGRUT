#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    m=json.loads((ROOT/'MATCHING_RESULTS.json').read_text())
    p=json.loads((ROOT/'PHASE_RESULTS.json').read_text())
    checks=[]
    def add(name,ok,evidence):checks.append(dict(name=name,passed=bool(ok),evidence=evidence))
    charges=dict(X=1,Y=3,Sigma=6)
    add('Renormalizable interactions and effective X6 respect the gauge charge',
        -charges['Y']+3*charges['X']==0 and -charges['Sigma']+2*charges['Y']==0 and -charges['Sigma']+6*charges['X']==0 and 5%6!=0,charges)
    q=m['quadratic_matching'];mch=m['matching']
    add('Exact heavy-field elimination agrees with a real-matrix solve and Fourier operator extraction',
        q['max_relative_error']<1e-12 and max(mch['Fourier_relative_errors'].values())<1e-12,
        dict(matrix=q['max_relative_error'],Fourier=mch['Fourier_relative_errors']))
    add('Full mediator quartic solution and gauge transformations pass separate controls',
        q['max_nonlinear_root_residual']<1e-10 and q['max_gauge_potential_invariance_error']<1e-12 and q['min_mediator_mass2_over_M2']>0,
        dict(root=q['max_nonlinear_root_residual'],gauge=q['max_gauge_potential_invariance_error'],min_mass2=q['min_mediator_mass2_over_M2']))
    s=m['stability'];z=np.geomspace(1,1e5,1001)
    bound=(z*z-.005)**2-(z-1)**2/2-.490025*z*z
    add('The supplied full tree potential has a positive global lower bound',
        s['benchmark_kappa']<s['sufficient_reserved_quartic_kappa_max'] and min(bound)>=-1e-12,
        dict(Young_limit=s['sufficient_reserved_quartic_kappa_max'],min_check=float(min(bound)),analytic_condition='z>=1; (z^2-.005)^2>=.990025*z^4 and (z-1)^2<=z^2'))
    c=m['compact_phase']
    add('Compact Gauss law and reduced phase inertia agree with the unreduced kinetic action',
        max(c[k] for k in ['max_Gauss_residual','max_kinetic_reduction_relative_error','max_A0_gauge_shift_error'])<1e-11,c)
    residual=max(r['source_integral_error'] for r in p['rows'])
    relation=max(r['reconstructed_phase_rate_error'] for r in p['rows'])
    bound_error=max(abs(r['p_final'])-abs(r['eta'])*(1-np.exp(-r['N_end'])) for r in p['rows'])
    add('Conditional phase histories satisfy integrated source and absolute production bound',
        residual<1e-8 and relation<1e-10 and bound_error<1e-12,dict(source=residual,relative_phase=relation,bound_excess=bound_error))
    mirror=[];zeros=[]
    for eta in [1e-6,1e-4,.01]:
        for ratio in [.01,1.,6.,30.]:
            rows=[r for r in p['rows'] if r['eta']==eta and r['f_initial_over_Higgs']==ratio]
            plus=next(r['p_final'] for r in rows if r['initial_relative_phase']>0)
            minus=next(r['p_final'] for r in rows if r['initial_relative_phase']<0)
            zeros.append(next(abs(r['p_final']) for r in rows if r['initial_relative_phase']==0))
            mirror.append(abs(plus+minus)/max(abs(plus),1e-100))
    add('Zero-phase and conjugate preparations retain their distinct charge predictions',max(zeros)==0 and max(mirror)<1e-10,
        dict(zero_phase=max(zeros),CP_pair_relative_error=max(mirror),net_gauge_charge=p['net_gauge_charge_generated']))
    add('The phase tolerance control passes and the unspecified cosmic inputs remain explicit',
        p['tolerance_absolute_p_shift']<1e-10 and not p['abundance_predicted'] and not p['entropy_history_selected'],
        dict(phase_tolerance_shift=p['tolerance_absolute_p_shift'],no_abundance_prediction=not p['abundance_predicted']))
    out=dict(scope='Author controls of a specified new tree-level EFT and conditional phase diagnostic, not cosmological, UV-quantum or independent Claude validation.',checks=checks,
             passed=sum(c['passed'] for c in checks),total=len(checks),conditional_ground_response_passes_engineering_target=mch['rate_times_reference_age']<.01,
             full_quantum_matching=False,cosmological_abundance_predicted=False,formation_or_fragmentation_rerun=False,independent_review_completed=False)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for c in checks:print('PASS' if c['passed'] else 'FAIL',c['name'],c['evidence'])
    assert out['passed']==out['total']


if __name__=='__main__':main()
