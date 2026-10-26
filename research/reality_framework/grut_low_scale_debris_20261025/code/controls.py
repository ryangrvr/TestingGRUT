"""Scientific cross-checks and physical verdicts kept in separate records."""
import ast,csv,json,math,re
import numpy as np
from thermo import BoseGas,adaptive
from relic import run
from run import ROOT,write


def main():
    modes=json.loads((ROOT/'MODE_RESULTS.json').read_text());rel=json.loads((ROOT/'RELEASE_RESULTS.json').read_text())
    dec=json.loads((ROOT/'DECAY_RESULTS.json').read_text());vec=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    bench=json.loads((ROOT/'BENCHMARK_RESULTS.json').read_text());th=json.loads((ROOT/'THRESHOLD_RESULTS.json').read_text())
    combined=json.loads((ROOT/'COMBINED_RESULTS.json').read_text());checks={};numbers={}
    err=max(v for r in rel.values() for v in r['C_vs_real_matrix_endpoint_relative_error'].values())
    fd=max(x['forced_derivative_FD_relative_error'] for x in modes['endpoint_rows'])
    checks['independent_real_matrix_projection']=err<1e-10 and fd<1e-6
    numbers['max_C_vs_real_matrix_error']=err;numbers['max_forced_derivative_FD_error']=fd
    drift=[];tol=[]
    for c in [10,17]:
        a=rel[f'{c}_long'];b=rel[f'{c}_tight']
        tol.append(abs(a['windows']['600-1000']['epsX']-b['windows']['600-1000']['epsX']))
        for k in ['NYH','NYL','NSIG']:
            early=a['windows']['600-1000'][k];late=a['windows']['2000-3000'][k]
            drift.append(abs(late/early-1))
    checks['release_tolerance_and_late_actions']=max(tol)<1e-8 and max(drift)<1e-4
    numbers['max_release_epsilon_tolerance_shift']=max(tol);numbers['max_particle_action_window_shift']=max(drift)
    checks['Gauss_and_work']=max(x['controls']['max_Gauss_over_initial_Ncirc'] for x in rel.values())<1e-8 and max(x['controls']['max_Hamiltonian_work_error_over_initial_energy'] for x in rel.values())<5e-6
    quad=[]
    for t,u in [(.003,.8611951308217225),(.01,.8611951308217225),(.1,.3),(1,.8611951308217225),(10,.8611951308217225),(30,.8611951308217225),(1,1.)]:
        ref=adaptive(t,u);a=BoseGas(64).thermo(t,u);b=BoseGas(128).thermo(t,u)
        quad.append(dict(T_over_m=t,mu_over_m=u,
            gauss64_vs_adaptive=max(abs(a[k]-ref[k])/max(abs(ref[k]),1e-290) for k in ref),
            gauss128_vs_adaptive=max(abs(b[k]-ref[k])/max(abs(ref[k]),1e-290) for k in ref)))
    checks['Bose_integrals_vs_adaptive']=max(x[k] for x in quad for k in ['gauss64_vs_adaptive','gauss128_vs_adaptive'])<1e-8
    eos=max(x['max_thermo_constraint_error'] for x in bench.values())
    en=max(x['initial_energy_relative_error'] for x in bench.values())
    checks['thermodynamic_constraints']=eos<1e-8 and en<1e-8
    numbers['max_thermo_constraint_error']=eos;numbers['max_initial_energy_error']=en
    steps=[];late=[];grids=[]
    for c in [10,17]:
        r=next(x for x in combined['rows'] if x['candidate']==c);ym=r['ymult_raw']
        fine=bench[f'{c}_{ym:.6g}_R_lower']['f_extra']
        steps.append(abs(bench[f'{c}_step520']['f_extra']/fine-1))
        grids.append(abs(bench[f'{c}_quad128']['f_extra']/bench[f'{c}_step520']['f_extra']-1))
        late.append(abs(bench[f'{c}_late']['f_extra']/fine-1))
    checks['relic_step_quadrature_and_late_cutoff']=max(steps)<.001 and max(grids)<1e-7 and max(late)<.002
    numbers['max_520_1040_relic_shift']=max(steps);numbers['max_96_128_relic_shift']=max(grids);numbers['max_late_endpoint_relic_shift']=max(late)
    rt=[]
    for r in th:
        a,b=r['R_L_520_bracket'],r['R_L_1040_bracket'];rt.append(abs(sum(a)/sum(b)-1))
    checks['threshold_resolution']=max(rt)<.003
    numbers['max_threshold_step_shift']=max(rt)
    checks['decay_time_dilation_energy_work']=max(x['quadrature_energy_work_relative_error'] for x in dec['vector_time_dilation_controls'])<2e-6
    core=vec['modes'][2]['changing']['occupation']
    pshift=abs(vec['tighter_start_and_tolerance']['occupation']/core-1)
    plast=abs(vec['shorter_late_endpoint']['occupation']/core-1)
    wr=max(x[k]['scaled_Wronskian_error'] for x in vec['modes'] for k in ['constant','changing'])
    checks['changing_mass_Proca_controls']=pshift<1e-5 and plast<.002 and wr<1e-8
    numbers.update(Proca_start_and_tolerance_shift=pshift,Proca_late_endpoint_shift=plast,max_scaled_Wronskian_error=wr)
    refpath=ROOT/'inputs/REVIEW_BENCHMARK_REPLAY.json'
    if refpath.exists():
        ref=json.loads(refpath.read_text())['result']
        # Review uses 260 sample points, i.e. 259 intervals.
        own=run(ref['eps'],ref['R'],ref['ymult'],steps=259,endT=1e-6,scheme='review_endpoint',table='review')
        shift=abs(own['f_extra']/ref['f_extra']-1)
        checks['one_reviewer_benchmark_replayed']=shift<1e-7
        numbers['reference_replay_relative_relic_shift']=shift
        write(ROOT/'REFERENCE_COMPARISON.json',dict(reference=ref,own_matched_grid=own,relative_f_shift=shift,
            scope='One rerun of pinned reviewer code vs original thermodynamics implementation, same 260 grid points; not outside replication.'))
    else:checks['one_reviewer_benchmark_replayed']=False
    checks={k:bool(v) for k,v in checks.items()}
    physical=[]
    for r in combined['rows']:
        thresholds=[q for q in th if q['candidate']==r['candidate']]
        best=max(q['R_L_1040_bracket'][1] for q in thresholds)
        conservative=bench[f"{r['candidate']}_{r['ymult_raw']:.6g}_R_lower"]
        physical.append(dict(candidate=r['candidate'],R_lower=r['R_lower'],largest_computed_RL_upper=best,
            lower_injection_exceeds_threshold=bool(r['R_lower']>best),
            conservative_f=conservative['f_extra'],f_after_max_future_erasure_bound=conservative['f_after_max_future_conversion_lower_bound'],
            case_L_passes=False,raw_yield_over_target=r['ymult_raw']))
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text())
    for p in ROOT.rglob('*.json'):
        json.loads(p.read_text(),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
    out=dict(author_control_groups=checks,passed=sum(checks.values()),total=len(checks),numerical_evidence=numbers,
        thermodynamic_quadratures=quad,physical_verdicts_separate=physical,
        universal_exclusion_of_all_continuous_phases=False,full_quantum_kinetic_evolution_completed=False)
    write(ROOT/'VERIFICATION_RESULTS.json',out)
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' not in link and not link.startswith('#'):assert (p.parent/link.split('#')[0]).exists(),(p,link)
    assert all(checks.values()),checks
    assert all(x['lower_injection_exceeds_threshold'] and x['f_after_max_future_erasure_bound']>.01 for x in physical)
    print('NUMERICAL_CONTROLS',sum(checks.values()),'/',len(checks),flush=True)
    print('PHYSICAL_VERDICTS',physical,flush=True)


if __name__=='__main__':main()
