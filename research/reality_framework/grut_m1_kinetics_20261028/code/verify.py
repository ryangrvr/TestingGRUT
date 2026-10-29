"""Same-author verification of the bounded kinetic packet; no empirical verdict."""
from setup import *
from transport import grid,core_operator,Model,initial
import ast,platform,subprocess,scipy,matplotlib

def main():
    checks=[]
    add=lambda name,passed,**detail:checks.append({'name':name,'pass':bool(passed),**detail})
    src=json.loads((ROOT/'SOURCE_CHECKS.json').read_text())
    add('exact_inherited_pulse_history',src['history_scaled_error']==0,error=src['history_scaled_error'])
    c=json.loads((ROOT/'CORE_CHECKS.json').read_text())
    err=max(x['max_probability_difference'] for x in c['same_operator_comparison'])
    add('inherited_scattering_and_reciprocity',err<1e-10 and c['max_reciprocity_error']<1e-9,
        probability_difference=err,reciprocity_error=c['max_reciprocity_error'],unitarity_error=c['max_unitarity_error'])
    ph=json.loads((ROOT/'PHYSICS_CHECKS.json').read_text())
    therm=max(abs(x['relative_error']) for x in ph['thermal_inverse_decay_normalization'])
    bw=max(abs(x['relative_difference']) for x in ph['finite_width_comparison'])
    add('independent_resonant_normalization',therm<1e-9 and bw<.001,thermal_error=therm,sampled_finite_width_difference=bw)
    e=grid(0,8,'translated',32);L,om,rr,err=core_operator(e);n=len(e)
    j=int(np.argmin(abs(e-1.38)));i=int(np.argmin(abs(e-(e[j]+2*OMEGA))))
    support=np.array([i,n+j]);outside=np.ones(2*n,bool);outside[support]=False
    leak=float(np.max(abs(L[outside][:,support])))
    column_error=float(np.max(abs(L[:,support].sum(axis=0))))
    add('closed_two_state_exchange_no_grid_diffusion',leak<1e-10 and column_error<1e-10,
        incident_E=e[j],matched_E=e[i],outside_pair_generator=leak,number_column_error=column_error)
    model=Model(e,'full',6);y,_=initial(e);rng=np.random.default_rng(239)
    delta=rng.normal(size=len(y));h=1e-4
    analytic=model.rates(y,True)@delta
    numeric=(model.rates(y+h*delta)[0]-model.rates(y-h*delta)[0])/(2*h)
    jerr=float(np.linalg.norm(analytic-numeric)/np.linalg.norm(numeric))
    add('analytic_jacobian_directional_check',jerr<1e-7,relative_error=jerr)
    primary=json.loads((ROOT/'translated_full_96.json').read_text())
    numerical=primary['checks'];groups=['resonant','elastic']
    local=max(numerical['core_column_invariant_error'],*[numerical[g]['energy_error'] for g in groups])
    drift=numerical['invariant_relative_drift'];minbin=numerical['min_bin_population']
    add('event_and_trajectory_ledgers',local<1e-9 and drift<1e-9 and minbin>-1e-6,
        event_energy_error=local,invariant_relative_drift=drift,minimum_population=minbin,
        charge_note='Qcore is reconstructed algebraically, so charge closure alone is not an independent ODE test')
    finals={}
    for name in ('translated_full_32','translated_full_64','translated_full_96','translated_full_64_range',
                 'translated_full_64_angles','translated_core_64','translated_resonant_64'):
        z=json.loads((ROOT/(name+'.json')).read_text());finals[name]=z['snapshots'][-1]
    base=finals['translated_full_64'];fine=finals['translated_full_96']
    refine={k:abs(fine[k]/base[k]-1) for k in ('Nplus','Nminus','gas_kinetic_mX')}
    add('graded_translated_grid_refinement',refine['Nplus']<.01 and refine['Nminus']<.03 and refine['gas_kinetic_mX']<.01,
        relative_changes_64_to_96=refine,coarse_32=finals['translated_full_32'])
    comparisons={}
    for tag in ('translated_full_64_range','translated_full_64_angles'):
        comparisons[tag]={k:abs(finals[tag][k]/base[k]-1) for k in ('Nplus','Nminus','gas_kinetic_mX')}
    add('energy_support_and_elastic_angle_controls',max(v for row in comparisons.values() for v in row.values())<.01,
        relative_changes=comparisons,primary_omitted_event_integral=fine['omitted_event_integral'])
    eq=primary['checks']['equilibrium']
    add('gas_only_thermal_distribution_resolution',max(x['L1_residual_per_event'] for x in eq)<.01,
        rows=eq,meaning='finite-grid residual per collision; not exact discrete detailed balance or a bound on every nonthermal history')
    q=json.loads((ROOT/'aligned_full_64.json').read_text())['snapshots'][-1]
    tight=json.loads((ROOT/'aligned_full_64_tight.json').read_text())['snapshots'][-1]
    tol=max(abs(tight[k]/q[k]-1) for k in ('Nplus','Nminus','gas_kinetic_mX'))
    add('ODE_tolerance_control_on_aligned_mesh',tol<1e-5,relative_change=tol)
    add('stationary_mass_BVP_and_initial_heat_remainder',max(x['BVP_residual'] for x in ph['stationary_mass_checks'])<3e-9,
        residual_energy_retained_mX=ph['source_energy_ledger']['unprocessed_core_and_matching_energy_mX'],
        frozen_profile_not_self_consistent=True)
    all_paths=list((ROOT/'code').glob('*.py'))
    for p in all_paths:ast.parse(p.read_text())
    finite=True
    for p in ROOT.glob('*.npz'):
        with np.load(p,allow_pickle=False) as d:
            for k in d.files:
                if np.issubdtype(d[k].dtype,np.number):finite &= bool(np.isfinite(d[k]).all())
    add('syntax_and_finite_saved_arrays',finite,python_files=len(all_paths))
    pins=['grut_m1_scattering_20261028/code/scattering.py','grut_m1_scattering_20261028/code/common.py',
          'grut_matched_candidate_20261027/PROFILE_M1.csv','grut_matched_candidate_20261027/CANDIDATE_M1.json',
          'grut_matched_candidate_20261027/DEFORMATION_RESULTS.json','grut_matched_candidate_20261027/code/deformation.py',
          'grut_m1_dynamics_20261027/code/radial.cpp','grut_m1_dynamics_20261027/EXCITATION_DATA.npz',
          'grut_m1_dynamics_20261027/EXCITATION_RESULTS.json']
    write_json('INPUT_PINS.json',{'parent_commit':'bf8ab3b96c07691410145c3b8c21ca18b009d66d',
        'sha256':{p:sha(RF/p) for p in pins}})
    write_json('ENVIRONMENT.json',{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'matplotlib':matplotlib.__version__,'platform':platform.platform(),
        'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0]})
    result={'verification_kind':'same-author, method/normalization/control checks',
        'checks':checks,'all_checks_pass':all(x['pass'] for x in checks),
        'selection_of_numerical_tolerances':'engineering checks, not preregistered empirical evidence',
        'independent_author_review':False,'complete_cosmological_cleanup':False}
    write_json('VERIFICATION.json',result)
    print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
