"""Numerical checks separated from the unresolved physical population verdict."""
from common import *
from analyze import positive_energy
from analysis_io import read_analysis

def main():
    engine=json.loads((ROOT/'ENGINE_CHECKS.json').read_text())
    radial=json.loads((ROOT/'RADIAL_CHECKS.json').read_text())
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases']
    result=json.loads((ROOT/'RESULTS.json').read_text())
    data=read_analysis()['cases']
    bulk=json.loads((ROOT/'BULK_RESULTS.json').read_text())
    groups=[]
    def group(name,ok,**details):groups.append({'name':name,'pass':bool(ok),**details})
    group('changed_3D_force_and_independent_evolution',engine['pass'],details=engine)
    emax=max(v['energy_relative_drift'] for v in runs.values());qmax=max(v['charge_relative_drift'] for v in runs.values())
    expected_hash=hashlib.sha256((ROOT/'code/field3d.cpp').read_bytes()+(ROOT/'code/common.py').read_bytes()).hexdigest()
    source_identity=all(v['code_sha256']==expected_hash for v in runs.values())
    point=json.loads((PARENT/'CANDIDATE_M1.json').read_text())
    fixed_inputs=(G==point['light_renormalized_targets']['gX'] and LX==point['light_renormalized_targets']['lambdaX']
                  and LP==point['light_renormalized_targets']['lambdaPhi'] and EPS==point['binding_potential_deformation']['epsilon'])
    group('all_3D_conservation',emax<1e-4 and qmax<1e-10 and source_identity and fixed_inputs,
          max_energy_drift=emax,max_charge_drift=qmax,histories=len(runs),source_hash_matched=source_identity,M1_couplings_unchanged=fixed_inputs)
    half=result['numerical_comparisons']['cold64_vs_cold64_half']
    ratio=runs['cold64']['energy_relative_drift']/runs['cold64_half']['energy_relative_drift']
    group('cold_time_refinement',half['final_field_RMS_relative_difference']<.001 and 14<ratio<18,
          field_RMS_difference=half['final_field_RMS_relative_difference'],energy_drift_ratio=ratio)
    core=data['object96_perturbed']['rows'][-1]['central_spheres'][1]
    oldcore=data['object64_perturbed']['rows'][-1]['central_spheres'][1]
    group('supplied_object_survival_and_grid_control',core['E_over_Q']<1 and core['Q']>199000 and abs(core['E_over_Q']/oldcore['E_over_Q']-1)<.001,
          fine_core=core,coarse_core=oldcore)
    roots=np.array(bulk['stationary_quartic_roots'])
    real_positive=roots[(abs(roots[:,1])<1e-10)&(roots[:,0]>0),0]
    group('bulk_minimum_and_bound',len(real_positive)==1 and 0<real_positive[0]<1 and bulk['direct_three_equation_residual']<1e-12 and bulk['sampled_potential_bound_min_margin']>0,
          bulk=bulk['bulk'],positive_stationary_roots=real_positive.tolist())
    group('independent_radial_evolution_and_receiver',radial['pass'],details=radial)
    ext=json.loads((ROOT/'EXTERIOR_RESULTS.json').read_text())['cases']['object96_perturbed']['rows']
    last=[v for v in ext if v['t']>=160 and v['cut_radius']>=16]
    cerr=max(v['Fourier_real_space_charge_relative_error'] for v in ext)
    group('exterior_charge_identity_and_weak_interaction',cerr<1e-12 and max(v['absolute_interaction_over_free_energy'] for v in last)<.001,
          max_charge_identity_error=cerr,interaction_max=max(v['absolute_interaction_over_free_energy'] for v in last))
    # Reconstruct final scalar energy/current independently from the full stored fields.
    reconstruction=[]
    for tag in ('cold64','dilute64','ellipse96','object96_perturbed'):
        N,L,t,a=snapshot(runs[tag]['snapshots'][-1]);q=np.asarray(a[:3]);p=np.asarray(a[3:])
        rho=positive_energy(q,p,L).mean();nq=np.mean(q[0]*p[1]-q[1]*p[0])
        reconstruction.append({'tag':tag,'energy':float(abs(rho/runs[tag]['final']['rho']-1)),
                               'charge':float(abs(nq/runs[tag]['final']['nQ']-1)),
                               'positive_vs_signed_gradient':float(abs(rho/local_energy(q,p,L).mean()-1))})
    group('archived_full_field_energy_charge',max(max(v[k] for k in ('energy','charge','positive_vs_signed_gradient')) for v in reconstruction)<1e-10,
          reconstruction=reconstruction)
    physical={'actual_cosmic_eta':None,'actual_cosmic_pair_ratio':None,'actual_final_abundance':None,
              'formed_isolated_object_population_certified':False,
              'hot_spatial_convergence_certified':False,
              'vacuum_density_predicted':False,
              'reason':'finite local scalar controls; source/expansion/quantum chemistry and core-wave matching remain'}
    out={'groups':groups,'passed':sum(v['pass'] for v in groups),'total':len(groups),
         'independent_author_review':False,'physical_status':physical}
    save('VERIFICATION.json',out)
    print(json.dumps({'passed':out['passed'],'total':out['total'],'physical_status':physical},indent=2))
    if out['passed']!=out['total']:raise RuntimeError('Numerical verification failed; inspect groups')

if __name__=='__main__':main()
