"""Mathematical checks of completed outputs; robustness values are reported separately."""
from pt_common import *
from pt_production import CASES

def main():
    groups={};engine=read('ENGINE_CHECKS.json');groups['discrete_engine_and_tangents']=dict(passed=engine['passed'],evidence=engine)
    prep=read('PREPARATION_RESULTS.json');canonical=0.
    J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])
    for L in [12,24]:
        d=np.load(ROOT/f'PREPARATION_L{L}.npz');lam=d['initial_mass_eigenvalues'];U=d['initial_mass_eigenvectors']
        for width in [1,3]:
          for k in d['ks']:
            w2=k*k+lam;w=np.sqrt(np.where(w2>0,w2,k*k+width*width))
            V=np.zeros((12,12));V[:6,:6]=(U/w)@U.T/2;V[6:,6:]=(U*w)@U.T/2
            canonical=max(canonical,float(np.max(abs(V@J@V-J/4))))
    sym=max(c['symplectic_absolute'] for row in prep['rows'] for c in row['checks']);det=max(c['phi_determinant_error'] for row in prep['rows'] for c in row['checks'])
    groups['Gaussian_and_transfer_canonical']=dict(passed=canonical<1e-11 and sym<1e-8 and det<1e-8,initial_pure_covariance_error=canonical,scalar_transfer_symplectic=sym,phi_determinant=det)
    result=read('RESULTS.json');con=result['conservation']
    groups['nonlinear_conservation']=dict(passed=con['maximum_energy_relative_drift']<1e-6 and con['maximum_charge_relative_drift']<1e-10,evidence=con)
    groups['averaged_export_reconstruction']=dict(passed=con['averaged_export_energy_relative_error']<1e-12 and con['averaged_export_charge_relative_error']<1e-12,
        energy_relative=con['averaged_export_energy_relative_error'],charge_relative=con['averaged_export_charge_relative_error'])
    groups['time_refinement_of_action']=dict(passed=result['time_refinement_maximum_total_action_absolute']<1e-5,maximum_action_absolute_difference=result['time_refinement_maximum_total_action_absolute'],
        scope='A numerical check of the stated classical statistic, not a particle-definition or state-selection test.')
    co=read('COHERENT_MATCHING.json');slow=[r for r in co['rows'] if r['duration']==80]
    charge=max(r['charge_relative_error'] for r in co['rows'] if r['charge_relative_error'] is not None)
    work=max(r['work_identity_relative'] for r in co['rows'] if r['work_identity_relative'] is not None)
    actionerr=max(abs(r['free_nminus_per_charge']-co['radial_action_per_charge']) for r in slow)
    fixed=max(abs(r['final_radial_action_per_charge']-co['radial_action_per_charge']) for r in co['fixed_coupling_mass_ramp_control'])
    review=Path(REPO.parent/'review_fd71d3a/tracks/en_particles/skeptic/results/sk_orbit_tau20.json')
    reference=json.loads(review.read_text());reproduction=abs(co['lab_negative_action_per_original_net_charge']-reference['n_minus_per_Qtotal'])
    groups['coherent_action_and_matching']=dict(passed=co['lab_charge_reconstruction_relative']<1e-9 and charge<1e-7 and work<1e-7 and actionerr<5e-5 and fixed<2e-6 and reproduction<1e-8,
        lab_charge_identity=co['lab_charge_reconstruction_relative'],review_action_reproduction_absolute=reproduction,charge_relative=charge,work_relative=work,
        slow_free_action_maximum_difference=actionerr,fixed_quartic_radial_action_maximum_difference=fixed,
        scope='Controlled homogeneous matching histories, not a cosmological transition.')
    free=read('RESPONSE_FREE_CHECK.json');groups['matrix_response_free_control']=dict(passed=free['passed'] and free['maximum_cross_response']==0 and free['maximum_anomalous_response']==0,evidence=free)
    groups['nonlinear_response_canonical_transform']=dict(passed=result['response']['maximum_canonical_sum_deviation']<.03,
        maximum_canonical_sum_deviation=result['response']['maximum_canonical_sum_deviation'],
        scope='Finite-window/probe canonical transform check; does not certify poles or populations.')
    assert len(groups)==8
    out=dict(passed=all(g['passed'] for g in groups.values()),groups=groups,
        infrared_robustness_observations=result['negative_action_ranges_by_window'],IR_spatial_refinement_observations=result['IR_spatial_refinement'],
        independent_reviewer_of_this_packet=False,physical_particle_population_certified=False)
    save('VERIFICATION_RESULTS.json',out);print(json.dumps(out,indent=2));assert out['passed']

if __name__=='__main__':main()
