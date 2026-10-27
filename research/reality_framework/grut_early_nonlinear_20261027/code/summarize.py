"""Consolidate the regulated pilot and retain absent physical outputs."""
from en_common import *

def main():
    runs=read('RUNS.json')['cases'];fields=read('FIELD_RESULTS.json')['cases'];groups=read('ENGINE_CHECKS.json')['groups']
    qe=max(r['max_recorded_charge_relative_drift'] for r in runs.values());ee=max(r['max_recorded_energy_relative_drift'] for r in runs.values())
    dt=max(r['maximum_recorded_tau_spacing'] for r in runs.values())
    groups['nonlinear_conservation_and_dense_clock']=dict(passed=qe<1e-10 and ee<1e-5 and dt<.05000001,
        charge_relative=qe,energy_relative=ee,maximum_tau_spacing=dt,
        halfstep_energy_deviation=runs['halfstep48']['max_recorded_energy_relative_drift'])
    ex=read('EXPORT_CHECKS.json');groups['field_covariance_export']=dict(passed=ex['passed'],
        maximum_charge_identity=max(r['charge_relative'] for r in ex['rows']),maximum_energy_identity=max(r['energy_relative'] for r in ex['rows']))
    work=max(r['scalar_exchange_work_identity_relative'] for r in fields.values())
    groups['scalar_exchange_work_ledger']=dict(passed=work<1e-4,maximum_relative_residual=work,
        scope='Trapezoid integral on recorded dense conformal times; active X/phi exchange only, not the omitted heavy/gauge/quantum sectors.')
    save('VERIFICATION_RESULTS.json',dict(groups=groups,independent_reviewer_ran=False,physical_f_sym=None))
    assert all(g['passed'] for g in groups.values())
    control=read('LINEAR_PHI_CONTROL.json');final={}
    det=control['transfer_determinant_error']
    groups['full_scalar_and_phi_linear_preparation']['linear_phi_control_determinant']=det
    groups['full_scalar_and_phi_linear_preparation']['passed'] &= det<1e-7
    save('VERIFICATION_RESULTS.json',dict(groups=groups,independent_reviewer_ran=False,physical_f_sym=None))
    assert all(g['passed'] for g in groups.values())
    for tag,f in fields.items():
        final[tag]=dict(first_charge_contrast_one=f['thresholds']['first_charge_contrast_one'],
            first_X_variance_half=f['thresholds']['first_X_variance_half'],
            Ephi=f['final']['Ephi_fraction'],Ephi_plus_mix=f['final']['phi_plus_interaction_energy_fraction'],
            energy_outside_spatial_mean=f['final']['energy_outside_spatial_mean_fraction'],
            spatial_mean_charge_over_net=f['final']['spatial_mean_charge_fraction'],
            boundary_shell_075=f['spectra']['boundary_shell_075'],
            raw_nonzero_reference_nminus_over_net=[r['raw_nminus_over_net_current'] for r in f['spectra']['reference_projections']])
    save('RESULTS.json',dict(parent='b73de649843204639dfa1965a6326d4fc60afe06',review='ec1eb7f860f98f8d889f90ff94ebd53396f740bb',
        dense_sampling_corrected=True,matched_full_linear_fields=['X','Y','Sigma','phi'],nonlinear_fields=['X','phi'],
        nonlinear_histories=7,end_tau=800,final=final,
        homogeneous_linear_phi_control_final=control['rows'][-1]['phi_plus_interaction_raw_fraction'],
        physical_nbar_per_charge=None,physical_f_sym=None,physical_threshold=None,
        physical_initial_state_selected=False,full_gauged_model_evolved=False,
        scope='Regulated scalar-limit nonlinear field evolution and covariance export. No physical particle population, bound objects or cosmological point inferred.'))
    lines=['| Control | Phi energy + interaction | Energy outside spatial mean | Mean charge / net charge | Boundary-shell energy |',
           '|---|---:|---:|---:|---:|']
    for tag,x in final.items():
        lines.append(f"| {tag} | {100*x['Ephi_plus_mix']:.3f}% | {100*x['energy_outside_spatial_mean']:.2f}% | {x['spatial_mean_charge_over_net']:.4f} | {100*x['boundary_shell_075']:.3f}% |")
    (ROOT/'RESULT_TABLE.md').write_text('\n'.join(lines)+'\n')
    print('SUMMARY',final,flush=True)

if __name__=='__main__':main()
