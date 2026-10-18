#!/usr/bin/env python3
"""Check the claimed leading-regime evidence, including independent sampled ODEs.

This checks a numerical implementation and its scope-specific acceptance rules;
it is not a proof of the proposed physics or a quantum/global stability test.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parent


def main():
    scalar = json.loads((ROOT/'results.json').read_text())
    reporter = json.loads((ROOT/'reporter_results.json').read_text())
    checks = []

    def require(name, passed, evidence):
        checks.append({'name': name, 'passed': bool(passed), 'evidence': evidence})

    obs = scalar['observables']; scheck = scalar['checks']
    require('four_converged_bvp_runs', len(scalar['attempts']) == 4 and
            all(a['status'] == 0 for a in scalar['attempts']),
            [a['status'] for a in scalar['attempts']])
    require('final_collocation_residual', scalar['attempts'][-1]['max_rms_residual'] < 1.1e-9,
            scalar['attempts'][-1]['max_rms_residual'])
    prev = scalar['convergence'][-2]
    for key in ['Qbar', 'Ebar']:
        difference = abs(obs[key]-prev[key])/abs(obs[key])
        require('last_domain_relative_'+key, difference < 1e-6, difference)
    require('bound_against_all_free_carriers', obs['Ebar_over_Qbar'] < 1., obs['Ebar_over_Qbar'])
    require('stress_conservation', scheck['stress_conservation_max_abs'] < 1e-7,
            scheck['stress_conservation_max_abs'])
    require('virial_identity', abs(scheck['virial_relative_to_E']) < 1e-7,
            scheck['virial_relative_to_E'])
    require('active_mass_equals_energy', abs(scheck['active_mass_over_energy']-1) < 1e-7,
            scheck['active_mass_over_energy'])
    require('weak_gravity', max(obs['Schwarzschild_radius_over_R90'],
            abs(scheck['central_Phi']), abs(scheck['central_Psi'])) < 1e-4,
            {'compactness':obs['Schwarzschild_radius_over_R90'],
             'Phi':scheck['central_Phi'],'Psi':scheck['central_Psi']})
    require('exterior_lensing_correspondence',
            abs(scheck['light_bending'][1]['ratio_to_point_mass']-1) < 1e-7,
            scheck['light_bending'][1]['ratio_to_point_mass'])

    # Recheck sampled field equations with finite differences, not BVP residuals.
    d = np.load(ROOT/'results_profile.npz')
    r=d['rho']; F,Fp,S,Sp=d['y']; omega=float(d['Omega'])
    a=scalar['configuration']['a']; b=scalar['configuration']['b']
    mask=(r>.02)&(r<40.)
    rf=np.gradient(Fp,r,edge_order=2)+np.divide(2*Fp,r,out=np.zeros_like(r),where=r>0)-(S*S+b*F*F-omega*omega)*F
    rs=np.gradient(Sp,r,edge_order=2)+np.divide(2*Sp,r,out=np.zeros_like(r),where=r>0)-(a*(S*S-1)+F*F)*S
    residual=max(float(np.max(abs(rf[mask]))),float(np.max(abs(rs[mask]))))
    require('independent_sampled_radial_equations', residual < 1e-6, residual)
    require('profile_is_nontrivial_and_admissible', F[0] > 1. and .0 < S[0] < .1 and
            float(F.min()) > -1e-7 and float(S.min()) > 0.,
            {'F0':float(F[0]),'S0':float(S[0]),'min_F':float(F.min())})
    q_sample=4*np.pi*omega*simpson(r*r*F*F,x=r)
    require('sampled_profile_charge_matches_record', abs(q_sample/obs['Qbar']-1) < 1e-8,
            float(q_sample/obs['Qbar']))

    spectra=scalar['spectra']
    finest=[x for x in spectra if x['N']==3200]
    require('positive_tested_nonsymmetry_amplitude_sectors',
            all(x['lowest_eigenvalues'][0]>0 for x in finest if x['ell'] in [0,2,3]),
            {str(x['ell']):x['lowest_eigenvalues'][0] for x in finest})
    trans=[x for x in spectra if x['ell']==1]
    scaled=[abs(x['lowest_eigenvalues'][0])/x['h']**2 for x in trans]
    require('translation_mode_overlap_and_refinement',
            all(x['lowest_translation_overlap_abs']>.9999 for x in trans) and
            max(scaled)/min(scaled)<1.05,
            {'absolute_eigenvalue_over_h_squared':scaled})
    phase=[x for x in scheck['phase_sector_spectrum'] if x['ell']==0]
    scaled_phase=[abs(x['lowest_eigenvalues'][0])/x['h']**2 for x in phase]
    require('phase_zero_mode_refines_quadratically',max(scaled_phase)/min(scaled_phase)<1.05,
            scaled_phase)

    for x in reporter['phase_results']:
        k=x['k_over_gv']
        require(f'phase_two_method_agreement_k_{k}',x['convergence'][-1]['method_difference'] < 1e-7,
                x['convergence'][-1]['method_difference'])
        require(f'analytic_born_bound_k_{k}',x['actual_born_difference'] <
                x['analytic_born_remainder_bound_for_exact_profile'],
                {'error':x['actual_born_difference'],
                 'bound':x['analytic_born_remainder_bound_for_exact_profile']})
    require('packet_norm_conservation',reporter['maximum_packet_norm_error'] < 1e-10,
            reporter['maximum_packet_norm_error'])
    probability_error=0.
    for run in reporter['packet_grid_runs']:
        for x in run['trajectory']:
            p=x['selected_probability']; plus=x['joint_plus_probability']; minus=x['joint_minus_probability']
            probability_error=max(probability_error,abs(plus+minus-p))
            assert -.0000000001<=plus<=p+1e-10 and -.0000000001<=minus<=p+1e-10
            assert 0.<=x['conditional_plus_probability']<=1. and x['visibility']<=1.+1e-10
    require('selected_packet_probability_sums',probability_error<1e-12,probability_error)
    final_packets=[x['trajectory'][-1] for x in reporter['packet_grid_runs']]
    require('packet_refinement_improves_conditional_output',
            abs(final_packets[2]['conditional_plus_probability']-final_packets[1]['conditional_plus_probability']) <
            abs(final_packets[1]['conditional_plus_probability']-final_packets[0]['conditional_plus_probability']),
            [x['conditional_plus_probability'] for x in final_packets])
    report={'scope':'internal leading-regime numerical evidence only',
            'passed':all(x['passed'] for x in checks),'checks':checks,
            'unverified':['global fission and nonlinear stability','quantum lifetime',
                          'microscopic apparatus','empirical realization','UV completion']}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'number_of_checks':len(checks),
                      'failed':[x['name'] for x in checks if not x['passed']]},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
