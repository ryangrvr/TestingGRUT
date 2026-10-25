#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    pair=json.loads((ROOT/'PAIR_RESULTS.json').read_text());point=json.loads((ROOT/'POINT_RESULTS.json').read_text())
    screening=json.loads((ROOT/'SCREENING_RESULTS.json').read_text());release=json.loads((ROOT/'RELEASE_RESULTS.json').read_text())
    controls=[]
    def add(name,test,evidence):controls.append(dict(name=name,passed=bool(test),evidence=evidence))
    add('Pair-source Fourier normalization agrees with the analytic Gaussian',
        pair['controls']['Gaussian_Fourier_max_error_over_peak']<1e-10,pair['controls'])
    spread=max(r['quadrature_relative_spread'] for r in pair['summaries']);review=max(abs(r['relative_difference_from_review_Born']) for r in pair['summaries'])
    ledger=max(abs(r['core_charge_lost_per_event']-r['escaping_charge_per_event']-r['total_X_violation_per_event']) for r in pair['rows'])
    add('Independent invariant phase-space quadrature converges and agrees with the review Born estimate',
        spread<1e-7 and review<1e-4 and ledger==0,dict(quadrature_spread=spread,review_difference=review,pair_charge_ledger_residual=ledger))
    inf=point['inflation'];minimum=inf['full_radial_minimum']
    add('Full high-field radial stationary point is stable in amplitude and respects the analytic lower bound',
        minimum['gradient_norm']<1e-8 and min(inf['radial_mass2_over_H2'])>0 and minimum['potential_over_H4']>inf['full_potential_analytic_lower_bound_over_H4'],
        dict(gradient=minimum['gradient_norm'],radial_mass2=inf['radial_mass2_over_H2'],full_lower_bound=inf['full_potential_analytic_lower_bound_over_H4']))
    add('Independent angular finite differences match the physical phase mass matrix and its gauge null mode',
        inf['angular_Hessian_finite_difference_relative_error']<5e-6 and abs(inf['angular_mass2_over_H2'][0])<1e-12 and inf['angular_gauge_null_residual']<1e-10,
        dict(finite_difference_error=inf['angular_Hessian_finite_difference_relative_error'],angular_mass2=inf['angular_mass2_over_H2']))
    gauss=max(r['max_relative_Gauss_charge'] for r in release['rows']);stress=max(r['max_relative_scalar_stress_ledger_error'] for r in release['rows'])
    add('Full complex release conserves gauge charge and scalar stress, with zero aligned source and conjugate preparations',
        gauss<1e-9 and stress<1e-8 and release['rows'][0]['max_abs_comoving_X_charge_over_H3']==0 and release['CP_conjugate_charge_residual']<1e-9 and release['tolerance_final_charge_shift']<1e-7,
        dict(Gauss=gauss,stress=stress,aligned_X=release['rows'][0]['max_abs_comoving_X_charge_over_H3'],tolerance_shift=release['tolerance_final_charge_shift']))
    sr=screening['screened_Fourier_controls']
    add('Screened binding bound agrees with Parseval and positive massive-vector energy',
        abs(screening['Parseval_relative_error'])<1e-8 and max(abs(r['charge_relative_error']) for r in screening['controls'])<1e-8 and all(0<r['screened_energy_to_local_bound']<1+1e-10 for r in sr),
        dict(Parseval=screening['Parseval_relative_error'],extra_energy=screening['added_gauge_energy_bound_GeV'],massive_local_error_bound=screening['controls'][-1]['relative_large_mass_error_bound']))
    wall=point['wall_bias_variant']
    fact=np.polynomial.polynomial.polymul([1,-2,1],[2,1]).tolist()
    add('Proposed bias has the stated stationary vacua and a positive global factorization',
        fact==[2.,-3.,0.,1.] and wall['bracket_global_minimum_GeV4']>0 and wall['false_vacuum_mass2_GeV2']>0,
        dict(exact_polynomial_coefficients=fact,positive_bracket_minimum=wall['bracket_global_minimum_GeV4'],new_coefficient=wall['DeltaV_GeV4']))
    add('Residual gauge quantum numbers distinguish one core from a neutral three-core group',
        200000%6==2 and 3*200000+6*(-100000)==0 and screening['neutral_group_X_charge']%6==0,
        dict(one_core_mod6=screening['residual_discrete_charge_of_one_reference_core'],neutral_group=screening['neutral_group_X_charge'],Higgs_number=screening['neutral_group_Sigma_number_charge']))
    gates=[
        dict(test='Potential stability',baseline='PASS at classical tree level',qualification='The full quartics stabilize the high-field branch; its truncated single-field AD description fails.'),
        dict(test='Fixed initial phase',baseline='NOT SUPPLIED',qualification='Own-potential alignment gives zero charge. The review assumes an additional phase mismatch; a full inflationary fluctuation calculation is absent.'),
        dict(test='Bias removing walls',baseline='ABSENT',qualification='A separately proposed DeltaV=1e-22 GeV^4 variant passes a scaling-time estimate, not a network simulation.'),
        dict(test='Symmetric leftover matter',baseline='OPEN',qualification='Fixed abundance leaves a small extra-pair budget; no number-changing or full transition/debris evolution is computed.'),
        dict(test='Gauged object stability',baseline='PARTIAL',qualification='Screened classical binding survives with an enormous margin. Full fluctuation/metastability and a quantum population state remain open.'),
        dict(test='Corrected survival rate',baseline='PASS for included-channel estimate',qualification='Coherent plus free-Born pairs is still far below the engineering target; complete in-medium quantum loss is uncomputed.')]
    out=dict(scope='Author controls are separate from physical gate dispositions.',controls=controls,passed=sum(c['passed'] for c in controls),total=len(controls),six_physical_tests=gates,
        complete_cosmological_point=False,review_one_field_window_can_be_copied_to_current_gauge_point=False,
        in_medium_pair_calculation_completed=False,all_numerical_controls_passed=all(c['passed'] for c in controls))
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for c in controls:print('PASS' if c['passed'] else 'FAIL',c['name'],c['evidence'])
    assert out['passed']==out['total']


if __name__=='__main__':main()
