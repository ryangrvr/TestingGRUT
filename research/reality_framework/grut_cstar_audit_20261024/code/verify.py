#!/usr/bin/env python3
import json
import numpy as np
from closure_checks import ROOT


def main():
    closure=json.loads((ROOT/'CLOSURE_CHECKS.json').read_text());rate=json.loads((ROOT/'CHARGE_LEAKAGE_RESULTS.json').read_text())
    checks=[]
    def check(name,ok,evidence):checks.append(dict(name=name,passed=bool(ok),evidence=evidence))
    check('Reduced-Planck de Sitter identity and compact positive-vacuum GR counterexample',closure['deSitter']['identity_relative_error']<1e-12 and closure['compact_GR_counterexample']['Friedmann_identity_max_error']<1e-12,
        dict(entropy_identity=closure['deSitter']['identity_relative_error'],Friedmann=closure['compact_GR_counterexample']['Friedmann_identity_max_error']))
    clock=closure['finite_closed_clock'];values=[r['conditional_sigma_x'] for r in clock['relational_clock_values']]
    check('Finite constrained clock has nontrivial correlations and an exact U1',clock['physical_H_constraint_dimension']==2 and clock['H_state_residual']<1e-12 and clock['exact_U1_commutator_norm']<1e-12 and np.max(abs(np.array(values)-[1,0,-1]))<1e-12,values)
    torque=max(r['torque_finite_difference_relative_error'] for r in closure['explicit_breaking_Ward_controls'])
    check('Charge-breaking Ward torque agrees with an independent angular finite difference',torque<1e-7,torque)
    scan=closure['operator_scan'];z5=scan['residual_Z5'];z6=scan['residual_Z6']
    check('Scalar operator scan identifies relevant breaking and residual-gauge selection rules',min(r['dimension'] for r in scan['without_X_protection'])==1 and len(z5)==len(z6)==1 and z5[0]['dimension']==5 and z6[0]['dimension']==6,dict(Z5=z5,Z6=z6))
    free=max(r['free_amplitude_relative_error'] for r in rate['free_controls']);other=max(r['free_other_channel_amplitude_max'] for r in rate['free_controls'])
    check('Forced free radiation matches an independent radial Green-function integral',free<1e-6 and other<1e-12,dict(relative_amplitude_error=free,other_channel_amplitude=other))
    work=max(r['source_work_relative_residual'] for r in rate['rows']+rate['free_controls']+rate['background_controls'])
    amplitude=max(r['maximum_dimensionless_field_response_bound_for_c1'] for r in rate['rows'])
    check('Outgoing flux matches source work and unit-coefficient perturbations remain weak',work<1e-6 and amplitude<1e-10,dict(max_work_residual=work,max_dimensionless_field_response_bound=amplitude))
    spread={};profile={}
    for power in [5,6]:
        rows=[r for r in rate['rows'] if r['power']==power];values=[r['initial_fractional_core_charge_loss_rate_per_second_for_c1'] for r in rows]
        spread[str(power)]=max(values)/min(values)-1
        coarse=next(r for r in rate['background_controls'] if r['power']==power)
        profile[str(power)]=abs(coarse['initial_fractional_core_charge_loss_rate_per_second_for_c1']/values[-1]-1)
    check('Forced rates are stable across domains/meshes/tolerances and two profile tables',max(spread.values())<1e-5 and max(profile.values())<1e-4,dict(solver_spreads=spread,profile_shifts=profile))
    errors=[]
    for r in rate['rows']:
        rates=np.array(r['number_rates_per_second_for_c1']);n=r['power'];nu=r['nu'];omega=nu/n
        core=np.dot(np.array(r['frequencies_over_mX'])/omega,rates)
        global_loss=core-(rates[0]-rates[1])
        errors.extend([abs(core/r['core_charge_loss_per_second_for_c1']-1),abs(global_loss/r['global_charge_loss_per_second_for_c1']-1),abs(r['max_abs_coefficient_for_one_percent_initial_rate_age_target']**2*r['frozen_profile_linear_charge_loss_fraction_over_reference_age_for_c1']/.01-1)])
    check('Emitted energy/charge and weak-loss coefficient conditions agree',max(errors)<1e-12,max(errors))
    best={str(n):next(r for r in reversed(rate['rows']) if r['power']==n) for n in [5,6]}
    out=dict(scope='Author mathematical and numerical controls of conditional C* claims and explicit EFT deformations, not verification of C* or a UV quantum-gravity completion.',
        checks=checks,passed=sum(r['passed'] for r in checks),total=len(checks),
        unit_coefficient_Planck_probe_passes_one_percent_reference_age_target={n:r['frozen_profile_linear_charge_loss_fraction_over_reference_age_for_c1']<=.01 for n,r in best.items()},
        independently_selected_vacuum_density=False,observer_encoding_constructed_for_GRUT=False,X_charge_UV_protection_constructed=False,quantum_gravity_Wilson_coefficients_predicted=False,
        independent_Claude_review_claimed=False)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in checks:print('PASS' if r['passed'] else 'FAIL',r['name'],r['evidence'])
    print('SPECIFIED_OPERATOR_BENCHMARK',out['unit_coefficient_Planck_probe_passes_one_percent_reference_age_target'])
    assert all(r['passed'] for r in checks)


if __name__=='__main__':main()
