#!/usr/bin/env python3
"""Cross-check exact RG evidence, independent numerical routes and claim limits."""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    checks = []
    def check(name, condition, detail=None):
        checks.append({'name': name, 'passed': bool(condition), 'detail': detail})
    rg = json.loads((HERE/'rg_results.json').read_text())
    pot = json.loads((HERE/'potential_results.json').read_text())
    flow = json.loads((HERE/'flow_results.json').read_text())
    rmod = module('grut_rg', 'rg_symbolic.py')
    pmod = module('grut_potential', 'effective_potential.py')
    fmod = module('grut_flow', 'flow_check.py')
    check('all exact algebraic checks passed', rg['checks_passed'] == 290 and all(c['passed'] for c in rg['checks']))
    check('full polynomial identity exactly zero', rg['full_identity_residual_monomials'] == 0)
    check('RG protocol identity', rg['protocol_sha256'] == hashlib.sha256((HERE/'rg_protocol.json').read_bytes()).hexdigest())
    check('RG source identity', rg['implementation_sha256'] == hashlib.sha256((HERE/'rg_symbolic.py').read_bytes()).hexdigest())
    check('flow protocol identity', flow['protocol_sha256'] == hashlib.sha256((HERE/'flow_protocol.json').read_bytes()).hexdigest())
    check('flow acceptance conditions passed', flow['all_checks_pass'] and all(flow['checks'].values()))
    check('potential matching agrees with flow matching', np.array_equal([pmod.PARAM[k] for k in fmod.NAMES], fmod.THETA0))
    translate = dict(zip(fmod.NAMES, ['lp','lx','c','d','ld','alpha','mp','mx','md','C']))
    rg_beta = np.array([rg['matching_beta'][translate[k]]['beta'] for k in fmod.NAMES])
    check('symbolic and numerical matching betas agree', np.max(abs(fmod.beta(fmod.THETA0)-rg_beta)) < 1e-18)
    check('generated alpha portal nonzero', rg['matching_beta']['alpha']['beta'] > 0)
    check('generated X mass coefficient nonzero', rg['matching_beta']['mx']['beta'] < 0)
    check('generated vacuum constant nonzero', rg['matching_beta']['C']['beta'] > 0)
    check('omitting each generated operator breaks identity', all(len(v['terms']) > 0 for v in rg['incorrect_omission_residuals'].values()))
    # General parameter values (including nonzero alpha and independent masses)
    # test full Hessian normalization and beta routines beyond matching inputs.
    th = np.array([.013,.002,.011,.0002,.0015,.00007,-.012,.0002,.008,.003])
    sym_sub = {translate[k]: rmod.F(str(value)) for k,value in zip(fmod.NAMES, th)}
    pol = json.loads((HERE/'rg_polynomials.json').read_text())
    for name in fmod.NAMES:
        expr = rmod.Poly({tuple(term['powers'].get(n, 0) for n in rmod.NAMES):rmod.F(term['coefficient'])
                          for term in pol['beta_B'][translate[name]]})
        expected = float(expr.substitute(sym_sub).constant())/(16*np.pi**2)
        i = fmod.NAMES.index(name)
        check('general-parameter beta '+name, abs(fmod.beta(th)[i]-expected) < 1e-18)
    herror = cerror = 0.
    for q in json.loads((HERE/'flow_protocol.json').read_text())['backgrounds_phi_x_y_D']:
        q = np.asarray(q)
        herror = max(herror, float(np.max(abs(pmod.tree_hessian(q)-fmod.mass_matrix(q, fmod.THETA0)))))
        cerror = max(cerror, abs(pmod.cw(q)[0].real-fmod.cw(q,fmod.THETA0,.1)))
    check('independent Hessian implementations agree', herror < 1e-14, herror)
    check('independent CW eigenvalue implementations agree', cerror < 1e-16, cerror)
    vac = pot['vacuum']
    bulk = pot['bulk']
    check('strict tadpole equation', np.max(abs(np.array(vac['CW_tadpole_at_tree_point'])+pmod.tree_hessian(np.array([1.,0,0,0])) @ np.array(vac['strict_one_loop_field_shift']))) < 1e-18)
    check('truncated vacuum is stationary', np.max(np.abs(vac['truncated_vacuum_gradient'])) < 1e-13)
    check('vacuum curvatures remain positive', min(vac['strict_one_loop_curvature_masses_squared']) > 0)
    check('strict vacuum shift small at benchmark', abs(vac['strict_one_loop_field_shift'][0]) < 1e-3)
    check('bulk strict and truncated calculations consistent', abs(bulk['strict_one_loop_ratio_mass_squared']-bulk['truncated_minimum_ratio_mass_squared']) < 1e-10)
    check('bulk all three bounded starts converge', all(o['success'] and abs(o['value']-bulk['truncated_minimum_ratio_mass_squared']) < 1e-12 for o in bulk['two_dimensional_optimizer_checks']))
    check('bulk rectangle tree Hessian lower bound positive', bulk['rectangle_tree_convexity_certificate']['determinant_lower'] > 0 and bulk['rectangle_tree_convexity_certificate']['phase_mass_lower'] > 0)
    check('bulk tree mass metadata evaluated at tree point', np.max(abs(np.array(bulk['tree_background_masses_squared'])-np.linalg.eigvalsh(fmod.mass_matrix([0.,10**.25,0,0],fmod.THETA0)))) < 1e-15)
    best_errors = []
    for point in pmod.PROTOCOL['derivatives']['points']:
        relevant = [r for r in pot['derivative_checks'] if r['point'] == point]
        best_errors.append(min(r['CW_hessian_from_gradients_max_abs_error'] for r in relevant))
    check('analytic CW derivatives checked by finite differences', max(best_errors) < 1e-12, best_errors)
    profile = np.genfromtxt(HERE/'potential_profile_flags.csv',delimiter=',',names=True)
    check('wall negative modes preserved in raw output', np.sum(profile['negative_mode_count'] > 0) == pot['classical_profile_flags']['negative_Hessian_points'] > 0)
    check('complex CW wall diagnostic preserved', np.max(abs(profile['CW_imaginary'])) > 0)
    check('no invalid quantum profile-energy integration', pot['classical_profile_flags']['CW_profile_energy_integral_performed'] is False)
    check('relative vacuum subtraction not gravitational calibration', pot['gravity_vacuum_shift']['gravity_calibration_applied'] is False)
    check('vacuum calibration correction has opposite sign to raw energy', abs(pot['gravity_vacuum_shift']['required_constant_shift_for_zero_vacuum_strict_over_v4']+vac['strict_one_loop_vacuum_energy']) < 1e-20)
    reference = HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz'
    check('classical profile input unchanged', hashlib.sha256(reference.read_bytes()).hexdigest() == '0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3')
    result = {'checks_passed':sum(c['passed'] for c in checks), 'total_checks':len(checks),
              'all_checks_pass':all(c['passed'] for c in checks), 'checks':checks,
              'scope':'Internal algebra, numeric and scope checks. No independent validation or quantum soliton stability follows.'}
    (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks_passed':result['checks_passed'],'total_checks':result['total_checks'],'all_checks_pass':result['all_checks_pass']}))
    for c in checks:
        if not c['passed']:
            print(json.dumps(c))
    if not result['all_checks_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
