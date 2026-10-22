"""Independent formulation checks of the 78352e2 handoff, not an empirical test.

No imports from any audited track or from the earlier Codex implementation.
Constant-K runs evolve (H, rho_V, rho_r, N) in cosmic time, with Hdot
from Raychaudhuri. The closed N-law and receiver quadrature are comparisons.
KSV static states and the K22 example are derived directly below.
Run: python tracks/codex_review/code/review.py
"""
from pathlib import Path
import json
import math

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
import sympy as sp

HERE = Path(__file__).resolve().parents[1]
I = json.loads((HERE / 'INPUTS.json').read_text())
CHECKS = {}


def check(name, condition):
    CHECKS[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def radiation_run(v0, r0, gamma, n_end, rtol=2e-11):
    """M=1. Includes the receiving radiation in both Einstein equations."""
    def rhs(t, y):
        h, v, r, n = y
        source = 4 * gamma * abs(h) * v * v / 9
        return [-2 * r / 3, -source, -4 * h * r + source, h]

    def done(t, y):
        return y[3] - n_end
    done.terminal = True
    done.direction = 1
    sol = solve_ivp(rhs, (0, 2e5), [math.sqrt((v0 + r0) / 3), v0, r0, 0],
                    method='DOP853', rtol=rtol, atol=1e-14,
                    events=done, dense_output=True)
    if sol.status != 1:
        raise RuntimeError(sol.message)
    tt = np.linspace(0, sol.t[-1], 501)
    h, v, r, n = sol.sol(tt)
    exact = v0 / (1 + 4 * gamma * v0 * n / 9)
    # Retarded radiation receiver, evaluated without e^(4N) overflow.
    rend = quad(lambda u: math.exp(-4 * (n_end-u)) * 4 * gamma / 9 *
                (v0 / (1 + 4 * gamma * v0 * u / 9))**2,
                0, n_end, epsabs=1e-16, epsrel=1e-11)[0] + r0 * math.exp(-4*n_end)
    scale = 3*h*h + np.abs(v) + np.abs(r)
    return {
        'v_initial': v0, 'r_initial': r0, 'N_end': float(n[-1]),
        'v_final': float(v[-1]), 'r_final': float(r[-1]),
        'max_relative_closed_v_error': float(np.max(abs(v/exact-1))),
        'max_relative_constraint': float(np.max(abs(3*h*h-v-r)/scale)),
        'receiver_relative_error': float(abs(r[-1]/rend-1)),
    }


def main():
    out = {'author': I['author'], 'versions': {'numpy': np.__version__,
           'scipy': scipy.__version__, 'sympy': sp.__version__}}
    gamma = 1/(32*math.pi**2)
    h, v, r, Q = sp.symbols('H v r Q', real=True)
    C = 3*h*h - v-r
    cdot = sp.diff(C,h)*(-2*r/3)+sp.diff(C,v)*(-Q)+sp.diff(C,r)*(-4*h*r+Q)
    check('Einstein_constraint_with_radiation_receiver', sp.simplify(cdot) == 0)
    out['constraint_identity'] = str(sp.simplify(cdot))
    # The boundary condition, rather than the equations of motion, moves mu.
    q, dq, mu, dmu, e1 = sp.symbols('q dq mu dmu e1', nonzero=True)
    dp = (e1*dq + 1 - mu*dq-q*dmu).subs(e1,mu)
    check('Minkowski_state_reselection_dmu_dc_equals_inverse_q',
          sp.simplify(dp.subs(dmu,1/q)) == 0)
    out['constant_shift'] = {'fixed_mu': 'additive c gravitates',
                            'reselected_equilibrium': 'dmu0/dc = 1/q0'}

    runs = [radiation_run(v0,r0,gamma,I['dissipation_cases']['N_end'])
            for v0 in I['dissipation_cases']['v_initial']
            for r0 in I['dissipation_cases']['radiation_initial']]
    check('radiation_closed_flow_all_runs',max(x['max_relative_closed_v_error'] for x in runs)<1e-9)
    check('radiation_receiver_all_runs',max(x['receiver_relative_error'] for x in runs)<2e-8)
    check('constraint_all_runs',max(x['max_relative_constraint'] for x in runs)<1e-9)
    check('vacuum_N_flow_independent_of_radiation_density',
          max(abs(runs[k]['v_final']/runs[k+1]['v_final']-1) for k in [0,2])<1e-9)
    fine = radiation_run(1,0.2,gamma,12,rtol=2e-12)
    check('selected_tolerance_refinement',abs(fine['v_final']/runs[2]['v_final']-1)<1e-9)
    out['radiation_runs'] = runs
    out['refinement'] = fine

    M = I['M_GeV']; rho = I['rho_Lambda_reference_GeV4']
    B = M**4; rL = rho/B; a = 4*gamma/9
    n_std = math.log(M/I['T0_reference_GeV'])+math.log(I['gs_hot']/I['gs_cold'])/3
    envelope = 1/(a*n_std)
    # For finite positive v_i the often-called "floor" is an UPPER envelope.
    end_from_one = 1/(1+a*n_std)
    check('claimed_universal_floor_is_upper_envelope',0 < end_from_one < 1 < envelope)
    t = I['age_reference_Gyr']*1e9*I['seconds_per_year']/I['hbar_GeV_s']
    rho_age = (2*gamma*t/(3*math.sqrt(3)*M**5))**(-2/3)
    req_n = (1/rL-1)/a
    rho_bbn = math.pi**2/30*I['BBN_gstar']*I['BBN_temperature_GeV']**4
    n_bbn = math.log(I['BBN_temperature_GeV']/I['T0_reference_GeV'])+math.log(I['BBN_gs']/I['gs_cold'])/3
    bbn_relative_depletion_bound = a*(rho_bbn/B)*n_bbn
    out['finite_history'] = {
        'gamma':gamma, 'rho_Lambda_in_M4':rL,
        'N_standard_entropy_conserving_span':n_std,
        'upper_envelope_in_M4':envelope,
        'final_in_M4_if_initial_M4':end_from_one,
        'final_over_rho_Lambda_if_initial_M4':end_from_one/rL,
        'N_M4_to_reference_rho_Lambda':req_n,
        'early_start_attractor_at_reference_age_GeV4':rho_age,
        'early_start_attractor_over_reference_rho_Lambda':rho_age/rho,
        'BBN_reference_radiation_GeV4':rho_bbn,
        'N_BBN_to_present_reference':n_bbn,
        'max_fractional_self_drain_after_BBN_if_initial_v_le_rho_r':bbn_relative_depletion_bound,
        'scope':'Radiation-only exact flow, constant K, fixed gamma. Age asymptote assumes early-start vacuum dominance; entropy span is a reference history, not a proof for arbitrary histories.'}
    check('finite_history_matches_both_builds',abs(req_n/9.913512607512736e122-1)<1e-12 and abs(rho_age/rho/8.902768918001714e81-1)<1e-12)

    # Negative branch: do not divide by H at the turnaround.
    def neg_rhs(t,y):
        h,v,r,n=y; src=a*abs(h)*v*v
        return [-2*r/3,-src,-4*h*r+src,h]
    def stop_h(t,y): return y[0]
    stop_h.terminal=True; stop_h.direction=-1
    neg=solve_ivp(neg_rhs,(0,100),[math.sqrt(.01/3),-.01,.02,0],
                  events=stop_h,method='DOP853',rtol=1e-12,atol=1e-15)
    check('negative_vacuum_not_restored',neg.status==1 and neg.y[1,-1]<-.01)
    out['negative_control']={'turnaround_time_M_inverse':float(neg.t[-1]),
                             'N':float(neg.y[3,-1]),'v_final':float(neg.y[1,-1])}

    # KSV16 static family, K=f/3 in the source's dimensionless conventions.
    f=sp.symbols('f',positive=True)
    eps=-f*f/2+f**4/6
    hs=sp.diff(eps,f)-eps/f
    us=sp.diff(eps,f)-2*hs
    static_res=sp.simplify(f*hs-(eps-us*f))
    check('KSV16_has_continuum_of_static_de_Sitter_states',static_res==0 and hs.subs(f,sp.Rational(11,10))>0)
    out['KSV16_static_family']={'h_squared':str(sp.factor(hs)),
                              'u':str(sp.factor(us)),
                              'source_at_each_static_state':0,
                              'meaning':'Quadratic time-derivative dissipation cannot select a unique Minkowski endpoint.'}

    # K22 missing verifier's table example, separate energy variables.
    # G=1, K=1/(8*pi), source=g_tilde H^5; xi reconstructed at lambda=mu^2/4.
    P=I['K22']; coeff=8*math.pi/3
    vv=P['lambda']+P['xi_initial']**2-P['mu']*P['xi_initial']
    rr=P['h_initial']**2/coeff-vv
    def k_rhs(t,y):
        hv,rv,rad=y; src=P['g_tilde']*abs(hv)*hv**4
        return [-16*math.pi*rad/3,-src,-4*hv*rad+src]
    ks=solve_ivp(k_rhs,(1,1e6),[P['h_initial'],vv,rr],method='Radau',
                 t_eval=[1e2,1e3,1e4,1e5,1e6],rtol=2e-11,atol=1e-17)
    check('K22_table_run_completed',ks.success)
    kv=[]
    for j,tt in enumerate(ks.t):
        hv,rv,rad=ks.y[:,j]
        kv.append({'tau':float(tt),'h':float(hv),'xi':float(2+math.sqrt(rv)),
                   'r_V':float(rv),'r_r':float(rad)})
    check('K22_TableIII_headline_reproduced',abs(kv[-1]['h']/.0430-1)<.001 and abs(kv[-1]['r_V']/.000221-1)<.002)
    out['K22_table_III']=kv
    x,lam=sp.symbols('xi lambda',real=True)
    xr=2+sp.sqrt(4-lam)
    resid=sp.simplify((2*x-4).subs(x,xr))
    out['K22_original_action_endpoint_residual']=str(resid)
    check('K22_endpoint_is_not_original_algebraic_field_equation',resid.subs(lam,1)!=0)
    check('K22_H5_source_does_not_stop_at_zero_vacuum_with_radiation',P['g_tilde']*(coeff*.1)**2.5>0)

    # Two vacua in the SAME fixed-mu action, not a separate mu for each phase.
    s,vs,ll,delta=sp.symbols('s v_s lambda delta',positive=True)
    pot=ll*(s*s-vs*vs)**2/4+delta*(2-3*s/vs+(s/vs)**3)/4
    mp=sp.simplify(sp.diff(pot,s,2).subs(s,vs))
    mm=sp.simplify(sp.diff(pot,s,2).subs(s,-vs))
    check('false_vacuum_gap_coexists_with_true_zero',pot.subs(s,vs)==0 and sp.simplify(pot.subs(s,-vs)-delta)==0)
    check('tree_mass_gap_relation',sp.simplify(vs*vs*(mp-mm)/3-delta)==0)
    out['multiple_vacua']={'same_mu_example':{'epsilon':'q^2/2 + 1/2','mu':1,'q':1,'V_true':0,'V_false':0.01,'rho_true':0,'rho_false':0.01},
         'biased_scalar_potential':str(pot),'m_plus_squared':str(mp),'m_minus_squared':str(mm),
         'local_minima_condition':'0 < delta < 4*lambda*v_s^4/3',
         'conclusion':'Final true-vacuum matching alone does not exclude a positive false phase today; its gap and preparation remain inputs.'}
    um=I['EotWash_lambda_digitized_um']
    out['unscreened_limit']={'source':'Lee et al. 2020 + Claude checkpoint SVG digitization',
        'lambda_um':um,'m_meV':.1973269804/um*1000,
        'scope':'alpha=1/3 unscreened, linear laboratory response; not a universal screened bound'}
    om=I['omega_m_reference'];orr=I['omega_r_reference']
    out['fixed_density_flat_zero_residual_background']={'H0_km_s_Mpc':100*math.sqrt(om+orr),
            'q_dec_today':(.5*om+orr)/(om+orr),
            'scope':'Pinned omega_m and omega_r, spatial flatness and heavy q. Not an optimized CMB likelihood.'}
    out['checks']=CHECKS
    out['passed']=all(CHECKS.values())
    (HERE/'results'/'review_results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'checks':len(CHECKS),'all_passed':out['passed'],
                      'N_required':req_n,'age_ratio':rho_age/rho,'upper_envelope_M4':envelope,
                      'actual_final_M4':end_from_one},indent=2))


if __name__=='__main__':
    main()
