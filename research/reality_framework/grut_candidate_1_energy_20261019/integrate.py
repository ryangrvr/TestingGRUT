#!/usr/bin/env python3
"""Integrate the audited energy, integer charge and like-order pole threshold."""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_bvp, simpson, cumulative_trapezoid
import spectral as sp

HERE=Path(__file__).resolve().parent


def load(name):return json.loads((HERE/name).read_text())


def integer_check(N, g):
    r=np.linspace(0,50,2001); y=sp.SPL(r); guess=np.vstack([y,cumulative_trapezoid(r*r*y[0]**2,r,initial=0)])
    def fun(r,y,p):
        F,Fp,S,Sp,J=y; O=p[0]
        return np.array([Fp,(S*S+.1*F*F-O*O)*F,Sp,(S*S-1+F*F)*S,r*r*F*F])
    def bc(a,b,p):return np.array([a[1],a[3],a[4],b[0],b[2]-1,4*np.pi*p[0]*b[4]-N*g*g])
    sol=solve_bvp(fun,bc,r,guess,p=[sp.OMEGA],S=np.diag([0,-2,0,-2,0]),tol=1e-10,bc_tol=1e-10,max_nodes=30000)
    rr=np.linspace(0,50,40001); F,Fp,S,Sp,J=sol.sol(rr); O=float(sol.p[0])
    E=4*np.pi*simpson(rr**2*(.5*(Fp*Fp+Sp*Sp)+.5*O*O*F*F+.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4),x=rr)
    Q=4*np.pi*O*simpson(rr**2*F*F,x=rr)
    def sens_fun(r,y):
        u,up,z,zp=y; F,_,S,_=sp.SPL(r)
        return np.array([up,(S*S+.3*F*F-sp.OMEGA**2)*u+2*F*S*z-2*sp.OMEGA*F,zp,(3*S*S-1+F*F)*z+2*F*S*u])
    sens=solve_bvp(sens_fun,lambda a,b:np.array([a[1],a[3],b[0],b[2]]),r,np.zeros((4,len(r))),S=np.diag([0,-2,0,-2]),tol=1e-10,bc_tol=1e-10,max_nodes=30000)
    F=sp.SPL(rr)[0]; derivative=4*np.pi*simpson(rr**2*(F*F+2*sp.OMEGA*F*sens.sol(rr)[0]),x=rr)
    return {'integer_BVP_success':bool(sol.success),'integer_BVP_max_residual':float(max(sol.rms_residuals)), 'Omega_at_N':O,'Ebar_at_N':float(E),'Qbar_at_N':float(Q),'sensitivity_success':bool(sens.success),'dQbar_dOmega_reference':float(derivative),'sensitivity_max_residual':float(max(sens.rms_residuals))}


def main():
    started=time.perf_counter(); graph=load('graphs_results.json'); spectral=load('spectral_results.json'); audit=load('audit_results.json'); tail=load('tail_refinement_results.json'); toy=load('theory_toy_results.json')
    if not audit.get('completed') or not audit['all_checks_pass'] or not graph['all_checks_pass'] or not tail['all_checks_pass']:raise RuntimeError('Incomplete or failed audit')
    original=json.loads((HERE.parent/'grut_candidate_1_20261018/code/results.json').read_text())['observables']
    pole_path=HERE.parent/'grut_candidate_1_fluctuations_20261019/pole_results.json'; pole=json.loads(pole_path.read_text())
    g=.1; Qref=original['Qbar']/g**2; N=round(Qref); delta=N-Qref
    integer=integer_check(N,g)
    Eref=original['Ebar']/g**2; Eint=Eref+sp.OMEGA*delta
    graphs_at_R25=next(x for x in audit['support_graphs'] if x['R']==25.)
    restored=graphs_at_R25['first']+graphs_at_R25['second']; K=restored+audit['energy_replay']
    uv=audit['tails'][1]['UV_bound']; angular=tail['angular_bound']; empirical=.02; graph_allowance=1e-6; classical_allowance=1e-5
    total_allowance=uv+angular+empirical+graph_allowance+classical_allowance
    # Round UP, never toward a favorable binding decision; not a rigorous exact-profile interval.
    quoted_allowance=np.ceil(total_allowance*100)/100
    mass=pole['matching_result']['strict_pole_mass_over_v']; threshold=N*mass/g
    energy=Eint+K; binding=energy-threshold
    quadratic=g*g*delta*delta/(2*integer['dQbar_dOmega_reference'])
    checks={'N_is_frozen_nearest_integer':N==39944,'all_inputs_same_profile':graph['profile_sha256']==audit['profile_sha256']==hashlib.sha256(sp.REF.read_bytes()).hexdigest(),
       'finite_cutoff_complete':spectral['complete_n6'] and spectral['complete_n8'] and len(spectral['skipped'])==0,
       'all_spectral_points_succeeded':all(x['success'] for x in spectral['points']),
       'no_discarded_large_determinant_phase':max(abs(x['logJost_phase_mod_2pi']) for x in spectral['points']+audit['replay_points'])<1e-7,
       'toy_exact_charge_identity':toy['exact_full_minus_zeta2_reduced_residual']=={},
       'toy_frequency_integral_agrees':abs(toy['integral_frequency_difference'])<1e-10,
       'pole_calculation_passed':pole['all_checks_pass'],
       'integer_and_sensitivity_BVP_converged':integer['integer_BVP_success'] and integer['sensitivity_success'],
       'integer_energy_expansion_matches_direct_BVP':abs(integer['Ebar_at_N']/g**2-Eint)<1e-5,
       'numerical_freeX_binding_interval_negative':binding+quoted_allowance<0}
    result={'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'parent_commit':'1f370987f2fa96a52426ec95680a4f3ad7cbd373','g':g,'Omega_reference':sp.OMEGA,'Q_reference_continuous':Qref,'N_integer':N,'delta_charge':delta,'Eclass_reference_over_gv':Eref,'integer_linear_adjustment_over_gv':sp.OMEGA*delta,'Eclass_integer_first_order_over_gv':Eint,'integer_check':integer,'integer_quadratic_estimate_over_gv':quadratic,'direct_integer_minus_linear_over_gv':integer['Ebar_at_N']/g**2-Eint,
      'spectral_remainder_initial_n6':spectral['energy_n6'],'spectral_remainder_initial_n8':spectral['energy_n8'],'spectral_remainder_replay':audit['energy_replay'],'spectral_replay_change':audit['replay_difference'],'graphs_matched_R25':graphs_at_R25,'renormalized_K1_over_gv':K,'renormalized_K1_over_gv_numerical_range':[K-quoted_allowance,K+quoted_allowance],
      'error_budget_over_gv':{'frequency_tail_analytical':uv,'angular_tail_analytical':angular,'empirical_numerical_allowance':empirical,'graph_allowance':graph_allowance,'classical_allowance':classical_allowance,'sum_before_roundup':total_allowance,'rounded_reported_allowance':float(quoted_allowance),'qualification':'Analytical truncation bounds conditional on supplied spline enclosures, plus an empirical numerical allowance. Not a certified exact-profile enclosure, statistical confidence interval, or higher-loop error bound.'},
      'X_strict_one_loop_pole_mass_over_v':mass,'free_NX_threshold_over_gv':threshold,'soliton_0plus1_energy_over_gv':energy,'DeltaE_freeX_over_gv':binding,'DeltaE_freeX_numerical_range':[binding-quoted_allowance,binding+quoted_allowance],'loop_shift_percent_of_classical_energy':100*K/Eint,'E_over_free_threshold':energy/threshold,
      'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':all(checks.values()),'claim':'Numerical one-scalar-loop support for binding against N separated free X carriers at the fixed original benchmark, conditional on the stable-saddle semiclassical construction. Not global quantum stability, lifetime, full-continuum stability proof, empirical viability or a universal theory.',
      'execution_seconds':time.perf_counter()-started}
    result['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [sp.REF,pole_path,HERE/'graphs_results.json',HERE/'spectral_results.json',HERE/'audit_results.json',HERE/'tail_refinement_results.json',HERE/'binding_protocol.json',HERE/'integrate.py']}
    (HERE/'binding_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
