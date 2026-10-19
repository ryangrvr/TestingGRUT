#!/usr/bin/env python3
"""Small physical calculations deciding premises of the proposed research fork.

No new one-loop spectrum is computed. Tree channel energies use tree masses.
The two-point test concerns a conformal scalar and fixed proper-time midpoint.
"""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

from decimal import Decimal as D, getcontext
import numpy as np
import scipy
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
REF=HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz'
G=.1
PI=D('3.141592653589793238462643383279502884197169399375105820974944592307816')


def sinh(x):return (x.exp()-(-x).exp())/2
def cosh(x):return (x.exp()+(-x).exp())/2
def coth(x):return cosh(x)/sinh(x)


def solve_charge(N,L,tol,spline,Qref):
    rr=np.linspace(0,L,int(L*40)+1)
    shrink=(N/Qref)**(1/3)
    y=spline(np.minimum(rr/shrink,50.))
    y[[1,3]]/=shrink
    y[:,rr/shrink>=50]=np.array([0.,0.,1.,0.])[:,None]
    guess=np.vstack([y,cumulative_trapezoid(rr**2*y[0]**2,rr,initial=0)])
    def fun(r,y,p):
        F,Fp,S,Sp,J=y; O=p[0]
        return np.array([Fp,(S*S+.1*F*F-O*O)*F,Sp,(S*S-1+F*F)*S,r*r*F*F])
    def bc(a,b,p):
        return np.array([a[1],a[3],a[4],b[0],b[2]-1,4*np.pi*p[0]*b[4]-N*G**2])
    start=time.perf_counter()
    result=solve_bvp(fun,bc,rr,guess,p=[.8 if N>30000 else .87],S=np.diag([0,-2,0,-2,0]),tol=tol,bc_tol=tol,max_nodes=40000)
    r=np.linspace(0,L,40001); F,Fp,S,Sp,J=result.sol(r); O=float(result.p[0])
    U=.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4
    grad=.5*(Fp*Fp+Sp*Sp); tim=.5*O*O*F*F
    integ=lambda a:4*np.pi*simpson(r*r*a,x=r)
    E=integ(U+grad+tim);Q=O*integ(F*F)
    return {'N':N,'L':L,'tol':tol,'success':bool(result.success),'status':result.status,'message':result.message,'mesh_nodes':len(result.x),'maximum_BVP_residual':float(max(result.rms_residuals)),
      'Omega':O,'Q':float(Q/G**2),'E_over_gv':float(E/G**2),'F0':float(F[0]),'S0':float(S[0]),'minF':float(min(F)),'minS':float(min(S)),
      'virial_relative':float((integ(grad)+3*integ(U-tim))/E),'seconds':time.perf_counter()-start}


def W(e,t1,t2):
    if e==0:
        a=lambda t:t.exp(); eta=lambda t:-(-t).exp()
    else:
        a=lambda t:(1+e*t)**(1/e)
        eta=lambda t:-(1+e*t)**((e-1)/e)/(1-e)
    return -1/(4*PI**2*a(t1)*a(t2)*(eta(t1)-eta(t2))**2)


def thermal(s,h):return -h*h/(16*PI**2*sinh(h*s/2)**2)


def main():
    started=time.perf_counter(); raw=(HERE/'protocol.json').read_bytes(); protocol=json.loads(raw)
    data=np.load(REF); spline=CubicSpline(data['rho'],data['y'],axis=1)
    ref=json.loads((HERE.parent/'grut_candidate_1_20261018/code/results.json').read_text())['observables']; Qref=ref['Qbar']/G**2
    checks={}; runs=[]
    for L,tol in zip(protocol['classical_channels']['domains'],protocol['classical_channels']['tolerances']):
        for N in protocol['classical_channels']['charges']:
            row=solve_charge(N,L,tol,spline,Qref);runs.append(row);print(json.dumps(row),flush=True)
    selected={x['N']:x for x in runs if x['L']==50}
    emission=selected[39944]['E_over_gv']-selected[39943]['E_over_gv']-1
    split=selected[39944]['E_over_gv']-2*selected[19972]['E_over_gv']
    coarse={x['N']:x for x in runs if x['L']==40}
    changes={str(N):abs(selected[N]['E_over_gv']-coarse[N]['E_over_gv']) for N in selected}
    checks['all_six_BVPs_converged']=all(x['success'] and 0<x['Omega']<1 and x['minF']>-1e-7 and x['minS']>0 for x in runs)
    checks['integer_charge_residuals_small']=max(abs(x['Q']-x['N']) for x in runs)<1e-5
    checks['classical_channel_energy_convergence']=max(changes.values())<1e-5
    checks['virial_identity']=max(abs(x['virial_relative']) for x in runs)<1e-7
    checks['tree_single_carrier_channel_closed']=emission<-.19
    checks['tree_equal_split_channel_closed']=split<0
    # The analytic identity is written in DERIVATION.md. Here test the exact
    # power-law expression with arbitrary-precision elementary arithmetic.
    refinement=json.loads((HERE/'refinement_protocol.json').read_text()) if (HERE/'refinement_protocol.json').exists() else None
    epsilon_grid=protocol['slow_roll_test']['epsilon']+(refinement['additional_epsilon'] if refinement else [])
    getcontext().prec=60; slow=[]
    for separation in protocol['slow_roll_test']['separations']:
        s=D(str(separation)); points=[]
        for eps in epsilon_grid:
            e=D(str(eps));mid=W(e,s/2,-s/2); shifted=thermal(s,1-e)
            points.append({'epsilon':eps,'relative_midpoint_residual':float(mid/shifted-1),
              'relative_unshifted_residual':float(mid/thermal(s,1)-1),
              'relative_endpoint_residual':float(W(e,0,-s)/shifted-1),
              'mean_time_translation_change':float(W(e,1+s/2,1-s/2)/mid-1)})
        errors=[abs(x['relative_midpoint_residual']) for x in points]
        ratios=[errors[i]/errors[i+1] for i in range(len(errors)-1)]
        endpoint_limit=s/2*(2-s*coth(s/2))
        calculated=points[-1]['relative_endpoint_residual']/points[-1]['epsilon']
        c2=s*s/4-s**3*coth(s/2)/12
        slow.append({'separation':separation,'points':points,'midpoint_error_halving_ratios':ratios,'endpoint_linear_coefficient_analytic':float(endpoint_limit),'endpoint_linear_coefficient_numerical':calculated,'midpoint_quadratic_coefficient_analytic':float(c2),'midpoint_quadratic_coefficient_numerical':points[-1]['relative_midpoint_residual']/points[-1]['epsilon']**2})
    checks['midpoint_residual_scales_quadratically']=all(abs(r-4)<.1 for row in slow for r in row['midpoint_error_halving_ratios'][-3:]) if refinement else all(3.7<r<4.3 for row in slow for r in row['midpoint_error_halving_ratios'])
    checks['midpoint_second_order_coefficient_matches']=all(abs(row['midpoint_quadratic_coefficient_numerical']/row['midpoint_quadratic_coefficient_analytic']-1)<.02 for row in slow)
    checks['endpoint_correction_is_linear_and_nonzero']=all(abs(row['endpoint_linear_coefficient_numerical']/row['endpoint_linear_coefficient_analytic']-1)<.02 for row in slow)
    checks['proper_time_translation_is_not_stationary_at_first_order']=all(abs(row['points'][-1]['mean_time_translation_change'])>row['points'][-1]['epsilon']*.01 for row in slow)
    # Closed de Sitter analytically continues from a regular Euclidean hemisphere.
    # This checks a sign/magnitude implication, not the observed inflation history.
    curvature=[]
    for N in [1,5,10,60]:
        n=D(N);a=n.exp();H=(1-1/a**2).sqrt()
        k=-1/(a*H)**2; identity=-1/((2*n).exp()-1)
        curvature.append({'efolds':N,'Omega_K':float(k),'identity_relative_residual':float(abs(k/identity-1)),
          'H_FLRW_over_H_Lambda':float(H),'epsilon_H':float(k),
          'imported_geodesic_temperature_over_HLambda_div_2pi':1.,
          'a_double_dot_over_a_dot_over_HLambda':float(1/H),
          'relative_temperature_discrepancy':float(1/H-1),
          'temperature_discrepancy_div_abs_Omega_K':float((1/H-1)/abs(k))})
    checks['closed_de_Sitter_curvature_identity']=max(row['identity_relative_residual'] for row in curvature)<1e-50
    checks['positive_curvature_can_be_below_1e_minus3']=curvature[1]['Omega_K']<0 and abs(curvature[1]['Omega_K'])<1e-3
    checks['closed_geometry_changes_flat_temperature_formula_at_curvature_order']=abs(curvature[1]['temperature_discrepancy_div_abs_Omega_K']-.5)<1e-4
    # The compact flat torus has constant spatial metric and periodic coordinates.
    result={'protocol_sha256':hashlib.sha256(raw).hexdigest(),'refinement_protocol_sha256':hashlib.sha256((HERE/'refinement_protocol.json').read_bytes()).hexdigest() if refinement else None,'curvature_transfer_protocol_sha256':hashlib.sha256((HERE/'curvature_transfer_protocol.json').read_bytes()).hexdigest(),'profile_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),
      'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'elementary_arithmetic':'Python Decimal, 60 digits','platform':platform.platform()},
      'classical_BVP_runs':runs,'classical_energy_changes_over_gv':changes,
      'classical_channels':{'N':39944,'tree_X_mass_over_gv':1.,'Delta_single_emission_over_gv':emission,'Delta_equal_split_over_gv':split,'Omega_at_N':selected[39944]['Omega'],'E_N_minus_E_Nminus1_over_gv':emission+1,'E_N_over_gv':selected[39944]['E_over_gv'],'E_half_over_gv':selected[19972]['E_over_gv'],'order':'Tree level, same original couplings. No quantum fission or lifetime claim.'},
      'analytical_midpoint_identity':'See DERIVATION.md; exact power-law predictions checked below, no computer-algebra claim.','slow_roll_points':slow,'curvature_counterexample':curvature,
      'flat_torus':{'metric':'dell^2=a(t)^2 (dx^2+dy^2+dz^2), x,y,z periodically identified with periods Lx,Ly,Lz','spatial_volume':'a(t)^3 Lx Ly Lz','spatial_Ricci_scalar':0,'purpose':'Topological compactness alone does not imply k=+1; this is not a solution of every no-boundary or sequestering boundary-value problem.'},
      'checks':{key:bool(value) for key,value in checks.items()},'all_checks_pass':all(checks.values()),'elapsed_seconds':time.perf_counter()-started,
      'scope':'Classical channel arithmetic, a conformal-field finite-separation expansion, and exact geometric counterexamples. No new quantum energy, no general no-boundary theorem, no empirical fit, no novelty certification.'}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
