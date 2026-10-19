#!/usr/bin/env python3
"""One-loop vacuum X pole threshold in the fixed matched scalar theory.

No finite counterterm is fitted. This is not a finite-charge quantum energy.
Run from any directory. Output files are created beside this source.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import root_scalar
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent/'grut_candidate_1_quantum_20261018'
LOOP = 16*np.pi**2


def load_flow():
    spec=importlib.util.spec_from_file_location('previous_grut_flow',PARENT/'flow_check.py')
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def B0(s,a,b,mu):
    return quad(lambda z:-np.log((z*a+(1-z)*b-z*(1-z)*s)/mu**2),0,1,
                epsabs=1e-12,epsrel=1e-12,limit=100)[0]


def B0_gauss(s,a,b,mu,n):
    x,w=leggauss(n);z=(1+x)/2
    return float(-.5*np.sum(w*np.log((z*a+(1-z)*b-z*(1-z)*s)/mu**2)))


def A(t,mu):
    return t*(np.log(t/mu**2)-1)/LOOP


def result_at(theta,mu):
    lp,lx,c,d,ld,alpha,mp,mx,md,C=theta
    vcl=np.sqrt(-mp/lp)
    a=-2*mp;b=mx+c*vcl**2;D2=md+d*vcl**2;g3=2*c*vcl
    derivatives=np.array([6*lp*vcl,2*c*vcl,2*c*vcl,2*d*vcl])
    masses=np.array([a,b,b,D2])
    tadpole=float(np.sum(masses*derivatives*(np.log(masses/mu**2)-1))/(2*LOOP))
    dv=-tadpole/a
    shift=2*c*vcl*dv
    tad=c*A(a,mu)+4*lx*A(b,mu)+alpha*A(D2,mu)
    bubble0=-g3*g3*B0(0,a,b,mu)/LOOP
    bubblepole=-g3*g3*B0(b,a,b,mu)/LOOP
    curvature=b+shift+tad+bubble0
    pole2=b+shift+tad+bubblepole
    root=root_scalar(lambda s:-s+b+shift+tad-g3*g3*B0(s,a,b,mu)/LOOP,
                     bracket=(.9*b,1.1*b),xtol=1e-15)
    dB=quad(lambda z:z*(1-z)/(z*a+(1-z)*b-z*(1-z)*b),0,1,
            epsabs=1e-12,epsrel=1e-12)[0]
    zmin=np.clip((2*b-a)/(2*b),0,1)
    denommin=zmin*a+(1-zmin)*b-zmin*(1-zmin)*b
    return {'mu_over_v':float(mu),'tree_vacuum_phi_over_v':float(vcl),
            'tree_internal_masses_squared_over_v2':{'phi':float(a),'X':float(b),'D':float(D2)},
            'cubic_phi_xx_over_v':float(g3),'tadpole_over_v3':tadpole,'strict_vacuum_shift_over_v':float(dv),
            'vacuum_displacement_mass_squared_correction_over_v2':float(shift),
            'quartic_tadpole_mass_squared_correction_over_v2':float(tad),
            'bubble_at_zero_over_v2':float(bubble0),'bubble_at_tree_shell_over_v2':float(bubblepole),
            'strict_curvature_squared_over_v2':float(curvature),
            'strict_pole_squared_over_v2':float(pole2),'strict_pole_mass_over_v':float(np.sqrt(b)+.5*(pole2-b)/np.sqrt(b)),
            'sqrt_of_truncated_strict_pole_squared_over_v':float(np.sqrt(pole2)),
            'pole_minus_curvature_squared_over_v2':float(pole2-curvature),
            'B0_zero':float(B0(0,a,b,mu)),'B0_tree_shell':float(B0(b,a,b,mu)),
            'B0_shell_minus_zero':float(B0(b,a,b,mu)-B0(0,a,b,mu)),
            'dB0_ds_at_tree_shell_times_v2':float(dB),
            'strict_pole_residue':float(1-g3*g3*dB/LOOP),
            'scalar_pair_cut_threshold_over_v2':float((np.sqrt(a)+np.sqrt(b))**2),
            'minimum_Feynman_denominator_at_tree_shell_over_v2':float(denommin),
            'implicit_root_of_truncated_one_loop_inverse_over_v2':float(root.root),
            'implicit_minus_strict_pole_squared_over_v2':float(root.root-pole2),
            'implicit_root_converged':bool(root.converged)}


def main():
    start=time.perf_counter()
    protocol_bytes=(HERE/'pole_protocol.json').read_bytes()
    protocol=json.loads(protocol_bytes)
    flow=load_flow();theta0=flow.THETA0.copy()
    assert np.array_equal(theta0,[protocol['parameters'][k] for k in flow.NAMES])
    res=result_at(theta0,.1)
    old=json.loads((PARENT/'potential_results.json').read_text())
    exact=1+np.log(2)-np.pi/2
    gauss=[]
    for n in [32,64,128,256]:
        for s in [0.,.01,res['strict_pole_squared_over_v2']]:
            bg=B0_gauss(s,.02,.01,.1,n);br=B0(s,.02,.01,.1)
            gauss.append({'n':n,'s':s,'B0':bg,'quad_minus_gauss':float(br-bg)})
    fd=[]
    for h in [1e-4,3e-5,1e-5,3e-6]:
        num=(B0(.01+h,.02,.01,.1)-B0(.01-h,.02,.01,.1))/(2*h)
        fd.append({'step':h,'derivative':float(num),'relative_error':float(abs(num/res['dB0_ds_at_tree_shell_times_v2']-1))})
    scales=[]
    for ratio in [.5,1.,2.]:
        if ratio==1:theta=theta0
        else:
            sol=solve_ivp(lambda t,th:flow.beta(th),(0,np.log(ratio)),theta0,method='DOP853',rtol=1e-12,atol=1e-16)
            assert sol.success
            theta=sol.y[:,-1]
        scales.append({'mu_over_mu0':ratio,'parameters':dict(zip(flow.NAMES,map(float,theta))),**result_at(theta,.1*ratio)})
    checks={
        'independent_tadpole_matches_previous_potential':abs(res['tadpole_over_v3']-old['vacuum']['CW_tadpole_at_tree_point'][0])<1e-14,
        'independent_zero_momentum_curvature_matches_previous_potential':abs(res['strict_curvature_squared_over_v2']-old['vacuum']['strict_one_loop_curvature_diagonal'][1])<1e-14,
        'analytic_B0_difference':abs(res['B0_shell_minus_zero']-exact)<1e-13,
        'adaptive_and_Gauss_quadratures_agree':max(abs(v['quad_minus_gauss']) for v in gauss)<1e-11,
        'analytic_momentum_derivative_agrees':min(v['relative_error'] for v in fd)<1e-7,
        'strict_and_implicit_poles_below_scalar_pair_cut':res['minimum_Feynman_denominator_at_tree_shell_over_v2']>0 and 0<res['implicit_root_of_truncated_one_loop_inverse_over_v2']<res['scalar_pair_cut_threshold_over_v2'],
        'implicit_root_converged':res['implicit_root_converged'],
        'pole_shift_below_curvature_shift':res['pole_minus_curvature_squared_over_v2']<0,
        'pole_residue_positive_at_this_order':res['strict_pole_residue']>0,
        'implicit_induced_higher_order_difference_small':abs(res['implicit_minus_strict_pole_squared_over_v2'])<1e-9,
    }
    checks={name:bool(value) for name,value in checks.items()}
    result={'protocol_sha256':hashlib.sha256(protocol_bytes).hexdigest(),
            'source_inputs':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PARENT/'flow_check.py',PARENT/'potential_results.json']},
            'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},
            'matching_result':res,'analytic_B0_difference':float(exact),'quadrature_checks':gauss,'derivative_checks':fd,
            'scale_results':scales,
            'pole_squared_scale_spread_over_v2':float(np.ptp([s['strict_pole_squared_over_v2'] for s in scales])),
            'checks':checks,'all_checks_pass':all(checks.values()),'elapsed_seconds':time.perf_counter()-start,
            'claim_boundary':'One-scalar-loop vacuum charged-particle threshold only; no finite-Q quantum binding/lifetime or gravitational vacuum matching.'}
    (HERE/'pole_results.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'pole_execution.jsonl').write_text(json.dumps({'command':'python pole.py','checks':checks,'elapsed_seconds':result['elapsed_seconds']})+'\n')
    print(json.dumps({'matching_result':res,'checks':checks,'scale_spread_over_v2':result['pole_squared_scale_spread_over_v2'],'elapsed_seconds':result['elapsed_seconds']},indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
