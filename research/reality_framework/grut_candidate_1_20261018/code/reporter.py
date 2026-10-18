#!/usr/bin/env python3
"""Weak reporter propagation on the numerically computed Candidate 1 background.

This is the linear test-field limit of the SAME action, not a full quantum
scattering calculation with soliton recoil, inelasticity or graviton loops.
Run after qball.py: python reporter.py
"""
from pathlib import Path
import json, hashlib, sys, time
import numpy as np
import scipy
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal

ROOT = Path(__file__).resolve().parent
CONFIG = {
    "a_D": .01, "m_D0_over_gv": 1.,
    "momenta_over_gv": [3., 5., 7.],
    "radial_domain": 50.,
    "solver_tolerances": [1e-8, 1e-10, 1e-12],
    "packet": {"k": 5., "center": 22., "sigma": 3.,
               "domain": 80., "grid_points": 1400,
               "refinement_grids": [700, 1400, 2800],
               "times": [0., 10., 20., 30., 40.]},
    "scope": "classical linear test field and its ideal elastic-mode interferometer",
}


def phase_difference(x):
    return float((x + np.pi/2) % np.pi - np.pi/2)


def main():
    t0 = time.perf_counter()
    profile_path = ROOT / 'results_profile.npz'
    if not profile_path.exists():
        raise SystemExit('Run qball.py successfully before reporter.py.')
    report = {'configuration': CONFIG,
              'profile_sha256': hashlib.sha256(profile_path.read_bytes()).hexdigest(),
              'versions': {'python': sys.version, 'numpy': np.__version__,
                           'scipy': scipy.__version__},
              'disclaimer': 'No full coupled quantum or empirical claim.'}
    (ROOT/'reporter_prescription.json').write_text(json.dumps(report, indent=2)+'\n')
    d = np.load(profile_path)
    r = d['rho']; S = d['y'][2]
    spline = CubicSpline(r, S)
    def sval(x):
        a = np.asarray(x)
        return np.where(a <= r[-1], spline(np.minimum(a, r[-1])), 1.)
    def V(x):
        return CONFIG['a_D']*(sval(x)**2-1.)
    R = CONFIG['radial_domain']
    rq = np.linspace(0, R, 40001)
    abs_integral = float(simpson(np.abs(V(rq)), x=rq))
    report['potential_integral_abs'] = abs_integral
    results = []
    for k in CONFIG['momenta_over_gv']:
        exact=[]
        for tol in CONFIG['solver_tolerances']:
            # Variable-phase equation: exact for this real s-wave test potential.
            phase = solve_ivp(lambda x,y: [-V(x)/k*np.sin(k*x+y[0])**2],
                              [0,R],[0.],method='DOP853',rtol=tol,atol=tol*.01,
                              max_step=.15)
            # Independent radial wave equation, fitted to A sin(kr+delta).
            wave = solve_ivp(lambda x,y:[y[1],(V(x)-k*k)*y[0]],
                             [0,R],[0.,1.],method='DOP853',rtol=tol,atol=tol*.01,
                             max_step=.15)
            assert phase.success and wave.success
            delta=float(phase.y[0,-1])
            delta_wave=phase_difference(np.arctan2(k*wave.y[0,-1],wave.y[1,-1])-k*R)
            exact.append({'tolerance':tol,'delta':delta,'delta_wave_fit':delta_wave,
                          'method_difference':abs(delta-delta_wave)})
        delta=exact[-1]['delta']
        born=float(-simpson(V(rq)*np.sin(k*rq)**2,x=rq)/k)
        remainder=abs_integral**2/(2*k*k)
        results.append({'k_over_gv':k,
                        'energy_over_gv':float(np.sqrt(k*k+1.+CONFIG['a_D'])),
                        'phase_delta_s_wave':delta,
                        'S_phase':2*delta,
                        'born_delta':born,
                        'actual_born_difference':abs(delta-born),
                        'analytic_born_remainder_bound_for_exact_profile':remainder,
                        'P_plus_reference_phase_pi_over_2':float((1-np.sin(2*delta))/2),
                        'convergence':exact})
        assert abs(delta-born) < remainder
        assert exact[-1]['method_difference'] < 5e-8
    report['phase_results']=results

    # An executed positive-frequency radial packet on the same fixed background.
    # sqrt(K) evolves the positive-frequency solution of the local KG equation.
    p=CONFIG['packet']
    packet_runs=[]
    for n in p['refinement_grids']:
        packet_runs.append(packet_run(sval, n))
    report['packet_grid_runs']=packet_runs
    trajectories=next(x['trajectory'] for x in packet_runs if x['N']==p['grid_points'])
    report['packet_trajectory']=trajectories
    report['maximum_packet_norm_error']=max(abs(x['norm']-1) for run in packet_runs for x in run['trajectory'])
    report['seconds']=time.perf_counter()-t0
    report['status']='test_field_calculation_complete'
    (ROOT/'reporter_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'phase_results':results,
                      'packet_refinements':packet_runs,'seconds':report['seconds']}))


def packet_run(sval, n):
    p=CONFIG['packet']; L=p['domain']; h=L/(n+1)
    rr=h*np.arange(1,n+1)
    off=-np.ones(n-1)/h**2
    diag=2/h**2 + 1. + CONFIG['a_D']*sval(rr)**2
    diag0=np.full(n,2/h**2+1.+CONFIG['a_D'])
    e,U=eigh_tridiagonal(diag,off)
    e0,U0=eigh_tridiagonal(diag0,off)
    assert np.min(e)>0 and np.min(e0)>0
    seed=np.exp(-.5*((rr-p['center'])/p['sigma'])**2)*np.exp(-1j*p['k']*rr)
    seed/=np.linalg.norm(seed)
    c=U.T@seed;c0=U0.T@seed
    trajectories=[]
    for t in p['times']:
        u=U@(np.exp(-1j*np.sqrt(e)*t)*c)
        u0=U0@(np.exp(-1j*np.sqrt(e0)*t)*c0)
        mask=(rr>=12)&(rr<=60)
        norm=float(np.vdot(u,u).real); norm0=float(np.vdot(u0,u0).real)
        z=np.vdot(u0[mask],u[mask])
        ps=float(np.sum(abs(u[mask])**2))
        pr=float(np.sum(abs(u0[mask])**2))
        selected=(ps+pr)/2
        # Arm state (|signal>u + |reference>u0)/sqrt(2), beta=pi/2.
        # + record projects (|signal> + exp(+i beta)|reference>)/sqrt(2).
        plus=(ps+pr-2*z.imag)/4
        trajectories.append({'time_in_inverse_gv':t,'norm':norm,'reference_norm':norm0,
                             'mean_radius':float(np.sum(rr*abs(u)**2)),
                             'outside_readout_probability':ps,
                             'reference_readout_probability':pr,
                             'overlap_with_reference':{'real':float(z.real),'imag':float(z.imag)},
                             'relative_phase_on_readout':float(np.angle(z)),
                             'selected_probability':selected,
                             'joint_plus_probability':float(plus),
                             'joint_minus_probability':float(selected-plus),
                             'conditional_plus_probability':float(plus/selected),
                             'visibility':float(2*abs(z)/(ps+pr))})
    return {'N':n,'h':h,'trajectory':trajectories}


if __name__=='__main__':
    main()
