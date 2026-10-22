#!/usr/bin/env python3
"""Conserved late-pole homogeneous-position probe; no old solver imports."""
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
I = json.loads((ROOT / 'INPUTS.json').read_text())
M, b = I['M_GeV'], I['b']
kappa = I['kappa_over_c_squared'] * I['c']**2
if I['c'] != 1:
    raise ValueError('This packet specifies c=1; do not silently change its exact-parent comparisons.')
V0 = I['lambda_H'] * I['v_GeV']**4 / 4
xend = math.log1p(math.sqrt(2/kappa))
xpivot = brentq(lambda x: kappa/2*(math.exp(x)-x-math.exp(xend)+xend)-I['N_pivot'], xend, 20)
epsilon = 2/(kappa*math.expm1(xpivot)**2)
alpha = math.sqrt(24*math.pi**2*I['A_s']*epsilon/(-math.expm1(-xpivot))**2)
A = alpha**2*M**4
s0 = V0/A
smin_exact = math.log1p(V0/A)
rho_star = V0**2/A
mI = math.sqrt(2/kappa)*alpha*M
rad_coeff = math.pi**2*I['g_radiation']/30
rho_stop = rad_coeff*I['T_stop_GeV']**4
Hstop2 = (rho_stop+V0)/(3*M*M)
minimum_mass2 = 2*rho_star/(b*b*M*M)
qstop = minimum_mass2/Hstop2
CHECKS = []


def check(name, passed, value=None, limit=None):
    CHECKS.append(dict(name=name, passed=bool(passed), value=value, limit=limit))


def potential(u):
    """Leading controlled small-positive-sigma potential; force kept separately."""
    eu = math.exp(u)
    post = rho_star*eu*eu
    total = V0 + post - 2*rho_star*eu
    force = 2*rho_star*eu*math.expm1(u)
    return total, force, post


def algebraic(u, vel, N, rhoi):
    U, force, post = potential(u)
    rad = rhoi*math.exp(-4*N)
    H2 = (rad+U)/(M*M*(3-.5*b*b*vel*vel))
    eps = .5*b*b*vel*vel + 2*rad/(3*M*M*H2)
    return U, force, post, rad, H2, eps


def trajectory(T, u0, vel0=0., refined=False):
    rhoi = rad_coeff*T**4
    Nend = math.log(T/I['T_stop_GeV'])
    Hini = math.sqrt(algebraic(u0, vel0, 0, rhoi)[4])

    def rhs(N, y):
        u, vel, _ = y
        _, force, _, _, H2, eps = algebraic(u, vel, N, rhoi)
        return vel, -(3-eps)*vel-force/(b*b*M*M*H2), -eps

    sol = solve_ivp(rhs, (0, Nend), [u0, vel0, math.log(Hini)], method='DOP853',
                    rtol=2e-12 if refined else 2e-10,
                    atol=[2e-12, 2e-13, 2e-13] if refined else [2e-10, 2e-11, 2e-11],
                    max_step=.025 if refined else .1, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    grid = np.linspace(0, Nend, I['raw_points'])
    state = sol.sol(grid)
    raw = []
    for n, (u, vel, lh) in zip(grid, state.T):
        U, force, post, rad, H2, eps = algebraic(u, vel, n, rhoi)
        kin = .5*b*b*M*M*H2*vel*vel
        total = 3*M*M*H2
        raw.append([n, u, vel, math.sqrt(H2), rad/total, kin/total, U/total,
                    eps, lh-.5*math.log(H2), post/rho_stop])
    raw = np.array(raw)
    U, force, post, rad, H2, eps = algebraic(state[0,-1], state[1,-1], Nend, rhoi)
    initial = algebraic(u0, vel0, 0, rhoi)
    row = dict(T_initial_GeV=T, initial_u=u0, initial_du_dN=vel0, N_available=Nend,
               final_u=float(state[0,-1]), final_du_dN=float(state[1,-1]),
               final_log10_position_ratio=float(state[0,-1]/math.log(10)),
               final_postEW_potential_GeV4=post,
               postEW_energy_ratio_to_prepared=float(math.exp(2*state[0,-1])),
               final_radiation_fraction=float(raw[-1,4]),
               final_scalar_kinetic_fraction=float(raw[-1,5]),
               final_epsilon=float(eps),
               initial_radiation_fraction=float(initial[3]/(3*M*M*initial[4])),
               max_abs_integrated_logH_error=float(np.max(np.abs(raw[:,8]))),
               max_abs_du_dN=float(np.max(np.abs(raw[:,2]))),
               distance_lower_bound_N=b*abs(u0)/math.sqrt(6),
               evaluation_count=sol.nfev)
    return row, raw, sol


def cosmic_time_control(T, u0, vel0, Nend):
    """Second-order cosmic-time KG, log of cosmic time as numerical clock."""
    rhoi = rad_coeff*T**4
    Href = math.sqrt(rhoi/(3*M*M))
    Hini = math.sqrt(algebraic(u0, vel0, 0, rhoi)[4])

    def rhs(z, y):
        n, u, W = y  # W=du/d(Href*t), tau=(exp(z)-1)/2
        U, force, _ = potential(u)
        h = math.sqrt(math.exp(-4*n)+U/rhoi+b*b*W*W/6)
        clock = .5*math.exp(z)
        return clock*h, clock*W, clock*(-3*h*W-3*force/(b*b*rhoi))

    def stop(z, y):
        return y[0]-Nend
    stop.terminal = True
    stop.direction = 1
    sol = solve_ivp(rhs, (0, 2*Nend+5), [0, u0, vel0*Hini/Href], events=stop,
                    method='DOP853', rtol=3e-12, atol=[2e-12, 2e-12, 1e-100],
                    max_step=.04)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError('Cosmic-time endpoint not reached: '+sol.message)
    n, u, W = sol.y_events[0][0]
    U, _, _ = potential(u)
    h = math.sqrt(math.exp(-4*n)+U/rhoi+b*b*W*W/6)
    return dict(final_u=float(u), final_du_dN=float(W/h), final_H_GeV=float(Href*h),
                clock_endpoint=float(sol.t_events[0][0]), evaluation_count=sol.nfev)


def minimum_response(T):
    """u/u_initial=1+qstop*R+O(qstop^2); resolves a ~1e-59 change."""
    rhoi = rad_coeff*T**4
    Nend = math.log(T/I['T_stop_GeV'])

    def rhs(N, y):
        rad = rhoi*math.exp(-4*N)
        eps = 2*rad/(rad+V0)
        Hratio = (rho_stop+V0)/(rad+V0)
        return y[1], -(3-eps)*y[1]-Hratio

    sol = solve_ivp(rhs, (0, Nend), [0., 0.], method='DOP853', rtol=2e-12,
                    atol=2e-14, max_step=.04)
    R, Rp = sol.y[:,-1]
    # Separate integrating-factor quadrature for the first-order Green function.
    # R(Nf)=-int dN a^3/H * m^2 * int_N^Nf dL/(a^3 H), divided by qstop.
    from scipy.integrate import quad

    def H_ratio_to_end(n):
        return math.sqrt((rhoi*math.exp(-4*n)+V0)/(rho_stop+V0))

    def outer(n):
        inner = quad(lambda ell: math.exp(-3*(ell-n))/H_ratio_to_end(ell), n, Nend,
                     epsabs=1e-13, epsrel=1e-11)[0]
        return inner/H_ratio_to_end(n)

    Rquad = -quad(outer, 0, Nend, epsabs=1e-13, epsrel=1e-11)[0]
    return dict(T_initial_GeV=T, scaled_displacement_response=float(R),
                scaled_velocity_response=float(Rp), quadrature_response=float(Rquad),
                leading_fractional_contraction=float(-qstop*R),
                radiation_only_asymptote=.05)


def phase_rhs(lam):
    def rhs(n, xy):
        X, Y = xy
        eps = 3*X*X+2*(1-X*X-Y*Y)
        return [-3*X+math.sqrt(1.5)*lam*Y*Y+eps*X,
                -math.sqrt(1.5)*lam*X*Y+eps*Y]
    return rhs


def phase_plane(bb):
    lam = 2/bb
    if bb < 1:
        X = math.sqrt(1.5)*(4/3)/lam
        Y = math.sqrt(4/(3*lam*lam))
        kind = 'radiation_scaling'
    else:
        X = lam/math.sqrt(6)
        Y = math.sqrt(1-lam*lam/6)
        kind = 'scalar_dominated'
    f = phase_rhs(lam)
    h = 1e-5
    J = np.column_stack([(np.array(f(0,[X+h,Y]))-np.array(f(0,[X-h,Y])))/(2*h),
                         (np.array(f(0,[X,Y+h]))-np.array(f(0,[X,Y-h])))/(2*h)])
    eigen = np.linalg.eigvals(J)
    rows = []
    for start in I['phase_plane_initial_X_Y']:
        sol = solve_ivp(f, (0, I['phase_plane_N_end']), start, method='DOP853',
                        rtol=2e-12, atol=2e-13, max_step=.1, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        grid = np.linspace(0, I['phase_plane_N_end'], 801)
        xy = sol.sol(grid)
        distance = float(np.max(np.abs(xy[:,-1]-[X,Y])))
        rows.append(dict(start=start, endpoint=xy[:,-1].tolist(), endpoint_error=distance))
    return dict(b=bb, lambda_exponential=lam, attractor=kind, fixed_point=[X,Y],
                scalar_fraction=X*X+Y*Y, total_w=-1+2/3*(3*X*X+2*(1-X*X-Y*Y)),
                jacobian_eigenvalues=[[float(e.real),float(e.imag)] for e in eigen],
                fixed_point_residual=float(np.max(np.abs(f(0,[X,Y])))), trajectories=rows)


def main():
    data = ROOT/'data'
    data.mkdir(exist_ok=True)
    cases = []
    raw_named = {}
    for T in I['T_initial_GeV']:
        for u in I['initial_log_position_over_minimum']:
            row, raw, _ = trajectory(T, u)
            cases.append(row)
            if T == I['T_initial_GeV'][0] and u in [-10.,0.,60.,90.,95.]:
                name = 'history_u'+str(int(u)).replace('-','minus')+'.csv'
                raw_named[name] = raw
    for u in I['velocity_controls_log_position']:
        for vel in I['velocity_controls_dlog_position_dN']:
            row, _, _ = trajectory(I['T_initial_GeV'][0], u, vel)
            cases.append(row)

    columns = ['N','u','du_dN','H_GeV','Omega_radiation','Omega_kinetic','Omega_potential',
               'epsilon','integrated_logH_error','postEW_potential_over_radiation_at_stop']
    for name, raw in raw_named.items():
        np.savetxt(data/name, raw, delimiter=',', header=','.join(columns), comments='', fmt='%.17e')
    simple_keys = [k for k in cases[0] if k != 'evaluation_count']
    with (data/'endpoint_map.csv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=simple_keys,lineterminator='\n')
        writer.writeheader()
        writer.writerows({k:v for k,v in row.items() if k in simple_keys} for row in cases)

    controls=[]
    for u,vel in [(60.,0.), (90.,0.), (95.,0.), (0.,.1)]:
        T=I['T_initial_GeV'][0]
        base, raw0, _=trajectory(T,u,vel)
        fine, raw1, _=trajectory(T,u,vel,refined=True)
        independent=cosmic_time_control(T,u,vel,base['N_available'])
        control=dict(initial_u=u,initial_du_dN=vel,
                     refinement_max_abs_u=float(np.max(np.abs(raw0[:,1]-raw1[:,1]))),
                     refinement_max_abs_velocity=float(np.max(np.abs(raw0[:,2]-raw1[:,2]))),
                     refinement_max_relative_H=float(np.max(np.abs(raw0[:,3]/raw1[:,3]-1))),
                     cosmic_abs_u_error=abs(independent['final_u']-fine['final_u']),
                     cosmic_abs_velocity_error=abs(independent['final_du_dN']-fine['final_du_dN']),
                     cosmic_relative_H_error=abs(independent['final_H_GeV']/raw1[-1,3]-1),
                     independent=independent)
        controls.append(control)

    responses=[minimum_response(T) for T in I['T_initial_GeV']]
    phases=[phase_plane(bb) for bb in I['phase_plane_b_values']]
    maxs=max(s0*math.exp(row['initial_u']+abs(row['initial_du_dN'])*row['N_available']) for row in cases)
    maxH=max(math.sqrt(algebraic(row['initial_u'],row['initial_du_dN'],0,
                                    rad_coeff*row['T_initial_GeV']**4)[4]) for row in cases)
    # Exact c=1 metric: K(s)=kappa*s^2+b^2[s/(2sinh(s/2))]^2.
    # |K/b^2-1| <= (kappa/b^2+1/12)*smax^2 to leading order; add factor two margin.
    metric_bound=2*(kappa/b**2+1/12)*maxs**2
    potential_bound=4*maxs
    heavy_bound=6/kappa*(maxH/mI)**2
    theory=dict(alpha=alpha,kappa=kappa,N_pivot=I['N_pivot'],x_pivot=xpivot,
                V0_GeV4=V0,A_GeV4=A,s0=s0,exact_parent_minimum=smin_exact,
                rho_star_GeV4=rho_star,m_inflaton_GeV=mI,
                rho_radiation_stop_GeV4=rho_stop,H_stop_GeV=math.sqrt(Hstop2),
                canonical_mass_squared_at_minimum_GeV2=minimum_mass2,
                mass_squared_over_Hstop_squared=qstop,
                max_sigma_over_M_conservative=maxs,
                kinetic_metric_relative_error_bound=metric_bound,
                small_sigma_component_error_bound=potential_bound,
                heavy_inflaton_adiabatic_error_scale=heavy_bound,
                radiation_scaling_b_requires_less_than=1.,
                scalar_acceleration_b_requires_greater_than=math.sqrt(2),
                max_trace_acceleration_if_abs_trace_le_4rho=6*maxs/b**2)
    check('all_initial_conditions_radiation_dominated', min(r['initial_radiation_fraction'] for r in cases)>.98,
          min(r['initial_radiation_fraction'] for r in cases),.98)
    check('energy_conservation_integrated_logH',max(r['max_abs_integrated_logH_error'] for r in cases)<2e-8,
          max(r['max_abs_integrated_logH_error'] for r in cases),2e-8)
    check('positive_Friedmann_kinetic_budget',max(r['max_abs_du_dN'] for r in cases)<math.sqrt(6)/b,
          max(r['max_abs_du_dN'] for r in cases),math.sqrt(6)/b)
    check('controlled_small_positive_sigma',maxs<1e-12,maxs,1e-12)
    check('heavy_inflaton_control',heavy_bound<1e-12,heavy_bound,1e-12)
    for j,r in enumerate(controls):
        check(f'refinement_u_{j}',r['refinement_max_abs_u']<2e-8,r['refinement_max_abs_u'],2e-8)
        check(f'refinement_H_{j}',r['refinement_max_relative_H']<2e-8,r['refinement_max_relative_H'],2e-8)
        check(f'cosmic_time_u_{j}',r['cosmic_abs_u_error']<2e-8,r['cosmic_abs_u_error'],2e-8)
        check(f'cosmic_time_H_{j}',r['cosmic_relative_H_error']<2e-8,r['cosmic_relative_H_error'],2e-8)
    for j,r in enumerate(responses):
        check(f'minimum_response_quadrature_{j}',abs(r['scaled_displacement_response']-r['quadrature_response'])<2e-10,
              abs(r['scaled_displacement_response']-r['quadrature_response']),2e-10)
        check(f'minimum_has_negligible_contraction_{j}',0<r['leading_fractional_contraction']<1e-55,
              r['leading_fractional_contraction'],1e-55)
    for j,p in enumerate(phases):
        err=max(r['endpoint_error'] for r in p['trajectories'])
        check(f'phase_fixed_point_{j}',p['fixed_point_residual']<1e-12,p['fixed_point_residual'],1e-12)
        check(f'phase_attractor_reached_{j}',err<2e-7,err,2e-7)
        check(f'phase_attractor_stable_{j}',max(e[0] for e in p['jacobian_eigenvalues'])<0,
              max(e[0] for e in p['jacobian_eigenvalues']),0)
    for T in I['T_initial_GeV']:
        low=next(r for r in cases if r['T_initial_GeV']==T and r['initial_u']==-10 and r['initial_du_dN']==0)
        high=next(r for r in cases if r['T_initial_GeV']==T and r['initial_u']==90 and r['initial_du_dN']==0)
        check(f'below_target_retains_preparation_{T:g}',abs(low['final_u']+10)<1e-10,low['final_u'],-10)
        check(f'high_start_not_selected_{T:g}',high['final_u']>60,high['final_u'],60)
    result=dict(classification='No generic post-inflaton position selection in this controlled late-pole interval; no all-actions no-go.',
                theory=theory,histories=cases,controls=controls,minimum_responses=responses,phase_plane=phases,
                checks=CHECKS,check_count=len(CHECKS),failed_checks=sum(not c['passed'] for c in CHECKS))
    (ROOT/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    env=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
             input_sha256=hashlib.sha256((ROOT/'INPUTS.json').read_bytes()).hexdigest())
    (ROOT/'ENVIRONMENT.json').write_text(json.dumps(env,indent=2)+'\n')
    print(json.dumps(dict(checks=len(CHECKS),failed=result['failed_checks'],theory=theory,
                          controls=controls,minimum_responses=responses),indent=2))
    if result['failed_checks']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
