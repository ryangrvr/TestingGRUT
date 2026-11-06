#!/usr/bin/env python3
"""Independent unrotated ODE and vacuum controls for the coupled solver."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

from run import compile_engine, run_engine


def model(**changes):
    base = dict(
        expand=True, subtract=True, dressed_start=False, nq=4, nsteps=160, stride=1,
        t0=0.0, t1=1.6, qmax=4.0, mp=4.0, **{'lambda': 0.08}, v=1.3, g2=0.22,
        sc=3.0, c1=0.12, m2=-0.25, mu2=1.0, cb=0.15, offset=2.4, gs2=0.37,
        gamma_s=0.13, a0=0.91, H0=0.23, rho_sm=0.06, rho_dr=0.03,
        fields=[0.43, 0.70, 0.28], velocities=[0.89, -0.30, 0.21],
        species=[(0, 0.85, 1.30, 0.2), (1, 0.60, 0.90, 0.1)],
    )
    return dict(base, **changes)


def inertias(cfg, a):
    """Independent quadrature of the two adiabatic energy coefficients."""
    coef = np.ones(3)
    slope = np.zeros(3)
    gravity = 0.0
    if cfg['subtract']:
        cutoff = cfg['qmax'] / a
        for parent, mass, f, _ in cfg['species']:
            if mass == 0:
                continue
            factor = mass**2 / (8 * math.pi**2)
            coeff = quad(lambda p: p*p/math.hypot(p, mass)**3, 0, cutoff, epsabs=2e-12)[0]
            gravity += factor * quad(lambda p: p**4/math.hypot(p, mass)**5, 0, cutoff, epsabs=2e-12)[0]
            coef[parent] -= factor * coeff / f**2
            slope[parent] += factor / f**2 * (cutoff/math.hypot(cutoff, mass))**3
    return coef, slope, gravity


def potential(cfg, field):
    sig, chi, spectator = field
    ds = sig - cfg['sc']
    v2 = cfg['v']**2
    lam = cfg['lambda']
    g2, gs2, m2, c1, cb = [cfg[k] for k in ('g2', 'gs2', 'm2', 'c1', 'cb')]
    assert sig < cfg['sc'] or cfg['g2'] == cfg['lambda'] == 0
    value = (lam*(chi*chi-v2)**2/4 + g2*sig**2*chi**2/2 + gs2*chi**2*spectator**2/2
             + c1*ds - m2*ds**2/2 - cb*chi + cfg['offset'])
    gradient = np.array([g2*sig*chi**2 + c1 - m2*ds,
                         lam*chi*(chi**2-v2) + g2*sig**2*chi + gs2*chi*spectator**2 - cb,
                         gs2*chi**2*spectator])
    return value, gradient


def unrotated_reference(cfg, times):
    parents, mass, f, width = (np.array(col) for col in zip(*cfg['species'], strict=True))
    parents = parents.astype(int)
    edges = np.linspace(0, cfg['qmax'], cfg['nq']+1)
    q = (edges[:-1] + edges[1:]) / 2
    q_signed = np.r_[-q, q]
    weight = np.tile(np.diff(edges**3) / (6 * math.pi**2), 2)

    def geometry(a, phi):
        angle = phi[parents] / f
        k = np.broadcast_to(q_signed / a, (len(mass), len(q_signed)))
        bx = np.broadcast_to((mass * np.cos(angle))[:, None], k.shape)
        by = np.broadcast_to((-mass * np.sin(angle))[:, None], k.shape)
        b = np.stack([bx, by, k], axis=-1)
        energy = np.hypot(k, mass[:, None])
        return b, energy, angle, k

    def quantities(y):
        a, phi, momentum = y[0], y[1:4], y[4:7]
        rsm, rdr = y[7:9] / a**4
        r = y[10:].reshape(len(mass), len(q_signed), 3)
        A, slope, G = inertias(cfg, a)
        vel = momentum/(a**3 * A)
        b, energy, angle, k = geometry(a, phi)
        axis = b/energy[..., None]
        projection = np.sum(r * axis, axis=-1)
        particles = (projection+1)/2
        species_rho = ((energy * particles) @ weight)/a**3
        species_p = (k*(r[..., 2]+k/energy)) @ weight/(6*a**3)
        V, dV = potential(cfg, phi)
        density = 0.5*(A @ vel**2) + V + sum(species_rho) + rsm + rdr
        pressure = 0.5*(A @ vel**2) + slope @ vel**2/6 - V + sum(species_p) + (rsm+rdr)/3
        H = math.sqrt(density/(3*cfg['mp']**2+G/2)) if cfg['expand'] else 0.0
        return a, r, vel, H, b, axis, projection, particles, energy, angle, dV, species_rho, pressure

    def rhs(t, y):
        a,r,vel,H,b,axis,pr,n,E,theta,dV,_,pressure = quantities(y)
        derivative = -mass[:, None]/f[:, None] * (
            np.sin(theta)[:, None]*r[..., 0] + np.cos(theta)[:, None]*r[..., 1])
        force = np.zeros(3)
        np.add.at(force, parents, (derivative @ weight)/(2*a**3))
        loss_rate = width[:, None]*mass[:, None]/E
        dr = 2*np.cross(b,r) - loss_rate[..., None] * (
            (pr+1)[..., None]*axis + (r-pr[..., None]*axis)/2)
        dmomentum = -a**3*(dV+force)
        dmomentum[2] -= cfg['gamma_s'] * y[6]
        dsm = a * np.sum((loss_rate * n * E) @ weight)
        ddr = a**4 * cfg['gamma_s'] * vel[2]**2
        return np.r_[a*H, vel, dmomentum, dsm, ddr, 3*a**3*H*pressure, dr.ravel()]

    a, phi = cfg['a0'], np.array(cfg['fields'])
    A, _, _ = inertias(cfg, a)
    b, energy, _, _ = geometry(a, phi)
    assert not cfg['dressed_start']
    y0 = np.r_[a, phi, a**3 * A * np.array(cfg['velocities']), a**4*cfg['rho_sm'],
               a**4*cfg['rho_dr'], 0.0, (-b/energy[..., None]).ravel()]
    ode = solve_ivp(rhs, (cfg['t0'], cfg['t1']), y0, t_eval=np.minimum(times, cfg['t1']), rtol=3e-13, atol=1e-14, method='DOP853')
    assert ode.success
    comparison = dict(a=ode.y[0], sigma=ode.y[1], chi=ode.y[2], s=ode.y[3],
                      rho_sm=ode.y[7]/ode.y[0]**4, rho_dr=ode.y[8]/ode.y[0]**4, work=ode.y[9])
    rows = [quantities(ode.y[:, i]) for i in range(len(times))]
    for j, key in enumerate(('dsigma', 'dchi', 'ds')):
        comparison[key] = np.array([r[2][j] for r in rows])
    for j in range(len(mass)):
        comparison[f'rho_{j+1}'] = np.array([r[11][j] for r in rows])
    return comparison


def independent_ode(binary, work):
    results = {}
    for expand in (False, True):
        errors = []
        for steps in (80, 160, 320):
            case = model(expand=expand, nsteps=steps, stride=max(1, steps//80))
            if not expand:
                case['gamma_s'] = 0.0
                case['species'] = [(i,m,f,0.0) for i,m,f,_ in case['species']]
                case['H0'] = 0.0
            h, _, _ = run_engine(binary, work, f'ode_{expand}_{steps}', case)
            ode = unrotated_reference(case, h['t'])
            per_var = {key: float(np.max(abs(h[key]-value))/max(np.max(abs(value)), 1e-5))
                       for key, value in ode.items()}
            errors.append(dict(steps=steps, relative_errors=per_var, max_relative=max(per_var.values()),
                               max_budget_error=float(max(abs(h['residual'])))))
        assert errors[-1]['max_relative'] < 1e-3
        assert errors[-1]['max_relative'] < errors[0]['max_relative'] / 10
        if not expand:
            assert max(row['max_budget_error'] for row in errors) < 1e-10
        results['expanding' if expand else 'flat'] = errors
    return results


def vacuum_controls(binary, work):
    free = dict(**{'lambda': 0.0}, g2=0.0, c1=0.0, m2=0.0, gs2=0.0, cb=0.0, gamma_s=0.0)
    static = model(expand=False, dressed_start=True, H0=0.0, velocities=[0.0,0.0,0.0],
                   species=[(0,0.85,1.3,0.0)], **free)
    h, _, _ = run_engine(binary, work, 'static', static)
    vacuum = float(max(abs(h['rho_1'])))
    assert vacuum < 2e-13
    cfg = dict(static, expand=True, H0=0.23, velocities=[4.0,0.1,0.0])
    empty, _, _ = run_engine(binary, work, 'massless_empty', dict(cfg, species=[]))
    zero, _, _ = run_engine(binary, work, 'massless', dict(cfg, species=[(0,0.0,1.0,0.0)]))
    comparison = max(float(max(abs(zero[key]-empty[key]))) for key in ('a','sigma','chi','rho'))
    assert comparison < 1e-11 and max(abs(zero['rho_1'])) < 1e-13

    rotor = dict(static, dressed_start=False, velocities=[0.89,0.0,0.0])
    h, _, _ = run_engine(binary, work, 'shift_current', rotor)
    shift_error = float(max(abs(h['Q_sigma']-h['Q_sigma'][0])))
    assert shift_error < 1e-11
    assert max(abs(h['dsigma']-h['dsigma'][0])) > 1e-5
    assert max(abs(h['residual'])) < 1e-10

    uv = []
    mass, f = static['species'][0][1:3]
    velocity = 0.9
    for cutoff in (16,32,64,128):
        run = dict(static, a0=1.0, qmax=cutoff, nq=8192, nsteps=4, stride=1,
                   t1=0.01, velocities=[velocity,0.0,0.0])
        h, _, _ = run_engine(binary, work, f'uv_{cutoff}', run)
        mu = velocity/(2*f)

        def integrand(p):
            old = np.array([mass, p]) / math.hypot(p, mass)
            value = 0.0
            for sign in (-1,1):
                axis = np.array([mass,p+sign*mu]) / math.hypot(p+sign*mu,mass)
                value += np.sum((axis-old)**2)/4
            return value * math.hypot(p,mass)*p*p/(2*math.pi**2)

        expected = quad(integrand, 0, cutoff, epsabs=1e-12)[0]
        got = float(h['rho_1'][-1])
        assert abs(got/expected-1)<2e-5
        uv.append(dict(qmax=cutoff, bare_energy=got, adiabatic_energy=float(h['scalar_ct'][-1]),
                       subtracted_energy=float(got-h['scalar_ct'][-1]), quadrature_error=abs(got/expected-1)))
    assert uv[-1]['bare_energy'] > uv[0]['bare_energy'] * 1.2
    assert abs(uv[-1]['subtracted_energy']-uv[-2]['subtracted_energy']) < 2e-6
    return dict(static_vacuum_energy=vacuum, massless_vs_no_species_max_error=comparison,
                canonical_shift_charge_max_error=shift_error, rotating_vacuum_cutoffs=uv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    engine = compile_engine(work)
    outcomes = dict(ode=independent_ode(engine,work), vacuum=vacuum_controls(engine,work))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(outcomes,indent=2,allow_nan=False)+'\n')
    print('Independent ODE, vacuum, massless, shift charge and cutoff controls passed.',flush=True)


if __name__ == '__main__':
    main()
