#!/usr/bin/env python3
"""Affleck-Dine charge generation for the Codex/GRUT complex scalar X (own code).

Conditional model calculation, not a GRUT prediction.

During the AD epoch |X| >> v, the binding field is pinned at phi=0 (its mass
sqrt(2) g|X| >> H, omega), so the X potential is
    V = lambda_X |X|^4 + (lambda_d X^d + h.c.),   lambda_d = c_d / M^(d-4),
plus a supplied Hubble-induced mass -c_H H^2 |X|^2 that only fixes the release
amplitude.  Conformal field chi = a X, conformal time tau:
    chi'' - (a''/a) chi + 2 lambda_X |chi|^2 chi + d lambda_d^* chi*^(d-1) a^(4-d) = 0.
Dimensionless: s = omega_0 tau, u = chi / X_0, omega_0 = sqrt(2 lambda_X) X_0, a=1 at release:
    u'' = (a''/a) u - |u|^2 u - (d/2) eta_d e^{-i delta} conj(u)^(d-1) a^(4-d),
    eta_d = |c_d| X_0^(d-2) / (lambda_X M^(d-4)).
Comoving charge (Codex sign convention Q = -2 Im X* Xdot):
    n a^3 = q X_0^2 omega_0,   q = -2 Im(conj(u) u').
Release parameter h = H_rel / omega_0 (h >> 1: frozen field; h = 1/sqrt(c_H):
release of a field held at the Hubble-induced minimum).  Initial radial
velocity Xdot/X = -beta H_rel.  The A-term is V_A = 2|lambda_d| R^d cos(phi),
phi = d*theta + delta; phi0 is the supplied initial angle (pi = valley).
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parents[1] / 'results'


def run(d, eta, phi0, h, beta=0.0, era='RD', a_end=None, rtol=1e-10, atol=1e-12, eta_sign=1.0):
    theta0 = phi0 / d
    if era == 'RD':
        s_r = 1.0 / h
        a = lambda s: s / s_r
        app = lambda s: 0.0
    else:  # MD: a = (s/s_r)^2, H_rel = 2/tau_r
        s_r = 2.0 / h
        a = lambda s: (s / s_r) ** 2
        app = lambda s: 2.0 / s ** 2
    if a_end is None:
        a_end = 1e3
    s_end = s_r * (a_end if era == 'RD' else np.sqrt(a_end))
    coef = 0.5 * d * eta * eta_sign

    def rhs(s, y):
        u = y[0] + 1j * y[1]
        aa = a(s)
        acc = app(s) * u - abs(u) ** 2 * u - coef * np.conj(u) ** (d - 1) * aa ** (4 - d)
        return [y[2], y[3], acc.real, acc.imag]

    u0 = np.exp(1j * theta0)
    du0 = h * (1.0 - beta) * u0
    track = {'max_X_over_X0': 0.0}

    def runaway(s, y):  # physical |X|/X0 > 3 means the A-term valley is running away
        return abs(y[0] + 1j * y[1]) / a(s) - 3.0
    runaway.terminal = True
    sol = solve_ivp(rhs, (s_r, s_end), [u0.real, u0.imag, du0.real, du0.imag], method='DOP853',
                    rtol=rtol, atol=atol, events=runaway, dense_output=False, max_step=np.inf)
    y = sol.y
    u = y[0] + 1j * y[1]
    du = y[2] + 1j * y[3]
    q = -2 * np.imag(np.conj(u) * du)
    av = np.array([a(s) for s in sol.t])
    xphys = np.abs(u) / av
    ran_away = sol.status == 1
    # late-time average over the last 20% of samples (the residual A-term oscillation decays)
    n = len(q)
    tail = slice(int(0.8 * n), n)
    q_inf = float(np.mean(q[tail]))
    q_spread = float(np.std(q[tail]))
    E = np.abs(du[-1]) ** 2 + 0.5 * np.abs(u[-1]) ** 4
    qmax = 2 * (2 * E / 3) ** 0.75
    return dict(d=d, eta=eta, phi0=phi0, h=h, beta=beta, era=era, a_end=a_end,
                q_inf=q_inf, q_tail_std=q_spread, epsilon=float(abs(q_inf) / qmax),
                max_X_over_X0_after_release=float(xphys.max()), ran_away=bool(ran_away),
                success=bool(sol.success or ran_away), nfev=int(sol.nfev))


def _job(args):
    return run(*args[:5], era=args[5])


def main():
    t0 = time.time()
    cases = [
        ('frozen_RD', 10.0, 0.0, 'RD'), ('frozen_RD_h30', 30.0, 0.0, 'RD'),
        ('frozen_MD', 10.0, 0.0, 'MD'), ('frozen_MD_h30', 30.0, 0.0, 'MD'),
        ('hubble_RD_cH1_b0', 1.0, 0.0, 'RD'), ('hubble_RD_cH1_b1.5', 1.0, 1.5, 'RD'),
        ('hubble_RD_cH10_b0', 10 ** -0.5, 0.0, 'RD'), ('hubble_RD_cH100_b0', 0.1, 0.0, 'RD'),
        ('hubble_MD_cH1_b0', 1.0, 0.0, 'MD'), ('hubble_MD_cH10_b0', 10 ** -0.5, 0.0, 'MD'),
        ('hubble_MD_cH100_b0', 0.1, 0.0, 'MD'),
    ]
    small = [1e-5, 1e-4, 1e-3]
    large = [0.01, 0.03, 0.1, 0.2, 0.3, 0.4]
    phis = [k * np.pi / 8 for k in range(1, 9)]  # pi/8 ... pi
    jobs, labels = [], []
    for d in (6, 5):
        for name, h, beta, era in cases:
            if name.endswith('h30'):
                grid = [(e, np.pi / 2) for e in small + [0.1]]
            else:
                grid = [(e, np.pi / 2) for e in small] + [(e, p) for e in large for p in phis]
            for eta, phi in grid:
                jobs.append((d, eta, phi, h, beta, era)); labels.append(name)
    with Pool(4) as pool:
        rows = pool.map(_job, jobs, chunksize=4)
    for r, lab in zip(rows, labels):
        r['case'] = lab
    # controls: eta=0 gives exactly zero charge; tolerance halving; longer integration
    controls = dict(
        zero_eta=run(6, 0.0, np.pi / 2, 1.0),
        tol_6=[run(6, 0.1, np.pi / 2, 1.0, rtol=r, atol=r * 1e-2) for r in (1e-9, 1e-11)],
        tol_5=[run(5, 0.1, np.pi / 2, 1.0, rtol=r, atol=r * 1e-2) for r in (1e-9, 1e-11)],
        long_6=run(6, 0.1, np.pi / 2, 1.0, a_end=1e4),
        long_5=run(5, 0.1, np.pi / 2, 1.0, a_end=1e4),
        mirror_6=run(6, 0.1, -np.pi / 2, 1.0),
    )
    out = dict(scope='Own AD condensate solver; conditional model calculation, not a GRUT prediction.',
               equation='u\'\' = (a\'\'/a)u - |u|^2 u - (d/2) eta e^{-i delta} conj(u)^(d-1) a^(4-d); q=-2Im(conj(u)u\')',
               rows=rows, controls=controls, runtime_s=time.time() - t0)
    OUT.mkdir(exist_ok=True)
    (OUT / 'ad_scan.json').write_text(json.dumps(out, indent=1))
    print('rows', len(rows), 'time', time.time() - t0)


if __name__ == '__main__':
    main()
