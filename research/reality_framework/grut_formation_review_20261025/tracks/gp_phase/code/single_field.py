#!/usr/bin/env python3
"""gp_phase (c), reduced comparator: ONE complex field, global U(1) (no Gauss constraint, no Y,
no Sigma), with two Z_6-class terms of different amplitude scaling and an irremovable relative
phase delta:

    V = lambda |X|^4 + 2 Re(a X^6) + 2 Re(b |X|^2 X^6)   [ - c_H H^2 |X|^2 before release ]

This is the moving-minimum hypothesis in its purest form, an UPPER-ENVELOPE comparator: it omits
the gauge-point physics (charge compensation by Y/Sigma, multi-field chaos). The field starts at
rest at its own inflationary minimum (c_H = 1); release into radiation with the Hubble mass off.
Dimensionless (same conventions as tracks/affleck_dine/code/ad_condensate.py):
    u'' = -|u|^2 u - 3 alpha* conj(u)^5 / a^2 - (beta u^7 + 7 beta* u conj(u)^6) / (2 a^4),
    alpha = a X0^2/lambda = eta (real > 0 by rephasing),  beta = r0 eta exp(i delta),
    X0 = sqrt(c_H/(2 lambda)) H,  s = omega0 tau,  a = s (h = 1).
At |u| = 1 the X^6 term has relative weight eta (AD-track eta_6) and the second term r0 times that.
For |X| decreasing the minimum of the angular potential moves from 6 theta = pi + Delta(r0, delta)
to 6 theta = pi. eps = |q| / (2 (2E/3)^(3/4)), q = -2 Im(conj(u) u'), E = |u'|^2 + |u|^4/2.
"""
import json, sys, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parents[1] / 'results'
EPS_REQ = (1 + 0.01 / 0.9) ** (-0.75)


def vinf(u, alpha, beta):
    r2 = abs(u) ** 2
    return r2 * r2 - 2 * r2 + 2 * np.real(alpha * u ** 6) + 2 * np.real(beta * r2 * u ** 6)


def initial_minimum(alpha, beta):
    """Lowest interior local minimum near |u| = 1 (the X^6 terms make V unbounded at large |u|).
    Grid search over R in [0.2, 1.8] and one full period of theta (pi/3), then a local polish."""
    R = np.linspace(0.2, 1.8, 321); T = np.linspace(-np.pi / 6, np.pi / 6, 240, endpoint=False)
    RR, TT = np.meshgrid(R, T, indexing='ij')
    V = vinf(RR * np.exp(1j * TT), alpha, beta)
    core = V[1:-1]
    nb = [np.roll(V, s, axis=1)[1 + dr:V.shape[0] - 1 + dr] for dr in (-1, 0, 1) for s in (-1, 0, 1) if (dr, s) != (0, 0)]
    mask = np.all([core < n for n in nb], axis=0)
    idx = np.argwhere(mask)
    if len(idx) == 0:
        return None
    best = None
    for i, j in idx:
        x0 = [R[i + 1], T[j]]
        r = minimize(lambda v: vinf(v[0] * np.exp(1j * v[1]), alpha, beta), x0, method='Nelder-Mead',
                     options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
        Rm, tm = r.x
        if 0.21 < Rm < 1.79 and abs(Rm - x0[0]) < 0.1 and (best is None or r.fun < best[1]):
            best = (Rm * np.exp(1j * tm), r.fun)
    return best[0] if best else None


def run(eta, r0, delta, a_end=1000.0, rtol=1e-10, atol=1e-12):
    alpha = eta; beta = r0 * eta * np.exp(1j * delta)
    u0 = initial_minimum(alpha, beta)
    if u0 is None:
        return dict(eta=eta, r0=r0, delta=delta, status='no local minimum')
    ac, bc = np.conj(alpha), np.conj(beta)

    def rhs(s, y):
        u = y[0] + 1j * y[1]; a = s
        acc = -abs(u) ** 2 * u - 3 * ac * np.conj(u) ** 5 / a ** 2 - (beta * u ** 7 + 7 * bc * u * np.conj(u) ** 6) / (2 * a ** 4)
        return [y[2], y[3], acc.real, acc.imag]

    def runaway(s, y):
        return abs(y[0] + 1j * y[1]) / s - 3.0
    runaway.terminal = True
    du0 = u0  # zero physical velocity: u' = h u with h = 1
    sol = solve_ivp(rhs, (1.0, a_end), [u0.real, u0.imag, du0.real, du0.imag], method='DOP853', rtol=rtol, atol=atol,
                    events=runaway, t_eval=np.geomspace(1.0, a_end, 400))
    u = sol.y[0] + 1j * sol.y[1]; du = sol.y[2] + 1j * sol.y[3]
    q = -2 * np.imag(np.conj(u) * du)
    E = np.abs(du) ** 2 + 0.5 * np.abs(u) ** 4
    eps = q / (2 * (2 * E / 3) ** 0.75)
    n = len(q); tail = slice(int(0.85 * n), n)
    return dict(eta=eta, r0=r0, delta=delta, status='runaway' if sol.status == 1 else 'ok',
                initial_angle_6theta=float(np.angle(u0 ** 6)), initial_abs_u=float(abs(u0)),
                eps_final=float(np.mean(eps[tail])), eps_tail_std=float(np.std(eps[tail])),
                eps_at_a100=float(np.interp(100.0, sol.t, eps)) if sol.t[-1] >= 100 else None)


def _job(a):
    return run(*a)


def main():
    etas = [0.005, 0.01, 0.02, 0.04, 0.06, 0.08]  # Hubble-held minimum exists only for eta < 1/12
    r0s = [0.03, 0.1, 0.3, 1.0, 3.0]
    deltas = list(np.linspace(np.pi / 16, np.pi, 16))
    jobs = [(e, r, d) for e in etas for r in r0s for d in deltas]
    with Pool(4) as pool:
        rows = pool.map(_job, jobs, chunksize=4)
    ok = [r for r in rows if r['status'] == 'ok']
    best = max(ok, key=lambda r: abs(r['eps_final']))
    # refine around the best point
    e0, r00, d0 = best['eta'], best['r0'], best['delta']
    fine = [(e, r, d) for e in np.linspace(max(0.005, e0 * 0.6), min(0.083, e0 * 1.4), 9)
            for r in np.geomspace(r00 / 2, r00 * 2, 9) for d in np.linspace(d0 - np.pi / 16, d0 + np.pi / 16, 9)]
    with Pool(4) as pool:
        frows = pool.map(_job, fine, chunksize=8)
    fok = [r for r in frows if r['status'] == 'ok']
    fbest = max(fok, key=lambda r: abs(r['eps_final']))
    # controls: tolerance and a_end
    ctrl = dict(tight=run(fbest['eta'], fbest['r0'], fbest['delta'], rtol=1e-12, atol=1e-14),
                long=run(fbest['eta'], fbest['r0'], fbest['delta'], a_end=1e4),
                cp_mirror=run(fbest['eta'], fbest['r0'], -fbest['delta']),
                no_second_term=run(fbest['eta'], 0.0, fbest['delta']))
    out = dict(scope='Reduced one-field global-U(1) comparator (upper envelope); omits gauge-point Y/Sigma physics.',
               eps_requirement=EPS_REQ, coarse=rows, n_coarse_runaway=sum(r['status'] == 'runaway' for r in rows),
               coarse_best=best, fine_best=fbest, n_fine=len(frows), n_fine_runaway=sum(r['status'] == 'runaway' for r in frows),
               controls=ctrl,
               max_abs_eps_by_eta={repr(e): max([abs(r['eps_final']) for r in ok if r['eta'] == e] or [float('nan')]) for e in etas})
    (OUT / 'single_field.json').write_text(json.dumps(out, indent=1, default=float))
    print('coarse best', best); print('fine best', fbest)
    print('controls', {k: (v.get('eps_final'), v.get('status')) for k, v in ctrl.items()})
    print('max|eps| by eta', out['max_abs_eps_by_eta'], 'runaways', out['n_coarse_runaway'])


if __name__ == '__main__':
    main()
