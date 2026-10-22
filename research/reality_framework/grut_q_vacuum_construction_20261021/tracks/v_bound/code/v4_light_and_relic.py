"""v_bound: independent re-runs of bound A.6 (regime (ii), light q: phi-MDE, K drift) and A.7 (regime (iii) relic toy).

Same formulation as v2_frw_adiabatic.py (N = ln a, H^2 from the constraint), written independently of bound/code.
Part II  family F1 (eps = (e/2)(1+q^2), K = 1 + 0.4 (q-1), m^2 = (2/3) e/0.16 at q0 = 1), mu = mu0 + dmu calibrated so
         that the static de Sitter root has H_s^2 = rho_V/(3 K(q_s)) = 0.685 (H0 = 1); q starts at q_s with dq/dN = 0 at
         a_i = 1e-4; rho_d0 = 3 Om K(q_s), rho_r0 = 3 Or K(q_s); a = 1 is defined by rho_d = rho_d0.
         Diagnostics: K(a_i)/K(a=1), H(a=1), Omega_m(a=1) = rho_d/(3 K_now H^2), p = dln a/dln t at z = 100 (t from a
         radiation-era origin t_i = 1/(2 H_i)).
Part III one fluid with trace fraction kappa (w_f = (1-kappa)/3), mu = mu0, q = 1 at H_i = 100 m or 1000 m, integrated
         to 4 e-folds of H/m^-1/2 after H = m; relic energy (3/4) K m^2 (K' A/K)^2 from the oscillation amplitude A of
         q about the quasi-static solution; scaled back to H = m.
Run: python3 v4_light_and_relic.py  (writes ../results/v4_light_and_relic.json)
"""
import json
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OM, OR = 0.315, 9.1e-5
OL = 1 - OM - OR


def F1(m):
    k = 0.4
    e = 1.5 * m**2 * k**2
    return dict(eps=lambda q: 0.5 * e * (1 + q * q), e1=lambda q: e * q, K=lambda q: 1 + k * (q - 1), K1=lambda q: k + 0 * q,
                K2=lambda q: 0 * q, mu0=e, m=m)


def static_roots(fam, mu):
    g = lambda q: (fam['e1'](q) - mu) * fam['K'](q) - 2 * fam['K1'](q) * (fam['eps'](q) - mu * q)
    grid = np.linspace(-1.4, 30, 200001)          # K > 0 requires q > -1.5
    v = g(grid)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    return [brentq(g, grid[i], grid[i + 1], xtol=1e-15) for i in idx]


def calibrate_light(fam):
    """dmu with the de Sitter root nearest q = 1 giving H_s^2 = OL; scan dmu < 0."""
    sols = []
    for dmu in -np.logspace(-4, 1.3, 600):
        mu = fam['mu0'] + dmu
        for qs in static_roots(fam, mu):
            H2 = (fam['eps'](qs) - mu * qs) / (3 * fam['K'](qs))
            sols.append((dmu, qs, H2))
    # follow the branch through q ~ 1 (smallest |dmu| region continuously connected), pick the first crossing of OL
    best = None
    for (d1, q1, h1), (d2, q2, h2) in zip(sols[:-1], sols[1:]):
        if abs(q2 - q1) < 0.2 and (h1 - OL) * (h2 - OL) < 0:
            best = (d1, d2, q1)
            break
    if best is None:
        raise RuntimeError('no calibration')

    def g(dmu):
        mu = fam['mu0'] + dmu
        r = static_roots(fam, mu)
        qs = min(r, key=lambda x: abs(x - best[2]))
        return (fam['eps'](qs) - mu * qs) / (3 * fam['K'](qs)) - OL
    dmu = brentq(g, best[0], best[1], xtol=1e-14)
    mu = fam['mu0'] + dmu
    qs = min(static_roots(fam, mu), key=lambda x: abs(x - best[2]))
    return dmu, mu, qs


def make_rhs(fam, mu, rd0, rr0, kappa=None, wf=None):
    def U(q):
        return fam['eps'](q) - mu * q

    def Rq(q):
        return 2 * (fam['e1'](q) - mu) / fam['K1'](q)

    def dens(N):
        if kappa is None:
            return rd0 * np.exp(-3 * N), rr0 * np.exp(-4 * N), rd0 * np.exp(-3 * N)      # total parts, trace source
        rf = rd0 * np.exp(-3 * (1 + wf) * N)
        return rf, 0.0, kappa * rf

    def H2(N, q, p):
        a_, b_, _ = dens(N)
        return (a_ + b_ + U(q)) / (3 * (fam['K'](q) + fam['K1'](q) * p))

    def rhs(N, y):
        q, p = y[0], y[1]
        h2 = H2(N, q, p)
        R = Rq(q)
        _, _, tr = dens(N)
        S = (tr + 4 * U(q) - fam['K'](q) * R) / 3
        pp = (S / h2 - fam['K2'](q) * p * p - fam['K1'](q) * p * (R / (6 * h2) + 1)) / fam['K1'](q)
        return [p, pp, 1.0 / np.sqrt(h2)]           # third component: t (dt/dN = 1/H)
    return rhs, H2, U, Rq


def light_run(m, a_i=1e-4):
    fam = F1(m)
    dmu, mu, qs = calibrate_light(fam)
    Ks = fam['K'](qs)
    rd0, rr0 = 3 * OM * Ks, 3 * OR * Ks
    rhs, H2, U, Rq = make_rhs(fam, mu, rd0, rr0)
    Ni = np.log(a_i)
    Hi = np.sqrt(H2(Ni, qs, 0.0))
    per_N = 2 * np.pi * 1.0 / m          # crude bound on oscillation scale in N at H ~ 1
    sol = solve_ivp(rhs, (Ni, 0.0), [qs, 0.0, 1 / (2 * Hi)], method='DOP853', rtol=1e-11, atol=1e-13,
                    max_step=min(0.01, per_N / 40), dense_output=True)
    q1, p1, t1 = sol.y[:, -1]
    Kn = fam['K'](q1)
    H1 = np.sqrt(H2(0.0, q1, p1))
    N100 = -np.log(101.0)
    q100, p100, t100 = sol.sol(N100)
    H100 = np.sqrt(H2(N100, q100, p100))
    return {'m/H0': m, 'dmu': dmu, 'q_s': qs, 'K(q_s)': Ks, 'H(a=1)': float(H1), 'K(a_i)/K(a=1)': float(Ks / Kn),
            'Omega_m(a=1) (K_now)': float(rd0 / (3 * Kn * H1**2)), 'p = H t at z=100': float(H100 * t100),
            'status': int(sol.status)}


def relic_run(m, kappa, Hi_over_m, n_after=4.0):
    fam = F1(m)
    mu = fam['mu0']
    wf = (1 - kappa) / 3
    rf_i = 3 * (Hi_over_m * m)**2            # K(q0) = 1, U(q0) = 0 at mu0
    rhs, H2, U, Rq = make_rhs(fam, mu, rf_i, 0.0, kappa=kappa, wf=wf)

    def ev(N, y):
        return np.sqrt(H2(N, y[0], y[1])) - m * np.exp(-2 * n_after)
    ev.terminal = True
    sol = solve_ivp(rhs, (0.0, 60.0), [1.0, 0.0, 1 / (2 * Hi_over_m * m)], method='LSODA', rtol=1e-10, atol=1e-13,
                    events=[ev], dense_output=True, max_step=0.01)
    Nend = sol.t[-1]
    # N where H = m
    Ng = np.linspace(0, Nend, 20001)
    qg, pg, _ = sol.sol(Ng)
    Hg = np.sqrt(H2(Ng, qg, pg))
    Nm = float(np.interp(-m, -Hg, Ng))
    # last stretch: about 3 oscillation periods in N
    Hend = np.sqrt(H2(Nend, *sol.y[:2, -1]))
    dNper = 2 * np.pi * Hend / m
    Nw = np.linspace(Nend - 3 * dNper, Nend, 3001)
    qw, pw, _ = sol.sol(Nw)
    qsl = np.empty_like(Nw)
    for i, N in enumerate(Nw):
        rf = rf_i * np.exp(-3 * (1 + wf) * N)
        qsl[i] = brentq(lambda x: kappa * rf + 4 * U(x) - fam['K'](x) * Rq(x), 0.5, 1.5, xtol=1e-15)
    A = 0.5 * (np.max(qw - qsl) - np.min(qw - qsl))
    qe = qw[-1]
    rf_end = rf_i * np.exp(-3 * (1 + wf) * Nend)
    E = 0.75 * fam['K'](qe) * m**2 * (fam['K1'](qe) * A / fam['K'](qe))**2
    frac_end = E / rf_end
    growth = np.exp((Nend - Nm) * 3 * wf)                # relic/fluid ~ a^-3 / a^-3(1+w_f) = a^(3 w_f)
    frac_on = frac_end / growth
    return {'m': m, 'kappa': kappa, 'H_i/m': Hi_over_m, 'N(H=m)': Nm, 'N_end': float(Nend),
            'E_osc/rho_f at end': float(frac_end), 'relic fraction at onset': float(frac_on),
            'coefficient = fraction/(kappa^2 ln^2(H_i/m))': float(frac_on / (kappa**2 * np.log(Hi_over_m)**2))}


def main():
    t0 = time.time()
    out = {}
    II = []
    for m in [0.1, 1.0, 10.0, 100.0]:
        r = light_run(m)
        II.append(r)
        print('II', r, flush=True)
    out['partII_light'] = II
    III = []
    for kap in [1e-3, 1e-2]:
        for hi in [1e2, 1e3]:
            r = relic_run(10.0, kap, hi)
            III.append(r)
            print('III', r, flush=True)
    out['partIII_relic'] = III
    out['runtime_s'] = time.time() - t0
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'v4_light_and_relic.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
