"""v_scenario / v2: independent full-equation integration of the compressed-hierarchy toy (scenario NOTES §2).

Independent of the scenario code in three ways:
  * H is NOT integrated from Hdot; it is solved from the Friedmann constraint at every step
        3 H^2 (1+y) + 3 H ydot = rho_m + rho_V(y)
    (so the result does not rely on Bianchi propagation, which v1 checks symbolically);
  * initial data and the quasi-static track used to measure the oscillation energy come from my own
    high-order adiabatic series (v2_qs_series.py, coupled y and H^2 expansion to O(m^-10)), not from the
    scenario's O(m^-4) track (a local least-squares fit of y(t) was tried first and found biased at the 1e-2
    level by the a^{-3/2} amplitude drift; it is kept as fit_osc for reference only);
  * the KV q-sector functions are taken from the sympy derivation of v1 (check C), K0 = 1, time unit 1/H_c.
Toy matter (same as the scenario's §1.5 toy so the numbers can be compared): radiation + dust + matter-vacuum step
DeltaV(N) = Delta (1 - tanh(N/w))/2, latent heat to radiation, Delta = 3 fD, initial densities at N_i = -1:
rho_r = 3(1-fdust-fD) e^{4}, rho_d = 3 fdust e^{3}; fD = 0.05, fdust = 0.02.

Outputs ../results/v2_toy_full.json:
  (1) reduction error eps0 = <H>_N/<H_R0>_N - 1 and eps1 (against R1) for m = 10..1000 (w = 1);
  (2) deposited oscillation energy / [(4/3) Delta^2/m^2 (x/sinh x)^2] for x = 0.1, 1, 2, 4 at m = 100, 300;
  (3) Delta = 0 controls.
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FD, FDUST, NI = 0.05, 0.02, -1.0
GLX, GLW = np.polynomial.legendre.leggauss(24)


class Toy:
    def __init__(self, m, w, Delta_on=True, fD=FD, fdust=FDUST, Ni=NI):
        self.m = float(m)
        self.m2 = self.m ** 2
        self.w = float(w)
        self.D = 3 * fD if Delta_on else 0.0
        self.Ni = Ni
        self.rr_i = 3 * (1 - fdust - fD) * math.exp(-4 * Ni)
        self.rd_i = 3 * fdust * math.exp(-3 * Ni)

    # matter vacuum step and its N-derivatives
    def dV(self, N):
        return self.D * 0.5 * (1 - math.tanh(N / self.w))

    def dV1(self, N):
        ex = math.exp(-2 * abs(N / self.w))
        return -self.D * 0.5 / self.w * 4 * ex / (1 + ex) ** 2

    def rho_d(self, N):
        return self.rd_i * math.exp(-3 * (N - self.Ni))

    # q-sector (KV quartic, from v1 check C): rho_V = (3/4) m^2 y^2 (1 + 4y/3 + y^2/3); R = 3 m^2 y (1 + 2y + 2y^2/3)
    def rhoV(self, y):
        return 0.75 * self.m2 * y * y * (1 + 4 * y / 3 + y * y / 3)

    def R(self, y):
        return 3 * self.m2 * y * (1 + 2 * y + 2 * y * y / 3)

    def H(self, y, u, rho):
        A = 3 * (1 + y)
        return (-3 * u + math.sqrt(9 * u * u + 4 * A * (rho + self.rhoV(y)))) / (2 * A)

    def rhs(self, t, s):
        y, u, N, rr, J = s
        dv = self.dV(N)
        rd = self.rho_d(N)
        rho = rr + rd + dv
        S = rd + 4 * dv
        H = self.H(y, u, rho)
        ud = (S + 4 * self.rhoV(y) - (1 + y) * self.R(y)) / 3 - 3 * H * u
        return [u, ud, H, H * (-4 * rr - self.dV1(N)), H * H]


def matter_of_N(toy, Ngrid):
    """rho_r(N) by integrating d rho_r/dN = -4 rho_r - DeltaV'(N) (model independent)."""
    s = solve_ivp(lambda N, r: [-4 * r[0] - toy.dV1(N)], (toy.Ni, Ngrid[-1]), [toy.rr_i], method='DOP853',
                  rtol=1e-13, atol=1e-30, dense_output=True, max_step=min(0.01, toy.w / 8))
    return s.sol


def H_red(toy, rr_sol, N, order):
    rr = float(rr_sol(N)[0])
    rd = toy.rho_d(N)
    dv = toy.dV(N)
    rho = rr + rd + dv
    if order == 0:
        return math.sqrt(rho / 3)
    S = rd + 4 * dv
    SN = -3 * rd + 4 * toy.dV1(N)
    return math.sqrt((rho + S * S / (12 * toy.m2)) / (3 * (1 + (S + SN) / (3 * toy.m2))))


_QS = None


def qs_series(order=5):
    global _QS
    if _QS is None:
        import v2_qs_series
        _QS = v2_qs_series.build(order)
    return _QS


_QS3 = None


def qs_series_post():
    global _QS3
    if _QS3 is None:
        import v2_qs_series
        _QS3 = v2_qs_series.build(3)
    return _QS3


def initial_state(toy, qs=None):
    """Quasi-static start from the O(m^-10) series of v2_qs_series (my construction), evaluated with mpmath.
    Matter derivatives at N_i: S = B e^{-3N} + 4 dV, rho = R + B e^{-3N} + dV, R' = -4R - dV'."""
    import mpmath as mp
    import sympy as sp
    mp.mp.dps = 40
    yk, ykN, pk, nder, _ = qs if qs is not None else qs_series()
    Nn = sp.symbols('N')
    dVs = toy.D * sp.Rational(1, 2) * (1 - sp.tanh(Nn / sp.Float(toy.w, 30)))
    Ni = mp.mpf(toy.Ni)
    dVd = [mp.mpf(sp.N(dVs.diff(Nn, j).subs(Nn, sp.Float(toy.Ni, 30)), 35)) for j in range(nder)]
    B = mp.mpf(toy.rd_i) * mp.exp(3 * Ni)
    Rd = [mp.mpf(toy.rr_i)]
    for j in range(nder - 1):
        Rd.append(-4 * Rd[j] - dVd[j + 1])
    Sd = [B * (-3) ** j * mp.exp(-3 * Ni) + 4 * dVd[j] for j in range(nder)]
    rd = [Rd[j] + B * (-3) ** j * mp.exp(-3 * Ni) + dVd[j] for j in range(nder)]
    e = 1 / mp.mpf(toy.m2)
    y = sum(e ** (k + 1) * yk[k](Sd, rd) for k in range(len(yk)))
    yN = sum(e ** (k + 1) * ykN[k](Sd, rd) for k in range(len(ykN)))
    P = sum(e ** k * pk[k](Sd, rd) for k in range(len(pk)))
    terms = [float(abs(e ** (k + 1) * yk[k](Sd, rd))) for k in range(len(yk))]
    u = mp.sqrt(P) * yN
    return [float(y), float(u), float(Ni), float(toy.rr_i), 0.0], terms


def fit_osc(ts, ys, m):
    s = ts - ts.mean()
    A = np.vstack([np.ones_like(s), s, s * s, s ** 3, np.cos(m * s), s * np.cos(m * s), np.sin(m * s),
                   s * np.sin(m * s)]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    return coef[4], coef[6]


def run(m, w, Delta_on=True, Ni=NI, Nf=1.5, Nmin=-0.25, rtol=1e-12, n_win=60, k=None, fit_periods=4,
        s0=None, Efit_from=None, fD=FD):
    t0 = time.time()
    toy = Toy(m, w, Delta_on, Ni=Ni, fD=fD)
    rr_sol = matter_of_N(toy, [Nf + 0.2])
    terms = None
    if s0 is None:
        s0, terms = initial_state(toy)
    ev = lambda t, s: s[2] - Nf
    ev.terminal = True
    sol = solve_ivp(toy.rhs, (0, 1e6), s0, method='DOP853', rtol=rtol, atol=1e-30, dense_output=True,
                    events=ev, max_step=min(0.25 * 2 * math.pi / m, 1.0))
    tf = sol.t[-1]
    period = 2 * math.pi / m
    if k is None:
        k = max(1, int(round(0.002 * m)))
    win = k * period
    # windows log-spaced in t
    centres = np.exp(np.linspace(math.log(3 * win), math.log(tf - win), n_win))
    rows = []
    last_hi = -1
    for c in centres:
        lo, hi = c - win / 2, c + win / 2
        if lo <= last_hi:
            continue
        last_hi = hi
        s_lo, s_hi = sol.sol(lo), sol.sol(hi)
        N1, N2 = s_lo[2], s_hi[2]
        Hbar = (s_hi[4] - s_lo[4]) / (N2 - N1)
        nodes = 0.5 * (N1 + N2) + 0.5 * (N2 - N1) * GLX
        H0 = 0.5 * sum(wi * H_red(toy, rr_sol, n, 0) for wi, n in zip(GLW, nodes))
        H1 = 0.5 * sum(wi * H_red(toy, rr_sol, n, 1) for wi, n in zip(GLW, nodes))
        rows.append(dict(N=0.5 * (N1 + N2), eps0=Hbar / H0 - 1, eps1=Hbar / H1 - 1, pred=H1 / H0 - 1))
    rows_sel = [r for r in rows if r['N'] >= Nmin]
    # oscillation energy E a^3 after the step: Einstein-frame energy about my own quasi-static series (O(m^-6)),
    # averaged over one period; post-step matter is pure radiation + dust (analytic N-derivatives)
    import mpmath as mp
    yk, ykN, pk, nder, _ = qs_series_post()
    B = toy.rd_i * math.exp(3 * toy.Ni)

    def qs_at(Nv, rr):
        Sd = [mp.mpf(B) * (-3) ** j * mp.exp(-3 * Nv) for j in range(nder)]
        rd_ = [mp.mpf(rr) * (-4) ** j + mp.mpf(B) * (-3) ** j * mp.exp(-3 * Nv) for j in range(nder)]
        e = 1 / mp.mpf(toy.m2)
        yv = sum(e ** (k + 1) * yk[k](Sd, rd_) for k in range(len(yk)))
        yN = sum(e ** (k + 1) * ykN[k](Sd, rd_) for k in range(len(ykN)))
        P = sum(e ** k * pk[k](Sd, rd_) for k in range(len(pk)))
        return float(yv), float(mp.sqrt(P) * yN)
    if Efit_from is None:
        Efit_from = 6 * w + 0.2
    tt = np.linspace(0, tf, 400001)
    Ns = sol.sol(tt)[2]
    Efits = []
    for Nt in (np.linspace(Efit_from, Nf - 0.05, 5) if Efit_from < Nf - 0.1 else []):
        tc = float(np.interp(Nt, Ns, tt))
        ts = tc + period * (np.arange(32) / 32 - 0.5)
        ss = sol.sol(ts)
        Es = []
        for j in range(len(ts)):
            yq, uq = qs_at(ss[2][j], ss[3][j])
            Es.append(0.75 * ((ss[1][j] - uq) ** 2 + toy.m2 * (ss[0][j] - yq) ** 2) * math.exp(3 * ss[2][j]))
        Efits.append(dict(N=float(Nt), E_a3=float(np.mean(Es))))
    # H at the step centre from the full solution (for the x convention)
    t0c = float(np.interp(0.0, Ns, tt))
    sc = sol.sol(t0c)
    H_full_c = toy.H(sc[0], sc[1], sc[3] + toy.rho_d(sc[2]) + toy.dV(sc[2]))
    Hc = H_red(toy, rr_sol, 0.0, 0)                 # R0 value at N = 0
    Hc1 = H_red(toy, rr_sol, 0.0, 1)                # R1 value (the scenario's convention)
    x = math.pi * m * (w / Hc) / 2                  # tau = w/H_c
    FF = lambda xx: (2 * xx * math.exp(-xx) / (1 - math.exp(-2 * xx))) ** 2
    x1 = math.pi * m * (w / Hc1) / 2
    xf = math.pi * m * (w / H_full_c) / 2
    Ef = (4.0 / 3.0) * toy.D ** 2 / toy.m2 * FF(x) if toy.D else float('nan')
    Ea3 = np.array([e['E_a3'] for e in Efits]) if Efits else np.array([float('nan')])
    res = dict(m=m, w=w, Delta=toy.D, H_c=Hc, m_over_Hc=m / Hc, x=x, k=k, rtol=rtol, nfev=int(sol.nfev),
               t_f=float(tf), wall_s=time.time() - t0,
               max_abs_eps0=float(max(abs(r['eps0']) for r in rows_sel)),
               max_abs_eps1=float(max(abs(r['eps1']) for r in rows_sel)),
               max_abs_pred=float(max(abs(r['pred']) for r in rows_sel)), Ni=Ni, Nf=Nf, Nmin=Nmin,
               qs_series_terms=terms,
               E_a3_mean=float(Ea3.mean()), E_a3_rel_spread=float(Ea3.std() / Ea3.mean()) if Ea3.mean() else None,
               E_formula=Ef, E_over_formula=float(Ea3.mean() / Ef) if toy.D else None,
               E_over_formula_Hc_R1=float(Ea3.mean() / Ef * FF(x) / FF(x1)) if toy.D else None,
               E_over_formula_Hc_full=float(Ea3.mean() / Ef * FF(x) / FF(xf)) if toy.D else None,
               H_c_R0=Hc, H_c_R1=Hc1, H_c_full=H_full_c,
               E_over_sudden=float(Ea3.mean() / ((4.0 / 3.0) * (3 * fD) ** 2 / toy.m2)),
               rows=rows, Efits=Efits)
    return res


def job(kw):
    return run(**kw)


if __name__ == '__main__':
    qs = qs_series()
    jobs = []
    # (1) reduction error, adiabatic step w = 1 (same N_i as the scenario: -0.25 (m=10), -0.6 (30), -1 otherwise)
    NImap = {10: -0.25, 30: -0.6}
    for m in (10, 30, 100, 300, 1000, 3000):
        jobs.append(dict(m=m, w=1.0, Ni=NImap.get(m, -1.0), Nf=1.5))
    # (2) deposit vs x at m = 100, 300 (tau = w/H_c, H_c = H(N=0) = 0.98742 for this toy)
    for m in (100, 300):
        for x in (0.1, 1.0, 2.0, 4.0, 8.0):
            w = 2 * x * 0.98742 / (math.pi * m)
            jobs.append(dict(m=m, w=w, Ni=-0.7, Nf=max(1.0, 6 * w + 0.9), Efit_from=max(0.4, 6 * w + 0.2)))
    # Delta^2 scaling at m = 100, x = 1
    for fD in (0.2, 0.0125):
        jobs.append(dict(m=100, w=2 * 0.98742 / (math.pi * 100), Ni=-0.7, Nf=1.0, Efit_from=0.4, fD=fD))
    # (3) controls (Delta = 0)
    for m in (100, 300):
        jobs.append(dict(m=m, w=2.0 / (math.pi * m), Ni=-0.7, Nf=1.0, Efit_from=0.4, Delta_on=False))
    for j in jobs:
        toy = Toy(j['m'], j['w'], j.get('Delta_on', True), Ni=j['Ni'], fD=j.get('fD', FD))
        j['s0'], _ = initial_state(toy, qs)
    jobs.sort(key=lambda j: -j['m'] * (2 if j['w'] == 1.0 else 1))
    t0 = time.time()
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, jobs))
    out = dict(script='v2_toy_full.py', wall_s=time.time() - t0, runs=res)
    sel = [r for r in res if r['w'] == 1.0 and r['Delta'] > 0 and r['m'] >= 100]
    lx = np.log([r['m_over_Hc'] for r in sel])
    ly = np.log([r['max_abs_eps0'] for r in sel])
    sl, ic = np.polyfit(lx, ly, 1)
    out['eps0_fit'] = dict(slope=float(sl), coeff=float(math.exp(ic)), points=len(sel))
    here = os.path.dirname(os.path.abspath(__file__))
    for r in res:
        r.pop('rows_full', None)
    with open(os.path.join(here, '..', 'results', 'v2_toy_full.json'), 'w') as f:
        json.dump(out, f, indent=1)
    for r in res:
        print('m=%6g w=%.4g D=%.4f m/Hc=%.1f x=%.3g  eps0=%.4g pred=%.4g eps1=%.3g  E/formula(R0,R1)=%s,%s  E/sudden=%.3g '
              'spread=%s nfev=%d %.0fs' % (r['m'], r['w'], r['Delta'], r['m_over_Hc'], r['x'], r['max_abs_eps0'],
                                          r['max_abs_pred'], r['max_abs_eps1'], r['E_over_formula'], r['E_over_formula_Hc_R1'],
                                          r['E_over_sudden'], r['E_a3_rel_spread'], r['nfev'], r['wall_s']))
    print(out['eps0_fit'])
