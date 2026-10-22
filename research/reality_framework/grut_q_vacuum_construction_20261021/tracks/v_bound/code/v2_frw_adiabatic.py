"""v_bound: independent full-FRW check of bound A.2/A.3 (regime (i) relation) with a different formulation.

Written from scratch (nothing imported from ../../bound/code).  Differences from the bound track's integrator:
  * independent variable N = ln a (not t); state (q, p = dq/dN) only;
  * H^2 is NOT evolved: it is solved from the Friedmann constraint at every step,
        H^2 = (rho_d + rho_r + U) / (3 (K + K' p)),   U = eps - mu q
    (equations from v1_symbolic.py V1a, derived from the minisuperspace action), so the constraint holds exactly;
  * the q-equation K' R = 2 U' is monitored instead (R_geom = 3 dH^2/dN + 12 H^2 vs R_q = 2 U'/K');
  * slow-manifold initial data by fixed-point iteration of the full equation (Chebyshev-smoothed derivatives),
    not by a second-order formula;
  * oscillation removed by a Hann window of two periods (same window for the run and for the formula).
Trace equation used for p' (V1a):  Kddot + 3 H Kdot = (rho_d + 4 U - K R_q)/3  =>
        p' = [ S/H^2 - K'' p^2 - K' p (R_q/(6 H^2) + 1) ] / K',   S = (rho_d + 4 U - K R_q)/3.
Quantity compared (bound A.2):  Delta = 3 K_lab H^2 - rho_d - rho_r - Lam_cal  vs
        P = [ (9/4) rho_d^2 + 2 rho_d rho_r - 4 Lam_cal rho_r ] / (3 K_lab m0^2),
K_lab = K at the R = 0 point of the calibrated mu, m0^2 = (2/3) K eps''/K'^2 there, Lam_cal = 3 K_lab H_s^2 (static root),
mu calibrated so that Lam_cal = 3 K_lab * 0.685 (H0 = 1 units), rho_d0 = 3 Om K_lab, rho_r0 = 3 Or K_lab.
Families: F1-F3 as described in bound NOTES A.3 (re-implemented from the NOTES text), plus F4 (new, not in bound).
Run: python3 v2_frw_adiabatic.py   (writes ../results/v2_frw_adiabatic.json)
"""
import json
import os
import time
import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OM, OR = 0.315, 9.1e-5
OL = 1.0 - OM - OR


def family(name, m):
    """returns dict of callables eps, e1, e2, K, K1, K2 and mu0 (Minkowski point q0 = 1, K(1) = 1, m(q0) = m)."""
    if name == 'F1':      # eps = (e/2)(1 + q^2), K = 1 + 0.4 (q - 1): f exactly quadratic
        k = 0.4
        e = 1.5 * m**2 * k**2
        return dict(eps=lambda q: 0.5 * e * (1 + q * q), e1=lambda q: e * q, e2=lambda q: e + 0 * q,
                    K=lambda q: 1 + k * (q - 1), K1=lambda q: k + 0 * q, K2=lambda q: 0 * q, mu0=e)
    if name == 'F2':      # KV 0806.2805: eps = (E/2)(-q^2 + q^4/3), K = q
        E = 1.5 * m**2
        return dict(eps=lambda q: 0.5 * E * (-q * q + q**4 / 3), e1=lambda q: E * (-q + 2 * q**3 / 3),
                    e2=lambda q: E * (-1 + 2 * q * q), K=lambda q: q, K1=lambda q: 1 + 0 * q, K2=lambda q: 0 * q, mu0=-E / 3)
    if name == 'F3':      # eps = e x^2/2 + 3e x^3/6 + 10e x^4/24 + mu0 q, K = exp(0.7 x)
        b = 0.7
        e = 1.5 * m**2 * b**2
        mu0 = 0.3 * e
        return dict(eps=lambda q: e * (q - 1)**2 / 2 + 3 * e * (q - 1)**3 / 6 + 10 * e * (q - 1)**4 / 24 + mu0 * q,
                    e1=lambda q: e * (q - 1) + 3 * e * (q - 1)**2 / 2 + 10 * e * (q - 1)**3 / 6 + mu0,
                    e2=lambda q: e + 3 * e * (q - 1) + 10 * e * (q - 1)**2 / 2,
                    K=lambda q: np.exp(b * (q - 1)), K1=lambda q: b * np.exp(b * (q - 1)), K2=lambda q: b * b * np.exp(b * (q - 1)),
                    mu0=mu0)
    if name == 'F4':      # new: eps = mu0 q + e [cosh(x) - 1] (x = q - 1), K = 1/(1 - 0.5 x)  (K'' != 0, eps''' = 0 at q0)
        g = 0.5
        e = 1.5 * m**2 * g**2
        mu0 = 2.0 * e
        return dict(eps=lambda q: mu0 * q + e * (np.cosh(q - 1) - 1), e1=lambda q: mu0 + e * np.sinh(q - 1),
                    e2=lambda q: e * np.cosh(q - 1), K=lambda q: 1 / (1 - g * (q - 1)), K1=lambda q: g / (1 - g * (q - 1))**2,
                    K2=lambda q: 2 * g * g / (1 - g * (q - 1))**3, mu0=mu0)
    raise ValueError(name)


def root_near(fun, x0, w0=1e-3, wmax=0.9):
    w = w0
    while w < wmax:
        lo, hi = x0 - w, x0 + w
        if np.sign(fun(lo)) != np.sign(fun(hi)):
            return brentq(fun, lo, hi, xtol=1e-15, rtol=1e-15)
        w *= 1.3
    raise RuntimeError('no root near %g' % x0)


class Model:
    def __init__(self, fam, mu, rd0, rr0):
        self.f, self.mu, self.rd0, self.rr0 = fam, mu, rd0, rr0

    def U(self, q):
        return self.f['eps'](q) - self.mu * q

    def Rq(self, q):
        return 2 * (self.f['e1'](q) - self.mu) / self.f['K1'](q)

    def H2(self, N, q, p):
        rd, rr = self.rd0 * np.exp(-3 * N), self.rr0 * np.exp(-4 * N)
        return (rd + rr + self.U(q)) / (3 * (self.f['K'](q) + self.f['K1'](q) * p))

    def pprime(self, N, q, p):
        f = self.f
        rd = self.rd0 * np.exp(-3 * N)
        H2 = self.H2(N, q, p)
        R = self.Rq(q)
        S = (rd + 4 * self.U(q) - f['K'](q) * R) / 3
        return (S / H2 - f['K2'](q) * p * p - f['K1'](q) * p * (R / (6 * H2) + 1)) / f['K1'](q)

    def rhs(self, N, y):
        q, p = y
        return [p, self.pprime(N, q, p)]

    def R_geom(self, N, q, p):
        """3 dH^2/dN + 12 H^2 with dH^2/dN by the chain rule (central differences in each argument)."""
        h = 1e-7
        dN = (self.H2(N + h, q, p) - self.H2(N - h, q, p)) / (2 * h)
        dq = (self.H2(N, q + h, p) - self.H2(N, q - h, p)) / (2 * h)
        dp = (self.H2(N, q, p + h) - self.H2(N, q, p - h)) / (2 * h)
        dH2 = dN + dq * p + dp * self.pprime(N, q, p)
        return 3 * dH2 + 12 * self.H2(N, q, p)


def calibrate(fam):
    """mu such that Lam_cal = 3 K_lab H_s^2 = 3 K_lab OL; K_lab at the R = 0 point (eps' = mu)."""
    def pieces(mu):
        qlab = root_near(lambda q: fam['e1'](q) - mu, 1.0)
        stat = lambda q: (fam['e1'](q) - mu) * fam['K'](q) - 2 * fam['K1'](q) * (fam['eps'](q) - mu * q)
        qs = root_near(stat, qlab, w0=1e-6)
        H2s = (fam['eps'](qs) - mu * qs) / (3 * fam['K'](qs))
        return qlab, qs, H2s

    def g(dmu):
        qlab, qs, H2s = pieces(fam['mu0'] + dmu)
        return H2s - OL
    # dmu < 0 gives dS (rho_V = -q0 dmu > 0)
    lo, hi = -1e-6, -0.2
    while np.sign(g(hi)) == np.sign(g(lo)):
        hi *= 1.4
    dmu = brentq(g, hi, lo, xtol=1e-15, rtol=1e-15)
    qlab, qs, H2s = pieces(fam['mu0'] + dmu)
    Klab = fam['K'](qlab)
    m0sq = (2.0 / 3.0) * Klab * fam['e2'](qlab) / fam['K1'](qlab)**2
    return dict(dmu=dmu, mu=fam['mu0'] + dmu, qlab=qlab, qs=qs, Klab=Klab, Lam=3 * Klab * H2s, m0=np.sqrt(m0sq), Hs=np.sqrt(H2s))


def slow_manifold(model, Nc, half=0.6, npts=1201, iters=6, deg=24):
    """fixed-point iteration q^(k+1)(N) solving the p'-equation for q with p, p' from q^(k) (Chebyshev-smoothed)."""
    Ng = np.linspace(Nc - half, Nc + half, npts)
    f = model.f
    # zeroth iterate: quasi-static, S = 0 at p = 0
    q = np.empty_like(Ng)
    guess = 1.0
    for i, N in enumerate(Ng):
        rd = model.rd0 * np.exp(-3 * N)
        q[i] = root_near(lambda x: rd + 4 * model.U(x) - f['K'](x) * model.Rq(x), guess, w0=1e-7)
        guess = q[i]
    hist = []
    for k in range(iters):
        cf = C.chebfit((Ng - Nc) / half, q, deg)
        p = C.chebval((Ng - Nc) / half, C.chebder(cf)) / half
        pp = C.chebval((Ng - Nc) / half, C.chebder(cf, 2)) / half**2
        qn = np.empty_like(q)
        for i, N in enumerate(Ng):
            def resid(x):
                return model.pprime(N, x, p[i]) - pp[i]
            qn[i] = root_near(resid, q[i], w0=1e-9)
        hist.append(float(np.max(np.abs(qn - q))))
        q = qn
    cf = C.chebfit((Ng - Nc) / half, q, deg)
    qc = float(C.chebval(0.0, cf))
    pc = float(C.chebval(0.0, C.chebder(cf)) / half)
    return qc, pc, hist


def run(name, m, a_evals=(1.0,), start_ratio=1.25, iters=6):
    fam = family(name, m)
    cal = calibrate(fam)
    rd0, rr0 = 3 * OM * cal['Klab'], 3 * OR * cal['Klab']
    model = Model(fam, cal['mu'], rd0, rr0)
    m0 = cal['m0']
    res = {'family': name, 'm_target': m, 'm0': m0, 'dmu': cal['dmu'], 'K_lab': cal['Klab'], 'Lam_cal': cal['Lam'],
           'H_static': cal['Hs'], 'evals': []}
    for ae in a_evals:
        Ne = np.log(ae)
        Ns = Ne - np.log(start_ratio)
        q0_, p0_, hist = slow_manifold(model, Ns, iters=iters)
        H_e_est = np.sqrt(model.H2(Ne, 1.0, 0.0))
        dNper = 2 * np.pi * H_e_est / m0
        Nend = Ne + 1.5 * dNper
        sol = solve_ivp(model.rhs, (Ns, Nend), [q0_, p0_], method='DOP853', rtol=1e-12, atol=1e-14,
                        max_step=dNper / 60, dense_output=True)
        Nw = np.linspace(Ne - dNper, Ne + dNper, 2001)
        w = np.sin(np.pi * (Nw - Nw[0]) / (Nw[-1] - Nw[0]))**2          # Hann window over two periods
        q, p = sol.sol(Nw)
        H2 = model.H2(Nw, q, p)
        rd, rr = rd0 * np.exp(-3 * Nw), rr0 * np.exp(-4 * Nw)
        Klab, Lam = cal['Klab'], cal['Lam']
        D = 3 * Klab * H2 - rd - rr - Lam
        P = (2.25 * rd**2 + 2 * rd * rr - 4 * Lam * rr) / (3 * Klab * m0**2)
        ratio = np.trapezoid(w * D, Nw) / np.trapezoid(w * P, Nw)
        Hm = np.sqrt(np.trapezoid(w * H2, Nw) / np.trapezoid(w, Nw))
        # q-equation (not imposed) check along the run
        Nck = np.linspace(Ns, Nend, 300)
        qq, pp = sol.sol(Nck)
        Rg = np.array([model.R_geom(N, a, b) for N, a, b in zip(Nck, qq, pp)])
        Rqv = model.Rq(qq)
        qeq = float(np.max(np.abs(Rg - Rqv) / np.abs(Rqv)))
        osc = float(np.std(D - P) / np.mean(P))
        res['evals'].append({'a': ae, 'a_start': float(np.exp(Ns)), 'H/m0 at a': float(Hm / m0),
                             'ratio <Delta>/<P> (Hann, 2 periods)': float(ratio),
                             'c = (1 - ratio)/(H/m0)^2': float((1 - ratio) / (Hm / m0)**2),
                             'Omega_dyn measured': float(np.trapezoid(w * D, Nw) / np.trapezoid(w * 3 * Klab * H2, Nw)),
                             'Omega_dyn formula (9/4) Om_run^2 (H/m0)^2': float(2.25 * (rd0 / (3 * Klab * Hm**2))**2 * (Hm / m0)**2),
                             'slow-manifold iteration max|dq| per iteration': hist,
                             'residual oscillation rms(Delta - P)/<P>': osc,
                             'max relative q-equation violation (not imposed)': qeq})
    return res


def main():
    t0 = time.time()
    out = {'inputs': {'Omega_m': OM, 'Omega_r': OR, 'Omega_L': OL, 'units': 'H0 = 1 (calibrated), K(q0) = 1'}}
    runs = []
    for name in ['F1', 'F2', 'F3', 'F4']:
        for m in [30.0, 100.0, 300.0]:
            r = run(name, m)
            runs.append(r)
            e = r['evals'][0]
            print(name, m, round(e['H/m0 at a'], 6), e['ratio <Delta>/<P> (Hann, 2 periods)'], e['c = (1 - ratio)/(H/m0)^2'],
                  e['Omega_dyn measured'], e['residual oscillation rms(Delta - P)/<P>'], e['max relative q-equation violation (not imposed)'],
                  flush=True)
    out['runs'] = runs
    # m = 10 (outside the regime, for the 'fails below ~10 H0' statement)
    r10 = []
    for name in ['F1', 'F2']:
        try:
            r = run(name, 10.0, iters=1)          # one fixed-point step only (the iteration diverges at H/m ~ 0.1)
            r['note'] = 'initial data: quasi-static + one fixed-point step; residual free oscillation larger'
            r10.append(r)
            e = r['evals'][0]
            print(name, 10, e['ratio <Delta>/<P> (Hann, 2 periods)'], e['c = (1 - ratio)/(H/m0)^2'], flush=True)
        except Exception as ex:          # slow-manifold iteration may not converge at H/m ~ 0.1
            r10.append({'family': name, 'm_target': 10.0, 'error': repr(ex)})
    out['runs_m10'] = r10
    out['runtime_s'] = time.time() - t0
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'v2_frw_adiabatic.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
