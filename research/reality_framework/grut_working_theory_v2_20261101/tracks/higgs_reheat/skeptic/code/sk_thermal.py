"""Skeptic of track higgs_reheat: thermal O(3) bounce re-computed with own thermal functions (Bessel-K series, not
quadrature) and own shooting, in two field normalisations.  Own code (Claude, skeptic pass, 2 Nov 2026).

'track': V/T^4 = lam_eff_flat(max(h,T)) u^4/4 + V_T(masses g^2 u^2/4 ...), u canonical  (lam_eff_flat contains
         e^{4 Gamma}: this is what thermal.py does).
'canon': the canonical field phi = e^{Gamma(T)} h0 carries both the quartic lam_eff (without e^{4 Gamma}) and the
         field-dependent masses, V/T^4 = lam_eff u^4/4 + V_T(g^2 u^2/4 ...), u = phi/T.
V_T = [6 J_B(m_W^2/T^2) + 3 J_B(m_Z^2/T^2) + 12 J_F(m_t^2/T^2)]/(2 pi^2) minus its h = 0 value,
J_B(y) = -sum_n y K_2(n sqrt y)/n^2, J_F(y) = sum_n (-1)^n y K_2(n sqrt y)/n^2 (same sign convention as thermal.py).
Writes ../results/sk_thermal.json.
"""
import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import kn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'code')))
import smpot as S  # noqa: E402

JB0, JF0 = -np.pi ** 4 / 45, -7 * np.pi ** 4 / 360


def Jseries(y, ferm, nmax=400):
    y = np.atleast_1d(np.asarray(y, float))
    out = np.empty_like(y)
    for i, yy in enumerate(y):
        if yy < 1e-6:
            out[i] = (JF0 + np.pi ** 2 * yy / 24) if ferm else (JB0 + np.pi ** 2 * yy / 12 - np.pi * yy ** 1.5 / 6)
            continue
        s = np.sqrt(yy)
        n = np.arange(1, max(nmax, int(40 / s)) + 1)
        terms = yy * kn(2, n * s) / n ** 2
        out[i] = np.sum(((-1.0) ** n) * terms) if ferm else -np.sum(terms)
    return out


class TV:
    def __init__(self, sm, T, norm):
        c = sm.cpl(np.pi * T)
        g2, gp2, yt2 = c['G2'], 0.6 * c['G1'], c['T']
        lh = np.linspace(np.log(1e10), np.log(1e19), 900)
        le = sm.lam_eff_flat(np.exp(lh))
        if norm == 'canon':
            le = le * np.exp(-4 * sm.cpl(T)['Gam'])
        lam = CubicSpline(lh, le)
        us = np.concatenate([[0.0], np.exp(np.linspace(np.log(1e-4), np.log(60.0), 3000))])
        hs = np.maximum(us * T, T)
        V0 = lam(np.log(hs)) * us ** 4 / 4
        yW, yZ, yt = g2 * us ** 2 / 4, (g2 + gp2) * us ** 2 / 4, yt2 * us ** 2 / 2
        VT = (6 * (Jseries(yW, False) - JB0) + 3 * (Jseries(yZ, False) - JB0) + 12 * (Jseries(yt, True) - JF0)) / (2 * np.pi ** 2)
        self.us, self.Vs = us, V0 + VT
        self.cs = CubicSpline(us, self.Vs)
        self.d1 = self.cs.derivative()
        self.kT = 3 * g2 / 16 + gp2 / 16 + yt2 / 4


def bounce(tv):
    us, Vs = tv.us, tv.Vs
    k = int(np.argmax(Vs[us < 20]))
    ub = us[k]
    neg = np.where((us > ub) & (Vs < 0))[0]
    if not len(neg):
        return dict(S3_over_T=np.inf, barrier_u=float(ub))
    ue = us[neg[0]]
    dV = lambda u: float(tv.d1(u))

    def shoot(u0, full=False):
        # y = (u, du/dr); start from small r with the regular expansion
        r0 = 1e-5
        a0 = dV(u0) / 3

        def rhs(r, y):
            return [y[1], dV(min(max(y[0], 0.0), 59.0)) - 2 * y[1] / r]

        def ev_over(r, y):
            return y[0]
        ev_over.terminal = True

        def ev_under(r, y):
            return y[1]
        ev_under.terminal = True
        ev_under.direction = 1
        sol = solve_ivp(rhs, (r0, 400.0), [u0 + 0.5 * a0 * r0 ** 2, a0 * r0], method='LSODA', rtol=1e-11,
                        atol=1e-13, events=[ev_over, ev_under], dense_output=full)
        return (1 if len(sol.t_events[0]) else -1), sol
    lo, hi = ue, 1.5 * ue
    while shoot(hi)[0] < 0:
        hi *= 1.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if shoot(mid)[0] > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-12 * hi:
            break
    _, sol = shoot(lo, full=True)
    # use the trajectory up to the turning point (undershoot): the bounce tail is within accuracy there
    r = np.linspace(sol.t[0], sol.t[-1], 20000)
    u, up = sol.sol(r)
    Vr = tv.cs(np.clip(u, 0, 60))
    K = 4 * np.pi * np.trapezoid(r ** 2 * 0.5 * up ** 2, r)
    P = 4 * np.pi * np.trapezoid(r ** 2 * Vr, r)
    return dict(S3_over_T=float(K + P), derrick_K_over_minus3P=float(K / (-3 * P)), u0=float(lo), barrier_u=float(ub),
                kappa_T=float(tv.kT), r_end=float(sol.t[-1]))


def main():
    out = {}
    for Mt in (172.87, 172.60):
        sm = S.SM(Mt)
        rows = {}
        for T in (1e13, 1.6e14, 3e14):
            r = {}
            for norm in ('track', 'canon'):
                r[norm] = bounce(TV(sm, T, norm))
            r['e4Gamma(T)'] = float(np.exp(4 * sm.cpl(T)['Gam']))
            r['ratio_canon_over_track'] = r['canon']['S3_over_T'] / r['track']['S3_over_T']
            rows[f'{T:.3g}'] = r
            print(Mt, T, json.dumps(r, default=float), flush=True)
        out[f'Mt={Mt}'] = rows
    # J check against quadrature
    from scipy.integrate import quad
    chk = {}
    for y in (0.01, 1.0, 10.0):
        qb = quad(lambda x: x * x * np.log1p(-np.exp(-np.sqrt(x * x + y))), 0, 60)[0]
        qf = quad(lambda x: -x * x * np.log1p(np.exp(-np.sqrt(x * x + y))), 0, 60)[0]
        chk[f'y={y}'] = dict(JB_series=float(Jseries(y, False)[0]), JB_quad=qb, JF_series=float(Jseries(y, True)[0]),
                             JF_quad=qf)
    out['J_check'] = chk
    print(chk)
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'sk_thermal.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
