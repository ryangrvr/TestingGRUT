#!/usr/bin/env python3
"""cf_scale_skeptic: (1) independent P4b re-check (closed forms + own integration with betas retyped from the
Machacek-Vaughn check); (2) nonminimal-coupling relevance at the two Eichhorn-Pauly fixed points (their eq. 32,
quoted) -> whether GRUT's xi_X = -5/6 survives P4-M; (3) missed postulates: total asymptotic freedom (minimal SU(2)
embedding of hidden charges 1, 3, 6) and classical scale invariance (Coleman-Weinberg breaking of the charge-6 Sigma).
Usage: python3 -B p4b_and_missed.py <out.json>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
K = 16 * np.pi**2
MP, MY, MX, H = 2.435e18, 4.0e12, 0.0924, 5.33e12


def betas(y, gauge=True):
    lX, lY, lS, k, xy, xs, ys, g2, pY, pS, lp, e = y
    e = e if gauge else 0.0
    e2, e4 = e * e, e**4
    b = [6 * e4 - 12 * e2 * lX + 2 * g2**2 + 18 * k**2 + xs**2 + xy**2 + 20 * lX**2,
         486 * e4 - 108 * e2 * lY + xy**2 + ys**2 + 2 * pY**2 + 20 * lY**2,
         7776 * e4 - 432 * e2 * lS + xs**2 + ys**2 + 2 * pS**2 + 20 * lS**2,
         -36 * e2 * k + 6 * k * xy + 12 * k * lX,
         108 * e4 - 60 * e2 * xy + 4 * g2 * pY + 36 * k**2 + 2 * xs * ys + 4 * xy**2 + 8 * xy * (lX + lY),
         432 * e4 - 222 * e2 * xs + 4 * g2 * pS + 4 * xs**2 + 8 * xs * (lS + lX) + 2 * xy * ys,
         3888 * e4 - 270 * e2 * ys + 2 * xs * xy + 4 * ys**2 + 8 * ys * (lS + lY) + 4 * pS * pY,
         -6 * e2 * g2 + 8 * g2**2 + 8 * g2 * lX + 6 * g2 * lp + 2 * xs * pS + 2 * xy * pY,
         -54 * e2 * pY + 2 * g2 * xy + 2 * ys * pS + 8 * pY**2 + 8 * pY * lY + 6 * pY * lp,
         -216 * e2 * pS + 2 * g2 * xs + 2 * ys * pY + 8 * pS**2 + 8 * pS * lS + 6 * pS * lp,
         4 * g2**2 + 4 * pS**2 + 4 * pY**2 + 18 * lp**2,
         (46 / 3) * e**3 if gauge else 0.0]
    return np.array(b) / K


def reaches(y0, mu0, mu1, thr, gauge_above=None):
    ev = lambda t, y: np.max(np.abs(y[:11])) - thr
    ev.terminal = True
    rhs = lambda t, y: betas(y, gauge=(gauge_above is None or mu0 * np.exp(t) > gauge_above))
    s = solve_ivp(rhs, (0, np.log(mu1 / mu0)), y0, method='RK45', rtol=1e-9, atol=1e-14, events=ev)
    return (mu0 * np.exp(s.t_events[0][0]) if len(s.t_events[0]) else None), s.y[:, -1]


def vec(**kw):
    names = ['lX', 'lY', 'lS', 'k', 'xy', 'xs', 'ys', 'g2', 'pY', 'pS', 'lp', 'e']
    return np.array([kw.get(n, 0.0) for n in names], float)


def maxval(make, mu0, thr, lo, hi, gauge_above=None, strict1=False):
    def ok(x):
        pole, yend = reaches(make(x), mu0, MP, thr, gauge_above)
        return pole is None and (not strict1 or np.max(np.abs(yend[:11])) <= 1.0)
    return brentq(lambda x: 1.0 if ok(x) else -1.0, lo, hi, xtol=1e-5) if (ok(lo) and not ok(hi)) else None


def main(out):
    res = dict(method=__doc__)
    # ---------------- (1) P4b
    tP = np.log(MP / MY)
    analytic = dict(pole_lam1_GeV=MY * np.exp((1 - 1 / (4 * np.pi)) * K / 20),
                    common_max_4pi=1 / (20 * tP / K + 1 / (4 * np.pi)), common_max_strict1=1 / (20 * tP / K + 1.0),
                    common_max_unitarity_8pi_over_3=1 / (20 * tP / K + 3 / (8 * np.pi)),
                    kappa_floor_lamY=(0.0034586847 / ((1e-3 / (1.5 * 2**(-1 / 3)))**0.75))**4)
    spot = dict(lX=1e-3, lY=1.0, lS=1.0, k=0.0034586847, g2=0.01, lp=0.03, e=1e-7)
    pole, _ = reaches(vec(**spot), MY, MP, 4 * np.pi)
    num = dict(pole_spot_GeV=pole)
    num['common_max_4pi'] = maxval(lambda x: vec(**{**spot, 'lY': x, 'lS': x}), MY, 4 * np.pi, 0.1, 1.0)
    num['common_max_strict1'] = maxval(lambda x: vec(**{**spot, 'lY': x, 'lS': x}), MY, 4 * np.pi, 0.1, 1.0, strict1=True)
    c1 = dict(lX=1e-3, g2=0.01, lp=0.03, e=1e-7)
    num['C1_lamX_max'] = maxval(lambda x: vec(**{**c1, 'lX': x}), MX, 4 * np.pi, 1e-3, 1.0, gauge_above=H)
    num['C1_lamphi_max'] = maxval(lambda x: vec(**{**c1, 'lp': x}), MX, 4 * np.pi, 1e-3, 1.0, gauge_above=H)
    g2max = maxval(lambda x: vec(**{**c1, 'g2': x}), MX, 4 * np.pi, 1e-3, 2.0, gauge_above=H)
    num['C1_g_max'] = np.sqrt(g2max)
    num['C1_ms_over_mX_max_at_g0.1'] = np.sqrt(2 * num['C1_lamphi_max']) / 0.1
    num['e_max_MY_lam0.3'] = maxval(lambda x: vec(**{**spot, 'lY': 0.3, 'lS': 0.3, 'e': x}), MY, 4 * np.pi, 1e-3, 1.0)
    res['P4b'] = dict(analytic=analytic, numeric=num,
                      track=dict(pole=5.73e15, common_max_4pi=0.566, strict=0.372, floor=0.241, lamX=0.174, lamphi=0.193,
                                 g=0.570, ms_ratio=6.2, e_max=0.25))
    print('P4b analytic:', {k: float('%.4g' % v) for k, v in analytic.items()})
    print('P4b numeric :', {k: (None if v is None else float('%.4g' % v)) for k, v in num.items()}, flush=True)
    # ---------------- (2) xi relevance at the Eichhorn-Pauly fixed points (eq. 32 gravitational part, quoted)
    def grav_xi(g, l):
        return g * (2 * (3 - l) / (3 * np.pi * (3 - 4 * l)**2) - 5 * (3 + 10 * l) / (18 * np.pi * (1 - 2 * l)**3))
    def theta_lam(g, l):
        return -2 * (51 - 144 * l + 104 * l**2) * g / (3 * np.pi * (3 - 4 * l)**2 * (1 - 2 * l)**2)
    fps = {'pure_gravity': (1.32, 0.13), 'SM_content': (7.63, -9.97), 'U1D_deBrito': (2.78, -3.58)}
    xi = {}
    for k, (g, l) in fps.items():
        c = grav_xi(g, l)
        xi[k] = dict(beta_xi_grav_coeff=c, theta_xi=-c, xi_relevant=bool(-c > 0), A_lam=-theta_lam(g, l),
                     eq41_linear_xixi_coeff_times_g2=g**2 * (160 / (1 - 2 * l)**3 + 864 / (3 - 4 * l)**3))
    res['xi_relevance'] = dict(source='Eichhorn-Pauly arXiv:2009.13543 eq. (32) gravitational part and eq. (41) [Q]; '
                                      'eq. (41) overall normalisation as rendered by ar5iv, not independently verified',
                               discriminant_of_51_144_104=144**2 - 4 * 104 * 51, rows=xi)
    print('xi relevance:', {k: (round(v['theta_xi'], 4), v['xi_relevant'], round(v['A_lam'], 4)) for k, v in xi.items()}, flush=True)
    # ---------------- (3a) TAF: minimal SU(2) embedding with q = 2 T3 (charges 1,3,6 -> j >= 1/2, 3/2, 3)
    T = lambda j: j * (j + 1) * (2 * j + 1) / 3
    b0 = 22 / 3 - (T(0.5) + T(1.5) + T(3.0)) / 3 - T(1.0) / 6        # complex scalars + one real adjoint
    res['TAF_SU2'] = dict(T_half=T(0.5), T_3half=T(1.5), T_3=T(3.0), b0=b0, asymptotically_free=bool(b0 > 0),
                          note='16pi^2 beta_g = -b0 g^3; q = 2T3 is the minimal normalisation admitting charge 1')
    print('TAF minimal SU(2) embedding: b0 = %.3f (AF needs > 0)' % b0)
    # ---------------- (3b) classical scale invariance: CW breaking by the charge-6 Sigma (Landau gauge, MS-bar)
    # V_eff(phi) = lam phi^4/4 + (3 M_A^4/64pi^2)(ln M_A^2/mu^2 - 5/6), M_A = q e phi, |Sigma| = phi/sqrt2 (lam|S|^4 = lam phi^4/4)
    q = 6
    def cw_quant(e, w):
        B = 3 * (q * e)**4 / (64 * np.pi**2)
        # dV/dphi = 0 at phi = w, mu = w: lam w^3 + B(4 w^3 (ln(q^2 e^2) - 5/6) + 2 w^3) = 0
        # choose mu = q e w (log vanishes at the minimum):  lam(mu) = -4B(-5/6) - 2B = (4/3) B
        lam_muMA = (10 / 3) * B - 2 * B
        mS2 = 8 * B * w**2                       # second derivative at the minimum: V'' = 8 B w^2 (scheme-independent)
        return dict(lam_at_mu_eq_MA=lam_muMA, m_Sigma=np.sqrt(mS2), m_A=q * e * w, ratio=np.sqrt(mS2) / (q * e * w))
    res['scale_invariance_CW'] = {('e=%g' % e): cw_quant(e, H) for e in (1e-7, 3e-6, 1e-2, 0.3)}
    res['scale_invariance_CW']['note'] = ('w = h = 5.33e12 GeV taken as the generated scale; heavy spot needs lambda_Sigma '
                                         '= O(0.24-1) at h, i.e. m_Sigma ~ sqrt(2 lam) h; Higgs portal needed for the EW '
                                         'scale ~ (m_h/h)^2 ~ 5e-22')
    for k, v in res['scale_invariance_CW'].items():
        if k != 'note':
            print('CW %s: lam(mu=M_A) = %.3e  m_Sigma = %.3e GeV  m_A = %.3e GeV  m_Sigma/m_A = %.3e' % (
                k, v['lam_at_mu_eq_MA'], v['m_Sigma'], v['m_A'], v['ratio']))
    res['scale_invariance_CW']['portal_for_EW_scale'] = (125.13 / H)**2
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
