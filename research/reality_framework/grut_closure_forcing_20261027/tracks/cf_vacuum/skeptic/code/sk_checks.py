"""Skeptic checks of the cf_vacuum track's hidden-sector, threshold and portal claims, by methods different from the
track's (direct multi-field minimisation, symbolic Hessians, own RG extension).  Nothing imported from ../../code.

Run: python3 -B sk_checks.py  (~1-2 min).  Writes ../results/sk_checks.json.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize
from scipy.integrate import solve_ivp

import sk_rge as R

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'sk_checks.json')
res = {}

# ------------------------------------------------------------------ GRUT heavy reference (04c53f5 THEORY_DATA.json)
H_I = 133349619398.33566
M_Y = 4000488581950.07
h = 5333984775933.427
mu = 954710753970.2344
lamY = 1.0
s_vac = h / np.sqrt(2)

# (1) hidden-sector P3a bound by direct 3-variable minimisation (|Sigma|, |Y|, relative phase), not via m^4/(4 lam_Y)
def V_hidden(s, y, ph, lamS, lamSY, lamY=1.0):
    # V = lamS(s^2 - h^2/2)^2 + M_Y^2 y^2 + lamY y^4 + lamSY s^2 y^2 + 2 mu s y^2 cos(ph)   (|Y| = y, |Sigma| = s)
    return lamS * (s ** 2 - s_vac ** 2) ** 2 + M_Y ** 2 * y ** 2 + lamY * y ** 4 + lamSY * s ** 2 * y ** 2 \
        + 2 * mu * s * y ** 2 * np.cos(ph)


S_GRID = np.concatenate([np.linspace(0.0, 10.0, 2001), np.logspace(1, 3.5, 1500)])   # s/s_vac
Y_GRID = np.concatenate([[0.0], np.logspace(-4, 3, 1400)])                             # y/s_vac


def global_min(lamS, lamSY, lamY=1.0):
    """Brute-force 2-D grid in (s, y) at the phase cos = -1 (which minimises the only phase-dependent term,
    2 mu s y^2 cos(ph), for s, y > 0 [D]), then Nelder-Mead refinement in all three variables."""
    S, Y = np.meshgrid(S_GRID * s_vac, Y_GRID * s_vac, indexing='ij')
    Vg = V_hidden(S, Y, np.pi, lamS, lamSY, lamY) / s_vac ** 4
    i, j = np.unravel_index(np.argmin(Vg), Vg.shape)
    f = lambda p: V_hidden(p[0] * s_vac, p[1] * s_vac, p[2], lamS, lamSY, lamY) / s_vac ** 4
    r = minimize(f, [S_GRID[i], Y_GRID[j], np.pi], method='Nelder-Mead',
                 options=dict(xatol=1e-12, fatol=1e-18, maxiter=20000))
    return (min(r.fun, Vg[i, j]), r.x)


def lamS_min(lamSY, lamY=1.0):
    g = lambda lS: global_min(lS, lamSY, lamY)[0]
    return brentq(lambda x: g(10 ** x) + 1e-14, -9, 0, xtol=1e-6)


hb = {}
hb['r_at_vacuum'] = 2 * mu * s_vac / M_Y ** 2
hb['tachyon_onset_s_over_svac'] = M_Y ** 2 / (2 * mu * s_vac)
hb['lamS_min_lamSY0_lamY1'] = 10 ** lamS_min(0.0)
hb['lamS_min_lamSY0.03_lamY1'] = 10 ** lamS_min(0.03)
hb['lamSY_threshold_mu2_over_MY2'] = (mu / M_Y) ** 2
gm = global_min(2.8e-5, 0.0)
hb['window_point_lamS2.8e-5_lamSY0: Vmin/s_vac^4, (s/svac, y/svac, phase)'] = [gm[0], list(gm[1])]
hb['window_point_lamS2.8e-5_lamSY0: depth GeV^4'] = gm[0] * s_vac ** 4
RHO_L = 2.5e-47  # GeV^4, rho_Lambda (Planck 2018, [R])
hb['window_point depth / rho_Lambda'] = -gm[0] * s_vac ** 4 / RHO_L
for lamS in (2.8e-5, 1.8e-10):
    f = lambda lSY: global_min(lamS, lSY)[0] + 1e-14
    hb[f'lamSY_needed_at_lamS={lamS}'] = brentq(f, 0.0, 0.0570)
res['hidden_sector_P3a'] = hb

# (2) one-loop Y threshold on lambda_Sigma and m_Sigma^2: symbolic CW expansion
s_, M2, m_, muR = sp.symbols('s M2 mu muR', positive=True)
f = lambda m2: m2 ** 2 * (sp.log(m2 / muR ** 2) - sp.Rational(3, 2)) / (64 * sp.pi ** 2)
Vcw = f(M2 + 2 * m_ * s_) + f(M2 - 2 * m_ * s_)
ser = sp.series(Vcw, s_, 0, 5).removeO()
c4 = sp.simplify(ser.coeff(s_, 4))
c2 = sp.simplify(ser.coeff(s_, 2))
res['CW_Y_threshold'] = {
    'coeff_s4 (= delta lambda_Sigma, |Sigma|=s)': str(c4),
    'numeric': float(c4.subs({m_: mu, M2: M_Y ** 2})),
    'track value': -1.37e-5,
    'coeff_s2 (= delta m_Sigma^2)': str(sp.simplify(c2)),
}

# (3) mixing identity: theta * M_S = sqrt(2 delta) v, exact 2x2 diagonalisation (EEGLS normalisation)
v = 246.22
rows = []
rng = np.random.default_rng(1)
for _ in range(6):
    lS = 10 ** rng.uniform(-14, -4)
    w = 10 ** rng.uniform(8, 13)
    delta = rng.uniform(0.007, 0.36)
    lHS = np.sqrt(delta * lS)
    lH = 0.13 + delta
    Mm = np.array([[2 * lH * v ** 2, 2 * lHS * v * w], [2 * lHS * v * w, 2 * lS * w ** 2]])
    ev, U = np.linalg.eigh(Mm)
    theta = abs(U[0, 1])
    MS = np.sqrt(ev[1])
    rows.append({'lam_S': lS, 'w': w, 'delta': delta, 'theta*M_S': theta * MS,
                 'sqrt(2 delta) v': np.sqrt(2 * delta) * v})
res['mixing_identity'] = rows
res['band_sqrt2delta_v'] = [np.sqrt(2 * 0.007) * v, np.sqrt(2 * 0.36) * v]

# (4) valley running identity (EEGLS one-loop scalar RGEs, SM pieces generic)
lH, lS, lHS, gH, Rr = sp.symbols('lH lS lHS gammaH R')
bH = 24 * lH ** 2 + 4 * lHS ** 2 + 4 * gH * lH + Rr        # 16pi^2 d lH / d ln mu
bS = 8 * lHS ** 2 + 20 * lS ** 2
bHS = 4 * lHS * (3 * lH + 2 * lHS + 2 * lS) + 2 * gH * lHS
delta = lHS ** 2 / lS
bdelta = sp.diff(delta, lHS) * bHS + sp.diff(delta, lS) * bS
lv = sp.Symbol('lv')
lhs = sp.expand((bH - bdelta).subs(lH, lv + delta))
rhs = sp.expand(24 * lv ** 2 + 4 * gH * lv + Rr + delta * (24 * lv + 8 * delta + 8 * lS - 16 * lHS))
res['valley_identity_holds'] = bool(sp.simplify(lhs - rhs) == 0)

# (5) one-loop scalar beta-function pieces for the Higgs portals, by 1/2 Tr M^4 (Coleman-Weinberg log coefficient)
hs = sp.symbols('h1:5', real=True)
ph, x1, x2 = sp.symbols('phi x1 x2', real=True)
lam, b, lphi, lHX, lX = sp.symbols('lam b lphi lHX lX')
HH = sum(t ** 2 for t in hs) / 2
XX = (x1 ** 2 + x2 ** 2) / 2
Vs = lam * HH ** 2 + b * ph ** 2 * HH + lphi / 4 * ph ** 4 + lHX * HH * XX + lX * XX ** 2
flds = list(hs) + [ph, x1, x2]
Hs = sp.hessian(Vs, flds)
trM4 = sp.expand((Hs * Hs).trace() / 2)
# read off coefficients of the |H|^4, phi^2|H|^2, |X|^2|H|^2 monomials at a point
pt = {hs[0]: sp.Symbol('a'), hs[1]: 0, hs[2]: 0, hs[3]: 0, ph: sp.Symbol('p'), x1: sp.Symbol('q'), x2: 0}
T = sp.expand(trM4.subs(pt))
a_, p_, q_ = sp.Symbol('a'), sp.Symbol('p'), sp.Symbol('q')
coef_H4 = sp.simplify(T.coeff(a_, 4).subs({p_: 0, q_: 0}) * 4)          # |H|^4 = a^4/4
coef_bp = sp.simplify(T.coeff(a_, 2).coeff(p_, 2) * 2)                  # phi^2|H|^2 = p^2 a^2/2
coef_hx = sp.simplify(T.coeff(a_, 2).coeff(q_, 2) * 4)                  # |H|^2|X|^2 = a^2 q^2/4
res['portal_betas_scalar_part'] = {
    '16pi^2 dlam/dlnmu (scalar part)': str(coef_H4),
    '16pi^2 db/dlnmu (scalar part)': str(coef_bp),
    '16pi^2 dlamHX/dlnmu (scalar part)': str(coef_hx)}

# (6) loop stabilisation by phi or X: SM 3-loop (Bednyakov matching) + 1-loop portal pieces, PDG 2026 central
MH, AS, MT = 125.13, 0.1180, 172.60
b_lam_extra = sp.lambdify((lam, b, lphi, lHX, lX), coef_H4 - (24 * lam ** 2), 'numpy')
bb = sp.lambdify((lam, b, lphi, lHX, lX), coef_bp, 'numpy')
bx = sp.lambdify((lam, b, lphi, lHX, lX), coef_hx, 'numpy')


def run_portal(b0, lhx0, lphi0=0.03, lX0=1e-3):
    mu0, c = R.match_B15(MT, MH, AS)
    y0 = [5 / 3 * c['gY'] ** 2, c['g2'] ** 2, c['g3'] ** 2, c['yt'] ** 2, c['lam'], c['yb'] ** 2, c['ytau'] ** 2,
          b0, lhx0, lphi0, lX0]

    def rhs(t, y):
        sm = R.betas(t, y[:7], 3)
        G1, G2, G3, YT, L = y[:5]
        B, LHX, LP, LX = y[7:]
        gam = 6 * YT - 9 * G2 / 2 - 9 * G1 / 10   # 2 gamma_H for one |H|^2 (per ln mu), x 16pi^2
        # d/d ln mu^2 = 1/2 d/d ln mu
        dL_extra = 0.5 * b_lam_extra(L, B, LP, LHX, LX) / R.K1
        dB = 0.5 * (bb(L, B, LP, LHX, LX) + gam * B) / R.K1
        dLHX = 0.5 * (bx(L, B, LP, LHX, LX) + gam * LHX) / R.K1
        dLP = 0.5 * (18 * LP ** 2 + 8 * B ** 2) / R.K1     # 1/2 Tr M^4 for lphi/4 phi^4 + b phi^2|H|^2 [D below]
        dLX = 0.5 * (20 * LX ** 2 + 2 * LHX ** 2) / R.K1
        sm[4] += dL_extra
        return sm + [dB, dLHX, dLP, dLX]
    return solve_ivp(rhs, (np.log(mu0 ** 2), np.log(R.MPL ** 2)), y0, method='DOP853', rtol=1e-9, atol=1e-12,
                     dense_output=True)


def minlam(sol):
    xs = np.linspace(np.log(1e3 ** 2), np.log(R.MPL ** 2), 800)
    return float(np.min(sol.sol(xs)[4]))


port = {}
port['b_needed'] = brentq(lambda b0: minlam(run_portal(b0, 0.0)), 0.01, 0.5, xtol=1e-5)
port['lamHX_needed'] = brentq(lambda l0: minlam(run_portal(0.0, l0)), 0.01, 0.5, xtol=1e-5)
port['track'] = {'b': 0.117, 'lamHX': 0.175}
# invisible width bounds
GSM = 4.07e-3       # GeV, SM Higgs width at ~125 GeV (LHC HXSWG, [R])
Binv = 0.107        # ATLAS 2301.10731 (95% CL)
Ginv = Binv / (1 - Binv) * GSM
port['Gamma_inv_max_GeV'] = Ginv
port['b_max'] = np.sqrt(Ginv * 8 * np.pi * MH / v ** 2)
port['lamHX_max'] = np.sqrt(Ginv * 16 * np.pi * MH / v ** 2)
port['Gamma(h->phiphi) at b_needed'] = port['b_needed'] ** 2 * v ** 2 / (8 * np.pi * MH)
port['ratio b_needed/b_max'] = port['b_needed'] / port['b_max']
port['ratio lamHX_needed/lamHX_max'] = port['lamHX_needed'] / port['lamHX_max']
res['portal_loop_stabilisation'] = port

# (7) GRUT reference hidden quartics: one-loop Landau pole (16pi^2 dlam/dlnmu = 20 lam^2, lam|Z|^4 normalisation)
lp = {}
for name, l0, m0 in (('lambda_Sigma=1 at h', 1.0, h), ('lambda_Y=1 at M_Y', 1.0, M_Y)):
    tpole = R.K1 / (20 * l0)
    t4pi = R.K1 / 20 * (1 / l0 - 1 / (4 * np.pi))
    lp[name] = {'pole_GeV': m0 * np.exp(tpole), 'lambda=4pi_GeV': m0 * np.exp(t4pi)}
lp['check 20 lam^2'] = 'O(2) model, V = lam (phi.phi/2)^2: 16pi^2 beta = 20 lam^2 [D: 1/2 Tr M^4 above with x1,x2]'
res['hidden_landau_pole'] = lp
res['xX_check_lX_beta_coefficient'] = str(sp.simplify(T.coeff(q_, 4).subs({a_: 0, p_: 0}) * 4))

# (8) inflation-scale numbers
hmax = 2.7e11
res['inflation'] = {'H_i/h_max': H_I / hmax, 'h_rms_60efolds': np.sqrt(60) * H_I / (2 * np.pi),
                    'P(|h|>h_max) per Hubble patch (Gaussian)': float(__import__('scipy.special', fromlist=['erfc'])
                                                                       .erfc(hmax / (np.sqrt(2) * np.sqrt(60) * H_I / (2 * np.pi)))),
                    'patches ~ e^(3N), N=60': 'e^180 = 1.5e78',
                    'xi_H window (quoted, abstracts)': 'lower ~0.05-0.07 (2011.03763); upper ~1 from post-inflation '
                                                      'curvature-induced Higgs production (1506.04065; 1602.00483)'}

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
