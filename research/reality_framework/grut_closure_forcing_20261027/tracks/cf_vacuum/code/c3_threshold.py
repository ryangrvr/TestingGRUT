"""Can new scalars (generic, or GRUT's Sigma, Y, phi, X) make the electroweak vacuum absolutely stable?

Track cf_vacuum (closure-forcing pass, 27 Oct 2026).  Run: python3 -B c3_threshold.py  (~3 min).
Writes ../results/threshold.json.

Physics (credited): tree-level threshold of a heavy singlet with a vev, Elias-Miro, Espinosa, Giudice, Lee, Strumia,
arXiv:1203.0237 (EEGLS, TeX read: eqs. pot, RG, match, condp, condm); Lebedev 1203.0156.  EEGLS normalisation:
    V = lam_H (H^+H - v^2/2)^2 + lam_S (S^+S - w^2/2)^2 + 2 lam_HS (H^+H - v^2/2)(S^+S - w^2/2),
    M_S^2 = 2 lam_S w^2,  delta = lam_HS^2 / lam_S,  lam_SM(M_S) = lam_H(M_S) - delta.
GRUT normalisation (THEORY.md / ACTION.md): lam_Sigma (|Sigma|^2 - h^2/2)^2 + lam_HSigma |H|^2 |Sigma|^2, so
    lam_Sigma = lam_S, h = w, lam_HSigma = 2 lam_HS.

Precision: SM couplings with 3-loop RGEs (sm_rge.py, from SMDR-matched values at M_t, PDG 2026 central point), singlet
couplings and the singlet contribution to beta(lam_H) at one loop, tree-level matching at mu = M_S (no one-loop
threshold pieces), RG-improved tree potential with mu = h, and the MSbar criterion lam >= 0 used for the SM line.
Expected accuracy: delta_min to ~0.002 absolute, scale boundaries to a factor ~2 (criterion choice: the lambda_eff
crossing is ~7x higher than the MSbar one, see sm_relation.json).
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import sm_rge as R

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
K16 = 16 * np.pi ** 2
MPL = 1.2209e19
VEW = 246.22
MT, MH0, AS0 = 172.60, 125.13, 0.1180
W_GRUT = 5.333984775933427e12      # Sigma vev scale h of the heavy reference (THEORY_DATA.json)
out = {}

# ============================================================ A. one-loop scalar betas from the general formula (sympy)
# V = (1/24) lam_abcd phi_a phi_b phi_c phi_d ; 16 pi^2 d lam_abcd / d ln mu = sum over the 3 channel pairings
def scalar_betas(V, fields, couplings):
    n = len(fields)
    L = [[[[sp.diff(V, fields[a], fields[b], fields[c], fields[d]) for d in range(n)] for c in range(n)]
          for b in range(n)] for a in range(n)]
    Vb = 0
    import itertools
    # beta-polynomial: (1/24) sum_abcd beta_abcd phi^4 with beta_abcd = sum_ef (L_abef L_efcd + L_acef L_efbd + L_adef L_efbc)
    # use the identity: (1/24) beta_abcd phi^4 = (1/8) * sum_ef (V_ef)^2 ... computed directly:
    He = sp.Matrix(n, n, lambda i, j: sp.diff(V, fields[i], fields[j]))     # field-dependent mass matrix
    # one-loop: dV/dln mu = (1/(32 pi^2)) Tr[(M^2)^2]  (the quartic part reproduces the beta functions)
    Vb = sp.expand((He * He).trace() / 2)          # times 1/(16 pi^2)
    return Vb


h1, h2, h3, h4, s1, s2 = sp.symbols('h1 h2 h3 h4 s1 s2', real=True)
lH, lS, lHS = sp.symbols('lam_H lam_S lam_HS', real=True)
HH = (h1 ** 2 + h2 ** 2 + h3 ** 2 + h4 ** 2) / 2
SS = (s1 ** 2 + s2 ** 2) / 2
Vc = lH * HH ** 2 + lS * SS ** 2 + 2 * lHS * HH * SS
Vb = scalar_betas(Vc, [h1, h2, h3, h4, s1, s2], [lH, lS, lHS])
# project on the monomials HH^2, SS^2, HH*SS at a test point
x, y = sp.symbols('x y', positive=True)
sub = {h1: sp.sqrt(2 * x), h2: 0, h3: 0, h4: 0, s1: sp.sqrt(2 * y), s2: 0}
Vb_xy = sp.expand(Vb.subs(sub))
bH = sp.simplify(Vb_xy.coeff(x, 2).coeff(y, 0))
bS = sp.simplify(Vb_xy.coeff(y, 2).coeff(x, 0))
bHS = sp.simplify(Vb_xy.coeff(x, 1).coeff(y, 1) / 2)
out['A_sympy_scalar_betas_complex_singlet'] = {'16pi^2 beta_lamH (scalar part)': str(bH), '16pi^2 beta_lamS': str(bS),
                                               '16pi^2 beta_lamHS (scalar part)': str(bHS),
                                               'EEGLS eq. RG': '24 lH^2 + 4 lHS^2 ; 8 lHS^2 + 20 lS^2 ; 4 lHS(3 lH + 2 lS) + 8 lHS^2'}
# real singlet, GRUT phi normalisation: V = lam (H^+H)^2 + b phi^2 H^+H + (lam_phi/4) phi^4
f_, b_, lp = sp.symbols('phi b lam_phi', real=True)
Vr = lH * HH ** 2 + b_ * f_ ** 2 * HH + lp / 4 * f_ ** 4
Vbr = sp.expand(scalar_betas(Vr, [h1, h2, h3, h4, f_], None).subs({h1: sp.sqrt(2 * x), h2: 0, h3: 0, h4: 0, f_: sp.sqrt(y)}))
out['A_sympy_scalar_betas_real_singlet'] = {'16pi^2 beta_lam (scalar part)': str(sp.simplify(Vbr.coeff(x, 2).coeff(y, 0))),
                                            '16pi^2 beta_b (scalar part)': str(sp.simplify(Vbr.coeff(x, 1).coeff(y, 1))),
                                            '16pi^2 beta_lam_phi': str(sp.simplify(4 * Vbr.coeff(y, 2).coeff(x, 0)))}
# complex singlet X with lam_HX |X|^2 H^+H (GRUT X normalisation, lam_X |X|^4)
X1, X2, lHX, lX = sp.symbols('X1 X2 lam_HX lam_X', real=True)
XX = (X1 ** 2 + X2 ** 2) / 2
Vx = lH * HH ** 2 + lHX * XX * HH + lX * XX ** 2
Vbx = sp.expand(scalar_betas(Vx, [h1, h2, h3, h4, X1, X2], None).subs({h1: sp.sqrt(2 * x), h2: 0, h3: 0, h4: 0,
                                                                       X1: sp.sqrt(2 * y), X2: 0}))
out['A_sympy_scalar_betas_complex_X'] = {'16pi^2 beta_lam (scalar part)': str(sp.simplify(Vbx.coeff(x, 2).coeff(y, 0))),
                                         '16pi^2 beta_lamHX (scalar part)': str(sp.simplify(Vbx.coeff(x, 1).coeff(y, 1))),
                                         '16pi^2 beta_lamX': str(sp.simplify(Vbx.coeff(y, 2).coeff(x, 0)))}
print(json.dumps({k: out[k] for k in out}, indent=1))

# ============================================================ B. SM running from the SMDR-matched PDG 2026 point
grid = json.load(open(os.path.join(RES, 'c1_fit_grid.json')))['fits']
c = [f for f in grid if np.allclose(f['inputs'], [MT, MH0, AS0])][0]
c0 = dict(g1=np.sqrt(5 / 3) * c['gp'], g2=c['g'], g3=c['g3'], yt=c['yt'], yb=c['yb'], ytau=c['ytau'], lam=c['k'])
smrun = R.run(c0, MT, mu_max=1.3e19, loops=3, n=2000)
Y_at = lambda mu: smrun['sol'].sol(np.log(mu / MT))
lam_sm = lambda mu: Y_at(mu)[6]
mus = np.logspace(np.log10(MT), np.log10(MPL), 3000)
lam_arr = np.array([lam_sm(m) for m in mus])
LAM_LAMBDA = float(np.exp(brentq(lambda t: lam_sm(np.exp(t)), np.log(1e9), np.log(1e12))))
out['B_SM'] = {'Lambda_lambda_GeV': LAM_LAMBDA, 'min_lambda': float(lam_arr.min()),
               'mu_at_min': float(mus[np.argmin(lam_arr)]), 'lambda(MPl)': float(lam_sm(MPL))}


# ============================================================ C. full theory above M_S (SM 3-loop + singlet 1-loop)
def rhs_full(t, Z, cplx=True):
    """t = ln mu. Z = (G1, G2, G3, T, B, L, lamH, Gamma, lamS, lamHS). EEGLS normalisation (complex S)."""
    y7 = Z[:7]
    d = 2.0 * R.betas(y7, 3)                       # d/d ln mu of G1..L, lam (SM 3-loop)
    G1, G2, G3, T = Z[0], Z[1], Z[2], Z[3]
    lH_, lS_, lHS_ = Z[6], Z[8], Z[9]
    gY2 = 3 / 5 * G1
    d[6] += 4 * lHS_ ** 2 / K16
    dS = (8 * lHS_ ** 2 + 20 * lS_ ** 2) / K16
    dHS = (0.5 * (12 * T - 3 * gY2 - 9 * G2) * lHS_ + 4 * lHS_ * (3 * lH_ + 2 * lS_) + 8 * lHS_ ** 2) / K16
    return np.concatenate([d, [R.gamma_h(y7, 3)], [dS, dHS]])


def run_full(MS, delta, lamS, sign=+1, mu_max=MPL):
    Y0 = Y_at(MS).copy()
    lHS0 = sign * np.sqrt(delta * lamS)
    Z0 = np.concatenate([Y0[:6], [Y0[6] + delta], [Y0[7]], [lamS, lHS0]])
    blow = lambda t, Z: 4 * np.pi * 2 - max(abs(Z[6]), abs(Z[8]), abs(Z[9]))
    blow.terminal = True
    sol = solve_ivp(rhs_full, (np.log(MS), np.log(mu_max)), Z0, method='DOP853', rtol=1e-9, atol=1e-14,
                    dense_output=True, events=blow)
    if sol.status == 1:
        sol.success = False
    return sol


def stable(MS, delta, lamS, w, sign=+1, nh=600, return_detail=False):
    """RG-improved tree-level test V(h, min_S) >= 0 for all h <= MPl, and perturbativity (all quartics < 4 pi)."""
    if MS > MPL:
        return False
    # below M_S: SM EFT lambda
    hb = np.logspace(np.log10(MT), np.log10(MS), 200)
    if np.min([lam_sm(m) for m in hb]) < 0:
        return (False, 'SM lambda < 0 below M_S') if return_detail else False
    sol = run_full(MS, delta, lamS, sign)
    if not sol.success:
        return (False, 'integration failed (Landau pole)') if return_detail else False
    hs = np.logspace(np.log10(MS), np.log10(MPL), nh)
    Z = sol.sol(np.log(hs))
    lH_, lS_, lHS_ = Z[6], Z[8], Z[9]
    if np.max(np.abs([lH_, lS_, lHS_])) > 4 * np.pi:
        return (False, 'non-perturbative (> 4 pi)') if return_detail else False
    if np.min(lS_) <= 0:
        return (False, 'lam_S < 0') if return_detail else False
    xh = hs ** 2 / 2
    ystar = w ** 2 / 2 - lHS_ * xh / lS_
    Vval = np.where(ystar > 0, (lH_ - lHS_ ** 2 / lS_) * xh ** 2, lH_ * xh ** 2 + lS_ * w ** 4 / 4 - lHS_ * xh * w ** 2)
    ok = bool(np.min(Vval) >= 0)
    if return_detail:
        k = int(np.argmin(Vval))
        return ok, {'min_V_over_h4': float(np.min(Vval / hs ** 4)), 'at_h': float(hs[k]),
                    'h_c_tree_GeV': float(w * np.sqrt(abs(lS_[0] / lHS_[0]))) if lHS_[0] > 0 else None,
                    'lamH_max': float(lH_.max()), 'lamS_MPl': float(lS_[-1]), 'lamHS_MPl': float(lHS_[-1])}
    return ok


def delta_min(MS, lamS_of_delta, w_of_delta, sign=+1, lo=1e-5, hi=2.0):
    """smallest delta in [lo, hi] giving stability (bisection on a monotone-in-delta assumption, checked on a grid)."""
    grid_d = np.logspace(np.log10(lo), np.log10(hi), 24)
    okv = [stable(MS, d, lamS_of_delta(d), w_of_delta(d), sign) for d in grid_d]
    if not any(okv):
        return None, None
    i = okv.index(True)
    a, b = (grid_d[i - 1] if i > 0 else lo), grid_d[i]
    for _ in range(16):
        m = np.sqrt(a * b)
        if stable(MS, m, lamS_of_delta(m), w_of_delta(m), sign):
            b = m
        else:
            a = m
    # largest stable delta (perturbativity upper edge), refined by bisection between grid points
    j = len(okv) - 1 - okv[::-1].index(True)
    if j == len(okv) - 1:
        return float(b), float(grid_d[j])
    a2, b2 = grid_d[j], grid_d[j + 1]
    for _ in range(16):
        m = np.sqrt(a2 * b2)
        if stable(MS, m, lamS_of_delta(m), w_of_delta(m), sign):
            a2 = m
        else:
            b2 = m
    return float(b), float(a2)


def mixing_angle(lamH, lamS, lamHS, w, v=VEW):
    M = 2 * np.array([[lamH * v ** 2, lamHS * v * w], [lamHS * v * w, lamS * w ** 2]])
    ev, U = np.linalg.eigh(M)
    return float(abs(U[0, 1]) if abs(U[0, 0]) > 0.5 else abs(U[0, 0])), np.sqrt(ev)


# ---------------- C1: generic EEGLS (free w), fixed lam_HS, positive and negative portal
gen = {}
for MS in (1e3, 1e8, 1e10, 3e10, 1e11, 3e11):
    row = {}
    for lHS_fix in (0.01, 0.1, 1.0):
        for sign in (+1, -1):
            lamS_of = lambda d, L=lHS_fix: L ** 2 / d
            w_of = lambda d, L=lHS_fix, M=MS: M / np.sqrt(2 * L ** 2 / d)
            dmin, dmax = delta_min(MS, lamS_of, w_of, sign)
            row[f'lamHS={lHS_fix:g}, sign={sign:+d}'] = {'delta_min': dmin, 'delta_max': dmax}
    gen[f'{MS:.0e}'] = row
    print('generic', MS, row)
out['C1_generic_EEGLS_free_vev'] = gen

# ---------------- C2: GRUT Sigma with the heavy-reference vev h = 5.334e12 GeV, M_Sigma free (lam_Sigma = M^2/(2 h^2))
grut = {}
for MS in (3e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 3e8, 1e9, 3e9, 1e10, 2e10, 3e10, 4e10, 1e11, 7.54e12):
    lamS = MS ** 2 / (2 * W_GRUT ** 2)
    row = {'lam_Sigma': lamS}
    for sign in (+1, -1):
        dmin, dmax = delta_min(MS, lambda d, L=lamS: L, lambda d: W_GRUT, sign)
        r = {'delta_min': dmin, 'delta_max': dmax}
        if dmin is not None:
            for lab, d in (('at_delta_min', dmin), ('at_delta_max', dmax)):
                lHS = sign * np.sqrt(d * lamS)
                th, masses = mixing_angle(lam_sm(MS) + d, lamS, lHS, W_GRUT)
                r[lab] = {'lam_HSigma_GRUTnorm': 2 * lHS, 'theta_mixing': th, 'theta_approx sqrt(2 delta) v/M': np.sqrt(2 * d) * VEW / MS,
                          'h_c_GeV': (W_GRUT * np.sqrt(lamS / lHS) if lHS > 0 else None),
                          'detail': stable(MS, d, lamS, W_GRUT, sign, return_detail=True)[1]}
        else:
            r['why_unstable_at_delta=0.02'] = stable(MS, 0.02, lamS, W_GRUT, sign, return_detail=True)[1]
            r['why_unstable_at_delta=0.5'] = stable(MS, 0.5, lamS, W_GRUT, sign, return_detail=True)[1]
        row[f'sign={sign:+d}'] = r
    grut[f'{MS:.3g}'] = row
    print('GRUT', MS, {k: (v.get('delta_min'), v.get('delta_max')) if isinstance(v, dict) else v for k, v in row.items()})
out['C2_GRUT_Sigma_fixed_vev'] = grut
# analytic edge for the positive branch: h_c = sqrt(w M_S) (2 delta)^(-1/4) <= Lambda  =>  M_S <= sqrt(2 delta) Lambda^2 / w
out['C2_tree_only_edge_estimate_SUPERSEDED_by_numerics_(ignores_8delta^2_running)'] = {f'M_S,max(delta={d})': np.sqrt(2 * d) * LAM_LAMBDA ** 2 / W_GRUT for d in (0.0105, 0.02, 0.1, 0.5, 1.0)}

# ---------------- C3: the heavy-reference point itself: lam_Sigma = 1 (M_Sigma = sqrt(2) h = 7.54e12 GeV)
MS_ref = np.sqrt(2) * W_GRUT
need = {}
for lHSig in (1e-10, 1e-3, 0.1, 1.0, 4 * np.pi):
    lHS = lHSig / 2
    need[f'lam_HSigma={lHSig:g}'] = {'h_c_GeV': W_GRUT * np.sqrt(1.0 / lHS), 'delta': lHS ** 2 / 1.0,
                                     'SM lambda < 0 between Lambda_lambda and M_Sigma': bool(MS_ref > LAM_LAMBDA)}
out['C3_heavy_reference_lamSigma=1'] = {'M_Sigma_GeV': MS_ref, 'Lambda_lambda_GeV': LAM_LAMBDA,
                                        'M_Sigma/Lambda_lambda': MS_ref / LAM_LAMBDA, 'portal_cases': need}


# ============================================================ D. light singlets: loop-only stabilisation (phi, X)
def rhs_light(t, Z, kind):
    """kind='phi': real singlet, V ⊃ b phi^2 H^+H + lam_phi/4 phi^4; kind='X': complex, lam_HX |X|^2 H^+H + lam_X |X|^4.
    Scalar parts from section A (sympy); H-leg anomalous dimension (6 y_t^2 - 9/2 g^2 - 3/2 g'^2) per H^+H."""
    y7 = Z[:7]
    d = 2.0 * R.betas(y7, 3)
    G1, G2, T, lam = Z[0], Z[1], Z[3], Z[6]
    gY2 = 3 / 5 * G1
    p, q = Z[8], Z[9]
    an = 6 * T - 4.5 * G2 - 1.5 * gY2
    if kind == 'phi':
        d[6] += 2 * p ** 2 / K16
        dp = (p * an + 12 * lam * p + 8 * p ** 2 + 6 * p * q) / K16
        dq = (8 * p ** 2 + 18 * q ** 2) / K16
    else:
        d[6] += p ** 2 / K16
        dp = (p * an + 12 * lam * p + 4 * p ** 2 + 8 * p * q) / K16
        dq = (2 * p ** 2 + 20 * q ** 2) / K16
    return np.concatenate([d, [R.gamma_h(y7, 3)], [dp, dq]])


def light_ok(p0, q0, kind):
    Y0 = Y_at(MT).copy()
    Z0 = np.concatenate([Y0[:8], [p0, q0]])
    blow = lambda t, Z, kind: 4 * np.pi * 2 - max(abs(Z[6]), abs(Z[8]), abs(Z[9]))
    blow.terminal = True
    sol = solve_ivp(rhs_light, (np.log(MT), np.log(MPL)), Z0, args=(kind,), method='DOP853', rtol=1e-9, atol=1e-14,
                    dense_output=True, events=blow)
    if not sol.success or sol.status == 1:
        return False
    Z = sol.sol(np.linspace(np.log(MT), np.log(MPL), 800))
    return bool(Z[6].min() >= 0 and np.abs(Z[6:]).max() < 4 * np.pi)


light = {}
for kind, q0, bmax, lab in (('phi', 0.03, 5.0501e-3, 'phi: b phi^2 H^+H, lam_phi = 0.03'),
                            ('X', 1e-3, 7.1419e-3, 'X: lam_HX |X|^2 H^+H, lam_X = 1e-3')):
    grid_p = np.linspace(0.0, 1.0, 101)
    okp = [light_ok(p, q0, kind) for p in grid_p]
    pmin = None
    if any(okp):
        i = okp.index(True)
        pmin = brentq(lambda p: (1 if light_ok(p, q0, kind) else -1), grid_p[i - 1], grid_p[i], xtol=1e-4) if i > 0 else 0.0
    light[lab] = {'portal_min_for_stability (one loop, singlet active from M_t)': pmin,
                  'invisible-width bound (ATLAS B_inv < 0.107, rp_portals [C])': bmax,
                  'ratio_needed_over_allowed': (pmin / bmax if pmin else None),
                  'GRUT detector-level portal (REVIEW_REHEATING)': '1e-10 to 1.5e-10'}
    # h -> ss width at the needed portal (b phi^2 H^+H: Gamma = b^2 v^2/(8 pi m_h); X complex: lam^2 v^2/(16 pi m_h))
    if pmin:
        G = (pmin ** 2 * VEW ** 2 / (8 * np.pi * MH0)) if kind == 'phi' else (pmin ** 2 * VEW ** 2 / (16 * np.pi * MH0))
        light[lab]['Gamma(h->ss)_GeV_at_needed_portal'] = G
        light[lab]['Gamma_h_SM_GeV'] = 4.1e-3
out['D_light_singlets'] = light
print(json.dumps(light, indent=1))

# ============================================================ E. Y (no vev, M_Y = 4e12 GeV)
out['E_Y'] = {'statement': 'Y has no vev (tree threshold needs a vev) and M_Y = 4.0e12 GeV > Lambda_lambda; its loops act only '
                           'above M_Y, so lambda_SM < 0 on (Lambda_lambda, M_Y) is untouched: cannot stabilise',
              'M_Y_over_Lambda_lambda': 4.0e12 / LAM_LAMBDA}

json.dump(out, open(os.path.join(RES, 'threshold.json'), 'w'), indent=1, default=float)
