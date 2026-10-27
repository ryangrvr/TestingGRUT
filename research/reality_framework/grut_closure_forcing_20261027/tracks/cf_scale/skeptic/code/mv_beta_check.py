#!/usr/bin/env python3
"""cf_scale_skeptic: independent re-derivation of the cf_scale one-loop beta functions (own code, no sympy, no
cf_scale modules imported).

Method (different from cf_scale's field-dependent-Hessian projection):
  * quartic tensor lambda_abcd = d^4 V4 / d phi_a d phi_b d phi_c d phi_d, obtained by exact 4-fold mixed forward
    differences of the homogeneous quartic V4 (exact for degree-4 polynomials);
  * Machacek-Vaughn / Luo-Wang-Xiao one-loop tensor formula for real scalars + one abelian gauge field
      16 pi^2 beta_abcd = Lambda^2_abcd - 3 Lambda^S_abcd + 3 A_abcd,
      Lambda^2 = lam_abef lam_efcd + 2 perms,  Lambda^S = (C_a + C_b + C_c + C_d) lam_abcd,  C_a = e^2 q_a^2,
      A_abcd = X_ab X_cd + 2 perms,  X_ab = {theta,theta}_ab = 2 e^2 q_a^2 delta_ab;
  * projection of the beta polynomial on the 11 operators by least squares at random field points (residual reported).
Compared with the closed forms printed in cf_scale results/beta_derive.json (read as data).
Also: (i) independent check of the one-loop monotonicity statement with a different minimiser and the envelope
formula; (ii) closed-form small-e critical-surface values of P4-M compared with cf_scale as_flow.json;
(iii) whether the P4-M charged-quartic values allow a Coleman-Weinberg minimum below M_P.
Usage: python3 -B mv_beta_check.py <out.json>
"""
import json, sys, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.dont_write_bytecode = True
TRACK = Path(__file__).resolve().parents[2]
NAMES = ['lam_X', 'lam_Y', 'lam_S', 'kap', 'l_XY', 'l_XS', 'l_YS', 'g2', 'l_pY', 'l_pS', 'lam_phi']
Q = np.array([1, 1, 3, 3, 6, 6, 0], float)       # real components x1 x2 y1 y2 s1 s2 phi
K16 = 16 * np.pi**2


def monos(f):
    x1, x2, y1, y2, s1, s2, p = f
    X = (x1 + 1j * x2) / np.sqrt(2); Y = (y1 + 1j * y2) / np.sqrt(2); S = (s1 + 1j * s2) / np.sqrt(2)
    aX, aY, aS = abs(X)**2, abs(Y)**2, abs(S)**2
    return np.array([aX**2, aY**2, aS**2, 2 * (np.conj(Y) * X**3).real, aX * aY, aX * aS, aY * aS,
                     p**2 * aX, p**2 * aY, p**2 * aS, p**4 / 4])


def V4(c, f):
    return float(np.dot(c, monos(f)))


def quartic_tensor(c):
    E = np.eye(7); T = np.zeros((7,) * 4)
    for a, b, cc, d in itertools.combinations_with_replacement(range(7), 4):
        s = 0.0
        for eps in itertools.product((0, 1), repeat=4):
            v = eps[0] * E[a] + eps[1] * E[b] + eps[2] * E[cc] + eps[3] * E[d]
            s += (-1)**(4 - sum(eps)) * V4(c, v)
        for perm in set(itertools.permutations((a, b, cc, d))):
            T[perm] = s
    return T


def beta_tensor(L, e):
    L2 = np.einsum('abef,efcd->abcd', L, L)
    Lam2 = L2 + np.einsum('acef,efbd->abcd', L, L) + np.einsum('adef,efbc->abcd', L, L)
    C = e**2 * Q**2
    LamS = (C[:, None, None, None] + C[None, :, None, None] + C[None, None, :, None] + C[None, None, None, :]) * L
    Xm = np.diag(2 * e**2 * Q**2)
    A = (np.einsum('ab,cd->abcd', Xm, Xm) + np.einsum('ac,bd->abcd', Xm, Xm) + np.einsum('ad,bc->abcd', Xm, Xm))
    return Lam2 - 3 * LamS + 3 * A          # = 16 pi^2 beta_abcd


def project(B, rng, npts=400):
    F = rng.normal(size=(npts, 7))
    M = np.array([monos(f) for f in F])
    rhs = np.einsum('abcd,na,nb,nc,nd->n', B, F, F, F, F) / 24.0
    coef, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    resid = float(np.max(np.abs(M @ coef - rhs)) / max(np.max(np.abs(rhs)), 1e-300))
    return coef, resid


def track_forms():
    js = json.loads((TRACK / 'results' / 'beta_derive.json').read_text())
    return js['beta_16pi2']


def eval_track(forms, c, e):
    ns = dict(zip(NAMES, c)); ns['e'] = e
    return np.array([eval(forms[n], {'__builtins__': {}}, ns) for n in NAMES])


def mv_beta(c, e, rng):
    return project(beta_tensor(quartic_tensor(c), e), rng)


def unit_min(c, rng, n=40, mask=None):
    mask = np.ones(7) if mask is None else mask
    best = None
    cons = {'type': 'eq', 'fun': lambda v: np.dot(v * mask, v * mask) - 1.0}
    for _ in range(n):
        v0 = rng.normal(size=7) * mask; v0 /= np.linalg.norm(v0)
        r = minimize(lambda v: V4(c, v * mask), v0, method='SLSQP', constraints=[cons], options=dict(ftol=1e-14, maxiter=500))
        if r.success and (best is None or r.fun < best.fun):
            best = r
    v = best.x * mask
    return best.fun, v / np.linalg.norm(v)


def main(out):
    rng = np.random.default_rng(20261027)
    res = dict(method=__doc__)
    forms = track_forms()
    # ---------------- 1. coefficient-by-coefficient comparison at random couplings
    worst = 0.0; rows = []
    for k in range(12):
        c = rng.normal(size=11) * 10**rng.uniform(-2, 0, size=11)
        e = 10**rng.uniform(-2, 0)
        mv, resid = mv_beta(c, e, rng)
        tr = eval_track(forms, c, e)
        rel = float(np.max(np.abs(mv - tr)) / np.max(np.abs(tr)))
        worst = max(worst, rel)
        rows.append(dict(e=e, max_rel_diff_vs_track=rel, projection_residual=resid))
    res['mv_vs_track'] = dict(n=len(rows), worst_rel_diff=worst, rows=rows)
    print('MV tensor betas vs cf_scale closed forms: worst rel diff = %.2e over %d random points' % (worst, len(rows)))
    # controls
    c = np.zeros(11); c[0] = 0.3
    mv, _ = mv_beta(c, 0.5, rng)
    res['control_scalarQED'] = dict(mv=float(mv[0]), expected=20 * 0.09 - 12 * 0.25 * 0.3 + 6 * 0.0625)
    # inhomogeneous gauge terms (all quartics zero): s_i
    mv0, _ = mv_beta(np.zeros(11), 1.0, rng)
    res['inhomogeneous_e4_coefficients'] = dict(zip(NAMES, [round(float(v), 6) for v in mv0]))
    print('inhomogeneous e^4 coefficients:', res['inhomogeneous_e4_coefficients'])
    # ---------------- 2. monotonicity, independent minimiser + envelope formula
    mono = []
    mask_all = np.ones(7)
    while len(mono) < 60:
        c = rng.normal(size=11) * 10**rng.uniform(-3, -0.5, size=11)
        e = 10**rng.uniform(-3, -0.3)
        lmin, n = unit_min(c, rng, n=30)
        if lmin > 0:
            continue
        mv, _ = mv_beta(c, e, rng)
        dl_env = V4(mv, n) / K16                     # envelope theorem: d lambda_min / dt = beta_V(n*)/16pi^2
        h = 1e-5
        lp, _ = unit_min(c + h * mv / K16, rng, n=30); lm, _ = unit_min(c - h * mv / K16, rng, n=30)
        mono.append(dict(lambda_min=lmin, e=e, d_envelope=dl_env, d_fd=(lp - lm) / (2 * h)))
    res['monotonicity'] = dict(n=len(mono), min_envelope=min(r['d_envelope'] for r in mono),
                               min_fd=min(r['d_fd'] for r in mono), rows=mono)
    print('monotonicity: %d cases with lambda_min<=0; min envelope derivative %.3e; min FD derivative %.3e' % (
        len(mono), res['monotonicity']['min_envelope'], res['monotonicity']['min_fd']))
    # ---------------- 3. P4-M closed-form small-e critical surface vs cf_scale as_flow.json
    s = mv0                                   # 16pi^2 beta_i at lambda = 0, e = 1
    fg = 9.83e-3; MP, MY = 2.435e18, 4.0e12; L = np.log(MP / MY)
    asf = json.loads((TRACK / 'results' / 'as_flow.json').read_text())
    comp = []
    for row in asf['C_P4M_critical_surface']['rows']:
        if row.get('diverged') or row['e_MP'] > 1e-3:
            continue
        A, eP = row['A_lam'], row['e_MP']
        pred_MP = -s * eP**4 / (K16 * (A + 4 * fg))
        pred_MY = pred_MP - s * eP**4 / K16 * L
        trMY = np.array([row['at_MY'][n] for n in NAMES])
        trMP = np.array([row['at_MP'][n] for n in NAMES])
        idx = [0, 1, 2, 4, 5, 6]
        comp.append(dict(A_lam=A, e_MP=eP, closed_form_MY={NAMES[i]: float(pred_MY[i]) for i in idx},
                         track_MY={NAMES[i]: float(trMY[i]) for i in idx},
                         max_rel_diff_MY=float(np.max(np.abs(pred_MY[idx] / trMY[idx] - 1))),
                         max_rel_diff_MP=float(np.max(np.abs(pred_MP[idx] / trMP[idx] - 1)))))
        print('P4-M A=%.4g e=%.1e: closed-form lam_S(M_Y)=%.3e  track=%.3e  max rel diff M_Y %.1e, M_P %.1e' % (
            A, eP, pred_MY[2], trMY[2], comp[-1]['max_rel_diff_MY'], comp[-1]['max_rel_diff_MP']))
    res['P4M_closed_form'] = dict(formula='lam_i(M_P) = -s_i e_P^4 / [16pi^2 (A_lam + 4 f_g)]; lam_i(mu) = lam_i(M_P) - '
                                  's_i e^4 ln(M_P/mu)/16pi^2 (small e, linear)', f_g=fg, rows=comp)
    # ---------------- 4. Coleman-Weinberg check along the pure-Sigma direction under P4-M
    # RG-improved V = lam(t)|S|^4: stationary point needs 4 lam + beta_lam = 0. With lam(M_P) = -s/(A+4fg) * e^4/16pi^2
    # and beta = s e^4/16pi^2 (small e): 4 lam + beta = (s e^4/16pi^2)[1 - 4/(A+4fg) - 4 ln(M_P/mu)] < 0 for all mu <= M_P
    # iff A + 4 f_g < 4.
    cw = {}
    for A in (0.0237, 0.056, 1.0, 2.83, 10.0):
        x = 1 - 4 / (A + 4 * fg)
        cw['A=%g' % A] = dict(bracket_at_MP=x, CW_minimum_below_MP=bool(x > 0),
                              mu_min_over_MP=float(np.exp(-x / 4)) if x > 0 else None)
    res['CW_pure_Sigma_P4M'] = dict(note='no stationary point of the RG-improved quartic potential below M_P unless '
                                         'A_lam + 4 f_g > 4; then only at mu ~ M_P, never near h = 5.33e12 GeV', rows=cw)
    print('CW:', cw)
    # vacuum-decay exponent for the P4-M quartic at the heavy-window e (illustrative, Fubini-Lipatov bounce for |S|^4)
    lamS = abs(comp[0]['closed_form_MY']['lam_S']) if comp else None
    if lamS:
        res['bounce_exponent_8pi2_over_3abslam'] = 8 * np.pi**2 / (3 * lamS)
        print('bounce exponent for |lam_S| = %.2e: S_E ~ %.2e (vacuum decay irrelevant)' % (lamS, 8 * np.pi**2 / (3 * lamS)))
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main(sys.argv[1])
