#!/usr/bin/env python3
"""cf_scale: sympy derivation of the one-loop beta functions of the GRUT hidden sector (own code).

Fields (real components): X = (x1 + i x2)/sqrt2, Y, Sigma likewise, phi real. Hidden U(1) charges 1, 3, 6, gauge e.
Tree quartic potential (GRUT THEORY.md sec. 2 normalisation, plus all gauge-invariant Z2_phi-even quartics):
  V4 = lam_X|X|^4 + lam_Y|Y|^4 + lam_S|S|^4 + 2 kap Re(Y* X^3) + l_XY|X|^2|Y|^2 + l_XS|X|^2|S|^2 + l_YS|Y|^2|S|^2
       + g2 phi^2|X|^2 + l_pY phi^2|Y|^2 + l_pS phi^2|S|^2 + (lam_phi/4) phi^4          (g2 = g^2 of Candidate 1)
(The tree term lam_S(|S|^2 - h^2/2)^2 has quartic lam_S|S|^4; kap is taken real, as in hw_gates gate 4.)

Method [D]: the one-loop RG of the quartic part of V in Landau gauge,
  16 pi^2 mu dV4/dmu = (1/2) Tr[(V4'')^2] + (3/2) Tr[(M_A^2)^2] - 3 e^2 sum_i q_i^2 phi_i dV4/dphi_i,
with V4'' the field-dependent scalar mass matrix (all real components, Goldstones included) and
M_A^2 = 2 e^2 sum_z q_z^2 |z|^2 the vector mass. Equivalent to the Machacek-Vaughn tensor formula
Lambda^2 - 3 g^2 Lambda^S + 3 g^4 A. Gauge: 16 pi^2 de/dmu = (e^3/3) sum_z q_z^2 (complex scalars).
The result is projected on the 11-operator basis (closure checked: zero residual).
Controls: scalar QED 16pi^2 beta = 20 lam^2 - 12 e^2 lam + 6 e^4; lone real phi^4/4: 18 lam^2; hw_gates gate 4 numbers.
Usage: python3 -B beta_derive.py <out.json>
"""
import json, sys
from pathlib import Path
import sympy as sp

sys.dont_write_bytecode = True

NAMES = ['lam_X', 'lam_Y', 'lam_S', 'kap', 'l_XY', 'l_XS', 'l_YS', 'g2', 'l_pY', 'l_pS', 'lam_phi']
x1, x2, y1, y2, s1, s2, p = F = sp.symbols('x1 x2 y1 y2 s1 s2 p', real=True)
CH = [1, 1, 3, 3, 6, 6, 0]
e = sp.Symbol('e', positive=True)
C = dict(zip(NAMES, sp.symbols(NAMES, real=True)))
X = (x1 + sp.I * x2) / sp.sqrt(2); Y = (y1 + sp.I * y2) / sp.sqrt(2); S = (s1 + sp.I * s2) / sp.sqrt(2)
aX, aY, aS = [sp.expand(z * sp.conjugate(z)) for z in (X, Y, S)]
MONO = dict(lam_X=aX**2, lam_Y=aY**2, lam_S=aS**2, kap=sp.expand(2 * sp.re(sp.expand(sp.conjugate(Y) * X**3))),
            l_XY=aX * aY, l_XS=aX * aS, l_YS=aY * aS, g2=p**2 * aX, l_pY=p**2 * aY, l_pS=p**2 * aS, lam_phi=p**4 / 4)
MONO = {k: sp.expand(v) for k, v in MONO.items()}


def beta_V(V4, gauge=True):
    Hm = sp.Matrix(7, 7, lambda i, j: sp.diff(V4, F[i], F[j]))
    tr = sp.expand((Hm * Hm).trace())
    out = tr / 2
    if gauge:
        MA2 = 2 * e**2 * sum(q**2 * a for q, a in zip((1, 3, 6), (aX, aY, aS)))
        out += sp.Rational(3, 2) * sp.expand(MA2**2)
        out -= 3 * e**2 * sum(CH[i]**2 * F[i] * sp.diff(V4, F[i]) for i in range(7))
    return sp.expand(out)


def project(Bpoly):
    """Solve Bpoly = sum_k b_k MONO[k] exactly; return dict and residual polynomial."""
    b = sp.symbols('b0:%d' % len(NAMES))
    resid = sp.expand(Bpoly - sum(bi * MONO[n] for bi, n in zip(b, NAMES)))
    eqs = sp.Poly(resid, *F).coeffs()
    sol = sp.solve(eqs, b, dict=True)
    assert len(sol) == 1, 'basis not closed'
    sol = sol[0]
    res = sp.expand(resid.subs(sol))
    return {n: sp.factor(sp.expand(sol[bi])) for bi, n in zip(b, NAMES)}, res


def main(out):
    V4 = sum(C[n] * MONO[n] for n in NAMES)
    B = beta_V(V4)
    betas, res = project(B)
    assert res == 0, res
    beta_e = sp.Rational(1, 3) * (1 + 9 + 36) * e**3
    res_json = dict(method=__doc__, convention='16 pi^2 beta_i = expression below; beta = d/d ln mu',
                    beta_16pi2={k: str(sp.expand(v)) for k, v in betas.items()}, beta_e_16pi2=str(beta_e))
    for k in NAMES:
        print('16pi^2 beta(%s) = %s' % (k, sp.expand(betas[k])))
    print('16pi^2 beta(e) =', beta_e)
    # ---- controls
    lam = sp.Symbol('lam')
    ctrl = beta_V(lam * aX**2, gauge=False) + 0  # scalar part only, X alone (Y,S,p fields still present but decouple)
    bq, _ = project(beta_V(lam * aX**2))
    res_json['control_scalarQED_lamX'] = str(sp.expand(bq['lam_X']))
    print('control scalar QED (charge 1):', sp.expand(bq['lam_X']), ' expected 20 lam^2 - 12 e^2 lam + 6 e^4')
    bp, _ = project(beta_V(lam * p**4 / 4))
    res_json['control_real_phi4'] = str(sp.expand(bp['lam_phi']))
    print('control lone real phi (V = lam phi^4/4):', bp['lam_phi'], ' expected 18 lam^2')
    # ---- X-phi-e sector closed forms (Y = Sigma = 0 couplings)
    sub = {C[n]: 0 for n in NAMES if n not in ('lam_X', 'g2', 'lam_phi')}
    xp = {k: sp.factor(sp.expand(betas[k].subs(sub))) for k in ('lam_X', 'g2', 'lam_phi')}
    # in the X-phi sector the U(1) still has the charge-3,6 fields as spectators for beta_e only
    res_json['X_phi_sector_16pi2'] = {k: str(sp.expand(v)) for k, v in xp.items()}
    print('X-phi-e sector:'); [print('   16pi^2 beta(%s) = %s' % (k, sp.expand(v))) for k, v in xp.items()]
    # ---- invariant surfaces [D]
    phi_surface = {C['g2']: 0, C['l_pY']: 0, C['l_pS']: 0}
    inv_phi = {k: sp.simplify(betas[k].subs(phi_surface)) for k in ('g2', 'l_pY', 'l_pS')}
    inv_phi_full = {k: sp.simplify(betas[k].subs(phi_surface).subs(C['lam_phi'], 0)) for k in ('lam_phi',)}
    inv_kap = sp.simplify(betas['kap'].subs(C['kap'], 0))
    res_json['invariant_surfaces'] = dict(
        phi_portals_zero=({k: str(v) for k, v in inv_phi.items()}),
        lam_phi_on_decoupled_surface=str(inv_phi_full['lam_phi']),
        kap_zero=str(inv_kap))
    print('beta of phi portals on {g2=l_pY=l_pS=0}:', inv_phi, '; beta(lam_phi) there with lam_phi=0:', inv_phi_full,
          '; beta(kap) at kap=0:', inv_kap)
    # ---- gauge-only inhomogeneous terms (all quartics zero) [D]
    zero = {C[n]: 0 for n in NAMES}
    inh = {k: sp.expand(betas[k].subs(zero)) for k in NAMES}
    res_json['inhomogeneous_at_zero_quartics_16pi2'] = {k: str(v) for k, v in inh.items()}
    print('gauge-only terms at lambda = 0:', inh)
    # positivity of the generated V4 at lambda = 0: beta_V = (3/2) Tr M_A^4 >= 0 in every direction [D]
    bV0 = sp.factor(B.subs(zero))
    res_json['beta_V_at_zero_quartics_16pi2'] = str(bV0)
    print('16pi^2 beta_V at lambda = 0 :', bV0)
    # ---- numeric comparison with hw_gates gate 4 (pure scalar, D couplings zero), e = 0
    start = dict(lam_X=1e-3, lam_Y=1.0, lam_S=1.0, kap=0.0034586847211042523, l_XY=0., l_XS=0., l_YS=0., g2=0.01,
                 l_pY=0., l_pS=0., lam_phi=0.03)
    num = {k: float(betas[k].subs({C[n]: start[n] for n in NAMES}).subs(e, 0) / (16 * sp.pi**2)) for k in NAMES}
    res_json['numeric_beta_at_spot_e0'] = num
    print('numeric betas at the spot point (e=0):', {k: '%.4e' % v for k, v in num.items()})
    Path(out).write_text(json.dumps(res_json, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
