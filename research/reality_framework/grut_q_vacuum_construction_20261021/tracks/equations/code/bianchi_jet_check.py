"""Exact (rational-arithmetic) check of the off-shell Noether/Bianchi identity of q-theory
on random, non-symmetric 4D configurations.

Identity to verify (derived in NOTES.md section 1.4):
  E_mn := K(q) G_mn + (eps(q) - mu(x) q) g_mn + (g_mn Box - nabla_m nabla_n) K(q),
  mu(x) := eps'(q) - K'(q) R / 2,
  nabla^m E_mn = - q d_n mu         (for ANY metric g and ANY scalar q(x)).
On shell E_mn = T^m_mn, so nabla^m T^m_mn = -q d_n mu: matter conservation <=> mu const.

Method: all fields are represented by their Taylor polynomials (degree <= 3) about the origin
with random rational coefficients; every intermediate quantity is truncated at the degree
that is still exact (g^-1: 3, Gamma: 2, Riemann/G/E: 1, divergence: 0). The identity is
then checked at the origin in exact rational arithmetic. This is equivalent to a symbolic
check at a generic point for generic fields up to the derivative order involved.
Negative controls: wrong sign of the mu q term; dropped (g Box - nabla nabla)K term;
mu replaced by a constant (then the residual must equal q d mu, nonzero).
"""
import json, os, random
from fractions import Fraction
from sympy import ring, QQ, Rational

Rg, T, X, Y, Z = ring("t,x,y,z", QQ)
V = [T, X, Y, Z]
gens = V


def trunc(p, d):
    return Rg({m: c for m, c in p.items() if sum(m) <= d})


def rnd(scale=1, den=7):
    return QQ(random.randint(-9 * scale, 9 * scale), random.randint(1, den) * 10)


def rand_poly(deg, const=None, scale=1):
    p = Rg(0)
    monos = []
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            for k in range(deg + 1 - i - j):
                for l in range(deg + 1 - i - j - k):
                    monos.append((i, j, k, l))
    for m in monos:
        if sum(m) == 0:
            continue
        p += Rg({m: rnd(scale)})
    if const is not None:
        p += const
    return p


def mat_inv_const(M):
    # exact inverse of 4x4 constant (QQ) matrix via sympy Matrix
    from sympy import Matrix
    A = Matrix(4, 4, lambda i, j: Rational(str(M[i][j])))
    Ai = A.inv()
    return [[QQ(int(Ai[i, j].p), int(Ai[i, j].q)) for j in range(4)] for i in range(4)]


def const_part(p):
    return p.coeff(1) if p != 0 else QQ(0)


def compose(f_coeffs, center, s):
    """f(s) = sum_n c_n (s - center)^n, s a polynomial; exact composition, then truncated later."""
    u = s - center
    res = Rg(0)
    pw = Rg(1)
    for c in f_coeffs:
        res += c * pw
        pw = trunc(pw * u, 3)
    return trunc(res, 3)


def dcoeffs(cs):
    return [cs[n] * n for n in range(1, len(cs))]


def run(metric, qpoly, Kc, ec, qc, variant='correct'):
    g = metric
    g0 = [[const_part(g[i][j]) for j in range(4)] for i in range(4)]
    g0i = mat_inv_const(g0)
    dlt = [[g[i][j] - g0[i][j] for j in range(4)] for i in range(4)]
    # g^-1 = sum_k (-g0^-1 dlt)^k g0^-1, truncated at degree 3
    Mneg = [[-sum(g0i[i][k] * dlt[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    gi = [[Rg(g0i[i][j]) for j in range(4)] for i in range(4)]
    term = [[Rg(g0i[i][j]) for j in range(4)] for i in range(4)]
    for _ in range(3):
        term = [[trunc(sum(Mneg[i][k] * term[k][j] for k in range(4)), 3) for j in range(4)] for i in range(4)]
        gi = [[gi[i][j] + term[i][j] for j in range(4)] for i in range(4)]
    # check inverse
    for i in range(4):
        for j in range(4):
            e = trunc(sum(g[i][k] * gi[k][j] for k in range(4)), 3)
            assert e == (1 if i == j else 0), 'inverse failed'
    dg = [[[g[i][j].diff(V[k]) for k in range(4)] for j in range(4)] for i in range(4)]
    Gam = [[[trunc(sum(gi[l][s] * (dg[s][n][m] + dg[s][m][n] - dg[m][n][s]) for s in range(4)) * QQ(1, 2), 2)
             for n in range(4)] for m in range(4)] for l in range(4)]
    # Riemann R^r_{s m n}, degree 1
    def Riem(r, s, m, n):
        e = Gam[r][n][s].diff(V[m]) - Gam[r][m][s].diff(V[n])
        e += sum(Gam[r][m][l] * Gam[l][n][s] - Gam[r][n][l] * Gam[l][m][s] for l in range(4))
        return trunc(e, 1)
    Ric = [[trunc(sum(Riem(r, s, r, n) for r in range(4)), 1) for n in range(4)] for s in range(4)]
    Rs = trunc(sum(gi[i][j] * Ric[i][j] for i in range(4) for j in range(4)), 1)
    G = [[trunc(Ric[i][j] - g[i][j] * Rs * QQ(1, 2), 1) for j in range(4)] for i in range(4)]
    q = qpoly
    K = compose(Kc, qc, q)
    K1 = compose(dcoeffs(Kc), qc, q)
    eps = compose(ec, qc, q)
    e1 = compose(dcoeffs(ec), qc, q)
    mu = trunc(e1 - K1 * Rs * QQ(1, 2), 1)
    if variant == 'mu_const':
        # constant mu (as in the field equation (g)) that differs from the local mu(x) at the origin by 1/3:
        # expected  nabla^m E_mn = (mu(x) - mu_const) d_n q = -(1/3) d_n q  at the origin
        mu = Rg(const_part(mu) + QQ(1, 3))
    dK = [K.diff(V[m]) for m in range(4)]
    DDK = [[trunc(K.diff(V[m]).diff(V[n]) - sum(Gam[l][m][n] * dK[l] for l in range(4)), 1) for n in range(4)] for m in range(4)]
    BoxK = trunc(sum(gi[m][n] * DDK[m][n] for m in range(4) for n in range(4)), 1)
    sgn = -1 if variant == 'wrong_sign_muq' else 1
    dterm = 0 if variant == 'drop_K_derivative_terms' else 1
    E = [[trunc(K * G[m][n] + (eps - sgn * mu * q) * g[m][n] + dterm * (g[m][n] * BoxK - DDK[m][n]), 1)
          for n in range(4)] for m in range(4)]
    # divergence at origin
    res = []
    for n in range(4):
        e = Rg(0)
        for m in range(4):
            for r in range(4):
                cd = E[r][n].diff(V[m]) - sum(Gam[l][m][r] * E[l][n] + Gam[l][m][n] * E[r][l] for l in range(4))
                e += gi[m][r] * cd
        e = const_part(trunc(e, 0))
        qdmu = const_part(trunc(q * mu.diff(V[n]), 0))
        if variant == 'mu_const':
            qdmu = QQ(1, 3) * const_part(trunc(q.diff(V[n]), 0))   # residual e + qdmu should vanish
        res.append((e, qdmu))
    return res


def random_metric(kind):
    eta = [-1, 1, 1, 1]
    if kind == 'generic':
        g = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i, 4):
                c = QQ(eta[i]) if i == j else QQ(0)
                c += rnd() / 5
                p = rand_poly(3, const=c, scale=1) * QQ(1, 3)
                p = p - const_part(p) + c
                g[i][j] = p
                g[j][i] = p
        return g
    if kind == 'frw':
        a = rand_poly_1d(3, const=QQ(1))
        a2 = a * a
        Zr = Rg(0)
        return [[Rg(-1), Zr, Zr, Zr], [Zr, a2, Zr, Zr], [Zr, Zr, a2, Zr], [Zr, Zr, Zr, a2]]
    raise ValueError


def rand_poly_1d(deg, const):
    p = Rg(const)
    for n in range(1, deg + 1):
        p += rnd() * T**n
    return p


def main():
    random.seed(20261021)
    out = {'description': __doc__.strip().split('\n')[0], 'cases': []}
    for kind in ['generic', 'generic', 'generic', 'frw', 'frw']:
        g = random_metric(kind)
        if kind == 'frw':
            q = rand_poly_1d(3, const=QQ(7, 5))
        else:
            q = rand_poly(3, const=QQ(7, 5))
        Kc = [QQ(1)] + [rnd() for _ in range(3)]      # K(q) = 1 + k1 s + k2 s^2 + k3 s^3, s = q - qc
        ec = [rnd() for _ in range(4)]               # eps(q) cubic in s
        qc = QQ(1, 2)
        for variant in ['correct', 'wrong_sign_muq', 'drop_K_derivative_terms', 'mu_const']:
            res = run(g, q, Kc, ec, qc, variant)
            if variant == 'mu_const':
                resid = [e + qd for e, qd in res]
                out['cases'].append({'metric': kind, 'variant': 'constant mu offset by 1/3 from local mu(x)',
                                     'check': 'div E_const = (mu(x)-mu_const) d q',
                                     'residual': [str(r) for r in resid],
                                     'div_E_const': [str(e) for e, _ in res],
                                     'identity_holds_exactly': all(r == 0 for r in resid)})
            else:
                resid = [e + qd for e, qd in res]
                out['cases'].append({'metric': kind, 'variant': variant,
                                     'residual_divE_plus_q_dmu': [str(r) for r in resid],
                                     'q_dmu': [str(qd) for _, qd in res],
                                     'identity_holds_exactly': all(r == 0 for r in resid)})
            print(out['cases'][-1])
    p = os.path.join(os.path.dirname(__file__), '..', 'results', 'bianchi_jet_check.json')
    with open(p, 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
