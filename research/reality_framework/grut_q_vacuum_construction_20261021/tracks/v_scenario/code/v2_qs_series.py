"""v_scenario / v2 helper: high-order quasi-static (adiabatic) series for the q-mode, my own derivation.

Unknowns y(N) = sum_k e^k y_k (e = 1/m^2, K0 = 1) and P(N) = H^2 = sum_k e^k p_k, with d/dt = H d/dN:
  (i)  P y'' + (P'/2) y' + 3 P y' + F(y)/e = S/3,          F(y) = y + 2y^2 + (4/3)y^3 + y^4/3 (KV quartic)
  (ii) 3 P (1 + y + y') = rho + (3/(4e)) y^2 (1 + 4y/3 + y^2/3)
solved order by order with generic S(N), rho(N). Returns lambdified y_k, p_k in terms of the derivatives
S^(j), rho^(j) at a point. Used for clean initial data (the scenario used an O(m^-4) track; this goes to O(m^-10)).
"""
import sympy as sp

N, e = sp.symbols('N e')
Sf = sp.Function('S')(N)
rf = sp.Function('rho')(N)


def build(order=5):
    ys = [sp.Function('y%d' % k)(N) for k in range(1, order + 1)]
    ps = [sp.Function('p%d' % k)(N) for k in range(0, order + 1)]
    y = sum(e ** (k + 1) * ys[k] for k in range(order))
    P = sum(e ** k * ps[k] for k in range(order + 1))
    F = y + 2 * y ** 2 + sp.Rational(4, 3) * y ** 3 + y ** 4 / 3
    eq1 = sp.expand(P * y.diff(N, 2) + P.diff(N) / 2 * y.diff(N) + 3 * P * y.diff(N) + F / e - Sf / 3)
    eq2 = sp.expand(3 * P * (1 + y + y.diff(N)) - rf - sp.Rational(3, 4) / e * y ** 2 * (1 + 4 * y / 3 + y ** 2 / 3))
    sol = {}
    # order e^0: eq1 -> y1 ; eq2 -> p0
    for k in range(order):
        c1 = eq1.coeff(e, k)
        c1 = c1.subs(sol).doit()
        yk = sp.solve(c1, ys[k])[0]
        sol[ys[k]] = sp.simplify(yk)
        # substitute derivatives of solved functions
        c2 = eq2.coeff(e, k).subs(sol).doit()
        pk = sp.solve(c2, ps[k])[0]
        sol[ps[k]] = sp.simplify(pk)
        # make derivative substitutions resolve: replace function by expression (doit evaluates derivatives)
        sol = {kk: vv.subs(sol).doit() for kk, vv in sol.items()}
    # collect derivative symbols
    nder = 2 * order + 2
    Sd = sp.symbols('S0:%d' % nder)
    Rd = sp.symbols('R0:%d' % nder)
    rep = {}
    for j in range(nder - 1, 0, -1):
        rep[Sf.diff(N, j)] = Sd[j]
        rep[rf.diff(N, j)] = Rd[j]
    rep[Sf] = Sd[0]
    rep[rf] = Rd[0]

    def sub(expr):
        for j in range(nder - 1, 0, -1):
            expr = expr.subs(Sf.diff(N, j), Sd[j]).subs(rf.diff(N, j), Rd[j])
        return expr.subs(Sf, Sd[0]).subs(rf, Rd[0])
    yk_l = [sp.lambdify((Sd, Rd), sub(sol[ys[k]]), 'mpmath') for k in range(order)]
    ykN_l = [sp.lambdify((Sd, Rd), sub(sol[ys[k]].diff(N)), 'mpmath') for k in range(order)]
    pk_l = [sp.lambdify((Sd, Rd), sub(sol[ps[k]]), 'mpmath') for k in range(order)]
    exprs = dict(y1=str(sol[ys[0]]), p0=str(sol[ps[0]]), y2=str(sol[ys[1]]) if order > 1 else None,
                 p1=str(sol[ps[1]]) if order > 1 else None)
    return yk_l, ykN_l, pk_l, nder, exprs


if __name__ == '__main__':
    yk, ykN, pk, nder, ex = build(3)
    print(ex)
