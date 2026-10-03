"""BRI1 analytic preflight -- EXACT SYMBOLIC checks only (sympy power series; no sampling, no Monte Carlo, no quadrature).
Independent code path, not independent reviewer.
(1) PF-1: C1 bath elimination identity  F = A'(q)[F_free - int_0^t gamma(t-s) d/ds A(q(s)) ds]  for one harmonic mode,
    A(q) = q + q^3/3, checked as an exact Taylor identity in t for a generic clamped polynomial q(t) with q(0)=q'(0)=0.
(2) PF-4: X1 single-Duffing first-order coefficient c(t) = Cov(x0(t)^2, y1(t)) for protocol P1 on [0, pi]
    (q = s(t/pi)), as an exact small-t series with Gibbs moments reduced to m2 = E[a^2] by the identity
    m_{k+1} + m_{k+3} = k m_{k-1} (integration by parts against exp(-a^2/2 - a^4/4)).
(3) PF-3/PF-5 parity facts are structural (x -> -x symmetry) and proved in the charter text.
"""
import sympy as sp

t, s_, w, c, T0 = sp.symbols('t s omega c T0', positive=True)
ORDER = 12

# ---------------- (1) C1 identity, one mode, exact series ----------------
q2, q3, q4 = sp.symbols('q2 q3 q4')
q = q3 * t**3 + q4 * t**4                      # clamped prefix with q(0)=q'(0)=0 (generic cubic/quartic start)
A = q + q**3 / 3; Ap = 1 + q**2
x0, p0 = sp.symbols('x0 p0')                   # bath initial data at q(0)=0 (shifted = unshifted since A(0)=0)
# exact bath solution by variation of constants
u = sp.symbols('u')
xs = x0 * sp.cos(w * t) + p0 / w * sp.sin(w * t) + sp.integrate((c * A.subs(t, u)) * sp.sin(w * (t - u)) / w, (u, 0, t))
force = c * Ap * (xs - c * A / w**2)            # F = -dH/dq for H_B = p^2/2 + w^2/2 (x - c A(q)/w^2)^2
Ffree = c * (x0 * sp.cos(w * t) + p0 / w * sp.sin(w * t))
gamma = lambda tau: c**2 / w**2 * sp.cos(w * tau)
mem = sp.integrate(gamma(t - u) * sp.diff(A, t).subs(t, u), (u, 0, t))
rhs = Ap * (Ffree - mem)
d1 = sp.series(sp.simplify(force - rhs), t, 0, ORDER).removeO()
print("(1) C1: force - A'(q)[F_free - int gamma dA] series through t^%d:" % (ORDER - 1), sp.simplify(sp.expand(d1)))

# ---------------- (2) X1 single-oscillator first-order coefficient, exact small-t series ----------------
a, b = sp.symbols('a b')                       # x0(0)=a ~ Gibbs exp(-a^2/2-a^4/4), p0(0)=b ~ N(0,1), independent
N_ = 14
# Taylor solution of x'' = -x - x^3
xc = [a, b]
for k in range(2, N_ + 1):
    X = sum(xc[i] * t**i for i in range(k))
    acc = sp.expand(-(X + X**3)).coeff(t, k - 2)
    xc.append(acc / (k * (k - 1)))
xser = sum(xc[i] * t**i for i in range(N_ + 1))
# P1 on [0,pi]: q = 10 u^3 - 15 u^4 + 6 u^5, u = t/pi
qq = 10 * (t / sp.pi)**3 - 15 * (t / sp.pi)**4 + 6 * (t / sp.pi)**5
# linearised response y'' = -(1+3 x0^2) y + q, y(0)=y'(0)=0
yc = [0, 0]
for k in range(2, N_ + 1):
    Y = sum(yc[i] * t**i for i in range(k))
    acc = sp.expand(-(1 + 3 * xser**2) * Y + qq)
    acc = sp.series(acc, t, 0, k - 1).removeO().coeff(t, k - 2) if k - 2 > 0 else acc.subs(t, 0)
    yc.append(sp.expand(acc) / (k * (k - 1)))
yser = sum(yc[i] * t**i for i in range(N_ + 1))
x2 = sp.expand(sp.series(sp.expand(xser**2), t, 0, N_ + 1).removeO())
prod = sp.expand(sp.series(sp.expand(x2 * yser), t, 0, N_ + 1).removeO())
# expectation operator: b ~ N(0,1); a moments reduced to m2
m2 = sp.symbols('m2', positive=True)
am = {0: sp.Integer(1), 2: m2}
for k in range(1, 40, 2):                      # m_{k+3} = k m_{k-1} - m_{k+1}
    am[k + 3] = sp.expand(k * am[k - 1] - am[k + 1])
def E(expr):
    expr = sp.expand(expr); out = 0
    for term in sp.Add.make_args(expr):
        pa = sp.Poly(term, a, b); val = 0
        for (ea, eb), co in pa.terms():
            if ea % 2 or eb % 2: continue
            val += co * am[ea] * sp.factorial2(eb - 1)
        out += val
    return sp.expand(out)
cser = sp.expand(E(prod) - E(x2) * E(yser))
lead = [sp.factor(cser.coeff(t, k)) for k in range(0, N_ + 1)]
for k, v in enumerate(lead):
    if v != 0:
        print(f"(2) X1/P1: c(t) coefficient of t^{k}: {v}")
print("(2) check Var(a^2) = m4 - m2^2 =", sp.factor(am[4] - m2**2), " (strictly positive for a non-degenerate a^2)")
