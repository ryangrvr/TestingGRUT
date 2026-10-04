"""VER0-B blind reproduction, R3: small-time series of K_{P1}(t) = 3 Cov_mu( X0(t)^2 , Y(t) ).

X0(t): unforced oscillator  x'' + x + x^3 = 0,  x(0)=x, x'(0)=p.
Y(t) : first epsilon-variation  Y'' + (1 + 3 X0(t)^2) Y = q(t),  Y(0)=Y'(0)=0,
       q(t) = 10 (t/pi)^3 - 15 (t/pi)^4 + 6 (t/pi)^5  (P1 on [0, pi]).
Gibbs mu: density prop. exp(-(p^2/2 + x^2/2 + x^4/4)); p ~ N(0,1) independent of x.
x-moments reduced by integration by parts:  (n-1) m_{n-2} = m_n + m_{n+2}  (all in terms of m2).

This script only produces the formal Taylor coefficients (exact polynomial algebra);
the analytic remainder bound is proved in the reproduction document, not here.
"""
import sympy as sp

t, x, p, m2 = sp.symbols('t x p m2')
ORDER = 11  # keep terms t^0 .. t^ORDER
PI = sp.pi


def trunc(expr, n=ORDER):
    expr = sp.expand(expr)
    return sum(expr.coeff(t, k) * t**k for k in range(n + 1))


def integrate0(expr):
    return trunc(sp.integrate(sp.expand(expr), (t, 0, t)))


# Picard iteration for X0 (each iteration gains >= 2 orders in t)
X0 = x + p * t
for _ in range(ORDER):
    V = p + integrate0(-(X0 + X0**3))
    X0 = trunc(x + integrate0(V))

q = 10 * (t / PI)**3 - 15 * (t / PI)**4 + 6 * (t / PI)**5
Y = sp.Integer(0)
for _ in range(ORDER):
    W = integrate0(q - (1 + 3 * X0**2) * Y)
    Y = trunc(integrate0(W))

# Gibbs moments
xmom = {0: sp.Integer(1), 2: m2}
for n in range(1, 40):
    # (n-1) m_{n-2} = m_n + m_{n+2}  for even n>=2  (n=1 trivially odd)
    if n % 2 == 0:
        xmom[n + 2] = sp.expand((n - 1) * xmom[n - 2] - xmom[n])


def pmom(k):
    return sp.Integer(0) if k % 2 else sp.factorial2(k - 1) if k > 0 else sp.Integer(1)


def xm(k):
    return sp.Integer(0) if k % 2 else xmom[k]


def E(poly):
    poly = sp.Poly(sp.expand(poly), x, p)
    out = sp.Integer(0)
    for (i, j), c in poly.terms():
        out += c * xm(i) * pmom(j)
    return sp.expand(out)


# sanity: invariance of Gibbs under unforced flow => E[X0(t)^2] = m2 to all computed orders
EX02 = trunc(E(trunc(X0**2)))
print("E[X0(t)^2] series            :", sp.simplify(EX02))
EX03 = trunc(E(trunc(X0**3)))
print("E[X0(t)^3] series (parity)   :", sp.simplify(EX03))

K = trunc(3 * (E(trunc(X0**2 * Y)) - m2 * E(Y)))
print("\nY(t) leading terms:")
for k in range(0, 9):
    c = sp.factor(sp.expand(Y).coeff(t, k))
    print(f"  t^{k}: {c}")

print("\nK(t) = 3 Cov(X0(t)^2, Y(t)) coefficients (m2 = E_mu[x^2]):")
for k in range(0, ORDER + 1):
    c = sp.factor(sp.expand(K).coeff(t, k))
    print(f"  t^{k}: {c}")

c7 = sp.expand(K).coeff(t, 7)
target = -sp.Rational(3, 28) / PI**3 * (1 - m2 - m2**2)
print("\nc7 - ( -3/(28 pi^3) * (m4 - m2^2) with m4 = 1 - m2 ) =", sp.simplify(c7 - target))
print("Var_mu(x^2) = m4 - m2^2 =", sp.expand(xm(4) - m2**2))
print("check moments m4, m6, m8:", xm(4), xm(6), xm(8))
