# P-17 C4: inside the white-noise C-B class (1 retained site, K11 = 23/10), which noise statistics reach the
# retained-mean Taylor coefficients?  Generator for additive white noise with cumulant rates Q (=2T), k3, k4
# (Kramers-Moyal / compound-Poisson + Gaussian): L = f d/dx + (Q/2) d^2 + (k3/6) d^3 + (k4/24) d^4.
# Exact rational arithmetic (sympy).
import sympy as sp
x, a, beta, T, k3, k4 = sp.symbols('x a beta T kappa3 kappa4')
K = sp.Rational(23, 10)
f = -K*x - 4*beta*x**3
def L(g, Q, c3, c4):
    return sp.expand(f*sp.diff(g, x) + Q/2*sp.diff(g, x, 2) + c3/6*sp.diff(g, x, 3) + c4/24*sp.diff(g, x, 4))
def coeffs(Q, c3, c4, n=4):
    out, g = [], x
    for _ in range(n+1):
        out.append(sp.expand(g.subs(x, a))); g = L(g, Q, c3, c4)
    return out
det = coeffs(0, 0, 0); gau = coeffs(2*T, 0, 0); gen = coeffs(2*T, k3, k4)
for n in (2, 3, 4):
    print(f"  Delta c_{n}  Gaussian white : {sp.factor(gau[n]-det[n])}")
    print(f"  Delta c_{n}  + k3, k4 extra  : {sp.factor(sp.expand(gen[n]-gau[n]))}")
print("  record check: Delta c2 = -24 beta T a :", sp.simplify(gau[2]-det[2] + 24*beta*T*a) == 0,
      "; Delta c3 = 24 beta T a (44 beta a^2 + 5 K11):", sp.simplify(gau[3]-det[3] - 24*beta*T*a*(44*beta*a**2 + 5*K)) == 0)
