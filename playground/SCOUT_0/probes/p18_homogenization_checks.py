# P-18 exact checks (sympy; no simulation, per charter). Fast deterministic chaos = doubling map T(y) = 2y mod 1
# (Lebesgue-invariant, mixing, uniformly expanding). Slow: x_{n+1} = x_n + eps^2 f(x_n) + eps h(y_n).
# Homogenization (imported: Melbourne-Stuart 2011; Gottwald-Melbourne 2013; Kelly-Melbourne 2016): x -> dX = f dt + Sigma dW,
# Sigma^2 = C_0 + 2 sum_{n>=1} C_n  (Green-Kubo), C_n = int h(y) h(T^n y) dy, provided int h = 0.
import sympy as sp
y = sp.symbols('y', real=True)
c = lambda m: sp.cos(2*sp.pi*m*y)
def corr(h, n):   # int_0^1 h(y) h(2^n y) dy, exact
    return sp.simplify(sp.integrate(sp.expand(sp.expand_trig(h*h.subs(y, 2**n*y))), (y, 0, 1)))
def green_kubo(h, nmax=6):
    C = [corr(h, n) for n in range(nmax+1)]
    return C, sp.simplify(C[0] + 2*sum(C[1:]))
cases = {
    "h = cos 2pi y            ": c(1),
    "h = cos 2pi y + cos 4pi y": c(1) + c(2),
    "h = chi o T - chi, chi = cos 2pi y (coboundary)": c(2) - c(1),
}
print("=== C1 Green-Kubo noise intensity of a deterministic chaotic forcing (exact) ===")
for name, h in cases.items():
    mean = sp.integrate(h, (y, 0, 1))
    C, S2 = green_kubo(h)
    print(f"  {name}: mean={mean}, C_0..C_3 = {C[:4]}, Sigma^2 = {S2}",
          "-> NO noise survives (coboundary)" if S2 == 0 else "-> Brownian limit, Q = Sigma^2")
print("\n=== C2 third cumulant: present at finite eps, O(eps) in the slow increment, gone in the limit ===")
h = c(1) + c(2)
k3 = sp.simplify(sp.integrate(sp.expand(h**3), (y, 0, 1)))
print(f"  E_mu[h^3] = {k3} (skewed forcing). Over slow time t: kappa3(eps S_N) ~ eps^3 * N * kbar3 = eps*t*kbar3 -> 0")
print("  O(beta) lemma: even-in-a term -4 beta int g kappa3 = O(eps); odd-in-a term covariance-only -> C-B with Q = Sigma^2")
print("\n=== C3 short-time order at finite eps (smooth forcing, first-order slow system) ===")
t, s, sig2, b, a = sp.symbols('t s sigma2 beta a', positive=True)
lead = sp.integrate(-12*b*a*sig2*s**2, (s, 0, t))
print(f"  v(s) ~ sigma_F^2 s^2 for s << correlation time  =>  S - D ~ {lead}  (t^3, vs C-B t^2 beyond the window)")
