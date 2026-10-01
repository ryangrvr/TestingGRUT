# P-23a / P-23 checks (requires sympy, numpy). See P23A_RESULT.md, P23_RESULT.md.
import sympy as sp, numpy as np

print("=== A1 (P-23a) the variance relation is the generic small-noise expansion, for ANY drift f ===")
x, a, T = sp.symbols('x a T'); f = sp.Function('f')
L = lambda g, Q: sp.expand(f(x)*sp.diff(g, x) + Q/2*sp.diff(g, x, 2))
d2 = sp.simplify((L(L(x, 2*T), 2*T) - L(L(x, 0), 0)).subs(x, a))
print("  Delta c2 for arbitrary drift f:", d2, " -> O(noise) mean shift = (1/2) int g f''(mu) v ds")
beta, k = sp.symbols('beta k')
fq = -k*x - 4*beta*x**3
print("  quartic C-B drift: T f''(a) =", sp.expand(T*sp.diff(fq, x, 2).subs(x, a)), "(record: -24 beta T a)")
fd = x - x**3                       # double well (Duffing-type), non-GRUT
print("  double-well x - x^3 (non-GRUT): T f''(a) =", sp.expand(T*sp.diff(fd, x, 2).subs(x, a)),
      " -> same relation form, (S-D) = -3 int g mu v")

print("\n=== B1 (P-23) the one law-level lock on record: interior family mu=1+x*al, eta=1/(1+x*al), Sigma=(mu(1+eta))/2 ===")
xx, al = sp.symbols('x alpha', positive=True)
mu = 1 + xx*al; eta = 1/(1 + xx*al); Sig = sp.simplify(mu*(1 + eta)/2)
print("  Sigma =", Sig, ";  Sigma0 - mu0/2 =", sp.simplify((Sig - 1) - (mu - 1)/2), " (x and alpha both eliminated)")
print("  mu*eta =", sp.simplify(mu*eta), " (curvature potential Phi unmodified; only Psi modified)")
# conformal scalar-tensor / f(R) quasi-static small-scale limit for comparison
mufR, etafR = sp.Rational(4, 3), sp.Rational(1, 2)
print("  f(R) small-scale: mu=4/3, eta=1/2 -> Sigma =", mufR*(1 + etafR)/2, "(Sigma0 = 0; off the GRUT line Sigma0 = mu0/2)")

print("\n=== B2 (P-23) CM of the retained response: symmetric K (reversible) vs cycle affinity (record E-4) ===")
Ksym = np.array([[2.3, -1, 0], [-1, 2.3, -1], [0, -1, 1.3]])
w = np.linalg.eigvals(Ksym)
print("  symmetric K_b-type block: eigenvalues real:", np.allclose(np.imag(w), 0), "-> k(t) = sum u^2 e^{-lam t} is CM (Debye/reversible: standard)")
eps = 0.8                            # driven ring: asymmetric hopping = nonzero cycle affinity
Kring = np.array([[2, -(1+eps), -(1-eps)], [-(1-eps), 2, -(1+eps)], [-(1+eps), -(1-eps), 2]])
w2 = np.linalg.eigvals(Kring)
print("  3-ring with affinity: eigenvalues", np.round(w2, 4), "-> complex pair, k(t) oscillates, not CM (standard: broken detailed balance)")

print("\n=== B3 (P-23 charter) decoherence-rate mass exponent moves with supplied geometry / correlation length ===")
# N point masses m0 on a cubic lattice (spacing 1), saturated regime: Gamma ~ sum_ij m0^2 C(r_i - r_j), C(r) = exp(-r^2/(4 l^2))
def gamma(n, l):
    g = np.arange(n); P = np.array(np.meshgrid(g, g, g)).reshape(3, -1).T.astype(float)
    D2 = ((P[:, None, :] - P[None, :, :])**2).sum(-1)
    return np.exp(-D2/(4*l*l)).sum()
for l in (50.0, 2.0, 0.3):
    ns = [2, 4, 8]; G = [gamma(n, l) for n in ns]; m = [n**3 for n in ns]
    ex = [np.log(G[i+1]/G[i])/np.log(m[i+1]/m[i]) for i in range(2)]
    print(f"  correlation length l={l:5}: d ln Gamma / d ln m over m=8->64->512 :", ", ".join(f"{e:.3f}" for e in ex))
print("  -> exponent 2 (coherent, R << l) ... 1 (incoherent, R >> l): set by supplied R/l, shared by DP/CSL/any mass-density coupling")
