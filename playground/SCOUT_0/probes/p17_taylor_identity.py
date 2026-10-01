# P-17 C2: exact Taylor coefficients of the retained mean E q(t) at t=0 for
#   A  = autonomous Caldeira-Leggett Hamiltonian (finite harmonic bath), bath IC thermal (conditional on q(0)),
#   B  = exogenous Gaussian forcing with the SAME covariance T*gamma(t-s), same memory kernel (GLE),
#   D  = same GLE with no fluctuating force (T = 0 bath: friction kept, forcing absent),
#   R  = autonomous bath, random-phase / fixed-energy preparation (same covariance, non-Gaussian).
# Exact rational arithmetic (sympy). Requires sympy.
import sympy as sp
from functools import lru_cache
from itertools import product

import os
NMAX = int(os.environ.get('P17_NMAX', 8))                                  # Taylor order in t
beta, T, a = sp.symbols('beta T a')
M, k = sp.Rational(1), sp.Rational(23, 10)   # system mass, retained stiffness K11 (record)
modes = [(sp.Rational(1), sp.Rational(1, 2)), (sp.Rational(2), sp.Rational(3, 4))][:int(os.environ.get('P17_MODES', 2))]   # (omega_j, c_j)
J = len(modes)
q, p = sp.symbols('q p')
xs = sp.symbols('x0:%d' % J); ps = sp.symbols('p0:%d' % J)
V = k*q**2/2 + beta*q**4
H = p**2/(2*M) + V + sum(ps[j]**2/2 + w**2*(xs[j] - c*q/w**2)**2/2 for j, (w, c) in enumerate(modes))
Z = [q, p] + list(xs) + list(ps)
flow = {q: sp.diff(H, p), p: -sp.diff(H, q)}
for j in range(J):
    flow[xs[j]] = sp.diff(H, ps[j]); flow[ps[j]] = -sp.diff(H, xs[j])
def lie(g): return sp.expand(sum(sp.diff(g, z)*flow[z] for z in Z))

# --- A: Hamiltonian Taylor coefficients q^(n)(0) as polynomials in the bath ICs
xi = sp.symbols('xi0:%d' % J); eta = sp.symbols('eta0:%d' % J)
ic = {q: a, p: 0}
for j, (w, c) in enumerate(modes):
    ic[xs[j]] = c*a/w**2 + xi[j]; ic[ps[j]] = eta[j]
derivA, g = [], q
for n in range(NMAX+1):
    derivA.append(sp.expand(g.subs(ic))); g = lie(g)

def gauss_mono(e_xi, e_eta):   # independent N(0,T/w^2), N(0,T)
    out = 1
    for j, (w, c) in enumerate(modes):
        for e, var in ((e_xi[j], T/w**2), (e_eta[j], T)):
            if e % 2: return 0
            out *= sp.factorial2(e-1)*var**(e//2) if e else 1
    return out
def phase_mono(e_xi, e_eta):   # xi = sqrt(2T)/w cos phi, eta = -sqrt(2T) sin phi, phi uniform
    out = 1
    for j, (w, c) in enumerate(modes):
        i, l = e_xi[j], e_eta[j]
        if i % 2 or l % 2: return 0
        m, n = i//2, l//2
        mom = sp.factorial(2*m)*sp.factorial(2*n)/(4**(m+n)*sp.factorial(m)*sp.factorial(n)*sp.factorial(m+n))
        out *= (2*T)**(m+n)*w**(-i)*(-1)**l*mom
    return out
def expect(poly, mono):
    P = sp.Poly(poly, *xi, *eta)
    return sp.expand(sum(coef*mono(mon[:J], mon[J:]) for mon, coef in P.terms()))

# --- B / D: GLE  M q'' = -V'(q) - int_0^t gamma(t-s) q'(s) ds + F(t), solved order by order
gam = lambda m: (0 if m % 2 else sum(c**2/w**2*(-1)**(m//2)*w**m for w, c in modes))
f = sp.symbols('f0:%d' % (NMAX+1))
def gle_derivs(fvals):
    d = [a, sp.Integer(0)]
    tt = sp.symbols('tt')
    for n in range(NMAX-1):            # determine d[n+2]
        qser = sum(d[i]*tt**i/sp.factorial(i) for i in range(len(d)))
        Vp = sp.diff(V, q).subs(q, qser)
        dnVp = sp.diff(Vp, tt, n).subs(tt, 0)
        mem = sum(gam(m)*d[n-m] for m in range(n)) if n >= 1 else 0
        d.append(sp.expand((-dnVp - mem + fvals[n])/M))
    return d
derivB = gle_derivs(f)
derivD = gle_derivs([0]*(NMAX+1))
covf = lambda m, l: T*(-1)**l*gam(m+l)       # E F^(m)(0) F^(l)(0) = T (-1)^l gamma^(m+l)(0)
@lru_cache(None)
def wick(idx):                               # E[f_i1 ... f_in] for zero-mean Gaussian, Isserlis
    if not idx: return sp.Integer(1)
    if len(idx) % 2: return sp.Integer(0)
    i0, rest = idx[0], idx[1:]
    return sp.expand(sum(covf(i0, rest[r])*wick(rest[:r]+rest[r+1:]) for r in range(len(rest))))
def expectB(poly):
    P = sp.Poly(poly, *f)
    tot = 0
    for mon, coef in P.terms():
        idx = tuple(i for i, e in enumerate(mon) for _ in range(e))
        tot += coef*wick(idx)
    return sp.expand(tot)

print("=== C2: E q^(n)(0), n = 0..%d ;  A (CL thermal) vs B (exogenous Gaussian, matched cov) ===" % NMAX)
EA = [expect(d, gauss_mono) for d in derivA]
EB = [expectB(d) for d in derivB]
ER = [expect(d, phase_mono) for d in derivA]
ED = [sp.expand(d) for d in derivD]
for n in range(NMAX+1):
    print(f"  n={n}:  A-B = {sp.simplify(EA[n]-EB[n])}   A-D = {sp.factor(EA[n]-ED[n])}   A-R = {sp.factor(EA[n]-ER[n])}")
g0 = gam(0)
print("  predicted leading S-D: t^6 coefficient -beta*a*T*gamma(0)/(10 M^3) -> 6! * that =",
      sp.factor(-sp.factorial(6)*beta*a*T*g0/(10*M**3)))
# controls
print("  control beta=0: A-D all zero:", all(sp.simplify((EA[n]-ED[n]).subs(beta, 0)) == 0 for n in range(NMAX+1)))
print("  control T=0 (bath at rest): A-D all zero:", all(sp.simplify((EA[n]-ED[n]).subs(T, 0)) == 0 for n in range(NMAX+1)))
# product (unshifted) preparation: bath Gibbs around x_j = 0 regardless of q(0) -> initial slip
ic2 = dict(ic)
for j, (w, c) in enumerate(modes): ic2[xs[j]] = xi[j]
g, derivP = q, []
for n in range(NMAX+1):
    derivP.append(sp.expand(g.subs(ic2))); g = lie(g)
EP = [expect(d, gauss_mono) for d in derivP]
slipF = [sp.expand(-a*sum(c**2/w**2*sp.diff(sp.cos(w*sp.Symbol('t')), sp.Symbol('t'), n).subs(sp.Symbol('t'), 0) for w, c in modes)) for n in range(NMAX+1)]
EBslip = [expectB(d) for d in gle_derivs([f[n] + slipF[n] for n in range(NMAX+1)])]
print("  product prep: A_prod - B(shift-matched) nonzero at n:", [n for n in range(NMAX+1) if sp.simplify(EP[n]-EB[n]) != 0],
      "; A_prod - B(+deterministic slip -gamma(t)q0) all zero:", all(sp.simplify(EP[n]-EBslip[n]) == 0 for n in range(NMAX+1)))

# --- B_R: EXOGENOUS random-phase forcing F(t) = sum_j A_j cos(w_j t + theta_j), A_j = c_j sqrt(2T)/w_j,
#     theta_j uniform (no bath, no bath variables; same law as R's free force) -> must equal R exactly
C = sp.symbols('C0:%d' % J); S = sp.symbols('S0:%d' % J)       # cos(theta_j), sin(theta_j)
fR = [sum(c*sp.sqrt(2*T)/w * w**n * (C[j]*sp.cos(sp.pi*n/2) - S[j]*sp.sin(sp.pi*n/2))
          for j, (w, c) in enumerate(modes)) for n in range(NMAX+1)]
def expect_theta(poly):
    P = sp.Poly(sp.expand(poly), *C, *S)
    tot = 0
    for mon, coef in P.terms():
        val = 1
        for j in range(J):
            i, l = mon[j], mon[J+j]
            if i % 2 or l % 2: val = 0; break
            m, n2 = i//2, l//2
            val *= sp.factorial(2*m)*sp.factorial(2*n2)/(4**(m+n2)*sp.factorial(m)*sp.factorial(n2)*sp.factorial(m+n2))
        tot += coef*val
    return sp.expand(tot)
EBR = [expect_theta(d) for d in gle_derivs(fR)]
print("  R (autonomous, random-phase) - B_R (exogenous random-phase, same law): all zero:",
      all(sp.simplify(ER[n]-EBR[n]) == 0 for n in range(NMAX+1)))
print("  B_R - B (exogenous, same covariance, different law): nonzero at n:",
      [n for n in range(NMAX+1) if sp.simplify(EBR[n]-EB[n]) != 0])
