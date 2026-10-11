#!/usr/bin/env python3
"""
V2 -- The Einstein equation as the equilibrium equation of state of the responsive medium (Jacobson 1995),
re-tiered: its two inputs are (i) KMS at beta = 2pi/kappa for the wedge (Bisognano-Wichmann: a THEOREM of the
covariant vacuum GRUT already pays for) and (ii) horizon entropy proportional to area with coefficient
eta = 1/(4 G hbar) -- area-proportionality is a QFT theorem (entanglement area law); the coefficient is the
DATUM G.  So GR costs GRUT III one datum, not the record's '+2 assumed'.

Checks:
 2.1  Raychaudhuri + Clausius bookkeeping (sympy): with theta(0)=sigma=0 at the point,
        T_kk = (hbar eta / 2 pi) R_kk  and with eta = 1/(4 G hbar):  T_kk = R_kk/(8 pi G).
 2.2  Tensor lemma (numerics): X_ab k^a k^b = 0 for all null k  =>  X_ab = f g_ab   (nullspace is 1-dim, = metric).
 2.3  Bianchi closure (sympy, FRW example): R_ab - 8 pi G T_ab = f g_ab with conserved T forces f = R/2 + Lambda,
      i.e. G_ab + Lambda g_ab = 8 pi G T_ab with Lambda an INTEGRATION CONSTANT (the medium does not fix it).
"""
import numpy as np
import sympy as sp
ok_all = True

print("== 2.1 Clausius bookkeeping ==")
hbar, eta, kappa, Lam, Rkk, Tkk, G = sp.symbols('hbar eta kappa Lambda R_kk T_kk G', positive=True)
lam = sp.symbols('lambda', positive=True)
theta = -Rkk*lam                       # Raychaudhuri, linearised about theta(0)=0, sigma=0 (no theta^2 term at this order)
dA_per_area = sp.integrate(theta, (lam, 0, Lam))            # fractional area change of the horizon patch
chi = -kappa*lam                        # boost Killing field along the generator: chi^a = -kappa lambda k^a
dQ_per_area = sp.integrate(Tkk*chi, (lam, 0, Lam))           # heat flux through the patch, per unit area
T_unruh = hbar*kappa/(2*sp.pi)
dS_per_area = eta*dA_per_area
clausius = sp.Eq(dQ_per_area, T_unruh*dS_per_area)
sol = sp.solve(clausius, Tkk)[0]
print("  Clausius dQ = T dS  =>  T_kk =", sp.simplify(sol))
want = hbar*eta*Rkk/(2*sp.pi)
ok = sp.simplify(sol - want) == 0
print("   matches (hbar eta/2pi) R_kk:", "PASS" if ok else "FAIL"); ok_all &= ok
sol_G = sp.simplify(sol.subs(eta, 1/(4*G*hbar)))
print("  with eta = 1/(4 G hbar):  T_kk =", sol_G, " ->", "PASS" if sp.simplify(sol_G - Rkk/(8*sp.pi*G)) == 0 else "FAIL")
ok_all &= sp.simplify(sol_G - Rkk/(8*sp.pi*G)) == 0

print("== 2.2 tensor lemma: null-cone data determine X up to the metric ==")
rng = np.random.default_rng(7)
g = np.diag([-1., 1., 1., 1.])
rows = []
for _ in range(60):
    v = rng.normal(size=3)
    k = np.concatenate([[np.linalg.norm(v)], v])      # null: -k0^2 + |v|^2 = 0
    # unknowns: 10 independent X_ab (symmetric); row = coefficients of X_ab in k^a k^b X_ab
    row = []
    for a in range(4):
        for b_ in range(a, 4):
            row.append(k[a]*k[b_]*(1 if a == b_ else 2))
    rows.append(row)
M = np.array(rows)
U, s, Vt = np.linalg.svd(M)
print("  singular values:", np.array2string(s, precision=3))
null = Vt[-1]
# reconstruct X from the 10-vector
X = np.zeros((4, 4)); idx = 0
for a in range(4):
    for b_ in range(a, 4):
        X[a, b_] = X[b_, a] = null[idx]; idx += 1
X = X/X[1, 1]
ok = (s[-1] < 1e-10*s[0]) and (s[-2] > 1e-3*s[0]) and np.allclose(X, g, atol=1e-9)
print("  one-dimensional nullspace, spanned by:\n", np.array2string(X, precision=6, suppress_small=True))
print("   = eta_ab:", "PASS" if ok else "FAIL"); ok_all &= ok

print("== 2.3 Bianchi closure on a spatially flat FRW example ==")
tt = sp.symbols('t', positive=True)
a = sp.Function('a')(tt)
rho, p = sp.Function('rho')(tt), sp.Function('p')(tt)
Hh = sp.diff(a, tt)/a
# Ricci tensor components (mostly plus): R_tt = -3 addot/a ; R_ij = (a addot + 2 adot^2) delta_ij ; R = 6(addot/a + adot^2/a^2)
addot = sp.diff(a, tt, 2); adot = sp.diff(a, tt)
R_tt = -3*addot/a
R_xx = a*addot + 2*adot**2
Rs = 6*(addot/a + adot**2/a**2)
# perfect fluid, g_tt=-1, g_xx=a^2:  T_tt = rho, T_xx = p a^2
Gn = sp.symbols('G', positive=True)
f = sp.Function('f')(tt)
# demand R_ab - 8 pi G T_ab = f g_ab  for both components
eq1 = sp.Eq(R_tt - 8*sp.pi*Gn*rho, f*(-1))
eq2 = sp.Eq(R_xx - 8*sp.pi*Gn*p*a**2, f*a**2)
f1 = sp.solve(eq1, f)[0]; f2 = sp.solve(eq2, f)[0]
cons = sp.Eq(sp.diff(rho, tt), -3*Hh*(rho + p))          # nabla_a T^ab = 0
# consistency f1 == f2 gives the Friedmann-type relation; check d/dt (f - R/2) = 0 using conservation
expr = sp.simplify(sp.diff(f1 - Rs/2, tt).subs(sp.diff(rho, tt), -3*Hh*(rho + p)).subs(p, sp.solve(sp.Eq(f1, f2), p)[0]))
print("  d/dt [ f - R/2 ] =", sp.simplify(expr), " -> f = R/2 + const:", "PASS" if sp.simplify(expr) == 0 else "FAIL")
ok_all &= sp.simplify(expr) == 0
print("  => G_ab + Lambda g_ab = 8 pi G T_ab, Lambda = integration constant (not a response of the medium).")
print("\nV2 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
