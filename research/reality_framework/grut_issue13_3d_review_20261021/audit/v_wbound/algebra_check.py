#!/usr/bin/env python3
"""Verifier v_wbound: independent symbolic checks of the wbound track (sympy).

Differences from wbound/algebra.py: X is kept complex (X, Xbar treated as independent fields), the stress tensor is
built from the Hilbert form T_mn = 2 dL/dg^mn - g_mn L evaluated on Minkowski, Euler's theorem is not used, the
FRW identity is checked symbolically with an arbitrary a(t), and the boosted-kink integral is done in the lab frame.
Writes results/algebra_check.json.
"""
import json
import os

import sympy as sp

lam, g, lX, hD, lD, mD, v = sp.symbols("lambda g lambda_X h_D lambda_D m_D0 v", positive=True)
ph, X, Xb, D = sp.symbols("phi X Xb D")
U = lam / 4 * (ph**2 - v**2) ** 2 + g**2 * ph**2 * X * Xb + lX * (X * Xb) ** 2 \
    + sp.Rational(1, 2) * (mD**2 + hD**2 * ph**2) * D**2 + lD / 4 * D**4
out = {}

# ---- (a) trace from the Hilbert stress tensor, using derivative symbols ----
# kinetic invariants: K_ab = g^mn d_m f_a d_n f_b.  L = -1/2 g^mn dphi dphi - g^mn dX dXb - 1/2 g^mn dD dD - U
# T_mn = d_m phi d_n phi + (d_m X d_n Xb + d_n X d_m Xb) + d_m D d_n D + g_mn L
d = {f: sp.symbols(f"d{n}0:4") for f, n in ((ph, "p"), (X, "x"), (Xb, "xb"), (D, "D"))}
eta = sp.diag(-1, 1, 1, 1)
inv = eta  # Minkowski inverse equals itself
L = -sp.Rational(1, 2) * sum(inv[m, m] * d[ph][m] ** 2 for m in range(4)) \
    - sum(inv[m, m] * d[X][m] * d[Xb][m] for m in range(4)) \
    - sp.Rational(1, 2) * sum(inv[m, m] * d[D][m] ** 2 for m in range(4)) - U
T = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        T[m, n] = d[ph][m] * d[ph][n] + d[X][m] * d[Xb][n] + d[X][n] * d[Xb][m] + d[D][m] * d[D][n] + eta[m, n] * L
rho = T[0, 0]
p = (T[1, 1] + T[2, 2] + T[3, 3]) / 3
kin = d[ph][0] ** 2 + 2 * d[X][0] * d[Xb][0] + d[D][0] ** 2
grad = sum(d[ph][i] ** 2 + 2 * d[X][i] * d[Xb][i] + d[D][i] ** 2 for i in (1, 2, 3))
out["a_residual"] = str(sp.expand(3 * p - rho - (kin - grad - 4 * U)))
out["a_trace_residual"] = str(sp.expand(sum(inv[m, m] * T[m, m] for m in range(4)) - (3 * p - rho)))
out["rho_check"] = str(sp.expand(rho - (sp.Rational(1, 2) * d[ph][0] ** 2 + d[X][0] * d[Xb][0] + sp.Rational(1, 2) * d[D][0] ** 2
                                       + sp.Rational(1, 2) * sum(d[ph][i] ** 2 + 2 * d[X][i] * d[Xb][i] + d[D][i] ** 2 for i in (1, 2, 3)) + U)))


# ---- (b) virial integrands, computed directly (complex X: (X - cX) dU/dX + (Xb - cXb) dU/dXb) ----
def Phi(cphi):
    return sp.expand((ph - cphi) * sp.diff(U, ph) + X * sp.diff(U, X) + Xb * sp.diff(U, Xb) + D * sp.diff(U, D) - 4 * U)


Pv = Phi(v)
P0 = Phi(0)
out["Phi_plus_v_residual"] = str(sp.expand(Pv - (-lam * v * (ph - v) ** 2 * (ph + v) - 2 * g**2 * v * ph * X * Xb
                                                 - D**2 * (mD**2 + hD**2 * v * ph))))
out["Phi_0_residual"] = str(sp.expand(P0 - (lam * v**2 * (ph**2 - v**2) - mD**2 * D**2)))
out["Phi_v_phi_only_factor"] = str(sp.factor(Pv.subs({X: 0, Xb: 0, D: 0})))
# pointwise condition with X for -v <= phi < 0 (D = 0): Phi_v <= 0  <=>  |X|^2 <= lam (phi-v)^2 (phi+v) / (2 g^2 |phi|)
Xa = sp.symbols("Xa", nonnegative=True)  # |X|^2
crit = sp.solve(sp.Eq(Pv.subs({X * Xb: Xa}).subs({X: sp.sqrt(Xa), Xb: sp.sqrt(Xa), D: 0}), 0), Xa)
out["Phi_v_zero_in_X2"] = [str(sp.simplify(c)) for c in crit]

# ---- FRW identity, homogeneous + gradient-free check with arbitrary a(t) (field eqs used explicitly) ----
t = sp.symbols("t")
a = sp.Function("a")(t)
H = sp.diff(a, t) / a
# use one real field with generic potential V(f) and generic gradient-energy density (not needed: spatial means of
# gradient terms are handled by integration by parts; the time-derivative algebra below is the nontrivial part)
f = sp.Function("f")(t)
V = sp.Function("V")
c = sp.symbols("c")
fdd = -3 * H * sp.diff(f, t) - sp.diff(V(f), f)  # homogeneous EOM
gc = (f - c) * sp.diff(f, t)
lhs = a**-3 * sp.diff(a**3 * gc, t)
lhs = lhs.subs(sp.Derivative(f, (t, 2)), fdd)
three_p_minus_rho = sp.diff(f, t) ** 2 - 4 * V(f)
Phi_c = (f - c) * sp.diff(V(f), f) - 4 * V(f)
out["frw_identity_residual_homog"] = str(sp.simplify(sp.expand(three_p_minus_rho - Phi_c - lhs)))
# continuity: d(rho a^4)/dt = -H a^4 (3p - rho)
rho_h = sp.diff(f, t) ** 2 / 2 + V(f)
p_h = sp.diff(f, t) ** 2 / 2 - V(f)
cont = sp.diff(rho_h * a**4, t).subs(sp.Derivative(f, (t, 2)), fdd) + H * a**4 * (3 * p_h - rho_h)
out["frw_continuity_residual_homog"] = str(sp.simplify(sp.expand(cont)))
# with gradients: 1D field f(t,z), EOM f_tt + 3H f_t - f_zz/a^2 + V' = 0; check pointwise identity
z = sp.symbols("z")
F = sp.Function("F")(t, z)
Ftt = -3 * H * sp.diff(F, t) + sp.diff(F, z, 2) / a**2 - sp.diff(V(F), F)
gloc = (F - c) * sp.diff(F, t)
lhs2 = (a**-3 * sp.diff(a**3 * gloc, t)).subs(sp.Derivative(F, (t, 2)), Ftt)
tpr = sp.diff(F, t) ** 2 - sp.diff(F, z) ** 2 / a**2 - 4 * V(F)
Phi2 = (F - c) * sp.diff(V(F), F) - 4 * V(F)
# pointwise: 3p - rho - Phi - a^-3 d(a^3 g)/dt = - d/dz[(F - c) F_z]/a^2  (a total spatial derivative -> zero mean)
resid2 = sp.simplify(sp.expand(tpr - Phi2 - lhs2 + sp.diff((F - c) * sp.diff(F, z), z) / a**2))
out["frw_identity_residual_with_gradient_up_to_divergence"] = str(resid2)

# ---- boosted static kink: integral of (3p - rho) dz in the lab frame ----
zz, tt = sp.symbols("zz tt", real=True)
u = sp.symbols("u", positive=True)
gam = 1 / sp.sqrt(1 - u**2)
m = v * sp.sqrt(lam / 2)
kink = v * sp.tanh(m * gam * (zz - u * tt))
Uk = lam / 4 * (kink**2 - v**2) ** 2
tr = sp.diff(kink, tt) ** 2 - sp.diff(kink, zz) ** 2 - 4 * Uk  # 3p - rho (only z dependence)
en = sp.diff(kink, tt) ** 2 / 2 + sp.diff(kink, zz) ** 2 / 2 + Uk
s = sp.symbols("s", real=True)
tr_s = sp.simplify(tr.subs(tt, 0).subs(zz, s / (m * gam)) / (m * gam))  # dz = ds/(m gamma)
en_s = sp.simplify(en.subs(tt, 0).subs(zz, s / (m * gam)) / (m * gam))
I_tr = sp.simplify(sp.integrate(sp.simplify(tr_s), (s, -sp.oo, sp.oo)))
I_en = sp.simplify(sp.integrate(sp.simplify(en_s), (s, -sp.oo, sp.oo)))
sigma = 2 * sp.sqrt(2 * lam) * v**3 / 3
out["kink_int_3p_minus_rho_over_sigma"] = str(sp.simplify(I_tr / sigma))
out["kink_int_rho_over_sigma"] = str(sp.simplify(I_en / sigma))
w_kink = sp.simplify((I_tr / I_en + 1) / 3)
out["kink_w"] = str(w_kink)
out["kink_w_minus_NG(u^2-2/3)"] = str(sp.simplify(w_kink - (u**2 - sp.Rational(2, 3))))

# ---- Cauchy-Schwarz cap: lam v^2 <phi^2 - v^2> <= lam v^2 sqrt(<(phi^2-v^2)^2>) <= lam v^2 sqrt(4 rho/lam) ----
rb = sp.symbols("rhobar", positive=True)
cap = sp.simplify(lam * v**2 * sp.sqrt(4 * rb / lam) / (3 * rb))
out["cs_cap_wbar_minus_third"] = str(cap)
out["cs_cap_matches_(2/3)sqrt(lam v^4/rho)"] = str(sp.simplify(cap - sp.Rational(2, 3) * sp.sqrt(lam * v**4 / rb)))

# ---- numbers ----
import math
r1 = 7 / 8 * (4 / 11) ** (4 / 3)
out["rho1nu_over_rhogamma"] = r1
out["DNeff_eq_Neff3.044"] = 1 / r1 + 3.044
out["DNeff_eq_Neff3"] = 1 / r1 + 3.0
MP = 2.435323e18
out["boundary_(2/sqrt3)(v/MP)_v1TeV"] = 2 / math.sqrt(3) * 1e3 / MP
out["eps_star_GeV4_if_0.58768"] = 0.58768 * 0.01 * 1e12
out["amp2_energy_lam_v4"] = (4 - 1) ** 2 / 4
out["amp2_energy_v4_at_lam0.01"] = 0.01 * (4 - 1) ** 2 / 4

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "results", "algebra_check.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
