#!/usr/bin/env python3
"""Track wbound, items (a), (b) algebra and (e) number.  Symbolic (sympy), deterministic.

(a) 3p - rho for L = -1/2 (dphi)^2 - |dX|^2 - 1/2 (dD)^2 - U, metric diag(-1,1,1,1), X = (x1 + i x2)/sqrt(2).
(b) Virial integrand Phi_c = sum_a (f_a - c_a) dU/df_a - 4U for origins c = (+v,0,0,0) and c = 0.
(e) Delta N_eff-equivalent of rho_hidden = rho_rad(standard, N_eff = 3.044).
Writes results/algebra.json.
"""
import json
import os

import sympy as sp

t, x, y, z = sp.symbols("t x y z", real=True)
lam, g, lX, hD, lD, mD0, v = sp.symbols("lambda g lambda_X h_D lambda_D m_D0 v", positive=True)
X = [t, x, y, z]
eta = sp.diag(-1, 1, 1, 1)

fphi, fx1, fx2, fD = [sp.Function(n)(t, x, y, z) for n in ("phi", "x1", "x2", "D")]
fields = [fphi, fx1, fx2, fD]
P, X1, X2, Dd = sp.symbols("P X1 X2 Dd", real=True)          # plain symbols for the potential
X2abs = (X1**2 + X2**2) / 2                                     # |X|^2
U = lam / 4 * (P**2 - v**2) ** 2 + g**2 * P**2 * X2abs + lX * X2abs**2 \
    + sp.Rational(1, 2) * (mD0**2 + hD**2 * P**2) * Dd**2 + lD / 4 * Dd**4
subs_f = {P: fphi, X1: fx1, X2: fx2, Dd: fD}
Uf = U.subs(subs_f)

# canonical stress tensor T_mn = sum_a d_m f_a d_n f_a + eta_mn L  (all four real fields canonically normalised)
def dmf(f, m):
    return sp.diff(f, X[m])

L = -sp.Rational(1, 2) * sum(eta[m, m] * dmf(f, m) ** 2 for f in fields for m in range(4)) - Uf
T = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        T[m, n] = sum(dmf(f, m) * dmf(f, n) for f in fields) + eta[m, n] * L
rho = sp.simplify(T[0, 0])
p = sp.simplify((T[1, 1] + T[2, 2] + T[3, 3]) / 3)
kin = sum(dmf(f, 0) ** 2 for f in fields)                   # phidot^2 + 2|Xdot|^2 + Ddot^2
grad2 = sum(dmf(f, m) ** 2 for f in fields for m in (1, 2, 3))
claim_a = kin - grad2 - 4 * Uf
res_a = sp.simplify(sp.expand(3 * p - rho - claim_a))
trace = sp.simplify(sum(eta[m, m] * T[m, m] for m in range(4)))   # T^mu_mu
res_trace = sp.simplify(sp.expand(trace - (3 * p - rho)))

# (b) virial integrands
def Phi(c):
    cP, c1, c2, cD = c
    return sp.expand((P - cP) * sp.diff(U, P) + (X1 - c1) * sp.diff(U, X1) + (X2 - c2) * sp.diff(U, X2)
                     + (Dd - cD) * sp.diff(U, Dd) - 4 * U)

s = sp.symbols("s", real=True)
Phi_v = Phi((v, 0, 0, 0))
claim_v = -lam * v * (P - v) ** 2 * (P + v) - 2 * g**2 * v * P * X2abs - Dd**2 * (mD0**2 + hD**2 * v * P)
res_v = sp.simplify(sp.expand(Phi_v - claim_v))
Phi_v_phi_only = sp.factor(Phi_v.subs({X1: 0, X2: 0, Dd: 0}))
Phi_0 = Phi((0, 0, 0, 0))
claim_0 = lam * v**2 * (P**2 - v**2) - mD0**2 * Dd**2
res_0 = sp.simplify(sp.expand(Phi_0 - claim_0))
# origin independence: Phi_c - Phi_0 = -c . grad U
res_shift = sp.simplify(sp.expand(Phi_v - Phi_0 + v * sp.diff(U, P)))
# degree decomposition about phi = +v: sum (k-4) U_k
U_s = sp.expand(U.subs(P, v + s))
poly = sp.Poly(U_s, s, X1, X2, Dd)
euler = sum((sum(mon) - 4) * coeff * s**mon[0] * X1**mon[1] * X2**mon[2] * Dd**mon[3]
            for mon, coeff in poly.terms())
res_euler = sp.simplify(sp.expand(euler - Phi_v.subs(P, v + s)))
# NG wall check: static tanh kink, integral of lam v^2 (phi^2 - v^2) dz = -3 sigma
zz = sp.symbols("zz", real=True)
kink = v * sp.tanh(v * sp.sqrt(lam / 2) * zz)
I_y = sp.integrate(lam * v**2 * (kink**2 - v**2), (zz, -sp.oo, sp.oo))
sigma = 2 * sp.sqrt(2 * lam) / 3 * v**3
res_kink = sp.simplify(I_y + 3 * sigma)

# (e) Delta N_eff-equivalent
Neff_std = 3.044
r1 = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)        # rho_1nu / rho_gamma
dN_eq = (1 + Neff_std * r1) / r1                       # rho_rad/rho_1nu
dN_eq_3 = (1 + 3.0 * r1) / r1

out = {
    "a_3p_minus_rho_residual": str(res_a),
    "a_3p_minus_rho_equals_trace_residual": str(res_trace),
    "a_statement": "3p - rho = (phidot^2 + 2|Xdot|^2 + Ddot^2) - ((grad phi)^2 + 2|grad X|^2 + (grad D)^2) - 4U = T^mu_mu",
    "b_Phi_plus_v_residual": str(res_v),
    "b_Phi_plus_v": "-lambda v (phi-v)^2 (phi+v) - 2 g^2 v phi |X|^2 - D^2 (m_D0^2 + h_D^2 v phi)",
    "b_Phi_plus_v_phi_only_factored": str(Phi_v_phi_only),
    "b_Phi_0_residual": str(res_0),
    "b_Phi_0": "lambda v^2 (phi^2 - v^2) - m_D0^2 D^2   (X drops out: all X terms are quartic about the origin)",
    "b_origin_shift_residual": str(res_shift),
    "b_euler_degree_decomposition_residual": str(res_euler),
    "b_static_kink_identity_residual (int lam v^2 (phi^2-v^2) dz + 3 sigma)": str(res_kink),
    "e_rho1nu_over_rhogamma": r1,
    "e_DeltaNeff_equivalent_of_rho_h_eq_rho_rad_Neff3.044": dN_eq,
    "e_DeltaNeff_equivalent_of_rho_h_eq_rho_rad_Neff3": dN_eq_3,
}
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "results", "algebra.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
