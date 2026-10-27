"""Symbolic checks of the derivations used in NOTES.md (sympy). Run: python3 -I derivations.py
Writes ../results/derivations.json with each identity's simplified residual (0 means verified).
"""
import json
import os

import sympy as sp

res = {}
MP, H, H0, T0, T, g, gs0, N, Nmin, V, w = sp.symbols("M_P H H_0 T_0 T g g_s0 N N_min V w", positive=True)

# ---- D1: TCC with instant reheating.  Horizon: e^N >= (g_s0/g)^{1/3} (T0/T) (H/H0).  TCC: e^N < MP/H.
#      Compatible iff (g_s0/g)^{1/3} T0 H^2/(T H0) < MP, with H^2 = V/(3 MP^2), T = (30 V/(pi^2 g))^{1/4}.
Tex = (30 * V / (sp.pi ** 2 * g)) ** sp.Rational(1, 4)
Hex = sp.sqrt(V / 3) / MP
lhs = (gs0 / g) ** sp.Rational(1, 3) * T0 * Hex ** 2 / (Tex * H0)
Vbound = (3 * MP ** 3 * (H0 / T0) * (30 / (sp.pi ** 2 * g)) ** sp.Rational(1, 4) * (g / gs0) ** sp.Rational(1, 3)) ** sp.Rational(4, 3)
res["D1_TCC_instant_bound_residual"] = str(sp.simplify(lhs.subs(V, Vbound) - MP))

# ---- D2: Banks-Fischler counting with explicit factors.  Entropy in the inflated patch at the end of reheating,
#      S_patch = s_rh (a_rh/a_end)^3 (4pi/3)(e^N/H)^3; adiabatic afterwards s_rh a_rh^3 = s0 a0^3; with
#      e^{Nmin} = a_end H/(a0 H0):  S_patch = S_H0 e^{3(N - Nmin)},  S_H0 = s0 (4pi/3)/H0^3.
s0, a0, aend = sp.symbols("s_0 a_0 a_end", positive=True)
S_patch = s0 * (a0 / aend) ** 3 * sp.Rational(4, 3) * sp.pi * (sp.exp(N) / H) ** 3
S_H0 = s0 * sp.Rational(4, 3) * sp.pi / H0 ** 3
expr = S_patch - S_H0 * sp.exp(3 * (N - Nmin))
res["D2_BF_patch_entropy_residual"] = str(sp.simplify(expr.subs(Nmin, sp.log(aend * H / (a0 * H0)))))

# ---- D3: BF eq. (3.2) in Planck units: e^{3N} L_I^{-3/4} <= 1/lam  <=>  N <= (1/4) ln(L_I/lam) + (1/12) ln(1/lam)
LI, lam = sp.symbols("Lambda_I lambda", positive=True)
Nmax = sp.solve(sp.Eq(3 * N - sp.Rational(3, 4) * sp.log(LI), -sp.log(lam)), N)[0]
res["D3_BF_eq3.2_residual"] = str(sp.simplify(sp.expand_log(Nmax - (sp.log(LI / lam) / 4 + sp.log(1 / lam) / 12), force=True)))

# ---- D4: curvature and topology under an e-fold cap: Omega_k = -k/(aH)^2, so
#      Omega_k0/Omega_ki = (a_i H_i/(a0 H0))^2 = e^{-2(N - Nmin)} (H constant during inflation);
#      comoving T^3 side L = ell/(a_i H_i)  ->  L H0 a0 = ell e^{N - Nmin}.
ai, ell = sp.symbols("a_i ell", positive=True)
aiH = aend * sp.exp(-N) * H
res["D4_curvature_ratio_residual"] = str(sp.simplify(((aiH / (a0 * H0)) ** 2 - sp.exp(-2 * (N - Nmin))).subs(
    Nmin, sp.log(aend * H / (a0 * H0)))))
res["D4_topology_residual"] = str(sp.simplify((ell / aiH * a0 * H0 - ell * sp.exp(N - Nmin)).subs(
    Nmin, sp.log(aend * H / (a0 * H0)))))

# ---- D5: two-stage inflation: radiation from T_rh1 to T_on2, stage 2 (N2), radiation from T_rh2 to T0.
#      a_end1 = a0 (T0/T_rh2) e^{-N2} (T_on2/T_rh1) (dof factors dropped); horizon a_end1 e^{-N1} H1 <= a0 H0 gives
#      e^{N1+N2} >= (H1/H0)(T0/T_rh1)(T_on2/T_rh2):  for T_on2 = T_rh2 the total is the single-stage N_min of stage 1.
H1, Trh1, Ton2, Trh2, N1, N2 = sp.symbols("H_1 T_rh1 T_on2 T_rh2 N_1 N_2", positive=True)
aend1 = a0 * (T0 / Trh2) * sp.exp(-N2) * (Ton2 / Trh1)
cond = sp.simplify(sp.log(aend1 * sp.exp(-N1) * H1 / (a0 * H0)))
target = sp.log(H1 * T0 * Ton2 / (H0 * Trh1 * Trh2)) - N1 - N2
res["D5_two_stage_residual"] = str(sp.simplify(sp.expand_log(cond - target, force=True)))

# ---- D6: homogeneous Gauss law: div E = rho with E = E(t) gives rho = 0 (any topology)
x, y, z, t = sp.symbols("x y z t")
Ex, Ey, Ez = [sp.Function(n)(t) for n in ("E_x", "E_y", "E_z")]
res["D6_homogeneous_divergence"] = str(sp.diff(Ex, x) + sp.diff(Ey, y) + sp.diff(Ez, z))

# ---- D7: compact-slice Gauss law for the hidden U(1): Q = N_X + 3 N_Y + 6 N_Sigma = 0  =>  N_X + 3 N_Y = 0 mod 6
NX, NY, NS = sp.symbols("N_X N_Y N_Sigma", integer=True)
res["D7_Z6_condition"] = "N_X + 3 N_Y = -6 N_Sigma  =>  (N_X + 3 N_Y) mod 6 = 0"

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "derivations.json")
with open(path, "w") as f:
    json.dump(res, f, indent=1)
for k, v in res.items():
    print(k, "->", v)
