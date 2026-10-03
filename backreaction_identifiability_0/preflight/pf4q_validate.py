"""BRI1-PF4Q validations V1-V3 (BRI1_PF4Q.md sec. 3).  Prints NO frozen-tau coefficient.  Evidence-grade except V1
(arb-certified enclosures).  No Monte Carlo."""
import sys, numpy as np, sympy as sp
from flint import arb, acb, ctx
import pf4q_core as C

ctx.prec = 200
# ---------------- V1: certified Gibbs moments (arb)
L = 8
def mom(k):
    f = lambda x, analytic: x**k * (-(x * x) / 2 - x**4 / 4).exp()
    core = 2 * acb.integral(f, 0, L).real
    tail = arb(0, rad=arb(10) ** (-400))          # int_L^inf x^k e^{-x^4/4} <= e^{-L^4/4 + ...}, L=8: < 1e-400
    return core + tail
Z, M2, M4 = mom(0), mom(2), mom(4)
m2, m4 = M2 / Z, M4 / Z
print("V1 certified: Z_x =", Z.str(20, radius=True))
print("V1 certified: m2 =", m2.str(20, radius=True), " m4 =", m4.str(20, radius=True))
s = m2 + m4; var = 1 - m2 - m2 * m2
print("V1 m2+m4 =", s.str(20, radius=True), " contains 1:", bool(s.contains(1)))
print("V1 Var(x0^2) = 1-m2-m2^2 =", var.str(20, radius=True), " excludes 0 and >0:", bool(var > 0))
m2f = float(m2.mid())
for n in (48, 64, 96):
    X, P, W, pr = C.rule_gh(n)
    print(f"V1 quadrature GH n={n}: m2 {W @ X**2:.15f}  m4 {W @ X**4:.15f}  (|dm2| {abs(W @ X**2 - m2f):.1e}; pruned mass {pr:.1e}; nodes {len(W)})")
for h in (0.12, 0.08):
    X, P, W = C.rule_trap(h)
    print(f"V1 quadrature trap h={h}: m2 {W @ X**2:.15f} (|dm2| {abs(W @ X**2 - m2f):.1e}; nodes {len(W)})")
sys.stdout.flush()
# ---------------- V2: P0 one-time variance time independence (Method A nodes, n=64)
X, P, W, _ = C.rule_gh(64)
sol = C.variational(X, P, [C.T_V] + list(C.TAU))
for t in [C.T_V] + list(C.TAU):
    print(f"V2 P0: E[x0(t)^2] at t={t:.6f}: {W @ sol[t][0]**2:.15f}  (minus m2: {W @ sol[t][0]**2 - m2f:.1e})")
# ---------------- V3: small validation time t_v = pi/8 vs exact series (P1)
mm = sp.symbols('m2')
coef = {7: (mm**2 + mm - 1) / (28 * sp.pi**3), 8: -3 * (mm**2 + mm - 1) / (112 * sp.pi**4),
        9: (sp.pi**2 * mm**2 + 12 * mm**2 + 12 * mm + 19 * sp.pi**2 * mm - 12 - sp.pi**2) / (2016 * sp.pi**5),
        10: -(mm**2 + 19 * mm - 1) / (3360 * sp.pi**4),
        11: (2 * mm**2 + 11 * sp.pi**2 * mm**2 + 38 * mm + 38 * sp.pi**2 * mm - 60 * sp.pi**2 - 2) / (36960 * sp.pi**5),
        12: -(11 * mm**2 + 38 * mm - 60) / (73920 * sp.pi**4),
        13: (792 * mm**2 + 617 * sp.pi**2 * mm**2 + 2736 * mm + 19607 * sp.pi**2 * mm - 725 * sp.pi**2 - 4320) / (34594560 * sp.pi**5),
        14: -(617 * mm**2 + 19607 * mm - 725) / (80720640 * sp.pi**4)}
tv = C.T_V
series = float(sum(v.subs(mm, m2f) * sp.pi**0 * (sp.pi / 8)**k for k, v in coef.items()))
last = float(coef[14].subs(mm, m2f) * (sp.pi / 8)**14)
xa = sol[tv][0]; ya = sol[tv][1]
cA = C.c_term(W, xa, xa, ya)
Xb, Pb, Wb = C.rule_trap(0.08)
cB = {}
for eps in (1e-3, 5e-4):
    s1 = C.nonlinear(Xb, Pb, eps, 1, [tv])
    x = s1[tv]; cB[eps] = C.k3(Wb, x, x, x) / eps / 3          # kappa3/eps = K = 3c at one time
cBR = (4 * cB[5e-4] - cB[1e-3]) / 3
print(f"V3 t_v=pi/8 (validation only): series c = {series:.6e} (last term {last:.1e}); method A c = {cA:.6e}; "
      f"method B c (Richardson) = {cBR:.6e}; sign negative: {series < 0 and cA < 0 and cBR < 0}")
