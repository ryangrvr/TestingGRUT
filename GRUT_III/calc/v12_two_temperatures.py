#!/usr/bin/env python3
"""
V12 -- Two first-order temperatures in slowly evolving FRW, and which one GRUT III's principle calls a temperature.

Question (reviewer, 2026-10-11): the comoving conformal-vacuum detector reads T_det = (a''/a')/2pi = (H/2pi)(1 - eps) (v11),
while the cosmological apparent horizon carries the Hayward-Kodama surface gravity kappa = -(1/R_A)(1 - R_A_dot/(2 H R_A))
(Hayward 1998; used with the full kappa by Akbar-Cai 2007), i.e. T_KH = |kappa|/2pi = (H/2pi)(1 - eps/2) in flat FRW.  ATTRIBUTION
(sixth-round audit): Cai-Kim 2005 and Cai-Cao-Hu 2009 use T0 = 1/(2 pi R_A) = H/2pi, NOT T_KH (Cai-Cao-Hu leave the |kappa|/2pi
relation open), and Tian-Booth 2014 compare the two and adopt the Cai-Kim T0; so three first-order candidates exist: T_det, T_KH, T0.  They agree in de Sitter and
differ at O(eps).  GRUT III's principle (P1: the bath is the causal complement, its temperature is the KMS parameter of the
restricted state's modular flow where that flow is geometric; P2: local equilibrium at LOCAL Rindler horizons) must say
which is the bath temperature.

Checks:
 12.1  (sympy) flat FRW: R_A = 1/H, R_A_dot = eps  =>  T_KH = (H/2pi)(1 - eps/2) exactly; T_det = (a''/a')/2pi = (H/2pi)(1 - eps);
       T_det/T_KH = 1 - eps/2 + O(eps^2); equal iff eps = 0 (de Sitter).
 12.2  (sympy) The Kodama vector K = d_t - H r d_r is null on the apparent horizon R_A = 1/H and its Lie derivative of the
       FRW 2-metric has the single component (L_K h)_tr = -a^2 r Hdot: K is Killing iff Hdot = 0 iff eps = 0.  So the apparent
       horizon is a Killing horizon only in de Sitter; for eps != 0 the mismatch T_det/T_KH - 1 = -eps/2 = -R_A_dot/2 is HALF the
       horizon's expansion rate (relative to T0 = 1/(2 pi R_A) it is the full rate, T_det/T0 - 1 = -eps exactly).  The step from
       'not Killing' to 'P1 assigns no temperature' is an application of P1's definition (geometric modular flow), stated in 09, not
       a computation here.
 12.3  P2 is used in D3 only at local Rindler horizons (Jacobson), whose boost kappa cancels; the Friedmann equations follow
       from D3 without any cosmological temperature.  Check: the Cai-Kim relation (flux = T0 dS with T0 = 1/(2 pi R_A)) and
       the Akbar-Cai unified first law (dE = (kappa/2pi) dS + W dV with the SIGNED kappa) are each identities of the Friedmann
       equations, not independent inputs -- verified symbolically for flat FRW with a perfect fluid.
 12.4  Lemma for the next computation (O5 item 1, other fields): a massless minimally coupled field in power-law FRW,
       a ∝ (-eta)^(1/2 - nu), obeys v'' + (k^2 - (nu^2 - 1/4)/eta^2) v = 0 -- the SAME mode equation as a field of mass m
       in de Sitter with nu^2 = 9/4 - m^2/H^2 (mode equations in CONFORMAL time only; worldline functions need the power-law
       t(eta) -- done in v14).  (sympy: a''/a = (nu^2 - 1/4)/eta^2 for a ∝ (-eta)^(1/2 - nu); the exponent -p/(p-1) for a ∝ t^p is
       derived from eta = int dt/a, not typed in.)
Conclusion printed at the end.
"""
import sympy as sp
ok_all = True
print("== 12.1 the two temperatures ==")
H, eps, t = sp.symbols('H epsilon t', positive=True)
R_A = 1/H                               # flat FRW apparent horizon radius
# R_A_dot = -Hdot/H^2 = eps
R_A_dot = eps
kappa = -(1/R_A)*(1 - R_A_dot/(2*H*R_A))
T_KH = sp.simplify(-kappa/(2*sp.pi))
T_det = sp.simplify(H*(1 - eps)/(2*sp.pi))          # (a''/a')/2pi = (H + Hdot/H)/2pi
print("  T_KH  =", T_KH, "   T_det =", T_det)
ok = sp.simplify(T_KH - H*(1-eps/2)/(2*sp.pi)) == 0 and sp.simplify(T_det - H*(1-eps)/(2*sp.pi)) == 0
print("  T_KH = (H/2pi)(1 - eps/2), T_det = (H/2pi)(1 - eps):", "PASS" if ok else "FAIL"); ok_all &= ok
ratio = sp.series(T_det/T_KH, eps, 0, 3).removeO()
print("  T_det/T_KH =", ratio, " -> equal iff eps = 0:", "PASS" if sp.simplify(ratio.subs(eps, 0) - 1) == 0 and ratio.coeff(eps, 1) == -sp.Rational(1, 2) else "FAIL")
ok_all &= (sp.simplify(ratio.subs(eps, 0) - 1) == 0 and ratio.coeff(eps, 1) == -sp.Rational(1, 2))
# a''/a' check from a general a(t): a''/a' = H + Hdot/H = H(1 - eps)
a = sp.Function('a')(t)
Hs = sp.diff(a, t)/a
expr = sp.simplify(sp.diff(a, t, 2)/sp.diff(a, t) - (Hs + sp.diff(Hs, t)/Hs))
print("  a''/a' - (H + Hdot/H) =", expr, "->", "PASS" if expr == 0 else "FAIL"); ok_all &= (expr == 0)

print("== 12.2 Killing character of the apparent horizon (Kodama vector; sympy) ==")
r = sp.symbols('r', positive=True)
coords = [t, r]; hmet = sp.Matrix([[-1, 0], [0, a**2]]); Kv = [sp.Integer(1), -Hs*r]      # K = d_t - H r d_r on the (t, r) 2-metric
LK = sp.Matrix(2, 2, lambda i, j: sp.simplify(sum(Kv[c]*sp.diff(hmet[i, j], coords[c]) + hmet[c, j]*sp.diff(Kv[c], coords[i]) + hmet[i, c]*sp.diff(Kv[c], coords[j]) for c in range(2))))
Knorm = sp.simplify(sum(hmet[i, j]*Kv[i]*Kv[j] for i in range(2) for j in range(2)))
print("  K.K =", Knorm, "-> null on a r = 1/H (the apparent horizon):", "PASS" if sp.simplify(Knorm.subs(r, 1/(a*Hs))) == 0 else "FAIL")
ok_all &= (sp.simplify(Knorm.subs(r, 1/(a*Hs))) == 0)
Hdot_expr = sp.simplify(sp.diff(Hs, t))
print("  (L_K h)_tt, (L_K h)_rr =", LK[0, 0], ",", LK[1, 1], ";  (L_K h)_tr =", sp.factor(LK[0, 1]))
ok = (LK[0, 0] == 0 and LK[1, 1] == 0 and sp.simplify(LK[0, 1] + a**2*r*Hdot_expr) == 0)
print("  (L_K h)_tr = -a^2 r Hdot: K is Killing iff Hdot = 0 iff eps = 0:", "PASS" if ok else "FAIL"); ok_all &= ok
print("  => the apparent horizon is a Killing horizon only in de Sitter; for eps != 0 it is a trapping horizon with a non-Killing")
print("     Kodama vector. P1 (01) assigns a KMS temperature only to regions whose modular flow is geometric; the inference that")
print("     none exists here is P1's definition applied (09), not computed.")

print("== 12.3 the apparent-horizon first law is an identity of the Einstein equations (flat FRW, perfect fluid) ==")
G = sp.symbols('G', positive=True)
rho, p, Hf = sp.Function('rho')(t), sp.Function('p')(t), sp.Function('H')(t)
Hdot = -4*sp.pi*G*(rho + p); rhodot = -3*Hf*(rho + p)
def on_shell(e):
    return sp.simplify(e.subs(sp.Derivative(Hf, t), Hdot).subs(sp.Derivative(rho, t), rhodot).subs(rho, 3*Hf**2/(8*sp.pi*G)))
RA = 1/Hf; RAdot = -sp.Derivative(Hf, t)/Hf**2
S = sp.pi*RA**2/G; Sdot = sp.diff(S, t)
kap = -(1/RA)*(1 - RAdot/(2*Hf*RA))          # Hayward-Kodama surface gravity (signed; negative for the inner apparent horizon)
T0 = 1/(2*sp.pi*RA)                           # Cai-Kim 'static' temperature
# (a) Cai-Kim 2005: energy flux through the apparent horizon in dt equals T0 dS
flux = 4*sp.pi*RA**2*(rho + p)*Hf*RA
ck = on_shell(flux - T0*Sdot)
print("  Cai-Kim:  A(rho+p)H R_A - T0 dS/dt =", ck, "->", "PASS (identity with T0 = 1/2 pi R_A)" if ck == 0 else "FAIL")
ok_all &= (ck == 0)
# (b) Akbar-Cai 2007 unified first law: dE = (kappa/2pi) dS + W dV, E = Misner-Sharp = rho V, W = (rho-p)/2, kappa SIGNED
V = sp.Rational(4, 3)*sp.pi*RA**3; E = rho*V; W = (rho - p)/2
ac = on_shell(sp.diff(E, t) - ((kap/(2*sp.pi))*Sdot + W*sp.diff(V, t)))
print("  Akbar-Cai: dE/dt - [(kappa/2pi) dS/dt + W dV/dt] =", ac, "->", "PASS (identity with the signed Hayward kappa, |kappa|/2pi = T_KH)" if ac == 0 else "FAIL")
ok_all &= (ac == 0)
print("  => both horizon 'first laws' are identities of the Einstein equations GRUT III already has from D3; T_KH (with its sign)")
print("     is the bookkeeping temperature that makes Hayward's unified first law close, not a KMS parameter supplied by P1.")
print("     (Tian-Booth 2014 compare the Cai-Kim T0 = 1/(2 pi R_A) with Hayward's kappa/2pi and adopt T0; their positive-heat-out")
print("      convention fixes the sign of the Gibbs relation, not the temperature. With +|kappa|/2pi the Akbar-Cai law does NOT close;")
print("      it closes with the signed kappa.)")

print("== 12.4 lemma: massless field in power-law FRW == massive field in de Sitter (same nu) ==")
eta, nu, k = sp.symbols('eta nu k', positive=True)
a_pl = (-eta)**(sp.Rational(1, 2) - nu)      # eta < 0; a ∝ (-eta)^(1/2-nu)
etan = sp.symbols('etan', positive=True)     # etan = -eta > 0
a_pl = etan**(sp.Rational(1, 2) - nu)
app_over_a = sp.simplify(sp.diff(a_pl, etan, 2)/a_pl)
target = (nu**2 - sp.Rational(1, 4))/etan**2
print("  a''/a =", app_over_a, " vs (nu^2-1/4)/eta^2:", "PASS" if sp.simplify(app_over_a - target) == 0 else "FAIL")
ok_all &= (sp.simplify(app_over_a - target) == 0)
# relation of nu to the slow-roll parameter, DERIVED: a ∝ t^p, eps = 1/p; eta = int dt/a = t^(1-p)/(1-p) (< 0 for p > 1)
pp = sp.symbols('p', positive=True); tt = sp.symbols('t_', positive=True)
eta_of_t = tt**(1 - pp)/(1 - pp)
expo = sp.simplify((sp.diff(sp.log(tt**pp), tt))/(sp.diff(sp.log(-eta_of_t), tt)))     # d ln a / d ln(-eta)
print("  a ∝ t^p: d ln a/d ln(-eta) =", expo, "= -p/(p-1):", "PASS" if sp.simplify(expo + pp/(pp - 1)) == 0 else "FAIL")
ok_all &= (sp.simplify(expo + pp/(pp - 1)) == 0)
nu_pl = sp.Rational(1, 2) - expo                     # from a ∝ (-eta)^(1/2 - nu)
nu_eps = sp.simplify(nu_pl.subs(pp, 1/eps))
print("  nu(eps) =", nu_eps, " = 3/2 + eps/(1-eps):", "PASS" if sp.simplify(nu_eps - (sp.Rational(3, 2) + eps/(1-eps))) == 0 else "FAIL")
ok_all &= (sp.simplify(nu_eps - (sp.Rational(3, 2) + eps/(1-eps))) == 0)
print("  dS analogue mass: m^2/H^2 = 9/4 - nu^2 =", sp.simplify(sp.Rational(9, 4) - nu_eps**2), " (< 0: the dS analogue is tachyonic; the position correlator is IR divergent for every eps > 0,")
print("   the velocity correlator is IR finite for nu < 7/2, i.e. eps < 2/3)")

print("\nCONCLUSION (GRUT III's answer to the two-temperature question):")
print("  (b') They differ because they are different objects. P1's bath temperature is the KMS parameter of the restricted")
print("  state where the modular flow is geometric; for the comoving observer in the conformal vacuum that is T_det = (a''/a')/2pi")
print("  at first order (derived, v11). The apparent horizon is not a Killing horizon for eps != 0, so P1 assigns it no")
print("  temperature; T_KH is Hayward's dynamical surface gravity, used with its sign by Akbar-Cai as the bookkeeping that turns the")
print("  unified first law into an identity of the Einstein equations (12.3); Cai-Kim and Cai-Cao-Hu use T0 = 1/(2 pi R_A) instead. No")
print("  heat-flow law is implied; the O(eps) mismatch is HALF the apparent horizon's expansion rate relative to T_KH (T_det/T_KH - 1 =")
print("  -R_A_dot/2) and the full rate relative to T0 (T_det/T0 - 1 = -eps exactly). Conventional statement (v11, fixed-mean-time")
print("  convention): the comoving worldline's two-point spectrum for a conformal field is KMS at first adiabatic order at (a''/a')/2pi,")
print("  not at T_KH or T0 -- conventional QFT, checkable by anyone.")
print("\nV12 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
