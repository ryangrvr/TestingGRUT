#!/usr/bin/env python3
"""
V12 -- Two first-order temperatures in slowly evolving FRW, and which one GRUT III's principle calls a temperature.

Question (reviewer, 2026-10-11): the comoving conformal-vacuum detector reads T_det = (a''/a')/2pi = (H/2pi)(1 - eps) (v11),
while the cosmological apparent horizon carries the Hayward-Kodama surface gravity kappa = -(1/R_A)(1 - R_A_dot/(2 H R_A))
(Hayward 1998; Cai-Kim 2005; Cai-Cao-Hu 2009), i.e. T_KH = (H/2pi)(1 - eps/2) in flat FRW.  They agree in de Sitter and
differ at O(eps).  GRUT III's principle (P1: the bath is the causal complement, its temperature is the KMS parameter of the
restricted state's modular flow where that flow is geometric; P2: local equilibrium at LOCAL Rindler horizons) must say
which is the bath temperature.

Checks:
 12.1  (sympy) flat FRW: R_A = 1/H, R_A_dot = eps  =>  T_KH = (H/2pi)(1 - eps/2) exactly; T_det = (a''/a')/2pi = (H/2pi)(1 - eps);
       T_det/T_KH = 1 - eps/2 + O(eps^2); equal iff eps = 0 (de Sitter).
 12.2  The apparent horizon is a Killing horizon iff R_A_dot = 0 iff eps = 0: the first-order mismatch is exactly the
       non-Killing character of the apparent horizon.  Where the horizon is not Killing, P1 assigns NO temperature to the
       region it bounds (the modular flow is not geometric), so T_KH is a geometric label, not a KMS parameter.
 12.3  P2 is used in D3 only at local Rindler horizons (Jacobson), whose boost kappa cancels; the Friedmann equations follow
       from D3 without any cosmological temperature.  Check: the Hayward unified first law at the apparent horizon with
       T_KH and S = A/4G reproduces the SAME Friedmann equation (Cai-Kim), i.e. it is an identity of the Einstein equations,
       not an independent input -- verified symbolically for flat FRW with a perfect fluid.
 12.4  Lemma for the next computation (O5 item 1, other fields): a massless minimally coupled field in power-law FRW,
       a ∝ (-eta)^(1/2 - nu), obeys v'' + (k^2 - (nu^2 - 1/4)/eta^2) v = 0 -- the SAME mode equation as a field of mass m
       in de Sitter with nu^2 = 9/4 - m^2/H^2.  Hence its worldline correlator is the dS massive one times a ratio of
       scale factors; the first-order temperature of its stationary (velocity) part is computable from known dS functions.
       (sympy: a''/a = (nu^2 - 1/4)/eta^2 for a ∝ (-eta)^(1/2 - nu).)
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

print("== 12.2 Killing character of the apparent horizon ==")
print("  R_A_dot = eps: the apparent horizon is stationary (Killing) iff eps = 0. For eps != 0 the modular flow of the")
print("  restricted state is not geometric, P1 assigns no temperature to the region inside R_A, and T_KH is a surface-gravity label.")

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
print("     (The sign/'heat-out' convention is the 'temperature confusion' of Tian-Booth 2014.)")

print("== 12.4 lemma: massless field in power-law FRW == massive field in de Sitter (same nu) ==")
eta, nu, k = sp.symbols('eta nu k', positive=True)
a_pl = (-eta)**(sp.Rational(1, 2) - nu)      # eta < 0; a ∝ (-eta)^(1/2-nu)
etan = sp.symbols('etan', positive=True)     # etan = -eta > 0
a_pl = etan**(sp.Rational(1, 2) - nu)
app_over_a = sp.simplify(sp.diff(a_pl, etan, 2)/a_pl)
target = (nu**2 - sp.Rational(1, 4))/etan**2
print("  a''/a =", app_over_a, " vs (nu^2-1/4)/eta^2:", "PASS" if sp.simplify(app_over_a - target) == 0 else "FAIL")
ok_all &= (sp.simplify(app_over_a - target) == 0)
# relation of nu to the slow-roll parameter: a ∝ t^p, eps = 1/p; conformal time: a ∝ (-eta)^(-p/(p-1)) => 1/2 - nu = -p/(p-1)
pp = sp.symbols('p', positive=True)
nu_pl = sp.Rational(1, 2) + pp/(pp - 1)
nu_eps = sp.simplify(nu_pl.subs(pp, 1/eps))
print("  nu(eps) =", nu_eps, " = 3/2 + eps/(1-eps):", "PASS" if sp.simplify(nu_eps - (sp.Rational(3, 2) + eps/(1-eps))) == 0 else "FAIL")
ok_all &= (sp.simplify(nu_eps - (sp.Rational(3, 2) + eps/(1-eps))) == 0)
print("  dS analogue mass: m^2/H^2 = 9/4 - nu^2 =", sp.simplify(sp.Rational(9, 4) - nu_eps**2), " (< 0: the dS analogue is tachyonic; the IR needs care, the velocity correlator does not)")

print("\nCONCLUSION (GRUT III's answer to the two-temperature question):")
print("  (b') They differ because they are different objects. P1's bath temperature is the KMS parameter of the restricted")
print("  state where the modular flow is geometric; for the comoving observer in the conformal vacuum that is T_det = (a''/a')/2pi")
print("  at first order (derived, v11). The apparent horizon is not a Killing horizon for eps != 0, so P1 assigns it no")
print("  temperature; T_KH is Hayward's dynamical surface gravity, used by Cai-Kim as the bookkeeping that turns the horizon")
print("  first law into an identity of the Einstein equations (12.3). No heat-flow law is implied; the O(eps) mismatch equals")
print("  the apparent horizon's expansion rate R_A_dot = eps. Testable statement: a comoving detector coupled to a conformal field")
print("  in the conformal vacuum thermalises (to first order) at (a''/a')/2pi, not at T_KH -- conventional QFT, checkable by anyone.")
print("\nV12 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
