#!/usr/bin/env python3
"""
V13 -- Does any horizon temperature fix the O(eps) Fokker-Planck coefficients of the massless zero mode in power-law FRW?
(The decisive cheap computation proposed by the O5-3 survey, `10`.)

Conventional input (Cho, Kim & Kitamoto 2015, arXiv:1508.07877, eqs. 2.5, 2.13, 3.14-3.16; read in full text):
  a = (1 + eps H0 t)^(1/eps),  nu0 = (3 - eps)/(2(1 - eps)) = 3/2 + eps/(1 - eps)   [= the v12 12.4 lemma index],
  Langevin noise  <phibar0_dot(t) phibar0_dot(t')> = A0 H^3 delta(t - t'),   A0 = (2 - 2 eps)^(2 nu0 - 3) Gamma^2(nu0) (1 - eps)^3 / pi^3,
  Fokker-Planck   rho_dot = (1/2) A0 H^3 d^2 rho/dphi^2 + (1/((3 - eps) H)) d/dphi (rho V').
At eps = 0: A0 = 1/(4 pi^2), i.e. Starobinsky's H^3/(8 pi^2) diffusion and 1/(3H) mobility.

The stationary solution is rho ∝ exp(-V/T_zm) with the ZERO-MODE TEMPERATURE  T_zm = (3 - eps) A0 H^4 / 2.
In exact de Sitter T_zm = 3 H^4/(8 pi^2) = 1/Vol(S^4): the Hartle-Hawking weight exp(-S_E), the Euclidean face of the KMS
state at beta = 2 pi/H.  Question: in power-law FRW, is T_zm(eps) the Euclidean weight of SOME horizon temperature, i.e.
T_zm = 3 H_X^4/(8 pi^2) with H_X = H(1 - c_X eps), for the detector value c = 1 (a''/a', D11), the Kodama-Hayward value
c = 1/2 (09), or the naive c = 0?

Results (exact, sympy):
 13.1  d ln A0/d eps |_0 = 1 - 2 gamma_E - 2 ln 2 = -1.5407...  (checked against a direct numerical evaluation of A0(eps)).
 13.2  T_zm(eps)/T_zm(0) = 1 + (c1 - 1/3) eps + O(eps^2) = 1 - 1.874 eps + ...
 13.3  The Euclidean-volume reading requires (H_X/H)^4 = T_zm(eps)/T_zm(0), i.e. H_X = H(1 - 0.4685 eps): c_X = (2 gamma_E + 2 ln 2 - 2/3)/4,
       which matches NEITHER the detector temperature (c = 1) NOR the Kodama-Hayward one (c = 1/2) NOR H itself (c = 0).
 13.4  The diffusion coefficient alone: D(eps)/D(0) = 1 + c1 eps with c1 = -1.54, vs (H_det/H)^3 = 1 - 3 eps and (H_KH/H)^3 = 1 - 1.5 eps:
       again no match (the KH value is numerically close, -1.50 vs -1.54, but differs by (2 gamma_E + 2 ln 2 - 5/2)).
Conclusion: the O(eps) coefficients of the massless zero mode are fixed by the Hankel-function normalisation (Gamma(nu0), 2^(2nu0-3)),
carrying Euler's gamma and ln 2; they are not the Euclidean weight of any geometric horizon temperature.  P1 is consistent with
this -- for eps != 0 the modular flow is not geometric and P1 assigns no temperature -- and it adds nothing: conventional
theory fixes the coefficients completely.  O5 items 1 and 3 close for GRUT III.
"""
import sympy as sp
import mpmath as mp
ok_all = True
eps = sp.symbols('epsilon', positive=True)
nu0 = (3 - eps)/(2*(1 - eps))
A0 = (2 - 2*eps)**(2*nu0 - 3)*sp.gamma(nu0)**2*(1 - eps)**3/sp.pi**3
print("== 13.0 checks of the input ==")
ok = sp.simplify(nu0 - (sp.Rational(3, 2) + eps/(1 - eps))) == 0
print("  nu0 = 3/2 + eps/(1-eps):", "PASS" if ok else "FAIL"); ok_all &= ok
A00 = sp.simplify(A0.subs(eps, 0))
print("  A0(0) =", A00, " = 1/(4 pi^2):", "PASS" if sp.simplify(A00 - 1/(4*sp.pi**2)) == 0 else "FAIL"); ok_all &= sp.simplify(A00 - 1/(4*sp.pi**2)) == 0
print("== 13.1 first-order coefficient of the noise amplitude ==")
c1 = sp.simplify(sp.diff(sp.log(A0), eps).subs(eps, 0))
c1_target = 1 - 2*sp.EulerGamma - 2*sp.log(2)
print("  c1 = d ln A0/d eps|0 =", c1, "=", sp.N(c1, 8), " vs 1 - 2 gamma_E - 2 ln 2:", "PASS" if sp.simplify(c1 - c1_target) == 0 else "FAIL")
ok_all &= sp.simplify(c1 - c1_target) == 0
# numerical cross-check of the derivative from A0 itself
mp.mp.dps = 30
def A0num(e):
    e = mp.mpf(e); nu = (3 - e)/(2*(1 - e))
    return (2 - 2*e)**(2*nu - 3)*mp.gamma(nu)**2*(1 - e)**3/mp.pi**3
h = mp.mpf('1e-6'); c1_num = (mp.log(A0num(h)) - mp.log(A0num(-h)))/(2*h)
ok = abs(c1_num - mp.mpf(sp.N(c1, 30))) < 1e-9
print(f"  numerical central difference: {mp.nstr(c1_num, 10)} ->", "PASS" if ok else "FAIL"); ok_all &= ok
print("== 13.2 zero-mode temperature T_zm = (3 - eps) A0 H^4 / 2 ==")
H = sp.symbols('H', positive=True)
Tzm = (3 - eps)*A0*H**4/2
T0 = sp.simplify(Tzm.subs(eps, 0))
print("  T_zm(0) =", T0, " = 3H^4/(8 pi^2) = 1/Vol(S^4):", "PASS" if sp.simplify(T0 - 3*H**4/(8*sp.pi**2)) == 0 else "FAIL"); ok_all &= sp.simplify(T0 - 3*H**4/(8*sp.pi**2)) == 0
k = sp.simplify(sp.diff(sp.log(Tzm), eps).subs(eps, 0))
print("  d ln T_zm/d eps|0 =", k, "=", sp.N(k, 8), " (= c1 - 1/3):", "PASS" if sp.simplify(k - (c1 - sp.Rational(1, 3))) == 0 else "FAIL")
ok_all &= sp.simplify(k - (c1 - sp.Rational(1, 3))) == 0
print("== 13.3 which horizon temperature would reproduce it as a Euclidean weight? ==")
cX = sp.simplify(-k/4)
print("  H_X = H(1 - c_X eps) with c_X = -(c1 - 1/3)/4 =", cX, "=", sp.N(cX, 6))
for name, c in [("detector (a''/a')", 1), ("Kodama-Hayward", sp.Rational(1, 2)), ("naive H", 0)]:
    match = sp.simplify(cX - c) == 0
    print(f"    vs {name}: c = {c}  ->", "MATCH" if match else f"no match (difference {sp.N(cX - c, 5)})")
    ok_all &= (not match)
print("== 13.4 diffusion coefficient alone ==")
for name, c in [("detector", 1), ("Kodama-Hayward", sp.Rational(1, 2))]:
    print(f"    (H_{name}/H)^3 first-order coefficient = {sp.N(-3*c, 5)}  vs c1 = {sp.N(c1, 5)}  -> no match (difference {sp.N(c1 + 3*c, 5)})")
    ok_all &= sp.simplify(c1 + 3*c) != 0
res = sp.simplify(c1 + sp.Rational(3, 2)); res_target = sp.Rational(5, 2) - 2*sp.EulerGamma - 2*sp.log(2)
ok = (sp.simplify(res - res_target) == 0) and abs(float(res)) > 1e-3
print(f"    Kodama-Hayward near-coincidence is not an identity: c1 + 3/2 = {res} = {sp.N(res, 5)} (transcendental, nonzero):", "PASS" if ok else "FAIL")
ok_all &= ok
print("\nCONCLUSION: the O(eps) Fokker-Planck coefficients of the massless zero mode carry gamma_E and ln 2 from the Hankel normalisation")
print("and are not the Euclidean weight of any geometric horizon temperature; P1 assigns none (non-geometric modular flow) and adds")
print("nothing: conventional theory (Cho-Kim-Kitamoto 2015) fixes them. O5 items 1 and 3 close for GRUT III.")
print("\nV13 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
