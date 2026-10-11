#!/usr/bin/env python3
"""
V13 (rev. 2) -- The O(eps) Fokker-Planck coefficients of the massless zero mode in power-law FRW: can they test a horizon
temperature?  ANSWER (rev. 2): NO -- not because no temperature matches, but because the first-order coefficient is not a
physical invariant: it depends logarithmically on the coarse-graining scale sigma of the stochastic formalism.

History.  Rev. 1 (revision 2.3, commit b6f0e19) took the published coefficients (Cho, Kim & Kitamoto 2015, arXiv:1508.07877,
eqs. 2.5, 2.13, 3.14-3.16), which are evaluated at the coarse-graining scale k = aH, found d ln A0/d eps = 1 - 2 gamma_E - 2 ln 2
and concluded that the implied zero-mode temperature matches no horizon temperature (detector c = 1, Kodama-Hayward c = 1/2,
naive c = 0; required c_X = 0.4685).  The arithmetic was right; the INFERENCE WAS NOT ESTABLISHED (reviewer, fifth round):
the stochastic noise is evaluated at k = sigma a H and sigma is a choice.  In exact dS the noise H^3/(8 pi^2) is
sigma-independent (why Starobinsky's formalism works); at O(eps) it is not.  This revision derives the sigma dependence from
the mode function itself and records the consequence.

Conventional input (CKK 2015):  a = (1 + eps H0 t)^(1/eps) = (-(1-eps) H0 eta)^(-1/(1-eps)),  nu0 = (3-eps)/(2(1-eps)),
  phi_k(eta) = (sqrt(pi)/2) (-eta)^(1/2) H^(1)_nu0(-k eta) / a(eta)   (Bunch-Davies),
  small-argument (frozen) amplitude  |phi_k|^2 = Gamma(nu0)^2 2^(2 nu0) ((1-eps)H0)^(2/(1-eps)) k^(-2 nu0) / (4 pi)   [= CKK (3.14) inner line],
  Fokker-Planck  rho_dot = (1/2) A0 H^3 rho'' + (1/((3-eps)H)) (rho V')',   A0 = (2-2eps)^(2nu0-3) Gamma^2(nu0) (1-eps)^3 / pi^3  [CKK (3.16), sigma = 1].

Derivation of the sigma dependence (13.1):  the coarse-grained field phibar collects modes k < sigma a H; its variance is
  <phibar^2> = (1/2pi^2) int^{sigma a H} k^2 dk |phi_k|^2 = C (sigma a H)^(3-2nu0) / (2 pi^2 (3-2nu0)) + const,
so the white-noise amplitude is d<phibar^2>/dt = C (sigma a H)^(3-2nu0) H (1-eps) / (2 pi^2) = A0 sigma^(3-2nu0) H^3,
i.e.  A0(sigma) = (2(1-eps)/sigma)^(2nu0-3) Gamma^2(nu0) (1-eps)^3 / pi^3,  with 3 - 2nu0 = -2 eps/(1-eps).
Consequences:  c1(sigma) = d ln A0/d eps|0 = 1 - 2 gamma_E - 2 ln(2 sigma);  the "zero-mode temperature" T_zm = (3-eps) A0 H^4/2
(stationary weight exp(-V/T_zm)) has d ln T_zm/d eps|0 = c1(sigma) - 1/3, and reading it as the Euclidean weight (H_X/H)^4 of a
rate H_X = H(1 - c_X eps) gives  c_X(sigma) = (2 gamma_E + 2 ln 2 - 2/3)/4 + (1/2) ln sigma = 0.4685 + 0.5 ln sigma.
So c_X = 0 at sigma = 0.39, = 1/2 (Kodama-Hayward) at sigma = 1.07, = 1 (detector) at sigma = 2.9, and c_X -> -infinity as
sigma -> 0, which is the regime the formalism needs.  The gamma_E and ln 2 of rev. 1 are the sigma = 1 value of gamma_E + ln(2 sigma):
a matching-scale artifact.  The O(eps) zero-mode coefficient cannot test horizon temperatures.

What stands:  O5 items 1 and 3 close on other grounds -- P1 assigns no temperature when the modular flow is not geometric,
and the sigma-independent correlators are fixed by conventional QFT.  The sigma-independent decisive object is the comoving
detector's frequency-resolved KMS ratio for the massless sector at O(eps): `v14`.
"""
import sympy as sp
import mpmath as mp
ok_all = True
def check(label, cond):
    global ok_all
    print(f"  {label}:", "PASS" if cond else "FAIL"); ok_all &= bool(cond)
def num_zero(expr, syms, pts=((sp.Rational(1, 10), sp.Rational(7, 10), sp.Rational(13, 10), sp.Rational(21, 10), sp.Rational(1, 2)),
                               (sp.Rational(3, 100), sp.Rational(3, 2), sp.Rational(4, 5), sp.Rational(1, 3), 2),
                               (sp.Rational(1, 4), sp.Rational(1, 5), 3, 5, sp.Rational(1, 7)))):
    """exact-in-principle identities that sympy cannot simplify through nested symbolic powers: evaluate at rational points to 40 digits"""
    return all(abs(sp.N(expr.subs(dict(zip(syms, pt))), 40)) < sp.Float('1e-30') for pt in pts)

eps, sig, H, H0, k, eta = sp.symbols('epsilon sigma H H_0 k eta', positive=True)
nu0 = (3 - eps)/(2*(1 - eps))
A0 = (2 - 2*eps)**(2*nu0 - 3)*sp.gamma(nu0)**2*(1 - eps)**3/sp.pi**3          # CKK (3.16), sigma = 1

print("== 13.0 checks of the input ==")
check("nu0 = 3/2 + eps/(1-eps)", sp.simplify(nu0 - (sp.Rational(3, 2) + eps/(1 - eps))) == 0)
check("A0(0) = 1/(4 pi^2) (Starobinsky)", sp.simplify(A0.subs(eps, 0) - 1/(4*sp.pi**2)) == 0)

print("== 13.1 the coarse-graining scale: noise amplitude at k = sigma a H, derived from the mode function ==")
# a(eta), H(eta) for CKK (2.5); frozen small-argument amplitude of the Bunch-Davies mode
a_eta = ((1 - eps)*H0*eta)**(-1/(1 - eps))               # eta := -(conformal time) > 0; CKK (2.5)
Heta = -sp.diff(a_eta, eta)/a_eta**2                       # H = (da/d tau)/a^2 with d/d tau = -d/d eta
check("a H (-eta) = 1/(1-eps)", sp.simplify(a_eta*Heta*eta - 1/(1 - eps)) == 0)
check("H = H0 a^(-eps)", num_zero(Heta - H0*a_eta**(-eps), (eps, eta, H0, k, sig)))
# |phi_k|^2 at small argument: (pi/4) (-eta) |H_nu|^2 / a^2 with |H_nu(x)|^2 -> (Gamma(nu)/pi)^2 (x/2)^(-2nu)
phi2 = sp.pi/4*eta*(sp.gamma(nu0)/sp.pi)**2*(k*eta/2)**(-2*nu0)/a_eta**2
expo = sp.simplify(sp.diff(sp.log(phi2), eta)*eta)          # power of eta
check("frozen: |phi_k|^2 is eta-independent (exponent 0)", sp.simplify(expo) == 0)
C = sp.simplify(phi2*k**(2*nu0))                            # the k-independent constant
C_ckk = 2**(2*nu0 - 2)*sp.gamma(nu0)**2*((1 - eps)*H0)**(2/(1 - eps))/sp.pi
check("|phi_k|^2 k^(2nu0) = CKK (3.14) inner line", num_zero(C/C_ckk - 1, (eps, eta, H0, k, sig)))
# variance of the coarse-grained field and its time derivative (the white-noise amplitude)
aH = sp.simplify(a_eta*Heta)
var_top = C*(sig*aH)**(3 - 2*nu0)/(2*sp.pi**2*(3 - 2*nu0))   # upper-limit contribution to <phibar^2>; the IR end is a constant
# d/dt = (1/a) d/d(conformal time) = -(1/a) d/d eta
noise = sp.simplify(-sp.diff(var_top, eta)/a_eta)
A0_sig = sp.simplify(noise/Heta**3)
A0_sig_target = (2*(1 - eps)/sig)**(2*nu0 - 3)*sp.gamma(nu0)**2*(1 - eps)**3/sp.pi**3
check("noise = A0(sigma) H^3 with A0(sigma) = (2(1-eps)/sigma)^(2nu0-3) Gamma^2(nu0) (1-eps)^3/pi^3 (eta-independent)", num_zero(A0_sig/A0_sig_target - 1, (eps, eta, H0, k, sig)))
check("sigma = 1 reproduces CKK (3.16)", num_zero(A0_sig.subs(sig, 1)/A0 - 1, (eps, eta, H0, k, sig)))
check("exact dS (eps = 0): sigma-independent, = 1/(4 pi^2)", sp.simplify(A0_sig.subs(eps, 0) - 1/(4*sp.pi**2)) == 0)
check("A0(sigma)/A0(1) = sigma^(-2 eps/(1-eps))", num_zero(A0_sig/A0_sig.subs(sig, 1) - sig**(-2*eps/(1 - eps)), (eps, eta, H0, k, sig)))

print("== 13.2 first-order coefficient c1(sigma) = d ln A0(sigma)/d eps |0 ==")
c1s = sp.simplify(sp.diff(sp.log(A0_sig_target), eps).subs(eps, 0))
c1_target = 1 - 2*sp.EulerGamma - 2*sp.log(2*sig)
check("c1(sigma) = 1 - 2 gamma_E - 2 ln(2 sigma)", sp.simplify(sp.expand_log(c1s - c1_target, force=True)) == 0)
c1 = sp.simplify(c1s.subs(sig, 1))
print("  c1(1) =", c1, "=", sp.N(c1, 8), " (rev. 1 value)")
mp.mp.dps = 30
def A0num(e, s=1):
    e = mp.mpf(e); nu = (3 - e)/(2*(1 - e))
    return (2*(1 - e)/s)**(2*nu - 3)*mp.gamma(nu)**2*(1 - e)**3/mp.pi**3
h = mp.mpf('1e-6')
for s in ['1', '0.1']:
    num = (mp.log(A0num(h, mp.mpf(s))) - mp.log(A0num(-h, mp.mpf(s))))/(2*h)
    exact = mp.mpf(sp.N(c1s.subs(sig, sp.Rational(s)), 30))
    check(f"numerical central difference at sigma = {s}: {mp.nstr(num, 10)}", abs(num - exact) < 1e-9)

print("== 13.3 the 'zero-mode temperature' and the horizon-temperature comparison, as a function of sigma ==")
Tzm = (3 - eps)*A0_sig_target*H**4/2
check("T_zm(0) = 3H^4/(8 pi^2) = 1/Vol(S^4), any sigma", sp.simplify(Tzm.subs(eps, 0) - 3*H**4/(8*sp.pi**2)) == 0)
ks = sp.simplify(sp.diff(sp.log(Tzm), eps).subs(eps, 0))
check("d ln T_zm/d eps|0 = c1(sigma) - 1/3", sp.simplify(sp.expand_log(ks - (c1s - sp.Rational(1, 3)), force=True)) == 0)
cXs = sp.simplify(-ks/4)
cX_target = (2*sp.EulerGamma + 2*sp.log(2) - sp.Rational(2, 3))/4 + sp.log(sig)/2
check("c_X(sigma) = 0.4685 + (1/2) ln sigma", sp.simplify(sp.expand_log(cXs - cX_target, force=True)) == 0)
cX1 = cXs.subs(sig, 1)
print(f"  c_X(1) = {sp.N(cX1, 6)} (rev. 1 value);  c_X(sigma) = c_X(1) + 0.5 ln sigma")
print("  sigma      c_X(sigma)   diffusion coeff. c1(sigma)")
for s in [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(39, 100), 1, sp.Rational(107, 100), 3]:
    print(f"  {float(s):7.3f}   {float(cXs.subs(sig, s)):+9.4f}   {float(c1s.subs(sig, s)):+9.4f}")
print("  sigma at which c_X equals each candidate:")
for name, c in [("naive H (c = 0)", 0), ("Kodama-Hayward (c = 1/2)", sp.Rational(1, 2)), ("detector a''/a' (c = 1)", 1)]:
    s_star = sp.exp(2*(c - cX1))
    print(f"    {name}: sigma = {sp.N(s_star, 4)}")
    check(f"    c_X(sigma) = {c} is attainable", sp.simplify(sp.expand_log(cXs.subs(sig, s_star) - c, force=True)) == 0)
check("no sigma -> 0 limit: c_X -> -infinity", sp.limit(cXs, sig, 0) == -sp.oo)

print("== 13.4 the rev.-1 near-coincidence, restated ==")
res = sp.simplify(c1 + sp.Rational(3, 2))
print(f"  at sigma = 1: c1 + 3/2 = {res} = {sp.N(res, 5)};  it vanishes at sigma = exp(5/4 - gamma_E - ln 2) = {sp.N(sp.exp(sp.Rational(5,4) - sp.EulerGamma - sp.log(2)), 4)}")
check("the 'transcendental gap' of rev. 1 is gamma_E + ln(2 sigma) at sigma = 1: a matching-scale artifact", sp.simplify(sp.expand_log(c1s - (1 - 2*(sp.EulerGamma + sp.log(2*sig))), force=True)) == 0)

print("\nCONCLUSION (rev. 2): the O(eps) Fokker-Planck coefficients of the massless zero mode depend on the coarse-graining scale sigma")
print("(c1 = 1 - 2 gamma_E - 2 ln(2 sigma); sigma-independent only at eps = 0). They are not a physical invariant and cannot test horizon")
print("temperatures; rev. 1's 'no temperature matches' is withdrawn as an inference (the arithmetic at sigma = 1 stands). O5 items 1 and 3")
print("still close: P1 assigns no temperature for non-geometric modular flow, and sigma-independent correlators are fixed by conventional")
print("QFT (CKK 2015). The sigma-independent test is the comoving detector's KMS ratio at omega != 0 (`v14`).")
print("\nV13 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
