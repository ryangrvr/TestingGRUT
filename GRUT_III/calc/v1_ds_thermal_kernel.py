#!/usr/bin/env python3
"""
V1 -- The de Sitter medium seen by a geodesic (comoving) observer: the noise kernel is EXACTLY thermal,
the local response is EXACTLY local (Huygens), and the asymptotic decay rate of the noise correlator is derived, not supplied.

Checks (all exact statements, verified numerically/symbolically here):
 1.1  a(eta)a(eta')(eta-eta')^2 = 4 sinh^2(H dt/2)/H^2  -- the identity that makes every conformal
      correlator stationary in cosmic proper time along a comoving worldline.
 1.2  Conformal-scalar Bunch-Davies Wightman function on the worldline,
         g(tau) = -(H^2/16 pi^2) / sinh^2(H (tau - i eps)/2),
      has Fourier transform  G+(w) = (w/2pi) / (1 - exp(-beta w)),  beta = 2 pi/H   (Planckian, exactly).
      Symmetrised spectrum S(w) = (w/2pi) coth(beta w/2);  S(0) = H/(2 pi^2)  (classical Nyquist floor).
      [This is the object the record's calc/RESULTS_worldline_reduction.md left 'UNRESOLVED'.]
 1.3  Image-sum identity  sum_n 1/(z + i n beta)^2 = (pi/beta)^2 / sinh^2(pi z/beta): the dS worldline
      correlator IS the flat thermal correlator at T = H/2pi.
 1.4  Flat-space <T00(x) T00(0)>_c for one real conformal scalar at Euclidean separation t along the
      time axis = (3/4) C_T / t^8 with C_T = 1/(3 pi^4), i.e. 1/(4 pi^4 t^8); checked by explicit Wick
      contraction (sympy), together with tracelessness <T_mumu T00> = 0.
      Hence on the dS worldline:  <T^00(t) T^00(0)>_c = (1/4pi^4) [H / (2 sinh(H dt/2))]^8.
 1.5  Closed form of its Wightman spectrum (derived here by contour shift + sech^2n recursion):
         G+_T(w) = w (w^2+H^2)(w^2+4H^2)(w^2+9H^2) / (10080 pi^3 (1 - exp(-2 pi w/H)))
      -> flat limit w^7/(10080 pi^3); KMS detailed balance G+(w)/G+(-w) = e^{beta w};
      dissipative part Im chi ~ odd POLYNOMIAL of degree 7 (=> local-in-time response, Huygens);
      classical floor G+_T(0) = H^7/(560 pi^4);  asymptotic time-domain decay rate 4H, i.e. e^{-4H|tau|} for
      H|tau| >> 1 (poles at +-4iH, +-5iH, ... after the polynomial cancels the Matsubara rungs n=1,2,3).
      4H is a large-|tau| decay RATE, not a 'memory time': the pointwise correlator diverges as tau^-8 at short
      times, so no finite normalised autocorrelation time exists without smearing (see V3, 3.2).
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
H = mp.mpf(1)
beta = 2*mp.pi/H
ok_all = True
def report(name, got, want, tol=1e-9):
    global ok_all
    err = abs(mp.mpf(got) - mp.mpf(want))
    rel = err/abs(mp.mpf(want)) if want != 0 else err
    ok = rel < tol
    ok_all &= bool(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {mp.nstr(got,12)}  want {mp.nstr(want,12)}  rel.err {mp.nstr(rel,3)}")
    return ok

print("== 1.1 stationarity identity (sympy) ==")
Hs, t1, t2 = sp.symbols('H t1 t2', positive=True)
eta1, eta2 = -sp.exp(-Hs*t1)/Hs, -sp.exp(-Hs*t2)/Hs
a1, a2 = -1/(Hs*eta1), -1/(Hs*eta2)
lhs = a1*a2*(eta1-eta2)**2
rhs = 4*sp.sinh(Hs*(t1-t2)/2)**2/Hs**2
d = sp.simplify((lhs-rhs).rewrite(sp.exp))
print("  a a' (eta-eta')^2 - 4 sinh^2(H dt/2)/H^2 =", d, " ->", "PASS" if d == 0 else "FAIL")
ok_all &= (d == 0)

print("== 1.2 worldline Wightman spectrum: contour-shifted quadrature vs Planck formula ==")
# Shift tau -> u - i beta/2: sinh(H(u - i pi/H)/2) = -i cosh(Hu/2), so g -> +(H^2/16pi^2) sech^2(Hu/2) (smooth).
def Gplus_scalar(w):
    f = lambda u: mp.e**(1j*w*u) * (H**2/(16*mp.pi**2)) / mp.cosh(H*u/2)**2
    val = mp.quad(f, [-mp.inf, -20, -8, -3, -1, 0, 1, 3, 8, 20, mp.inf])
    return mp.re(mp.e**(beta*w/2) * val)
def planck_scalar(w):
    w = mp.mpf(w)
    return (w/(2*mp.pi)) / (1 - mp.e**(-beta*w))
for w in [0.1, 0.5, 1, 2, 5, -0.5, -2]:
    report(f"G+(w={w})", Gplus_scalar(mp.mpf(w)), planck_scalar(w))
S01 = planck_scalar(0.1) + planck_scalar(-0.1)
Gp01 = planck_scalar(0.1)
record_val = mp.mpf('0.03393')
gap = (Gp01 - record_val)/Gp01
print(f"  symmetrised S(0.1) [flat line w/2pi] = {mp.nstr(S01,8)}  (= (w/2pi) coth(pi w/H));  exact floor S(0)=H/(2pi^2)={mp.nstr(H/(2*mp.pi**2),8)}")
print(f"  one-sided G+(0.1) = {mp.nstr(Gp01,8)};  the record's calc/worldline_reduction.py 'S(w)' is 2 Re Int_0^inf e^(iw tau) g = G+(w),")
print(f"  so its quadrature value 0.03393 is to be compared with G+(0.1): relative gap (G+ - record)/G+ = {mp.nstr(100*gap,3)}%.")
print(f"  The identification 'record value = exact Planck G+' holds only if the record's eps-extrapolated quadrature precision")
print(f"  is coarser than {mp.nstr(100*gap,3)}%; the record itself quotes <3% validation against the flat half-line result, so it is.")

print("== 1.3 image-sum identity ==")
z, b = mp.mpf('0.7'), beta
lhs13 = mp.nsum(lambda n: 1/(z + 1j*n*b)**2, [-mp.inf, mp.inf])
rhs13 = (mp.pi/b)**2 / mp.sinh(mp.pi*z/b)**2
report("sum_n 1/(z+i n beta)^2 vs (pi/beta)^2/sinh^2(pi z/beta)", mp.re(lhs13), rhs13)

print("== 1.4 Wick contraction: <T00 T00>_c and tracelessness for the improved conformal scalar ==")
t, x1, x2, x3 = sp.symbols('t x1 x2 x3', real=True)
X = [t, x1, x2, x3]
G = 1/(4*sp.pi**2*(t**2+x1**2+x2**2+x3**2))
xi = sp.Rational(1, 6)
def T_terms(mu, nu):
    """T_{mu nu} = d_mu phi d_nu phi - 1/2 delta (dphi)^2 + xi (delta d^2 - d_mu d_nu) phi^2
       as a list of (coef, multi-index on first phi, multi-index on second phi)."""
    terms = [(sp.Integer(1), (mu,), (nu,))]
    if mu == nu:
        for al in range(4):
            terms.append((-sp.Rational(1, 2), (al,), (al,)))
            terms.append((2*xi, (al,), (al,)))
            terms.append((2*xi, (), (al, al)))
    terms.append((-2*xi, (mu,), (nu,)))
    terms.append((-2*xi, (), (mu, nu)))
    return terms
def Dx(expr, idx):
    for i in idx:
        expr = sp.diff(expr, X[i])
    return expr
cache = {}
def corr(a, bidx):
    key = tuple(sorted(a + bidx)), len(bidx) % 2
    if key not in cache:
        cache[key] = (-1)**len(bidx) * Dx(G, a + bidx)
    return cache[key]
def two_point(A, B):
    tot = 0
    for ca, a1, a2 in A:
        for cb, b1, b2 in B:
            tot += ca*cb*(corr(a1, b1)*corr(a2, b2) + corr(a1, b2)*corr(a2, b1))
    return tot
T00 = T_terms(0, 0)
res = sp.simplify(two_point(T00, T00).subs({x1: 0, x2: 0, x3: 0}))
want = 1/(4*sp.pi**4*t**8)
print("  <T00(t) T00(0)>_c =", res, "  expected (3/4)C_T/t^8 = 1/(4 pi^4 t^8):", "PASS" if sp.simplify(res-want) == 0 else "FAIL")
ok_all &= (sp.simplify(res-want) == 0)
trace = 0
for mu in range(4):
    trace += two_point(T_terms(mu, mu), T00)
trace = sp.simplify(trace.subs({x1: 0, x2: 0, x3: 0}))
print("  sum_mu <T_mumu(t) T00(0)>_c =", trace, " ->", "PASS" if trace == 0 else "FAIL")
ok_all &= (trace == 0)
# generic point (not on the axis) tracelessness, as a stronger check
tr2 = sp.simplify(sum(two_point(T_terms(m, m), T00) for m in range(4)).subs({x1: sp.Rational(1,3), x2: sp.Rational(-2,7), x3: sp.Rational(5,11), t: sp.Rational(2,5)}))
print("  tracelessness at a generic point:", tr2, " ->", "PASS" if tr2 == 0 else "FAIL")
ok_all &= (tr2 == 0)

print("== 1.5 stress-tensor Wightman spectrum on the dS worldline: quadrature vs closed form ==")
def Gplus_T(w):
    f = lambda u: mp.e**(1j*w*u) * (1/(4*mp.pi**4)) * (H/2)**8 / mp.cosh(H*u/2)**8
    val = mp.quad(f, [-mp.inf, -20, -8, -3, -1, 0, 1, 3, 8, 20, mp.inf])
    return mp.re(mp.e**(beta*w/2) * val)
def closed_T(w):
    w = mp.mpf(w)
    return w*(w**2+H**2)*(w**2+4*H**2)*(w**2+9*H**2) / (10080*mp.pi**3*(1 - mp.e**(-2*mp.pi*w/H)))
for w in [0.05, 0.5, 1, 3, 8, -1, -3]:
    report(f"G+_T(w={w})", Gplus_T(mp.mpf(w)), closed_T(w))
for w in [1, 3]:
    report(f"KMS ratio G+_T({w})/G+_T({-w}) = e^(beta w)", closed_T(w)/closed_T(-w), mp.e**(beta*w))
for w in [50, 200]:
    report(f"flat limit G+_T({w}) / (w^7/10080 pi^3) -> 1", closed_T(w)/(mp.mpf(w)**7/(10080*mp.pi**3)), 1, tol=1e-2 if w == 50 else 1e-3)
floor = mp.limit(lambda w: closed_T(w), 0)
report("classical floor G+_T(0) = H^7/(560 pi^4)", floor, H**7/(560*mp.pi**4))
# asymptotic time-domain decay rate of the noise kernel: K8(tau) ~ e^{-4 H tau} for H tau >> 1
K8 = lambda tau: (1/(4*mp.pi**4))*(H/(2*mp.sinh(H*tau/2)))**8
rate = -mp.log(K8(12)/K8(10))/2
report("noise-kernel asymptotic decay rate (tau=10..12) -> 4H", rate, 4*H, tol=1e-3)
print("  dissipative part Im chi_T(w) propto G+(w)-G+(-w) = w(w^2+H^2)(w^2+4H^2)(w^2+9H^2)/(10080 pi^3): an odd polynomial")
print("  -> the retarded response at a point is a finite sum of delta-function derivatives (no memory: Huygens);")
print("  -> all memory sits in the NOISE, with derived ASYMPTOTIC DECAY RATE 4H (large-|tau| tail e^{-4H|tau|}) and Planckian occupation at T=H/2pi.")
print("  NOTE: 4H is a decay rate, not a 'memory time'. The pointwise <T00 T00> correlator diverges as tau^-8 as tau -> 0,")
print("  so no finite normalised autocorrelation time Int K dtau / K(0) exists without smearing the operator over a finite region.")
print("\nV1 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
