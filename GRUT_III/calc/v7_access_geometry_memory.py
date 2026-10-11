#!/usr/bin/env python3
"""
V7 -- 'Memory' of a response is a property of (observable, state), not of the vacuum:
access geometry (the spatial smearing that defines what is measured) alone produces apparent memory,
and the state alone (via KMS) fixes the noise once the commutator is known.  Flat space, free fields.

WHAT IS ESTABLISHED HERE (each item is a PASS/FAIL check below; exit status 0 iff all pass)

 (A) Free massless scalar, G_R(t,r) = theta(t) delta(t-r)/(4 pi r), both endpoints smeared with the
     normalised Gaussian f_L(x) = exp(-|x|^2/2L^2)/(2 pi L^2)^{3/2}:
   A.1 (sympy, exact)  f_L * f_L = f_{sqrt2 L}  and  chi_L(t) = theta(t) t exp(-t^2/4L^2) / (8 pi^{3/2} L^3).
   A.2 (sympy, exact)  int_0^oo chi_L dt = 1/(4 pi^{3/2} L) = (1/(4 pi L^2)) * (L/sqrt(pi)); it equals the
       smeared static Coulomb potential int d^3r f_{sqrt2 L}(r)/(4 pi r) (independent integral).
       Mean memory time <t> = int t chi / int chi = sqrt(pi) L  (apparent memory ~ L/c, from geometry alone).
       The commutator spectral density is rho_L(w) = (w/2pi) e^{-w^2 L^2}: the pointwise flat line w/2pi
       (v1) times the form factor |f~_L(k=w)|^2.
   A.3 Compact support: for a profile supported in |x|<R, supp(f*f) is in |x|<2R and chi vanishes for t>2R
       (support of the convolution; G_R lives on r=t).  Demonstrated with a top hat: numerical 3D convolution
       agrees with the analytic lens volume (pi/12)(4R+d)(2R-d)^2 (3/(4 pi R^3))^2 for d<2R, is 0 for d>=2R,
       and int_0^{2R} chi dt = 3/(10 pi R) (twice the Coulomb self-energy of a uniform ball).

 (B) Stress tensor of one real conformal scalar (xi=1/6), flat space.
   B.0 (sympy, exact Wick contraction)  <T00(x) T00(0)> = C_T (I_00^2 - 1/4)/(x^2)^4 with C_T = 1/(3 pi^4),
       I_00 = 1 - 2 tau^2/x^2 (Euclidean), at EVERY point, not only on the time axis (v1 checked the axis).
       Lorentzian (mostly plus): W(t,r) = C_T [(1+2t^2/x^2)^2 - 1/4]/(x^2)^4, x^2 = r^2 - (t - i eps)^2.
   B.1 (sympy, exact) Closed form of the smeared retarded kernel chi_L(t) = -2 theta(t) Im W_L(t),
       W_L(t) = int d^3r f_{sqrt2 L}(r) W(t,r), obtained two independent ways which agree identically:
        (a) position space, Im 1/(u - t^2 + i0)^n = -pi (-1)^{n-1} delta^{(n-1)}(u-t^2)/(n-1)!, after writing
            the numerator as 3/(4a^4) - 4u/a^5 + 4u^2/a^6 (a = x^2, u = r^2);
        (b) momentum space: the CFT spectral function rho_{T00T00}(w,K) = K^4 theta(w-K)/(720 pi) (spin-2
            projector Pi_{00,00} = (2/3) K^4/(p^2)^2 times (p^2)^2; normalisation fixed by v1's w^7/(10080 pi^3)),
            so rho_L(w) = (1/(1440 pi^3)) int_0^w K^6 e^{-K^2 L^2} dK and chi_L(t) = (1/pi) int_0^oo rho_L sin(wt) dw.
       Result:  chi_L(t) = (120 L^6 - 180 L^4 t^2 + 30 L^2 t^4 - t^6) e^{-t^2/4L^2} / (184320 pi^{7/2} L^13 t).
       Small t: chi_L -> 1/(1536 pi^{7/2} L^7 t) = rho_L(oo)/(pi t): spatial smearing alone leaves a 1/t edge
       (the dissipative spectrum tends to a constant, not zero, at high frequency: a spatially smeared T00 still
       creates pairs of arbitrarily high energy with bounded total momentum).  Large t: Gaussian e^{-t^2/4L^2}.
   B.2 (mpmath) Direct quadrature of W_L(t) on the rotated radial contour r = s e^{i pi/8} (exact eps->0 limit)
       reproduces the closed form to 1e-12 at t/L = 0.05 ... 4; the rotation is validated against real-axis
       quadrature at finite eps.
   B.3 Scaling: L^8 chi_L(t; L) at fixed t/L is L-independent for L = 1/2, 1, 2 (numerically, not just by formula).
       NOTE: the task text said L^{-9}; the correct count for normalised smearing is L^{-8}, because
       F_L = int d^3x f_L T00 has mass dimension 4, so <F_L F_L> has dimension 8.
   B.4 (numpy) Fourier transform W~_L(w) = int dt e^{iwt} W_L(t), computed exactly on the shifted line
       Im t = -L (analyticity of the Wightman function in the lower half plane): W~_L(w) = rho_L^{(b)}(w) for w>0
       to 1e-10, |W~_L(-w)|/W~_L(w) < 1e-9 (positive frequency, T=0), small-w limit -> v1's w^7/(10080 pi^3),
       large-w plateau sqrt(pi)/(1536 pi^3 L^7).

 (C) KMS for smeared observables.  Thermal Wightman function of the massless scalar FIELD by image sum,
     G_beta(t,r) = (1/(8 pi r beta)) [coth(pi(r-t)/beta) + coth(pi(r+t)/beta)] (checked against the raw sum),
     smeared in space with f_{sqrt2 L}, Fourier transformed on the generic shifted line Im t = -beta/4 (a
     symmetric shift beta/2 would make detailed balance tautological; beta/4 does not).  For
     (beta, L) = (2 pi, 1), (2 pi, 1/2), (1, 0.3):  G+_L(w)/G+_L(-w) = e^{beta w},
     rho_L = G+_L(w) - G+_L(-w) = (w/2pi) e^{-w^2 L^2}  (the state-independent commutator of (A)),
     N_L = (G+_L(w)+G+_L(-w))/2 = (1/2) coth(beta w/2) rho_L(w),  all to better than 1e-6.

 CONCLUSION.  Given the commutator (fixed by the observable's access geometry: the smearing), the state fixes the
 noise for ANY smearing (KMS, (C)); the smearing fixes the apparent memory ((A): <t> = sqrt(pi) L, support < 2R for
 compact profiles; (B): Gaussian tail e^{-t^2/4L^2} of the T00 kernel).  Neither is a free 'vacuum memory'.

WHAT IS NOT ESTABLISHED HERE
 - Nothing about gravity, de Sitter or horizons: flat space, free (Gaussian) fields only.
 - The T00 kernel with spatial smearing alone is NOT an integrable function: the 1/t edge means the static
   susceptibility int chi_L dt diverges logarithmically in the UV; a temporal smearing (detector time resolution)
   is also needed.  The task's phrasing 'chi_L is a smooth function with support scale ~L' holds for t>0 but
   not uniformly down to t=0.
 - rho_{T00T00}(w,K) = K^4 theta(w-K)/(720 pi) is derived from conformal/Ward structure plus v1's normalisation
   and then confirmed by the independent position-space closed form and by quadrature; no interacting corrections.
"""
import sys, time
import numpy as np
import sympy as sp
import mpmath as mp
from numpy.polynomial.legendre import leggauss

mp.mp.dps = 30
ok_all = True
T0 = time.time()

def report(name, got, want, tol=1e-9, absolute=False):
    global ok_all
    got, want = mp.mpf(got), mp.mpf(want)
    err = abs(got - want)
    rel = err if (absolute or want == 0) else err/abs(want)
    ok = bool(rel < tol)
    ok_all &= ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {mp.nstr(got,12)}  want {mp.nstr(want,12)}  "
          f"{'abs' if (absolute or want==0) else 'rel'}.err {mp.nstr(rel,3)}")
    return ok

def report_bool(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{(': ' + detail) if detail else ''}")
    return bool(cond)

def gl(a, b, n):
    x, w = leggauss(n)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

# ====================================================================================================
print("== A. free scalar: smeared retarded kernel is fixed by the access geometry ==")
L, t, r, d, R = sp.symbols('L t r d R', positive=True)
xs = sp.symbols('xs', real=True)
f1 = lambda y, s: sp.exp(-y**2/(2*s**2))/sp.sqrt(2*sp.pi*s**2)                 # 1D factor of f_s
f3 = lambda rr, s: sp.exp(-rr**2/(2*s**2))/(2*sp.pi*s**2)**sp.Rational(3, 2)   # 3D normalised Gaussian

print("-- A.1 convolution of Gaussians and radial integral against delta(t-r) (sympy) --")
conv1 = sp.simplify(sp.integrate(f1(xs, L)*f1(d - xs, L), (xs, -sp.oo, sp.oo)))
report_bool("f_L * f_L = f_{sqrt2 L} (per Cartesian factor; 3D Gaussian factorises)",
            sp.simplify(conv1 - f1(d, sp.sqrt(2)*L)) == 0, f"(f_L*f_L)(d) = {conv1}")
chi_sym = sp.simplify(sp.integrate(4*sp.pi*r**2*f3(r, sp.sqrt(2)*L)*sp.DiracDelta(t - r)/(4*sp.pi*r), (r, 0, sp.oo)))
chi_want = t*sp.exp(-t**2/(4*L**2))/(8*sp.pi**sp.Rational(3, 2)*L**3)
report_bool("chi_L(t) = int d^3r f_{sqrt2 L}(r) delta(t-r)/(4 pi r) = t exp(-t^2/4L^2)/(8 pi^{3/2} L^3)",
            sp.simplify(chi_sym - chi_want) == 0, f"chi_L = {chi_sym}")
# numerical cross-check of the Gaussian convolution by direct 3D quadrature (same routine used for the top hat below)
def conv3d_num(fr, dv, breaks):
    """(f*f)(d) = 2 pi int r^2 dr int_{-1}^{1} dmu f(r) f(sqrt(r^2+d^2-2 r d mu)) for a radial profile f."""
    def inner(rr):
        return fr(rr)*mp.quad(lambda mu: fr(mp.sqrt(rr**2 + dv**2 - 2*rr*dv*mu)), [-1, 1])
    return 2*mp.pi*mp.quad(lambda rr: rr**2*inner(rr), breaks)
mp.mp.dps = 15
Lv = mp.mpf(1)
fG = lambda rr: mp.e**(-rr**2/(2*Lv**2))/(2*mp.pi*Lv**2)**mp.mpf(1.5)
dv = mp.mpf('1.3')
report("direct 3D quadrature (f_L*f_L)(1.3 L) vs f_{sqrt2 L}(1.3 L)", conv3d_num(fG, dv, [0, dv, 2*dv, 8]),
       mp.e**(-dv**2/(4*Lv**2))/(4*mp.pi*Lv**2)**mp.mpf(1.5), tol=1e-8)
mp.mp.dps = 30

print("-- A.2 normalisation, static limit, mean memory time (sympy) --")
norm = sp.simplify(sp.integrate(chi_want, (t, 0, sp.oo)))
report_bool("int_0^oo chi_L dt = 1/(4 pi^{3/2} L)", sp.simplify(norm - 1/(4*sp.pi**sp.Rational(3, 2)*L)) == 0,
            f"= {norm} = (1/(4 pi L^2)) * ({sp.simplify(norm*4*sp.pi*L**2)})")
static = sp.simplify(sp.integrate(4*sp.pi*r**2*f3(r, sp.sqrt(2)*L)/(4*sp.pi*r), (r, 0, sp.oo)))
report_bool("= smeared static Coulomb potential int d^3r f_{sqrt2 L}/(4 pi r) (independent integral)",
            sp.simplify(static - norm) == 0, f"{static}")
mean_t = sp.simplify(sp.integrate(t*chi_want, (t, 0, sp.oo))/norm)
report_bool("mean memory time <t> = sqrt(pi) L", sp.simplify(mean_t - sp.sqrt(sp.pi)*L) == 0, f"<t> = {mean_t}")
w = sp.symbols('w', real=True)
# commutator function <[F(t),F(0)]> = -i sgn(t) chi_L(|t|); its FT is the spectral density rho_L(w)
rho_A = sp.simplify(sp.integrate(sp.exp(sp.I*w*t)*(-sp.I)*t*sp.exp(-t**2/(4*L**2))/(8*sp.pi**sp.Rational(3, 2)*L**3),
                                 (t, -sp.oo, sp.oo)))
report_bool("rho_L(w) = FT of <[F_L(t),F_L(0)]> = (w/2pi) exp(-w^2 L^2)  [= pointwise w/2pi x form factor]",
            sp.simplify(rho_A - w*sp.exp(-w**2*L**2)/(2*sp.pi)) == 0, f"rho_L = {rho_A}")

print("-- A.3 compact support: top hat of radius R, chi vanishes for t > 2R --")
Vlens = sp.pi*(4*R + d)*(2*R - d)**2/12                         # overlap volume of two balls at distance d<2R
ff_top = Vlens*(3/(4*sp.pi*R**3))**2                             # (f_R * f_R)(d)
report_bool("int d^3x (f_R*f_R) = 1 (lens-volume formula correctly normalised)",
            sp.simplify(sp.integrate(4*sp.pi*d**2*ff_top, (d, 0, 2*R)) - 1) == 0)
selfE = sp.simplify(sp.integrate(d*ff_top, (d, 0, 2*R)))
report_bool("int_0^{2R} chi_R dt = int_0^{2R} t (f_R*f_R)(t) dt = 3/(10 pi R)", sp.simplify(selfE - 3/(10*sp.pi*R)) == 0,
            f"= {selfE}")
mp.mp.dps = 15
Rv = mp.mpf(1)
fT = lambda rr: 3/(4*mp.pi*Rv**3) if rr < Rv else mp.mpf(0)
def conv_top(dv):
    # mu-integral done exactly: f(sqrt(r^2+d^2-2rd mu)) != 0 iff mu > (r^2+d^2-R^2)/(2 r d)
    def inner(rr):
        mu0 = (rr**2 + dv**2 - Rv**2)/(2*rr*dv)
        if mu0 >= 1: return mp.mpf(0)
        return (1 - max(mu0, mp.mpf(-1)))*3/(4*mp.pi*Rv**3)
    br = sorted({mp.mpf(0), abs(Rv - dv), Rv} if abs(Rv-dv) < Rv else {mp.mpf(0), Rv})
    return 2*mp.pi*3/(4*mp.pi*Rv**3)*mp.quad(lambda rr: rr**2*inner(rr), br)
print("     d/R     numeric (f*f)(d)      analytic lens formula     chi_R(t=d) = d (f*f)(d)")
for dv in ['0.5', '1.0', '1.5', '1.9', '1.99', '2.0', '2.1', '3.0']:
    dv = mp.mpf(dv)
    num = conv_top(dv)
    an = (mp.pi*(4*Rv+dv)*(2*Rv-dv)**2/12)*(3/(4*mp.pi*Rv**3))**2 if dv < 2*Rv else mp.mpf(0)
    print(f"    {mp.nstr(dv,4):>5}   {mp.nstr(num,12):>18}   {mp.nstr(an,12):>18}   {mp.nstr(dv*num,12):>18}")
    if dv < 2*Rv: report(f"top hat (f*f)(d={dv}R) vs lens volume", num, an, tol=1e-9)
    else:         report(f"top hat chi_R(t={dv}R) = 0 (t > 2R)", dv*num, 0, tol=1e-30, absolute=True)
mp.mp.dps = 30

# ====================================================================================================
print("\n== B. conformal scalar stress tensor: smeared response kernel ==")
print("-- B.0 conformal structure of <T00 T00> at a generic point (exact Wick contraction, xi=1/6) --")
tE, x1, x2, x3 = sp.symbols('tE x1 x2 x3', real=True)
X = [tE, x1, x2, x3]
G = 1/(4*sp.pi**2*(tE**2 + x1**2 + x2**2 + x3**2))
xi = sp.Rational(1, 6)
def T_terms(mu, nu):
    terms = [(sp.Integer(1), (mu,), (nu,))]
    if mu == nu:
        for al in range(4):
            terms += [(-sp.Rational(1, 2), (al,), (al,)), (2*xi, (al,), (al,)), (2*xi, (), (al, al))]
    terms += [(-2*xi, (mu,), (nu,)), (-2*xi, (), (mu, nu))]
    return terms
def Dx(expr, idx):
    for i in idx: expr = sp.diff(expr, X[i])
    return expr
_cache = {}
def corr(a, bidx):
    key = tuple(sorted(a + bidx)), len(bidx) % 2
    if key not in _cache: _cache[key] = (-1)**len(bidx)*Dx(G, a + bidx)
    return _cache[key]
def two_point(A, B):
    return sum(ca*cb*(corr(a1, b1)*corr(a2, b2) + corr(a1, b2)*corr(a2, b1)) for ca, a1, a2 in A for cb, b1, b2 in B)
full = two_point(T_terms(0, 0), T_terms(0, 0))
xE2 = tE**2 + x1**2 + x2**2 + x3**2
CT = 1/(3*sp.pi**4)
I00E = 1 - 2*tE**2/xE2
report_bool("<T00(x)T00(0)>_Wick = C_T (I_00^2 - 1/4)/(x^2)^4, C_T = 1/(3 pi^4), I_00 = 1 - 2 tau^2/x^2, all x",
            sp.simplify(full - CT*(I00E**2 - sp.Rational(1, 4))/xE2**4) == 0)
report_bool("time axis: I_00 = 1, I_{00,00} = 3/4, <T00 T00> = 1/(4 pi^4 t^8) (v1)",
            sp.simplify(full.subs({x1: 0, x2: 0, x3: 0}) - 1/(4*sp.pi**4*tE**8)) == 0)

print("-- B.1 closed form of chi_L(t) two ways (sympy) --")
u = sp.symbols('u', positive=True)
g_u = 2*sp.pi*sp.sqrt(u)*sp.exp(-u/(4*L**2))/(4*sp.pi*L**2)**sp.Rational(3, 2)      # 4 pi r^2 f dr = g(u) du, u=r^2
# W = C_T int g(u) [3/(4a^4) - 4u/a^5 + 4u^2/a^6] du, a = u - t^2 + i0 (t>0);  Im 1/(a+i0)^n -> -pi(-1)^{n-1}delta^{(n-1)}/(n-1)!
chi_pos = 2*sp.pi*CT*(sp.Rational(3, 4)*sp.diff(g_u, u, 3)/6 - sp.diff(u*g_u, u, 4)/6 + sp.diff(u**2*g_u, u, 5)/30)
chi_pos = sp.simplify(chi_pos.subs(u, t**2))
K, a = sp.symbols('K a', positive=True)
I0 = sp.integrate(sp.exp(-a*K**2)*sp.cos(K*t), (K, 0, sp.oo))
I6 = -sp.diff(I0, a, 3)                                                             # int_0^oo K^6 e^{-aK^2} cos(Kt) dK
chi_mom = sp.simplify(I6.subs(a, L**2)/(1440*sp.pi**4*t))                           # (1/pi) int rho_L sin(wt), rho_L = int_0^w K^6 e^{-K^2L^2}/(1440 pi^3)
print(f"  chi_L(t) = {sp.factor(chi_pos)}")
report_bool("position-space (distributional) closed form == momentum-space (Ward/CFT) closed form",
            sp.simplify(chi_pos - chi_mom) == 0)
ser = sp.series(chi_pos, t, 0, 2).removeO()
lead = sp.simplify(sp.limit(chi_pos*t, t, 0))
rho_inf = sp.sqrt(sp.pi)/(1536*sp.pi**3*L**7)                                        # (1/(1440 pi^3)) int_0^oo K^6 e^{-K^2L^2} dK
report_bool("small-t edge: t chi_L(t) -> 1/(1536 pi^{7/2} L^7) = rho_L(oo)/pi", sp.simplify(lead - rho_inf/sp.pi) == 0,
            f"series: {ser}")
report_bool("rho_L(oo) = (1/(1440 pi^3)) int_0^oo K^6 e^{-K^2 L^2} dK = sqrt(pi)/(1536 pi^3 L^7)",
            sp.simplify(sp.integrate(K**6*sp.exp(-K**2*L**2), (K, 0, sp.oo))/(1440*sp.pi**3) - rho_inf) == 0)
chi_cf = sp.lambdify((t, L), chi_pos, 'mpmath')

print("-- B.2 quadrature on the rotated contour r = s e^{i pi/8} (eps -> 0 exactly) vs closed form --")
CTm = 1/(3*mp.pi**4)
def W_pt(tt, rr):
    x2 = rr**2 - tt**2
    return CTm*((1 + 2*tt**2/x2)**2 - mp.mpf(1)/4)/x2**4
def fs2(rr, Lv): return mp.e**(-rr**2/(4*Lv**2))/(4*mp.pi*Lv**2)**mp.mpf(1.5)
def W_L_rot(tt, Lv, theta=mp.pi/8, eps=0):
    ph = mp.e**(1j*theta)
    tc = tt - 1j*eps
    f = lambda s: ph*4*mp.pi*(s*ph)**2*fs2(s*ph, Lv)*W_pt(tc, s*ph)
    return mp.quad(f, [0, tt/2, tt, 2*tt, 2*tt + 4*Lv, mp.inf])
def W_L_real(tt, Lv, eps):
    tc = tt - 1j*eps
    f = lambda rr: 4*mp.pi*rr**2*fs2(rr, Lv)*W_pt(tc, rr)
    pts = [0] + [tt + k*eps for k in (-6, -3, -1, 0, 1, 3, 6)] + [2*tt + 4*Lv, mp.inf]
    return mp.quad(f, pts)
Lv = mp.mpf(1)
Wa, Wb = W_L_rot(mp.mpf('0.8'), Lv, eps=mp.mpf('0.2')), W_L_real(mp.mpf('0.8'), Lv, mp.mpf('0.2'))
report("contour rotation valid at finite eps=0.2L, t=0.8L: Re W (rotated vs real axis)", mp.re(Wa), mp.re(Wb), tol=1e-8)
report("contour rotation valid at finite eps=0.2L, t=0.8L: Im W (rotated vs real axis)", mp.im(Wa), mp.im(Wb), tol=1e-8)
print("     t/L      L^8 chi_L (quadrature)      L^8 chi_L (closed form)    L^8 Re W_L     chi_L(t)/chi_L(L)")
tab = {}
for x in ['0.05', '0.25', '0.5', '0.75', '1', '1.5', '2', '3', '4']:
    tt = mp.mpf(x)*Lv
    W = W_L_rot(tt, Lv)
    chi_num, chi_an = -2*mp.im(W), chi_cf(tt, Lv)
    tab[x] = chi_num
    print(f"    {x:>5}   {mp.nstr(Lv**8*chi_num,16):>24}   {mp.nstr(Lv**8*chi_an,16):>24}   {mp.nstr(Lv**8*mp.re(W),8):>12}")
    report(f"chi_L(t={x} L) quadrature vs closed form", chi_num, chi_an, tol=1e-12)
for x in ['2', '3', '4']:
    print(f"    chi_L({x}L)/chi_L(L) = {mp.nstr(tab[x]/tab['1'],6)};  Gaussian tail exp(-t^2/4L^2) x polynomial")
report("edge: t chi_L(t) at t = 0.05 L vs rho_L(oo)/pi (1/t law, 1 - 7 t^2/(4 L^2) correction)",
       mp.mpf('0.05')*tab['0.05'], (1/(1536*mp.pi**3.5))*(1 - 7*mp.mpf('0.05')**2/4), tol=2e-5)

print("-- B.3 scaling chi_L(t) = L^{-8} g(t/L) (numerical, L = 1/2, 1, 2 at t/L = 0.5 and 2) --")
for x in ['0.5', '2']:
    vals = [Lx**8*(-2*mp.im(W_L_rot(mp.mpf(x)*Lx, Lx))) for Lx in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2))]
    report(f"L^8 chi_L at t/L={x}: L=1/2 vs L=1", vals[0], vals[1], tol=1e-12)
    report(f"L^8 chi_L at t/L={x}: L=2 vs L=1", vals[2], vals[1], tol=1e-12)
print("  (F_L = int d^3x f_L T00 has mass dimension 4 => <F_L F_L> ~ L^-8; the task text's L^-9 is corrected to L^-8)")

print("-- B.4 Fourier transform of the smeared Wightman function at T=0 (shifted line Im t = -L, numpy) --")
CTf = 1/(3*np.pi**4)
def W_pt_np(tt, rr):
    x2 = rr**2 - tt**2
    return CTf*((1 + 2*tt**2/x2)**2 - 0.25)/x2**4
def WL_line(s, Lf, delta, nr=800, rmax=14.0):
    rr, wr = gl(0.0, rmax*Lf, nr)
    f = np.exp(-rr**2/(4*Lf**2))/(4*np.pi*Lf**2)**1.5
    return np.sum(wr*4*np.pi*rr**2*f*W_pt_np(s[:, None] - 1j*delta, rr[None, :]), axis=1)
def Wtilde(omegas, Lf=1.0, delta=1.0, ns=4000, smax=200.0, nr=800):
    s, ws = gl(-smax*Lf, smax*Lf, ns)
    WL = WL_line(s, Lf, delta, nr=nr)
    return np.array([np.exp(om*delta)*np.sum(ws*np.exp(1j*om*s)*WL) for om in omegas])
def rho_ward(om, Lf=1.0):
    return (1/(1440*mp.pi**3))*mp.quad(lambda k: k**6*mp.e**(-k**2*Lf**2), [0, om])
oms = [0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0]
Wt = Wtilde(oms + [-o for o in oms])
print("      w L      W~_L(w) [L^7]        rho_L^{(b)}(w)        |W~_L(-w)|/W~_L(w)    W~_L/(w^7/10080 pi^3)")
for i, o in enumerate(oms):
    ex = rho_ward(o)
    print(f"    {o:>5}   {Wt[i].real:.10e}   {mp.nstr(ex,11):>14}   {abs(Wt[i+len(oms)])/abs(Wt[i]):.2e}            "
          f"{mp.nstr(Wt[i].real/(o**7/(10080*mp.pi**3)),8)}")
    report(f"W~_L(w={o}/L) = rho_L^{{(b)}}(w) (Ward/CFT form)", Wt[i].real, ex, tol=1e-10)
    report_bool(f"positive frequency: |W~_L(-{o}/L)|/W~_L({o}/L) < 1e-9", abs(Wt[i+len(oms)])/abs(Wt[i]) < 1e-9,
                f"{abs(Wt[i+len(oms)])/abs(Wt[i]):.2e}")
report("small-w limit W~_L(0.25/L)/(w^7/10080 pi^3) -> 1 - (7/9)(wL)^2 (pointwise v1 spectrum recovered)",
       Wt[0].real/(0.25**7/(10080*mp.pi**3)), 1 - mp.mpf(7)/9*0.25**2, tol=2e-3)
report("large-w plateau W~_L(6/L) -> sqrt(pi)/(1536 pi^3 L^7)", Wt[len(oms)-1].real, mp.sqrt(mp.pi)/(1536*mp.pi**3), tol=1e-6)

# ====================================================================================================
print("\n== C. KMS for smeared observables: state x commutator = noise, for any smearing ==")
print("-- C.1 image sum closed form --")
beta_m, tt_c, r_c = 2*mp.pi, mp.mpf('0.7') - 0.3j, mp.mpf('1.3')
raw = mp.nsum(lambda n: 1/(4*mp.pi**2*(r_c**2 - (tt_c + 1j*n*beta_m)**2)), [-mp.inf, mp.inf])
cf = (1/(8*mp.pi*r_c*beta_m))*(mp.coth(mp.pi*(r_c - tt_c)/beta_m) + mp.coth(mp.pi*(r_c + tt_c)/beta_m))
report("sum_n 1/(4pi^2 (r^2-(t+i n beta)^2)) = [coth(pi(r-t)/beta)+coth(pi(r+t)/beta)]/(8 pi r beta) at a complex point",
       abs(raw - cf), 0, tol=1e-25, absolute=True)

print("-- C.2 smeared thermal Wightman spectrum, shifted line Im t = -beta/4 (generic: not the KMS-symmetric beta/2) --")
def coth_np(z): return 1/np.tanh(z)
def GL_line(s, Lf, beta, delta, nr=700, rmax=14.0):
    rr, wr = gl(0.0, rmax*Lf, nr)
    f = np.exp(-rr**2/(4*Lf**2))/(4*np.pi*Lf**2)**1.5
    tt = s[:, None] - 1j*delta
    integ = 4*np.pi*rr**2*f*(1/(8*np.pi*rr*beta))*(coth_np(np.pi*(rr - tt)/beta) + coth_np(np.pi*(rr + tt)/beta))
    return np.sum(wr*integ, axis=1)
def Gplus(omegas, Lf, beta, delta, ns=5000, nr=700):
    smax = 25*beta + 20*Lf
    s, ws = gl(-smax, smax, ns)
    Gs = GL_line(s, Lf, beta, delta, nr=nr)
    return np.array([np.exp(om*delta)*np.sum(ws*np.exp(1j*om*s)*Gs) for om in omegas])
for beta, Lf, oms_c in [(2*np.pi, 1.0, [0.25, 0.5, 1.0, 2.0]), (2*np.pi, 0.5, [0.25, 0.5, 1.0, 2.0]),
                        (1.0, 0.3, [0.25, 0.5, 1.0, 2.0, 3.0])]:
    Gp = Gplus(np.array(oms_c + [-o for o in oms_c]), Lf, beta, beta/4)
    n = len(oms_c)
    print(f"  beta = {beta:.6f}, L = {Lf}:")
    print("      w        G+_L(w)            G+_L(-w)          G+/G- e^-bw - 1    rho_L/rho_A - 1    N_L/(coth rho/2) - 1")
    for i, o in enumerate(oms_c):
        gp, gm = Gp[i].real, Gp[i+n].real
        rhoA = o*np.exp(-o**2*Lf**2)/(2*np.pi)
        r1, r2, r3 = gp/gm/np.exp(beta*o) - 1, (gp - gm)/rhoA - 1, 0.5*(gp + gm)/(0.5/np.tanh(beta*o/2)*rhoA) - 1
        print(f"    {o:>5}   {gp:.10e}   {gm:.6e}   {r1:+.2e}          {r2:+.2e}          {r3:+.2e}")
        report_bool(f"KMS  G+_L({o})/G+_L({-o}) = e^(beta w)  [beta={beta:.4g}, L={Lf}]", abs(r1) < 1e-6, f"{r1:.2e}")
        report_bool(f"rho_L({o}) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A)", abs(r2) < 1e-6, f"{r2:.2e}")
        report_bool(f"N_L({o}) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L", abs(r3) < 1e-6, f"{r3:.2e}")
        report_bool(f"Im G+_L({o}) = 0", abs(Gp[i].imag) < 1e-12*abs(gp), f"{abs(Gp[i].imag):.1e}")

print(f"\n  (runtime {time.time()-T0:.1f} s)")
print("V7 OVERALL:", "PASS" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
