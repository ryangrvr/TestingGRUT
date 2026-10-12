#!/usr/bin/env python3
"""
V11 (revision 3) -- The conventional KMS structure of the comoving conformal vacuum in FRW: exact first-order thermality,
the true second-order defect, and the present-day (non-perturbative) baseline.

CORRECTION OF REVISION 2.  The earlier version of this script reported a first-order temperature-shift coefficient
c ~ 1.014 and a '1%-of-epsilon frequency dependence' from finite-epsilon tables (p = 25..200) without extrapolating.
Two independent reviewers showed analytically that at fixed x = H_eff*Delta the exact power-law expression
    a1 a2 (eta1-eta2)^2 = (tbar/q)^2 (1-u^2/4)^(1-q) [(1+u/2)^q - (1-u/2)^q]^2 ,   q = p-1,  u = Delta/tbar,
has NO 1/q term: to first order the state is EXACTLY de Sitter-thermal at H_eff = q/tbar = Hbar(1-eps) = a''/a',
and the first non-thermal terms are O(eps^2).  Richardson extrapolation of the old numbers gives c = 1 to 3e-5.
Revision 2's headline was finite-epsilon contamination.  This file proves the corrected statements and computes
what is actually true.

 A. (sympy, exact) Power law: ln R := ln[a1 a2 (Delta eta)^2 / (4 sinh^2(x/2)/H_eff^2)] = [-x^2/4 + (x^3/12) coth(x/2)]/q^2 + O(q^-3)
    => W_eps / W_dS(H_eff) = 1 + eps^2 [x^2/4 - (x^3/12) coth(x/2)] + O(eps^3)   (reviewers' coefficient confirmed).
 B. (sympy, exact) GENERAL a(t) to second adiabatic order, with eps = -Hdot/H^2 and eta = eps_dot/(H eps):
    the first-order correction is a pure shift H -> H_eff = H + Hdot/H = H(1-eps) = (d/dt) ln(a_dot)  (THEOREM, any a(t) with eps AND eta = eps_dot/(H eps) first-order small; the eps*eta term below is why);
    the second-order non-thermal remainder depends on eps^2 and eps*eta with explicit coefficient functions of x (complete at
    that order only if the third flow parameter eps3 = eta_dot/(H eta) is also first-order small -- the full adiabatic hierarchy;
    the eps*eta coefficient is checked below against its closed form g(x) = 2 + x^2/4 - (x + x^3/24) coth(x/2)).
    SCOPE: the conformally coupled massless scalar's worldline Wightman function in spatially flat FRW (other scalar conformal
    primaries follow by Weyl covariance; composite, higher-spin and smeared observables are not computed).
 C. (mpmath) Numerics with Richardson extrapolation: D/eps -> 2 pi w/Hbar exactly (c = 1); the TRUE second-order
    KMS defect D2(w) := lim [ln(G+(w)/G+(-w)) - 2 pi w/H_eff]/eps^2, tabulated; its non-linearity in w is the
    convention-free violation.  It is O(eps^2), not O(eps).
 D. (mpmath) The present universe is NOT slow-roll: eps_0 = (3/2) Omega_m = 0.4725 (Omega_m = 0.315), eta_0 = -3 Omega_L = -2.06.
    Exact worldline Wightman function of a conformal scalar in flat LambdaCDM, a(t) = (Om/OL)^(1/3) sinh^(2/3)(3 sqrt(OL) H0 t/2),
    for pairs symmetric about t0, compared in the TIME domain with the first-order thermal form at H_eff(t0) = a''/a'|_0 and
    with the T = 0 (Minkowski) form.  No Fourier spectrum and no KMS test: the comoving worldline exists only for |Delta| < 2 t0.
    What this part shows (audit, 11_VERIFICATION_ROUND.md): the present-day kernel is not described by the first-order thermal
    form at H_eff -- the small-Delta agreement is the universal Hadamard term (the local effective rate in the Delta^2 correction
    is H0 sqrt(1 - 2 eps_0) = 0.23 H0, not H_eff = 0.53 H0), the T = 0 form is closer at every grid point, and the ratio to the
    truncated second-order form depends on an O(eps^3) convention (x = H_eff Delta vs x = H Delta) by a factor 1/(1-eps_0)^2 = 3.6,
    so the eps-expansion is not quantitatively usable at eps_0 = 0.47.  The pairs are symmetric about t0, so Delta = 1.5/H0 pairs
    z = 2.25 (eps = 1.4 there) with a point 10.9 Gyr in the future; the deviation grows without bound as Delta -> 2 t0.
    (Formal baseline for the vacuum medium's KMS structure; physically the CMB bath at 2.7 K dominates real photon correlators.)
Conclusions: the conformally coupled scalar's worldline function stays KMS through O(eps) (fixed-mean-time convention; omega >> eps H)
-- the FDT boundary of 07 applies through first order; the first conventional KMS violation is O(eps^2); the present-day comoving
vacuum is not described by the first-order thermal form, and the adiabatic expansion is not usable at eps_0 ~ 0.5, eta_0 ~ -2.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
ok_all = True
def report(name, got, want, tol):
    global ok_all
    err = abs(mp.mpf(got) - mp.mpf(want)); rel = err/abs(mp.mpf(want)) if want != 0 else err
    ok = rel < tol; ok_all &= bool(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {mp.nstr(got, 10)} want {mp.nstr(want, 10)} rel {mp.nstr(rel, 3)}")

print("== A. power law: exact first-order thermality and the second-order coefficient (sympy, explicit expansion) ==")
x, sq = sp.symbols('x s', positive=True)        # s = 1/q
# E_pm = (1/s) ln(1 +- x s/2) expanded to O(s^3):  +-x/2 - x^2 s/8 +- x^3 s^2/24 - x^4 s^3/64
Ep = x/2 - x**2*sq/8 + x**3*sq**2/24 - x**4*sq**3/64
Em = -x/2 - x**2*sq/8 - x**3*sq**2/24 - x**4*sq**3/64
# (1+u/2)^q - (1-u/2)^q = e^{Ep} - e^{Em} = 2 e^{(Ep+Em)/2} sinh((Ep-Em)/2)
lnR = (1 - 1/sq)*sp.log(1 - x**2*sq**2/4) + 2*((Ep+Em)/2 + sp.log(2*sp.sinh((Ep-Em)/2))) - sp.log(4*sp.sinh(x/2)**2)
ser = sp.series(lnR, sq, 0, 3).removeO()
c0 = sp.simplify(ser.coeff(sq, 0)); c1 = sp.simplify(ser.coeff(sq, 1)); c2 = sp.simplify(ser.coeff(sq, 2))
print("  O(1):", c0, "  O(1/q):", c1, "->", "PASS (no first-order term)" if (c0 == 0 and c1 == 0) else "FAIL")
ok_all &= (c0 == 0 and c1 == 0)
target_c2 = -x**2/4 + x**3*sp.coth(x/2)/12
d = sp.simplify((c2 - target_c2).rewrite(sp.exp))
print("  O(1/q^2):", c2, " vs -x^2/4 + (x^3/12)coth(x/2):", "PASS" if d == 0 else f"FAIL ({d})")
ok_all &= (d == 0)
# independent high-precision numerical check of the exact expression (no expansion at all)
def lnR_exact(xv, qv):
    xv = mp.mpf(xv); qv = mp.mpf(qv); uu = xv/qv
    return (1-qv)*mp.log(1-uu**2/4) + 2*mp.log((1+uu/2)**qv - (1-uu/2)**qv) - mp.log(4*mp.sinh(xv/2)**2)
for xv in ['0.7', '2.0']:
    v = lnR_exact(xv, 10**5)*mp.mpf(10**5)**2
    tgt = -mp.mpf(xv)**2/4 + mp.mpf(xv)**3*mp.coth(mp.mpf(xv)/2)/12
    report(f"  exact ln R * q^2 at q=1e5, x={xv}", v, tgt, tol=1e-4)
print("  => W_eps/W_dS(H_eff) = 1 + eps^2 [x^2/4 - (x^3/12) coth(x/2)] + O(eps^3), H_eff = Hbar(1-eps).  Reviewers confirmed.")
f2 = x**2/4 - x**3*sp.coth(x/2)/12
kms_break = sp.simplify((f2.subs(x, x - 2*sp.pi*sp.I) - f2.subs(x, -x)).rewrite(sp.exp))
print("  f2(x - 2 pi i) - f2(-x) =", kms_break, " (nonzero => KMS broken at O(eps^2))")
ok_all &= (kms_break != 0)

print("== B. general a(t) to second adiabatic order (sympy construction + numerical verification) ==")
H, eps, eta, D, dd, lam = sp.symbols('H epsilon eta Delta delta lambda', positive=True)
Hddot = eps*H**3*(2*eps - eta)                       # d^2H/dt^2 with eps=-Hdot/H^2, eta=eps_dot/(H eps); BOTH eps and eta are first-order small
lna = H*dd - (eps*H**2/2)*dd**2 + (Hddot/6)*dd**3    # ln a(tbar+delta) - ln a(tbar)
expo = sp.expand(((eps*H**2/2)*dd**2 - (Hddot/6)*dd**3).subs({eps: lam*eps, eta: lam*eta}))
pert = sp.expand(1 + expo + expo**2/2)
pert = sum(pert.coeff(lam, k)*lam**k for k in range(3))
deta = sp.integrate(sp.expand(pert*sp.exp(-H*dd)), (dd, -D/2, D/2))        # conformal-time difference (abar = 1)
aa = sp.exp((lna.subs(dd, D/2) + lna.subs(dd, -D/2)).subs({eps: lam*eps, eta: lam*eta}))
ratio = sp.series(sp.expand(aa*deta**2), lam, 0, 3).removeO()
Heff = H*(1 - lam*eps)
thermal = sp.series(4*sp.sinh(Heff*D/2)**2/Heff**2, lam, 0, 3).removeO()
diff = sp.expand(ratio - thermal)
first = diff.coeff(lam, 1); second = diff.coeff(lam, 2)
first_sym = sp.simplify(first)
print("  first-order remainder after the shift H -> H(1 - eps), simplified symbolically:", first_sym, "->", "PASS (identically zero)" if first_sym == 0 else "FAIL")
ok_all &= (first_sym == 0)
# numerical verification at random points (avoids heavy symbolic simplification)
import random
random.seed(3)
mx1 = 0; mx2 = 0
for _ in range(6):
    vals = {H: random.uniform(0.5, 2.0), D: random.uniform(0.3, 3.0), eps: random.uniform(0.1, 0.9), eta: random.uniform(-1, 1)}
    f1 = abs(sp.N(first.subs(vals), 30)); mx1 = max(mx1, f1)
    # power-law check: eta = 0 must give second/(4 sinh^2(HD/2)/H^2) = -eps^2 f2(x)... (ratio of a a' deta^2 to thermal form is R = 1 + c2/q^2, c2 = -x^2/4 + (x^3/12)coth(x/2))
    v0 = dict(vals); v0[eta] = 0
    sec = sp.N((second/(4*sp.sinh(H*D/2)**2/H**2)).subs(v0), 30)
    xv = v0[H]*v0[D]
    tgt = v0[eps]**2*(-xv**2/4 + xv**3/mp.tanh(xv/2)/12)
    mx2 = max(mx2, abs(float(sec) - float(tgt)))
print(f"  max |first-order difference from thermal form at H_eff = H(1-eps)| over random points: {mx1:.2e} (float inputs; zero to machine precision) ->", "PASS" if mx1 < 1e-12 else "FAIL")
ok_all &= (mx1 < 1e-12)
print(f"  max |second-order (eta=0) - power-law coefficient| over random points: {mx2:.2e} ->", "PASS" if mx2 < 1e-10 else "FAIL")
ok_all &= (mx2 < 1e-10)
sec_eta = sp.simplify(sp.expand(second).coeff(eta, 1)/(4*sp.sinh(H*D/2)**2/H**2))
print("  eps*eta term of the second-order relative remainder (general a(t)):", sp.simplify(sec_eta.subs(D, x/H)))
g_closed = lambda xv: 2 + xv**2/4 - (xv + xv**3/24)*mp.coth(xv/2)
mxg = max(abs(mp.mpf(sp.N((sec_eta/eps).subs({D: xv, H: 1}), 30)) - g_closed(mp.mpf(xv))) for xv in ['0.5', '1', '2', '3'])   # sec_eta = eps * g(x)
print(f"  eps*eta coefficient vs closed form g(x) = 2 + x^2/4 - (x + x^3/24) coth(x/2) at x = 0.5, 1, 2, 3: max diff {mp.nstr(mxg, 3)} ->", "PASS" if mxg < 1e-12 else "FAIL")
ok_all &= (mxg < 1e-12)
print("  THEOREM (any a(t) with eps and eta first-order small): to first adiabatic order the comoving conformal vacuum's fixed-mean-time")
print("  spectrum is thermal at T_eff = (H + Hdot/H)/2pi = (d/dt ln a_dot)/2pi (omega >> eps H; the eps^2 remainder is not uniform in Delta).")

print("== C. numerics: Richardson-extrapolated first- and second-order KMS defects (power law) ==")
def make_W(p, tbar):
    def a(t): return (t/tbar)**p
    def eta_c(t): return tbar**p * t**(1-p)/(1-p)
    def W(Delta):
        t1 = tbar + Delta/2; t2 = tbar - Delta/2
        return -1/(4*mp.pi**2 * a(t1)*a(t2) * (eta_c(t1)-eta_c(t2))**2)
    return W
def Gplus(W, Hshift, w, L=45):
    delta = mp.pi/Hshift
    f = lambda uu: mp.e**(1j*w*uu) * W(uu - 1j*delta)
    pts = [-L, -20, -8, -3, -1, 0, 1, 3, 8, 20, L]
    return mp.re(mp.e**(w*delta) * mp.quad(f, [xx/Hshift for xx in pts]))
Hb = mp.mpf(1); ws = [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(3)]
ps = [200, 400, 800, 1600]
D1 = {}; D2 = {}
for p in ps:
    e = 1/mp.mpf(p); tbar = p/Hb; W = make_W(p, tbar); Heff = (p-1)/tbar
    D1[p] = []; D2[p] = []
    for w in ws:
        gp = Gplus(W, Heff, w); gm = Gplus(W, Heff, -w)
        lr = mp.log(gp/gm)
        D1[p].append((lr - 2*mp.pi*w/Hb)/e)            # first-order defect relative to Hbar
        D2[p].append((lr - 2*mp.pi*w/Heff)/e**2)       # second-order defect relative to the exact H_eff
def richardson(seq):  # seq at eps, eps/2, eps/4, eps/8 ; error C1 eps + C2 eps^2 + C3 eps^3 -> eliminate successively with weights 2, 4/3, 8/7
    s = list(seq); k = 1                              # (rev. 3 used weight 2 at every level, leaving an eps^2/8 artefact ~3e-5; audit 11)
    while len(s) > 1:
        s = [(2**k*s[i+1] - s[i])/(2**k - 1) for i in range(len(s)-1)]; k += 1
    return s[0]
print("  D/eps relative to Hbar, Richardson over p =", ps, " (expected 2 pi w/Hbar, i.e. c = 1 exactly):")
for i, w in enumerate(ws):
    val = richardson([D1[p][i] for p in ps])
    report(f"   lim D/eps at w/H={mp.nstr(w,2)} = 2 pi w/H", val, 2*mp.pi*w/Hb, tol=1e-4)
print("  second-order defect D2(w) = lim [ln(G+/G-) - 2 pi w/H_eff]/eps^2  (the genuine non-thermal content):")
D2lim = []
for i, w in enumerate(ws):
    val = richardson([D2[p][i] for p in ps]); D2lim.append(val)
    print(f"    w/H = {mp.nstr(w,2)}:  D2 = {mp.nstr(val, 7)}   (p=1600 raw {mp.nstr(D2[1600][i],7)})   D2/(2 pi w/H) = {mp.nstr(val/(2*mp.pi*w/Hb), 6)}")
slope_var = D2lim[-1]/(2*mp.pi*ws[-1]/Hb) - D2lim[0]/(2*mp.pi*ws[0]/Hb)
print(f"  non-linearity of D2 in w (beta_eff variation per eps^2 between w=0.5H and 3H): {mp.nstr(slope_var, 6)}  -> frequency-dependent KMS violation is O(eps^2)")
ok_all &= abs(slope_var) > 1e-3

print("== D. the present universe: exact LambdaCDM worldline correlator vs the first-order thermal form (time domain) ==")
Om, OL = mp.mpf('0.315'), mp.mpf('0.685'); H0 = mp.mpf(1)
B = mp.mpf(3)/2*mp.sqrt(OL)*H0; A = (Om/OL)**(mp.mpf(1)/3)
a_L = lambda t: A*mp.sinh(B*t)**(mp.mpf(2)/3)
t0 = mp.asinh(mp.sqrt(OL/Om))/B                       # a(t0) = 1
eps0 = mp.mpf(3)/2*Om; Heff0 = H0*(1 - eps0)
print(f"  t0 = {mp.nstr(t0,6)}/H0 (a worldline exists only for |Delta| < 2 t0: no full Fourier spectrum at t0 without a window)")
eta0 = -3*OL
print(f"  eps_0 = {mp.nstr(eps0,4)}, eta_0 = -3 Omega_L = {mp.nstr(eta0,4)}, H_eff(t0) = H0 + Hdot/H0 = a''/a' = {mp.nstr(Heff0,4)} H0 -> T_eff = H_eff/2pi = {mp.nstr(Heff0/(2*mp.pi),4)} H0")
print(f"  local effective rate in the Delta^2 term of a1 a2 deta^2 = Delta^2 [1 + H^2 (1 - 2 eps) Delta^2/12 + ...]: H0 sqrt(1 - 2 eps_0) = {mp.nstr(mp.sqrt(1 - 2*eps0),4)} H0 (vs H_eff = {mp.nstr(Heff0,4)} H0)")
Hnum = lambda t: mp.diff(lambda s_: mp.log(a_L(s_)), t)
report("  eps_0 from a(t) directly", -mp.diff(Hnum, t0)/Hnum(t0)**2, eps0, tol=1e-8)
def aadeta2(Delta):
    t1 = t0 + Delta/2; t2 = t0 - Delta/2
    de = mp.quad(lambda t: 1/a_L(t), [t2, t1])
    return a_L(t1)*a_L(t2)*de**2
print("  Delta*H0   W_dS(H_eff)/W_exact - 1   W_Mink/W_exact - 1   ratio to eps0^2 f2 at x = H_eff Delta   at x = H0 Delta")
devs = []
for Dv in ['0.25', '0.5', '1.0', '1.5']:
    Dv = mp.mpf(Dv)
    exact = aadeta2(Dv)
    therm = 4*mp.sinh(Heff0*Dv/2)**2/Heff0**2
    dev = therm/exact - 1; dev0 = Dv**2/exact - 1
    f2 = lambda xv: xv**2/4 - xv**3*mp.coth(xv/2)/12
    devs.append((Dv, dev, dev0))
    print(f"   {mp.nstr(Dv,3):>5}      {mp.nstr(dev, 5):>12}          {mp.nstr(dev0, 5):>12}            {mp.nstr(dev/(eps0**2*f2(Heff0*Dv)), 4):>8}                      {mp.nstr(dev/(eps0**2*f2(H0*Dv)), 4)}")
print(f"  => deviation from the first-order thermal form at H_eff: {mp.nstr(100*devs[0][1],2)}% (Delta = 0.25/H0) to {mp.nstr(100*devs[-1][1],3)}% (1.5/H0);")
print(f"     the T = 0 (Minkowski) form is closer at every grid point ({mp.nstr(100*abs(devs[0][2]),2)}% to {mp.nstr(100*abs(devs[-1][2]),3)}%): the small-Delta agreement is the")
print("     universal Hadamard term, not thermality at H_eff. The ratio to the truncated second-order form depends on the x convention by")
print("     1/(1-eps_0)^2 = 3.6, so the eps-expansion is not quantitatively usable today (eps_0 = 0.47, eta_0 = -2.06). The exact kernel is needed.")
print("     (Pairs are symmetric about t0: Delta = 1.5/H0 pairs z = 2.25 with a point 10.9 Gyr in the future; the deviation is unbounded as Delta -> 2 t0.)")
print("  (formal baseline for the vacuum medium; the physical photon field is dominated by the 2.7 K CMB bath at these scales)")

print("\nV11 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
