#!/usr/bin/env python3
"""
V14 -- The massless minimally coupled sector in power-law FRW at first order in eps: the sigma-free computation.
(The decisive object the fifth-round reviewer asked for after v13's zero-mode coefficient turned out to depend on the
stochastic coarse-graining scale sigma.)

QUESTION.  Along a comoving worldline in power-law inflation a = (1 + eps H0 t)^(1/eps), is the Wightman function of a
massless minimally coupled scalar KMS at first order in eps, and at which temperature?  CONVENTION (sixth-round audit): the
statement is about the fixed-mean-time (Wigner) spectrum, W(tbar, Delta) Fourier-transformed in Delta at omega != 0 -- the
same convention as v11/D11 for the conformal sector.  An instantaneous one-sided Unruh-DeWitt rate at detector time tau uses
W(tau, tau - s) = W(tbar = tau - s/2, s) and differs at O(eps) by a frequency-dependent term from Hdot (it is not a pure
temperature shift in either sector); for THIS sector a linearly coupled detector is non-thermal already at eps = 0 because
of the secular zero-mode growth (Garbrecht-Prokopec 2004), which the Wigner convention places at omega = 0.  The conformal sector (v11) is
exactly thermal at T = (a''/a')/2pi = (H/2pi)(1 - eps) to O(eps).  Candidates for the massless sector: the same detector
rate (c = 1), the Kodama-Hayward apparent-horizon rate (c = 1/2), or neither (a frequency-dependent first-order defect).

METHOD (all closed form; every link checked numerically below).
 (L) v12 12.4 lemma: the power-law massless modes are the de Sitter modes of index nu = 3/2 + eps/(1-eps) =: 3/2 + d
     up to a ratio of scale factors.  Hence, with w = 1 + (tau - tau')^2/(4 tau tau'), tau = -1/((1-eps) a H):
        W(t,t') = (1-eps)^2 H(t) H(t') Gamma(3+d) Gamma(-d) 2F1(3+d, -d; 2; w) / (16 pi^2)  +  C_IR,
     where the hypergeometric is the Bunch-Davies de Sitter two-point function of index nu continued to nu > 3/2; the
     IR-cutoff Bunch-Davies state differs from it by the constant C_IR (frozen super-horizon modes).  [Part B checks
     this against a direct mode integral.]
 (A) Expansion around the pole at nu = 3/2:  2F1(3+d,-d;2;w) = 1 - d S1(w) - d^2 S2(w) + O(d^3) with
        S1 = w/(2(1-w)) - ln(1-w),
        S2 = -ln(1-w)/2 + Li2(w) + ln(1-w)/(2w) + 1/2 - (3/4) w/(1-w),
     and Gamma(3+d)Gamma(-d) = -(2/d)(1 + g1 d + g2 d^2 + ...), g1 = 3/2, g2 = 1/2 + pi^2/6.  Zeroth order:
     (H^2/16pi^2) 2 S1(w) = (H^2/8pi^2)[1/(1-Z) - ln(1-Z)] + const, the Allen-Folacci massless dS function.
 (C) First-order theorem (verified by Richardson in eps against the exact closed form):
        W(tbar, Delta) = (1 + g1 eps) W_mmc^dS(Delta; H_eff) + eps (Hbar^2/8pi^2) S2(cosh^2(Hbar Delta/2))
                         + c0 + c2 Delta^2 + O(eps^2),      H_eff = a''/a' = Hbar(1 - eps),  c2 = -(1-eps)^3 eps Hbar^4/(32 pi^2)
     (= -eps Hbar^4/32pi^2 at leading order; c0 depends on tbar -- it is the O(1) secular growth of the equal-time variance,
     (Hbar^3/4pi^2) tbar + C_IR -- and is Delta-independent at fixed tbar).
     The conformal-factor combination X/(w-1), X = (1-eps)^2 H(t)H(t'), is the v11 first-order theorem (pure shift to
     H_eff); the only new first-order term is the dilogarithm S2 evaluated on the zeroth-order w0 = cosh^2(Hbar Delta/2).
     Using Hbar or H_KH = Hbar(1 - eps/2) instead of H_eff, or dropping S2, leaves an O(eps) residual.
 (D) KMS.  Any function of w0 = cosh^2(Hbar Delta/2) is symmetric about the midline Im Delta = -beta/2, beta = 2pi/Hbar
     (w0(s - i beta/2) = -sinh^2(Hbar s/2) is even in s), i.e. KMS at beta to all orders by itself.  So the first-order
     KMS defect of the full function at beta_eff = 2pi/H_eff is O(eps^2): the massless sector is thermal at
     T = (a''/a')/2pi + O(eps^2), like the conformal sector; Kodama-Hayward is excluded at O(eps) (c2 Delta^2 and c0 are
     omega = 0 terms).  [14.D'' tests this on the EXACT closed form: the KMS defect with the omega = 0 polynomial removed scales as eps^2 at
     beta_eff and as eps at beta_KH and beta_bar.  14.D's midline checks are implementation checks only: any function of w0 is
     midline-symmetric by construction.]  The O(eps) physics is a SHAPE change of the spectral density, not a temperature change:
        rho(omega) = rho_mmc(omega; H_eff) [1 + eps r(omega)],  rho_mmc = (omega^2 + H^2)/(2 pi omega),
     with r(0) = 2 - 2 ln 2 = 0.614, r(Hbar) = 0.461, r -> 0 as omega -> inf (the S2 pole 3/(8pi^2 Delta^2) cancels the
     g1 amplitude: the Hadamard coefficient is unchanged, as it must be); large-omega asymptote r -> 3(2w^4 - w^2 - 1)/(4w^4(w^2+1))
     ~ 3/(2 w^2) from the short-distance expansion of S2 (checked at omega = 8).  The r table stops at omega = 8 because at
     mp.dps = 30 the e^(beta omega/2) factor amplifies the cancellation in h1 beyond that (precision, not physics; dps >= 50 fixes it).
 (E) Byproduct, sigma-free: the coefficient of |Delta| in the large-Delta expansion of W(tbar, Delta) (equivalently the
     1/omega^2 Brownian tail of the Wigner spectrum), W = const - (D/2)|Delta| - (3 eps Hbar^4/32 pi^2) Delta^2 + ..., the Delta^2
     term being the expansion of the secular variance growth (omega = 0).  Reading D as the zero-mode diffusion coefficient is
     the usual stochastic identification (an interpretive step; a Markovian walk with slowly varying D has no such Delta^2 term):
        D(tbar) = (Hbar^3/4pi^2) [1 - (1 + 2 ln 2) eps] + O(eps^2) = (H_eff^3/4pi^2)[1 + (2 - 2 ln 2) eps].
     The stochastic value A0(sigma) H^3 (v13) reproduces it only at sigma = exp(1 - gamma_E) = 1.526: the exact two-time
     function calibrates the coarse-graining scale for this observable.  (No gamma_E survives; one ln 2 does.)
Units: Hbar = H(tbar) = 1, a(tbar) = 1, tbar = 0.  Audit (11_VERIFICATION_ROUND.md): three lenses; mathematics and numerics
CONFIRMED (mutations of c2, g1 break 14.C; exact-strip KMS test passes; r asymptote matches), physics lens PARTIAL on the
convention point recorded above.  Part C's threshold is tuned to Im Delta = -0.3 (at Im Delta = -1 the ratios converge to 4 more slowly).
"""
import mpmath as mp, sympy as sp
ok_all = True
def check(label, cond):
    global ok_all
    print(f"  {label}:", "PASS" if cond else "FAIL"); ok_all &= bool(cond)
mp.mp.dps = 30
g1 = mp.mpf(3)/2; beta = 2*mp.pi; L2 = mp.log(2)
def S1(w): return w/(2*(1-w)) - mp.log(1-w)
def S2(w): return -mp.log(1-w)/2 + mp.polylog(2, w) + mp.log(1-w)/(2*w) + mp.mpf(1)/2 - mp.mpf(3)/4*w/(1-w)

print("== 14.A expansion of 2F1(3+d, -d; 2; w) around the massless pole ==")
for w in [mp.mpf('-0.7'), mp.mpf('0.3'), mp.mpc('0.5', '0.2'), mp.mpc('2.0', '-0.3')]:
    F = lambda d: mp.hyp2f1(3+d, -d, 2, w)
    d1 = mp.diff(F, 0, 1); d2 = mp.diff(F, 0, 2)/2
    check(f"w={mp.nstr(w,3)}: -F'(0) = S1 and -F''(0)/2 = S2", abs(-d1 - S1(w)) < 1e-15 and abs(-d2 - S2(w)) < 1e-15)
dd = sp.symbols('delta')
ser = sp.series(sp.gamma(3+dd)*sp.gamma(-dd)*(-dd/2), dd, 0, 3).removeO()
check("Gamma(3+d)Gamma(-d) = -(2/d)(1 + (3/2) d + (1/2 + pi^2/6) d^2 + ...)", sp.simplify(ser.coeff(dd, 1) - sp.Rational(3, 2)) == 0 and sp.simplify(ser.coeff(dd, 2) - (sp.Rational(1, 2) + sp.pi**2/6)) == 0)
Z = mp.mpf('-0.3'); w = (1+Z)/2; Zb = mp.mpf('0.6'); wb = (1+Zb)/2
check("zeroth order = Allen-Folacci (H^2/8pi^2)[1/(1-Z) - ln(1-Z)] up to a constant", abs((2*S1(w) - (2/(1-Z) - 2*mp.log(1-Z))) - (2*S1(wb) - (2/(1-Zb) - 2*mp.log(1-Zb)))) < 1e-25)

print("== 14.B the closed form vs a direct Bunch-Davies mode integral ==")
def setup(eps):
    nu = (3-eps)/(2*(1-eps)); d = nu - mp.mpf(3)/2
    a = lambda t: (1+eps*t)**(1/eps); H = lambda t: 1/(1+eps*t); tau = lambda t: -1/((1-eps)*a(t)*H(t))
    return nu, d, a, H, tau
def W_closed(eps, t1, t2):
    nu, d, a, H, tau = setup(eps)
    w = 1 + (tau(t1)-tau(t2))**2/(4*tau(t1)*tau(t2))
    return (1-eps)**2*H(t1)*H(t2)*mp.gamma(3+d)*mp.gamma(-d)*mp.hyp2f1(3+d, -d, 2, w)/(16*mp.pi**2)
def mode_integrand(eps, t1, t2):
    nu, d, a, H, tau = setup(eps); T1, T2 = tau(t1), tau(t2)
    a_of_tau = lambda T: (-(1-eps)*T)**(-1/(1-eps))
    return lambda k: k**2*(mp.pi/4)*mp.sqrt(-T1)*mp.sqrt(-T2)*mp.hankel1(nu, -k*T1)*mp.hankel2(nu, -k*T2)/(a_of_tau(T1)*a_of_tau(T2))/(2*mp.pi**2)
def W_modes(eps, t1, t2):
    return mp.quad(mode_integrand(eps, t1, t2), [0, 0.5, 2, 8, 30, 120, 500, mp.inf])
def W_modes_diff(eps, pair1, pair2, kmin):
    # the IR-divergent (frozen-mode) parts cancel pointwise in the integrand difference, so one well-conditioned integral
    f1, f2 = mode_integrand(eps, *pair1), mode_integrand(eps, *pair2)
    return mp.quad(lambda k: f1(k) - f2(k), [kmin, 1e-4, 1e-2, 0.1, 0.5, 2, 8, 30, 120, 300, 600, 1200, 2500, mp.inf])
ireg = mp.mpc(0, '-0.08')                       # the Wightman i-epsilon kept finite: both sides evaluated at the same complex t1
eps = mp.mpf('-0.25')                           # nu = 1.3 < 3/2: the mode integral converges in the IR
for (t1, t2) in [(mp.mpf('0.6'), mp.mpf('0.0')), (mp.mpf('0.3'), mp.mpf('-0.9'))]:
    wc = W_closed(eps, t1+ireg, t2); wm = W_modes(eps, t1+ireg, t2)
    check(f"eps={eps} (nu=1.3), t1={t1}, t2={t2}: closed form = mode integral (rel. {mp.nstr(abs(wc-wm)/abs(wc),2)})", abs(wc-wm)/abs(wc) < 1e-8)
eps = mp.mpf('0.2'); kIR = mp.mpf('1e-6')       # nu = 1.75 > 3/2: IR divergent; differences between two pairs are finite
ireg2 = mp.mpc(0, '-0.4')                     # larger regulator here: Im tau shrinks with t for eps > 0, so the tail would decay slowly
pairs = [(mp.mpf('0.5')+ireg2, mp.mpf('0.0')), (mp.mpf('1.2')+ireg2, mp.mpf('0.3'))]
wc = W_closed(eps, *pairs[0]) - W_closed(eps, *pairs[1]); wm = W_modes_diff(eps, pairs[0], pairs[1], kIR)
check(f"eps=0.2 (nu=1.75): closed-form difference = cutoff mode-integral difference (rel. {mp.nstr(abs(wc-wm)/abs(wc),2)}; C_IR drops out; the remainder scales as k_IR^(5-2nu))", abs(wc-wm)/abs(wc) < 1e-7)

print("== 14.C first-order theorem: W_exact - [(1 + g1 eps) W_mmc(Delta; H_X) + eps (1/8pi^2) S2(w0) + c2 Delta^2] = O(eps^2) iff H_X = H_eff ==")
mp.mp.dps = 40
def W_mmc(D, Hx): return Hx**2/(16*mp.pi**2)*2*S1(mp.cosh(Hx*D/2)**2)
def W_first(eps, D, Hx, with_S2=True):
    w0 = mp.cosh(D/2)**2
    return W_mmc(D, Hx) + eps/(16*mp.pi**2)*((2*S2(w0) if with_S2 else 0) + 2*g1*S1(w0)) - (1-eps)**3*eps/(32*mp.pi**2)*D**2
Dref = mp.mpc('0.5', '-0.3'); Ds = [mp.mpc(x, '-0.3') for x in ['0.2', '0.9', '1.6', '2.5', '3.5']]
epss = [mp.mpf('8e-3'), mp.mpf('4e-3'), mp.mpf('2e-3'), mp.mpf('1e-3')]
def resid(eps, D, Hx_of_eps, with_S2=True):
    Wf = lambda DD: W_first(eps, DD, Hx_of_eps(eps), with_S2)
    return (W_closed(eps, D/2, -D/2) - W_closed(eps, Dref/2, -Dref/2)) - (Wf(D) - Wf(Dref))
variants = [("H_eff = 1 - eps, with S2", lambda e: 1 - e, True, True), ("H_eff, S2 dropped", lambda e: 1 - e, False, False),
            ("Hbar (c = 0)", lambda e: mp.mpf(1), True, False), ("Kodama-Hayward 1 - eps/2", lambda e: 1 - e/2, True, False)]
for name, Hx, wS2, expect_second_order in variants:
    worst_ratio = 0; worst_lin = 1e9
    for D in Ds:
        R = [abs(resid(e, D, Hx, wS2)) for e in epss]
        ratios = [R[i]/R[i+1] for i in range(3)]            # ~4 for O(eps^2), ~2 for O(eps)
        worst_ratio = max(worst_ratio, max(abs(r - 4) for r in ratios)) if expect_second_order else worst_ratio
        worst_lin = min(worst_lin, min(ratios))
        print(f"    {name}: D={mp.nstr(D,3)}  |R|/eps^2 = " + ", ".join(mp.nstr(R[i]/epss[i]**2, 5) for i in range(4)) + f"   halving ratios {', '.join(mp.nstr(r,4) for r in ratios)}")
    if expect_second_order: check(f"{name}: residual is O(eps^2) (halving ratios -> 4; worst |ratio-4| = {mp.nstr(worst_ratio,3)})", worst_ratio < 0.1)
    else: check(f"{name}: residual is O(eps), NOT O(eps^2) (ratios ~2; min ratio {mp.nstr(worst_lin,4)})", worst_lin < 2.3)
mp.mp.dps = 30

print("== 14.D KMS, implementation check: the first-order shape f = 2 S2(w0) + 2 g1 S1(w0) is midline-symmetric (automatic for any function of w0) ==")
f = lambda D: 2*S2(mp.cosh(D/2)**2) + 2*g1*S1(mp.cosh(D/2)**2)
worst = 0
for s in ['0.3', '1.1', '2.7', '5.0']:
    s = mp.mpf(s); a = f(s - 1j*beta/2); b = f(-s - 1j*beta/2); worst = max(worst, abs(a-b), abs(mp.im(a)))
check(f"f(s - i beta/2) real and even in s (worst deviation {mp.nstr(worst,2)})", worst < 1e-25)
worst = 0
for D in ['0.4', '1.3', '3.0']:
    D = mp.mpf(D); e = mp.mpf('1e-8'); worst = max(worst, abs(f(D - 1j*beta + 1j*e) - f(-D - 1j*e)))
check(f"boundary values: f(Delta - i beta) = f(-Delta) (worst {mp.nstr(worst,2)})", worst < 1e-20)
print("  (these two checks verify the implementation of f, not the physics; the evidence is 14.C and 14.D'' below)")
print("== 14.D'' KMS of the EXACT closed form: K(Delta) = [W(Delta - i beta + i d) - c2 (.)^2] - [W(-Delta - i d) - c2 (.)^2], d = 0.1 ==")
mp.mp.dps = 40
def Kdef(eps, D, beta, dl=mp.mpf('0.1')):
    c2 = -(1-eps)**3*eps/(32*mp.pi**2); z1 = D - 1j*beta + 1j*dl; z2 = -D - 1j*dl
    return (W_closed(eps, z1/2, -z1/2) - c2*z1**2) - (W_closed(eps, z2/2, -z2/2) - c2*z2**2)
epsK = [mp.mpf('2e-3'), mp.mpf('1e-3'), mp.mpf('5e-4')]
for name, bfun, lo, hi in [("beta_eff = 2pi/(Hbar(1-eps))", lambda e: 2*mp.pi/(1-e), 3.5, 4.5), ("beta_KH = 2pi/(Hbar(1-eps/2))", lambda e: 2*mp.pi/(1-e/2), 1.7, 2.4), ("beta_bar = 2pi/Hbar", lambda e: 2*mp.pi, 1.7, 2.4)]:
    allok = True
    for D in [mp.mpf('0.8'), mp.mpf('2')]:
        K = [abs(Kdef(e, D, bfun(e))) for e in epsK]; ratios = [K[0]/K[1], K[1]/K[2]]
        print(f"    {name}, Delta={D}: |K| = " + ", ".join(mp.nstr(k, 3) for k in K) + f"   halving ratios {mp.nstr(ratios[0],4)}, {mp.nstr(ratios[1],4)}")
        allok &= all(lo < r < hi for r in ratios)
    check(f"{name}: exact-function KMS defect is O(eps^{2 if lo > 3 else 1})", allok)
mp.mp.dps = 30
print("  => the first-order KMS defect of the exact function (omega != 0 part) vanishes at beta_eff and not at beta_KH or beta_bar;")
print("     beta_KH - beta_eff = -pi eps/Hbar + O(eps^2), so the KH defect would be -(pi eps omega/Hbar) != 0.")

print("== 14.D' spectra via the midline, with the asymptotic polynomial removed analytically ==")
fmid = lambda s: f(s - 1j*beta/2).real; f0mid = lambda s: (2*S1(mp.cosh((s - 1j*beta/2)/2)**2)).real
A0, B0, C0 = mp.mpf(0), mp.mpf(-2), 4*L2 - 1
A1, B1, C1 = mp.mpf(-1), 4*L2 - 4, -4*L2**2 + 8*L2 + 1 - mp.pi**2/3
h0 = lambda s: f0mid(s) - (A0*s**2 + B0*abs(s) + C0); h1 = lambda s: fmid(s) - (A1*s**2 + B1*abs(s) + C1)
check(f"asymptotic polynomials correct: h0(30), h1(30) decay (|h0|={mp.nstr(abs(h0(mp.mpf(30))),2)}, |h1|={mp.nstr(abs(h1(mp.mpf(30))),2)})", abs(h0(mp.mpf(30))) < 1e-10 and abs(h1(mp.mpf(30))) < 1e-10)
def G(h, B, w):
    I = 2*mp.quad(lambda x: mp.cos(w*x)*h(x), [0, 0.5, 1, 2, 3, 5, 8, 12, 20, 30, 45, 70])
    return mp.e**(beta*w/2)*(I - 2*B/w**2)/(16*mp.pi**2)
G0x = lambda w: (w + 1/w)/(2*mp.pi*(1-mp.e**(-beta*w)))
worst = 0
for w in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('1'), mp.mpf('3'), mp.mpf('8')]:
    worst = max(worst, abs(G(h0, B0, w) - G0x(w))/G0x(w))
check(f"zeroth order: G0(omega) = (omega + Hbar^2/omega)/(2pi(1 - e^(-beta omega))) [massless dS detector spectrum] (worst rel. {mp.nstr(worst,2)})", worst < 1e-12)
print("  r(omega): rho(omega) = rho_mmc(omega; H_eff) [1 + eps r(omega)]")
rvals = {}
for w in ['0.05', '0.1', '0.25', '0.5', '1', '2', '4', '8']:
    w = mp.mpf(w); rvals[str(w)] = G(h1, B1, w)/G0x(w); print(f"    omega/Hbar = {mp.nstr(w,4)}:  r = {mp.nstr(rvals[str(w)], 8)}")
check("r(omega -> 0) -> B1/B0 = 2 - 2 ln 2 = 0.6137", abs(rvals['0.05'] - (2 - 2*L2)) < 2e-3 and abs(B1/B0 - (2 - 2*L2)) < 1e-25)
check("r(omega) -> 0 in the UV: the S2 pole coefficient -3/2 cancels g1 = +3/2 (Hadamard coefficient unchanged)", (-mp.mpf(3)/2 + g1) == 0 and abs(rvals['8.0']) < 0.05 and abs(rvals['4.0']) < 0.1)
check("r(omega) is monotone decreasing from the IR to the UV", all(rvals[a] > rvals[b] for a, b in zip(['0.05','0.1','0.25','0.5','1.0','2.0','4.0'], ['0.1','0.25','0.5','1.0','2.0','4.0','8.0'])))
r_asym = lambda w: 3*(2*w**4 - w**2 - 1)/(4*w**4*(w**2 + 1))
check(f"large-omega asymptote r -> 3(2w^4 - w^2 - 1)/(4 w^4 (w^2+1)) ~ 3/(2w^2): r(8) = {mp.nstr(rvals['8.0'],8)} vs {mp.nstr(r_asym(mp.mpf(8)),8)}", abs(rvals['8.0'] - r_asym(mp.mpf(8)))/r_asym(mp.mpf(8)) < 1e-3)

print("== 14.E the sigma-free first-order diffusion coefficient of the zero mode (two-time function at fixed mean time) ==")
Dcoef = -3 + B1/B0
check("D(tbar) = (Hbar^3/4pi^2)[1 + c_D eps], c_D = -3 + (2 - 2 ln 2) = -(1 + 2 ln 2) = -2.386", abs(Dcoef + (1 + 2*L2)) < 1e-25)
sig_star = mp.e**(1 - mp.euler)
c1 = lambda s: 1 - 2*mp.euler - 2*mp.log(2*s)
check(f"equals the stochastic c1(sigma) = 1 - 2 gamma_E - 2 ln(2 sigma) at sigma* = e^(1 - gamma_E) = {mp.nstr(sig_star,6)} (and at no other sigma)", abs(c1(sig_star) - Dcoef) < 1e-25)
print(f"  compare: naive H_eff^3 -> -3; stochastic at sigma = 1 -> {mp.nstr(c1(mp.mpf(1)),5)}; exact -> {mp.nstr(Dcoef,6)}")

print("\nCONCLUSION: at first order in eps the fixed-mean-time (Wigner) two-point spectrum of the massless minimally coupled sector along a")
print("comoving worldline is, at omega != 0, KMS at T = (a''/a')/2pi + O(eps^2) -- the same statement, in the same convention, as v11's for the")
print("conformal sector; the Kodama-Hayward rate is excluded at O(eps). (An instantaneous one-sided detector rate is a different convention;")
print("see the docstring.)")
print("the result is sigma-free. The O(eps) effect is a shape change of the spectral density, rho = rho_mmc(omega; H_eff)[1 + eps r(omega)],")
print("r(0) = 2 - 2 ln 2, r(inf) = 0, and the zero-mode diffusion coefficient is D = (Hbar^3/4pi^2)[1 - (1 + 2 ln 2) eps].")
print("This is conventional QFT; it establishes nothing about P1 beyond consistency. O5 item 1 is now computed, not argued.")
print("\nV14 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
