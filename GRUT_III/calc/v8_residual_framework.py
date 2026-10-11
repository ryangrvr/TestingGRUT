#!/usr/bin/env python3
"""
V8 -- The residual framework: what any 'GRUT memory' claim must actually exhibit.

THESIS.  For a fixed observable O and a fixed state, let G+(w) be the Wightman spectrum
    G+(w) = Int dt e^{i w t} <O(t) O(0)>        (w real, both signs)
and define
    rho(w) = G+(w) - G+(-w)                     commutator (dissipative / response) spectrum, odd,
    N(w)   = (G+(w) + G+(-w)) / 2               symmetrised noise spectrum, even.
If the state is KMS at inverse temperature beta  (G+(w)/G+(-w) = e^{beta w})  then IDENTICALLY
    N(w) = (1/2) coth(beta w / 2) rho(w).                                             (FDT)
So in any KMS state the noise carries NO information beyond rho: the only object that can
distinguish a 'GRUT' medium from the conventional one, for a given (observable, state), is the residual
    Delta rho(w) = rho_GRUT(w) - rho_conventional(w),
and a claim of 'extra memory' with Delta rho = 0 is empty, while a claim with Delta rho != 0 must name
the operator content and state that produce it (otherwise it is SUPPLIED, in SCOUT_0's sense).

CONVENTION NOTE.  v1_ds_thermal_kernel.py prints the 'symmetrised spectrum S(w) = (w/2pi) coth(beta w/2)',
computed there as G+(w)+G+(-w); in the present convention that is 2N(w).  The record's
calc/finite_T_exponent.py writes S(w,T) = J(w) coth(w/2T) (hbar = k_B = omega_c = 1); with T = H/2pi
that is coth(beta w/2) J(w) = 2N(w) for rho = J.  Amplitudes are free in both records; only shapes
are compared here.

CHECKS (PASS/FAIL printed per check; exit status 0 iff all pass):
 (A) helpers.  rho_from_Gplus, N_from_Gplus, kms_defect = |N - (1/2)coth(beta w/2) rho|, residual
     Delta rho between two models, loglog_slope and shape_residual (difference of local log-log slopes).
     Tested on the exact Planckian G+(w) = (w/2pi)/(1 - e^{-beta w}) (rho = w/2pi, N = (w/4pi)coth(beta w/2),
     kms_defect = 0 to 1e-25), on a planted NON-KMS G+ (defect must be detected: mutation test), on the
     zero residual of a model against itself, and on power laws w^a vs w^b (shape residual must equal a-b).
 (B) THEOREM CHECK (FDT cancellation for the pointwise conformal <T00 T00> correlator on the dS worldline).
     With the v1 closed form
        G+_T(w) = w (w^2+H^2)(w^2+4H^2)(w^2+9H^2) / (10080 pi^3 (1 - exp(-2 pi w/H))),   beta = 2 pi/H,
     sympy shows EXACTLY
        rho_T(w) = G+_T(w) - G+_T(-w) = P7(w)/(10080 pi^3),  P7 = w^7 + 14 H^2 w^5 + 49 H^4 w^3 + 36 H^6 w,
     i.e. the thermal denominator cancels (1/(1-e^{-x}) + 1/(1-e^{x}) = 1): rho_T is an odd POLYNOMIAL, so the
     dissipative part of the retarded kernel is a finite sum of delta-function derivatives (local in time),
     and N_T(w) = (1/2) coth(pi w/H) rho_T(w) exactly (symbolically and numerically), with N_T(0) = H^7/(560 pi^4).
     Also recorded: the local log-log slope of rho_T is 7 for w >> H (flat limit) and 1 for w << H (the 36 H^6 w
     term): the dS conformal medium's commutator spectrum is Ohmic-like in the deep IR and w^7 in the UV; it is
     NEVER w^3 (super-Ohmic s=3) and NEVER flat (s=0).
 (C) The OLD GRUT kernels as examples of Delta rho != 0 by fiat (beta = 2 pi/H, H = 1):
     (C1) calc/finite_T_exponent.py rung-3 noise model: rho_old1(w) = J(w) = w^3 exp(-(w/omega_c)^2),
          omega_c = 1 (record's unit; the cutoff location relative to H is itself a free choice).
     (C2) SCOUT_0 S5-1 Lorentzian-exponent family (P06/P06B): S_p(x) ∝ (x^2+kappa^2)^{-p}, kappa = 1/(2 sqrt 2.3),
          p = 1 the recorded kernel, read as an odd commutator spectrum rho_old2(w) = sgn(w) (w^2+kappa^2)^{-p}
          (the kernel f_p(t) ∝ t^{p-1/2} K_{p-1/2}(kappa t) is its Fourier transform).  A Drude-friction reading,
          rho = w (w^2+kappa^2)^{-p}, is tabulated as well.
     For each: normalise to agree with rho_T at w = 1, tabulate rho_old / rho_T and the local slopes at
     w = 0.1, 0.3, 1, 3, 10, report the IR exponent (3; 0 (all p) resp. 1 for the Drude reading) and the UV exponent
     (Gaussian cutoff -> no power law; -2p resp. 1-2p) against the conformal medium's (1 in the IR, 7 in the UV),
     and the shape residual.  PASS criteria: |Delta rho|/rho_T >= 10% at every tabulated w != 1 and
     |shape residual| >= 0.5 at w = 0.1 and w = 10 for C1 and C2; for the Drude reading, whose IR exponent 1
     coincides with the medium's, at w = 1 and w = 10 (its shape residual at w = 0.1 is only -0.20, stated as such).
 (D) residual_report(model_rho, observable_label, state_label, ...) -- the template a future GRUT memory claim
     must fill in (operator, state, reference medium, mechanism); demonstrated on (C1), (C2) and on the
     conformal medium itself (which must return Delta rho = 0).

Everything numerical is mpmath at 30 digits; numpy is used only for the grid bookkeeping in (D).
"""
import sys
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
H = mp.mpf(1)
BETA = 2 * mp.pi / H                    # de Sitter: T = H / 2 pi
OMEGAS = [mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf(1), mp.mpf(3), mp.mpf(10)]

ok_all = True
def check(name, cond, detail=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return bool(cond)
def close(a, b, tol=mp.mpf('1e-20'), rel=True):
    a, b = mp.mpf(a), mp.mpf(b)
    d = abs(a - b)
    return d <= tol * (abs(b) if (rel and b != 0) else 1)

# ---------------------------------------------------------------------------------------------
# (A) helpers
# ---------------------------------------------------------------------------------------------
def rho_from_Gplus(Gp):
    """Commutator spectrum rho(w) = G+(w) - G+(-w) (odd)."""
    return lambda w: Gp(mp.mpf(w)) - Gp(-mp.mpf(w))

def N_from_Gplus(Gp):
    """Symmetrised noise N(w) = (G+(w) + G+(-w))/2 (even)."""
    return lambda w: (Gp(mp.mpf(w)) + Gp(-mp.mpf(w))) / 2

def coth(x):
    x = mp.mpf(x)
    return mp.cosh(x) / mp.sinh(x)

def kms_defect(N, rho, beta, w):
    """|N(w) - (1/2) coth(beta w/2) rho(w)|: zero in any KMS state at inverse temperature beta."""
    w = mp.mpf(w)
    return abs(N(w) - coth(beta * w / 2) * rho(w) / 2)

def N_from_rho_kms(rho, beta):
    """The noise a KMS state at inverse temperature beta MUST carry for a given rho."""
    return lambda w: coth(beta * mp.mpf(w) / 2) * rho(w) / 2

def residual(rho_model, rho_ref):
    """Delta rho(w) = rho_model(w) - rho_ref(w)."""
    return lambda w: rho_model(w) - rho_ref(w)

def loglog_slope(f, w, h=mp.mpf('1e-8')):
    """Local exponent d ln|f| / d ln w at w (central difference in log w)."""
    w = mp.mpf(w)
    return (mp.log(abs(f(w * mp.e**h))) - mp.log(abs(f(w * mp.e**(-h))))) / (2 * h)

def shape_residual(rho_model, rho_ref, w):
    """Difference of local log-log slopes: zero iff the two spectra have the same local power law."""
    return loglog_slope(rho_model, w) - loglog_slope(rho_ref, w)

def normalise_at(rho, rho_ref, w0=mp.mpf(1)):
    """Rescale rho by a positive constant so that rho(w0) = rho_ref(w0) (amplitudes are free in the record)."""
    c = rho_ref(w0) / rho(w0)
    return (lambda w: c * rho(w)), c

print("=" * 96)
print("V8  residual framework:  Delta rho = rho_GRUT - rho_conventional is the decisive object   (H = 1, beta = 2 pi/H)")
print("=" * 96)
print("\n== (A) helpers, tested on exact objects ==")
planck = lambda w: (mp.mpf(w) / (2 * mp.pi)) / (1 - mp.e**(-BETA * mp.mpf(w)))
rho_pl, N_pl = rho_from_Gplus(planck), N_from_Gplus(planck)
check("Planckian G+: rho(w) = w/2pi exactly at w in {0.1, 1, 3}",
      all(close(rho_pl(w), w / (2 * mp.pi)) for w in (mp.mpf('0.1'), mp.mpf(1), mp.mpf(3))),
      f"rho(1) = {mp.nstr(rho_pl(1), 15)} vs 1/2pi = {mp.nstr(1/(2*mp.pi), 15)}")
check("Planckian G+: N(w) = (w/4pi) coth(beta w/2)",
      all(close(N_pl(w), w / (4 * mp.pi) * coth(BETA * w / 2)) for w in OMEGAS))
dmax = max(kms_defect(N_pl, rho_pl, BETA, w) for w in OMEGAS)
check("KMS defect of the Planckian spectrum = 0", dmax < mp.mpf('1e-25'), f"max defect {mp.nstr(dmax, 3)}")
# mutation: a non-KMS G+ (wrong temperature) must be flagged
bad = lambda w: (mp.mpf(w) / (2 * mp.pi)) / (1 - mp.e**(-2 * BETA * mp.mpf(w)))
dbad = kms_defect(N_from_Gplus(bad), rho_from_Gplus(bad), BETA, mp.mpf('0.5'))
check("mutation: G+ at temperature T/2 is flagged (defect > 1e-3)", dbad > mp.mpf('1e-3'),
      f"defect at w=0.5: {mp.nstr(dbad, 6)}  -- tester live")
check("residual of a model against itself vanishes",
      all(residual(rho_pl, rho_pl)(w) == 0 for w in OMEGAS))
pa, pb = (lambda w: mp.mpf(w)**3), (lambda w: mp.mpf(w)**7)
sr = shape_residual(pa, pb, mp.mpf('0.7'))
check("shape residual of w^3 against w^7 = -4", close(sr, -4, tol=mp.mpf('1e-10')), f"got {mp.nstr(sr, 12)}")

# ---------------------------------------------------------------------------------------------
# (B) theorem check: FDT cancellation for the conformal <T00 T00> correlator on the dS worldline
# ---------------------------------------------------------------------------------------------
print("\n== (B) FDT cancellation for G+_T (v1 closed form): rho_T is an odd polynomial ==")
w, Hs = sp.symbols('omega H', positive=True)
GpT_s = w * (w**2 + Hs**2) * (w**2 + 4*Hs**2) * (w**2 + 9*Hs**2) / (10080 * sp.pi**3 * (1 - sp.exp(-2*sp.pi*w/Hs)))
P7 = sp.expand(w * (w**2 + Hs**2) * (w**2 + 4*Hs**2) * (w**2 + 9*Hs**2))
rhoT_s = sp.simplify((GpT_s - GpT_s.subs(w, -w)).rewrite(sp.exp))
diff_rho = sp.simplify(sp.expand(rhoT_s - P7 / (10080 * sp.pi**3)))
check("sympy: G+_T(w) - G+_T(-w) == P7(w)/(10080 pi^3)", diff_rho == 0,
      f"P7 = {P7}")
coeffs = sp.Poly(P7, w).all_coeffs()
check("P7 coefficients (1, 0, 14H^2, 0, 49H^4, 0, 36H^6, 0): odd polynomial of degree 7",
      coeffs == [1, 0, 14*Hs**2, 0, 49*Hs**4, 0, 36*Hs**6, 0], f"coeffs = {coeffs}")
# the identity behind it
x = sp.symbols('x', real=True)
ident = sp.simplify(1/(1 - sp.exp(-x)) + 1/(1 - sp.exp(x)) - 1)
check("thermal denominators cancel: 1/(1-e^-x) + 1/(1-e^x) = 1", ident == 0)
NT_s = (GpT_s + GpT_s.subs(w, -w)) / 2
fdt_s = sp.simplify((NT_s - sp.coth(sp.pi*w/Hs) * P7 / (10080*sp.pi**3) / 2).rewrite(sp.exp))
check("sympy: N_T(w) == (1/2) coth(pi w/H) rho_T(w)", fdt_s == 0)
floor_s = sp.limit(NT_s, w, 0)
check("N_T(0) = G+_T(0) = H^7/(560 pi^4)  (classical Nyquist floor)", sp.simplify(floor_s - Hs**7/(560*sp.pi**4)) == 0,
      f"limit = {floor_s}")

def GpT(w):
    w = mp.mpf(w)
    return w*(w**2+H**2)*(w**2+4*H**2)*(w**2+9*H**2) / (10080*mp.pi**3*(1 - mp.e**(-2*mp.pi*w/H)))
def rho_T(w):
    w = mp.mpf(w)
    return (w**7 + 14*H**2*w**5 + 49*H**4*w**3 + 36*H**6*w) / (10080*mp.pi**3)
N_T = N_from_Gplus(GpT)
check("numeric: rho_from_Gplus(G+_T) == P7/(10080 pi^3) at the table frequencies",
      all(close(rho_from_Gplus(GpT)(w), rho_T(w)) for w in OMEGAS))
dT = max(kms_defect(N_T, rho_T, BETA, w) for w in OMEGAS)
check("numeric: KMS defect of G+_T = 0", dT < mp.mpf('1e-25'), f"max defect {mp.nstr(dT, 3)}")
sIR, sUV = loglog_slope(rho_T, mp.mpf('1e-3')), loglog_slope(rho_T, mp.mpf('1e3'))
check("rho_T local exponent: 1 at w = 1e-3 H (IR, from 36 H^6 w) and 7 at w = 1e3 H (flat limit)",
      close(sIR, 1, tol=mp.mpf('1e-4')) and close(sUV, 7, tol=mp.mpf('1e-4')),
      f"slope(1e-3) = {mp.nstr(sIR, 10)}, slope(1e3) = {mp.nstr(sUV, 10)}")
print("  -> rho_T(t) = sum_k c_k delta^(k)(t): the dissipative kernel at a point is LOCAL (Huygens); all memory sits in")
print("     N_T = (1/2) coth(pi w/H) rho_T, which is fixed by rho_T and the temperature T = H/2pi -- nothing is left to supply.")

# ---------------------------------------------------------------------------------------------
# (C) the old GRUT kernels: Delta rho != 0 by fiat
# ---------------------------------------------------------------------------------------------
print("\n== (C) old GRUT kernels against the conformal medium  (normalised to agree with rho_T at w = 1) ==")
OMEGA_C = mp.mpf(1)                                  # calc/finite_T_exponent.py: WC = 1.0, ETA = 1.0
def J_old(w):                                        # J(w) = ETA w^3 exp(-(w/WC)^2); S(w,T) = J coth(w/2T)
    w = mp.mpf(w)
    return w**3 * mp.e**(-(w/OMEGA_C)**2)
KAPPA = 1 / (2 * mp.sqrt(mp.mpf('2.3')))             # S5_OWNER_RULING_03.md:49 ; P06: S_p(x) ∝ (x^2+kappa^2)^-p
def lorentz_odd(p):
    return lambda w: mp.sign(mp.mpf(w)) * (mp.mpf(w)**2 + KAPPA**2)**(-mp.mpf(p))
def lorentz_drude(p):
    return lambda w: mp.mpf(w) * (mp.mpf(w)**2 + KAPPA**2)**(-mp.mpf(p))

def table(label, rho_old, ir_expected, ir_tol=mp.mpf('2e-2'), shape_ws=(OMEGAS[0], OMEGAS[-1])):
    """Normalise, tabulate, and test Delta rho != 0 and the shape residual (|slope difference| >= 0.5 at shape_ws)."""
    rho_n, c = normalise_at(rho_old, rho_T)
    print(f"\n  --- {label} ---   normalisation constant c = {mp.nstr(c, 8)} (the record left this amplitude free)")
    print(f"  {'w':>6} | {'rho_old/rho_T':>14} | {'Delta rho / rho_T':>17} | {'slope_old':>10} | {'slope_T':>8} | {'shape resid':>11}")
    rel_res, shapes = {}, {}
    for w_ in OMEGAS:
        r = rho_n(w_) / rho_T(w_)
        so, sT = loglog_slope(rho_n, w_), loglog_slope(rho_T, w_)
        rel_res[w_], shapes[w_] = r - 1, so - sT
        print(f"  {mp.nstr(w_, 3):>6} | {mp.nstr(r, 8):>14} | {mp.nstr(r - 1, 8):>17} | {mp.nstr(so, 6):>10} | {mp.nstr(sT, 6):>8} | {mp.nstr(so - sT, 6):>11}")
    ir_old = loglog_slope(rho_n, mp.mpf('1e-3'))
    ir_T = loglog_slope(rho_T, mp.mpf('1e-3'))
    print(f"  IR exponent (w = 1e-3): old = {mp.nstr(ir_old, 8)}  vs conformal medium = {mp.nstr(ir_T, 8)};  "
          f"UV exponent of the medium (w = 1e3) = {mp.nstr(loglog_slope(rho_T, mp.mpf('1e3')), 8)}")
    check(f"{label}: IR exponent = {ir_expected}", close(ir_old, ir_expected, tol=ir_tol, rel=False),
          f"got {mp.nstr(ir_old, 8)}")
    check(f"{label}: Delta rho != 0 -- |Delta rho|/rho_T >= 10% at every tabulated w != 1",
          all(abs(rel_res[w_]) >= mp.mpf('0.1') for w_ in OMEGAS if w_ != 1),
          "min |Delta rho|/rho_T off w=1: " + mp.nstr(min(abs(rel_res[w_]) for w_ in OMEGAS if w_ != 1), 6))
    check(f"{label}: shape residual |slope_old - slope_T| >= 0.5 at w = " + " and ".join(mp.nstr(w_, 3) for w_ in shape_ws),
          all(abs(shapes[w_]) >= mp.mpf('0.5') for w_ in shape_ws),
          "shape residual: " + ", ".join(f"at {mp.nstr(w_, 3)}: {mp.nstr(shapes[w_], 6)}" for w_ in shape_ws))
    return rho_n, rel_res, shapes

print("  Conformal medium reference: rho_T(w) = (w^7 + 14 w^5 + 49 w^3 + 36 w)/(10080 pi^3), slope 1 (w<<H) -> 7 (w>>H).")
rho1, res1, sh1 = table("C1  rung-3 noise model  rho = J = w^3 exp(-(w/omega_c)^2), omega_c = 1  [calc/finite_T_exponent.py]",
                        J_old, 3)
print("  (the record's S(w,T) = J coth(w/2T) is 2 N(w) for rho = J: KMS-locked by construction, so N adds no information)")
uv1 = loglog_slope(J_old, mp.mpf(3))
print(f"  UV: the Gaussian cutoff gives slope 3 - 2 (w/omega_c)^2 = {mp.nstr(uv1, 6)} already at w = 3: no power law, "
      f"vs 7 for the medium; the cutoff scale omega_c/H is a second free choice")
rho2, res2, sh2 = table("C2  S5-1 Lorentzian, odd extension  rho = sgn(w) (w^2+kappa^2)^-1, kappa = 1/(2 sqrt 2.3)  [P06/P06B, p = 1]",
                        lorentz_odd(1), 0)
print(f"  UV exponent of (w^2+kappa^2)^-p is -2p = -2 (p = 1): the family value;  kappa = {mp.nstr(KAPPA, 8)}")
for p in (mp.mpf('0.5'), mp.mpf(2)):
    f = lorentz_odd(p)
    print(f"  family member p = {mp.nstr(p, 3)}: IR exponent {mp.nstr(loglog_slope(f, mp.mpf('1e-3')), 6)} (0 for every p), "
          f"UV exponent {mp.nstr(loglog_slope(f, mp.mpf('1e3')), 6)} (= -2p)")
rho2d, res2d, sh2d = table("C2' S5-1 Lorentzian, Drude-friction reading  rho = w (w^2+kappa^2)^-1  (Ohmic-Drude)",
                           lorentz_drude(1), 1, shape_ws=(OMEGAS[2], OMEGAS[-1]))
print(f"  NOTE: in the Drude reading the IR exponent 1 coincides with the medium's IR exponent 1 (w << H), so the shape residual")
print(f"        at w = 0.1 is only {mp.nstr(sh2d[OMEGAS[0]], 4)} (tested at w = 1 and 10 instead); the amplitude residual there is still")
print(f"        {mp.nstr(res2d[OMEGAS[0]], 4)} x rho_T, and for w >= 1 Drude falls as 1/w while the medium rises as w^7: Delta rho != 0 regardless.")
print("\n  CONCLUSION (C): every registered GRUT kernel asserts Delta rho != 0 against the conformal vacuum medium of the")
print("  same (observable, state), with exponents 3 / 0 (or 1) in place of the medium's 1 -> 7, and gives no operator or state")
print("  that would produce the difference.  That is exactly what SCOUT_0 recorded as 'supplied' (noise profile, kernel class).")

# ---------------------------------------------------------------------------------------------
# (D) the template a future claim must fill in
# ---------------------------------------------------------------------------------------------
print("\n== (D) residual_report template, demonstrated ==")
def residual_report(model_rho, observable_label, state_label, *, reference_rho=rho_T, beta=BETA,
                    omegas=OMEGAS, mechanism=None, normalise_w0=mp.mpf(1), tol=mp.mpf('1e-6')):
    """TEMPLATE for any future 'GRUT memory' claim.

    A claim must supply:
      model_rho        : callable w -> rho_GRUT(w), the commutator spectrum of the SAME observable in the SAME state
                         as the reference (odd in w); the noise then follows from KMS and is NOT an independent input.
      observable_label : which operator (e.g. 'T^00 at a point on a comoving worldline').
      state_label      : which state (e.g. 'Bunch-Davies, geodesic observer, T = H/2pi').
      reference_rho    : rho_conventional for that observable/state (default: the dS conformal medium rho_T of v1).
      mechanism        : the operator/state content that produces Delta rho != 0, or None.
    Returns a dict with the residual table, local exponents, and a verdict:
      'EMPTY'    if Delta rho = 0 (nothing beyond the conventional medium is being claimed),
      'SUPPLIED' if Delta rho != 0 and mechanism is None (the record's situation),
      'TESTABLE' if Delta rho != 0 and a mechanism is named (then Delta rho is the prediction to falsify).
    """
    rho_n, c = normalise_at(model_rho, reference_rho, normalise_w0)
    N_model = N_from_rho_kms(rho_n, beta)                   # the only noise a KMS state allows for this rho
    om = np.array([float(w_) for w_ in omegas])
    rel = np.array([float(rho_n(w_) / reference_rho(w_) - 1) for w_ in omegas])
    shape = np.array([float(shape_residual(rho_n, reference_rho, w_)) for w_ in omegas])
    kms = float(max(kms_defect(N_model, rho_n, beta, w_) for w_ in omegas))
    nonzero = bool(np.max(np.abs(rel)) > float(tol))
    verdict = "EMPTY (Delta rho = 0)" if not nonzero else (
        "SUPPLIED (Delta rho != 0, no mechanism named)" if mechanism is None else
        f"TESTABLE (Delta rho != 0; mechanism: {mechanism})")
    rep = dict(observable=observable_label, state=state_label, beta=float(beta), normalisation=float(c),
               omegas=om, rel_residual=rel, shape_residual=shape, kms_defect=kms,
               ir_exponent=float(loglog_slope(rho_n, mp.mpf('1e-3'))), ir_exponent_ref=float(loglog_slope(reference_rho, mp.mpf('1e-3'))),
               uv_exponent_ref=float(loglog_slope(reference_rho, mp.mpf('1e3'))), mechanism=mechanism, verdict=verdict)
    print(f"  residual_report: observable = {observable_label}; state = {state_label}")
    print(f"    max |Delta rho|/rho_ref = {np.max(np.abs(rel)):.4g}; max |shape residual| = {np.max(np.abs(shape)):.4g}; "
          f"IR exponent model/ref = {rep['ir_exponent']:.4f}/{rep['ir_exponent_ref']:.4f}; KMS defect of implied N = {kms:.1e}")
    print(f"    VERDICT: {verdict}")
    return rep

r0 = residual_report(rho_T, "T^00 at a point on a comoving worldline (one conformal scalar)",
                     "Bunch-Davies / geodesic observer, KMS at T = H/2pi")
check("template on the conformal medium itself returns EMPTY (Delta rho = 0)", r0['verdict'].startswith("EMPTY"))
r1 = residual_report(J_old, "unspecified ('fast modes', rung 3)", "unspecified thermal state, T free")
check("template on C1 returns SUPPLIED", r1['verdict'].startswith("SUPPLIED"),
      f"IR exponent {r1['ir_exponent']:.4f} vs ref {r1['ir_exponent_ref']:.4f}")
r2 = residual_report(lorentz_odd(1), "unspecified (S5-1 scaled-limit resolvent density)", "unspecified, kappa = 1/(2 sqrt 2.3)")
check("template on C2 returns SUPPLIED", r2['verdict'].startswith("SUPPLIED"),
      f"IR exponent {r2['ir_exponent']:.4f} vs ref {r2['ir_exponent_ref']:.4f}")
r3 = residual_report(J_old, "unspecified", "unspecified", mechanism="(hypothetical) named operator content")
check("template on C1 with a named mechanism returns TESTABLE", r3['verdict'].startswith("TESTABLE"))
check("template: implied noise of every model is KMS-consistent by construction (defect < 1e-25)",
      max(r['kms_defect'] for r in (r0, r1, r2, r3)) < 1e-25)

print("\nV8 OVERALL:", "PASS" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
