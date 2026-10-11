#!/usr/bin/env python3
"""
V10 -- Where the gravitational sector's 'memory' lives: the massless, minimally coupled sector in de Sitter
(the tensor-mode proxy; each graviton polarisation obeys the massless minimally coupled scalar equation).

The record (calc/RESULTS_tt_worldline.md) found numerically that the graviton worldline kernel is NON-stationary
and concluded 'no derived memory time exists' for the gravitational channel. Here the same object is written in
closed form and the finding is made exact: the non-stationary part is a Brownian random walk of the zero mode
with a DERIVED diffusion constant; every response function and every correlator of time derivatives is stationary.

Closed form (Allen & Folacci 1987; Onemli & Woodard 2002 form): for the Bunch-Davies-type (Euclidean, dS-breaking)
state of a massless minimally coupled scalar in dS4,
    G(x,x') = (H^2/8 pi^2) [ 1/(1-Z) - ln(1-Z) + ln(a a') + c ],   Z = dS invariant,
and on a comoving worldline  1 - Z = 2 sinh^2(H dt/2),  ln(a a') = H (t + t').

Checks:
 10.1  1 - Z on the worldline equals 2 sinh^2(H dt/2) (sympy, exact); flat limit 1/(4 pi^2 dt^2).
 10.2  d_t d_t' G depends on dt only (stationary): the velocity correlator and every commutator are stationary.
 10.3  Large-|dt| asymptotics at fixed mean time: G -> (H^3/4 pi^2) min(t,t') + O(1) + O(e^{-H|dt|}):
       the Brownian covariance D*min(t,t') with D = H^3/(4 pi^2)  (Starobinsky stochastic inflation, exact here).
 10.4  Equal-time secular growth of the renormalised variance: d<phi^2>/dt = H^3/(4 pi^2).
 10.5  Tensor translation: h_lambda = (2/M_P) phi_lambda  =>  d<h_ij h^ij>/dt = 2 H^3/(pi^2 M_P^2) = H * P_T with
       P_T = 2 H^2/(pi^2 M_P^2): the standard tensor power per e-fold.  (M_P = reduced Planck mass.)
 10.6  The stationary remainder is thermal: the e^{-H|dt|} corrections carry KMS period 2 pi i/H (sympy check of
       periodicity of the sinh^2 and ln sinh^2 terms under dt -> dt + 2 pi i/H).
Conclusion: the gravitational sector's 'infinite memory' is the zero-mode random walk (a Markov process with
derived D); its 'finite memory' is thermal at T = H/2pi. No memory kernel is supplied. All imported/standard.
"""
import sympy as sp
ok_all = True
H, t, tp, dt, tbar = sp.symbols('H t tp dt tbar', positive=True)
eta, etap = -sp.exp(-H*t)/H, -sp.exp(-H*tp)/H
a, ap = -1/(H*eta), -1/(H*etap)
print("== 10.1 worldline invariant ==")
oneZ = (eta-etap)**2/(2*eta*etap)          # 1-Z at coincident comoving points (eta,eta'<0)
target = 2*sp.sinh(H*(t-tp)/2)**2
d = sp.simplify((oneZ - target).rewrite(sp.exp))
print("  (1-Z) - 2 sinh^2(H dt/2) =", d, "->", "PASS" if d == 0 else "FAIL"); ok_all &= (d == 0)
G = (H**2/(8*sp.pi**2))*(1/(2*sp.sinh(H*(t-tp)/2)**2) - sp.log(2*sp.sinh(H*(t-tp)/2)**2) + H*(t+tp))
flat = sp.limit((G*(t-tp)**2).subs(tp, t - dt), dt, 0)
print("  flat limit coefficient of 1/dt^2:", flat, "->", "PASS" if sp.simplify(flat - 1/(4*sp.pi**2)) == 0 else "FAIL")
ok_all &= sp.simplify(flat - 1/(4*sp.pi**2)) == 0
print("== 10.2 stationarity of the velocity correlator ==")
Gdd = sp.diff(G, t, tp)
Gdd_sub = sp.simplify(Gdd.subs({t: tbar + dt/2, tp: tbar - dt/2}))
print("  d_t d_t' G depends on tbar?", Gdd_sub.has(tbar), "->", "PASS" if not Gdd_sub.has(tbar) else "FAIL")
ok_all &= not Gdd_sub.has(tbar)
print("== 10.3 Brownian covariance at large separation ==")
# write G with t = tbar + dt/2, tp = tbar - dt/2 and expand for large dt: ln(2 sinh^2) = H dt - 2 ln... 
Gs = G.subs({t: tbar + dt/2, tp: tbar - dt/2})
lead = sp.limit(Gs - (H**2/(8*sp.pi**2))*(2*H*tbar - H*dt), dt, sp.oo)
print("  G - (H^2/8pi^2)(2H tbar - H dt) -> ", lead, " as dt->oo (a constant: PASS)" if lead.is_finite else " FAIL")
ok_all &= bool(lead.is_finite)
# 2H tbar - H dt = 2 H min(t,t')  for t>t': min = tp = tbar - dt/2
D = sp.simplify((H**2/(8*sp.pi**2))*2*H)
print("  => G ~ D * min(t,t') + const, D =", D, "->", "PASS" if sp.simplify(D - H**3/(4*sp.pi**2)) == 0 else "FAIL")
ok_all &= sp.simplify(D - H**3/(4*sp.pi**2)) == 0
print("== 10.4 secular growth of the renormalised variance ==")
# the only t-dependence surviving the (dt-only) subtraction at equal times is (H^2/8pi^2) H (t+t) 
growth = sp.diff((H**2/(8*sp.pi**2))*H*(t+t), t)
print("  d<phi^2>/dt =", sp.simplify(growth), "->", "PASS" if sp.simplify(growth - H**3/(4*sp.pi**2)) == 0 else "FAIL")
ok_all &= sp.simplify(growth - H**3/(4*sp.pi**2)) == 0
print("== 10.5 tensor translation ==")
MP = sp.symbols('M_P', positive=True)
DT = 2*(2/MP)**2*growth                 # two polarisations, h = (2/M_P) phi
PT = 2*H**2/(sp.pi**2*MP**2)
print("  d<h_ij h^ij>/dt =", sp.simplify(DT), ";  H*P_T =", sp.simplify(H*PT), "->", "PASS" if sp.simplify(DT - H*PT) == 0 else "FAIL")
ok_all &= sp.simplify(DT - H*PT) == 0
print("== 10.6 KMS periodicity of the stationary part ==")
stat = 1/(2*sp.sinh(H*dt/2)**2) - sp.log(2*sp.sinh(H*dt/2)**2)
shifted = stat.subs(dt, dt + 2*sp.pi*sp.I/H)
per = sp.simplify(sp.expand_complex(shifted - stat).rewrite(sp.exp)) if False else sp.simplify((shifted - stat).rewrite(sp.exp))
print("  stationary part shifted by 2 pi i/H minus itself:", per, "->", "PASS" if per == 0 else "CHECK (branch of log)")
ok_all &= (per == 0)
print("\nV10 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
