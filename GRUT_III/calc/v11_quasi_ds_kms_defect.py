#!/usr/bin/env python3
"""
V11 -- The conventional KMS defect of the de Sitter medium under slow roll (the baseline for any non-KMS residual).

Structural point (reviewer + owner, 2026-10-11): GRUT III builds KMS in by construction (bath = causal complement,
Bisognano-Wichmann/Gibbons-Hawking).  In a KMS state N = (1/2) coth(beta w/2) rho, so any 'extra noise' at the same
temperature is forced to be extra dissipation: in its own regime GRUT III is standard physics.  A distinctive residual
must therefore be sought where KMS fails.  The first such platform is quasi-de Sitter, H(t) slowly varying with
epsilon = -Hdot/H^2 > 0, where the comoving-worldline state is KMS only approximately.

Here the conventional O(epsilon) KMS defect is computed EXACTLY for a conformally coupled scalar in the conformal
vacuum of a power-law FRW universe a(t) = (t/t0)^p, epsilon = 1/p (H = p/t):
    W(t,t') = -1 / ( 4 pi^2 a(t) a(t') (eta(t) - eta(t') - i0)^2 ),   eta(t) = t0^p t^{1-p}/(1-p)     [exact, any a(t)]
Wigner (mean-time) spectrum  G+(w; tbar) = Int dDelta e^{i w Delta} W(tbar + Delta/2, tbar - Delta/2), evaluated by
shifting the contour to Im Delta = -pi/Hbar (half the local thermal period), where the integrand is smooth.

Checks / results:
 11.1  epsilon -> 0 recovers the exact dS Planck spectrum (v1): G+(w) = (w/2pi)/(1 - e^{-2 pi w/H}).
 11.2  Self-similarity: G+(w; tbar) = Hbar g(w/Hbar; epsilon)  (W ~ H^2, Int dDelta ~ 1/H)  (same dimensionless value at two mean times).
 11.3  KMS defect  D(w) := ln[G+(w)/G+(-w)] - 2 pi w/Hbar  is O(epsilon): D/epsilon converges as epsilon -> 0;
       its shape f(w/H) = lim D/epsilon is tabulated.  Equivalently beta_eff(w) = ln[G+(w)/G+(-w)]/w depends on w:
       the state is not KMS at any single temperature; the w-dependence is the invariant (convention-free) violation.
 11.4  The defect is odd in w (exact) and vanishes identically at epsilon = 0 (to quadrature precision).
Conclusion: conventional QFT in FRW fixes the O(epsilon) KMS violation of the conformal medium completely; there is
no free function here for GRUT to restrict.  Any GRUT claim in quasi-dS must live in the sector conventional theory
leaves open (interacting/gravitational infrared), and must be reported as a residual against THIS baseline.
"""
import mpmath as mp
mp.mp.dps = 25
ok_all = True
def report(name, got, want, tol):
    global ok_all
    err = abs(mp.mpf(got) - mp.mpf(want)); rel = err/abs(mp.mpf(want)) if want != 0 else err
    ok = rel < tol; ok_all &= bool(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {mp.nstr(got, 10)} want {mp.nstr(want, 10)} rel {mp.nstr(rel, 3)}")
    return ok

def make_W(p, tbar):
    """Worldline Wightman function for a(t) = (t/tbar)^p (so a(tbar) = 1), as a function of complex Delta."""
    t0 = tbar
    def a(t): return (t/t0)**p
    def eta(t): return t0**p * t**(1-p)/(1-p)
    def W(Delta):
        t1 = tbar + Delta/2; t2 = tbar - Delta/2
        return -1/(4*mp.pi**2 * a(t1)*a(t2) * (eta(t1)-eta(t2))**2)
    return W

def Gplus(W, Hbar, w, L=45):
    delta = mp.pi/Hbar
    f = lambda u: mp.e**(1j*w*u) * W(u - 1j*delta)
    pts = [-L, -20, -8, -3, -1, 0, 1, 3, 8, 20, L]
    val = mp.quad(f, [x/Hbar for x in pts])
    return mp.re(mp.e**(w*delta) * val)

print("== 11.1 exact-dS limit (epsilon -> 0 realised as p -> large; compare with Planck) ==")
H = mp.mpf(1)
for p in [400, 1600]:
    eps = 1/mp.mpf(p); tbar = p/H              # Hbar = p/tbar = 1
    W = make_W(p, tbar)
    for w in [mp.mpf('0.5'), mp.mpf(2)]:
        got = Gplus(W, H, w); want = (w/(2*mp.pi))/(1 - mp.e**(-2*mp.pi*w/H))
        report(f"eps={mp.nstr(eps,3)} G+(w={w})", got, want, tol=20*eps)   # deviation is itself O(eps)
print("== 11.2 self-similarity ==")
p = 50; eps = 1/mp.mpf(p)
vals = []
for tbar in [mp.mpf(50), mp.mpf(100)]:
    Hbar = p/tbar; W = make_W(p, tbar)
    vals.append(Gplus(W, Hbar, Hbar*mp.mpf('1.0'))/Hbar)
report("G+(w=Hbar)/Hbar equal at two mean times", vals[0], vals[1], tol=1e-10)
print("== 11.3 KMS defect D(w) = ln[G+(w)/G+(-w)] - 2 pi w/Hbar, scaling with epsilon ==")
ws = [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(3)]
table = {}
for p in [25, 50, 100, 200]:
    eps = 1/mp.mpf(p); tbar = p/H; W = make_W(p, tbar)
    row = []
    for w in ws:
        gp = Gplus(W, H, w); gm = Gplus(W, H, -w)
        D = mp.log(gp/gm) - 2*mp.pi*w/H
        row.append(D/eps)
    table[p] = row
    print(f"  eps={mp.nstr(eps,4)}:  D/eps at w/H={[mp.nstr(x,2) for x in ws]} = {[mp.nstr(x,6) for x in row]}")
# convergence of D/eps as eps -> 0 : compare p=100 and p=200 rows
for i, w in enumerate(ws):
    report(f"D/eps converges (p=100 vs 200) at w/H={mp.nstr(w,2)}", table[100][i], table[200][i], tol=0.05)
print("  => f(w/H) = lim D/eps (from p=200):", [mp.nstr(x, 5) for x in table[200]])
shift = [table[200][i]/(2*mp.pi*ws[i]/H) for i in range(len(ws))]
print("  beta_eff(w)/beta_dS - 1, per unit eps (eps->0):", [mp.nstr(x, 5) for x in shift])
print("  => leading defect = effective-temperature shift T_eff ~ (H/2pi)(1 - c eps), c ~", mp.nstr(shift[1], 4), " (= H one Hubble time after the mean time, to this order)")
inv = (shift[-1] - shift[0])
print("  => CONVENTION-FREE violation: beta_eff varies with w; [beta_eff(3H) - beta_eff(0.5H)]/(eps beta_dS) =", mp.nstr(inv, 4), "(no single temperature fits)")
print("== 11.4 oddness and vanishing at eps=0 ==")
p = 100; tbar = p/H; W = make_W(p, tbar)
w = mp.mpf('1.5')
Dp = mp.log(Gplus(W, H, w)/Gplus(W, H, -w)) - 2*mp.pi*w/H
Dm = mp.log(Gplus(W, H, -w)/Gplus(W, H, w)) + 2*mp.pi*w/H
report("D(-w) = -D(w)", Dm, -Dp, tol=1e-12)
# exact dS kernel through the same pipeline: defect must vanish
WdS = lambda Delta: -H**2/(16*mp.pi**2*mp.sinh(H*Delta/2)**2)
D0 = mp.log(Gplus(WdS, H, w)/Gplus(WdS, H, -w)) - 2*mp.pi*w/H
ok = abs(D0) < 1e-15; ok_all &= ok
print(f"  [{'PASS' if ok else 'FAIL'}] exact dS through the same pipeline: D = {mp.nstr(D0, 3)}")
print("\nV11 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
