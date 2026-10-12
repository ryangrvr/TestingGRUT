#!/usr/bin/env python3
"""
V15 -- The Sorkin-Johnston (SJ) state in an inflationary de Sitter slab: the one known parameter-free, temperature-free
state-selection rule, tested against GRUT III's own P2 and against the CMB.

WHY.  The baseline (12) says a distinctive GRUT claim needs a rule that selects a non-equilibrium state without assuming a
temperature.  The SJ prescription (Johnston 2009; Sorkin; Afshordi, Aslanbeigi & Sorkin 2012) is exactly such a rule: in a
bounded region O, take the Pauli-Jordan function iDelta as an operator on L^2(O, sqrt(-g) d^4x) and keep its positive
spectral part as the two-point function.  No Killing field, no temperature, no free parameter beyond the region.

MODE REDUCTION (spatially flat FRW, region = slab eta_i < eta < eta_e).  Per comoving mode k, with u_k the Bunch-Davies
mode (Wronskian-normalised) and v_k = alpha u_k + beta u_k* the SJ mode, the SJ conditions are
    |alpha|^2 - |beta|^2 = 1   and   int_O v_k^2 a^4 d eta = 0     (positive and negative eigenvectors L^2-orthogonal),
which with A = int |u|^2 a^4 d eta, B = int u^2 a^4 d eta give
    beta/alpha = -(A - lam)/conj(B),   lam = sqrt(A^2 - |B|^2),   |beta|^2 = (A - lam)/(2 lam).
de Sitter: a = -1/(H eta), H = 1; the slab ends at s_e = -eta_e (so a_e H = 1/s_e) and starts at s_i = -eta_i.
Fields: conformally coupled massless (u = -eta e^{-ik eta}/sqrt(2k)) and massless minimally coupled
(u = (1 + i k eta) e^{-ik eta}/sqrt(2k^3)), the latter standing in for the inflaton / curvature perturbation.

RESULTS (checked below).
 15.1 closed forms of A, B (incomplete gamma) agree with quadrature; Wronskians; the SJ conditions hold to 1e-150; the
      mode-reduced SJ two-point function equals the positive part of the discretised Pauli-Jordan operator (independent check).
 15.2 ULTRAVIOLET: |beta_k|^2 -> (a_e H)^2/(16 k^2) for both fields -- a power law, not faster than any power, so the SJ state
      in a sharp slab is not Hadamard (Fewster-Verch 2012), and the |beta|^2 part of its energy density grows like
      Lambda^2 with the cutoff: infinite.  This violates GRUT III's P2 (local equilibrium = Minkowski-like short distances).
 15.3 INFRARED: modes outside the horizon at the slab's END are suppressed: R_k = |v_k|^2/|u_k|^2 at eta_e behaves as
      (k/a_e H)^{1/2} (conformal) and (k/a_e H)^{3/2} (massless minimally coupled).  If the SJ region is the inflationary
      slab, every observable mode (exited ~50-60 e-folds before the end) has its power multiplied by ~(k/a_e H)^{3/2}: a
      spectral-tilt shift of +3/2 against Planck 2018's n_s - 1 = -0.035 +- 0.004.  Excluded.
 15.4 TELEOLOGY: at fixed comoving k the state depends on where the region ends in the FUTURE (through k s_e): moving the end
      of the slab changes the state of modes long before it.
CONCLUSION: the only known temperature-free, parameter-free selection rule fails GRUT III's P2 in a sharp region and, over
the inflationary slab, contradicts the CMB.  Smoothing it (Brum-Fredenhagen) restores Hadamard at the price of a free
function; choosing a region that ends today moves the IR suppression to the largest observable scales (k <~ a_0 H_0), where
cosmic variance hides it, but leaves the region choice and the teleology as inputs.  A parameter-free, Hadamard,
causal state-selection rule is exactly what is missing (12 Sec. 3, 13).
"""
import mpmath as mp
import numpy as np
ok_all = True
def check(label, cond):
    global ok_all
    print(f"  {label}:", "PASS" if cond else "FAIL"); ok_all &= bool(cond)
mp.mp.dps = 160

def In(n, k, s_e, s_i):
    """int_{eta_i}^{eta_e} e^{-2ik eta} eta^{-n} d eta, eta = -s, via the generalised incomplete gamma along -i R+"""
    return (-1)**n * (-2j*k)**(n - 1) * mp.gammainc(1 - n, -2j*k*s_e, -2j*k*s_i)
def AB(field, k, s_e, s_i):
    if field == 'conf':
        return (1/(2*k))*(1/s_e - 1/s_i), (1/(2*k))*In(2, k, s_e, s_i)
    A = (1/(2*k**3))*(mp.mpf(1)/3*(1/s_e**3 - 1/s_i**3) + k**2*(1/s_e - 1/s_i))
    B = (1/(2*k**3))*(In(4, k, s_e, s_i) + 2j*k*In(3, k, s_e, s_i) - k**2*In(2, k, s_e, s_i))
    return A, B
def u(field, k, eta):
    return (-eta*mp.exp(-1j*k*eta)/mp.sqrt(2*k)) if field == 'conf' else ((1 + 1j*k*eta)*mp.exp(-1j*k*eta)/mp.sqrt(2*k**3))
def sj(field, k, s_e, s_i):
    A, B = AB(field, k, s_e, s_i)
    lam = mp.sqrt(A**2 - abs(B)**2); r = -(A - lam)/mp.conj(B)
    al = 1/mp.sqrt(1 - abs(r)**2)
    return A, B, lam, al, r*al
def R_end(field, k, s_e, s_i):
    A, B, lam, al, be = sj(field, k, s_e, s_i); ue = u(field, k, -s_e)
    return abs(al*ue + be*mp.conj(ue))**2/abs(ue)**2, abs(be)**2

print("== 15.1 construction checks ==")
mp.mp.dps = 40
worst = 0
for n in [2, 3, 4]:
    for (k, se, si) in [(mp.mpf('0.7'), 1, 30), (mp.mpf(5), mp.mpf('0.5'), 8)]:
        quad = mp.quad(lambda eta: mp.exp(-2j*k*eta)*eta**(-n), mp.linspace(-si, -se, 400))
        worst = max(worst, abs(In(n, k, se, si) - quad)/abs(quad))
check(f"incomplete-gamma closed forms = quadrature (worst rel. {mp.nstr(worst, 2)})", worst < 1e-15)
for field in ['conf', 'mmc']:
    k, eta = mp.mpf('1.3'), mp.mpf('-0.7')
    du = mp.diff(lambda e: u(field, k, e), eta); a = -1/eta
    wr = a**2*(u(field, k, eta)*mp.conj(du) - mp.conj(u(field, k, eta))*du)
    check(f"{field}: Wronskian a^2 (u u*' - u* u') = i", abs(wr - 1j) < 1e-30)
mp.mp.dps = 160
worst = 0
for field in ['conf', 'mmc']:
    for k in ['1e-4', '0.1', '1', '100']:
        A, B, lam, al, be = sj(field, mp.mpf(k), mp.mpf(1), mp.mpf('1e8'))
        cond = al**2*B + 2*al*be*A + be**2*mp.conj(B)
        worst = max(worst, abs(cond)/A, abs(abs(al)**2 - abs(be)**2 - 1))
check(f"SJ conditions |alpha|^2 - |beta|^2 = 1 and int v^2 a^4 = 0 (worst {mp.nstr(worst, 2)})", worst < 1e-140)
# independent check: positive part of the discretised Pauli-Jordan operator equals v v^dagger
k0, se, si = 1.0, 1.0, 20.0
x, w = np.polynomial.legendre.leggauss(1600)
eta = -si + (x + 1)*(si - se)/2; wq = w*(si - se)/2*(1/eta**4)          # weights include a^4 = eta^{-4}
uu = -eta*np.exp(-1j*k0*eta)/np.sqrt(2*k0)
iD = np.outer(uu, uu.conj()) - np.outer(uu.conj(), uu)                   # i Delta_k(eta, eta')
S = np.sqrt(wq); M = S[:, None]*iD*S[None, :]                            # Hermitian (iDelta is i x real antisymmetric)
ev, V = np.linalg.eigh(M)
Pp = (V[:, ev > 0]*ev[ev > 0]) @ V[:, ev > 0].conj().T
Wsj = Pp/(S[:, None]*S[None, :])
mp.mp.dps = 40
A, B, lam, al, be = sj('conf', mp.mpf(k0), mp.mpf(se), mp.mpf(si))
vv = complex(al)*uu + complex(be)*uu.conj()
err = np.linalg.norm(Wsj - np.outer(vv, vv.conj()))/np.linalg.norm(np.outer(vv, vv.conj()))
check(f"mode-reduced SJ two-point function = positive part of the discretised Pauli-Jordan operator (rel. {err:.1e})", err < 1e-8)
mp.mp.dps = 160

print("== 15.2 ultraviolet: |beta_k|^2 tail and the Hadamard condition ==")
s_e, s_i = mp.mpf(1), mp.mpf('1e8')
for field in ['conf', 'mmc']:
    vals = [(k, R_end(field, mp.mpf(k), s_e, s_i)[1]) for k in ['100', '1000', '10000']]
    for k, b2 in vals:
        print(f"    {field}: k/(a_e H) = {k:>6}:  |beta|^2 = {mp.nstr(b2, 8)}   16 k^2 |beta|^2 (a_e H)^-2 = {mp.nstr(16*mp.mpf(k)**2*b2, 8)}")
    check(f"{field}: |beta_k|^2 -> (a_e H)^2/(16 k^2) (power law: not Hadamard)", abs(16*mp.mpf(10000)**2*vals[-1][1] - 1) < 1e-6)
# energy density of the |beta|^2 part at eta_e (conformal field): rho ~ a^-4 int d^3k/(2pi)^3 k |beta_k|^2 -> grows like Lambda^2
def rho_part(L):
    f = lambda k: k**3*R_end('conf', k, s_e, s_i)[1]/(2*mp.pi**2)
    return mp.quad(f, [mp.mpf(10), 30, 100, 300, L])
mp.mp.dps = 50
r1, r2 = rho_part(mp.mpf(1000)), rho_part(mp.mpf(4000))
slope = mp.log(r2/r1)/mp.log(4)
print(f"    conformal: int_10^Lambda k^3 |beta|^2 dk/2pi^2: Lambda = 1000 -> {mp.nstr(r1, 6)}, 4000 -> {mp.nstr(r2, 6)}; growth exponent {mp.nstr(slope, 4)}")
check("the |beta|^2 part of the SJ energy density diverges like Lambda^2 (exponent -> 2)", abs(slope - 2) < 0.05)
mp.mp.dps = 160

print("== 15.3 infrared: modes outside the horizon at the slab's end ==")
for field, expo in [('conf', 0.5), ('mmc', 1.5)]:
    ks = ['1e-12', '1e-9', '1e-6', '1e-3']
    Rs = [R_end(field, mp.mpf(k), s_e, mp.mpf('1e20'))[0] for k in ks]
    for k, R in zip(ks, Rs):
        print(f"    {field}: k/(a_e H) = {k:>6}:  R_k = |v_k|^2/|u_k|^2 at eta_e = {mp.nstr(R, 6)}")
    sl = [float(mp.log(Rs[i + 1]/Rs[i])/mp.log(mp.mpf(1000))) for i in range(len(Rs) - 1)]
    print(f"    {field}: local slopes d ln R/d ln k = {', '.join(f'{x:.4f}' for x in sl)}")
    check(f"{field}: R_k ~ (k/a_e H)^{expo} for k << a_e H", all(abs(x - expo) < 0.01 for x in sl))
print("    => SJ over the inflationary slab multiplies the power of every observable mode (k/a_e H ~ e^-50..e^-60) by ~(k/a_e H)^{3/2}:")
print("       a tilt shift of +1.5 against Planck 2018 n_s - 1 = -0.035 +- 0.004 (some 370 sigma). Excluded.")
check("tilt shift 1.5 is incompatible with |n_s - 1 - (-0.035)| < 5 x 0.0042", 1.5 + 0.035 > 5*0.0042)

print("== 15.4 dependence on the region's future boundary ==")
k = mp.mpf('1e-3')
for se in ['1', '1e-1', '1e-2']:
    R, b2 = R_end('mmc', k, mp.mpf(se), mp.mpf('1e20'))
    Ru = abs(sj('mmc', k, mp.mpf(se), mp.mpf('1e20'))[3]*u('mmc', k, mp.mpf(-1)) + sj('mmc', k, mp.mpf(se), mp.mpf('1e20'))[4]*mp.conj(u('mmc', k, mp.mpf(-1))))**2/abs(u('mmc', k, mp.mpf(-1)))**2
    print(f"    mmc, k = 1e-3: slab ends at s_e = {se:>5} (a_e H = {mp.nstr(1/mp.mpf(se),3)}): R at the FIXED time eta = -1 is {mp.nstr(Ru, 6)}")
Rs = [abs(sj('mmc', k, mp.mpf(se), mp.mpf('1e20'))[3]*u('mmc', k, mp.mpf(-1)) + sj('mmc', k, mp.mpf(se), mp.mpf('1e20'))[4]*mp.conj(u('mmc', k, mp.mpf(-1))))**2/abs(u('mmc', k, mp.mpf(-1)))**2 for se in ['1', '1e-2']]
check("extending the region into the future changes the state at an earlier fixed time (teleology)", abs(Rs[0]/Rs[1] - 1) > 0.5)

print("\nCONCLUSION: the SJ state -- the only known parameter-free, temperature-free selection rule -- is not Hadamard in a sharp")
print("region (|beta_k|^2 ~ (a_e H)^2/16k^2, divergent energy density: violates P2), suppresses modes that are superhorizon at the")
print("region's end by (k/a_e H)^{3/2} (massless field; over the inflationary slab this contradicts n_s), and depends on the")
print("region's future boundary. A parameter-free, Hadamard, causal state-selection rule is what is missing.")
print("\nV15 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
