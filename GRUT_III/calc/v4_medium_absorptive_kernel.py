#!/usr/bin/env python3
"""
V4 -- The vacuum medium's absorptive (dissipative) response is DERIVED from the particle spectrum:
by the optical theorem, Im Pi(s) for the graviton self-energy equals the total rate at which a
timelike graviton of invariant mass sqrt(s) pair-produces Standard-Model quanta.  Nothing is supplied:
the 'spectral edge data' (edge locations 2 m_i, edge exponents 5/2, 3/2, 1/2 for scalars, fermions,
vectors; weights 1 : 3 : 12 per real scalar : Weyl fermion : vector) are fixed by masses, spins and statistics.

Conventions: g = eta + kappa h, kappa^2 = 32 pi G, canonically normalised TT graviton, coupling -(kappa/2) h_mn T^mn,
propagator i/(k^2 - Pi).  |Im Pi(s)| = (1/2) S_f Int dPhi_2 |M|^2 (S_f = 1/2 for identical particles), averaged over
the 5 spin-2 polarisations (isotropic).  Polarisation sum: Pi_{mn,rs} = 1/2(P_mr P_ns + P_ms P_nr) - 1/3 P_mn P_rs,
P_mn = -g_mn + q_m q_n/M^2 (mostly-minus metric).

Checks:
 4.1 real scalar:   Sum_pol |eps.T|^2 = 8 p^4/3 exactly  =>  threshold law beta^5;  Im Pi = G s^2 beta^5/120.
 4.2 Dirac fermion: Sum_spins,pol |eps.T|^2 from explicit Dirac algebra;  massless ratio Dirac/real-scalar = 6
     (the C_T / Weyl-anomaly ratio);  threshold law beta^3 (1 + 8 r/3), r = m^2/s, extracted numerically.
 4.3 massless vector (two photons, physical polarisations): ratio vector/real-scalar = 12.
 4.4 Hence Im Pi(s) = c_eff(s) G s^2 with c_eff = (1/120) Sum_i w_i F_i(4 m_i^2/s): the Weyl-anomaly coefficient,
     resolved in the timelike invariant mass sqrt(s) (NOT in the frequency of a propagating GW: a free on-shell GW
     has s = 0, so no species is open regardless of its frequency).  Standard-Model staircase tabulated; it is
     ILLUSTRATIVE: the massive-vector threshold polynomial is NOT computed (leading beta only), and the hadronic region
     0.3-2 GeV is schematic.  In FRW the off-shell effect is O(c G H^2 omega), i.e. far below anything observable.
"""
import numpy as np
from math import pi, sqrt
ok_all = True
g = np.diag([1., -1., -1., -1.])
I2 = np.eye(2, dtype=complex); Z2 = np.zeros((2, 2), dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
gam = [np.block([[I2, Z2], [Z2, -I2]])] + [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
def slash(p):  # p: upper-index components
    pl = g @ p
    return sum(gam[m]*pl[m] for m in range(4))
def pol_sum_tensor(M):
    q = np.array([M, 0, 0, 0.]); ql = g @ q
    P = -g + np.outer(ql, ql)/M**2
    return 0.5*(np.einsum('mr,ns->mnrs', P, P) + np.einsum('ms,nr->mnrs', P, P)) - np.einsum('mn,rs->mnrs', P, P)/3
def contract(Pi, T):
    return np.einsum('mnrs,mn,rs->', Pi, T, np.conj(T)).real

def scalar_amp2(M, m):
    E = M/2; p = sqrt(max(E*E - m*m, 0.))
    p1 = np.array([E, 0, 0, p]); p2 = np.array([E, 0, 0, -p])
    T = np.outer(p1, p2) + np.outer(p2, p1)          # g^{mn} term drops against traceless eps
    return contract(pol_sum_tensor(M), T), p

def fermion_amp2(M, m):
    E = M/2; p = sqrt(max(E*E - m*m, 0.))
    p1 = np.array([E, 0, 0, p]); p2 = np.array([E, 0, 0, -p]); k = p1 - p2
    V = [[0.25*(gam[mu]*k[nu] + gam[nu]*k[mu]) for nu in range(4)] for mu in range(4)]
    A = slash(p1) + m*np.eye(4); B = slash(p2) - m*np.eye(4)
    Vbar = [[gam[0] @ V[r][s].conj().T @ gam[0] for s in range(4)] for r in range(4)]
    Ttens = np.zeros((4, 4, 4, 4), complex)
    for mu in range(4):
        for nu in range(4):
            AV = A @ V[mu][nu] @ B
            for r in range(4):
                for s in range(4):
                    Ttens[mu, nu, r, s] = np.trace(AV @ Vbar[r][s])
    Pi = pol_sum_tensor(M)
    return np.einsum('mnrs,mnrs->', Pi, Ttens).real, p

def photon_amp2(M):
    E = M/2
    k1 = np.array([E, 0, 0, E]); k2 = np.array([E, 0, 0, -E])
    pols = [np.array([0, 1, 0, 0.]), np.array([0, 0, 1, 0.])]
    Pi = pol_sum_tensor(M); tot = 0.
    for e1 in pols:
        for e2 in pols:
            F1 = 1j*(np.outer(k1, e1) - np.outer(e1, k1))   # F^{mu nu} (upper)
            F2 = 1j*(np.outer(k2, e2) - np.outer(e2, k2))
            # T^{mn} = F1^{m a} F2^{n}_a + F2^{m a} F1^{n}_a  (index lowered with g); trace term drops
            T = F1 @ g @ F2.T + F2 @ g @ F1.T
            tot += contract(Pi, T)
    return tot

print("== 4.1 real scalar ==")
M = 1.0
for m in [0.0, 0.2, 0.4]:
    val, p = scalar_amp2(M, m)
    ok = abs(val - 8*p**4/3) < 1e-12
    ok_all &= ok
    print(f"  m={m}: Sum_pol|eps.T|^2 = {val:.10f}  vs 8p^4/3 = {8*p**4/3:.10f}  [{'PASS' if ok else 'FAIL'}]")
print("  => threshold law beta^5 (beta = 2p/M); Im Pi_S = (1/2)(1/2)(kappa^2/4)(beta/8pi)(8p^4/3)/5 = kappa^2 beta^5 s^2/(3840 pi) = G beta^5 s^2/120")

print("== 4.2 Dirac fermion ==")
S0, _ = scalar_amp2(M, 0.0)
F0, _ = fermion_amp2(M, 0.0)
ratio = F0/(0.5*S0)            # identical-particle 1/2 for the real scalar
ok = abs(ratio - 6) < 1e-9; ok_all &= ok
print(f"  massless: Sum|M|^2_F / (1/2 Sum|M|^2_S) = {ratio:.10f}  (expected 6 = C_T^Dirac/C_T^scalar)  [{'PASS' if ok else 'FAIL'}]")
print("  threshold shape: F_F(r) := beta * Sum|M|^2_F(m)/Sum|M|^2_F(0)  vs  beta^3 (1 + 8r/3), r = m^2/s")
maxdev = 0.
for m in [0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]:
    val, p = fermion_amp2(M, m); beta = 2*p/M; r = m*m/M**2
    FF = beta*val/F0; pred = beta**3*(1 + 8*r/3)
    maxdev = max(maxdev, abs(FF - pred))
    print(f"   r={r:.4f}: numeric {FF:.10f}  formula {pred:.10f}")
ok = maxdev < 1e-9; ok_all &= ok
print(f"  max |numeric - beta^3(1+8r/3)| = {maxdev:.2e}  [{'PASS' if ok else 'FAIL'}]  (HLZ 1999 form reproduced from scratch)")

print("== 4.3 massless vector (photon pair) ==")
G0 = photon_amp2(M)
ratio_v = (0.5*G0)/(0.5*S0)
ok = abs(ratio_v - 12) < 1e-9; ok_all &= ok
print(f"  (1/2)Sum|M|^2_gamma / (1/2)Sum|M|^2_S = {ratio_v:.10f}  (expected 12)  [{'PASS' if ok else 'FAIL'}]")

print("== 4.4 Im Pi(s) = c_eff(s) G s^2 : Standard-Model staircase in the timelike invariant mass sqrt(s) ==")
print("  ILLUSTRATIVE ONLY: massive-vector (W, Z) threshold polynomial NOT computed (leading beta only);")
print("  hadronic region 0.3-2 GeV schematic (3 real pions as proxy, perturbative quarks+gluons above 2 GeV).")
GeV = 1.0
species = [  # name, weight w (120 c), mass [GeV], threshold exponent nu (F ~ beta^(2nu)), kind
    ("photon", 12, 0.0, 0.5), ("nu1", 3, 1e-12, 1.5), ("nu2", 3, 8.7e-12, 1.5), ("nu3", 3, 5.0e-11, 1.5),
    ("e", 6, 0.000511, 1.5), ("mu", 6, 0.10566, 1.5), ("tau", 6, 1.77686, 1.5),
    ("pi (3 real, hadronic proxy <2 GeV)", 3, 0.1396, 2.5),
    ("u,d,s + gluons (perturbative above 2 GeV)", 18*3 + 96, 1.0, 0.5),
    ("c", 18, 1.27, 1.5), ("b", 18, 4.18, 1.5), ("t", 18, 172.7, 1.5),
    ("W+,W- (12+1 each)", 26, 80.38, 0.5), ("Z (12+1)", 13, 91.19, 0.5), ("H", 1, 125.25, 2.5)]
def F(kind_nu, r):
    if 4*r >= 1: return 0.0
    beta = sqrt(1 - 4*r)
    if kind_nu == 2.5: return beta**5
    if kind_nu == 1.5: return beta**3*(1 + 8*r/3)
    return beta  # vector: leading threshold power only (massive-vector polynomial not computed here)
def c_eff(sqrt_s):
    s = sqrt_s**2; tot = 0.
    for name, w, m, nu in species:
        if name.startswith("pi") and sqrt_s > 2.0: continue
        tot += w*F(nu, m*m/s)
    return tot/120.
print("  sqrt(s)           c_eff      120 c_eff   (c_SM(45 Weyl) = 283/120 = 2.358 at high energy)")
for sq in [1e-12, 1e-10, 1e-9, 1e-6, 1e-3, 2e-3, 1e-2, 0.3, 1.0, 3.0, 10., 200., 1e3, 1e4]:
    c = c_eff(sq); print(f"  {sq:9.1e} GeV   {c:8.4f}   {120*c:8.2f}")
# The threshold variable is the TIMELIKE invariant mass sqrt(s) of an off-shell graviton, not a GW frequency.
# For scale: hbar * 2 pi * 100 Hz = 4.1e-13 eV = 4.1e-22 GeV.
c_ligo = c_eff(4.1e-22)
print(f"  for a timelike graviton of invariant mass sqrt(s) ~ 4e-22 GeV: c_eff = {c_ligo:.4f} (photon only; neutrinos closed by their masses)")
print("  a free on-shell GW has s = 0, so no species is open regardless of frequency: the above is NOT a statement about LIGO-band GWs.")
G_SI = 6.67430e-11; hbar = 1.054571817e-34; c_l = 299792458.0
H0 = 2.18e-18
om = 2*pi*100.
# order-of-magnitude FRW absorption: s ~ H0^2, Gamma ~ Im Pi/(2 omega) = c G hbar H0^4/(2 omega c^5)   [SI: G hbar/c^5 = t_P^2]
tP2 = G_SI*hbar/c_l**5
Gam = c_ligo*tP2*H0**4/(2*om)
print(f"  on-shell flat-space absorption of a GW: Im Pi(s=0) = 0 exactly.  FRW order-of-magnitude estimate (s ~ H0^2 off-shellness,")
print(f"  photon channel open) Gamma ~ c_eff t_P^2 H0^4/(2 omega) = {Gam:.1e} s^-1 at omega = 2 pi x 100 Hz")
print("  => the medium's absorptive response in sqrt(s) is derived (thresholds fixed by masses/spins), and unobservably weak for GWs: consistent with GW170817 trivially.")
print("\nV4 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
