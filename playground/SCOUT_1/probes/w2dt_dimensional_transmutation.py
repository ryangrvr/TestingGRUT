"""SCOUT-1 W2-DT: dimensional transmutation vs TC-1.

 1. large-N Gross-Neveu (2D, classically scale-invariant, dimensionless coupling).
    Euclidean effective potential per flavour with hard cutoff Lam:
      V(s) = s^2/(2 lam) - (1/4pi) [ (Lam^2+s^2) ln(Lam^2+s^2) - s^2 ln s^2 - Lam^2 ln Lam^2 ]
    gap: 1/lam = (1/2pi) ln(1 + Lam^2/m^2)  =>  m = Lam / sqrt(exp(2pi/lam) - 1)
    (a) m is RG-invariant: many (Lam, lam(Lam)) pairs on one trajectory give the same m;
    (b) dimensionless observables in units of m are cutoff/coupling independent (no dimensionless input left);
    (c) the absolute m needs one boundary datum (lam at a reference Lam) -- the G-orbit coordinate.
 2. one-loop QCD: Lambda = mu exp(-2 pi / (b0 alpha_s(mu))); input alpha_s(M_Z) = 0.118 at M_Z = 91.19 GeV;
    hierarchy Lambda / M_Pl from a dimensionless UV boundary coupling (exponential sensitivity).
 3. Coleman-Weinberg (massless scalar QED): <phi> traded for lambda; m_S^2 / m_V^2 = (3 e^2) / (8 pi^2)
    fixed by the remaining dimensionless coupling; absolute <phi> = renormalization-point datum.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def V(s, lam, Lam):
    s2 = s * s
    L2 = Lam * Lam
    # stable form of (L2+s2)ln(L2+s2) - s2 ln s2 - L2 ln L2 = L2 ln(1+s2/L2) + s2 ln(1+L2/s2)
    t = L2 * np.log1p(s2 / L2) + (s2 * np.log1p(L2 / s2) if s2 > 0 else 0.0)
    return s2 / (2 * lam) - t / (4 * np.pi)


def gap_m(lam, Lam):
    return Lam / np.sqrt(np.expm1(2 * np.pi / lam))


print("=== 1a. Gross-Neveu: one RG trajectory, many (cutoff, coupling) pairs ===")
m_target = 1.0
print("   Lam        lam(Lam)     m(gap eq)    m(minimize V)   Delta V / m^2")
for Lam in (1e1, 1e2, 1e4, 1e8, 1e16):
    lam = 2 * np.pi / np.log1p(Lam ** 2 / m_target ** 2)     # coupling that puts this cutoff on the m=1 trajectory
    m_gap = gap_m(lam, Lam)
    res = minimize_scalar(lambda s: V(s, lam, Lam), bounds=(1e-3, 10.0), method="bounded", options={"xatol": 1e-12})
    dV = (V(res.x, lam, Lam) - V(0.0, lam, Lam)) / res.x ** 2
    print(f"  {Lam:8.0e}   {lam:.6f}   {m_gap:.10f}   {res.x:.8f}      {dV:+.8f}")
print("   continuum value of Delta V/m^2 = -1/(4 pi) = %.8f  (pure number: no dimensionless input survives)" % (-1 / (4 * np.pi)))

print("\n=== 1b. the coupling is NOT a free dimensionless parameter: every lam is the same theory in other units ===")
Lam = 1e6
for lam in (0.2, 0.3, 0.5, 1.0):
    m = gap_m(lam, Lam)
    res = minimize_scalar(lambda s: V(s, lam, Lam), bounds=(m * 1e-3, m * 10), method="bounded", options={"xatol": 1e-14 * m})
    dV = (V(res.x, lam, Lam) - V(0.0, lam, Lam)) / res.x ** 2
    print(f"   lam={lam:4.2f}: m/Lam = {m / Lam:.3e}   Delta V/m^2 = {dV:+.8f}")
print("   classical theory: one dimensionless parameter (lam). quantum theory: zero dimensionless parameters,")
print("   one scale m = the coordinate along the (anomalous) scaling orbit. Absolute m requires lam at a reference Lam.")

print("\n=== 2. one-loop QCD ===")
MZ, aZ = 91.1876, 0.118
for nf in (5,):
    b0 = 11 - 2 * nf / 3
    Lam = MZ * np.exp(-2 * np.pi / (b0 * aZ))
    print(f"   n_f={nf}: b0={b0:.3f}, Lambda_1loop = {Lam * 1e3:.1f} MeV  (inputs: alpha_s = {aZ} AT mu = M_Z = {MZ} GeV)")
MPl = 1.22e19
b0 = 11 - 2 * 5 / 3
print("   hierarchy from a dimensionless UV boundary coupling at M_Pl (n_f = 5 throughout, illustrative):")
for aPl in (0.020, 0.025, 0.030):
    print(f"     alpha_s(M_Pl) = {aPl}: Lambda/M_Pl = exp(-2pi/(b0 alpha)) = {np.exp(-2 * np.pi / (b0 * aPl)):.3e}")
print("   d ln(Lambda/M_Pl) / d alpha = 2 pi / (b0 alpha^2): a 1% change in alpha_UV moves the hierarchy by a factor",
      f"{np.exp(2 * np.pi / (b0 * 0.025) * 0.01):.2f}")

print("\n=== 3. Coleman-Weinberg (massless scalar QED, 1-loop, lambda ~ e^4) ===")
for e in (0.1, 0.3, 0.5):
    print(f"   e = {e}: m_S^2/m_V^2 = 3 e^2/(8 pi^2) = {3 * e ** 2 / (8 * np.pi ** 2):.5f}   (lambda(<phi>) = 33 e^4/(8 pi^2) eliminated)")
print("   <phi> itself = the renormalization point at which lambda takes its CW value: a boundary datum.")
