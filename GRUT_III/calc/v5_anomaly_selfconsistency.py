#!/usr/bin/env python3
"""
V5 -- The trace anomaly as the medium's 'alpha leg', made quantitative with Standard-Model content.
 5.1 Anomaly coefficients a = (N_s + 11 N_D + 62 N_V)/360, c = (N_s + 6 N_D + 12 N_V)/120 (N_D = Dirac = 2 Weyl).
     SM: N_s = 4, N_W = 45 (48 with nu_R), N_V = 12.  Check: conformal scalar <T^mu_mu>_dS = -3 a H^4/(2 pi^2) = -H^4/(240 pi^2).
     C_T = 40 c/pi^4 -> scalar 1/(3 pi^4) (consistent with V1).
 5.2 Starobinsky (1980) self-consistent de Sitter: trace of G_mn = 8 pi G <T_mn> with W = 0, E_4 = 24 H^4 gives
     H^2 = pi M_P^2 / a  (M_P^2 = 1/G).  With SM content H ~ 1.07 M_P: Planckian, NOT the observed inflation.
     Species count needed for H_inf = 1e13 GeV: a ~ 3e12.  (The 'closed' anomaly-driven GRUT is falsified at SM content;
     Starobinsky's R^2 coefficient is a supplied scale.)
"""
from math import pi, sqrt
ok_all = True
def anomaly(Ns, NW, NV):
    ND = NW/2
    return (Ns + 11*ND + 62*NV)/360, (Ns + 6*ND + 12*NV)/120
for label, Ns, NW, NV in [("SM (45 Weyl)", 4, 45, 12), ("SM + 3 nu_R (48 Weyl)", 4, 48, 12), ("one real conformal scalar", 1, 0, 0)]:
    a, c = anomaly(Ns, NW, NV)
    print(f"  {label:28s}: a = {a:.5f}  c = {c:.5f}  (120c = {120*c:.1f}, 360a = {360*a:.1f})  C_T = 40c/pi^4 = {40*c/pi**4:.6f}")
a_s, c_s = anomaly(1, 0, 0)
ok = abs(3*a_s/(2*pi**2) - 1/(240*pi**2)) < 1e-15 and abs(40*c_s/pi**4 - 1/(3*pi**4)) < 1e-15
ok_all &= ok
print(f"  scalar checks: -3aH^4/2pi^2 = -H^4/240pi^2 and C_T = 1/(3 pi^4): [{'PASS' if ok else 'FAIL'}]")
print("== 5.2 self-consistent de Sitter from the anomaly alone ==")
MP = 1.220890e19  # GeV, non-reduced Planck mass
for label, Ns, NW, NV in [("SM (45 Weyl)", 4, 45, 12), ("SM + 3 nu_R", 4, 48, 12)]:
    a, c = anomaly(Ns, NW, NV)
    H = sqrt(pi/a)*MP
    print(f"  {label}: H_sc = sqrt(pi/a) M_P = {H/MP:.3f} M_P = {H:.3e} GeV")
H_inf_max = 4.7e13  # GeV, from r < 0.036
a_needed = pi*(MP/H_inf_max)**2
print(f"  observed bound H_inf < {H_inf_max:.1e} GeV (r<0.036)  =>  anomaly coefficient needed a > {a_needed:.2e}  (~{360*a_needed/62:.1e} vector species)")
print("  => anomaly-driven self-consistent dS is Planckian at SM content: GRUT III inherits Starobinsky's conclusion that the R^2 scale is supplied.")
print("\nV5 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
