#!/usr/bin/env python3
"""
V6 -- CONJECTURE-GRADE (not derived): the one place the responsive-horizon picture produces a scale nature displays.
Jacobson's local-equilibrium construction needs a local Rindler horizon (distance c^2/a) that fits inside the
cosmological horizon (c/H); for a < c H the two horizons compete and the Unruh temperature in de Sitter becomes
T = (hbar/2 pi k_B c) sqrt(a^2 + c^2 H^2) (Narnhofer-Peter-Thirring 1996; Deser-Levin 1997).  Milgrom (1999) and
Verlinde (2016) built MOND-like phenomenology on exactly this.  Here we only tabulate the dimensional facts and the
DISCRIMINATING prediction (a_0 tracks H(z) if the apparent horizon sets it; constant if the Lambda-horizon does).
Nothing here is a derivation of a_0's coefficient.
"""
from math import pi, sqrt
c = 299792458.0; Mpc = 3.0856775814913673e22
a0_obs, a0_err = 1.20e-10, 0.02e-10   # McGaugh, Lelli, Schombert 2016 (random error; systematic ~0.24)
Om, OL = 0.315, 0.685
for H0km in [67.4, 73.0]:
    H0 = H0km*1e3/Mpc
    cH = c*H0
    print(f"H0 = {H0km} km/s/Mpc: cH0 = {cH:.3e} m/s^2;  a0_obs/cH0 = {a0_obs/cH:.3f}")
    for name, val in [("cH0/(2 pi)", cH/(2*pi)), ("cH0/6 (Verlinde)", cH/6), ("c sqrt(Lambda/3)/(2 pi)", cH*sqrt(OL)/(2*pi)), ("cH0/(2 pi) * sqrt(Omega_L)^-1", cH/(2*pi)/sqrt(OL))]:
        print(f"   {name:32s} = {val:.3e}  ratio to a0_obs = {val/a0_obs:.3f}")
print("Discriminating prediction if a0 is set by the apparent-horizon temperature: a0(z)/a0(0) = H(z)/H0")
for z in [0.5, 1.0, 2.0, 3.0]:
    Hz = sqrt(Om*(1+z)**3 + OL)
    print(f"   z={z}: H(z)/H0 = {Hz:.3f}   (Lambda-horizon alternative: 1.000)")
print("Status: CONJECTURE. Must survive: CMB acoustic peaks, cluster masses, lensing, BBN (unaffected), the RAR scatter. Not derived here.")
print("\nV6 OVERALL: INFORMATIONAL (no pass/fail)")
