#!/usr/bin/env python3
"""
V6 -- NUMERICAL COINCIDENCE ONLY (not derived; no mechanism within GRUT III).
The MOND acceleration scale a_0 ~ 1.2e-10 m/s^2 is numerically close to c H_0/(2 pi) ~ 1.1e-10 m/s^2.  This
coincidence was noted by Milgrom (1999) and used by Verlinde (2016); GRUT III supplies NO mechanism for it.
In particular, Jacobson's local-equilibrium derivation does NOT fail for a_galaxy < c H: the local Rindler
construction uses an auxiliary boost at an arbitrary point whose kappa cancels out of the Clausius relation,
and is unrelated to a star's orbital acceleration.  Nothing below bears on that derivation.

Here we only tabulate (i) the dimensional-coincidence table and (ii) the one DISCRIMINATING observable that
separates the two ways the coincidence could be realised: a_0 ∝ H(z) (apparent horizon) vs a_0 = const
(Lambda-horizon).  The a_0 ∝ H(z) test has been proposed before (Milgrom 1999 and later authors); it is
listed here as a discriminator, not as a GRUT III prediction.  Nothing here is a derivation of a_0's coefficient.
"""
from math import pi, sqrt
c = 299792458.0; Mpc = 3.0856775814913673e22
a0_obs, a0_err = 1.20e-10, 0.02e-10   # McGaugh, Lelli, Schombert 2016 (random error; systematic ~0.24)
Om, OL = 0.315, 0.685
print("Status: numerical coincidence (Milgrom 1999, Verlinde 2016); no mechanism within GRUT III.")
print("== dimensional-coincidence table ==")
for H0km in [67.4, 73.0]:
    H0 = H0km*1e3/Mpc
    cH = c*H0
    print(f"H0 = {H0km} km/s/Mpc: cH0 = {cH:.3e} m/s^2;  a0_obs/cH0 = {a0_obs/cH:.3f}")
    for name, val in [("cH0/(2 pi)", cH/(2*pi)), ("cH0/6 (Verlinde)", cH/6), ("c sqrt(Lambda/3)/(2 pi)", cH*sqrt(OL)/(2*pi)), ("cH0/(2 pi) * sqrt(Omega_L)^-1", cH/(2*pi)/sqrt(OL))]:
        print(f"   {name:32s} = {val:.3e}  ratio to a0_obs = {val/a0_obs:.3f}")
print("== discriminator: a0 ∝ H(z) (apparent horizon) vs a0 = const (Lambda-horizon) ==")
print("  (this test has been proposed before; it is not a GRUT III prediction)")
for z in [0.5, 1.0, 2.0, 3.0]:
    Hz = sqrt(Om*(1+z)**3 + OL)
    print(f"   z={z}: H(z)/H0 = {Hz:.3f}   (Lambda-horizon alternative: 1.000)")
print("Status: numerical coincidence (Milgrom 1999, Verlinde 2016); no mechanism within GRUT III; the a0 ∝ H(z) test has been proposed before.")
print("\nV6 OVERALL: INFORMATIONAL (no pass/fail)")
