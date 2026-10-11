#!/usr/bin/env python3
"""
V3 -- The responsive medium's transport data are fixed by (hbar, c, G, k_B): nothing is supplied.
 3.1  Membrane-paradigm SURFACE shear viscosity of the horizon eta = c^3/(16 pi G) [units kg/s = Pa s m: a 2D membrane
      viscosity, NOT a bulk Pa s], surface entropy density s = k_B c^3/(4 G hbar) [J K^-1 m^-2]  =>  eta/s = hbar/(4 pi k_B)
      (the KSS value).  The horizon is a 2D fluid (Damour 1979; Thorne-Price-Macdonald 1986); this does NOT say that empty
      space is a bulk fluid with an ordinary shear viscosity.
 3.2  Thermal (beta/2pi) time scale tau_thermal = hbar/(2 pi k_B T) = 1/kappa:  for de Sitter T = hbar H/(2 pi k_B c)...
      tau_thermal = 1/H; for a solar-mass black hole tau_thermal = 4GM/c^3; for Rindler tau_thermal = c/a.  This is a DIFFERENT
      quantity from the 1/(4H) of V1, which is the asymptotic decay rate of the pointwise <T T> correlator; that correlator
      diverges as t^-8 at short times, so no finite normalised autocorrelation time exists without smearing.
 3.3  de Sitter static patch identities (exact rationals): E = rho_Lambda V = 1/(2GH) = T S  =>  free energy F = 0;
      |E_Komar| = 1/(GH) = (1/2) N k_B T with N = A/l_P^2  (holographic equipartition).
"""
import sympy as sp
from math import pi
ok_all = True
# CODATA 2018
c = 299792458.0; G = 6.67430e-11; hbar = 1.054571817e-34; kB = 1.380649e-23; Msun = 1.98847e30
H0 = 67.4e3/3.0856775814913673e22   # s^-1 (Planck 2018)
print("== 3.1 horizon (membrane) fluid transport ==")
eta = c**3/(16*pi*G); s_area = kB*c**3/(4*G*hbar)
print(f"  eta  = c^3/(16 pi G) = {eta:.4e} kg s^-1  (= Pa s m: membrane-paradigm SURFACE shear viscosity of the horizon, a 2D fluid;")
print(f"         not a bulk Pa s viscosity of empty space)")
print(f"  s    = k_B c^3/(4 G hbar) = {s_area:.4e} J K^-1 m^-2  (surface entropy density)")
print(f"  eta/s = {eta/s_area:.4e} K s   vs hbar/(4 pi k_B) = {hbar/(4*pi*kB):.4e} K s  ->",
      "PASS" if abs(eta/s_area - hbar/(4*pi*kB)) < 1e-9*eta/s_area else "FAIL")
ok_all &= abs(eta/s_area - hbar/(4*pi*kB)) < 1e-9*eta/s_area
print("== 3.2 thermal (beta/2pi) time scale tau_thermal = hbar/(2 pi k_B T) = 1/kappa ==")
T_dS = hbar*H0/(2*pi*kB)
tau = hbar/(2*pi*kB*T_dS)
print(f"  T_dS(H0) = {T_dS:.3e} K ;  tau_thermal = hbar/(2 pi k_B T) = {tau:.4e} s = 1/H0 = {1/H0:.4e} s ->", "PASS" if abs(tau*H0-1) < 1e-12 else "FAIL")
ok_all &= abs(tau*H0-1) < 1e-12
for M in [1, 10, 1e6]:
    TH = hbar*c**3/(8*pi*G*M*Msun*kB); tauBH = hbar/(2*pi*kB*TH)
    print(f"  BH M={M:g} Msun: T_H={TH:.3e} K, tau_thermal {tauBH:.4e} s = 4GM/c^3 = {4*G*M*Msun/c**3:.4e} s")
print("  NOTE: tau_thermal = 1/kappa is the beta/2pi scale of the Planckian occupation. It is a DIFFERENT quantity from the")
print("  1/(4H) of V1, which is the asymptotic (large-|t|) decay rate of the pointwise <T00(t) T00(0)> correlator. That pointwise")
print("  correlator diverges as t^-8 as t -> 0, so no finite normalised autocorrelation time exists without smearing;")
print("  neither number is a 'memory time' of the medium in the usual sense.")
print("== 3.3 de Sitter static-patch identities (hbar=c=k_B=1, exact) ==")
Gs, Hs = sp.symbols('G H', positive=True)
A = 4*sp.pi/Hs**2; S = A/(4*Gs); T = Hs/(2*sp.pi)
rhoL = 3*Hs**2/(8*sp.pi*Gs); V = sp.Rational(4,3)*sp.pi/Hs**3
E = sp.simplify(rhoL*V); TS = sp.simplify(T*S)
print(f"  E = rho_Lambda V = {E},  T S = {TS},  F = E - TS = {sp.simplify(E-TS)} ->", "PASS" if sp.simplify(E-TS) == 0 else "FAIL")
ok_all &= sp.simplify(E-TS) == 0
N = A/Gs; EK = sp.simplify(2*rhoL*V)   # |Komar| for Lambda: (rho+3p) = -2 rho
print(f"  |E_Komar| = {EK},  (1/2) N T = {sp.simplify(N*T/2)} ->", "PASS" if sp.simplify(EK - N*T/2) == 0 else "FAIL")
ok_all &= sp.simplify(EK - N*T/2) == 0
print("\nV3 OVERALL:", "PASS" if ok_all else "FAIL")
raise SystemExit(0 if ok_all else 1)
