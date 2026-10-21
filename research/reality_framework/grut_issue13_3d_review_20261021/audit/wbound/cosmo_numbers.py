#!/usr/bin/env python3
"""Track wbound: magnitudes of the FRW corrections to the virial bound at wall-radiation parity.  [E]
Inputs from Claude's packet (c7b8480): v = 1 TeV, lambda_phi = 0.01, T_par = 5.8-6.4 keV,
rho_tot(T_dom) = 3.7e-21 GeV^4 (verification/v_cosmo/VERIFY.md).  M_P = 2.435323e18 GeV.
Writes results/cosmo_numbers.json.
"""
import json
import math
import os

MP = 2.435323e18
v, lam = 1e3, 0.01
m_phi = math.sqrt(2 * lam) * v
g_rad = 2 + 7 / 8 * 2 * 3.044 * (4 / 11) ** (4 / 3)     # photons + neutrinos (g_rho after e+e-)
out = {"g_rho_late": g_rad}
for Tpar_keV in (5.8, 6.4):
    T = Tpar_keV * 1e-6
    rho_r = math.pi**2 / 30 * g_rad * T**4
    rho_tot = 2 * rho_r                      # walls = radiation at parity
    H = math.sqrt(rho_tot / 3) / MP
    rho_h = rho_r
    out[f"T_par={Tpar_keV}keV"] = dict(
        rho_rad_GeV4=rho_r, rho_tot_GeV4=rho_tot, H_GeV=H, H_over_m_phi=H / m_phi,
        boundary_bound_H_v_sqrt2_over_rho=H * v * math.sqrt(2 / rho_h),
        barrier_lam_v4_over_4_over_mean=lam * v**4 / 4 / rho_h,
        eps_star_over_mean=0.5876775 * lam * v**4 / rho_h,
        wall_core_density_lam_v4_over_2_over_mean=lam * v**4 / 2 / rho_h,
        Hinv_over_wall_width=(1 / H) / (math.sqrt(2 / lam) / v),
        xi_R_v2_over_rho_estimate=0.033 / 12 * (rho_tot / MP**2) * v**2 / rho_h)
out["note"] = ("boundary_bound: |[H a^4 g0]|/(rho a^4) <= H v sqrt(2/rho_h) per e-fold (g0 = <phi phidot>, |phi| ~ v); "
               "the addot term is of the same order per e-fold.  eps_star = 0.5877 lam v^4 is the homogeneous "
               "threshold for wbar > 1/3 (zero_d.py).  xi term: |xi_phi| <= 0.033/12 (packet), R ~ rho_tot/M_P^2.")
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "results", "cosmo_numbers.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
