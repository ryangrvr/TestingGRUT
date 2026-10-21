#!/usr/bin/env python3
"""Independent cosmology mapping checks (verifier v_charge): Y_B, T_n(r), H/m_X, KZ length, n_Q,max, xi_BEC.

Inputs (quoted, with provenance): omega_b = 0.02233 (Mossa et al. 2020 via BRIEF), T_CMB = 2.7255 K,
G_N = 6.67430e-11 (CODATA 2018), 1 GeV/c^2 = 1.78266192e-27 kg, hbar c = 0.1973269804 GeV fm,
m_p = 0.93827209 GeV, g_*s = 106.75 (SM, no nu_R, no hidden dof), reduced M_P = 2.435323e18 GeV (BRIEF).
"""
import json, os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_cosmo.json")
g, lp, lX, v = 0.1, 0.01, 0.001, 1000.0
MP = 2.435323e18
gs = 106.75

# ---- Y_B from omega_b, computed from scratch ----
H100 = 100 * 1e3 / 3.0856775814913673e22  # s^-1
G = 6.67430e-11
rho_c_h2_kg_m3 = 3 * H100**2 / (8 * np.pi * G)
GeV_kg = 1.78266192e-27
rho_c_h2_GeV_cm3 = rho_c_h2_kg_m3 / GeV_kg * 1e-6
kB_GeV = 8.617333262e-14  # GeV/K
hbarc_cm = 0.1973269804e-13  # GeV cm
T0 = 2.7255 * kB_GeV
s0 = (2 * np.pi**2 / 45) * 3.90909 * (T0 / hbarc_cm) ** 3  # g_s0 = 43/11
omega_b = 0.02233
res_YB = {}
for name, mb in [("m_p", 0.93827209), ("m_b=0.938 (Codex)", 0.938), ("m_u", 0.93149410)]:
    nb = omega_b * rho_c_h2_GeV_cm3 / mb
    res_YB[name] = dict(n_b_cm3=nb, Y_B=nb / s0)
# using module's quoted inputs exactly
YB_module_inputs = 0.02233 * 1.05375e-5 / 0.938272 / 2891.2

# ---- T_n(r) ----
YB = res_YB["m_p"]["Y_B"]
s_coef = 2 * np.pi**2 / 45 * gs
n_c = np.sqrt(8 * lX) * (lp / (2 * g**2)) ** 1.5         # v^3
thresholds = dict(n_c=n_c, n_t=0.04010261923326426, n_s1=0.04276668660663895,
                  n_star=np.sqrt(lp) * (lp / lX) ** 0.25)
Tn = {}
for r in [1.0, 0.64, 0.1, 0.01]:
    Tn[f"r={r}"] = {k: (val / (r * YB * s_coef)) ** (1 / 3) * v / 1e3 for k, val in thresholds.items()}  # TeV
T_prod = (0.0220085210770737 / (YB * s_coef)) ** (1 / 3) * v / 1e3


def Hubble(T):
    return np.sqrt(np.pi**2 * gs / 90) * T**2 / MP


Tn1 = Tn["r=1.0"]["n_c"] * 1e3  # GeV
H_over_mX = Hubble(Tn1) / (g * v)
# KZ estimate (mean field nu=1/2, z=1): xi = (m tau_Q)^{1/3}/m with m = sqrt(lp) v, tau_Q = 1/(3H)
m = np.sqrt(lp) * v
xiKZ = {}
for r in [1.0, 0.01]:
    T = Tn[f"r={r}"]["n_c"] * 1e3
    tau = 1 / (3 * Hubble(T))
    xi = (m * tau) ** (1 / 3)  # in ell = 1/m
    RH = 1 / Hubble(T) * m      # Hubble radius in ell
    xiKZ[f"r={r}"] = dict(xi_KZ_ell=xi, R_H_ell=RH, domains_per_hubble_volume=(RH / xi) ** 3,
                          charge_per_KZ_volume=n_c * v**3 * (xi / m) ** 3)


# ---- Bose gas, complex scalar (particle + antiparticle, 1 dof each) ----
def nQ(mu, M, T=1.0):
    def integ(p):
        E = np.sqrt(p * p + M * M)
        a = (p * p + M * M - mu * mu) / (E + mu)  # E - mu without cancellation
        return p * p * (1 / np.expm1(a / T) - 1 / np.expm1((E + mu) / T))
    val, err = quad(integ, 0, 60 * T, limit=400, epsabs=0, epsrel=1e-12)
    return val / (2 * np.pi**2)


PiX = g**2 / 12 + lX / 3
Mt = np.sqrt(PiX)
nmax = nQ(Mt, Mt)
# cross-check normalisation: massless, small mu -> mu T^2/3
chk_massless = nQ(1e-3, 0.0) / (1e-3 / 3)
xi_BEC_coef = (YB * s_coef / nmax) ** (1 / 3)
# relic_yield table spot check: xi=7.58e-4, r=1 condensed fraction
nreq = YB * s_coef / (7.58e-4) ** 3
cond_frac = 1 - nmax / nreq
# T_c0/T_n ratio
ratio = 1.5575 / (Tn["r=1.0"]["n_c"])

out = dict(rho_c_h2_GeV_cm3=rho_c_h2_GeV_cm3, s0_cm3=s0, Y_B=res_YB, Y_B_with_module_inputs=YB_module_inputs,
           s_SM_coef=s_coef, nQ_over_T3_coef=YB * s_coef, n_c=n_c, T_n_TeV=Tn, T_prod_TeV=T_prod,
           H_over_mX_at_Tn_r1=H_over_mX, KZ=xiKZ, PiX=PiX, M_over_T=Mt, nQmax_over_T3=nmax,
           leading_MT2_over_3=Mt / 3, rel_diff_exact_vs_leading=(nmax - Mt / 3) / (Mt / 3),
           massless_norm_check_ratio=chk_massless, xi_BEC_coefficient=xi_BEC_coef,
           cond_fraction_xi7p58em4_r1=cond_frac, Tc0_over_Tn_r1=ratio)
with open(OUT, "w") as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
