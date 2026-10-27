"""P2: finite state capacity. Heavy-window inflation numbers and the capacity-type bounds.

Run:  python3 -I p2_capacity.py   (writes ../results/p2_capacity.json)
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cfs_common import (A_S, CAND1, G_GRUT, G_S0, G_SM, GAUGEPOINT_H, H0_DESI_CMB, H0_PLANCK, HEAVY,
                        K_PIVOT_MPC, M_P, M_PL, OMEGA_L_PLANCK, R_BK18, R_TRISTRAM, SCANPOINT_H, SEC_PER_YR,
                        T0, T_BBN_MIN, GEV_PER_INV_S, H0_GeV, H_lambda, H_of_r, H_of_V, N_of_k, S_dS,
                        T_of_rho, V_of_H, ln_aH_end_over_a0H0, r_of_H, rho_lambda, rho_of_T, s_of_T)

out = {}
L10 = math.log(10.0)

# ============================================================== 1. heavy-window inflation numbers
H = HEAVY["H_i"]
V = V_of_H(H)
heavy = {
    "H_i_GeV": H,
    "V_GeV4": V,
    "V_quarter_GeV": V ** 0.25,
    "r_if_H_i_is_H_at_pivot_exit": r_of_H(H),
    "epsilon_equiv_r_over_16": r_of_H(H) / 16.0,
    "T_instant_SM_GeV": T_of_rho(V, G_SM),
    "T_instant_GRUTdof_GeV": T_of_rho(V, G_GRUT),
    "H_over_M_P": H / M_P,
}
# minimal e-folds: today's Hubble scale inside the horizon at the start of inflation
cases = {}
for label, kw in {
    "instant_SM": dict(g_rh=G_SM),
    "instant_GRUTdof": dict(g_rh=G_GRUT),
    "TRH_2e14_matter_era": dict(T_rh=HEAVY["T_RH_lo"], w=0.0, g_rh=G_GRUT),
    "TRH_5e14_cap_at_instant": dict(T_rh=min(HEAVY["T_RH_hi"], T_of_rho(V, G_GRUT)), w=0.0, g_rh=G_GRUT),
}.items():
    n0 = ln_aH_end_over_a0H0(H, **kw)
    nstar = N_of_k(K_PIVOT_MPC, H, **kw)
    cases[label] = {"N_min_today_Hubble": n0, "N_star_k0.05": nstar,
                    "N_flat_Ok_1e-3_from_Oki_1": n0 + 0.5 * math.log(1e3),
                    "N_flat_Ok_1e-4_from_Oki_1": n0 + 0.5 * math.log(1e4)}
heavy["efolds"] = cases
heavy["N_min_central"] = cases["instant_GRUTdof"]["N_min_today_Hubble"]
heavy["N_star_central"] = cases["instant_GRUTdof"]["N_star_k0.05"]
# r against data (retrospective): BK18 r<0.036 -> H_max
heavy["H_max_from_BK18_r0.036"] = H_of_r(R_BK18)
heavy["H_max_from_r0.032"] = H_of_r(R_TRISTRAM)
heavy["r_gaugepoint_H4.7e13"] = r_of_H(GAUGEPOINT_H)
heavy["r_scanpoint_range"] = [r_of_H(SCANPOINT_H[0]), r_of_H(SCANPOINT_H[1])]
out["heavy_window"] = heavy

# ============================================================== 2. TCC (Bedroya-Vafa 1909.11063)
# TCC for one expanding phase ending at H_f:  a_f/a_i < M_P/H_f  (reduced Planck mass).
tcc = {}
N_tcc_heavy = math.log(M_P / H)
tcc["N_TCC_max_heavy"] = N_tcc_heavy
tcc["N_min_heavy"] = heavy["N_min_central"]
tcc["deficit_efolds_heavy"] = heavy["N_min_central"] - N_tcc_heavy


def tcc_margin_instant(Hx, g=G_SM):
    return math.log(M_P / Hx) - ln_aH_end_over_a0H0(Hx, g_rh=g)


H_tcc_std = brentq(lambda lx: tcc_margin_instant(10 ** lx), -10, 18)
H_tcc_std = 10 ** H_tcc_std
tcc["standard_instant_reheating"] = {
    "H_max_GeV": H_tcc_std,
    "V_quarter_max_GeV": V_of_H(H_tcc_std) ** 0.25,
    "r_max": r_of_H(H_tcc_std),
    "T_rh_max_GeV": T_of_rho(V_of_H(H_tcc_std)),
    "N_at_bound": math.log(M_P / H_tcc_std),
    "note": "BBLV 1909.11106 quote V^{1/4} < 6e8 GeV and r < 1e-30 [Q]; this is our own reproduction",
}
H_tcc_grut = 10 ** brentq(lambda lx: tcc_margin_instant(10 ** lx, G_GRUT), -10, 18)
tcc["standard_instant_reheating_GRUTdof"] = {"H_max_GeV": H_tcc_grut,
                                             "V_quarter_max_GeV": V_of_H(H_tcc_grut) ** 0.25}
# analytic form of the standard bound (derived in NOTES): V^{3/4} < 3 M_P^3 (H0/T0) (30/(pi^2 g))^{1/4} (g/g_s0)^{1/3}
Van = (3 * M_P ** 3 * (H0_GeV() / T0) * (30 / (math.pi ** 2 * G_SM)) ** 0.25 * (G_SM / G_S0) ** (1 / 3)) ** (4 / 3)
tcc["standard_analytic_V_quarter_GeV"] = Van ** 0.25

# most relaxed post-inflationary history: constant w in [-1/3, 1] until T_rh in [T_BBN, T_inst]


def nmin_relaxed(Hx):
    """minimum over (w, T_rh) of ln(a_end H_end/a0H0); deceleration requires w >= -1/3."""
    best = np.inf
    Tinst = T_of_rho(V_of_H(Hx))
    for w in np.linspace(-1 / 3 + 1e-9, 1.0, 61):
        for lT in np.linspace(math.log10(T_BBN_MIN), math.log10(Tinst), 81):
            try:
                n = ln_aH_end_over_a0H0(Hx, T_rh=10 ** lT, w=w, g_rh=G_SM if 10 ** lT > 200 else 10.75)
            except ValueError:
                continue
            best = min(best, n)
    return best


H_tcc_relaxed = 10 ** brentq(lambda lx: math.log(M_P / 10 ** lx) - nmin_relaxed(10 ** lx), 0, 18, xtol=1e-4)
tcc["most_relaxed_w-1/3_until_TBBN"] = {
    "H_max_GeV": H_tcc_relaxed,
    "V_quarter_max_GeV": V_of_H(H_tcc_relaxed) ** 0.25,
    "r_max": r_of_H(H_tcc_relaxed),
    "N_min_at_bound": nmin_relaxed(H_tcc_relaxed),
    "note": "minimised over w in [-1/3,1], T_rh in [5 MeV, T_inst]; optimum is w=-1/3, T_rh=T_BBN "
            "(Mizuno-Mukohyama-Pi-Zhang 1910.02979 give r <~ 1e-8 at the BBN reheating bound [Q via search])",
}
tcc["heavy_over_relaxed_H"] = H / H_tcc_relaxed
tcc["heavy_over_standard_H"] = H / H_tcc_std
tcc["heavy_over_standard_Vquarter"] = V ** 0.25 / V_of_H(H_tcc_std) ** 0.25

# heavy window with its own required reheating (T_RH 2-5e14, within ~ one Hubble time): strict line
tcc["heavy_own_reheating_N_min"] = [cases["TRH_2e14_matter_era"]["N_min_today_Hubble"],
                                    cases["TRH_5e14_cap_at_instant"]["N_min_today_Hubble"]]
# GRUT inflationary branches vs TCC
branches = {"heavy_spot": H, "scanpoint_lo": SCANPOINT_H[0], "scanpoint_hi": SCANPOINT_H[1],
            "gaugepoint": GAUGEPOINT_H}
tcc["branches"] = {k: {"H": v, "N_TCC": math.log(M_P / v), "N_min_instant": ln_aH_end_over_a0H0(v, g_rh=G_GRUT),
                       "H_over_H_TCC_std": v / H_tcc_std, "H_over_H_TCC_relaxed": v / H_tcc_relaxed}
                   for k, v in branches.items()}
# family floors inside the uniform scaling family (ls_scaling: yield ~ eps c_H^{3/2} H^{3/2}/lambda_X; c6 ~ 1/H^2)
c6_floor_H = H * math.sqrt(HEAVY["c6"] / HEAVY["c6_ceiling"])
yield_floor_H_eps1 = H * (HEAVY["eps"] / 1.0) ** (2 / 3)
tcc["family_floors"] = {
    "H_min_from_c6_survival_fixed_eps": c6_floor_H,
    "H_min_from_yield_with_eps_le_1": yield_floor_H_eps1,
    "required_eps_cH32_over_lambdaX_boost_to_reach_H_TCC_std": (H / H_tcc_std) ** 1.5,
    "note": "uniform scaling (M_Y,h,mu,H)->s(...) at fixed kappa, lambdas, e, c_H, delta [Q ls_scaling]",
}
# multistage rescue: stage 1 = heavy window with N1 <= N_TCC; stage 2 supplies the rest and dilutes the charge
N2_min = heavy["N_min_central"] - N_tcc_heavy
tcc["multistage"] = {"N2_min": N2_min, "log10_dilution_min_e^{3N2}": 3 * N2_min / L10,
                     "log10_H_boost_needed_via_Y~H^1.5": (3 * N2_min / L10) / 1.5}
# late-time TCC: lifetime of the present acceleration
HL = H_lambda()
t_tcc = math.log(M_P / HL) / HL * GEV_PER_INV_S / SEC_PER_YR
tcc["late_dS_lifetime_bound"] = {"ln_MP_over_HL": math.log(M_P / HL), "t_max_yr": t_tcc,
                                 "one_over_HL_yr": 1 / HL * GEV_PER_INV_S / SEC_PER_YR}
out["TCC"] = tcc

# ============================================================== 3. Banks-Fischler (astro-ph/0307459)
bf = {}
rhoL = rho_lambda()
SdS = S_dS(HL)
bf["rho_Lambda_GeV4"] = rhoL
bf["H_Lambda_GeV"] = HL
bf["S_dS_Lambda"] = SdS
bf["log10_dimH_eq_S_over_ln10"] = SdS / L10
# eq. (3.2) N <= (1/4) ln(Lambda_I/lambda) + (1/12) ln(M^4/lambda)  ["Planck units": M ambiguous]
for tag, M in (("reduced", M_P), ("nonreduced", M_PL)):
    bf[f"eq3.2_heavy_{tag}"] = 0.25 * math.log(V / rhoL) + math.log(M ** 4 / rhoL) / 12.0
bf["eq2.8_kappa1/3_heavy"] = 0.25 * math.log(V / rhoL)
bf["eq2.8_kappa1_heavy"] = math.log(V / rhoL) / 3.0
# reproduce BF's illustrative numbers: Lambda_I = (1e16 GeV)^4, lambda = (1e-3 eV)^4 = (1e-12 GeV)^4
LI, lam = 1e64, 1e-48
bf["BF_check_eq2.9_kappa1"] = math.log(LI / lam) / 3.0
bf["BF_check_eq2.8_kappa1/3"] = 0.25 * math.log(LI / lam)
bf["BF_check_eq3.2_nonreduced"] = 0.25 * math.log(LI / lam) + math.log(M_PL ** 4 / lam) / 12.0
bf["BF_check_eq3.2_reduced"] = 0.25 * math.log(LI / lam) + math.log(M_P ** 4 / lam) / 12.0
# own explicit-factor version: S_patch = S_H0 * exp(3 (N - N_min))  <=  S_dS  ->  Delta = (1/3) ln(S_dS/S_H0)
s0 = s_of_T(T0, G_S0)
S_H0 = s0 * 4 * math.pi / 3 / H0_GeV() ** 3
Delta = math.log(SdS / S_H0) / 3.0
bf["s0_GeV3"] = s0
bf["S_rad_in_today_Hubble_volume"] = S_H0
bf["Delta_Nmax_minus_Nmin"] = Delta
bf["own_N_max_heavy"] = heavy["N_min_central"] + Delta
bf["curvature_floor_factor_e^{-2Delta}"] = math.exp(-2 * Delta)
bf["topology_ceiling_factor_e^{Delta}"] = math.exp(Delta)
bf["heavy_window_compatible"] = bool(heavy["N_min_central"] <= bf["own_N_max_heavy"])
out["Banks_Fischler"] = bf

# ============================================================== 4. weaker capacity bounds (not binding)
out["weak_bounds"] = {
    "ADNTV_0704.1814_S_dS_inflation_heavy": S_dS(H),
    "DvaliGomezZell_1701.08776_MP2_over_H2_heavy": M_P ** 2 / H ** 2,
    "note": "N_e bounded by ~S_dS(H_inf) or ~M_P^2/H^2 [R]; ~1e14-1e16 e-folds at the heavy point: not binding",
}

# ============================================================== 5. Banks finite-N Lambda relation, GKS, Candidate 1
G = 1.0 / (8 * math.pi * M_P ** 2)
Lam_geom = 8 * math.pi * G * rhoL  # Lambda in GeV^2
out["finite_N_Lambda"] = {
    "Lambda_GeV2": Lam_geom,
    "S_from_3pi_over_G_Lambda": 3 * math.pi / (G * Lam_geom),
    "consistency_with_S_dS": 3 * math.pi / (G * Lam_geom) / SdS,
    "log10_recurrence_time_in_Hubble_units": SdS / L10,
    "note": "Lambda = 3 pi/(G ln N) relabels Lambda by N: N supplied (Banks hep-th/0007146, Fischler 2000) [R]",
}
v4 = CAND1["v"] ** 4
out["candidate1"] = {
    "V1loop_GeV4": CAND1["V1loop_over_v4"] * v4,
    "V1loop_over_rho_Lambda": CAND1["V1loop_over_v4"] * v4 / rhoL,
    "T_rh_max_under_TCC_std_GeV": tcc["standard_instant_reheating"]["T_rh_max_GeV"],
    "mX_GeV": CAND1["mX"],
    "TCC_reheating_room_over_mX": tcc["standard_instant_reheating"]["T_rh_max_GeV"] / CAND1["mX"],
}

# ============================================================== 6. sensitivity of the heavy-window numbers to H0
out["H0_sensitivity"] = {
    "N_min_heavy_H0_DESI_CMB": ln_aH_end_over_a0H0(H, g_rh=G_GRUT, h0_kms=H0_DESI_CMB),
    "N_min_heavy_H0_Planck": ln_aH_end_over_a0H0(H, g_rh=G_GRUT, h0_kms=H0_PLANCK),
}

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "p2_capacity.json")
with open(path, "w") as f:
    json.dump(out, f, indent=1)


def show(d, pre=""):
    for k, v in d.items():
        if isinstance(v, dict):
            show(v, pre + k + ".")
        else:
            print(f"{pre}{k} = {v:.6g}" if isinstance(v, float) else f"{pre}{k} = {v}")


show(out)
