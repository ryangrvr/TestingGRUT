#!/usr/bin/env python3
"""history.py -- assemble the explicit hidden histories of the issue #13 wall calculation.

Recomputes, from the unchanged action and stated preparation inputs, the numbers that the report quotes for
the reference histories H-star (hot hidden sector), H-star-Q (charged variant) and C-star (cold hidden sector),
and places them next to the values computed by the modules (read from modules/*/results). Every input that the
action does not supply is marked "P" (preparation, outside the action).

Usage: python3 code/history.py            -> writes results/history.json and prints a summary
Requires numpy and scipy. Runs in about one second.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

PKT = Path(__file__).resolve().parent.parent
MOD = PKT / "modules"

# ---------------------------------------------------------------- action (unchanged benchmark)
v, gX, lphi, lX, hD, lD = 1000.0, 0.1, 0.01, 0.001, 0.01, 0.001
mX, mphi, mD = gX * v, math.sqrt(2 * lphi) * v, math.sqrt((gX * v) ** 2 + (hD * v) ** 2)
MP = 2.435323e18                       # reduced Planck mass, GeV
sigma = 2 * math.sqrt(2 * lphi) / 3 * v ** 3
mstar = v * min(gX, (lphi * lX) ** 0.25)

# ---------------------------------------------------------------- constants (calibration data, retrospective)
T0 = 2.7255 * 8.617333262e-14         # GeV
s0_cm3, rhoc_h2 = 2891.2, 1.05375e-5   # cm^-3, GeV cm^-3
m_p = 0.938272
omega_b = 0.02233                      # Mossa et al. 2020 (standard BBN)
Y_B = omega_b * rhoc_h2 / (m_p * s0_cm3)
H100 = 100 / 3.0856775814913673e19 * 6.582119569e-25  # GeV
gs_hi = 106.75
g_rho_late = 2 + 7 / 8 * 2 * 3.044 * (4 / 11) ** (4 / 3)
OMEGA_DM_H2 = 0.12                     # Planck 2018, retrospective


def s_SM(T, g=gs_hi):
    return 2 * math.pi ** 2 / 45 * g * T ** 3


def H_rad(T, g=gs_hi):
    return math.sqrt(math.pi ** 2 * g / 90) * T ** 2 / MP


def load(rel):
    p = MOD / rel
    return json.loads(p.read_text()) if p.exists() else None


# ---------------------------------------------------------------- 1. zero-charge transition (leading check)
Pi_phi = lphi / 4 + gX ** 2 / 6 + hD ** 2 / 12
Tc_lead = math.sqrt(lphi / Pi_phi) * v
tp = load("thermal_potential/results/summary.json") or {}
Tc_range = tp.get("Tc_adopted_range_over_v", [1.557140, 1.564147])

# ---------------------------------------------------------------- 2. charged phase structure at T = 0 (D = 0)
def rho_sym(n):          # phi = 0 branch, |X|^2 = (n^2/(8 lX))^(1/3)
    X2 = (n * n / (8 * lX)) ** (1 / 3)
    return n * n / (4 * X2) + lphi / 4 * v ** 4 + lX * X2 * X2


def rho_brk(n):          # phi != 0 branch, phi^2 = v^2 - 2 gX^2 X2/lphi
    def f(X2):
        p2 = max(v * v - 2 * gX * gX * X2 / lphi, 0.0)
        return n * n / (4 * X2) + lphi / 4 * (p2 - v * v) ** 2 + gX * gX * p2 * X2 + lX * X2 * X2
    r = minimize_scalar(f, bounds=(1e-12 * v * v, lphi * v * v / (2 * gX * gX)), method="bounded",
                        options=dict(xatol=1e-14 * v * v))
    return r.fun


n_c = math.sqrt(8 * lX) * (lphi * v * v / (2 * gX * gX)) ** 1.5
n_t = brentq(lambda n: rho_sym(n) - rho_brk(n), 0.035 * v ** 3, 0.05 * v ** 3, xtol=1e-12 * v ** 3)
# phi != 0 branch at fixed n: rho = n^2/(4 X2) + gX^2 v^2 X2 - c X2^2, c = gX^4/lphi - lX.
# Stationary: n^2 = 4 X2^2 (gX^2 v^2 - 2 c X2); spinodal: n^2 = 4 c X2^3  =>  X2 = gX^2 v^2 / (3 c).
c_ = gX ** 4 / lphi - lX
X2_s = gX ** 2 * v ** 2 / (3 * c_)
n_s1 = math.sqrt(4 * c_ * X2_s ** 3)
ns = np.linspace(0.02, 1.0, 200001) * v ** 3
rn = np.array([rho_sym(n) for n in ns]) / ns
n_star, mu_star = ns[np.argmin(rn)], rn.min()


def T_supply(n_thr, r):   # visible temperature at which a homogeneous Y_Q = r Y_B reaches n_thr
    return (n_thr / (r * Y_B * 2 * math.pi ** 2 / 45 * gs_hi)) ** (1 / 3)


# ---------------------------------------------------------------- 3. wall-network cosmology (simple cross-check)
def T_parity(A):        # rho_w = A sigma/t = 2 A sigma H (radiation), rho_r = 3 MP^2 H^2, late g_rho
    H = 2 * A * sigma / (3 * MP ** 2)
    rho = 3 * MP ** 2 * H ** 2
    return (rho / (math.pi ** 2 / 30 * g_rho_late)) ** 0.25


def H0_after_parity(A, n):   # rho_w ∝ a^-n after parity, photons T ∝ 1/a, matter neglected (cross-check only)
    Td = T_parity(A)
    rho_d = math.pi ** 2 / 30 * g_rho_late * Td ** 4
    return math.sqrt(rho_d * (T0 / Td) ** n / (3 * MP ** 2)) / H100 * 100   # km/s/Mpc


rho_nu1 = 7 / 8 * 2 * (4 / 11) ** (4 / 3)
dNeff_conservative = g_rho_late / rho_nu1   # rho_hidden >= rho_rad at recombination (w <= 1/3 after parity)

# ---------------------------------------------------------------- 4. relic bookkeeping
Ymax_strict = OMEGA_DM_H2 * rhoc_h2 / (s0_cm3 * mX)


def Y_thermal(xi, g_h=4, gs=gs_hi):     # mu_N = 0 Bose number per visible entropy
    return (1.2020569 / math.pi ** 2) * g_h * xi ** 3 / (2 * math.pi ** 2 / 45 * gs)


def omega_from_Y(Y, m):
    return Y * m * s0_cm3 / rhoc_h2


ry = load("relic_yield/results/headline.json") or {}
floor_coef = ry.get("headline", {}).get("floor", {}).get("coef_GeV_per_omega_sym", 79.28029870506704)

# ---------------------------------------------------------------- 5. reference histories
out = {"inputs": {"v": v, "g_X": gX, "lambda_phi": lphi, "lambda_X": lX, "h_D": hD, "lambda_D": lD, "m_D0": gX * v,
                  "M_P_reduced": MP, "sigma_GeV3": sigma, "m_star": mstar, "Y_B": Y_B},
       "phase": {"Tc_leading_over_v": Tc_lead / v, "Tc_resummed_range_over_v": Tc_range,
                 "n_c_over_v3": n_c / v ** 3, "n_t_over_v3": n_t / v ** 3, "n_s1_over_v3": n_s1 / v ** 3,
                 "n_star_over_v3": n_star / v ** 3, "Qmatter_E_per_Q_GeV": mu_star,
                 "T_n_r1_TeV": T_supply(n_c, 1) / 1e3, "T_s1_r1_TeV": T_supply(n_s1, 1) / 1e3},
       "network": {"T_parity_keV_A": {str(A): T_parity(A) * 1e6 for A in (0.7, 0.8, 0.9)},
                   "H0_pred_n1_kmsMpc_A0.8_simple": H0_after_parity(0.8, 1),
                   "H0_pred_n2_kmsMpc_A0.8_simple": H0_after_parity(0.8, 2),
                   "DeltaNeff_conservative_min": dNeff_conservative},
       "relic": {"Y_N_max_strict": Ymax_strict, "floor_coef_GeV_per_omega_sym": floor_coef}}

# H-star: hot, relic-compatible hidden sector, negligible charge
xi = 5e-4
Tc_mid = 0.5 * (Tc_range[0] + Tc_range[1]) * v
T_tr = Tc_mid / xi
YN = Y_thermal(xi)
om_sym = omega_from_Y(YN, 101.7)   # thermal-reference mean mass 101.7 GeV (relic_yield)
r_H = 1e-4
Hstar = {
    "inputs": {"H_inf_GeV [P]": 1e8, "reheating [P]": "instant, visible T_RH ~ 8e12 GeV",
               "hidden bath [P]": "P-invariant, all 4 dof, kinetically thermal, mu_N = 0", "xi [P]": xi,
               "r = Y_Q/Y_B [P]": r_H, "rho_Lambda [matching, free]": 0.0},
    "T_SM_at_transition_GeV": T_tr, "H_over_mX_at_transition": H_rad(T_tr) / mX,
    "Hinv_over_wall_width": 1 / H_rad(T_tr) / (math.sqrt(2) / (math.sqrt(lphi) * v)),
    "Y_N_mu0": YN, "omega_sym_estimate": om_sym, "omega_asym": r_H * omega_b * mX / m_p,
    "maser_product_GeV": mX * r_H + floor_coef * om_sym,
    "network": "forms (thermal exit, p = 1/2) and persists (neutral lattice + scaling literature)",
    "T_parity_keV": T_parity(0.8) * 1e6, "DeltaNeff_min_at_recombination": dNeff_conservative,
}
# C-star: cold hidden sector, low-scale inflation
H_c = 120.0
T_RH_inst = (90 / (math.pi ** 2 * gs_hi)) ** 0.25 * math.sqrt(H_c * MP)
r_C, T_Q = 0.05, 1e5
n_TQ = r_C * Y_B * s_SM(T_Q)
pm = load("production_map/results/headline.json") or {}
cB = pm.get("B", {}).get("Y_N floor/J coefficients (instant), Y = c (H/M_P)^1.5", [0.184, 2.83])
Cstar = {
    "inputs": {"H_inf_GeV [P]": H_c, "N_* [P]": "~48 (instant reheating)", "start [P]": "+v at start of last N_* e-folds",
               "reheating [P]": "instant, visible only, B_h = 0", "r [P]": r_C, "T_Q_GeV [P]": T_Q,
               "charge form [P]": "homogeneous condensate on the +v branch"},
    "T_RH_instant_GeV": T_RH_inst,
    "below regime-A boundary 145.6 GeV (minority_relic)": H_c < 145.6,
    "gravity-only heating below nonthermal restoration threshold 9.07e12 GeV": T_RH_inst < 9.07e12,
    "spectator Y_N range": [c * (H_c / MP) ** 1.5 for c in cB],
    "n_Q(T_Q)/n_c": n_TQ / n_c,
    "maser_product_GeV": mX * r_C,
    "network": "none (no restoration, no opposite-vacuum region expected, charge sub-critical on +v)",
}
# fate of a formed network is set by the systematic exit bias m0 (bias_threshold, residual_bias)
bt = load("bias_threshold/results/thresholds.json") or {}
rb = load("residual_bias/results/headline.json") or {}
hs_bt = bt.get("histories", {}).get("H_star", {}).get("central", {})
eps_p, eps_s = hs_bt.get("eps_persist"), hs_bt.get("eps_safe")
rA = rb.get("routes", {}).get("A_Hstar_thermal [P: xi=5e-4, bath from T_SM=1e8 GeV, phi_0 = v]", {})
log_m0 = rA.get("log10_abs_m0_central")


def fate(log10_m0):
    if None in (log10_m0, eps_p, eps_s):
        return None
    if log10_m0 < math.log10(eps_p):
        return "persists to parity (rejected by parity / DeltaN_eff / R_H)"
    if log10_m0 < math.log10(eps_s):
        return "decays after T_safe (stable-relic overproduction)"
    return "decays early (not rejected by the wall tests computed here)"


Hstar["bias"] = {"eps_persist": eps_p, "eps_safe": eps_s, "log10_m0_thermal_erasure": log_m0,
                 "fate": fate(log_m0)}
out["histories"] = {"H_star": Hstar, "C_star": Cstar}
out["bias_routes"] = {
    "thermal bath (H-star)": fate(log_m0),
    "Hubble friction only (no bath)": fate(0.0),
    "inflation trapping, H_inf = 5 TeV, fresh start (|m0| = 2.18e-3)": fate(math.log10(2.18e-3)),
    "charge forced to 2 n_c (|m0| ~ 0.5)": fate(math.log10(0.5)),
    "exact Z2 ensemble, no prior (m0 = 0)": "persists to parity (rejected by parity / DeltaN_eff / R_H)",
}

# ---------------------------------------------------------------- 6. module values placed side by side
co = load("cosmo_observable/results/headline.json") or {}
ce = load("charged_exit_sim/results/headline.json") or {}
ns_ = load("network_sim/results/headline.json") or {}
mr = load("minority_relic/results/headline.json") or {}
out["module_values"] = {
    "T_dom_keV_backreaction (cosmo_observable)": co.get("T_dom_keV (rho_w=rho_rad, A=0.8, backreaction) [C]"),
    "DeltaNeff stress n=4 (cosmo_observable)": co.get("stress n=4 (conversion to radiation at T_dom): DeltaNeff at T_rec [C]"),
    "A_3D_PRS_256 (network_sim)": ns_.get("A_3D_PRS_256_box_safe"),
    "kappa_E (minority_relic)": mr.get("collapse_energy", {}).get("kappa_E_RD"),
    "charged exit p(+) (charged_exit_sim)": ce.get("sign_statistics", {}).get("charged_all"),
}

(PKT / "results").mkdir(exist_ok=True)
(PKT / "results" / "history.json").write_text(json.dumps(out, indent=1, default=str) + "\n")

if __name__ == "__main__":
    ph = out["phase"]
    print(f"sigma = {sigma:.6e} GeV^3, m* = {mstar:.6f} GeV, Y_B = {Y_B:.4e}")
    print(f"T_c leading {Tc_lead / v:.6f} v; resummed {Tc_range}")
    print(f"n_c {ph['n_c_over_v3']:.6f}, n_t {ph['n_t_over_v3']:.6f}, n_s1 {ph['n_s1_over_v3']:.6f}, "
          f"n_* {ph['n_star_over_v3']:.4f} v^3; Q-matter E/Q {mu_star:.3f} GeV")
    print(f"T_n(r=1) {ph['T_n_r1_TeV']:.1f} TeV, T_s1(r=1) {ph['T_s1_r1_TeV']:.1f} TeV")
    print("T_parity keV", out["network"]["T_parity_keV_A"], " DeltaNeff_min", round(dNeff_conservative, 3))
    print("simple H0_pred n=1, n=2 (km/s/Mpc):", f"{H0_after_parity(0.8, 1):.3e}", f"{H0_after_parity(0.8, 2):.3e}")
    print("H*:", {k: Hstar[k] for k in ("T_SM_at_transition_GeV", "H_over_mX_at_transition", "Hinv_over_wall_width",
                                         "Y_N_mu0", "omega_sym_estimate", "maser_product_GeV")})
    print("C*:", {k: Cstar[k] for k in ("T_RH_instant_GeV", "spectator Y_N range", "n_Q(T_Q)/n_c", "maser_product_GeV")})
    print("module values:", out["module_values"])
