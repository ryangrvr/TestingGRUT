"""Part 5: map the charged phase structure onto a cosmological charge yield Y_Q = r Y_B.

Inputs (quoted, BRIEF.md): omega_b = 0.02233, s_0 = 2891.2 cm^-3, rho_crit/h^2 = 1.05375e-5 GeV cm^-3,
m_p = 0.938272 GeV, g_*s = 106.75 (SM, T > 200 GeV), reduced Planck mass 2.435323e18 GeV.
A HOMOGENEOUS hidden charge density n_Q(T) = Y_Q s_SM(T), s_SM = (2 pi^2/45) g_*s T^3, T = visible temperature.
Hidden temperature T_h = xi * T is an OUTSIDE-ACTION preparation input (xi constant is an idealisation).
Run after t0_phase.py and thermal_phase.py:  python3 cosmo_map.py
"""
import json
import os
import numpy as np
from common import V_GEV, LAM_PHI, RESULTS_DIR, save_json

OMEGA_B = 0.02233
S0 = 2891.2
RHOC_H2 = 1.05375e-5
M_P_GEV = 0.938272
GSTAR = 106.75
MPL_RED = 2.435323e18

n_b = OMEGA_B * RHOC_H2 / M_P_GEV          # cm^-3
Y_B = n_b / S0
S_COEF = 2 * np.pi**2 / 45 * GSTAR          # s = S_COEF T^3

t0 = json.load(open(os.path.join(RESULTS_DIR, "t0_phase.json")))
th = json.load(open(os.path.join(RESULTS_DIR, "thermal_phase.json")))

thresholds = {
    "n_c (phi=0 locally stable; T_h=0)": t0["n_c_over_v3"],
    "n_t (phi=0 global among homogeneous; T_h=0)": t0["homogeneous_global_switch_n_over_v3"],
    "n_s1 (phi!=0 branch disappears; T_h=0)": t0["branch1_spinodal_n_over_v3"],
    "n_star (binodal: homogeneous phi=0 is the ground state; T_h=0)": t0["binodal_n_star_over_v3"],
    "max_Th n_local(T_h) for T_h<Tc0 (leading order)": float(max(th["local_line_dense"]["n_local_over_v3"])),
}


def T_of_n(n_over_v3, r, v=V_GEV):
    """visible T (GeV) at which r Y_B s_SM(T) = n v^3."""
    return (n_over_v3 / (r * Y_B * S_COEF)) ** (1 / 3) * v


def hubble(T):
    return np.pi * np.sqrt(GSTAR / 90) * T**2 / MPL_RED


rs = [1.0, 0.64, 0.1, 0.01]
table = []
for name, nthr in thresholds.items():
    for r in rs:
        T = T_of_n(nthr, r)
        table.append(dict(threshold=name, n_thr_over_v3=nthr, r=r, T_vis_GeV=T, T_vis_over_v=T / V_GEV,
                          H_GeV=hubble(T), H_over_mX=hubble(T) / (0.1 * V_GEV)))
scaling = ("T_n(r, v) = (n_thr/(r Y_B s_coef))^{1/3} v = %.4g r^{-1/3} v  for n_thr = n_c (fixed couplings)"
           % ((t0["n_c_over_v3"] / (Y_B * S_COEF)) ** (1 / 3)))

# Kibble-Zurek style freeze-out length for a homogeneous phi=0 phase diluted through n_c by expansion
# (mean-field exponents nu=1/2, z=1: xi ~ m^-1 (m tau_Q)^{1/3}; tau_Q = 1/(3H) dilution time).  ESTIMATE.
kz = []
for r in rs:
    T = T_of_n(t0["n_c_over_v3"], r)
    H = hubble(T)
    m = np.sqrt(LAM_PHI) * V_GEV
    tauQ = 1 / (3 * H)
    xi_kz = (m * tauQ) ** (1 / 3) / m
    kz.append(dict(r=r, T_vis_GeV=T, m_unit_GeV=m, tauQ_GeV_inv=tauQ, xi_KZ_GeV_inv=xi_kz,
                   xi_KZ_in_ell=xi_kz * m, hubble_radius_in_ell=m / H,
                   domains_per_hubble_volume=(1 / (H * xi_kz)) ** 3))

# trajectories with a hidden/visible temperature ratio xi
Tl = np.array(th["local_line_dense"]["T_over_v"])
nl = np.array(th["local_line_dense"]["n_local_over_v3"])
Tc0 = th["Tc_zero_charge_over_v"]


def n_local(Th):
    if Th >= Tc0:
        return 0.0
    return float(np.interp(Th, Tl, nl))


traj = []
for r in rs:
    for xi in [0.0, 1e-3, 3e-3, 5e-3, 7.5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.3, 1.0]:
        Tgrid = np.geomspace(1e5 * V_GEV, 1e-2 * V_GEV, 200001) / V_GEV     # T_vis in units of v, decreasing
        n = r * Y_B * S_COEF * Tgrid**3
        Th = xi * Tgrid
        sym = np.array([(t >= Tc0) or (nn >= n_local(t)) for t, nn in zip(Th, n)]) if xi > 0 else n >= nl[0]
        if not sym[0]:
            traj.append(dict(r=r, xi=xi, note="never symmetric on grid"))
            continue
        k = int(np.argmax(~sym))
        Tex = Tgrid[k]
        route = "charge-supported (T_h,exit < T_cp~1.54 v)" if xi * Tex < 1.54 else "thermal (T_h,exit ~ Tc0)"
        traj.append(dict(r=r, xi=xi, T_vis_exit_GeV=Tex * V_GEV, T_h_exit_GeV=xi * Tex * V_GEV,
                         n_exit_over_v3=r * Y_B * S_COEF * Tex**3, route=route))

out = dict(
    status="derived/computed; Y_B quoted inputs; xi and the homogeneity of the charge are preparation inputs",
    n_b_cm3=n_b, Y_B=Y_B, s_coef=S_COEF, thresholds_over_v3=thresholds, table=table, scaling=scaling,
    kibble_zurek_estimate=kz, xi_trajectories=traj, Tc0_over_v=Tc0,
    codex_production_check=dict(n_over_v3=0.022008521077073767,
                                T_vis_GeV_r1=T_of_n(0.022008521077073767, 1.0), codex_reported_T_vis_GeV=1.7563e5),
)
save_json("cosmo_map.json", out)
if __name__ == "__main__":
    print("Y_B", Y_B, "n_b", n_b)
    for row in table:
        print(row["threshold"][:40], row["r"], "%.4g GeV" % row["T_vis_GeV"], "H/mX %.3g" % row["H_over_mX"])
    print(scaling)
    for k in kz:
        print(k)
    for t in traj:
        print(t)
    print(out["codex_production_check"])
