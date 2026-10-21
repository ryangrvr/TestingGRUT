"""Items 4 and 5: direct inflaton -> hidden branching B_h (OUTSIDE-ACTION preparation data) and nonthermal hidden
populations (induced curvature Pi_phi and thermalisation rates).

Pi_phi = sum_i c_i <chi_i^2>,  <chi_i^2> = int d^3p f_i(p) / ((2 pi)^3 E_p)  (thermal part of the coincident
two-point function of a real scalar; = T^2/12 for a massless Bose gas), c_phi = 3 lambda_phi, c_q1 = c_q2 = g_X^2,
c_D = h_D^2.  Restoration criterion used: Pi_phi > lambda_phi v^2 (requested), and the Hubble-damped variant
Pi_phi - lambda_phi v^2 > H^2 at the same time (the zero mode must respond).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import zeta
from common import (M_P, V_BENCH, LAM_PHI, G_X, H_D, G_H, G_STAR_SM, C_AVG, C_I, c_H, tc, T_of_rho, dump, TC_BAND)
import c3_freezein_history as c3

LV2 = LAM_PHI * V_BENCH**2
SIG_M2 = 36 * LAM_PHI**2      # largest contact |M|^2 among hidden 2->2 channels (phi phi -> phi phi)


def thermal_Bh_bound(T_max, f=1.0):
    T_c = tc(V_BENCH, f)
    x = (G_H / G_STAR_SM) * (T_c / T_max)**4
    return x / (1 + x)


def gamma22(n, E):
    s = 2 * E**2
    return n * SIG_M2 / (32 * np.pi * s)


def part_thermal():
    out = {}
    rows = []
    for Tm in [1e3, 3e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e10, 1e12, 1e14, 1e16]:
        row = {"T_max_visible": Tm}
        for lab, f in [("nominal", 1.0), ("minus30", 1 - TC_BAND), ("plus30", 1 + TC_BAND)]:
            row[f"B_h_max_{lab}"] = thermal_Bh_bound(Tm, f) if Tm > tc(V_BENCH, f) * (G_H / G_STAR_SM)**0.25 else 1.0
        row["xi_at_boundary"] = min(tc(V_BENCH) / Tm, 1.0)
        rows.append(row)
    out["rows"] = rows
    out["formula"] = "T_h,max = (B_h g*/((1-B_h) g_h))^(1/4) T_max; B_h,max = (g_h/g*)(T_c/T_max)^4 for B_h<<1"
    # numerical check with the reheating ODE (no gravitational sources): rho_h/rho_R -> B_h/(1-B_h)
    rr = c3.rates()
    B = 1e-12
    r = c3.reheat(1e9, 1e-2, 0.0, m_I=3e9, B_h=B, N_after=3.0, npts=4000)
    ratio = r["rho_h"][-1] / r["rho_R"][-1]
    out["ode_check"] = {"B_h": B, "rho_h/rho_R_end": ratio, "expected": B / (1 - B),
                        "T_h,max/T_max_num": r["Th_max"] / r["T_max"],
                        "expected_(B g*/g_h)^1/4": (B * G_STAR_SM / G_H)**0.25}
    # hidden thermal equilibrium requirement at T_h = T_c given SM at T_SM (estimate)
    T_c = tc(V_BENCH)
    n_eq = G_H * zeta(3) * T_c**3 / np.pi**2
    s_mean = 2 * (2.701 * T_c)**2
    Gam = n_eq * SIG_M2 / (32 * np.pi * s_mean)
    T_SM_max = np.sqrt(Gam * M_P / c_H())
    out["kinetic_equilibrium_at_Tc"] = {"Gamma_22_at_Tc": Gam, "T_SM_max_for_Gamma>H": T_SM_max,
                                        "xi_min": T_c / T_SM_max,
                                        "thermal_potential_module_xi_min_kinetic": 6.6e-5,
                                        "thermal_potential_module_xi_min_number_changing": 0.068,
                                        "consequence": "thermal B_h bound is self-consistent only if T_max <~ T_c/xi_min"}
    out["T_max_selfconsistency_kinetic"] = T_c / 6.6e-5
    out["T_max_selfconsistency_number_changing"] = T_c / 0.068
    return out


def part_nonthermal():
    """Direct branching into hidden quanta of E = m_I/2 that do not thermalise."""
    out = {}
    rows = []
    for H_e, T_RH_target, mH in [(1e5, 1e8, 3.0), (1e9, 1e8, 3.0), (1e13, 1e8, 3.0), (1e13, 1e12, 3.0), (1e9, 1e6, 1.0)]:
        Gam = c_H() * T_RH_target**2 / M_P
        m_I = mH * H_e
        r1 = c3.reheat(H_e, Gam, 0.0, m_I=m_I, B_h=1.0e-30, N_after=3.0, npts=8000, c_eff=1.0)
        chi1 = r1["chi2"] / 1.0e-30      # chi2 per unit B_h (B_h << 1 so background unchanged)
        H2 = r1["H"]**2
        res = {"H_e": H_e, "T_RH": r1["T_RH"], "T_max": r1["T_max"], "m_I": m_I}
        for lab, ce in [("equal_share", C_AVG), ("phi_only", C_I["phi"]), ("D_only", C_I["D"])]:
            Bnh = LV2 / (ce * chi1.max())
            def g(lB):
                return np.log(max(np.max((10**lB * ce * chi1 - LV2) / H2), 1e-300))
            try:
                Bh = 10**brentq(g, np.log10(Bnh) - 1e-9, 40.0)
            except ValueError:
                Bh = None
            res[f"B_h_max_noHubble_{lab}"] = Bnh
            res[f"B_h_max_withHubble_{lab}"] = Bh
        res["B_h_max_thermal_same_T_max"] = thermal_Bh_bound(r1["T_max"])
        # thermalisation of the branching gas at the no-Hubble boundary (equal share), at a = a_e
        Bn = res["B_h_max_noHubble_equal_share"]
        n0 = 2 * Bn * Gam * 3 * M_P**2 * H_e**2 / m_I / H_e      # ~ produced in one Hubble time at a_e
        res["Gamma22_over_H_at_a_e_noHubble_boundary"] = gamma22(n0, m_I / 2) / H_e
        rows.append(res)
    out["rows"] = rows
    out["note"] = ("B_h needed for nonthermal restoration is many orders larger than the thermal bound; the thermal "
                   "bound needs hidden self-thermalisation, which itself needs xi >~ 6.6e-5")
    return out


def part_freezein_gas():
    rr = c3.rates()
    out = {"Pi_formula": "Pi = 4 C_avg gamma_chi (T_RH^3 - T^3) T^2/(3 c_H M_P^3)", "gamma_chi_1dof": rr["gamma_chi_1dof"],
           "mean_T_over_E": rr["TE"], "C_avg": C_AVG}
    TA, TB = c3.inst_Pi_bounds(rr["gamma_chi_1dof"])
    out["T_RH_threshold_noHubble"] = TA
    out["T_RH_threshold_withHubble"] = TB
    out["required_H_inf_for_instant_TRH_noHubble"] = c_H() * TA**2 / M_P
    return out


def part_condensate_gas():
    out = {}
    rows = []
    for H_e in [1e3, 1e5, 1e7, 1e9, 1e11, 1e13]:
        for mH in [1.0, 3.0]:
            m = mH * H_e
            rho_pk = 7.76e-3 * m * H_e**3
            n = rho_pk / m
            Gam = gamma22(n, m)
            H_pk = H_e * 1.27**-1.5
            rows.append({"H_e": H_e, "m_I/H_e": mH, "Gamma22_over_H_at_peak": Gam / H_pk,
                         "a_over_a_e_for_Gamma=H_MD": (H_pk / Gam)**2 * 1.27 if Gam > 0 else None})
    out["rows"] = rows
    out["scaling"] = "Gamma/H ∝ a^(1/2) in MD, ∝ a in RD after production"
    return out


def run():
    out = {"thermal_branching": part_thermal(), "nonthermal_branching": part_nonthermal(),
           "freezein_gas": part_freezein_gas(), "condensate_gas_thermalisation": part_condensate_gas(),
           "Pi_definition": "Pi_phi = 3 lambda_phi <dphi^2> + g_X^2 (<q1^2>+<q2^2>) + h_D^2 <D^2>, <chi^2> = int d^3p f/((2pi)^3 E)",
           "thermal_check": {"Pi/T^2 for thermal 4-dof bath": (3 * LAM_PHI + 2 * G_X**2 + H_D**2) / 12.0,
                             "brief_value": 0.004175}}
    dump("c5_branching_nonthermal.json", out)
    return out


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps(r, indent=1, default=float)[:7000])
