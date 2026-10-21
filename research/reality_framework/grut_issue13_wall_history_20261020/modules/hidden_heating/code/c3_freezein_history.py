"""Item 2b/2c and the bath part of item 5: integrate gravitational freeze-in of hidden energy from the SM bath.

Inputs from c2 (computed here): beta_E (energy rate into the 4 hidden real dof, T^9/M_P^4), pair rate,
<T/E> of produced quanta.  Outside-action preparation data: T_RH, H_e (end of inflation), inflaton decay
width Gamma_I (sets T_RH), inflaton mass m_I (condensate scattering, item 3), direct hidden branching B_h
(item 4).  The SM bath is assumed to thermalise instantly with g* = 106.75.

(a) instantaneous reheating (analytic, checked by ODE)
(b) non-instantaneous reheating: inflaton matter domination with perturbative decay (Giudice-Kolb-Riotto
    hep-ph/0005123 scaling; Mambrini-Olive 2102.06214 eqs.(10)-(21) for T_max/T_RH), integrated numerically.
(c) formal hidden temperature if the hidden gas self-thermalised instantly (T_h = (30 rho_h/pi^2 g_h)^(1/4)),
    and the nonthermal induced curvature Pi_phi = sum_i c_i <chi_i^2> with <chi^2> = sum n <1/E>.
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq
from common import (M_P, V_BENCH, LAM_PHI, G_H, G_STAR_SM, C_AVG, c_H, rho_rad, T_of_rho, tc, dump,
                    TC_BAND)
import c2_graviton_rates

R2 = None


def rates():
    global R2
    if R2 is None:
        R2 = c2_graviton_rates.run()
    r = R2
    return dict(beta_E=r["beta_E_4dof"], pair_1dof=r["per_hidden_real_dof"]["pair_rate_coeff_T8_over_MP4"],
                TE=r["mc_mean_T_over_E_weighted"], gamma_chi_1dof=r["gamma_chi_per_dof_T7_over_MP4"])


# ------------------------------------------------------------------ (a) instantaneous
def inst_Th_max(TRH, beta_E):
    cH = c_H()
    K = 10.0 * beta_E / (np.pi**2 * G_H * cH * M_P**3)   # T_h^4 = K T^4 (TRH^3 - T^3)
    C7 = (3.0 / 7.0) * (4.0 / 7.0)**(4.0 / 3.0)
    return (K * C7)**0.25 * TRH**1.75


def inst_TRH_max(Tc, beta_E):
    cH = c_H()
    K = 10.0 * beta_E / (np.pi**2 * G_H * cH * M_P**3)
    C7 = (3.0 / 7.0) * (4.0 / 7.0)**(4.0 / 3.0)
    return (Tc / (K * C7)**0.25)**(4.0 / 7.0)


def inst_xi_inf(TRH, beta_E):
    cH = c_H()
    return (30.0 / (np.pi**2 * G_H) * beta_E * TRH**3 / (3 * cH * M_P**3))**0.25


def inst_Pi(T, TRH, gamma_chi_1dof):
    """Nonthermal induced curvature from the freeze-in gas (equal sharing over 4 dof)."""
    cH = c_H()
    chi2_1dof = gamma_chi_1dof * (TRH**3 - T**3) * T**2 / (3 * cH * M_P**3)
    return 4 * C_AVG * chi2_1dof


def inst_Pi_bounds(gamma_chi_1dof, v=V_BENCH):
    """Return T_RH thresholds for (A) max Pi > lambda v^2, (B) Pi - lambda v^2 > H^2 at some T."""
    lv2 = LAM_PHI * v**2
    cH = c_H()
    # (A): Pi_max = 4 C gamma 0.6 (0.4)^(2/3) TRH^5/(3 cH MP^3)
    coefA = 4 * C_AVG * gamma_chi_1dof * 0.6 * 0.4**(2.0 / 3.0) / (3 * cH * M_P**3)
    TRH_A = (lv2 / coefA)**0.2

    def margin(TRH):
        Ts = np.logspace(np.log10(TRH) - 12, np.log10(TRH), 4000)
        Pi = inst_Pi(Ts, TRH, gamma_chi_1dof)
        H2 = (cH * Ts**2 / M_P)**2
        return np.max((Pi - lv2) / H2)
    TRH_B = 10**brentq(lambda l: np.log(max(margin(10**l), 1e-300)), np.log10(TRH_A), 19.5)
    return TRH_A, TRH_B


# ------------------------------------------------------------------ (b) reheating ODE
def reheat(H_e, Gam, beta_E, pair_1dof=None, TE=None, m_I=None, cond_scat=False, B_h=0.0, higgs_cond=True,
           N_after=6.0, c_eff=None, npts=20000):
    """Background (inflaton + SM radiation) by ODE in log variables; hidden energy and <chi^2> by quadrature
    (the hidden sector does not backreact: rho_h << rho_tot is checked).  N = ln(a/a_e).
    Hidden sources: bath freeze-in beta_E T^9/M_P^4; condensate scattering 4 m_I rho_I^2/(128 pi M_P^4)
    (quadratic inflaton minimum; outside-action); direct branching B_h Gamma rho_I.
    chi2 = sum over the 4 dof of <chi^2>: bath 4*2*pair_1dof*TE*T^7/M_P^4; condensate 4 rho_I^2/(128 pi M_P^4 m_I);
    branching 4 B_h Gamma rho_I/m_I^2 (two quanta of E = m_I/2 per decay)."""
    if c_eff is None:
        c_eff = C_AVG
    rho_e = 3 * M_P**2 * H_e**2
    ag = np.pi**2 * G_STAR_SM / 30.0
    kc = 4 * m_I / (128 * np.pi * M_P**4) if (cond_scat and higgs_cond) else 0.0

    def back(N, y):
        lY, lX = y
        rI = rho_e * np.exp(lY - 3 * N)
        rR = rho_e * np.exp(lX - 4 * N)
        H = np.sqrt((rI + rR) / 3.0) / M_P
        sR = (1 - B_h) * Gam * rI + kc * rI**2
        return [-min(Gam / H, 1e4), np.exp(4 * N - lX) * sR / (H * rho_e)]

    N0 = 1e-6
    X0 = ((1 - B_h) * Gam / H_e + kc * rho_e / H_e) * N0
    N_RH_est = max((2.0 / 3.0) * np.log(max(H_e / Gam, 1.0)), 0.0)
    N_end = N_RH_est + N_after
    Ngrid = np.concatenate([[N0], np.linspace(1e-4, N_end, npts)])
    sol = solve_ivp(back, (N0, N_end), [0.0, np.log(X0)], method="LSODA", rtol=1e-9, atol=1e-11, t_eval=Ngrid)
    N = sol.t
    rI = rho_e * np.exp(sol.y[0] - 3 * N)
    rR = rho_e * np.exp(sol.y[1] - 4 * N)
    H = np.sqrt((rI + rR) / 3.0) / M_P
    T = (rR / ag)**0.25
    sh = beta_E * T**9 / M_P**4
    sc = np.zeros_like(T)
    if pair_1dof is not None:
        sc += 4 * 2 * pair_1dof * TE * T**7 / M_P**4
    if cond_scat:
        sh = sh + 4 * m_I * rI**2 / (128 * np.pi * M_P**4)
        sc = sc + 4 * rI**2 / (128 * np.pi * M_P**4 * m_I)
    if B_h > 0:
        sh = sh + B_h * Gam * rI
        sc = sc + 4 * B_h * Gam * rI / m_I**2
    Xh = cumulative_trapezoid(np.exp(4 * N) * sh / H, N, initial=0.0)
    Xc = cumulative_trapezoid(np.exp(2 * N) * sc / H, N, initial=0.0)
    rh = Xh * np.exp(-4 * N)
    rc = Xc * np.exp(-2 * N)
    i_eq = int(np.argmax(rR >= rI)) if np.any(rR >= rI) else len(rR) - 1
    Th_formal = T_of_rho(rh, G_H)
    Pi = c_eff * rc
    return dict(N=N, backreaction_max=float(np.max(rh / (rI + rR))), rho_I=rI, rho_R=rR, rho_h=rh, chi2=rc, H=H, T=T, Th_formal=Th_formal, Pi=Pi,
                T_max=T.max(), T_RH=T[i_eq], N_RH=sol.t[i_eq], Th_max=Th_formal.max(), ok=sol.success)


def Gamma_for_TRH(H_e, TRH):
    """Perturbative-decay width giving T_RH by the standard estimate Gamma = H(T_RH) (refined numerically later)."""
    return c_H() * TRH**2 / M_P


def run():
    rr = rates()
    bE = rr["beta_E"]
    out = {"inputs": {"beta_E_4dof": bE, "pair_rate_1dof": rr["pair_1dof"], "mean_T_over_E": rr["TE"],
                      "gamma_chi_1dof": rr["gamma_chi_1dof"], "g_h": G_H, "g_star_SM": G_STAR_SM,
                      "C_avg": C_AVG, "M_P": M_P}}
    # (a) instantaneous
    v = V_BENCH
    inst = {}
    for lab, f in [("Tc_nominal", 1.0), ("Tc_minus30", 1 - TC_BAND), ("Tc_plus30", 1 + TC_BAND)]:
        T_c = tc(v, f)
        TRHm = inst_TRH_max(T_c, bE)
        inst[lab] = {"T_c": T_c, "T_RH_max": TRHm, "xi_inf_at_T_RH_max": inst_xi_inf(TRHm, bE),
                     "T_h_max_check": inst_Th_max(TRHm, bE)}
    inst["scaling_law"] = "T_RH,max = 1.0e{:.4f} GeV (v/TeV)^(4/7)".format(np.log10(inst["Tc_nominal"]["T_RH_max"]))
    inst["T_RH_max_prefactor_GeV"] = inst["Tc_nominal"]["T_RH_max"]
    inst["T_h_max_formula"] = "T_h,max = [10 beta_E C7/(pi^2 g_h c_H M_P^3)]^(1/4) T_RH^(7/4), C7=(3/7)(4/7)^(4/3)"
    inst["T_h_max_coefficient"] = inst_Th_max(1.0, bE)
    # ODE check of instantaneous: pure radiation from T_RH, no inflaton
    TRH = inst["Tc_nominal"]["T_RH_max"]
    ag = np.pi**2 * G_STAR_SM / 30
    def rhs(N, y):
        rR, rh = np.exp(y)
        T = (rR / ag)**0.25
        H = np.sqrt((rR + rh) / 3) / M_P
        return [-4.0, -4.0 + bE * T**9 / (M_P**4 * H * rh)]
    sol = solve_ivp(rhs, (0, 8), [np.log(ag * TRH**4), np.log(ag * TRH**4 * 1e-60)], rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0, 8, 4001))
    Th = T_of_rho(np.exp(sol.y[1]), G_H)
    inst["ode_check_T_h_max"] = Th.max()
    inst["ode_check_rel_diff"] = abs(Th.max() - inst["Tc_nominal"]["T_c"]) / inst["Tc_nominal"]["T_c"]
    Ts = np.exp(-sol.t) * TRH
    inst["xi_of_T_rows_at_T_RH_max"] = [{"T_SM": float(Ts[i]), "T_h_formal": float(Th[i]), "xi": float(Th[i] / Ts[i])}
                                        for i in range(0, 4001, 400)]
    # nonthermal Pi bounds (instantaneous)
    TA, TB = inst_Pi_bounds(rr["gamma_chi_1dof"])
    inst["nonthermal_Pi_T_RH_threshold_noHubble"] = TA
    inst["nonthermal_Pi_T_RH_threshold_withHubble"] = TB
    inst["Pi_formula"] = "Pi(T) = 4 C_avg gamma_chi (T_RH^3 - T^3) T^2 / (3 c_H M_P^3), gamma_chi = 2 pair_1dof <T/E>"
    # thermalisation check at the boundary T_RH (dilute gas; 2->2 contact scattering, |M|^2 <= 36 lambda_phi^2)
    rho_h = (bE * TRH**3 / (3 * c_H() * M_P**3)) * TRH**4
    n_h = rho_h / (3.936 * TRH)
    E = 3.936 * TRH
    sig_v = 36 * LAM_PHI**2 / (32 * np.pi * 2 * E**2)
    Gam22 = n_h * sig_v
    inst["thermalisation_at_T_RH_max"] = {"rho_h": rho_h, "n_h": n_h, "Gamma_22": Gam22, "H": c_H() * TRH**2 / M_P,
                                          "Gamma_22_over_H": Gam22 / (c_H() * TRH**2 / M_P),
                                          "note": "Gamma/H scales as 1/T in RD afterwards; self-thermalisation fails"}
    out["instantaneous"] = inst
    # (b) non-instantaneous
    rows = []
    for H_e in [1e3, 1e5, 1e7, 1e9, 1e11, 1e13]:
        # find Gamma such that formal T_h,max (bath freeze-in only) = T_c; scan in Gamma
        def f(lg):
            r = reheat(H_e, 10**lg, bE, N_after=4.0, npts=6000)
            return np.log(r["Th_max"] / tc(v)), r
        lg_hi = np.log10(H_e) - 1e-3   # Gamma ~ H_e: near-instantaneous
        lg_lo = max(np.log10(H_e) - 40, np.log10(c_H() * 10.0**2 / M_P))   # T_RH >= ~10 GeV
        fhi, rhi = f(lg_hi)
        flo, rlo = f(lg_lo)
        if flo > 0:
            rows.append({"H_e": H_e, "binding": "all T_RH >= 10 GeV exceed", "T_h_max_at_slowest": rlo["Th_max"]})
            continue
        row = {"H_e": H_e, "T_RH_inst_limit": T_of_rho(3 * M_P**2 * H_e**2, G_STAR_SM)}
        if fhi < 0:
            row.update({"binding": False, "T_h_max_at_fastest_decay": rhi["Th_max"], "T_RH_at_fastest": rhi["T_RH"],
                        "T_max_at_fastest": rhi["T_max"]})
        else:
            lg = brentq(lambda x: f(x)[0], lg_lo, lg_hi, xtol=1e-4)
            r = f(lg)[1]
            # Pi criterion along the same history
            lv2 = LAM_PHI * v**2
            r2 = reheat(H_e, 10**lg, bE, pair_1dof=rr["pair_1dof"], TE=rr["TE"], N_after=4.0, npts=6000)
            row.update({"binding": True, "Gamma_I": 10**lg, "T_RH_max": r["T_RH"], "T_max_at_boundary": r["T_max"],
                        "T_max_over_T_RH": r["T_max"] / r["T_RH"], "T_h_max": r["Th_max"],
                        "Pi_max_over_lambda_v2_on_boundary": float(np.max(r2["Pi"]) / lv2),
                        "max_(Pi-lv2)/H2_on_boundary": float(np.max((r2["Pi"] - lv2) / r2["H"]**2))})
        rows.append(row)
    out["noninstantaneous_bath_only"] = rows
    # analytic cross-check of T_max/T_RH (Mambrini-Olive eq.21) for one case
    H_e, TRHt = 1e11, 1e8
    Gam = Gamma_for_TRH(H_e, TRHt)
    r = reheat(H_e, Gam, bE, N_after=3.0)
    rho_e = 3 * M_P**2 * H_e**2
    mo = (25 * np.sqrt(rho_e) / (16 * np.sqrt(3) * Gam * M_P) * (3 / 8)**0.6)**0.25
    out["MO_eq21_check"] = {"H_e": H_e, "Gamma_I": Gam, "T_max_num": r["T_max"], "T_RH_num(rho_R=rho_I)": r["T_RH"],
                            "T_max/T_RH_num": r["T_max"] / r["T_RH"], "MO_eq21_T_max/T_RH": mo,
                            "note": "MO define T_RH by rho_R(a_RH)=rho_I(a_RH) up to the bracket; O(10%) agreement expected"}
    # non-instantaneous: +-30 % T_c band and v scaling of the T_max boundary at fixed H_e (bath freeze-in only)
    def boundary(H_e, T_c):
        def f(lg):
            r = reheat(H_e, 10**lg, bE, N_after=4.0, npts=6000)
            return np.log(r["Th_max"] / T_c)
        lg_hi = np.log10(H_e) - 1e-3
        lg_lo = max(np.log10(H_e) - 40, np.log10(c_H() * 10.0**2 / M_P))
        if f(lg_hi) < 0 or f(lg_lo) > 0:
            return None
        lg = brentq(f, lg_lo, lg_hi, xtol=1e-4)
        r = reheat(H_e, 10**lg, bE, N_after=4.0, npts=6000)
        return {"T_RH_max": r["T_RH"], "T_max_max": r["T_max"]}
    band = {}
    for H_e in [1e9, 1e13]:
        band[f"H_e={H_e:.0e}"] = {lab: boundary(H_e, tc(v, f)) for lab, f in
                                  [("nominal", 1.0), ("minus30", 1 - TC_BAND), ("plus30", 1 + TC_BAND)]}
    out["noninstantaneous_Tc_band"] = band
    vs = {}
    for vv in [1e2, 1e3, 1e4, 1e5]:
        vs[f"{vv:.0e}"] = boundary(1e13, tc(vv))
    lv = np.log10([1e2, 1e5])
    lt = np.log10([vs["1e+02"]["T_max_max"], vs["1e+05"]["T_max_max"]])
    out["noninstantaneous_v_scaling_H_e_1e13"] = {"rows": vs, "slope_dlnTmax_dlnv": (lt[1] - lt[0]) / (lv[1] - lv[0]),
                                                  "expected_slope": 4.0 / 9.0}
    # condensate-scattering numeric check of the analytic peak (c4): bath negligible, slow decay
    H_e, mI = 1e9, 3e9
    rc = reheat(H_e, 1e-6, 0.0, m_I=mI, cond_scat=True, higgs_cond=True, N_after=2.0, npts=20000)
    pk = 0.0448 * 4 * mI * H_e**3 * 0.0433
    out["condensate_numeric_check"] = {"H_e": H_e, "m_I": mI, "rho_h_peak_numeric": float(rc["rho_h"].max()),
                                       "rho_h_peak_analytic": pk, "ratio": float(rc["rho_h"].max() / pk),
                                       "T_h_formal_peak_numeric": rc["Th_max"],
                                       "T_h_formal_peak_analytic": 0.277 * (mI * H_e**3)**0.25}
    # v scaling of the instantaneous bound
    out["v_scaling"] = {f"{vv:.0e}": inst_TRH_max(tc(vv), bE) for vv in [1e2, 1e3, 1e4, 1e5]}
    dump("c3_freezein_history.json", out)
    return out


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps({k: r[k] for k in ["inputs", "instantaneous", "MO_eq21_check", "v_scaling"]}, indent=1, default=float)[:6000])
    for row in r["noninstantaneous_bath_only"]:
        print(row)
