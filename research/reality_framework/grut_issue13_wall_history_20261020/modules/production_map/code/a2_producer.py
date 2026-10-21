"""production_map / a2: heavy-inflaton direct branching into the hidden quartet -> (Y_N, eps, xi, mu_N), Pi(t), B_h window.

Output: results/a2_producer.json.  Deterministic, < 1 min.  Needs results/a1_restoration.json and
results/b1_spectator.json (run a1, b1 first).

Producer [P, outside the action]: inflaton of mass m_I, decay rate Gamma_I = H(T_RH) = H_COEF T_RH^2/M_P,
branching B_h into hidden pairs (phi phi, X Xbar, D D with equal sharing among the 4 real dof), each quantum injected
with E = m_I/2.  End-of-inflation Hubble H_e is a further free input.  The hidden quanta free-stream (no portal; their
own 2->2 rate is computed and is far too slow), so the curvature they induce at phi = 0 is the free-streaming
Pi = c_eff J, J = int d^3k f/E (summed over dof) [D].

Comoving bookkeeping in N = ln a (a_e = 1), per unit B_h (B_h << 1, no back-reaction):
  P = rho_I a^3/(He^2 MP^2):      dP/dN = -(gam/h) P
  R = rho_r a^4/(He^2 MP^2):      dR/dN = (gam/h) P a
  Nh = n_h a^3 m_I/(He^2 MP^2 B): dNh/dN = 2 (gam/h) P
  Eh = rho_h a^4/(He^2 MP^2 B):   dEh/dN = (gam/h) P a
  Jh = J a^2 m_I^2/(He^2 MP^2 B): dJh/dN = 4 (gam/h) P / a
  h^2 = (P/a^3 + R/a^4)/3,  gam = Gamma_I/H_e.
"""
import json
import math

import numpy as np
from scipy.integrate import solve_ivp

from common import (v, lphi, MP, gstar, gh, S_COEF, R_COEF, H_COEF, Y_MAX_STRICT, C_EFF_EQUAL, COUP_PHI,
                    gamma_el, H_of_T, T_of_H, dump, load)

OUT = {"module": "production_map", "part": "A2 heavy-inflaton producer map and B_h window",
       "labels": "[D] derived, [C] computed, [Q] quoted, [E] estimate, [P] preparation datum"}
b1 = load("b1_spectator.json")
inst = b1["instant reheating, minimal couplings [C]"]
YF = inst["Y_floor = coef (H/M_P)^(3/2) (mu=0 number at the same energy) [C]"]
YJ = inst["Y_J = coef (H/M_P)^(3/2) (condensate quanta conserved) [C]"]
md = b1["non-instant reheating (MD + decay, Gamma_I = gamma H_e) [C]"]
g_arr = np.array([r["gamma"] for r in md["rows"]])
fits = md["small-gamma fits R = A gamma^p"]
TRH_COEF = (90 / (math.pi ** 2 * gstar)) ** 0.25


def R_of(gam, which):
    arr = np.array([r["R_J" if which == "J" else "R_floor"] for r in md["rows"]])
    if gam >= g_arr.max():
        return 1.0
    if gam <= g_arr.min():
        f = fits["J" if which == "J" else "floor"]
        return f["A"] * gam ** f["p"]
    return float(np.exp(np.interp(math.log(gam), np.log(g_arr[::-1]), np.log(arr[::-1]))))


def Y_spec(He, TRH, which):
    gam = (TRH / TRH_COEF) ** 2 / (He * MP)
    return (YJ if which == "J" else YF) * (He / MP) ** 1.5 * R_of(gam, which)


LV2 = lphi * v ** 2


def producer(mI, TRH, He, Nmax_extra=6.0):
    Gam = H_COEF * TRH ** 2 / MP
    gam = Gam / He
    if gam > 1:
        return None

    def rhs(N, y):
        P, R, Nh, Eh, Jh = y
        a = math.exp(N)
        h = math.sqrt(max(P / a ** 3 + R / a ** 4, 1e-300) / 3)
        s = gam / h * P
        return [-s, s * a, 2 * s, s * a, 4 * s / a]

    # integrate until the inflaton is gone (gam/h >> 1 and P tiny) plus a few e-folds
    N_rh = (2 / 3) * math.log(1 / gam)
    Nend = N_rh + 12 + Nmax_extra
    Ngrid = np.linspace(0, Nend, 4000)
    sol = solve_ivp(rhs, (0, Nend), [3.0, 0.0, 0.0, 0.0, 0.0], method="LSODA", t_eval=Ngrid, rtol=1e-10,
                    atol=1e-30)
    P, R, Nh, Eh, Jh = sol.y
    a = np.exp(sol.t)
    scale = He ** 2 * MP ** 2
    rho_r = R * scale / a ** 4
    T = (rho_r / R_COEF) ** 0.25
    H = He * np.sqrt((P / a ** 3 + R / a ** 4) / 3)
    n_h = Nh * scale / (mI * a ** 3)            # per unit B
    rho_h = Eh * scale / a ** 4
    J = Jh * scale / (mI ** 2 * a ** 2)
    Pi_free = C_EFF_EQUAL * J                    # per unit B
    with np.errstate(divide="ignore", invalid="ignore"):
        Pi_th = np.where(rho_h > 0, 1.5 * C_EFF_EQUAL * n_h ** 2 / np.where(rho_h > 0, rho_h, 1.0), 0.0)
    s_sm = S_COEF * T ** 3
    i_late = -1
    Y_per_B = n_h[i_late] / s_sm[i_late]
    eps = rho_h[i_late] / rho_r[i_late]
    E_mean_over_T = (rho_h[i_late] / n_h[i_late]) / T[i_late]
    C = Y_per_B * mI / (1.5 * TRH)
    # late RD constants
    A_free_per_B = Pi_free[i_late] / H[i_late]
    A_th_per_B = Pi_th[i_late] / H[i_late]
    A_MD_formula_per_B = 24 * C_EFF_EQUAL * Gam * MP ** 2 / mI ** 2
    # extend in RD analytically: Pi ~ a^-2, H ~ a^-2 -> criterion (LV2 + H^2)/(Pi/B); include grid + RD tail
    ext = np.logspace(0, 12, 600)
    H_all = np.concatenate([H, H[i_late] / ext ** 2])
    Pf_all = np.concatenate([Pi_free, Pi_free[i_late] / ext ** 2])
    Pt_all = np.concatenate([Pi_th, Pi_th[i_late] / ext ** 2])
    T_all = np.concatenate([T, T[i_late] / ext])
    with np.errstate(divide="ignore", invalid="ignore"):
        Bf_eff = np.nanmin(np.where(Pf_all > 0, (LV2 + H_all ** 2) / Pf_all, np.inf))
        Bt_eff = np.nanmin(np.where(Pt_all > 0, (LV2 + H_all ** 2) / Pt_all, np.inf))
    Bf_noH = LV2 / np.nanmax(Pf_all)
    with np.errstate(divide="ignore", invalid="ignore"):
        i_opt = int(np.nanargmin(np.where(Pf_all > 0, (LV2 + H_all ** 2) / Pf_all, np.inf)))
    T_max = float(np.nanmax(T))
    return dict(mI=mI, TRH=TRH, He=He, gamma=gam, Gamma_I=Gam, T_max=T_max, Y_per_B=Y_per_B, C=C, eps_per_B=eps,
                E_mean_over_TSM=E_mean_over_T, xi_formal=E_mean_over_T / 3, A_free_per_B=A_free_per_B,
                A_th_per_B=A_th_per_B, A_MD_formula_per_B=A_MD_formula_per_B,
                Jfree_over_Jth_late=Pi_free[i_late] / Pi_th[i_late],
                B_min_eff_free=Bf_eff, B_min_eff_thermalized=Bt_eff, B_min_noHubble_free=Bf_noH,
                H_at_optimum=float(H_all[i_opt]), T_SM_at_optimum=float(T_all[i_opt]),
                B_max_relic=Y_MAX_STRICT / Y_per_B,
                _arrays=dict(T=T_all, H=H_all, Pf=Pf_all, n_h=np.concatenate([n_h, n_h[i_late] / ext ** 3]),
                             rho_h=np.concatenate([rho_h, rho_h[i_late] / ext ** 4])))


def edge_state(res, B):
    """State of the hidden gas for branching B: exit (last time Pi_free = lphi v^2), kinetic rate there."""
    arr = res["_arrays"]
    Pi = B * arr["Pf"]
    above = np.where(Pi > LV2)[0]
    if above.size == 0:
        return dict(restored_ever=False)
    i = above[-1]
    T_end = float(arr["T"][i])
    n = B * arr["n_h"][i]
    rho = B * arr["rho_h"][i]
    E = rho / n
    Th = E / 3
    G = gamma_el(n / 4, E, COUP_PHI) / arr["H"][i]
    # RD: Gamma/H ~ a -> kinetic at T_kin = T_end * G (if G < 1); hidden kinetic T then Th * G
    T_kin = T_end * G if G < 1 else T_end
    z = math.pi ** 2 * n / (4 * Th ** 3)
    A = B * res["A_free_per_B"]
    return dict(restored_ever=True, T_SM_end=T_end, T_h_kinetic_formal_at_end=Th, xi_formal=Th / T_end,
                z_formal=z, mu_over_T=math.log(z), Gamma_el_over_H_at_end=G, T_SM_at_kinetic_eq=T_kin,
                T_h_at_kinetic_eq=Th * T_kin / T_end, ever_kinetic_while_relativistic=bool(Th * T_kin / T_end > 100.0),
                A_GeV=A, zero_mode_responds=bool(A > 2 * math.sqrt(lphi) * v),
                cold_start_memory_erased_by_redshift_E=bool(A > 15 * v))


rows = []
for mI in (1e10, 1e11, 1e12, 1e13):
    for TRH in (1e7, 1e8, 1e9, 1e10, 1e11, 1e12):
        Gam = H_COEF * TRH ** 2 / MP
        for He_label, He in (("H_e = m_I/3", mI / 3), ("H_e = max(1e3 GeV, 30 Gamma_I)", max(1e3, 30 * Gam))):
            res = producer(mI, TRH, He)
            if res is None:
                continue
            Bmin, Bmax = res["B_min_eff_free"], res["B_max_relic"]
            row = {k: vv for k, vv in res.items() if not k.startswith("_")}
            row["H_e_choice"] = He_label
            row["validity_m_I_gt_2T_max (hard quanta, no thermal blocking)"] = bool(mI > 2 * res["T_max"])
            row["window_exists_free"] = bool(Bmin < Bmax)
            row["B_min_memory_E (A > 15 v, cold start, no damping)"] = 15 * v / res["A_free_per_B"]
            row["window_with_memory_E"] = ([15 * v / res["A_free_per_B"], Bmax]
                                           if 15 * v / res["A_free_per_B"] < Bmax else None)
            row["window_free"] = [Bmin, Bmax] if Bmin < Bmax else None
            row["edge_at_B_min"] = edge_state(res, Bmin) if Bmin < 1 else None
            row["edge_at_B_max"] = edge_state(res, Bmax) if Bmax < 1 else None
            row["Y_N_at_B_max"] = Y_MAX_STRICT
            row["spectator_Y_over_Ymax (floor, J) for this H_e [C b1]"] = [Y_spec(He, TRH, "floor") / Y_MAX_STRICT,
                                                                           Y_spec(He, TRH, "J") / Y_MAX_STRICT]
            rows.append(row)
OUT["producer grid [C]"] = rows

# analytic summary [D]: B_eff = 2 sqrt(lphi) v m_I^2/(24 c_eff k Gamma_I M_P^2), B_max = Y_max m_I/(1.5 C T_RH)
ks = [r["B_min_eff_free"] / (2 * math.sqrt(lphi) * v * r["mI"] ** 2 / (24 * C_EFF_EQUAL * r["Gamma_I"] * MP ** 2))
      for r in rows]
Cs = [r["C"] for r in rows]
k_inv = float(np.median(ks))
C_med = float(np.median(Cs))
ratio_coef = Y_MAX_STRICT * 24 * C_EFF_EQUAL * H_COEF * MP / (1.5 * C_med * 2 * math.sqrt(lphi) * v * k_inv)
OUT["analytic window [D, C fitted O(1)]"] = {
    "Y_N": "Y_N = 1.5 C B_h T_RH/m_I, C = %.4f (median; range %.4f-%.4f)" % (C_med, min(Cs), max(Cs)),
    "eps": "eps = B_h/(1-B_h) exactly",
    "xi_formal": "xi = <E>/(3 T_SM) = eps/(4 Y_N) = m_I/(6 C T_RH) (indep. of B_h); z = pi^2 n_h/(4 T_h^3) ~ 115.5 Y_N/xi^3",
    "B_min (zero mode responds)": "B_min = k * 2 sqrt(lphi) v m_I^2/(24 c_eff Gamma_I M_P^2), k = %.3f (median; range %.3f-%.3f)"
                                  % (k_inv, min(ks), max(ks)),
    "B_min numeric": "= %.3g (m_I/T_RH)^2 (v/TeV)" % (k_inv * 2 * math.sqrt(lphi) * v / (24 * C_EFF_EQUAL * H_COEF * MP)),
    "B_max numeric": "= %.3g (m_I/T_RH) (TeV/v)" % (Y_MAX_STRICT / (1.5 * C_med)),
    "window exists iff m_I/T_RH < %.3g (v/TeV)^-2" % ratio_coef: True,
    "comparison hidden_heating c5 'withHubble' B_h,max [Q]": "7.2e-14 at (H_e, m_I) = (1e9, 3e9) GeV, T_RH = 6.4e7 GeV",
}
chk = producer(3e9, 1e8, 1e9)
OUT["check vs hidden_heating c5 row (H_e=1e9, m_I=3e9, Gamma_I = H(1e8 GeV)) [C]"] = {
    "B_min_eff_free": chk["B_min_eff_free"], "B_min_noHubble_free": chk["B_min_noHubble_free"],
    "hidden_heating withHubble / noHubble [Q]": [7.219962118404513e-14, 3.406914385782353e-20],
    "note": "hidden_heating fixes Gamma_I = H(T = 1e8 GeV) and reports a rho_r-based T_RH = 6.4e7 GeV for the same run; "
            "this module defines T_RH by Gamma_I = H(T_RH)"}

OUT["headline"] = {
    "C (Y_N = 1.5 C B T_RH/m_I)": C_med,
    "k (B_min fit)": k_inv,
    "window existence m_I/T_RH <": ratio_coef,
    "examples (m_I, T_RH, H_e=m_I/3): window": [(r["mI"], r["TRH"], r["window_free"]) for r in rows
                                               if r["H_e_choice"] == "H_e = m_I/3"],
}
dump("a2_producer.json", OUT)
if __name__ == "__main__":
    print(json.dumps(OUT["analytic window [D, C fitted O(1)]"], indent=1))
    print(json.dumps(OUT["check vs hidden_heating c5 row (H_e=1e9, m_I=3e9, Gamma_I = H(1e8 GeV)) [C]"], indent=1, default=float))
    for r in rows:
        print("%.0e %.0e %-28s gam=%.2e C=%.3f xi=%.3g Tmax=%.2e Bmin=%.2e Bmax=%.2e win=%s spec=%s" % (
            r["mI"], r["TRH"], r["H_e_choice"], r["gamma"], r["C"], r["xi_formal"], r["T_max"],
            r["B_min_eff_free"], r["B_max_relic"], r["window_exists_free"],
            ["%.2g" % q for q in r["spectator_Y_over_Ymax (floor, J) for this H_e [C b1]"]]))
