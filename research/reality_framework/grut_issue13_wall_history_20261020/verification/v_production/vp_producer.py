#!/usr/bin/env python3
"""v_production independent check, part (4): heavy-inflaton producer map (outside-action preparation data [P]).

Perturbative decay of a matter-dominated inflaton, Gamma = H(T_RH) = sqrt(pi^2 g*/90) T_RH^2/M_P, g* = 106.75,
branching B into hidden pairs, each quantum injected with E = m_I/2, no hidden back-reaction on the decay.
Dimensionless: rho in units 3 M_P^2 Gamma^2, time via N = ln a, h = H/Gamma.
  r_I' = -3 r_I - r_I/h ; r_r' = -4 r_r + (1-B) r_I/h ; r_h' = -4 r_h + B r_I/h
  n~'  = -3 n~ + 2B r_I/h            (n_h = n~ 3 M_P^2 Gamma^2/m_I)
  J~'  = -2 J~ + 4B r_I/h            (J = sum_dof int d^3k f/E/(2pi)^3 = J~ 3 M_P^2 Gamma^2/m_I^2)
Late time: T/T_RH = r_r^(1/4) exactly (from the T_RH definition), so
  Y m_I/(B T_RH) = (3/4) n~/(B r_r^(3/4))  (g* cancels), C = Y/(1.5 B T_RH/m_I).
Curvature at phi = 0 of the free-streaming gas: Pi = c_eff J, A = Pi/H.
Written from scratch by the verifier.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp

v = 1000.0; lphi = 0.01; gX = 0.1; hD = 0.01
c_eff = (3 * lphi + 2 * gX**2 + hD**2) / 4
MP = 2.435323e18; gstar = 106.75
Ymax = 4.373e-12
sq = math.sqrt(math.pi**2 * gstar / 90)


def run(he, B=1e-6, Nmax=None):
    # rescaled comoving variables to avoid stiffness in magnitude
    def rhs(N, y):
        RI, Rr, Rh, n, J = y
        rI = RI * math.exp(-3 * N)
        h = math.sqrt(max(rI + (Rr + Rh) * math.exp(-4 * N), 1e-300))
        src = rI / h
        return [-RI / h,
                (1 - B) * src * math.exp(4 * N),
                B * src * math.exp(4 * N),
                2 * B * src * math.exp(3 * N),
                4 * B * src * math.exp(2 * N)]
    if Nmax is None:
        Nmax = (2 / 3) * math.log(he) + 4.5
    y0 = [he**2, 0.0, 0.0, 0.0, 0.0]
    # absolute tolerances scaled to each variable's natural size (R_I decays to ~0 after Gamma t >> 1)
    atol = [he**2 * 1e-22, he**2 * 1e-22 * math.exp(4 * Nmax), B * he**2 * 1e-22 * math.exp(4 * Nmax),
            B * he**2 * 1e-22 * math.exp(3 * Nmax), B * he**2 * 1e-22 * math.exp(2 * Nmax)]
    sol = solve_ivp(rhs, [0, Nmax], y0, method="LSODA", rtol=1e-11, atol=atol, dense_output=True)
    return sol


def late(he, B=1e-6):
    sol = run(he, B)
    N = sol.t[-1]
    RI, Rr, Rh, n, J = sol.y[:, -1]
    rr = Rr * math.exp(-4 * N); rh = Rh * math.exp(-4 * N); nt = n * math.exp(-3 * N); Jt = J * math.exp(-2 * N)
    h = math.sqrt(RI * math.exp(-3 * N) + rr + rh)
    C = 0.75 * nt / (B * rr**0.75) / 1.5
    kRD = (Jt / h) / (8 * B)
    EinvE = rh * Jt / nt**2
    return dict(he=he, B=B, C=C, coef_1p5C=1.5 * C, eps=rh / rr, eps_expected=B / (1 - B),
                A_RD_over_A_MD=kRD, E_times_invE=EinvE, E_over_TSM_per_mI_over_TRH=(rh / nt) / rr**0.25,
                residual_inflaton_fraction=RI * math.exp(-3 * N) / (rr + rh))


out = {"c_eff": c_eff, "sqrt_pi2gstar_over_90": sq}
out["late_time"] = [late(he) for he in [30.0, 1e2, 1e4, 1e8, 1e12]]
out["late_time_B0.1"] = late(1e8, B=0.1)

# MD analytic A = 24 c_eff B Gamma M_P^2/m_I^2 and the zero-mode criterion A > 2 sqrt(lam) v (max_H [A H - lam v^2 - H^2])
Bmin_coef = 2 * math.sqrt(lphi) * v / (24 * c_eff * sq * MP)       # B_min = coef (m_I/T_RH)^2  (A_MD)
C_inf = out["late_time"][3]["C"]
Bmax_coef = Ymax / (1.5 * C_inf)                                       # B_max = coef (m_I/T_RH)
out["Bmin_coef_MD"] = Bmin_coef
out["Bmin_coef_RD"] = Bmin_coef / out["late_time"][3]["A_RD_over_A_MD"]
out["Bmax_coef"] = Bmax_coef
out["Bmax_coef_C_at_he30"] = Ymax / (1.5 * out["late_time"][0]["C"])
out["window_closes_at_mI_over_TRH_MD"] = Bmax_coef / Bmin_coef
out["window_closes_at_mI_over_TRH_RD"] = Bmax_coef / out["Bmin_coef_RD"]
out["H_at_which_H_eq_100GeV_T_GeV"] = math.sqrt(100 * MP / sq)


# full numerical criterion for physical cells: max over t of [c_eff J - lam v^2 - H^2] > 0, H_e = m_I/3
def crit_cell(mI, TRH, He=None, crit_kappa=1.0):
    Gam = sq * TRH**2 / MP
    if He is None:
        He = mI / 3
    he = He / Gam
    N0 = (2 / 3) * math.log(he) + 4.5
    Hend0 = Gam * math.exp(-2 * 4.5)                       # rough H at N0 (RD, h ~ e^{-2 dN})
    Nmax = N0 + max(0.0, 0.5 * math.log(Hend0 / 5.0))      # continue until H ~ 5 GeV < lam v^2/A_opt scale
    sol = run(he, B=1.0e-12, Nmax=Nmax)  # J linear in B (no backreaction); B tiny
    Ns = np.linspace(0, sol.t[-1], 200001)
    Y = sol.sol(Ns)
    RI, Rr, Rh, n, J = Y
    rI = RI * np.exp(-3 * Ns); rr = Rr * np.exp(-4 * Ns); rh = Rh * np.exp(-4 * Ns)
    h = np.sqrt(rI + rr + rh)
    H = h * Gam
    Jphys = J * np.exp(-2 * Ns) * 3 * MP**2 * Gam**2 / mI**2   # for B=1e-12
    Pi1 = c_eff * Jphys / 1e-12                                 # per unit B
    # need B such that max_t [B Pi1 - lam v^2 - kappa H^2] = 0
    # For each t: B_t = (lam v^2 + kappa H^2)/Pi1 ; B_min = min_t B_t
    with np.errstate(divide="ignore", invalid="ignore"):
        Bt = (lphi * v * v + crit_kappa * H**2) / Pi1
    Bt[~np.isfinite(Bt)] = np.inf
    i = int(np.argmin(Bt))
    return dict(mI=mI, TRH=TRH, He=He, Bmin=float(Bt[i]), Bmin_over_x2=float(Bt[i]) / (mI / TRH)**2,
                H_at_opt=float(H[i]), phase="MD" if rI[i] > rr[i] else "RD")


cells = []
for mI, TRH in [(1e10, 1e9), (1e10, 1e7), (1e11, 1e10), (1e12, 1e11), (1e13, 1e12), (1e13, 1e9), (1e12, 1e8)]:
    cells.append(crit_cell(mI, TRH))
out["Bmin_cells_Hubble_damped_kappa1"] = cells
out["Bmin_cells_kappa_9over4"] = [crit_cell(1e10, 1e9, crit_kappa=2.25), crit_cell(1e13, 1e12, crit_kappa=2.25)]
# hidden_heating c5 cross-check quoted by the module: Gamma = H(1e8), He = 1e9, mI = 3e9
c5 = crit_cell(3e9, 1e8, He=1e9)
out["c5_check"] = c5

# Y_N, xi_formal, z in a sample cell, B at B_max
lt = out["late_time"][3]
out["xi_formal_coef"] = 1 / (6 * lt["C"])   # xi = <E>/(3 T_SM) = m_I/(6 C T_RH) for B<<1
out["z_coef"] = math.pi**2 * (2 * math.pi**2 * gstar / 45) / 4   # z = coef Y/xi^3 (MB, 4 dof)
mu = []
for x in [10, 100, 1e3, 1e4, 3.8e4]:
    Y = Ymax; xi = out["xi_formal_coef"] * x
    z = out["z_coef"] * Y / xi**3
    mu.append(dict(mI_over_TRH=x, xi_formal=xi, z_at_Bmax=z, mu_over_T_at_Bmax=math.log(z)))
out["mu_at_Bmax"] = mu

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vp_producer.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
