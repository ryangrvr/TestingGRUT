"""production_map / a1: Z2 restoration by a kinetically thermal hidden bath at FIXED number (mu_N <= 0).

Output: results/a1_restoration.json.  Deterministic, < 1 min.

Physics [D]:
  * Hartree (one-loop gap) curvature at phi = 0, common fugacity z = exp(mu_N/T_h) <= 1 for phi, q1, q2, D,
    Bose-Einstein distributions with quasi-particle masses (m_phi = 0 at the boundary, m_q and m_D self-consistent,
    m_D0 = 100 GeV kept).  Restoration iff Pi_phi(T_h, z) > lphi v^2.
  * Conserved hidden yield Y_N = n_h/s_SM and hidden energy fraction eps = rho_h/rho_SM (both constant while the
    bath is relativistic and the SM g* is constant).  Given (Y_N, eps) the bath state is (T_h, z) and the exit is
    at T_SM,end where Pi_phi = lphi v^2.
  * In the massless limit eps/Y_N^(4/3) = F(z) and T_SM,end = G(z) Y_N^(-1/3); F decreases monotonically with z.
"""
import math

import numpy as np
from scipy.optimize import brentq

from common import (v, lphi, MP, S_COEF, R_COEF, H_COEF, Y_MAX_STRICT, Y_MAX_THERMAL_COMPOSITION, TC_RESUMMED,
                    TC_LEADING, C_EFF_EQUAL, COUP_PHI, COUP_X, gh, gstar, zeta3, hartree_phi0, Tc_of_z, gamma_el,
                    H_of_T, T_of_H, dump)

OUT = {"module": "production_map", "part": "A1 restoration criterion for a fixed-number bath",
       "labels": "[D] derived, [C] computed, [Q] quoted, [E] estimate, [P] preparation datum"}

# ------------------------------------------------------------------ 1. T_c(z) and the per-z invariants
zs = np.concatenate([[1.0, 0.99, 0.97, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03,
                      0.02, 0.01], 10.0 ** np.arange(-2.5, -16.01, -0.5)])
rows = []
for z in zs:
    T = Tc_of_z(float(z))
    h = hartree_phi0(T, float(z))
    n, rho = h["n"], h["rho"]
    F = rho / R_COEF / (n / S_COEF) ** (4 / 3)          # eps / Y_N^(4/3)
    G = (n / S_COEF) ** (1 / 3)                           # T_SM,end * Y_N^(1/3)   [GeV]
    E_mean = rho / n
    n_dof_phi = h["n_phi"]
    rows.append(dict(z=float(z), mu_over_T=float(math.log(z)), Tc_over_v=T / v, mq_over_T=h["mq"] / T,
                     mD_over_T=h["mD"] / T, n_over_T3=n / T ** 3, rho_over_T4=rho / T ** 4, F=F, G=G,
                     E_mean_over_T=E_mean / T, n_phi_over_T3=n_dof_phi / T ** 3))
OUT["Tc_of_z [C] (Hartree, m_phi=0 at boundary, Bose, thermal masses)"] = rows

lnz = np.log([r["z"] for r in rows])
lnF = np.log([r["F"] for r in rows])
lnG = np.log([r["G"] for r in rows])
order = np.argsort(lnF)


def z_of(Y, eps):
    """fugacity of a kinetically thermal bath with conserved (Y_N, eps); None if overpopulated (eps < eps_th)."""
    F = eps / Y ** (4 / 3)
    if F < rows[0]["F"] * (1 - 1e-12):
        return None
    if F > rows[-1]["F"]:
        # MB asymptote: F ~ z^(-1/3) (rho/n^(4/3) ~ z^-1/3) -> extrapolate
        return math.exp(lnz[-1] - 3 * (math.log(F) - lnF[-1]))
    return float(np.exp(np.interp(math.log(F), lnF[order], lnz[order])))


def exit_state(Y, eps):
    """(T_SM,end, T_h,end, z) for a relativistic bath with conserved (Y_N, eps).  Direct root solve, not table."""
    z0 = z_of(Y, eps)
    if z0 is None:
        return None

    def mismatch(lz):
        z = math.exp(lz)
        T = Tc_of_z(z)
        h = hartree_phi0(T, z)
        Tsm = (h["n"] / (Y * S_COEF)) ** (1 / 3)
        return math.log(h["rho"] / (eps * R_COEF * Tsm ** 4))

    lz0 = math.log(z0)
    if z0 > 0.999999:
        lz = 0.0 if abs(mismatch(0.0)) < 1e-6 else brentq(mismatch, lz0 - 0.5, 0.0, xtol=1e-12)
    else:
        lo, hi = lz0 - 0.3, min(lz0 + 0.3, 0.0)
        lz = brentq(mismatch, lo, hi, xtol=1e-12)
    z = math.exp(lz)
    T = Tc_of_z(z)
    h = hartree_phi0(T, z)
    Tsm = (h["n"] / (Y * S_COEF)) ** (1 / 3)
    return dict(T_SM_end=Tsm, T_h_end=T, z=z, mu_over_T=lz, xi_end=T / Tsm, n=h["n"], rho=h["rho"],
                n_phi=h["n_phi"])


def eps_thermal(Y):
    return rows[0]["F"] * Y ** (4 / 3)


# ------------------------------------------------------------------ 2. checks
K_MB = 1.5 * C_EFF_EQUAL * S_COEF ** 2 / R_COEF        # Pi = K_MB (Y^2/eps) T_SM^2  (massless MB) [D]
checks = {
    "Tc(z=1)/v Hartree [C]": rows[0]["Tc_over_v"],
    "Tc(z=1) vs resummed 1.5571-1.5641 [Q thermal_potential] and charge_phase T_c0 = 1.5575 [Q]":
        "inside the resummed range" if TC_RESUMMED[0] <= rows[0]["Tc_over_v"] <= TC_RESUMMED[1] else "outside",
    "Tc(z)->Tc_leading*sqrt(zeta2/Li2(z)) as z->0 (thermal masses vanish) [C]":
        rows[-1]["Tc_over_v"] / math.sqrt((math.pi ** 2 / 6) / rows[-1]["z"]) / TC_LEADING,
    "thermal_potential quoted scaling Tc ~ sqrt(zeta2/Li2): z=0.5 -> x1.681, 0.1 -> x4.00 [Q]; Hartree here [C]": {
        str(r["z"]): r["Tc_over_v"] / rows[0]["Tc_over_v"] for r in rows if r["z"] in (0.9, 0.5, 0.1)},
    "MB closed form K_MB = (3/2) c_eff s^2/rho_SM (Pi = K (Y^2/eps) T_SM^2) [D]": K_MB,
    "c_eff (equal sharing) [D]": C_EFF_EQUAL,
}

# ------------------------------------------------------------------ 3. minimum T_SM for a relic-compatible restored bath
res_min = {}
for label, Y in (("strict 4.373e-12 [Q v_relic]", Y_MAX_STRICT),
                 ("thermal composition 4.30e-12 [Q relic_yield]", Y_MAX_THERMAL_COMPOSITION)):
    st = exit_state(Y, eps_thermal(Y))
    scale = [tc / st["T_h_end"] * v for tc in TC_RESUMMED]   # T_SM,end ∝ T_c at fixed z
    res_min[label] = dict(T_SM_min_Hartree_GeV=st["T_SM_end"], xi_at_exit=st["xi_end"],
                          T_SM_min_rescaled_to_resummed_Tc_GeV=[st["T_SM_end"] * s for s in scale],
                          eps_thermal=eps_thermal(Y))
OUT["minimum T_SM of a relic-compatible restored phase (z = 1) [C]"] = res_min
OUT["comparison relic_yield mu_N=0 [Q]"] = ("relic_yield: T_SM >= 2.04-2.10e6 GeV (v/TeV)^(4/3); v_relic route B with "
                                            "the strict Y: 2.066-2.088e6 GeV")

# ------------------------------------------------------------------ 4. T_SM,end(Y_N, eps) grid, kinetic and Hubble tests
T50 = T_of_H(math.sqrt(lphi) * v / 2)          # Pi = A H with A const: zero mode can respond iff A > 2 sqrt(l) v
Tmem = T_of_H(lphi * v ** 2 / (15 * v))        # [E] cold-start coherent amplitude below vacuum soft fluctuations
grid = []
for Y in (1e-16, 1e-15, 1e-14, 1e-13, 1e-12, Y_MAX_STRICT):
    e_th = eps_thermal(Y)
    for ratio in (1.0, 1.01, 1.1, 2.0, 10.0, 1e2, 1e3, 1e4, 1e6, 1e8, 1e10, 1e12):
        eps = e_th * ratio
        st = exit_state(Y, eps)
        Hend = H_of_T(st["T_SM_end"])
        Eh = st["rho"] / st["n"]
        nd_phi = st["n_phi"]
        g_phi = gamma_el(nd_phi, Eh, COUP_PHI) / Hend
        g_X = gamma_el(nd_phi, Eh, COUP_X) / Hend
        A = lphi * v ** 2 / Hend
        mb = math.sqrt(lphi * v ** 2 * eps / (K_MB * Y ** 2))
        grid.append(dict(Y_N=Y, eps=eps, eps_over_eps_th=ratio, z=st["z"], mu_over_T=st["mu_over_T"],
                         T_h_end_GeV=st["T_h_end"], xi_end=st["xi_end"], T_SM_end_GeV=st["T_SM_end"],
                         T_SM_end_MB_closed_form=mb, Gamma_phi_el_over_H_at_exit=g_phi,
                         Gamma_X_el_over_H_at_exit=g_X, A_equals_Pi_over_H_GeV=A,
                         effective_restoration_A_gt_2sqrt_lphi_v=A > 2 * math.sqrt(lphi) * v,
                         memory_erased_by_redshift_E=A > 15 * v))
OUT["exit grid T_SM,end(Y_N, eps) [C]"] = grid
OUT["Hubble-damped criterion [D]"] = {
    "statement": "In RD with conserved (Y_N, eps) and a relativistic bath, Pi_phi = A H with A constant. The zero "
                 "mode responds (Pi - lphi v^2 > H^2 at some time) iff A > 2 sqrt(lphi) v = 200 GeV, i.e. iff "
                 "T_SM,end < T(H = sqrt(lphi) v/2), and the bath existed at H = A/2.",
    "T_SM_end_max_for_effective_restoration_GeV": T50,
    "memory [E]": "cold +v start, no damping: coherent amplitude at the exit ~ 1.7 v sqrt(lphi) v/A; below the vacuum "
                  "soft-mode level 0.11 sqrt(lphi) v needs A > ~15 v, i.e. T_SM,end < %.3g GeV" % Tmem,
}

# ------------------------------------------------------------------ 5. region: kinetic + relic + restoring
# For each Y_N find the largest eps for which Gamma_phi,el/H >= 1 at the exit (kinetically thermal exit).
region = []
for Y in (1e-14, 1e-13, 3e-13, 1e-12, 2e-12, Y_MAX_STRICT):
    e_th = eps_thermal(Y)

    def g(lr):
        st = exit_state(Y, e_th * math.exp(lr))
        Hend = H_of_T(st["T_SM_end"])
        return math.log(gamma_el(st["n_phi"], st["rho"] / st["n"], COUP_PHI) / Hend)

    g0 = g(0.0)
    if g0 < 0:
        region.append(dict(Y_N=Y, kinetic_at_exit_even_for_mu0=False, Gamma_over_H_mu0=math.exp(g0)))
        continue
    lr = brentq(g, 0.0, 60.0, xtol=1e-8)
    st = exit_state(Y, e_th * math.exp(lr))
    st0 = exit_state(Y, e_th)
    region.append(dict(Y_N=Y, Gamma_over_H_mu0=math.exp(g0), eps_max_kinetic=e_th * math.exp(lr),
                       eps_max_over_eps_th=math.exp(lr), z_min=st["z"], mu_over_T_min=st["mu_over_T"],
                       xi_end_max=st["xi_end"], T_h_end_max_GeV=st["T_h_end"],
                       T_SM_end_range_GeV=[st0["T_SM_end"], st["T_SM_end"]]))
OUT["kinetically thermal exit region (Gamma_phi,el/H >= 1 at the exit) [C, E rates]"] = region

# xi_kin check at z = 1 against thermal_potential's 6.6e-5
def gk(lxi):
    xi = math.exp(lxi)
    T = rows[0]["Tc_over_v"] * v
    nphi = rows[0]["n_phi_over_T3"] * T ** 3
    E = rows[0]["E_mean_over_T"] * T
    return math.log(gamma_el(nphi, E, COUP_PHI) / H_of_T(T / xi))


checks["xi_kin(phi) at z=1 [C] vs 6.6e-5 [Q thermal_potential]"] = math.exp(brentq(gk, -20, -2))
OUT["checks"] = checks

# ------------------------------------------------------------------ 6. T_SM,end in terms of (Y_N, xi) [D]
xi_rows = []
for Y in (1e-14, 1e-13, 1e-12, Y_MAX_STRICT):
    for xi in (1e-4, 3e-4, 7.5e-4, 1e-3, 3e-3, 1e-2, 0.1, 1.0):
        # massless relation n_h/T_h^3 = Y S_COEF / xi^3 -> z; bracket via the table
        target = Y * S_COEF / xi ** 3
        nT3 = np.array([r["n_over_T3"] for r in rows])
        if target > nT3[0]:
            xi_rows.append(dict(Y_N=Y, xi=xi, state="overpopulated (needs condensate) for this (Y_N, xi)"))
            continue
        lzz = float(np.interp(math.log(target), np.log(nT3[::-1]), lnz[::-1])) if target > nT3[-1] else \
            lnz[-1] + math.log(target / nT3[-1])
        z = math.exp(lzz)
        T = Tc_of_z(z)
        xi_rows.append(dict(Y_N=Y, xi=xi, z=z, mu_over_T=lzz, T_h_end_GeV=T, T_SM_end_GeV=T / xi))
OUT["T_SM,end(Y_N, xi) [C]"] = xi_rows

# ------------------------------------------------------------------ 7. overpopulated bath (eps < eps_th) [D, E]
st1 = exit_state(Y_MAX_STRICT, eps_thermal(Y_MAX_STRICT))
lX_ = 0.001
A2 = lphi * v ** 2 / (0.457 * 0.01)              # g_X^2 <q1^2> = 0.457 gX2 A^2 > lphi v^2   (quartic oscillation)
J_cond = (4 / 3) * (lX_ * A2 ** 2 / 4) / (0.8472 * math.sqrt(lX_) * math.sqrt(A2))
T_cond = (J_cond / (Y_MAX_STRICT * S_COEF)) ** (1 / 3)
OUT["overpopulated bath [D, E]"] = {
    "gas bound [D]": "eps < eps_th(Y_N): kinetic equilibrium is a mu = m_min gas at the energy temperature plus a "
                     "condensate; the gas alone gives Pi >= Pi_th(T_E), so restoration lasts at least to "
                     "T_SM = T_c(1)/xi_E, xi_E = (eps g*/g_h)^(1/4).",
    "uncharged X-condensate placeholder [E]": dict(n_cond_min_over_v3=J_cond / v ** 3, T_SM_min_GeV=T_cond,
                                                   note="cold, uncharged, linearly oscillating X condensate at phi=0 "
                                                        "with g_X^2<q1^2> = lphi v^2; unstable to fragmentation; not a "
                                                        "computed bound"),
}

OUT["headline"] = {
    "Tc(z=1)/v Hartree": rows[0]["Tc_over_v"],
    "T_SM_min relic-compatible restored phase, strict Y (Hartree; rescaled resummed range)":
        [res_min["strict 4.373e-12 [Q v_relic]"]["T_SM_min_Hartree_GeV"]] +
        res_min["strict 4.373e-12 [Q v_relic]"]["T_SM_min_rescaled_to_resummed_Tc_GeV"],
    "T_SM,end scaling, MB limit": "T_SM,end = sqrt(lphi v^2 eps/(K Y_N^2)), K = %.4f -> %.4g GeV sqrt(eps)/Y_N (v/TeV)"
                                  % (K_MB, math.sqrt(lphi * v ** 2 / K_MB)),
    "underpopulated bath at fixed Y_N exits EARLIER (higher T_SM) than mu_N = 0": True,
    "T_SM,end max for zero-mode response (Hubble-damped)": T50,
}

dump("a1_restoration.json", OUT)
if __name__ == "__main__":
    import json
    print(json.dumps(OUT["headline"], indent=1))
    print(json.dumps(OUT["checks"], indent=1, default=float))
    for r in region:
        print(r)
