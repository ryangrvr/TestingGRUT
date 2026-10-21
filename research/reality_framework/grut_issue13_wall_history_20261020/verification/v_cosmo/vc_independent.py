#!/usr/bin/env python3
"""Independent verification (v_cosmo) of the cosmo_observable module.

Written from scratch; does not import or copy module code.
numpy/scipy only, deterministic (no random numbers).

Units: natural (hbar = c = k_B = 1), energies in GeV.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vc_results.json")

# ---------------------------------------------------------------- constants (CODATA 2018 / PDG)
HBAR_GeVs = 6.582119569e-25          # GeV s
HBARC_GeVm = 1.973269804e-16         # GeV m
C_ms = 2.99792458e8                  # m/s
GeV_J = 1.602176634e-10              # J
KB_GeVK = 8.617333262e-14            # GeV/K
MPC_m = 3.0856775814913673e22        # m
G_SI = 6.67430e-11                   # m^3 kg^-1 s^-2
YR_s = 365.25 * 86400.0
ME = 0.51099895e-3                   # GeV

MP = 2.435323e18                     # reduced Planck mass (BRIEF)
SIGMA_BENCH = 9.428090416e7          # GeV^3 (BRIEF)
T0 = 2.7255 * KB_GeVK                # GeV
KMSMPC_GeV = 1e3 / MPC_m * HBAR_GeVs  # 1 km/s/Mpc in GeV
RHO_C100 = 3 * MP**2 * (100 * KMSMPC_GeV) ** 2
OMEGA_M_H2 = 0.02237 + 0.1200        # baryons + CDM (Planck 2018 base), massless nu

# ---------------------------------------------------------------- SM thermodynamics below ~10 MeV
# photons + e+e- (exact FD integrals, mu=0) + 3 nu decoupled instantaneously at T_dec


def _fd_int(x, kind):
    """Dimensionless e+- integrals per internal d.o.f., x = m/T. kind: 'rho' or 'p'."""
    if kind == "rho":
        f = lambda u: u * u * math.sqrt(u * u + x * x) / (math.exp(math.sqrt(u * u + x * x)) + 1.0)
    else:
        f = lambda u: u**4 / math.sqrt(u * u + x * x) / (math.exp(math.sqrt(u * u + x * x)) + 1.0) / 3.0
    val, _ = quad(f, 0, 60 + x, limit=200, epsabs=0, epsrel=1e-12)
    return val / (2 * math.pi**2)  # rho/T^4 or p/T^4 per d.o.f.


def e_pm_rho_p_s(T):
    x = ME / T
    if x > 200:
        return 0.0, 0.0, 0.0
    r = 4 * _fd_int(x, "rho") * T**4
    p = 4 * _fd_int(x, "p") * T**4
    return r, p, (r + p) / T


def build_sm_tables(T_dec=None, Neff=3.0, Tmax=5.0e-3, Tmin=1.0e-14, npts=1200):
    """Return splines in lnT of: ln g_s_total, rho_r(T)/T^4 (gamma+e+nu), (T_nu/T)."""
    lnT = np.linspace(math.log(Tmin), math.log(Tmax), npts)
    T = np.exp(lnT)
    gsge = np.empty_like(T)   # entropy dof of gamma+e plasma
    rho_ge = np.empty_like(T)
    for i, Ti in enumerate(T):
        re, pe, se = e_pm_rho_p_s(Ti)
        sg = 4 * (math.pi**2 / 45) * Ti**3  # photon entropy: (2pi^2/45)*2*T^3
        rg = 2 * (math.pi**2 / 30) * Ti**4
        gsge[i] = (sg + se) / ((2 * math.pi**2 / 45) * Ti**3)
        rho_ge[i] = rg + re
    # value at decoupling
    if T_dec is not None:
        re, pe, se = e_pm_rho_p_s(T_dec)
        gsge_dec = 2 + se / ((2 * math.pi**2 / 45) * T_dec**3)
    if T_dec is None:
        # instantaneous decoupling while e+- fully relativistic: (T_nu/T)^3 = g_s,ge(T)/(11/2)
        ratio3 = np.minimum(1.0, gsge / 5.5)
    else:
        ratio3 = np.where(T >= T_dec, 1.0, gsge / gsge_dec)  # (T_nu/T)^3
    Tnu_T = ratio3 ** (1.0 / 3.0)
    gs_nu = 3 * 2 * (7.0 / 8.0) * ratio3                 # 3 species x (nu, nubar), entropy (not rescaled)
    gs_tot = gsge + gs_nu
    rho_nu = Neff / 3.0 * 3 * 2 * (7.0 / 8.0) * (math.pi**2 / 30) * (Tnu_T * T) ** 4
    rho_r = rho_ge + rho_nu
    return dict(
        lnT=lnT,
        ln_gs=CubicSpline(lnT, np.log(gs_tot)),
        ln_rhoT4=CubicSpline(lnT, np.log(rho_r / T**4)),
        TnuT=CubicSpline(lnT, Tnu_T),
        gs0=float(gs_tot[0]),
        grho0=float(rho_r[0] / T[0] ** 4 / (math.pi**2 / 30)),
    )


SM = build_sm_tables()
SM44 = build_sm_tables(Neff=3.044)
SM_dec2 = build_sm_tables(T_dec=2.0e-3)  # variant: instantaneous decoupling at 2 MeV
# rho_Lambda closing a no-wall universe at h = 0.6736 (Planck-like)
OMEGA_L_H2_PLANCK = 0.6736**2 - OMEGA_M_H2 - SM44["grho0"] * (math.pi**2 / 30) * T0**4 / RHO_C100


def rho_r(T, sm=SM):
    return math.exp(sm["ln_rhoT4"](math.log(T))) * T**4


def ln_a(T, sm=SM):
    """ln(a/a0) from entropy conservation (exact below Tmax of the table)."""
    return (math.log(sm["gs0"]) - sm["ln_gs"](math.log(T))) / 3.0 + math.log(T0 / T)


def dlna_dlnT(T, sm=SM):
    return -1.0 - sm["ln_gs"](math.log(T), 1) / 3.0


def rho_1nu(T, sm=SM):
    Tnu = sm["TnuT"](math.log(T)) * T
    return 2 * (7.0 / 8.0) * (math.pi**2 / 30) * Tnu**4


def rho_m(T, sm=SM):
    return OMEGA_M_H2 * RHO_C100 * math.exp(-3 * ln_a(T, sm))


# ---------------------------------------------------------------- background with walls
def evolve(sigma, A=0.8, conv="t", backreact=True, rhoL=0.0, T_start=5.0e-3, T_end=None,
           switch=None, n_post=1, sm=SM, rtol=1e-11, dense=False):
    """Integrate t(lnT).  conv='t': rho_w = A sigma / t ; conv='H': rho_w = 2 A sigma H (self-consistent).
    switch: None (scale to T_end) or ('rad', f) -> switch when rho_w = f rho_r, or ('tot', f) -> rho_w = f*(rho_r+rho_m+rhoL).
    After switch rho_w = rho_sw (a_sw/a)^n_post.  Returns dict.
    """
    if T_end is None:
        T_end = T0

    def H_of(T, t, rw_override=None):
        rr = rho_r(T, sm)
        rm = rho_m(T, sm)
        if not backreact:
            rho_other = rr + rm + rhoL
            return math.sqrt(rho_other / 3.0) / MP, rr, rm
        return None, rr, rm

    def rho_w_scaling(T, t):
        rr = rho_r(T, sm)
        rm = rho_m(T, sm)
        rest = rr + rm + rhoL
        if conv == "t":
            rw = A * sigma / t
            H = math.sqrt((rest + (rw if backreact else 0.0)) / 3.0) / MP
        else:  # rho_w = 2 A sigma H ; with backreaction solve 3 MP^2 H^2 = rest + 2 A sigma H
            if backreact:
                b = 2 * A * sigma
                H = (b + math.sqrt(b * b + 12 * MP**2 * rest)) / (6 * MP**2)
            else:
                H = math.sqrt(rest / 3.0) / MP
            rw = 2 * A * sigma * H
        return rw, H, rr, rm

    # initial t: radiation era t = 1/(2H) at T_start (error decays as (a_start/a)^2)
    rr0 = rho_r(T_start, sm)
    t_init = 1.0 / (2.0 * math.sqrt(rr0 / 3.0) / MP)

    def rhs(lnT, y):
        T = math.exp(lnT)
        t = math.exp(y[0])
        rw, H, rr, rm = rho_w_scaling(T, t)
        return [dlna_dlnT(T, sm) / (H * t)]

    events = []
    if switch is not None:
        kind, f = switch

        def ev(lnT, y):
            T = math.exp(lnT)
            rw, H, rr, rm = rho_w_scaling(T, math.exp(y[0]))
            ref = rr if kind == "rad" else (rr + rm + rhoL)
            return math.log(rw) - math.log(f * ref)
        ev.terminal = True
        ev.direction = 1
        events = [ev]

    sol = solve_ivp(rhs, [math.log(T_start), math.log(T_end)], [math.log(t_init)], method="DOP853",
                    rtol=rtol, atol=1e-14, events=events, dense_output=True)
    out = dict(status=sol.status)
    if switch is not None and sol.t_events[0].size:
        lnTs = sol.t_events[0][0]
        Ts = math.exp(lnTs)
        ts = math.exp(sol.y_events[0][0][0])
        rws, Hs, rrs, rms = rho_w_scaling(Ts, ts)
        out.update(T_sw=Ts, t_sw=ts, rho_w_sw=rws, H_sw=Hs, lna_sw=ln_a(Ts, sm))
        # post-switch: analytic rho_w(a); integrate t to T_end for the age
        lna_sw = ln_a(Ts, sm)

        def rw_post(T):
            return rws * math.exp(-n_post * (ln_a(T, sm) - lna_sw))

        def H_post(T):
            return math.sqrt((rho_r(T, sm) + rho_m(T, sm) + rhoL + rw_post(T)) / 3.0) / MP

        def rhs2(lnT, y):
            T = math.exp(lnT)
            return [dlna_dlnT(T, sm) / H_post(T)]
        sol2 = solve_ivp(rhs2, [lnTs, math.log(T_end)], [ts], method="DOP853", rtol=1e-11, atol=1e-30)
        out.update(rho_w_end=rw_post(T_end), H_end=H_post(T_end), t_end=float(sol2.y[0][-1]),
                   rho_r_end=rho_r(T_end, sm), rho_m_end=rho_m(T_end, sm))
    else:
        tE = math.exp(sol.y[0][-1])
        rw, H, rr, rm = rho_w_scaling(T_end, tE)
        out.update(rho_w_end=rw, H_end=H, t_end=tE, rho_r_end=rr, rho_m_end=rm)
    if dense:
        out["sol"] = sol
        out["rho_w_scaling"] = rho_w_scaling
    return out


def T_dom_numeric(sigma, A=0.8, conv="t", backreact=True, sm=SM, rtol=1e-11):
    r = evolve(sigma, A=A, conv=conv, backreact=backreact, switch=("rad", 1.0), sm=sm, rtol=rtol,
               T_end=1e-9)
    return r["T_sw"], r["t_sw"], r


def T_dom_closed(sigma, A, g):
    return (40 * A**2 * sigma**2 / (math.pi**2 * g * MP**2)) ** 0.25


def main():
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    res = {"inputs": dict(sigma=SIGMA_BENCH, M_P=MP, T0_GeV=T0, A=0.8, KMSMPC_GeV=KMSMPC_GeV,
                          rho_c100=RHO_C100, omega_m=OMEGA_M_H2,
                          g_rho0_Neff3=SM["grho0"], g_s0_Neff3=SM["gs0"],
                          g_rho0_Neff3044=SM44["grho0"], g_s0_Neff3044=SM44["gs0"])}
    sig = SIGMA_BENCH
    # ---------------- (1) T_dom
    g3 = 2 + 3 * 2 * (7 / 8) * (4 / 11) ** (4 / 3)
    g3044 = 2 + 3.044 * 2 * (7 / 8) * (4 / 11) ** (4 / 3)
    d1 = {}
    d1["g_rho_analytic_Neff3"] = g3
    d1["g_rho_analytic_Neff3044"] = g3044
    d1["Tdom_closed_noBR_g3_keV"] = T_dom_closed(sig, 0.8, g3) * 1e6
    d1["Tdom_closed_noBR_g3044_keV"] = T_dom_closed(sig, 0.8, g3044) * 1e6
    d1["Tdom_closed_noBR_g3383_keV"] = T_dom_closed(sig, 0.8, 3.383) * 1e6
    d1["Tdom_closed_2AsigmaH_BR_g3044_keV"] = T_dom_closed(sig, 0.8, g3044) * 1e6 * 2**0.25
    # wrong-g sanity: what would one get with g = 10.75 or g = 2
    d1["Tdom_closed_noBR_g10.75_keV(wrong)"] = T_dom_closed(sig, 0.8, 10.75) * 1e6
    d1["Tdom_num_Tdec2MeV_BR_keV"] = T_dom_numeric(sig, sm=SM_dec2)[0] * 1e6
    d1["g_rho0_Tdec2MeV"] = SM_dec2["grho0"]
    for label, sm in (("Neff3", SM), ("Neff3044", SM44)):
        for br in (True, False):
            Td, td, _ = T_dom_numeric(sig, backreact=br, sm=sm)
            d1[f"Tdom_num_{label}_{'BR' if br else 'noBR'}_keV"] = Td * 1e6
            d1[f"tdom_num_{label}_{'BR' if br else 'noBR'}_s"] = td * HBAR_GeVs
            d1[f"zdom_num_{label}_{'BR' if br else 'noBR'}"] = math.exp(-ln_a(Td, sm)) - 1
        Td, td, _ = T_dom_numeric(sig, conv="H", backreact=True, sm=sm)
        d1[f"Tdom_num_{label}_2AsigmaH_BR_keV"] = Td * 1e6
    # A variants (Neff 3.044)
    for A in (0.7, 0.9):
        d1[f"Tdom_num_Neff3044_BR_A{A}_keV"] = T_dom_numeric(sig, A=A, sm=SM44)[0] * 1e6
        d1[f"Tdom_num_Neff3044_noBR_A{A}_keV"] = T_dom_numeric(sig, A=A, backreact=False, sm=SM44)[0] * 1e6
        d1[f"Tdom_num_Neff3044_2AsigmaH_BR_A{A}_keV"] = T_dom_numeric(sig, A=A, conv="H", sm=SM44)[0] * 1e6
    # Saikawa eq. 2.19 t_dom = 3 MP^2/(4 A sigma)
    d1["tdom_Saikawa219_s"] = 3 * MP**2 / (4 * 0.8 * sig) * HBAR_GeVs
    # convergence
    d1["Tdom_rtol1e-8_vs_1e-12_rel"] = abs(T_dom_numeric(sig, sm=SM44, rtol=1e-8)[0] /
                                           T_dom_numeric(sig, sm=SM44, rtol=1e-12)[0] - 1)
    res["1_Tdom"] = d1
    print("(1)", json.dumps(d1, indent=1))

    # ---------------- (2) R_H
    d2 = {}
    H70 = 70 * KMSMPC_GeV
    for label, sm in (("Neff3", SM), ("Neff3044", SM44)):
        for n in (1, 2):
            for f in (1.0,):
                r = evolve(sig, switch=("rad", f), n_post=n, sm=sm)
                H0p = r["H_end"]
                d2[f"{label}_n{n}_rho_w0_GeV4"] = r["rho_w_end"]
                d2[f"{label}_n{n}_H0_kmsMpc"] = H0p / KMSMPC_GeV
                d2[f"{label}_n{n}_RH_70"] = H0p / H70
                d2[f"{label}_n{n}_RH_67"] = H0p / (67 * KMSMPC_GeV)
                d2[f"{label}_n{n}_RH_74"] = H0p / (74 * KMSMPC_GeV)
                d2[f"{label}_n{n}_age_yr"] = r["t_end"] * HBAR_GeVs / YR_s
    # analytic n=1: rho_w0 = rho_r(Tdom) (T0/Tdom)^n
    Td = d1["Tdom_num_Neff3044_BR_keV"] * 1e-6
    for n in (1, 2):
        rw0 = rho_r(Td, SM44) * (T0 / Td) ** n
        d2[f"analytic_n{n}_rho_w0"] = rw0
        d2[f"analytic_n{n}_RH_70"] = math.sqrt((rw0 + rho_r(T0, SM44) + OMEGA_M_H2 * RHO_C100) / 3) / MP / H70
    # f_switch and A/convention scan, Neff3044
    scan = {}
    for n in (1, 2):
        vals = []
        for A in (0.7, 0.8, 0.9):
            for conv in ("t", "H"):
                for f in (0.1, 0.3, 1.0, 3.0):
                    r = evolve(sig, A=A, conv=conv, switch=("rad", f), n_post=n, sm=SM44)
                    vals.append(r["H_end"] / H70)
        scan[f"n{n}_RH70_min"] = min(vals)
        scan[f"n{n}_RH70_max"] = max(vals)
    d2["scan_A_conv_fswitch"] = scan
    # wall fraction at recombination T_rec = T0 (1 + 1089.92)
    Trec = T0 * (1 + 1089.92)
    for n in (1, 2):
        r = evolve(sig, switch=("rad", 1.0), n_post=n, sm=SM44, T_end=Trec)
        tot = r["rho_w_end"] + r["rho_r_end"] + r["rho_m_end"]
        d2[f"n{n}_1-Omega_w_Trec"] = (r["rho_r_end"] + r["rho_m_end"]) / tot
        d2[f"n{n}_rho_w/rho_r_Trec"] = r["rho_w_end"] / r["rho_r_end"]
    res["2_RH"] = d2
    print("(2)", json.dumps(d2, indent=1))

    # ---------------- (3) stress tests
    d3 = {}
    a = (7 / 8) * (4 / 11) ** (4 / 3)
    for label, Neff in (("Neff3", 3.0), ("Neff3044", 3.044)):
        # all wall energy at T_dom (rho_w = rho_r) into dark radiation: rho_dr/rho_1nu at T_rec
        d3[f"DeltaNeff_rad_{label}"] = (1 + a * Neff) / a
    # into stable matter: Omega h^2 = rho_m0 / rho_c100, rho_m0 = rho_r(Tdom) (T0/Tdom)^3 (g_s const)
    for label, sm in (("Neff3", SM), ("Neff3044", SM44)):
        Tdx = T_dom_numeric(sig, sm=sm)[0]
        d3[f"Omega_h2_matter_{label}"] = rho_r(Tdx, sm) * (T0 / Tdx) ** 3 / RHO_C100
        # via rho/s
        s_dom = (2 * math.pi**2 / 45) * math.exp(sm["ln_gs"](math.log(Tdx))) * Tdx**3
        s0 = (2 * math.pi**2 / 45) * sm["gs0"] * T0**3
        d3[f"Omega_h2_matter_via_rho_over_s_{label}"] = rho_r(Tdx, sm) / s_dom * s0 / RHO_C100
        d3[f"rho_over_s_dom_GeV_{label}"] = rho_r(Tdx, sm) / s_dom
    d3["rho_over_s_bound_Oh2_0.12_GeV"] = 0.12 * RHO_C100 / ((2 * math.pi**2 / 45) * SM44["gs0"] * T0**3)
    res["3_stress"] = d3
    print("(3)", json.dumps(d3, indent=1))

    # ---------------- (4) LMS unit conversion
    d4 = {}
    GeV3_Jm2 = GeV_J / HBARC_GeVm**2
    GeV3_kgm2 = GeV3_Jm2 / C_ms**2
    d4["1GeV3_in_J_per_m2"] = GeV3_Jm2
    d4["1GeV3_in_kg_per_m2"] = GeV3_kgm2
    d4["1kg_per_m2_in_GeV3"] = 1 / GeV3_kgm2
    d4["benchmark_sigma_kg_per_m2"] = sig * GeV3_kgm2
    t0_s = 13.797e9 * YR_s
    for lab, val, label_MeV in (("RM", 4.22e-9, 0.96), ("RML", 3.85e-9, 0.93)):
        for unit, fac in (("kg/m2_literal", 1.0), ("kg/cm2", 1e4), ("g/cm2", 10.0)):
            s_kgm2 = val * fac
            s_GeV3 = s_kgm2 / GeV3_kgm2
            key = f"{lab}_{unit}"
            d4[key + "_sigma_GeV3"] = s_GeV3
            d4[key + "_sigma13_MeV"] = s_GeV3 ** (1 / 3) * 1e3
            d4[key + "_ratio_label_over_sigma13"] = label_MeV / (s_GeV3 ** (1 / 3) * 1e3)
            d4[key + "_G_sigma_t0_over_c"] = G_SI * s_kgm2 * t0_s / C_ms
            d4[key + "_6piG_sigma_t0_over_c"] = 6 * math.pi * G_SI * s_kgm2 * t0_s / C_ms
            d4[key + "_vmax_MeV_lambda0.01"] = (s_GeV3 / (2 * math.sqrt(2 * 0.01) / 3)) ** (1 / 3) * 1e3
        # what sigma would the label imply if sigma = (label)^3
        s_lab = (label_MeV * 1e-3) ** 3
        d4[f"{lab}_label_cubed_kg_per_m2"] = s_lab * GeV3_kgm2
        d4[f"{lab}_label_cubed_over_printed"] = s_lab * GeV3_kgm2 / val
        # alternative: label = vev eta with sigma = (2 sqrt(2 lam)/3) eta^3 ; lambda needed for literal value
        s_lit = val / GeV3_kgm2
        lam_needed = (s_lit / ((2 * math.sqrt(2) / 3) * (label_MeV * 1e-3) ** 3)) ** 2
        d4[f"{lab}_lambda_needed_if_label_is_vev_literal"] = lam_needed
    # cross-check G sigma t0 in natural units for sigma = (1 MeV)^3
    d4["G_sigma_t0_for_1MeV3"] = (1e-9) * (t0_s / HBAR_GeVs) / (8 * math.pi * MP**2)
    d4["G_SI_from_MP"] = HBARC_GeVm * C_ms**2 / (8 * math.pi * (MP * GeV_J / C_ms**2) ** 2) * 1.0  # check
    # G = hbar c / (8 pi M_P^2) with M_P in kg:  M_P[kg] = MP*GeV_J/c^2
    MPkg = MP * GeV_J / C_ms**2
    d4["G_SI_from_MP"] = (HBAR_GeVs * GeV_J) * C_ms / (8 * math.pi * MPkg**2)
    res["4_LMS"] = d4
    print("(4)", json.dumps(d4, indent=1))

    # ---------------- (5) v_max
    d5 = {}
    kappa = 2 * math.sqrt(2 * 0.01) / 3

    def vmax_from_sigma(s):
        return (s / kappa) ** (1 / 3)

    def resid(lns, A, rhoL, conv="t", frozen_after_acc=False):
        s = math.exp(lns)
        if not frozen_after_acc:
            r = evolve(s, A=A, conv=conv, rhoL=rhoL, sm=SM44)
            return math.log(r["rho_w_end"]) - math.log(r["rho_r_end"] + r["rho_m_end"] + rhoL)
        raise NotImplementedError

    rhoLP = OMEGA_L_H2_PLANCK * RHO_C100
    for lab, A, rhoL in (("A0.8_PlanckL", 0.8, rhoLP), ("A0.643_PlanckL", 0.643, rhoLP),
                         ("A0.7_PlanckL", 0.7, rhoLP), ("A0.9_PlanckL", 0.9, rhoLP),
                         ("A0.8_L0", 0.8, 0.0)):
        lns = brentq(resid, math.log(1e-7), math.log(1e-3), args=(A, rhoL), xtol=1e-10)
        s = math.exp(lns)
        d5[f"{lab}_sigma_max_GeV3"] = s
        d5[f"{lab}_vmax_MeV"] = vmax_from_sigma(s) * 1e3
    # analytic no-backreaction: sigma = 3 MP^2 H0^2 t0 / A, Planck H0=67.36, t0=13.797 Gyr
    H0P = 67.36 * KMSMPC_GeV
    sa = 3 * MP**2 * H0P**2 * (13.797e9 * YR_s / HBAR_GeVs) / 0.8
    d5["analytic_noBR_sigma_max"] = sa
    d5["analytic_noBR_vmax_MeV"] = vmax_from_sigma(sa) * 1e3
    # CMB (reconciled LMS R+M+Lambda) and literal
    d5["vmax_CMB_reconciled_RML_MeV"] = d4["RML_kg/cm2_vmax_MeV_lambda0.01"]
    d5["vmax_CMB_literal_RML_MeV"] = d4["RML_kg/m2_literal_vmax_MeV_lambda0.01"]
    d5["vmax_CMB_reconciled_RM_MeV"] = d4["RM_kg/cm2_vmax_MeV_lambda0.01"]
    res["5_vmax"] = d5
    print("(5)", json.dumps(d5, indent=1))

    # ---------------- (6) BBN
    d6 = {}
    r = evolve(sig, sm=SM44, T_end=1e-5, dense=True)
    sol = r["sol"]
    rws = r["rho_w_scaling"]
    for TT in (1e-3, 0.7e-3, 1e-4, 7e-5, 5e-5, 3e-5):
        t = math.exp(sol.sol(math.log(TT))[0])
        rw, H, rr, rm = rws(TT, t)
        d6[f"T{TT*1e3:g}MeV_t_s"] = t * HBAR_GeVs
        d6[f"T{TT*1e3:g}MeV_Omega_w"] = rw / (rw + rr + rm)
        d6[f"T{TT*1e3:g}MeV_rho_w_over_rho_r"] = rw / rr
        d6[f"T{TT*1e3:g}MeV_DeltaNeff_equiv"] = rw / rho_1nu(TT, SM44)
        d6[f"T{TT*1e3:g}MeV_TnuT"] = float(SM44["TnuT"](math.log(TT)))
    # A range at 0.1 MeV and 1 MeV
    for A in (0.7, 0.9):
        r = evolve(sig, A=A, sm=SM44, T_end=1e-5, dense=True)
        for TT in (1e-3, 1e-4):
            t = math.exp(r["sol"].sol(math.log(TT))[0])
            rw, H, rr, rm = r["rho_w_scaling"](TT, t)
            d6[f"A{A}_T{TT*1e3:g}MeV_DeltaNeff_equiv"] = rw / rho_1nu(TT, SM44)
    res["6_BBN"] = d6
    print("(6)", json.dumps(d6, indent=1))

    # ---------------- (7) conservative-argument helpers
    d7 = {}
    # e+- residual energy fraction at 20 keV and 6 keV
    for TT in (2e-5, 6.4e-6):
        re, pe, se = e_pm_rho_p_s(TT)
        d7[f"rho_e/rho_r_at_{TT*1e6:g}keV"] = re / rho_r(TT, SM44)
    # maximal drop of rho_X a^4 relative factor for SM radiation from 100 GeV / 10 MeV to keV: g_rho/g_s^{4/3}
    d7["SM_rho_a4_gain_10.75_to_late"] = (SM44["grho0"] / SM44["gs0"] ** (4 / 3)) / (10.75 / 10.75 ** (4 / 3))
    d7["SM_rho_a4_gain_106.75_to_late"] = (SM44["grho0"] / SM44["gs0"] ** (4 / 3)) / (106.75 / 106.75 ** (4 / 3))
    # barrier energy density vs total density at T_dom
    d7["barrier_lambda_v4_over_4_GeV4"] = 0.01 * 1000.0**4 / 4
    d7["rho_tot_at_Tdom_GeV4"] = 2 * rho_r(d1["Tdom_num_Neff3044_BR_keV"] * 1e-6, SM44)
    # PBH Hawking lifetime ~ 5120 pi G^2 M^3/(hbar c^4): mass evaporating by recombination (t ~ 1.2e13 s)
    coef = 5120 * math.pi * G_SI**2 / ((HBAR_GeVs * GeV_J) * C_ms**4)
    d7["Hawking_coef_s_per_kg3"] = coef
    d7["M_evap_by_t_rec_kg"] = (1.2e13 / coef) ** (1 / 3)
    sig_kgm2 = sig * GeV3_kgm2
    d7["closed_wall_radius_with_that_mass_m"] = math.sqrt(d7["M_evap_by_t_rec_kg"] / (4 * math.pi * sig_kgm2))
    # graviton-mediated h decay estimate Gamma ~ m^5/MP^4 (order of magnitude)
    d7["h_grav_decay_lifetime_s_est"] = HBAR_GeVs / (141.42**5 / MP**4)
    res["7_conservative"] = d7
    print("(7)", json.dumps(d7, indent=1))

    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
