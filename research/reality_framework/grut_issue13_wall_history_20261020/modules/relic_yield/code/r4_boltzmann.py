"""Item 2-3: symmetric relic abundance of a secluded hidden bath with T_h = xi T_SM (Boltzmann equations),
xi_max for Omega_sym h^2 < 0.12 (retrospective known datum) and the minimum visible temperature of a
Z2-restoring thermal history.

Model (stated approximations):
  * species X (X + Xbar, g = 2, common number chemical potential), h (g = 1), D (g = 1);
  * one common kinetic temperature T_h (kinetic equilibrium assumed), Bose-Einstein n, rho, P;
  * collision terms: every tree-level 2->2 conversion and 3->2 process open at threshold (rates.json),
    detailed balance through chemical potentials, MB collision kernels;
  * the hidden sector exchanges no energy with the SM (only gravity: negligible, hidden_heating module);
  * SM radiation with constant g* = g*s = 106.75; H includes the hidden energy density;
  * initial state at T_h = 4 m_X (x = 0.25) with mu = 0 (thermal reference); xi = T_h/T_SM there.
Writes results/boltzmann.json.
"""
import json
import math
import os
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from common import (GSTAR, M_D, M_H, M_X, MP, OMEGA_DM_KNOWN, OMEGA_PER_YIELD, RESULTS, RHO_SM_COEF, S_SM_COEF,
                    TC_LEADING, TC_RESUMMED, V, XI_CHEM_24, XI_KIN_PHI, ZETA3, be_thermo_psi, dump, log_n0_mb)

T0 = time.time()
RATES = json.load(open(os.path.join(RESULTS, "rates.json")))
SPECIES = ["X", "h", "D"]
G = np.array([2.0, 1.0, 1.0])
MS = np.array([M_X, M_H, M_D])
IDX = {"X": 0, "Xb": 0, "h": 1, "D": 2}
PM = {"X": M_X, "Xb": M_X, "h": M_H, "D": M_D}
X_INIT = 0.25
X_END = 2000.0


def _build_processes(include_D=True):
    procs = []
    for t in RATES["three_to_two"]:
        if not include_D and ("D" in t["initial"] or "D" in t["final"]):
            continue
        xs = np.array(t["x_grid"])
        sv = np.array(t["sigma_v2_thermal"])
        procs.append(dict(kind="32", ini=t["initial"], fin=t["final"], xs=np.log(xs), lsv=np.log(np.maximum(sv, 1e-300)),
                          sv_thr=t["sigma_v2_thr"], x_last=xs[-1], sv_last=sv[-1], S_i=t["S_i"]))
    for c in RATES["conversions_2to2"]:
        if not include_D and ("D" in c["initial"] or "D" in c["final"]):
            continue
        procs.append(dict(kind="22", ini=c["initial"], fin=c["final"], xg=np.array(c["x_grid"]), lR=np.array(c["log_R_eq"])))
    for p in procs:
        p["nu"] = np.zeros(3)
        for t in p["fin"]:
            p["nu"][IDX[t]] += 1
        for t in p["ini"]:
            p["nu"][IDX[t]] -= 1
        p["iA"] = [IDX[t] for t in p["ini"]]
        p["iB"] = [IDX[t] for t in p["fin"]]
    return procs


def _log_neq(T, m):
    return log_n0_mb(T, m) - m / T


def log_Req(p, T):
    x = M_X / T
    if p["kind"] == "22":
        xg, lR = p["xg"], p["lR"]
        if x <= xg[-1]:
            return float(np.interp(x, xg, lR))
        return float(lR[-1] + (lR[-1] - lR[-2]) / (xg[-1] - xg[-2]) * (x - xg[-1]))
    # 3->2: thermal <sigma v^2>(x) interpolated (log-log); beyond x_last, linear in 1/x towards threshold value
    if x <= p["x_last"]:
        lsv = float(np.interp(math.log(max(x, X_INIT)), p["xs"], p["lsv"]))
    else:
        f = p["x_last"] / x
        lsv = math.log(max(p["sv_thr"] + (p["sv_last"] - p["sv_thr"]) * f, 1e-300))
    return lsv - math.log(p["S_i"]) + sum(_log_neq(T, PM[t]) for t in p["ini"])


class Hidden:
    def __init__(self, xi, include_D=True, collisions=True):
        self.xi = xi
        self.include_D = include_D
        self.collisions = collisions
        self.procs = _build_processes(include_D) if collisions else []
        self.active = [0, 1, 2] if include_D else [0, 1]
        self.Th_i = M_X / X_INIT
        self.Tsm_i = self.Th_i / xi
        self.diag = dict(max_G32_over_H=0.0, max_Gconv_over_H=0.0)

    def thermo(self, lnT, psi):
        T = math.exp(lnT)
        out = []
        for s in range(3):
            if s in self.active:
                out.append(be_thermo_psi(T, psi[s], MS[s], G[s]))
            else:
                out.append(dict(n=0.0, rho=0.0, P=0.0, dn_dpsi=1.0, drho_dpsi=0.0, dn_dlnT=0.0, drho_dlnT=0.0))
        return T, out

    def rhs(self, N, y, record=False):
        """State y = [ln T_h, psi_X, psi_h, psi_D], psi = (mu - m)/T; independent variable N = ln a."""
        lnT, psi = y[0], y[1:]
        T, th = self.thermo(lnT, psi)
        eta = [psi[s] + MS[s] / T for s in range(3)]
        rho_h = sum(t["rho"] for t in th)
        P_h = sum(t["P"] for t in th)
        Tsm = self.Tsm_i * math.exp(-N)
        H = math.sqrt((RHO_SM_COEF * Tsm ** 4 + rho_h) / 3.0) / MP
        C = np.zeros(3)
        g32 = gcv = 0.0
        ntot = sum(th[s]["n"] for s in self.active)
        if self.collisions:
            lne = [_log_neq(T, MS[s]) for s in range(3)]
            for p in self.procs:
                lR = log_Req(p, T)
                lf = lR + sum(eta[i] for i in p["iA"])
                lb = lR + sum(eta[i] for i in p["iB"])
                r = -math.exp(lf) * math.expm1(lb - lf) if lf > lb else math.exp(lb) * math.expm1(lf - lb)
                C += r * p["nu"]
                if record:
                    if p["kind"] == "32":
                        g32 += math.exp(lf)
                    else:
                        gcv += math.exp(lf)
        if record:
            return dict(T=T, H=H, Tsm=Tsm, G32_over_H=g32 / max(ntot, 1e-300) / H, Gconv_over_H=gcv / max(ntot, 1e-300) / H,
                        n=[t["n"] for t in th], rho=rho_h, P=P_h)
        dn = np.array([-3 * th[s]["n"] + C[s] / H for s in range(3)])
        drho = -3 * (rho_h + P_h)
        a = np.array([t["dn_dlnT"] for t in th])
        b = np.array([t["dn_dpsi"] for t in th])
        c = np.array([t["drho_dlnT"] for t in th])
        e = np.array([t["drho_dpsi"] for t in th])
        act = self.active
        num = drho - sum(e[s] * dn[s] / b[s] for s in act)
        den = sum(c[s] - e[s] * a[s] / b[s] for s in act)
        tau = num / den
        deta = np.zeros(3)
        for s in act:
            deta[s] = (dn[s] - a[s] * tau) / b[s]
        return np.concatenate([[tau], deta])

    def run(self, nsave=400):
        y0 = np.array([math.log(self.Th_i)] + [-m / self.Th_i for m in MS])

        def ev(N, y):
            return math.log(M_X / X_END) - y[0]
        ev.terminal = True
        sol = solve_ivp(self.rhs, (0.0, 80.0), y0, method="LSODA", rtol=1e-9, atol=1e-11, events=ev,
                        dense_output=True)
        self.sol = sol
        Nend = sol.t[-1]
        yend = sol.y[:, -1]
        T, th = self.thermo(yend[0], yend[1:])
        Tsm = self.Tsm_i * math.exp(-Nend)
        s_sm = S_SM_COEF * Tsm ** 3
        Y = [t["n"] / s_sm for t in th]
        Om = [MS[s] * Y[s] * OMEGA_PER_YIELD for s in range(3)]
        # diagnostics along the trajectory
        Ns = np.linspace(0, Nend, nsave)
        g32, gcv, xs, ratio = [], [], [], []
        for N in Ns:
            d = self.rhs(N, sol.sol(N), record=True)
            g32.append(d["G32_over_H"])
            gcv.append(d["Gconv_over_H"])
            xs.append(M_X / d["T"])
            ratio.append(d["T"] / d["Tsm"])
        g32 = np.array(g32)
        xs = np.array(xs)
        above = np.where(g32 > 1)[0]
        x_fo32 = float(xs[above[-1]]) if len(above) else None
        return dict(xi=self.xi, include_D=self.include_D, collisions=self.collisions, nfev=int(sol.nfev), status=int(sol.status),
                    N_end=float(Nend), x_end=float(M_X / T), Y=Y, Omega_h2_species=Om, Omega_h2=float(sum(Om)),
                    max_G32_over_H=float(g32.max()), max_Gconv_over_H=float(np.max(gcv)), x_at_last_G32_gt_1=x_fo32,
                    Th_over_Tsm_at_x=[[float(xv), float(rv)] for xv, rv in zip(xs[::40], ratio[::40])])


def hot_relic_omega(xi, masses_g):
    """Massless-BE mu=0 yield frozen: Y_s = 45 zeta3/(2 pi^4) g_s xi^3/g*s (derived)."""
    c = 45 * ZETA3 / (2 * math.pi ** 4) / GSTAR
    return sum(m * g for m, g in masses_g) * c * xi ** 3 * OMEGA_PER_YIELD


def entropy_check():
    """Collisionless run: comoving hidden entropy and particle numbers must be conserved."""
    hs = Hidden(1e-3, collisions=False)
    y0 = np.array([math.log(hs.Th_i)] + [-m / hs.Th_i for m in MS])
    sol = solve_ivp(hs.rhs, (0.0, 4.0), y0, method="LSODA", rtol=1e-10, atol=1e-12, dense_output=True)
    vals = []
    for N in (0.0, 1.0, 2.0, 4.0):
        y = sol.sol(N)
        T, th = hs.thermo(y[0], y[1:])
        S = sum((t["rho"] + t["P"] - (y[1 + i] + MS[i] / T) * T * t["n"]) / T for i, t in enumerate(th)) * math.exp(3 * N)
        Ntot = [t["n"] * math.exp(3 * N) for t in th]
        vals.append([S] + Ntot)
    vals = np.array(vals)
    return dict(rel_change_entropy=float(abs(vals[-1, 0] / vals[0, 0] - 1)),
                rel_change_numbers=[float(abs(vals[-1, i] / vals[0, i] - 1)) for i in (1, 2, 3)],
                x_at_N4=float(M_X / math.exp(sol.sol(4.0)[0])))


if __name__ == "__main__":
    out = dict(status="computed (Boltzmann ODE) + derived (hot-relic formula); thermal reference initial state is a "
                      "preparation input", omega_bound=OMEGA_DM_KNOWN)
    out["entropy_check_collisionless"] = entropy_check()
    print("entropy check", out["entropy_check_collisionless"], time.time() - T0)

    mg_all = [(M_X, 2), (M_H, 1), (M_D, 1)]
    mg_conv = [(M_X, 2), (M_X, 1), (M_D, 1)]
    mg_noD = [(M_X, 2), (M_X, 1)]
    c_hot = hot_relic_omega(1.0, mg_all)
    out["hot_relic"] = dict(
        formula="Omega h^2 = sum_s m_s g_s * [45 zeta(3)/(2 pi^4 g*s)] * xi^3 * s0/(rho_c/h^2)",
        coef_all_species_GeV=c_hot,
        coef_h_converted_to_X=hot_relic_omega(1.0, mg_conv),
        coef_no_D_h_converted=hot_relic_omega(1.0, mg_noD),
        xi_max_all=(OMEGA_DM_KNOWN / c_hot) ** (1 / 3),
        xi_max_h_converted=(OMEGA_DM_KNOWN / hot_relic_omega(1.0, mg_conv)) ** (1 / 3),
        xi_max_no_D=(OMEGA_DM_KNOWN / hot_relic_omega(1.0, mg_noD)) ** (1 / 3),
        Y_N_max_for_mass_100GeV=OMEGA_DM_KNOWN / (100.0 * OMEGA_PER_YIELD))

    scan = []
    for xi in np.logspace(-5, 0, 11):
        r = Hidden(float(xi)).run()
        r["Omega_hot_all"] = hot_relic_omega(xi, mg_all)
        hs0 = Hidden(float(xi), collisions=False)
        Ti, th0 = hs0.thermo(math.log(hs0.Th_i), [-m / hs0.Th_i for m in MS])
        r["Omega_frozen_initial_numbers"] = sum(MS[s] * th0[s]["n"] for s in range(3)) / (S_SM_COEF * hs0.Tsm_i ** 3) * OMEGA_PER_YIELD
        r["ratio_to_hot_relic"] = r["Omega_h2"] / r["Omega_hot_all"]
        scan.append(r)
        print("xi=%.3g Omega=%.4g ratio_hot=%.4f maxG32/H=%.3g maxGconv/H=%.3g nfev=%d t=%.0f" % (
            xi, r["Omega_h2"], r["ratio_to_hot_relic"], r["max_G32_over_H"], r["max_Gconv_over_H"], r["nfev"], time.time() - T0))
    out["scan_all_species"] = scan

    lg = np.log([r["max_G32_over_H"] for r in scan])
    lxs = np.log([r["xi"] for r in scan])
    out["xi_32_equilibrium"] = dict(value=float(np.exp(np.interp(0.0, lg, lxs))),
                                    meaning="xi above which the broken-phase 3->2 rate exceeds H somewhere on the trajectory "
                                            "(cannibal/chemical equilibrium); below it the hidden particle number is frozen")

    def logratio(lxi, include_D=True):
        return math.log(Hidden(10 ** lxi, include_D=include_D).run()["Omega_h2"] / OMEGA_DM_KNOWN)

    lx = brentq(logratio, -4.0, -2.5, xtol=1e-5)
    xi_max = 10 ** lx
    rmax = Hidden(xi_max).run()
    lxn = brentq(lambda l: logratio(l, False), -4.0, -2.5, xtol=1e-5)
    xi_max_noD = 10 ** lxn
    rmaxn = Hidden(xi_max_noD, include_D=False).run()
    out["xi_max"] = dict(all_species=xi_max, run_at_xi_max=rmax, no_D_population=xi_max_noD, run_no_D=rmaxn,
                         definition="xi = T_h/T_SM at T_h = 4 m_X with mu = 0 (thermal reference), g*s = 106.75",
                         scaling="Omega ~ v xi^3 at fixed couplings -> xi_max ~ (v/TeV)^(-1/3)")
    print("xi_max", xi_max, xi_max_noD, time.time() - T0)

    # item 3: restoring history
    rest = []
    for label, Tc in [("leading", TC_LEADING), ("resummed_low", TC_RESUMMED[0]), ("resummed_high", TC_RESUMMED[1])]:
        for xlabel, xm in [("all_species", xi_max), ("no_D", xi_max_noD), ("hot_relic_all", out["hot_relic"]["xi_max_all"])]:
            rest.append(dict(Tc_label=label, Tc_over_v=Tc, xi_label=xlabel, xi_max=xm,
                             T_SM_min_GeV=Tc * V / xm, H_at_restoration_GeV=math.sqrt(RHO_SM_COEF / 3) * (Tc * V / xm) ** 2 / MP))
    out["restoration"] = dict(
        rows=rest,
        window=dict(xi_kinetic_min_phi=XI_KIN_PHI, xi_max=xi_max, xi_chem_2to4=XI_CHEM_24[0],
                    kinetic_window_exists=XI_KIN_PHI < xi_max, chemical_equilibrium_compatible=XI_CHEM_24[1] < xi_max),
        scaling="T_SM,min = T_c/xi_max ~ (v/TeV)^(4/3); xi_kin ~ (v/TeV)^(1/2) (thermal_potential)")
    out["runtime_s"] = time.time() - T0
    dump("boltzmann.json", out)
    print("done", time.time() - T0)
