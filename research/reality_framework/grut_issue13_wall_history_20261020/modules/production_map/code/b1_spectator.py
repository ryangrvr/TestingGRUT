"""production_map / b1: hidden relic from de Sitter spectator condensates for H_inf >> v.

Output: results/b1_spectator.json.  Deterministic (numpy default_rng seed 20261021).  ~2-4 min on 2 threads.

Chain [D unless noted]:
 1. End of inflation [P]: phi, q1, q2, D are light spectators.  Starobinsky-Yokoyama equilibrium
    P ~ exp(-8 pi^2 U/(3 H^4)) (quartic U; variant with the running xi Hubble masses and couplings at mu = H
    quoted from inflation_spectator/results/xi_running.json [Q]).
 2. Instant reheating into a radiation era [P].  With chi = a phi and conformal time, a'' = 0 in RD and U is a
    homogeneous quartic once |phi| >> v, so chi'' = -grad V4(chi) is autonomous and the conformal energy
    E_c = |chi'|^2/2 + V4(chi) is EXACTLY conserved.  The initial data chi = phi_0, chi' = a' phi_0 (+ slow roll)
    give  eps = rho_h/rho_SM = sum<phi_i^2>/(6 M_P^2) + <U>/(3 H^2 M_P^2) + O(slow roll).
    The frozen-field term phi_0^2/(6 M_P^2) dominates <U>/(3H^2M_P^2) by ~ 7/sqrt(lambda).
 3. Number: homogeneous multi-field conformal oscillation integrated (velocity Verlet), number estimator
    J = sum_i <chi_i'^2>/omega_i (omega_i from zero crossings; exact adiabatic invariant for one quartic field).
    Bracket for the final hidden number after fragmentation and kinetic equilibration:
      floor  : mu_N = 0 thermal number at the same energy (number-changing fully efficient, then frozen),
      upper  : J (condensate quanta conserved through fragmentation; no number-changing).
 4. Y_N = n/s_SM; threshold H_inf where Y_N = 4.373e-12 (strict).  Critic's estimate checked.
 5. Non-instant reheating (MD inflaton era, decay rate Gamma_I = gamma H_e) [P]: same fields in the MD+decay
    background, suppression factor R(gamma) of Y.
"""
import json
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from common import (v, gX2, lphi, lX, hD2, lD, mD0, MP, gstar, gh, zeta3, S_COEF, Y_MAX_STRICT,
                    Y_MAX_THERMAL_COMPOSITION, PKT, dump)
import os

T0 = time.time()
OUT = {"module": "production_map", "part": "B spectator-condensate relic for H_inf >> v",
       "labels": "[D] derived, [C] computed, [Q] quoted, [E] estimate, [P] preparation datum"}
rng = np.random.default_rng(20261021)

K_S = S_COEF * (90 / (math.pi ** 2 * gstar)) ** 0.75     # s a^3 = K_S (a' M_P)^(3/2) in RD   [D]
TRH_COEF = (90 / (math.pi ** 2 * gstar)) ** 0.25         # T_RH,inst = TRH_COEF sqrt(H M_P)
Y_TH_COEF = 45 * zeta3 * gh / (2 * math.pi ** 4 * gstar)  # Y = Y_TH_COEF xi^3  (mu = 0 massless)

# ------------------------------------------------------------------ couplings at mu = H [Q inflation_spectator]
xr = json.load(open(os.path.join(PKT, "modules", "inflation_spectator", "results", "xi_running.json")))
tabH = np.array([r["H_GeV"] for r in xr["table"]])


def couplings_at(H, variant):
    if variant == "minimal":
        return dict(lphi=lphi, lX=lX, gX2=gX2, hD2=hD2, lD=lD, m2phi=0.0, m2X=0.0, m2D=0.0)
    lH = math.log(H)
    def ip(key):
        return float(np.interp(lH, np.log(tabH), [r["muH"][key] for r in xr["table"]]))
    return dict(lphi=ip("lam_phi"), lX=ip("lam_X"), gX2=ip("gX2"), hD2=hD2, lD=lD,
                m2phi=ip("hubble_mass2_over_H2_phi") - lphi * v ** 2 / H ** 2,
                m2X=ip("hubble_mass2_over_H2_X"),
                m2D=ip("hubble_mass2_over_H2_D") + mD0 ** 2 / H ** 2)


def u_field(x, r, d, c):
    """U/H^4 with x = phi/H, r^2 = (q1^2+q2^2)/H^2, d = D/H."""
    return (c["lphi"] / 4 * x ** 4 + c["gX2"] / 2 * x ** 2 * r ** 2 + c["lX"] / 4 * r ** 4
            + c["hD2"] / 2 * x ** 2 * d ** 2 + c["lD"] / 4 * d ** 4
            + 0.5 * c["m2phi"] * x ** 2 + 0.5 * c["m2X"] * r ** 2 + 0.5 * c["m2D"] * d ** 2)


def u4_only(x, r, d, c):
    return (c["lphi"] / 4 * x ** 4 + c["gX2"] / 2 * x ** 2 * r ** 2 + c["lX"] / 4 * r ** 4
            + c["hD2"] / 2 * x ** 2 * d ** 2 + c["lD"] / 4 * d ** 4)


NX, NR, ND = 140, 220, 220
XG = np.linspace(0, 9.0, NX)
RG = np.linspace(0, 17.0, NR)
DG = np.linspace(0, 17.0, ND)


def sy_equilibrium(H, variant):
    c = couplings_at(H, variant)
    X, R, Dd = XG[:, None, None], RG[None, :, None], DG[None, None, :]
    u = u_field(X, R, Dd, c)
    w = np.exp(-(8 * math.pi ** 2 / 3) * (u - u.min())) * R      # 2D measure for (q1, q2): r dr dtheta
    wx = np.trapezoid(np.trapezoid(np.trapezoid(w, DG, axis=2), RG, axis=1), XG)
    def mom(f):
        return float(np.trapezoid(np.trapezoid(np.trapezoid(w * f, DG, axis=2), RG, axis=1), XG) / wx)
    res = dict(phi2=mom(X ** 2), r2=mom(R ** 2), D2=mom(Dd ** 2), U4=mom(u4_only(X, R, Dd, c)),
               U_total=mom(u), couplings=c)
    res["sum_chi2"] = res["phi2"] + res["r2"] + res["D2"]
    return res, w, c


# ------------------------------------------------------------------ 1. equilibrium moments
eq = {}
eq["minimal (mu0 couplings, xi = 0, exact quartic, H-independent)"], W_min, C_min = sy_equilibrium(1e10, "minimal")
for H in (1e8, 1e10, 1e11, 1e12, 1e13):
    eq[f"running xi and couplings at mu=H, H={H:g}"] = sy_equilibrium(H, "running")[0]
# analytic single-field checks [D]
chk = {"single quartic <x^2> = Gamma(3/4)/Gamma(1/4) sqrt(3/(2 pi^2 lambda)) [D]":
           {"lphi": math.gamma(0.75) / math.gamma(0.25) * math.sqrt(3 / (2 * math.pi ** 2 * lphi)),
            "lD": math.gamma(0.75) / math.gamma(0.25) * math.sqrt(3 / (2 * math.pi ** 2 * lD))},
       "X alone <r^2> = sqrt(3/(2 pi^3 lX)) = 6.955 [Q hidden_heating]": math.sqrt(3 / (2 * math.pi ** 3 * lX)),
       "inflation_spectator 2D FP at 1e13 GeV: <phi^2> = 0.97 H^2, <2|X|^2> = 8.0 H^2 [Q]": None,
       "hidden_heating <U_true> = c H^4, c = 0.042-0.094 [Q]": None}
OUT["SY equilibrium moments (units H^2, H^4) [C]"] = {k: {kk: vv for kk, vv in val.items()} for k, val in eq.items()}
OUT["checks"] = chk


# ------------------------------------------------------------------ 2-3. conformal dynamics
def grad_v4(x, r, d, c):
    gx = c["lphi"] * x ** 3 + c["gX2"] * x * r ** 2 + c["hD2"] * x * d ** 2
    gr = c["lX"] * r ** 3 + c["gX2"] * x ** 2 * r
    gd = c["lD"] * d ** 3 + c["hD2"] * x ** 2 * d
    return gx, gr, gd


def curv_v4(x, r, d, c):
    return (3 * c["lphi"] * x ** 2 + c["gX2"] * r ** 2 + c["hD2"] * d ** 2,
            3 * c["lX"] * r ** 2 + c["gX2"] * x ** 2,
            3 * c["lD"] * d ** 2 + c["hD2"] * x ** 2)


def sample_initial(w, n):
    p = w.ravel() / w.sum()
    idx = rng.choice(p.size, size=n, p=p)
    ix, ir, idd = np.unravel_index(idx, w.shape)
    dx, dr, dd = XG[1] - XG[0], RG[1] - RG[0], DG[1] - DG[0]
    x = (XG[ix] + (rng.random(n) - 0.5) * dx) * rng.choice([-1, 1], n)
    r = np.abs(RG[ir] + (rng.random(n) - 0.5) * dr) * rng.choice([-1, 1], n)   # signed coordinate on the X ray
    d = (DG[idd] + (rng.random(n) - 0.5) * dd) * rng.choice([-1, 1], n)
    return x, r, d


def run_conformal(x, r, d, vx, vr, vd, c, tau_max=3000.0, dt=0.02, t_avg0=1000.0, background=None):
    """Velocity Verlet for chi'' = -grad V4(chi) + (a''/a)(tau) chi.  background: callable tau -> a''/a or None."""
    nsteps = int(round(tau_max / dt))
    i0 = int(round(t_avg0 / dt))
    acc = lambda x, r, d, t: tuple(-g + (background(t) * q if background else 0.0)
                                   for g, q in zip(grad_v4(x, r, d, c), (x, r, d)))
    ax, ar, ad = acc(x, r, d, 0.0)
    sv = [np.zeros_like(x) for _ in range(3)]
    sc = [np.zeros_like(x) for _ in range(3)]
    nz = [np.zeros_like(x) for _ in range(3)]
    sE = np.zeros_like(x)
    prev = [np.sign(x), np.sign(r), np.sign(d)]
    cnt = 0
    for i in range(nsteps):
        t = i * dt
        vx += 0.5 * dt * ax; vr += 0.5 * dt * ar; vd += 0.5 * dt * ad
        x += dt * vx; r += dt * vr; d += dt * vd
        ax, ar, ad = acc(x, r, d, t + dt)
        vx += 0.5 * dt * ax; vr += 0.5 * dt * ar; vd += 0.5 * dt * ad
        if i >= i0:
            cnt += 1
            for k, (q, vq) in enumerate(((x, vx), (r, vr), (d, vd))):
                sv[k] += vq * vq
                s = np.sign(q)
                nz[k] += (s != prev[k]) & (s != 0)
                prev[k] = s
            cv = curv_v4(x, r, d, c)
            for k in range(3):
                sc[k] += cv[k]
            sE += 0.5 * (vx * vx + vr * vr + vd * vd) + (c["lphi"] / 4 * x ** 4 + c["gX2"] / 2 * x ** 2 * r ** 2
                                                         + c["lX"] / 4 * r ** 4 + c["hD2"] / 2 * x ** 2 * d ** 2
                                                         + c["lD"] / 4 * d ** 4)
        else:
            prev = [np.sign(x), np.sign(r), np.sign(d)]
    Tw = cnt * dt
    v2 = [s / cnt for s in sv]
    om_zc = [np.pi * np.maximum(n, 1) / Tw for n in nz]
    om_cv = [np.sqrt(s / cnt) for s in sc]
    J_zc = sum(a / b for a, b in zip(v2, om_zc))
    J_cv = sum(a / b for a, b in zip(v2, om_cv))
    return dict(J_zc=J_zc, J_cv=J_cv, E=sE / cnt, v2=v2, om_zc=om_zc, nz=nz)


# calibration: single quartic field (phi only)
xs = np.array([0.5, 1.0, 2.0, 4.0])
cal = run_conformal(xs.copy(), np.zeros(4), np.zeros(4), xs.copy(), np.zeros(4), np.zeros(4), C_min,
                    tau_max=2000.0, t_avg0=200.0)
Ecal = 0.5 * xs ** 2 + lphi / 4 * xs ** 4
amp = (4 * Ecal / lphi) ** 0.25
J_exact = (4 / 3) * Ecal / (0.847213 * math.sqrt(lphi) * amp)
OUT["calibration single quartic field [C]"] = dict(J_zc_over_exact=(cal["J_zc"] / J_exact).tolist(),
                                                   J_curv_over_exact=(cal["J_cv"] / J_exact).tolist(),
                                                   E_drift=(cal["E"] / Ecal - 1).tolist(),
                                                   omega_exact_formula="0.847213 sqrt(lambda) amplitude [D]")

NS = 1500
x0, r0, d0 = sample_initial(W_min, NS)
# slow-roll velocity at the end of inflation (H = 1 units): phi_dot = -dU/dphi / 3   [D]
gx, gr, gd = grad_v4(x0, r0, d0, C_min)
dyn = {}
for label, sr in (("no slow-roll velocity", 0.0), ("with slow-roll velocity", 1.0)):
    vx0, vr0, vd0 = x0 - sr * gx / 3, r0 - sr * gr / 3, d0 - sr * gd / 3
    t1 = time.time()
    res = run_conformal(x0.copy(), r0.copy(), d0.copy(), vx0.copy(), vr0.copy(), vd0.copy(), C_min)
    Ec0 = 0.5 * (vx0 ** 2 + vr0 ** 2 + vd0 ** 2) + u4_only(x0, r0, d0, C_min)
    dyn[label] = dict(mean_E_c=float(Ec0.mean()), mean_E_c_late=float(res["E"].mean()),
                      energy_conservation=float(abs(res["E"].mean() / Ec0.mean() - 1)),
                      mean_J_zc=float(res["J_zc"].mean()), mean_J_curv=float(res["J_cv"].mean()),
                      J_zc_sem=float(res["J_zc"].std() / math.sqrt(NS)),
                      share_J_by_field=[float((res["v2"][k] / res["om_zc"][k]).mean() / res["J_zc"].mean())
                                        for k in range(3)],
                      min_zero_crossings=[int(n.min()) for n in res["nz"]], runtime_s=time.time() - t1)
OUT["conformal dynamics, instant reheating, minimal couplings (units H_e) [C]"] = dyn
base = dyn["with slow-roll velocity"]

e_min = eq["minimal (mu0 couplings, xi = 0, exact quartic, H-independent)"]
eps_coef_analytic = e_min["sum_chi2"] / 6 + e_min["U4"] / 3          # eps = coef (H/M_P)^2  [D]
eps_coef = base["mean_E_c"] / 3
J = base["mean_J_zc"]
Ycoef_J = J / K_S                                                      # Y_J = Ycoef_J (H/M_P)^(3/2)
xiE_coef = (eps_coef * gstar / gh) ** 0.25                             # xi_E = xiE_coef (H/M_P)^(1/2)
Ycoef_floor = Y_TH_COEF * xiE_coef ** 3
n_eq_coef = (zeta3 / math.pi ** 2) * gh * (30 * base["mean_E_c"] / (math.pi ** 2 * gh)) ** 0.75
crit_coef = 0.07 / 3 * 0.75 * TRH_COEF / 0.08                        # critic: Y = crit_coef (H/M_P)^(3/2) [Q critic]


def Hthr(coef, Y=Y_MAX_STRICT):
    return MP * (Y / coef) ** (2 / 3)


inst = {
    "eps = coef (H/M_P)^2: analytic sum<chi^2>/6 + <U4>/3 [D]": eps_coef_analytic,
    "eps coef from sampled initial data incl. slow roll [C]": eps_coef,
    "critic's eps coef (c/3, c = 0.07) [Q critic]": 0.07 / 3,
    "frozen-field term / potential term": (e_min["sum_chi2"] / 6) / (e_min["U4"] / 3),
    "xi_E = coef (H/M_P)^(1/2) [C]": xiE_coef,
    "Y_floor = coef (H/M_P)^(3/2) (mu=0 number at the same energy) [C]": Ycoef_floor,
    "Y_J = coef (H/M_P)^(3/2) (condensate quanta conserved) [C]": Ycoef_J,
    "critic Y = coef (H/M_P)^(3/2), coef = (c/3)(3/4)(T_RH/sqrt(H M_P))/0.08 [Q critic, re-derived]": crit_coef,
    "overpopulation J/n_eq(rho) [C]": J / n_eq_coef,
    "Y_N at H = 1e13 GeV: floor, J, critic": [Ycoef_floor * (1e13 / MP) ** 1.5, Ycoef_J * (1e13 / MP) ** 1.5,
                                              crit_coef * (1e13 / MP) ** 1.5],
    "H_inf threshold (Y = 4.373e-12): floor, J, critic-formula [C]": [Hthr(Ycoef_floor), Hthr(Ycoef_J),
                                                                     Hthr(crit_coef)],
    "H_inf threshold with Y = 4.30e-12 (critic's bound): critic formula": Hthr(crit_coef, Y_MAX_THERMAL_COMPOSITION),
    "variance budget for the floor: sqrt(sum<chi^2>)_max = sqrt(6 eps_th) M_P with eps_th(Y_max) [D]":
        math.sqrt(6 * gh / gstar * (Y_MAX_STRICT / Y_TH_COEF) ** (4 / 3)) * MP,
}
OUT["instant reheating, minimal couplings [C]"] = inst

# ------------------------------------------------------------------ running-xi variant (moments only; dynamics scale)
run_rows = []
for k, val in eq.items():
    if not k.startswith("running"):
        continue
    H = float(k.split("H=")[1])
    ratio_chi2 = val["sum_chi2"] / e_min["sum_chi2"]
    # J scales as (amplitude)^(3/2) at fixed direction distribution: J ~ <chi^2>^(3/4)  [D, scaling]
    eps_c = val["sum_chi2"] / 6 + val["U4"] / 3      # R = 0 in RD: the xi R term does not enter E_c
    xiE = (eps_c * gstar / gh) ** 0.25
    run_rows.append(dict(H=H, sum_chi2=val["sum_chi2"], phi2=val["phi2"], r2=val["r2"], D2=val["D2"],
                         ratio_to_minimal=ratio_chi2,
                         Y_floor=Y_TH_COEF * xiE ** 3 * (H / MP) ** 1.5,
                         Y_J=Ycoef_J * ratio_chi2 ** 0.75 * (H / MP) ** 1.5))
OUT["running-xi variant (couplings and Hubble masses at mu=H) [C]"] = run_rows
lHr = np.log([r["H"] for r in run_rows])


def Y_running(H, which):
    rr = float(np.interp(math.log(H), lHr, [r["ratio_to_minimal"] for r in run_rows]))
    s2 = float(np.interp(math.log(H), lHr, [r["sum_chi2"] for r in run_rows]))
    u4 = float(np.interp(math.log(H), lHr, [eq[f"running xi and couplings at mu=H, H={r['H']:g}"]["U4"]
                                            for r in run_rows]))
    if which == "J":
        return Ycoef_J * rr ** 0.75 * (H / MP) ** 1.5
    return Y_TH_COEF * ((s2 / 6 + u4 / 3) * gstar / gh) ** 0.75 * (H / MP) ** 1.5


from scipy.optimize import brentq as _brentq
OUT["running-xi thresholds H_inf (Y = 4.373e-12) [C]"] = {
    w: math.exp(_brentq(lambda lh: math.log(Y_running(math.exp(lh), w) / Y_MAX_STRICT), math.log(1e8),
                        math.log(1e13))) for w in ("floor", "J")}

# finite duration from zero: Gaussian, variance N/(4 pi^2) per real field (free diffusion) [E variant]
fin = {}
for N in (60.0,):
    s2 = N / (4 * math.pi ** 2)
    xg, rg, dg = (rng.normal(0, math.sqrt(s2), 600) for _ in range(3))
    # X: two real components; the conformal motion stays on the initial ray, use signed radius
    rg = np.sign(rng.random(600) - 0.5) * np.sqrt(rg ** 2 + rng.normal(0, math.sqrt(s2), 600) ** 2)
    resN = run_conformal(xg.copy(), rg.copy(), dg.copy(), xg.copy(), rg.copy(), dg.copy(), C_min)
    EcN = 0.5 * (xg ** 2 + rg ** 2 + dg ** 2) + u4_only(xg, rg, dg, C_min)
    epsN = EcN.mean() / 3
    fin[f"N={N:g} e-folds from zero"] = dict(sum_chi2=4 * s2, eps_coef=epsN,
                                             Y_floor_coef=Y_TH_COEF * (epsN * gstar / gh) ** 0.75,
                                             Y_J_coef=float(resN["J_zc"].mean()) / K_S,
                                             H_thr_floor=Hthr(Y_TH_COEF * (epsN * gstar / gh) ** 0.75),
                                             H_thr_J=Hthr(float(resN["J_zc"].mean()) / K_S))
OUT["finite-duration variant [C, E]"] = fin

# ------------------------------------------------------------------ 5. non-instant reheating: MD + decay
NSM = 300
xm, rm, dm = x0[:NSM].copy(), r0[:NSM].copy(), d0[:NSM].copy()
gxm, grm, gdm = grad_v4(xm, rm, dm, C_min)
md_rows = []
for gam in (10.0, 1.0, 0.1, 1e-2, 1e-3, 1e-4):
    # background in conformal time: y = (ln a, rhoI, rhoR) with rho in units H_e^2 M_P^2
    def bg(t, y):
        a, rI, rR = math.exp(y[0]), y[1], y[2]
        H = math.sqrt(max(rI + rR, 0) / 3)
        return [a * H, a * (-3 * H * rI - gam * rI), a * (-4 * H * rR + gam * rI)]
    tau_end_bg = 4000.0
    sol = solve_ivp(bg, (0, tau_end_bg), [0.0, 3.0, 0.0], method="LSODA", rtol=1e-10, atol=1e-14,
                    dense_output=True)
    def appa(t):
        y = sol.sol(min(t, tau_end_bg))
        return math.exp(2 * y[0]) * y[1] / 6
    t1 = time.time()
    resm = run_conformal(xm.copy(), rm.copy(), dm.copy(), xm - gxm / 3, rm - grm / 3, dm - gdm / 3, C_min,
                         background=appa)
    # late entropy: continue background until inflaton energy negligible (in N variable)
    def bgN(N, y):
        rI, rR = math.exp(y[0]), math.exp(y[1])
        H = math.sqrt((rI + rR) / 3)
        return [-3 - gam / H, -4 + gam * rI / (H * rR)]
    yend = sol.y[:, -1]
    N0 = yend[0]
    rI0, rR0 = max(yend[1], 1e-300), yend[2]
    solN = solve_ivp(bgN, (N0, N0 + 60), [math.log(rI0), math.log(rR0)], method="LSODA", rtol=1e-10,
                     atol=1e-12)
    rRa4 = math.exp(solN.y[1, -1] + 4 * solN.t[-1])
    frac_I = math.exp(solN.y[0, -1] - solN.y[1, -1])
    Ycoef = float(resm["J_zc"].mean()) / (S_COEF * (30 / (math.pi ** 2 * gstar)) ** 0.75 * rRa4 ** 0.75)
    eps_late = float(resm["E"].mean()) / rRa4
    md_rows.append(dict(gamma=gam, T_RH_over_T_RH_inst=math.sqrt(gam), Y_J_coef=Ycoef,
                        eps_late_coef=eps_late,
                        Y_floor_coef=Y_TH_COEF * (eps_late * gstar / gh) ** 0.75,
                        rho_R_a4_late=rRa4, residual_inflaton_fraction=frac_I, runtime_s=time.time() - t1))
# normalise R with the same 300-sample subset run in exact RD
res_sub = run_conformal(xm.copy(), rm.copy(), dm.copy(), xm - gxm / 3, rm - grm / 3, dm - gdm / 3, C_min)
Ysub = float(res_sub["J_zc"].mean()) / K_S
Esub = float(res_sub["E"].mean()) / 3
for rw in md_rows:
    rw["R_J"] = rw["Y_J_coef"] / Ysub
    rw["R_floor"] = (rw["eps_late_coef"] / Esub) ** 0.75
OUT["non-instant reheating (MD + decay, Gamma_I = gamma H_e) [C]"] = dict(
    rows=md_rows, subset_instant_Y_J_coef=Ysub, subset_instant_eps_coef=Esub,
    note="T_RH = TRH_COEF sqrt(gamma H_e M_P); R = Y/Y_instant at the same H_e")

# small-gamma power law fit for R (last three points)
g_arr = np.array([r["gamma"] for r in md_rows])
RJ = np.array([r["R_J"] for r in md_rows])
RF = np.array([r["R_floor"] for r in md_rows])
sel = g_arr <= 1e-2
pJ = np.polyfit(np.log(g_arr[sel]), np.log(RJ[sel]), 1)
pF = np.polyfit(np.log(g_arr[sel]), np.log(RF[sel]), 1)
OUT["non-instant reheating (MD + decay, Gamma_I = gamma H_e) [C]"]["small-gamma fits R = A gamma^p"] = dict(
    J=dict(p=pJ[0], A=math.exp(pJ[1])), floor=dict(p=pF[0], A=math.exp(pF[1])))


def R_of_gamma(gam, which):
    if gam >= g_arr.max():
        return 1.0 if which == "J" else 1.0
    arr = RJ if which == "J" else RF
    if gam <= g_arr.min():
        p = pJ if which == "J" else pF
        return math.exp(p[1]) * gam ** p[0]
    return float(np.exp(np.interp(math.log(gam), np.log(g_arr[::-1]), np.log(arr[::-1]))))


def Y_spectator(H_e, T_RH, which):
    """Y_N of the spectator relic for end-of-inflation Hubble H_e and reheating temperature T_RH [C, interpolated]."""
    gam = (T_RH / TRH_COEF) ** 2 / (H_e * MP)
    coef = Ycoef_J if which == "J" else Ycoef_floor
    return coef * (H_e / MP) ** 1.5 * (R_of_gamma(gam, which) if gam < 10 else 1.0)


tab = []
for TRH in (1e7, 1e8, 1e9, 1e10, 1e11, 1e12):
    row = dict(T_RH=TRH)
    for which in ("floor", "J"):
        # largest H_e with Y <= Y_max (H_e >= Gamma_I/10 so that gamma <= 10 is meaningful)
        Hgrid = np.logspace(3, 16, 1301)
        ok = [H for H in Hgrid if (TRH / TRH_COEF) ** 2 / (H * MP) <= 10 and
              Y_spectator(H, TRH, which) <= Y_MAX_STRICT]
        row[f"H_e_max_{which}"] = max(ok) if ok else None
    tab.append(row)
OUT["spectator bound H_e,max(T_RH) [C]"] = tab
for which, coef, fit in (("J", Ycoef_J, pJ), ("floor", Ycoef_floor, pF)):
    # gamma << 1:  Y = coef A_R gamma^p (H/M_P)^1.5, gamma = (T_RH/TRH_COEF)^2/(H M_P); for p = 1/2:
    # Y = (coef A_R/TRH_COEF) H T_RH / M_P^2   [D, C]
    kk = coef * math.exp(fit[1]) / TRH_COEF
    OUT[f"MD-reheating law ({which}): Y = k H_e T_RH/M_P^2"] = dict(
        k=kk, fitted_p=fit[0], H_e_times_T_RH_max_GeV2=Y_MAX_STRICT * MP ** 2 / kk,
        validity="gamma = Gamma_I/H_e <~ 1e-2 (T_RH <~ 0.1 T_RH,instant); fit over gamma = 1e-4..1e-2, "
                 "extrapolated below with the derived p = 1/2")
OUT["Y_spectator function"] = "Y = coef (H_e/M_P)^1.5 R(gamma), gamma = (T_RH/TRH_COEF)^2/(H_e M_P); coefs and R above"
OUT["headline"] = {
    "eps coefficient (instant): eps = c_eps (H/M_P)^2": eps_coef,
    "Y_N floor/J coefficients (instant), Y = c (H/M_P)^1.5": [Ycoef_floor, Ycoef_J],
    "H_inf threshold for Y_N = 4.373e-12 (floor, J)": inst["H_inf threshold (Y = 4.373e-12): floor, J, critic-formula [C]"][:2],
    "critic formula coefficient and threshold": [crit_coef, Hthr(crit_coef)],
    "overpopulation J/n_eq": J / n_eq_coef,
    "running-xi thresholds (floor, J)": OUT["running-xi thresholds H_inf (Y = 4.373e-12) [C]"],
    "N=60-from-zero thresholds (floor, J)": [fin["N=60 e-folds from zero"]["H_thr_floor"],
                                            fin["N=60 e-folds from zero"]["H_thr_J"]],
    "MD law H_e*T_RH max GeV^2 (floor, J)": [OUT["MD-reheating law (floor): Y = k H_e T_RH/M_P^2"]["H_e_times_T_RH_max_GeV2"],
                                             OUT["MD-reheating law (J): Y = k H_e T_RH/M_P^2"]["H_e_times_T_RH_max_GeV2"]],
    "runtime_s": time.time() - T0,
}
dump("b1_spectator.json", OUT)
if __name__ == "__main__":
    print(json.dumps(OUT["headline"], indent=1, default=float))
    print(json.dumps(OUT["instant reheating, minimal couplings [C]"], indent=1, default=float))
    print(json.dumps(OUT["calibration single quartic field [C]"], indent=1, default=float))
