"""Verifier v_residual_bias, part 2: route A integrals, coherent-remnant -> m0 conversion, xi scans, dilute gas,
Hubble-only counterfactual, route C noise.  Uses the verifier's own width (vrb_width.py), not module code.
Quoted inputs [Q]: thermal_potential kibble_zurek.json (xi_KZ, eps_freeze vs xi), T_c = 1557.5 GeV, t_G = 6.917e-4,
hard-particle rate 4.98e-7 T, xi_kin = 6.6e-5; bias_threshold eps_persist = 6e-24, eps_safe = 1.5e-18.
"""
import json
import math
import os
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import vrb_width as W

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
KZ = os.path.join(HERE, "..", "..", "modules", "thermal_potential", "results", "kibble_zurek.json")
MP = 2.435323e18
GSTAR = 106.75
TC = W.TC
EPS_G = 6.917e-4
LAM, V = W.lam, W.v
HARD = 4.981371739469606e-07


def hub(T):
    return math.sqrt(math.pi ** 2 * GSTAR / 90) * T * T / MP


def mphi(Th):
    return math.sqrt(max(W.PI_P * (Th * Th - TC * TC), 0.0))


# ---------------------------------------------------------------- width table (own)
def build_table():
    """full and nsv on a dense grid; classical (slow power-law tail) interpolated from width.json's 14-point table."""
    eps = np.concatenate([[EPS_G], np.geomspace(1e-3, 2e3, 46)])
    with open(os.path.join(OUT, "width.json")) as fh:
        wt = json.load(fh)["table"]
    ce = np.log([r["eps"] for r in wt])
    cg = np.log([r["classical"] for r in wt])
    rows = []
    for e in eps:
        T = TC * (1 + e)
        gcl = float(np.exp(np.interp(math.log(e), ce, cg)))
        rows.append((e, W.width_over_T(T, "full"), gcl, W.width_over_T(T, "nsv")))
    return np.array(rows)


class Interp:
    def __init__(self, tab, col):
        self.le = np.log(tab[:, 0])
        self.lg = np.log(tab[:, col])

    def __call__(self, e):
        return float(np.exp(np.interp(math.log(max(e, 1e-12)), self.le, self.lg)))


def D_int(xi, Tprod_SM, gfun, scale=1.0, eps_cut=EPS_G):
    Tend = TC * (1 + eps_cut) / xi
    if Tprod_SM <= Tend:
        return 0.0

    def f(lT):
        T = math.exp(lT)
        Th = xi * T
        G = scale * gfun(Th / TC - 1) * Th
        m = mphi(Th)
        r = G / 2 if G < 2 * m else G / 2 - math.sqrt(G * G / 4 - m * m)
        return r / hub(T)
    # dense near the exit where Gamma/T rises
    lo, hi = math.log(Tend), math.log(Tprod_SM)
    br = [lo] + [math.log(TC * (1 + e) / xi) for e in (2e-3, 1e-2, 0.1, 1.0) if lo < math.log(TC * (1 + e) / xi) < hi] + [hi]
    return sum(quad(f, a, b, limit=400, epsrel=1e-9)[0] for a, b in zip(br[:-1], br[1:]))


# ---------------------------------------------------------------- KZ and soft rms
with open(KZ) as fh:
    _kz = json.load(fh)
_kx = np.log([r["xi"] for r in _kz["rows"]])
_kL = np.log([r["z1_ising_ginzburg_matched"]["xi_KZ_GeV^-1"] for r in _kz["rows"]])
_kE = np.log([r["z1_ising_ginzburg_matched"]["eps_freeze"] for r in _kz["rows"]])


def kz(xi):
    l = math.log(xi)
    return math.exp(np.interp(l, _kx, _kL)), math.exp(np.interp(l, _kx, _kE))


# c_w by direct numerical integration: sigma^2 = T/(2 pi^2) Int k^2 e^{-k^2 R^2}/(k^2 + m^2) dk, R = 1/m
CW_NUM = quad(lambda x: x * x * math.exp(-x * x) / (1 + x * x), 0, np.inf)[0] / (2 * math.pi ** 2)


def coherent(xi, Tprod_SM, D, phi0=V):
    xk, ehat = kz(xi)
    mhat = 1 / xk
    Th0 = xi * Tprod_SM
    m0 = mphi(Th0)
    Tdec = TC * (1 + ehat)
    A_lin = phi0 * (Tdec / Th0) ** 1.5 * math.sqrt(m0 / mhat)
    sig = math.sqrt(CW_NUM * Tdec * mhat)
    l10r = math.log10(A_lin / sig) - D / math.log(10)
    l10m = math.log10(math.erf(10 ** l10r / math.sqrt(2))) if l10r > -3 else l10r + math.log10(math.sqrt(2 / math.pi))
    return dict(xi_KZ=xk, eps_hat=ehat, m_hat=mhat, sigma=sig, A_lin=A_lin, log10_r=l10r, log10_m0=l10m)


# ---------------------------------------------------------------- Hubble-only via the exact adiabatic invariant
def J_action(A, m2):
    """J = 4 Int_0^A sqrt(2 (V(A) - V(x))) dx, V = m2 x^2/2 + lam x^4/4 (m2 >= 0)."""
    VA = 0.5 * m2 * A * A + 0.25 * LAM * A ** 4
    f = lambda t: math.sqrt(max(2 * (VA - 0.5 * m2 * (A * t) ** 2 - 0.25 * LAM * (A * t) ** 4), 0.0)) * A
    return 4 * quad(f, 0, 1, limit=200)[0]


def hubble_only(xi=5e-4, Tprod_SM=1e8):
    Th0 = xi * Tprod_SM
    J0 = J_action(V, mphi(Th0) ** 2) * 1.0          # a0^3 = 1 ; a ~ 1/Th
    xk, ehat = kz(xi)
    Tdec = TC * (1 + ehat)
    target = J0 * (Th0 / Tdec) ** -3                # J(a) a^3 const -> J_dec = J0 (a0/adec)^3 = J0 (Tdec/Th0)^3
    target = J0 * (Tdec / Th0) ** 3
    m2dec = mphi(Tdec) ** 2
    A_dec = brentq(lambda A: J_action(A, m2dec) - target, 1e-6, 1e4)
    # onset of nonlinearity (lam A^2 = m^2) from the same invariant
    def gapf(lT):
        Th = math.exp(lT)
        m2 = mphi(Th) ** 2
        A = brentq(lambda A: J_action(A, m2) - J0 * (Th / Th0) ** 3, 1e-6, 1e5)
        return LAM * A * A - m2, A
    lo, hi = math.log(Tdec * 1.0000001), math.log(Th0)
    l1 = brentq(lambda l: gapf(l)[0], lo, hi)
    Tnl = math.exp(l1)
    A_nl = gapf(l1)[1]
    sig = math.sqrt(CW_NUM * Tdec / xk)
    # oscillation count between nonlinear onset and trapping (amplitude falls below v(T) below T_c), quartic period
    H = hub(TC / xi)
    # amplitude in the quartic regime ~ A_dec * (Tdec/Th); trapping when lam A^2 = Pi (Tc^2 - T^2)
    eps_trap = brentq(lambda e: LAM * (A_dec * (1 - e)) ** 2 - W.PI_P * TC ** 2 * (1 - (1 - e) ** 2), 1e-9, 0.5)
    om = 2 * math.pi * math.sqrt(LAM) * A_dec / 7.4163
    dt_before = math.log(Tnl / TC) / H
    dt_after = eps_trap / H
    N_osc = om * (dt_before + dt_after) / (2 * math.pi)
    return dict(T_nonlinear=Tnl, A_at_onset=A_nl, A_dec=A_dec, sigma=sig, r=A_dec / sig,
                m0_if_coherent=math.erf(A_dec / sig / math.sqrt(2)), eps_trap_below_Tc=eps_trap,
                quartic_omega_GeV=om, N_osc_onset_to_trap=N_osc,
                phase_spread_rad_for_dA_over_A_sigma=2 * math.pi * N_osc * sig / A_dec,
                lame_floquet_mu=0.036, resonance_efolds=0.036 * math.sqrt(LAM) * A_dec * (dt_before + dt_after))


# ---------------------------------------------------------------- dilute gas: direct numerical integral
CD = (18 * LAM ** 2 + 12 * W.g2 ** 2 + 6 * W.h2 ** 2) / (128 * math.pi)
CE = (3 * LAM + 2 * W.g2 + W.h2) / 4


def dilute_numeric(A):
    """D = Int (Gamma_d/2) dt, Gamma_d = (C/c_eff) Pi/m, Pi = A H, m^2 = Pi - lam v^2, RD t = 1/(2H).
    Integrate in x = ln H from H_exit = lam v^2/A upward (zero mode responding from early on)."""
    Hx = LAM * V * V / A
    f = lambda y: (CD / CE) * A * math.exp(y) / math.sqrt(max(A * math.exp(y) - LAM * V * V, 1e-300)) / 2 * \
        (1 / (2 * math.exp(y)))            # dt = dH/(2 H^2) -> dt = dy/(2H)
    # substitution near the exit singularity (integrable 1/sqrt)
    g = lambda w: f(math.log(Hx) + w * w) * 2 * w
    return quad(g, 0, math.sqrt(math.log(1e30)), limit=400)[0]


def main():
    t0 = time.time()
    out = {}
    tab = build_table()
    out["table_eps_full_cl_nsv"] = tab.tolist()
    gf, gc, gn = Interp(tab, 1), Interp(tab, 2), Interp(tab, 3)
    xi, Tp = 5e-4, 1e8
    Tex = TC / xi
    Ds = dict(full=D_int(xi, Tp, gf), classical=D_int(xi, Tp, gc), nsv=D_int(xi, Tp, gn),
              hard=D_int(xi, Tp, lambda e: HARD), full_half=D_int(xi, Tp, gf, 0.5), full_x2=D_int(xi, Tp, gf, 2.0),
              full_cut_at_eps1em2=D_int(xi, Tp, gf, eps_cut=1e-2))
    # crude analytic check: constant Gamma/T = 3.84e-5 (high-T value) all the way down
    Ds["analytic_const_GammaT"] = 3.8424e-5 * xi * MP / (2 * math.sqrt(math.pi ** 2 * GSTAR / 90)) * (1 / Tex - 1 / Tp)
    star = dict(A_GeV=W.PI_P * xi ** 2 * MP / math.sqrt(math.pi ** 2 * GSTAR / 90), A_from_lamv2_over_Hexit=LAM * V * V / hub(Tex),
                T_exit=Tex, H_exit=hub(Tex), Gamma_over_H_eps7em3=gf(7e-3) * TC * 1.007 / hub(Tex * 1.007), D=Ds)
    for k, D in Ds.items():
        star["m0_" + k] = coherent(xi, Tp, D)
    star["c_w_numeric"] = CW_NUM
    star["c_w_closed_form_module"] = (math.sqrt(math.pi) / 2 - (math.pi / 2) * math.e * math.erfc(1.0)) / (2 * math.pi ** 2)
    star["hubble_only"] = hubble_only(xi, Tp)
    # production_map WKB [Q formula] arithmetic
    Apm = LAM * V * V / hub(Tex)
    star["production_map_ratio"] = (1.7 * V * math.sqrt(LAM) * V / Apm) / (0.11 * math.sqrt(LAM) * V)
    star["dilute_formula_D_sameA"] = CD / (2 * CE) * math.pi * Apm / (2 * math.sqrt(LAM) * V)
    star["dilute_numeric_D_sameA"] = dilute_numeric(Apm)
    out["H_star"] = star
    print(json.dumps(star, indent=1, default=float), flush=True)

    # xi_min tables
    def lm0(x, f, scale=1.0):
        Tp_ = TC * f / x
        return coherent(x, Tp_, D_int(x, Tp_, gf, scale))["log10_m0"]
    xm = {}
    for f in (32.1, 10.0, 2.0, 1.1, 1.01):
        row = {}
        for tgt in (-17, -20, -30, math.log10(6e-24), math.log10(1.5e-18)):
            try:
                row["%.2f" % tgt] = brentq(lambda l: lm0(10 ** l, f) - tgt, -5.5, -2.0, xtol=1e-5)
                row["%.2f" % tgt] = 10 ** row["%.2f" % tgt]
            except ValueError:
                row["%.2f" % tgt] = None
        xm["f=%g" % f] = row
    out["xi_min"] = xm
    # band sensitivity at f = 32.1: Gamma x 1/2 and x 2
    out["xi_min_f32_Gamma_half"] = {"%.2f" % t: 10 ** brentq(lambda l: lm0(10 ** l, 32.1, 0.5) - t, -5.5, -2.0, xtol=1e-5)
                                    for t in (-20, math.log10(6e-24))}
    out["xi_min_f32_Gamma_x2"] = {"%.2f" % t: 10 ** brentq(lambda l: lm0(10 ** l, 32.1, 2.0) - t, -5.5, -2.0, xtol=1e-5)
                                  for t in (-20, math.log10(6e-24))}
    out["log10m0_at_xi_kin_6.6e-5"] = dict(f32=lm0(6.6e-5, 32.1), f32_half=lm0(6.6e-5, 32.1, 0.5),
                                          f10=lm0(6.6e-5, 10.0), f1p1=lm0(6.6e-5, 1.1))
    out["scan_selected"] = {"%g_f%g" % (x, f): lm0(x, f) for x in (1e-5, 3e-5, 4.5e-5, 1e-4, 2e-4, 1e-3)
                            for f in (32.1, 1.1, 1.01)}
    # dilute gas thresholds
    dd = {}
    for tgt in (-20,):
        for rH in (1.0, 100.0):
            Dn = math.log(rH * math.sqrt(2 / math.pi)) + abs(tgt) * math.log(10)
            A = brentq(lambda A: dilute_numeric(A) - Dn, 1e3, 1e12)
            Hx = LAM * V * V / A
            dd["m0=1e%d_rH=%g" % (tgt, rH)] = dict(A=A, T_SM_end=math.sqrt(Hx * MP / math.sqrt(math.pi ** 2 * GSTAR / 90)))
    for Tend in (2.6e7, 1e8, 1e9, 6e9):
        A = LAM * V * V / hub(Tend)
        dd["D_at_Tend_%g" % Tend] = dilute_numeric(A)
    dd["D_over_A_formula"] = CD / (2 * CE) * math.pi / (2 * math.sqrt(LAM) * V)
    dd["D_over_A_numeric_A1e7"] = dilute_numeric(1e7) / 1e7
    dd["C"], dd["c_eff"] = CD, CE
    out["dilute"] = dd
    # route C white-noise estimate: Hubble volumes in today's horizon at T = T_exit and KZ domains per volume
    T0 = 2.7255 * 8.617333e-14
    H0 = 0.674 * 2.1332e-42          # 100 h km/s/Mpc = 2.1332e-42 h GeV
    gs0 = 3.909
    a_ratio = Tex / T0 * (GSTAR / gs0) ** (1 / 3)      # a0/a
    NH = (hub(Tex) / (a_ratio * H0)) ** 3
    ndom = (1 / hub(Tex) / kz(xi)[0]) ** 3
    out["routeC"] = dict(N_Hubble=NH, dom_per_H=ndom, noise=(NH * ndom) ** -0.5)
    out["runtime_s"] = time.time() - t0
    with open(os.path.join(OUT, "damping.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ("table_eps_full_cl_nsv", "H_star")}, indent=1, default=float))


if __name__ == "__main__":
    main()
