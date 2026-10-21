#!/usr/bin/env python3
"""Verifier v_wbound: homogeneous phi-only orbits by ODE period integration (independent of wbound/zero_d.py, which
uses theta-substitution quadrature).  Units v = 1, lambda = 0.01; energies reported in lambda v^4.

Method: start at the outer turning point phi = A, phidot = 0; integrate phi'' = -lam phi (phi^2-1) together with the
accumulators int U dt and int phi^2 dt until phidot changes sign the second time (one full period for over-barrier
orbits; for in-well orbits the period).  wbar = 1 - 2<U>/eps.  Also checks the virial identity
(3 wbar - 1) eps = lam (<phi^2> - 1) at the end of a full period (boundary term exactly zero there).

Also: (i) a small-|X| perturbation scan of the phi-X homogeneous system to test whether X can raise wbar above the
phi-only maximum at fixed energy; (ii) a homogeneous FRW run (prescribed radiation a ~ t^1/2) checking the exact budget
Delta(rho a^4) = -int H a^4 Phi_0 dt - [adot a^3 g0] + int addot a^3 g0 dt to integration accuracy.
Writes results/homog_check.json.
"""
import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

LAM = 0.01
G2, LX = 0.01, 0.001


def U(p):
    return LAM / 4 * (p * p - 1) ** 2


def period_avgs(A):
    """A = outer turning point (> 1 over-barrier start at +A swings to -A; 1 < A < sqrt2 in-well)."""
    eps = U(A)

    def rhs(t, y):
        p, pd, iu, ip2 = y
        return [pd, -LAM * p * (p * p - 1), U(p), p * p]

    # event: phidot crosses zero from below to above (return to the outer turning point)
    def ev(t, y):
        return y[1]
    ev.direction = 1
    ev.terminal = True
    # tiny kick inward to leave the turning point, then integrate to the next upward zero crossing of phidot
    sol = solve_ivp(rhs, [0, 1e6], [A, 0.0, 0.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev,
                    first_step=1e-6)
    # first event may be at t ~ 0; take the last
    te = sol.t_events[0]
    ye = sol.y_events[0]
    k = -1
    T = te[k]
    iu, ip2 = ye[k][2], ye[k][3]
    Ubar, p2 = iu / T, ip2 / T
    wbar = 1 - 2 * Ubar / eps
    return dict(A=A, eps_lam=eps / LAM, period=T, wbar=wbar, phi2bar=p2,
                identity_residual=(3 * wbar - 1) * eps - LAM * (p2 - 1), min_phi=float(sol.y[0].min()))


def A_of_eps(eps_lam):
    return math.sqrt(1 + math.sqrt(4 * eps_lam))


def part_quadrature():
    rows = []
    for e in (0.0025, 0.025, 0.125, 0.225, 0.2475, 0.25025, 0.2525, 0.275, 0.375, 0.5, 0.75, 1.0, 1.5, 2.0,
              2.127, 2.25, 2.5, 3.0, 5.0, 16.0, 25.0, 250.0, 2500.0):
        A = A_of_eps(e)
        rows.append(period_avgs(A))
    estar = brentq(lambda e: period_avgs(A_of_eps(e))["wbar"] - 1 / 3, 0.3, 2.0, xtol=1e-11)
    r = minimize_scalar(lambda le: -period_avgs(A_of_eps(math.exp(le)))["wbar"], bracket=(math.log(1.5), math.log(2.1), math.log(3.0)),
                        tol=1e-9)
    emax = math.exp(r.x)
    mx = period_avgs(A_of_eps(emax))
    return dict(rows=rows, eps_star_lam=estar, amp_star=A_of_eps(estar), eps_max_lam=emax, amp_max=A_of_eps(emax),
                wbar_max=mx["wbar"], max_identity_residual=max(abs(r["identity_residual"]) for r in rows))


# ---------- phi + complex X (homogeneous), fixed energy, wbar vs X share ----------
def U2(p, x1, x2):
    X2 = 0.5 * (x1 * x1 + x2 * x2)
    return LAM / 4 * (p * p - 1) ** 2 + G2 * p * p * X2 + LX * X2 * X2


def rhs2(t, y):
    p, x1, x2 = y[:3]
    X2 = 0.5 * (x1 * x1 + x2 * x2)
    fp = LAM * p * (p * p - 1) + 2 * G2 * p * X2
    fx = G2 * p * p + 2 * LX * X2
    return [y[3], y[4], y[5], -fp, -fx * x1, -fx * x2]


def wbar_phiX(eps_lam, xamp, charge_frac, T=4000.0, seed_phase=0.0):
    """phi starts at its turning point with share of energy; X amplitude xamp (real comp), rotating with fraction
    charge_frac of the circular-orbit angular velocity at phi = A (0 = no charge)."""
    eps = eps_lam * LAM
    # find phi start so that U2(phi, x, 0) + X kinetic = eps
    om = math.sqrt(max(G2 * 0 + 2 * LX * 0.5 * xamp**2, 1e-30))
    xd = charge_frac * xamp * math.sqrt(G2 * 1.0 + LX * xamp**2)
    kinX = 0.5 * xd**2
    f = lambda p: U2(p, xamp, 0.0) + kinX - eps
    if f(1.0) > 0:
        return None
    pa = brentq(f, 1.0, 50.0)
    y0 = [pa, xamp * math.cos(seed_phase), xamp * math.sin(seed_phase), 0.0, -xd * math.sin(seed_phase), xd * math.cos(seed_phase)]
    ts = np.linspace(0, T, 200001)
    sol = solve_ivp(rhs2, [0, T], y0, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=ts)
    Y = sol.y
    K = 0.5 * (Y[3:] ** 2).sum(axis=0)
    Uu = U2(Y[0], Y[1], Y[2])
    p = K - Uu
    rho = K + Uu
    w = np.trapezoid(p, ts) / np.trapezoid(rho, ts)
    p2 = np.trapezoid(Y[0] ** 2, ts) / T
    return dict(eps_lam=eps_lam, xamp=xamp, charge_frac=charge_frac, wbar=float(w), phi2bar=float(p2),
                rho_drift=float(np.ptp(rho) / rho.mean()))


def part_X():
    out = []
    for e in (2.127, 5.0, 16.0):
        for xa in (0.0, 0.01, 0.05, 0.1, 0.3, 0.6, 1.0):
            for cf in (0.0, 0.5, 1.0):
                if xa == 0.0 and cf > 0:
                    continue
                r = wbar_phiX(e, xa, cf)
                if r:
                    out.append(r)
    return out


# ---------- homogeneous FRW: exact budget closure ----------
def part_frw(A0=2.2, H0=0.02, tmax=None):
    """Prescribed radiation background a = (t/t0)^1/2 with H(t0) = H0 -> t0 = 1/(2 H0).  phi homogeneous."""
    t0 = 1 / (2 * H0)
    if tmax is None:
        tmax = t0 * 16  # a: 1 -> 4
    a = lambda t: math.sqrt(t / t0)
    Hf = lambda t: 0.5 / t
    adot = lambda t: a(t) * Hf(t)
    addot = lambda t: -0.25 * a(t) / t**2

    def rhs(t, y):
        p, pd, R_int, B_int = y
        H = Hf(t)
        rho = 0.5 * pd * pd + U(p)
        Phi0 = LAM * (p * p - 1)
        g0 = p * pd
        return [pd, -3 * H * pd - LAM * p * (p * p - 1), -H * a(t) ** 4 * Phi0, addot(t) * a(t) ** 3 * g0]

    sol = solve_ivp(rhs, [t0, tmax], [A0, 0.0, 0.0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    p, pd, Rint, Bint = sol.y[:, -1]
    rho_end = 0.5 * pd**2 + U(p)
    rho0 = U(A0)
    lhs = rho_end * a(tmax) ** 4 - rho0 * a(t0) ** 4
    bnd = -(adot(tmax) * a(tmax) ** 3 * p * pd - adot(t0) * a(t0) ** 3 * A0 * 0.0)
    rhs_ = Rint + bnd + Bint
    return dict(A0=A0, H0=H0, a_end=a(tmax), rho_a4_ratio=lhs / (rho0) + 1, closure_rel=abs(lhs - rhs_) / (rho0),
                boundary_rel=abs(bnd) / (rho_end * a(tmax) ** 4), addot_rel=abs(Bint) / (rho_end * a(tmax) ** 4),
                main_term_rel=Rint / rho0)


if __name__ == "__main__":
    res = {"units": "v = 1, lambda = 0.01; energies eps_lam = eps/(lambda v^4)"}
    res["quadrature"] = part_quadrature()
    q = {k: v for k, v in res["quadrature"].items() if k != "rows"}
    print(json.dumps(q, indent=1))
    for r in res["quadrature"]["rows"]:
        print(f"eps={r['eps_lam']:9.5f} A={r['A']:.4f} wbar={r['wbar']:+.5f} minphi={r['min_phi']:+.3f} res={r['identity_residual']:.1e}")
    res["phiX"] = part_X()
    for r in res["phiX"]:
        print(f"phiX eps={r['eps_lam']} xamp={r['xamp']} cf={r['charge_frac']} wbar={r['wbar']:.5f} drift={r['rho_drift']:.1e}")
    res["frw"] = [part_frw(A0, H0) for A0 in (1.2, 2.2, 3.0) for H0 in (0.02, 0.005)]
    for r in res["frw"]:
        print("frw", json.dumps(r))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "homog_check.json"), "w") as fh:
        json.dump(res, fh, indent=1)
