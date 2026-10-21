#!/usr/bin/env python3
"""Track wbound, item (b): homogeneous (0D) orbits of the unchanged potential.  Units v = 1.

Part 1 (phi only, exact period quadrature).  For U = lam/4 (phi^2-1)^2 the period-averaged w depends only on
e = eps/(lam v^4/4).  Computes wbar(e), <phi^2>, the threshold e* where wbar = 1/3 (equivalently <phi^2> = v^2),
the maximum wbar and checks the virial identity 3pbar - rho = lam v^2 (<phi^2> - v^2).

Part 2 (phi, complex X, D, benchmark couplings, ODE).  Random homogeneous initial data at fixed energy; long-time
averages.  Checks the identity <3p - rho> = lam v^2 <phi^2 - 1> - m_D0^2 <D^2>, that orbits that keep one sign of
phi have wbar <= 1/3, and whether adding X (or D) raises the maximum wbar at fixed energy.
Writes results/zero_d.json.
"""
import json
import math
import os

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

LAM, G, LX, HD, LD, MD0 = 0.01, 0.1, 0.001, 0.01, 0.001, 0.1
BARRIER = LAM / 4


# ---------------- Part 1: phi only, exact quadrature ----------------
def U1(p):
    return LAM / 4 * (p * p - 1) ** 2


def orbit_avgs(e):
    """Period averages for energy eps = e*BARRIER.  Returns dict with period, <U>, <phi^2>, wbar."""
    eps = e * BARRIER
    if e > 1:   # over the barrier: symmetric orbit between -pm and pm
        pm = math.sqrt(1 + math.sqrt(4 * eps / LAM))
        c, d = 0.0, pm
    else:       # inside the +v well between pmin and pmax
        r = math.sqrt(4 * eps / LAM)
        pmax, pmin = math.sqrt(1 + r), math.sqrt(max(1 - r, 0.0))
        c, d = 0.5 * (pmax + pmin), 0.5 * (pmax - pmin)

    def dt(th):
        ph = c + d * math.sin(th)
        k = 2 * (eps - U1(ph))
        return d * math.cos(th) / math.sqrt(k) if k > 0 else 0.0

    opts = dict(limit=400, epsabs=0, epsrel=1e-11)
    pts = None
    if e > 1 and e < 1.5:     # slow passage near phi = 0 (theta = 0)
        pts = [0.0]
    T = quad(dt, -math.pi / 2, math.pi / 2, points=pts, **opts)[0]
    IU = quad(lambda th: dt(th) * U1(c + d * math.sin(th)), -math.pi / 2, math.pi / 2, points=pts, **opts)[0]
    IP2 = quad(lambda th: dt(th) * (c + d * math.sin(th)) ** 2, -math.pi / 2, math.pi / 2, points=pts, **opts)[0]
    Ubar, p2 = IU / T, IP2 / T
    wbar = 1 - 2 * Ubar / eps      # <p> = <K> - <U> = eps - 2<U>
    ident = (3 * wbar - 1) * eps - LAM * (p2 - 1)
    return dict(e=e, eps_v4=eps, amplitude=c + d, half_period=T, Ubar=Ubar, phi2bar=p2, wbar=wbar,
                identity_residual=ident)


def part1():
    es_in = [0.01, 0.1, 0.5, 0.9, 0.99]
    es_over = [1.001, 1.01, 1.1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 100, 300, 1e3, 1e4]
    rows = [orbit_avgs(e) for e in es_in + es_over]
    estar = brentq(lambda e: orbit_avgs(e)["wbar"] - 1 / 3, 1.01, 100, xtol=1e-12)
    st = orbit_avgs(estar)
    # maximum wbar over e > 1
    grid = np.geomspace(1.5, 1e4, 300)
    ws = np.array([orbit_avgs(e)["wbar"] for e in grid])
    i = int(np.argmax(ws))
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda le: -orbit_avgs(math.exp(le))["wbar"],
                        bracket=(math.log(grid[max(i - 1, 0)]), math.log(grid[i]), math.log(grid[min(i + 1, len(grid) - 1)])),
                        tol=1e-10)
    emax = math.exp(r.x)
    mx = orbit_avgs(emax)
    # asymptotics: pure quartic limit, wbar -> 1/3 + (4/3) c2 / A^2 ...
    return dict(rows=rows,
                threshold=dict(e_star=estar, eps_star_v4=estar * BARRIER, eps_star_over_lam_v4=estar / 4,
                               amplitude_star=st["amplitude"], phi2bar_at_star=st["phi2bar"],
                               rho_star_GeV4_at_v_1TeV=estar * BARRIER * 1e12),
                maximum=dict(e_max=emax, eps_max_v4=emax * BARRIER, amplitude=mx["amplitude"], wbar_max=mx["wbar"]),
                max_identity_residual=max(abs(rw["identity_residual"]) for rw in rows))


# ---------------- Part 2: phi, X, D homogeneous ODE ----------------
def U4(y):
    ph, x1, x2, D = y[0], y[1], y[2], y[3]
    X2 = 0.5 * (x1 * x1 + x2 * x2)
    return LAM / 4 * (ph * ph - 1) ** 2 + G**2 * ph * ph * X2 + LX * X2 * X2 \
        + 0.5 * (MD0**2 + HD**2 * ph * ph) * D * D + LD / 4 * D**4


def rhs(t, y):
    ph, x1, x2, D = y[:4]
    X2 = 0.5 * (x1 * x1 + x2 * x2)
    dph = LAM * ph * (ph * ph - 1) + 2 * G**2 * ph * X2 + HD**2 * ph * D * D
    dX2 = G**2 * ph * ph + 2 * LX * X2
    dD = (MD0**2 + HD**2 * ph * ph) * D + LD * D**3
    return np.concatenate([y[4:], [-dph, -dX2 * x1, -dX2 * x2, -dD]])


def run_state(y0, T=6000.0, n=120001):
    ts = np.linspace(0, T, n)
    sol = solve_ivp(rhs, [0, T], y0, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=ts)
    Y = sol.y
    K = 0.5 * (Y[4:] ** 2).sum(axis=0)
    Uu = U4(Y)
    rho = K + Uu
    p = K - Uu
    # trapezoid time averages
    av = lambda f: float(np.trapezoid(f, ts) / T)
    wbar = av(p) / av(rho)
    lhs = av(3 * p - rho)
    rhs_id = LAM * av(Y[0] ** 2 - 1) - MD0**2 * av(Y[3] ** 2)
    Gend = float(Y[0, -1] * Y[4, -1] + Y[1, -1] * Y[5, -1] + Y[2, -1] * Y[6, -1] + Y[3, -1] * Y[7, -1])
    G0 = float(y0[0] * y0[4] + y0[1] * y0[5] + y0[2] * y0[6] + y0[3] * y0[7])
    return dict(wbar=wbar, rho=av(rho), rho_drift=float(rho.std() / rho.mean()), lhs=lhs, rhs=rhs_id,
                boundary=(Gend - G0) / T, phi_min=float(Y[0].min()), phi_max=float(Y[0].max()),
                one_sign=bool(Y[0].min() > 0 or Y[0].max() < 0), absphi_le_1=bool(np.abs(Y[0]).max() <= 1),
                phi2bar=av(Y[0] ** 2), D2bar=av(Y[3] ** 2), X2bar=av(0.5 * (Y[1] ** 2 + Y[2] ** 2)))


def random_state(rng, eps, fX, fD, start_phi=None):
    """Random homogeneous state with total energy eps: field values random, then momenta rescaled."""
    for _ in range(1000):
        ph = start_phi if start_phi is not None else rng.uniform(-2.5, 2.5)
        x = rng.normal(size=2) * (rng.uniform(0, 3) if fX > 0 else 0)
        D = rng.normal() * (rng.uniform(0, 3) if fD > 0 else 0)
        y = np.array([ph, x[0], x[1], D, 0, 0, 0, 0], float)
        u = U4(y)
        if u >= eps:
            continue
        mom = rng.normal(size=4) * np.array([1, fX, fX, fD])
        if not np.any(mom):
            mom[0] = 1.0
        mom *= math.sqrt(2 * (eps - u)) / np.linalg.norm(mom)
        y[4:] = mom
        return y
    return None


def part2(seed=20261021, nper=16):
    rng = np.random.default_rng(seed)
    out = []
    for e in (0.5, 0.95, 1.5, 3.0, 6.0, 9.0, 20.0, 100.0):
        eps = e * BARRIER
        for label, fX, fD in (("phi_only", 0.0, 0.0), ("phi_X", 1.0, 0.0), ("phi_X_D", 1.0, 1.0)):
            for k in range(nper):
                y0 = random_state(rng, eps, fX, fD)
                if y0 is None:
                    continue
                r = run_state(y0)
                r.update(e=e, family=label, k=k)
                out.append(r)
    # summaries
    summ = {}
    for e in sorted({r["e"] for r in out}):
        for fam in ("phi_only", "phi_X", "phi_X_D"):
            rs = [r for r in out if r["e"] == e and r["family"] == fam]
            if not rs:
                continue
            summ[f"e={e}_{fam}"] = dict(n=len(rs), wbar_max=max(r["wbar"] for r in rs),
                                        wbar_min=min(r["wbar"] for r in rs),
                                        n_wbar_gt_third=sum(r["wbar"] > 1 / 3 for r in rs),
                                        n_one_sign=sum(r["one_sign"] for r in rs),
                                        max_wbar_one_sign=max([r["wbar"] for r in rs if r["one_sign"]], default=None),
                                        max_identity_err_rel=max(abs(r["lhs"] - r["rhs"] - r["boundary"]) / r["rho"] for r in rs),
                                        max_boundary_rel=max(abs(r["boundary"]) / r["rho"] for r in rs))
    viol_one_sign = [r for r in out if r["one_sign"] and r["wbar"] > 1 / 3 + 1e-6]
    viol_crit = [r for r in out if (r["wbar"] > 1 / 3) != (LAM * (r["phi2bar"] - 1) - MD0**2 * r["D2bar"] > 0)]
    return dict(summary=summ, n_states=len(out), one_sign_violations=len(viol_one_sign),
                criterion_mismatches=len(viol_crit),
                max_wbar_any=max(r["wbar"] for r in out),
                min_e_with_wbar_gt_third=min([r["e"] for r in out if r["wbar"] > 1 / 3], default=None),
                states=out)


if __name__ == "__main__":
    res = dict(units="v = 1, lambda_phi = 0.01, g_X = 0.1, lambda_X = 0.001, h_D = 0.01, lambda_D = 0.001, m_D0 = 0.1",
               part1_phi_only=part1())
    print(json.dumps({k: v for k, v in res["part1_phi_only"].items() if k != "rows"}, indent=1))
    res["part2_multifield"] = part2()
    p2 = {k: v for k, v in res["part2_multifield"].items() if k != "states"}
    print(json.dumps(p2, indent=1))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "results", "zero_d.json"), "w") as fh:
        json.dump(res, fh, indent=1)
