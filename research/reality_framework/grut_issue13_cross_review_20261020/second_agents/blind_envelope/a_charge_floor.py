#!/usr/bin/env python3
"""(a) Blind check of the classical charge-energy floor rho >= m* |n_Q| for GRUT Candidate 1.

Written without reading Codex's audit.py/report.py/DERIVATION.md or any Claude script.

Potential (metric -+++, canonical kinetic terms):
  U = (lp/4)(phi^2 - v^2)^2 + gX^2 phi^2 |X|^2 + lX |X|^4 + (1/2)(mD0^2 + hD^2 phi^2) D^2 + (lD/4) D^4

Derivation (x = |X| > 0, s = phi^2 >= 0, D-part >= 0 for mD0^2 >= 0, lD >= 0, so D = 0 minimises U/x^2):
  U/x^2 >= (lp/4)(s - v^2)^2 / x^2 + lX x^2 + gX^2 s
  min over x of A/x^2 + lX x^2 (A >= 0) is 2 sqrt(A lX), at x^4 = A/lX
  -> f(s) = gX^2 s + sqrt(lp lX) |v^2 - s|
  piecewise linear in s; on [0, v^2] slope gX^2 - sqrt(lp lX), on [v^2, inf) slope gX^2 + sqrt(lp lX) > 0
  -> min f = v^2 min(gX^2, sqrt(lp lX)) = m*^2,  m* = v min(gX, (lp lX)^(1/4)).
Energy density T_00 = |Xdot|^2 + |grad X|^2 + (1/2)(phidot^2 + |grad phi|^2 + Ddot^2 + |grad D|^2) + U
Charge density n_Q = j^0 = i(X* Xdot - X Xdot*) = -2 Im(X* Xdot), so |n_Q| <= 2 |X| |Xdot|.
  T_00 - m*|n_Q| >= |Xdot|^2 + m*^2 |X|^2 - m*|n_Q|  = |Xdot + i sgn(n_Q) m* X|^2  >= 0   (completing the square)
  (and U - m*^2|X|^2 >= 0 by the minimisation above).  Saturation: grad = 0, phidot = 0, D = 0,
  (phi, |X|) at the minimiser of U/|X|^2, Xdot = +-i m* X.
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(20261020)


def mstar(v, gX, lp, lX):
    return v * min(gX, (lp * lX) ** 0.25)


def U(phi, Xabs, D, v, gX, lp, lX, mD0, hD, lD):
    return (lp / 4) * (phi**2 - v**2) ** 2 + gX**2 * phi**2 * Xabs**2 + lX * Xabs**4 \
        + 0.5 * (mD0**2 + hD**2 * phi**2) * D**2 + (lD / 4) * D**4


def check_completing_square(n=200000):
    """Numerically confirm |Xdot|^2 + m^2|X|^2 - m|n_Q| == |Xdot + i sgn(n_Q) m X|^2 identically (n_Q = -2 Im(X* Xdot))."""
    X = rng.normal(size=n) + 1j * rng.normal(size=n)
    Xd = rng.normal(size=n) + 1j * rng.normal(size=n)
    m = rng.uniform(0.1, 10, size=n)
    nQ = -2 * np.imag(np.conj(X) * Xd)
    lhs = np.abs(Xd) ** 2 + m**2 * np.abs(X) ** 2 - m * np.abs(nQ)
    rhs = np.abs(Xd + 1j * np.sign(nQ) * m * X) ** 2
    return float(np.max(np.abs(lhs - rhs) / (1 + np.abs(lhs))))


def min_U_over_X2(v, gX, lp, lX, mD0, hD, lD, nstart=400):
    """Direct numerical minimisation of U/|X|^2 over (phi, |X|, D), many random starts, in units of v."""
    best = (np.inf, None)
    # dimensionless: phi = v p, |X| = v x, D = v d ; U/|X|^2 / v^2
    def F(z):
        p, lx, d = z
        x = np.exp(lx)
        return U(p, x, d, 1.0, gX, lp, lX, mD0 / v, hD, lD) / x**2
    for _ in range(nstart):
        z0 = [rng.uniform(-2, 2), rng.uniform(-6, 3), rng.uniform(-1, 1)]
        r = minimize(F, z0, method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000, maxfev=40000))
        if r.fun < best[0]:
            best = (r.fun, r.x)
    # brute 2D grid (D = 0) as a non-optimiser cross-check
    p = np.linspace(-2, 2, 2001)
    lx = np.linspace(-8, 4, 2001)
    P, LX = np.meshgrid(p, lx, indexing="ij")
    Xg = np.exp(LX)
    G = U(P, Xg, 0.0, 1.0, gX, lp, lX, 0.0, 0.0, 0.0) / Xg**2
    gi = np.unravel_index(np.argmin(G), G.shape)
    return best[0] * v**2, best[1], float(G[gi]) * v**2, (float(P[gi]), float(np.exp(LX[gi])))


def random_pointwise(v, gX, lp, lX, mD0, hD, lD, n=2_000_000, adversarial=True):
    """Sample field values/velocities/gradients; return min of (rho - m*|nQ|)/(m*|nQ|) and min rho/|nQ|."""
    ms = mstar(v, gX, lp, lX)
    out = {}
    # generic sampling with scales spanning the problem
    sc = v * 10 ** rng.uniform(-3, 1, size=n)
    phi = rng.normal(size=n) * sc
    X = (rng.normal(size=n) + 1j * rng.normal(size=n)) * sc
    D = rng.normal(size=n) * sc * 0.3
    Xd = (rng.normal(size=n) + 1j * rng.normal(size=n)) * sc * gX
    phid = rng.normal(size=n) * sc * gX
    Dd = rng.normal(size=n) * sc * gX
    gradX2 = np.sum(np.abs(rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))) ** 2, axis=1) * (sc * gX) ** 2 * rng.uniform(0, 1, size=n)
    gradphi2 = np.sum(rng.normal(size=(n, 3)) ** 2, axis=1) * (sc * gX) ** 2 * rng.uniform(0, 1, size=n)
    gradD2 = np.sum(rng.normal(size=(n, 3)) ** 2, axis=1) * (sc * gX) ** 2 * rng.uniform(0, 1, size=n)
    rho = np.abs(Xd) ** 2 + gradX2 + 0.5 * (phid**2 + gradphi2 + Dd**2 + gradD2) + U(phi, np.abs(X), D, v, gX, lp, lX, mD0, hD, lD)
    nQ = -2 * np.imag(np.conj(X) * Xd)
    ratio = rho / np.maximum(np.abs(nQ), 1e-300)
    out["generic_n"] = n
    out["generic_min_rho_over_absnQ_over_mstar"] = float(np.min(ratio) / ms)
    out["generic_violations"] = int(np.sum(rho < ms * np.abs(nQ) * (1 - 1e-12)))
    if adversarial:
        # sample near the saturating configuration
        s_star = 0.0 if gX**4 > lp * lX else v**2
        A = (lp / 4) * (s_star - v**2) ** 2
        x_star = (A / lX) ** 0.25 if A > 0 else v * 0.1
        m = n // 2
        eps = 10 ** rng.uniform(-6, -1, size=m)
        phi = np.sqrt(s_star) + rng.normal(size=m) * eps * v
        xa = x_star * (1 + rng.normal(size=m) * eps)
        th = rng.uniform(0, 2 * np.pi, size=m)
        X = xa * np.exp(1j * th)
        sgn = rng.choice([-1, 1], size=m)
        Xd = 1j * sgn * ms * X * (1 + rng.normal(size=m) * eps) + (rng.normal(size=m) + 1j * rng.normal(size=m)) * eps * ms * x_star
        D = rng.normal(size=m) * eps * v
        rho = np.abs(Xd) ** 2 + U(phi, np.abs(X), D, v, gX, lp, lX, mD0, hD, lD)
        nQ = -2 * np.imag(np.conj(X) * Xd)
        ratio = rho / np.abs(nQ)
        out["adversarial_n"] = m
        out["adversarial_min_rho_over_absnQ_over_mstar"] = float(np.min(ratio) / ms)
        out["adversarial_violations"] = int(np.sum(rho < ms * np.abs(nQ) * (1 - 1e-10)))
        out["saturation_point_phi_over_v"] = float(np.sqrt(s_star) / v)
        out["saturation_point_absX_over_v"] = float(x_star / v)
    # optimiser on rho/|nQ| over homogeneous configs (gradients only add positive terms)
    def R(z):
        p, lx, d, a, b = z  # Xdot = (a + i b) * v
        x = np.exp(lx) * v
        Xd = (a + 1j * b) * v
        rho = abs(Xd) ** 2 + U(p * v, x, d * v, v, gX, lp, lX, mD0, hD, lD)
        nQ = abs(-2 * np.imag(x * Xd))  # X real positive wlog (global phase)
        return rho / max(nQ, 1e-300) / ms
    best = np.inf
    for _ in range(200):
        z0 = [rng.uniform(-1.5, 1.5), rng.uniform(-5, 2), rng.uniform(-0.5, 0.5), rng.uniform(-0.2, 0.2), rng.uniform(-0.2, 0.2)]
        r = minimize(R, z0, method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-15, maxiter=40000, maxfev=80000))
        best = min(best, r.fun)
    out["optimizer_min_rho_over_absnQ_over_mstar"] = float(best)
    return out


def main():
    res = {}
    res["completing_square_identity_max_rel_err"] = check_completing_square()
    bench = dict(v=1000.0, gX=0.1, lp=0.01, lX=0.001)
    bench["mD0"] = bench["gX"] * bench["v"]
    bench.update(hD=0.01, lD=0.001)
    cases = {"benchmark": bench}
    # extra coupling sets, both classes, to test the min() formula
    for i in range(6):
        gX = 10 ** rng.uniform(-2, 0)
        lp = 10 ** rng.uniform(-3, 0)
        lX = 10 ** rng.uniform(-4, 0)
        cases[f"random{i}"] = dict(v=1000.0, gX=gX, lp=lp, lX=lX, mD0=gX * 1000.0 * rng.uniform(0.1, 2), hD=10 ** rng.uniform(-3, -1), lD=10 ** rng.uniform(-3, -1))
    for name, c in cases.items():
        ms = mstar(c["v"], c["gX"], c["lp"], c["lX"])
        mX = c["gX"] * c["v"]
        num, arg, grid, gridarg = min_U_over_X2(**c, nstart=150 if name != "benchmark" else 600)
        row = dict(params=c, binding_class=bool(c["gX"] ** 4 > c["lp"] * c["lX"]),
                   mstar_analytic_GeV=ms, mX_GeV=mX, mu=ms / mX,
                   sqrt_min_U_over_X2_optimizer_GeV=float(np.sqrt(num)),
                   optimizer_argmin_phi_over_v_absX_over_v_D_over_v=[float(arg[0]), float(np.exp(arg[1])), float(arg[2])],
                   sqrt_min_U_over_X2_grid_GeV=float(np.sqrt(grid)), grid_argmin=gridarg,
                   rel_diff_optimizer=float(np.sqrt(num) / ms - 1), rel_diff_grid=float(np.sqrt(grid) / ms - 1))
        row.update(random_pointwise(**c, n=2_000_000 if name == "benchmark" else 400_000))
        res[name] = row
    # premise check: mD0^2 < 0 (tachyonic bare D mass with positive m_D) breaks U >= 0 and the floor
    c = dict(bench)
    c["mD0"] = 1j * 5.0  # mD0^2 = -25 GeV^2 ; m_D^2 = -25 + (0.01*1000)^2 = 75 > 0
    phi0, X0 = 0.0, None
    s_star = 0.0
    A = (c["lp"] / 4) * c["v"] ** 4
    x_star = (A / c["lX"]) ** 0.25
    Dmin = np.sqrt(25.0 / c["lD"])
    UD = 0.5 * (-25.0) * Dmin**2 + c["lD"] / 4 * Dmin**4
    res["premise_mD0sq_negative_example"] = dict(mD0sq=-25.0, mD_vac_sq=-25 + (c["hD"] * c["v"]) ** 2,
                                                  U_D_at_phi0_min_GeV4=UD,
                                                  note="U_D<0 at phi=0 shifts U/|X|^2 below m*^2 for large enough D-domain; floor needs mD0^2>=0 (benchmark mD0=gX v satisfies it)")
    json.dump(res, sys.stdout, indent=1, default=float)
    print()


if __name__ == "__main__":
    main()
