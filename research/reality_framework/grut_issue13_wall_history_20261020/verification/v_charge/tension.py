#!/usr/bin/env python3
"""Independent Q-matter / +v-vacuum planar interface tension at mu = m* (verifier v_charge).

Grand-potential density in the frame co-rotating at mu, with F = sqrt2 |X| (canonical real amplitude):
  omega = 1/2 phi'^2 + 1/2 F'^2 + V(phi,F),
  V = lp/4 (phi^2-1)^2 + g^2 phi^2 F^2/2 + lX F^4/4 - mu^2 F^2/2,   (v = 1)
since |dX|^2 = F'^2/2, U contains g^2 phi^2 |X|^2 = g^2 phi^2 F^2/2 and lX |X|^4 = lX F^4/4, and the
rotation contributes -mu^2|X|^2 = -mu^2 F^2/2 to rho - mu n.  At mu = (lp lX)^{1/4} both bulk phases
(phi=1,F=0) and (phi=0, F^2 = mu^2/lX) have omega = 0, so sigma = integral of omega.
Method A: scipy solve_bvp on phi''=dV/dphi, F''=dV/dF with Dirichlet ends (shooting-type BVP); first integral check.
Method B: gradient-flow relaxation of a 1D profile with Dirichlet ends, then sigma = int omega dz, virial check.
Result reported in units of sqrt(lp) v^3 (=0.1 v^3) and GeV^3 at v = 1 TeV.
"""
import json, os
import numpy as np
from scipy.optimize import minimize

g, lp, lX = 0.1, 0.01, 0.001
mu = (lp * lX) ** 0.25
F2s = mu**2 / lX
Fs = np.sqrt(F2s)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_tension.json")


def V(p, F):
    return lp / 4 * (p**2 - 1) ** 2 + g**2 * p**2 * F**2 / 2 + lX * F**4 / 4 - mu**2 * F**2 / 2


def dV(p, F):
    dp = lp * (p**2 - 1) * p + g**2 * p * F**2
    dF = g**2 * p**2 * F + lX * F**3 - mu**2 * F
    return dp, dF


# check V >= 0 on a grid (coexistence requires both minima global)
P, FF = np.meshgrid(np.linspace(-1.5, 1.5, 601), np.linspace(0, 3, 601), indexing="ij")
Vmin_grid = float(V(P, FF).min())

# ---- Method A: boundary-value solution of the static field equations (shooting-type BVP) ----
from scipy.integrate import solve_bvp


def bvp(L=150.0, w=10.0):
    z = np.linspace(-L, L, 3001)
    p = 0.5 * (1 + np.tanh(z / w)); F = Fs * 0.5 * (1 - np.tanh(z / w))
    dp = np.gradient(p, z); dF = np.gradient(F, z)

    def rhs(zz, Y):
        a, b = dV(Y[0], Y[2])
        return np.vstack([Y[1], a, Y[3], b])

    def bc(Ya, Yb):
        return np.array([Ya[0] - 0.0, Ya[2] - Fs, Yb[0] - 1.0, Yb[2] - 0.0])

    sol = solve_bvp(rhs, bc, z, np.vstack([p, dp, F, dF]), tol=1e-9, max_nodes=200000)
    zz = np.linspace(-L, L, 200001)
    Y = sol.sol(zz)
    grad = 0.5 * (Y[1]**2 + Y[3]**2)
    pot = V(Y[0], Y[2])
    G = float(np.trapezoid(grad, zz)); P_ = float(np.trapezoid(pot, zz))
    return G + P_, G, P_, bool(sol.success), float(np.max(np.abs(grad - pot)))


sigA, GA, PA, okA, firstint = bvp()
sigA2, GA2, PA2, okA2, firstint2 = bvp(L=250.0, w=15.0)

# ---- Method B: direct minimisation of the discretised 1D grand-potential functional ----
def relax(L=150.0, h=0.25):
    z = np.arange(-L, L + h / 2, h)
    n = z.size
    w = 10.0
    p_init = 0.5 * (1 + np.tanh(z / w))           # p -> 1 at +L (vacuum), 0 at -L (Q-matter)
    F_init = Fs * 0.5 * (1 - np.tanh(z / w))

    def unpack(u):
        p = np.concatenate([[0.0], u[:n - 2], [1.0]])
        F = np.concatenate([[Fs], u[n - 2:], [0.0]])
        return p, F

    def E(u):
        p, F = unpack(u)
        dp = np.diff(p); dF = np.diff(F)
        e = 0.5 * np.sum(dp**2 + dF**2) / h + h * np.sum(V(p, F)[1:-1]) + 0.5 * h * (V(p[0], F[0]) + V(p[-1], F[-1]))
        gp = (2 * p[1:-1] - p[:-2] - p[2:]) / h
        gF = (2 * F[1:-1] - F[:-2] - F[2:]) / h
        vp, vF = dV(p[1:-1], F[1:-1])
        return e, np.concatenate([gp + h * vp, gF + h * vF])

    u0 = np.concatenate([p_init[1:-1], F_init[1:-1]])
    res = minimize(E, u0, jac=True, method="L-BFGS-B", options=dict(maxiter=200000, maxfun=10**7, ftol=1e-16, gtol=1e-13))
    p, F = unpack(res.x)
    G = float(0.5 * np.sum(np.diff(p)**2 + np.diff(F)**2) / h)
    Pot = float(res.fun - G)
    centre = float(np.interp(0.5, p, z))
    return float(res.fun), G, Pot, centre


sigB, GB, PB, cB = relax()
sigB2, GB2, PB2, cB2 = relax(L=250.0, h=0.125)

unit = np.sqrt(lp)  # sqrt(lp) v^3 with v = 1
out = dict(
    mu_v=mu, Fstar=Fs, F2star=F2s, Vmin_on_grid=Vmin_grid,
    sigma_bvp_v3=sigA, sigma_bvp_units=sigA / unit, bvp_gradient=GA / unit, bvp_potential=PA / unit, bvp_ok=okA,
    bvp_first_integral_maxdev=firstint, sigma_bvp2_units=sigA2 / unit, bvp2_ok=okA2,
    sigma_relax_v3=sigB, sigma_relax_units=sigB / unit, relax_gradient=GB / unit, relax_potential=PB / unit,
    relax_centre=cB,
    sigma_relax2_units=sigB2 / unit, relax2_gradient=GB2 / unit, relax2_potential=PB2 / unit, relax2_centre=cB2,
    sigma_GeV3_at_v1TeV=sigB2 * 1e9, sigma_cuberoot_GeV=(sigB2 * 1e9) ** (1 / 3),
    tree_wall_units=(2 * np.sqrt(2 * lp) / 3) / unit,
)
with open(OUT, "w") as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
