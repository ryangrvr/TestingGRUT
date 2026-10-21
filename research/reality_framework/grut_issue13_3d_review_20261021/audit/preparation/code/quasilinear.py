"""Quasi-linear preparation diagnostics for Codex a69bd86 (issue #13 audit, track 'preparation').

Model (conformal variables of expanding.cpp, radiation a = 1 + tau/40, D = 0):
  * linear S fluctuations:  chi_k'' = -(khat^2 + F^2 - a^2) chi_k,  chi_k(0)=A, chi_k'(0)=A/40
  * linear X fluctuations about the exact rotating orbit U=F cos(w t), W=-F sin(w t):
        du'' = -khat^2 du - 0.1[F^2 du + 2U(U du + W dw)]   (and u<->w)
    with the preparation's rigid-rotation momenta du'(0)=w dw(0), dw'(0)=-w du(0);
    dF2 = 2(U du + W dw)
  * the S zero mode (absent in every preparation, <S>=0 exactly on the grid) is
    sourced at second order by the box-level S-X cross-correlation:
        chi0'' = (a^2 - F^2) chi0 - <dF2 chi_S>_box  [- <chi^3>_box optionally]
All growth uses the lattice eigenvalue khat^2 = sum 4 sin^2(pi n_i/N)/dx^2.
This is an estimate [E] built from the preparation inputs [P]; it is compared with
Codex's recorded nonlinear mean_S series, not substituted for them.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

TAUI = 40.0
RATE = 1.0 / TAUI
EPS_S = 0.01
EPS_X = 0.001


def a_of(t):
    return 1.0 + RATE * t


def khat2_int(n, N, dx):
    n = np.atleast_2d(n).astype(float)
    return np.sum(4 * np.sin(np.pi * n / N) ** 2, axis=1) / dx ** 2


def s_growth(k2, kappa, tgrid):
    """g(k2,t): chi_S(t)/chi_S(0) for chi'(0)=chi(0)/40. Vectorized over k2."""
    k2 = np.asarray(k2, float)
    F2 = kappa ** (2 / 3)
    m = len(k2)

    def rhs(t, y):
        a = a_of(t)
        x, v = y[:m], y[m:]
        return np.concatenate([v, -(k2 + F2 - a * a) * x])

    y0 = np.concatenate([np.ones(m), RATE * np.ones(m)])
    sol = solve_ivp(rhs, [tgrid[0], tgrid[-1]], y0, t_eval=tgrid, method='DOP853', rtol=1e-10, atol=1e-13)
    return sol.y[:m].T, sol.y[m:].T  # shape (nt, m)


def x_response(k2, kappa, tgrid):
    """Return Au(t,k2), Aw(t,k2): dF2 per unit initial du, dw (rigid-rotation momenta)."""
    k2 = np.asarray(k2, float)
    F = kappa ** (1 / 3)
    F2 = F * F
    w = np.sqrt(0.1) * F
    m = len(k2)
    out = []
    for init in ('u', 'w'):
        du0 = np.ones(m) if init == 'u' else np.zeros(m)
        dw0 = np.zeros(m) if init == 'u' else np.ones(m)
        y0 = np.concatenate([du0, dw0, w * dw0, -w * du0])

        def rhs(t, y):
            du, dw, pu, pw = y[:m], y[m:2 * m], y[2 * m:3 * m], y[3 * m:]
            U, W = F * np.cos(w * t), -F * np.sin(w * t)
            proj = U * du + W * dw
            fu = -k2 * du - 0.1 * (F2 * du + 2 * U * proj)
            fw = -k2 * dw - 0.1 * (F2 * dw + 2 * W * proj)
            return np.concatenate([pu, pw, fu, fw])

        sol = solve_ivp(rhs, [tgrid[0], tgrid[-1]], y0, t_eval=tgrid, method='DOP853', rtol=1e-10, atol=1e-13)
        du, dw = sol.y[:m].T, sol.y[m:2 * m].T
        U = (F * np.cos(w * tgrid))[:, None]
        W = (-F * np.sin(w * tgrid))[:, None]
        out.append(2 * (U * du + W * dw))
    return out[0], out[1]


def zero_mode(tgrid, J, kappa, cubic=None):
    """chi0'' = (a^2-F^2) chi0 + J(t) [+ cubic(t, chi0)], chi0(0)=chi0'(0)=0."""
    F2 = kappa ** (2 / 3)

    def Jt(t):
        return np.interp(t, tgrid, J)

    def rhs(t, y):
        a = a_of(t)
        f = (a * a - F2) * y[0] + Jt(t)
        if cubic is not None:
            f += cubic(t, y[0])
        return [y[1], f]

    sol = solve_ivp(rhs, [tgrid[0], tgrid[-1]], [0.0, 0.0], t_eval=tgrid, method='DOP853', rtol=1e-10, atol=1e-16,
                    max_step=0.05)
    return sol.y[0]


def lam_exponent(k2, kappa, a):
    """Leading WKB growth exponent 40*int_{a_k}^{a} sqrt(a'^2-a_k^2) da' (0 before crossing)."""
    c2 = kappa ** (2 / 3) + np.asarray(k2, float)
    c = np.sqrt(c2)
    a = np.asarray(a, float)
    s = np.sqrt(np.maximum(a * a - c2, 0))
    val = 0.5 * (a * s - c2 * np.log(np.maximum((a + s) / c, 1.0)))
    return TAUI * np.where(a > c, val, 0.0)
