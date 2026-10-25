"""Stationary charged soliton of the packet's two-field radial action (own code).

X = f(r) exp(-i w0 t), S = s(r); in U-form U = r f, Z = r (s - 1):
   -L U + (s^2 + 0.1 f^2) U = w0^2 U
   -L Z + r (3 (s^2 - 1) + f^2) s = 0
with the same seven-point L and odd closures. Newton with a dense Jacobian.
"""
import numpy as np
from lattice import L_matrix


def residual(U, Z, r, L, w0):
    f = U / r
    s = 1 + Z / r
    F1 = -L @ U + (s * s + .1 * f * f) * U - w0 * w0 * U
    F2 = -L @ Z + r * (3 * (s * s - 1) + f * f) * s
    return np.r_[F1, F2]


def jacobian(U, Z, r, L, w0):
    f = U / r
    s = 1 + Z / r
    N = len(r)
    J = np.zeros((2 * N, 2 * N))
    J[:N, :N] = -L + np.diag(s * s + .3 * f * f - w0 * w0)
    J[:N, N:] = np.diag(2 * s * f)                 # d/dZ of s^2 U = 2 s (1/r) U = 2 s f
    J[N:, :N] = np.diag(2 * f * s)                 # d/dU of r f^2 s = 2 f s
    J[N:, N:] = -L + np.diag(9 * s * s - 3 + f * f)
    return J


def solve(w0, h=.25, Rb=40., amp=2.0, width=4.0, tol=1e-11, maxit=60, guess=None):
    n = int(round(Rb / h))
    r = np.arange(1, n) * h
    L = L_matrix(n, h)
    if guess is None:
        f = amp * np.exp(-r ** 2 / (2 * width ** 2))
        s = np.sqrt(np.clip(1 - f * f / 3, 0.0, 1))
        U, Z = r * f, r * (s - 1)
    else:
        U, Z = guess
    for it in range(maxit):
        F = residual(U, Z, r, L, w0)
        nr = np.max(abs(F))
        if nr < tol:
            break
        dx = np.linalg.solve(jacobian(U, Z, r, L, w0), -F)
        lam = 1.0
        while lam > 1e-3:
            U2, Z2 = U + lam * dx[:n - 1], Z + lam * dx[n - 1:]
            if np.max(abs(residual(U2, Z2, r, L, w0))) < (1 - .3 * lam) * nr:
                break
            lam /= 2
        U, Z = U2, Z2
    F = residual(U, Z, r, L, w0)
    return dict(r=r, U=U, Z=Z, h=h, n=n, w0=w0, residual=float(np.max(abs(F))), iterations=it,
                charge=float(4 * np.pi * h * np.sum(w0 * U * U)), fmax=float(np.max(abs(U / r))))


def solve_fixed_norm(Nnorm, h=.25, Rb=40., amp=3.0, width=6.4, tol=1e-10, maxit=200, w0_guess=.85):
    """Same equations with lambda = w0^2 unknown and h*sum(U^2) = Nnorm fixed (excludes U = 0)."""
    n = int(round(Rb / h))
    r = np.arange(1, n) * h
    L = L_matrix(n, h)
    N = n - 1
    f = amp * np.exp(-r ** 2 / (2 * width ** 2))
    s = np.sqrt(np.clip(1 - f * f / 3, 0.0, 1))
    U, Z = r * f, r * (s - 1)
    U *= np.sqrt(Nnorm / (h * np.sum(U * U)))
    lam = w0_guess ** 2

    def full(U, Z, lam):
        return np.r_[residual(U, Z, r, L, np.sqrt(lam)), h * np.sum(U * U) - Nnorm]

    for it in range(maxit):
        F = full(U, Z, lam)
        nr = np.max(abs(F[:-1])) + abs(F[-1]) / Nnorm
        if nr < tol:
            break
        J = np.zeros((2 * N + 1, 2 * N + 1))
        J[:2 * N, :2 * N] = jacobian(U, Z, r, L, np.sqrt(lam))
        J[:N, -1] = -U
        J[-1, :N] = 2 * h * U
        dx = np.linalg.solve(J, -F)
        step = 1.0
        while step > 1e-4:
            U2, Z2, l2 = U + step * dx[:N], Z + step * dx[N:2 * N], lam + step * dx[-1]
            F2 = full(U2, Z2, l2)
            if np.max(abs(F2[:-1])) + abs(F2[-1]) / Nnorm < (1 - .2 * step) * nr:
                break
            step /= 2
        U, Z, lam = U2, Z2, l2
    w0 = float(np.sqrt(lam))
    return dict(r=r, U=U, Z=Z, h=h, n=n, w0=w0, residual=float(np.max(abs(residual(U, Z, r, L, w0)))), iterations=it,
                charge=float(4 * np.pi * h * np.sum(w0 * U * U)), fmax=float(np.max(abs(U / r))),
                smin=float(np.min(1 + Z / r)))
