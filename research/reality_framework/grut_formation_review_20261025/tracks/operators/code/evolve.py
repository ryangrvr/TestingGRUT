"""Own fast radial evolution (numpy, shift stencils) of the packet's three-field action with
the causal CFS memory, plus my re-implementation of the frozen-snapshot bound-mode split.

State arrays: U, V of shape (3, n+1) on nodes 0..n (U[:,0]=U[:,n]=0); P (3, n) on faces;
Qm (3, n+1) on nodes. Fields: 0 = r x, 1 = r y, 2 = r (S-1).
"""
import numpy as np
from scipy.linalg import eig_banded
from lattice import A0, A1, A2, A3, C1, C3


class Evolver:
    def __init__(self, n, h, R, smax=20., alpha=.1, nonlinear=True, mass=1.0):
        self.n, self.h, self.alpha = n, h, alpha
        self.r = np.arange(n + 1) * h
        self.rs = self.r.copy(); self.rs[0] = 1.0            # avoid 0/0 (node 0 is never updated)
        W = n * h - R
        self.sn = smax * np.clip((self.r - R) / W, 0, 1) ** 3
        self.sf = smax * np.clip(((np.arange(n) + .5) * h - R) / W, 0, 1) ** 3
        self.sn[0] = self.sn[-1] = 0.0
        self.nonlinear = nonlinear
        self.mass = mass
        self.w = 4 * np.pi * h
        self.inner = slice(1, n)

    def ext_nodes(self, U):
        n = self.n
        return np.concatenate([-U[:, 3:0:-1], U, -U[:, n - 1:n - 4:-1]], axis=1)  # index k <-> node k-3

    def ext_faces(self, P):
        n = self.n
        return np.concatenate([P[:, 1::-1], P, P[:, n - 1:n - 3:-1]], axis=1)  # index k <-> face k-2

    def grad(self, U):
        e = self.ext_nodes(U); n = self.n
        # face f uses nodes f-1, f, f+1, f+2  -> ext indices f+2 .. f+5
        return (C1 * (e[:, 4:n + 4] - e[:, 3:n + 3]) + C3 * (e[:, 5:n + 5] - e[:, 2:n + 2])) / self.h

    def div(self, P):
        e = self.ext_faces(P); n = self.n
        out = np.zeros((P.shape[0], n + 1))
        i = np.arange(1, n)
        # node i uses faces i-2, i-1, i, i+1 -> ext indices i, i+1, i+2, i+3
        out[:, 1:n] = (C1 * (e[:, i + 2] - e[:, i + 1]) + C3 * (e[:, i + 3] - e[:, i])) / self.h
        return out

    def lap(self, U):
        e = self.ext_nodes(U); n = self.n
        out = np.zeros_like(U)
        c = lambda d: e[:, 3 + d + 1:3 + d + n]   # nodes 1+d .. n-1+d
        out[:, 1:n] = (A0 * c(0) + A1 * (c(1) + c(-1)) + A2 * (c(2) + c(-2)) + A3 * (c(3) + c(-3))) / self.h ** 2
        return out

    def force(self, U):
        r = self.rs
        if not self.nonlinear:
            F = np.empty_like(U)
            F[:2] = self.mass ** 2 * U[:2]; F[2] = 6 * U[2]
            return F
        x, y, z = U
        rho = (x * x + y * y) / (r * r)
        S = 1 + z / r
        g = S * S + .1 * rho
        F = np.array([g * x, g * y, r * (3 * (S * S - 1) + rho) * S])
        F[:, 0] = 0; F[:, -1] = 0
        return F

    def rhs(self, Y):
        U, V, P, Qm = Y
        LU = self.lap(U)
        DP = self.div(P)
        B = DP + Qm
        dV = LU - self.force(U) + B
        dV[:, 0] = 0; dV[:, -1] = 0
        dP = -(self.alpha + self.sf) * P - self.sf * self.grad(U)
        dQ = -(self.alpha + self.sn) * Qm - self.sn * (LU + DP)
        dQ[:, 0] = 0; dQ[:, -1] = 0
        return [V, dV, dP, dQ]

    def step(self, Y, dt):
        k1 = self.rhs(Y)
        k2 = self.rhs([y + .5 * dt * k for y, k in zip(Y, k1)])
        k3 = self.rhs([y + .5 * dt * k for y, k in zip(Y, k2)])
        k4 = self.rhs([y + dt * k for y, k in zip(Y, k3)])
        return [y + dt / 6 * (a + 2 * b + 2 * c + d) for y, a, b, c, d in zip(Y, k1, k2, k3, k4)]

    def charge_density(self, U, V):
        return U[1] * V[0] - U[0] * V[1]

    def Q(self, U, V):
        return self.w * np.sum(self.charge_density(U, V)[1:self.n])

    def Q_within(self, U, V, R):
        q = self.charge_density(U, V)
        r = self.r
        m = (r > 0) & (r < R - self.h * .1)
        edge = abs(r - R) < self.h * .1
        return self.w * (np.sum(q[m]) + .5 * np.sum(q[edge]))

    def energy(self, U, V):
        LU = self.lap(U)
        r = self.rs
        x, y, z = U
        rho = (x * x + y * y) / (r * r)
        S = 1 + z / r
        pot = r * r * (.75 * (S * S - 1) ** 2 + .5 * S * S * rho + .025 * rho * rho)
        e = np.sum(.5 * V * V - .5 * U * LU, axis=0) + pot
        return self.w * np.sum(e[1:self.n])


def frozen_bound_split(ev, U, V, upper=1 - 1e-9, rmax_core=60.):
    """My re-implementation of the packet's frozen instantaneous charged-operator split:
    H = -L + S^2 + 0.1 |X|^2 on nodes 1..n-1 (odd ghosts, wall), eigenvalues in (0, upper)
    counted as bound; canonical charge per real orthonormal eigenvector."""
    n, h = ev.n, ev.h
    r = ev.r[1:n]
    x, y, z = U[0, 1:n], U[1, 1:n], U[2, 1:n]
    S = 1 + z / r
    pot = S * S + .1 * (x * x + y * y) / (r * r)
    N = n - 1
    band = np.zeros((4, N))
    band[0] = pot - A0 / h ** 2
    band[0, 0] += A2 / h ** 2; band[0, -1] += A2 / h ** 2
    band[1, :-1] = -A1 / h ** 2; band[1, 0] += A3 / h ** 2; band[1, -2] += A3 / h ** 2
    band[2, :-2] = -A2 / h ** 2
    band[3, :-3] = -A3 / h ** 2
    vals, vecs = eig_banded(band, lower=True, select='v', select_range=(0., upper))
    u = x + 1j * y
    p = V[0, 1:n] + 1j * V[1, 1:n]
    Uj, Pj = vecs.T @ u, vecs.T @ p
    q = -ev.w * np.imag(Uj.conj() * Pj)
    inside = np.array([np.sum(vecs[r < rmax_core, j] ** 2) for j in range(len(vals))])
    return dict(omegas=np.sqrt(vals).tolist(), mode_charges=q.tolist(), Qbound=float(np.sum(q)),
                inside60=inside.tolist())
