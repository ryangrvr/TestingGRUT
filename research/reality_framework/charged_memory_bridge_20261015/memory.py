"""Sourced amplitude memory equation and a coupling-rank Volterra reduction.

All spectral preprocessing is finite dimensional. The equation, not the
preprocessing cost, is reduced to the numerical rank of the inter-sector
coupling. There is no rate fit, bath, measurement or Markov assumption.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh, svd


class MemoryBridge:
    def __init__(self, A, C, D, z0, relative_rank_tolerance=1e-12):
        self.A = np.asarray(A)
        self.C = np.asarray(C)
        self.D = np.asarray(D)
        self.z0 = np.asarray(z0)
        self.alpha, self.UA = eigh(self.A)
        self.delta, self.UD = eigh(self.D)
        U, sigma, Vh = svd(self.C, full_matrices=False)
        keep = sigma > relative_rank_tolerance * max(1., sigma[0])
        self.rank = int(keep.sum())
        self.U = U[:, keep]
        self.sigma = sigma[keep]
        self.V = Vh[keep].conj().T
        self.coupling_residual = float(np.linalg.norm(
            self.C - (self.U * self.sigma) @ self.V.conj().T, ord=2))
        self.Ceig = self.UA.conj().T @ self.C @ self.UD
        self.zeig = self.UD.conj().T @ self.z0
        self.return_leg = self.UD.conj().T @ (self.V * self.sigma)
        self.input_overlap = self.U.conj().T @ self.UA

    def integrated_blocks(self, t):
        """Return first-transfer g(t) and J(t) with L(t)=J(t) U^dagger.

        The oscillatory divided difference is evaluated with sinc so exact
        or near eigenvalue degeneracies need no special subtraction.
        """
        a, d = self.alpha[:, None], self.delta[None, :]
        I = t * np.exp(-.5j * (a + d) * t) * np.sinc((d - a) * t / (2*np.pi))
        F = I * self.Ceig
        g = -1j * self.UA @ (F @ self.zeig)
        J = self.UA @ (F @ self.return_leg)
        return g, J

    def memory_kernel(self, t):
        """K_D(t) on the coupling interface; K(t)=U K_D(t) U^dagger."""
        V = self.return_leg
        return V.conj().T @ (np.exp(-1j*self.delta*t)[:, None] * V)

    def solve(self, end, steps, readout_times):
        """Composite trapezoid for a Volterra equation of the second kind.

        Both y(0) and the integrated kernel at zero vanish. Thus this
        otherwise implicit rule is an explicit causal convolution. The
        grid study in run.py measures its error; this routine does not
        assert an interval-certified discretization error.
        """
        grid = np.linspace(0., end, steps+1)
        dt = end/steps
        g = np.empty((steps+1, self.A.shape[0]), dtype=complex)
        J = np.empty((steps+1, self.A.shape[0], self.rank), dtype=complex)
        for i, t in enumerate(grid):
            g[i], J[i] = self.integrated_blocks(float(t))
        h = g @ self.U.conj()
        L = np.einsum('ia,taj->tij', self.U.conj().T, J, optimize=True)
        y = np.zeros((steps+1, self.rank), dtype=complex)
        for n in range(1, steps+1):
            y[n] = h[n] - dt*np.einsum(
                'kij,kj->i', L[1:n], y[n-1:0:-1], optimize=True)
        readouts = []
        for t in readout_times:
            n = int(round(t/dt))
            assert abs(grid[n]-t) < 1e-12
            x = g[n] - dt*np.einsum(
                'kaj,kj->a', J[1:n], y[n-1:0:-1], optimize=True)
            readouts.append(dict(time=float(t), x=x,
                probability=float(np.vdot(x, x).real),
                first_transfer_probability=float(np.vdot(g[n], g[n]).real),
                interface_consistency_error=float(np.linalg.norm(self.U.conj().T@x-y[n]))))
        return readouts


def orthogonal_columns(X):
    U, s, _ = svd(X, full_matrices=False)
    return U[:, s > 1e-11 * max(1., s[0])]


def reservoir_krylov(D, C, z0, maximum_depth):
    """Block Krylov on the complement, seeded by source and return ports."""
    V = orthogonal_columns(np.column_stack([np.real_if_close(z0), C.conj().T]))
    block = V.copy()
    yield 1, V.copy()
    for depth in range(2, maximum_depth+1):
        W = D @ block
        for _ in range(2):
            W -= V @ (V.conj().T @ W)
        block = orthogonal_columns(W)
        if not block.shape[1]:
            break
        V = np.column_stack([V, block])
        yield depth, V.copy()
