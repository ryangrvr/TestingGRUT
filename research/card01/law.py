"""Frozen Card 1: a single operator supplies phase and decoherence geometry.

Only these equations are the candidate. Other generators belong to hostile
comparators. Matrix units, preparation, measurement and the clock are supplied.
"""
from __future__ import annotations

import numpy as np


def commutator(a, x):
    return a @ x - x @ a


def generator(a, x, tau):
    """K=0 means dA/dt=0 and d rho/dt=generator(A,rho,tau)."""
    return -1j * commutator(a, x) - tau / 2 * commutator(a, commutator(a, x))


def channel(a, x, tau, t):
    """Exact finite-dimensional solution, also acting on non-density matrices."""
    if tau <= 0 or t < 0:
        raise ValueError("tau must be positive and t nonnegative")
    a = np.asarray(a, dtype=complex)
    x = np.asarray(x, dtype=complex)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or x.shape != a.shape:
        raise ValueError("square matrices of equal size required")
    if not np.allclose(a, a.conj().T, atol=1e-13, rtol=0):
        raise ValueError("A must be Hermitian")
    values, vectors = np.linalg.eigh(a)
    gaps = values[:, None] - values[None, :]
    factors = np.exp((-1j * gaps - tau / 2 * gaps**2) * t)
    return vectors @ (factors * (vectors.conj().T @ x @ vectors)) @ vectors.conj().T


def spectral_blocks(a):
    """Exact interpretation groups equal eigenvalues; no fitted slow-mode cutoff.

    This implementation uses the eigensolver's equality, so floating accidental
    near-degeneracies are not certified as exact identities. Tests use exact
    diagonal degeneracies. No tolerance is a physical selection rule.
    """
    values, vectors = np.linalg.eigh(a)
    out = []
    for value in np.unique(values):
        v = vectors[:, values == value]
        out.append((float(value), v @ v.conj().T))
    return out
