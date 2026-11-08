"""Lattice Schwinger model (1+1D QED) in the Q=0 sector, Gauss law solved exactly.

Dimensionless W = 2H/(a g^2); physical H/g = W / (2 sqrt(x)).
"""
import itertools
import numpy as np
import scipy.sparse as sp


def basis(N):
    states = []
    for ups in itertools.combinations(range(N), N // 2):
        s = 0
        for u in ups:
            s |= 1 << u
        states.append(s)
    return np.array(sorted(states), dtype=np.int64)


def bits(states, N):
    return ((states[:, None] >> np.arange(N)) & 1).astype(np.int8)  # (D,N), 1 = up


def electric_fields(B):
    N = B.shape[1]
    sz = 2 * B - 1
    stag = (-1) ** np.arange(N)
    q = (sz + stag) / 2.0                 # staggered charge per site
    L = np.cumsum(q, axis=1)[:, :-1]      # links 0..N-2
    return L


def hamiltonian(N, x, m_over_g):
    st = basis(N)
    B = bits(st, N)
    D = len(st)
    sz = 2 * B - 1
    mu = 2 * np.sqrt(x) * m_over_g
    L = electric_fields(B)
    diag = (L ** 2).sum(1) + 0.5 * mu * (1 + ((-1) ** np.arange(N)) * sz).sum(1)
    rows, cols, vals = [np.arange(D)], [np.arange(D)], [diag]
    for n in range(N - 1):
        differ = B[:, n] != B[:, n + 1]
        idx = np.nonzero(differ)[0]
        new = st[idx] ^ ((1 << n) | (1 << (n + 1)))
        j = np.searchsorted(st, new)
        rows.append(idx); cols.append(j); vals.append(np.full(len(idx), x))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
    return H, L, st
