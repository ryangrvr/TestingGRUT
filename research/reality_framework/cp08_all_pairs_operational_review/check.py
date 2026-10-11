"""Independent exact CP-08 branch-resolved causal response; no new Hamiltonian.

This independently implements frozen CP-08 LAW.md instead of importing
its numerical module. Requires numpy and scipy; run python check.py.
"""
import itertools
import math
import numpy as np
from scipy.linalg import expm

J = 1.0
OMEGA = 0.25
GAMMA = math.pi / 4
PAULI = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.diag([1., -1.]).astype(complex)
}


def site_operator(n, changes):
    op = np.array([[1.]], complex)
    for k in range(n):
        op = np.kron(op, PAULI[changes.get(k, 'I')])
    return op


def build_independent(n, gamma=GAMMA):
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    posets = []
    for bits in range(1 << (n * (n - 1))):
        C = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            C[i, j] = (bits >> k) & 1
        if np.any(C * C.T):
            continue
        if not np.all((C @ C > 0) <= (C > 0)):
            continue
        posets.append(bits)
    lookup = {v: i for i, v in enumerate(posets)}
    matter_dim = 1 << n
    geom_dim = len(posets)
    H = np.zeros((geom_dim * matter_dim, geom_dim * matter_dim), complex)
    Xsum = sum(site_operator(n, {i: 'X'}) for i in range(n))
    transport = {}
    for i, j in pairs:
        P = site_operator(n, {i: 'Z', j: 'X'})
        transport[i, j] = (np.cos(gamma) * np.eye(matter_dim) -
                           1j * np.sin(gamma) * P)
    for k, mask in enumerate(posets):
        sl = slice(k * matter_dim, (k + 1) * matter_dim)
        H[sl, sl] = OMEGA * Xsum
        for p, (i, j) in enumerate(pairs):
            new = mask | (1 << p)
            if new == mask or new not in lookup:
                continue
            l = lookup[new]
            sl2 = slice(l * matter_dim, (l + 1) * matter_dim)
            H[sl2, sl] += -J * transport[i, j]
            H[sl, sl2] += -J * transport[i, j].conj().T
    assert np.linalg.norm(H - H.conj().T) < 1e-10
    return H, posets, pairs, lookup


def comm(A, B):
    return A @ B - B @ A


def coefficient_by_geometry(H, n, lookup, i, j):
    m = 1 << n
    src = site_operator(n, {i: 'X'})
    dst = site_operator(n, {j: 'Z'})
    Ypair = site_operator(n, {i: 'Y', j: 'Y'})
    B = np.kron(np.eye(len(lookup)), src)
    total = np.kron(np.eye(len(lookup)), dst)
    pieces = {}
    for mask, k in lookup.items():
        M = np.zeros_like(total)
        sl = slice(k * m, (k + 1) * m)
        M[sl, sl] = dst
        pieces[mask] = (-.5j * comm(comm(H, comm(H, M)), B))[:m, :m]
    qtotal = (-.5j * comm(comm(H, comm(H, total)), B))[:m, :m]
    assert np.linalg.norm(sum(pieces.values()) - qtotal) < 1e-11
    expected = -2 * J**2 * math.sin(2 * GAMMA) * Ypair
    assert np.linalg.norm(qtotal - expected) < 1e-11
    pairs = list(itertools.permutations(range(n), 2))
    forward = 1 << pairs.index((i, j))
    assert np.linalg.norm(pieces[forward] - expected) < 1e-11
    assert all(np.linalg.norm(Q) < 1e-11
               for mask, Q in pieces.items() if mask != forward)
    return forward


def branch_responses(H, n, lookup, i, j, t):
    m = 1 << n
    plus = np.array([1, 1j], complex) / math.sqrt(2)
    minus = np.array([1, -1j], complex) / math.sqrt(2)
    zero = np.array([1, 0], complex)
    state = np.array([1 + 0j])
    for k in range(n):
        state = np.kron(state, plus if k == i else minus if k == j else zero)
    psi = np.zeros(len(lookup) * m, complex)
    psi[:m] = state
    B = np.kron(np.eye(len(lookup)), site_operator(n, {i: 'X'}))
    U = expm(-1j * H * t)
    w, wB = U @ psi, U @ (B @ psi)
    vals = {}
    for mask, k in lookup.items():
        sl = slice(k * m, (k + 1) * m)
        D = site_operator(n, {j: 'Z'})
        c = np.vdot(w[sl], D @ wB[sl])
        vals[mask] = float(-2 * c.imag)
    return vals


def empty_survival(H, n, t):
    m = 1 << n
    state = np.eye(1, m, 0, dtype=complex).reshape(-1)
    psi = np.zeros(H.shape[0], complex)
    psi[:m] = state
    w = expm(-1j * H * t) @ psi
    return float(np.vdot(w[:m], w[:m]).real)


def main():
    for n, required in [(2, 3), (3, 19)]:
        H, posets, pairs, lookup = build_independent(n)
        assert len(posets) == required
        for i, j in itertools.permutations(range(n), 2):
            forward = coefficient_by_geometry(H, n, lookup, i, j)
            if (i, j) == (0, 1):
                for t in (0.04, 0.02, 0.01):
                    vals = branch_responses(H, n, lookup, i, j, t)
                    print(f'N={n}, t={t}, '
                          f'forward_response_over_t2={vals[forward]/t**2:.10f}, '
                          f'empty_response_over_t2={vals[0]/t**2:.10f}, '
                          f'other_branches_over_t2='
                          f'{sum(v for k,v in vals.items() if k!=forward)/t**2:.10f}')
                    assert abs(vals[forward] / t**2 - 2) < 0.02
                    assert abs(vals[0] / t**2) < 0.03
                    assert abs(sum(vals.values()) / t**2 - 2) < 0.02
        for t in (0.04, 0.02, 0.01):
            rate = (1 - empty_survival(H, n, t)) / t**2
            print(f'N={n}, t={t}: empty_departure_over_t2={rate:.10f}, '
                  f'exact_coefficient={J**2*n*(n-1):.0f}')
            assert abs(rate - J**2 * n * (n - 1)) < 0.1
        print(f'N={n}: all {n*(n-1)} ordered source/readout pairs '
              f'have the same operator coefficient; only their corresponding '
              f'one-edge final relation carries the t^2 term')
    ypair = site_operator(2, {0: 'Y', 1: 'Y'})
    for gamma in (0.0, math.pi/8, math.pi/4, math.pi/2):
        H, states, pairs, lookup = build_independent(2, gamma)
        A = np.kron(np.eye(len(states)), site_operator(2, {1: 'Z'}))
        B = np.kron(np.eye(len(states)), site_operator(2, {0: 'X'}))
        Q = (-.5j * comm(comm(H, comm(H, A)), B))[:4, :4]
        expected = -2 * J**2 * math.sin(2 * gamma) * ypair
        assert np.linalg.norm(Q - expected) < 1e-11
        print(f'gamma={gamma:.6f}: '
              f'exact_cross_coefficient={-2*math.sin(2*gamma):.10f}')
    print('PASS: uniform all-pairs short-time response, branch-resolved '
          'support, and empty-order quadratic escape, N=2,3')


if __name__ == '__main__':
    main()
