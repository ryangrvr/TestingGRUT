"""Fresh adjacent-charge-block checks; imports no earlier campaign code.

Finite small receivers are diagnostic models, not replacements for L=64N^2.
All photon occupations in each declared conserved-charge block are included.
"""
from itertools import permutations
from functools import lru_cache
from math import sqrt, factorial
from fractions import Fraction as Q
import numpy as np
from scipy.special import j1


def compositions(total, modes):
    if modes == 1:
        if total >= 0:
            yield (total,)
        return
    for n in range(total + 1):
        for rest in compositions(total - n, modes - 1):
            yield (n,) + rest


class Blocks:
    def __init__(self, M, L, delta, g):
        self.M, self.L, self.delta, self.lam = M, L, np.array(delta), g / (M + 1)
        self.words = [w for k in range(M + 1) for w in permutations(range(M), k)]

    @lru_cache(None)
    def basis(self, c):
        states = [(w, n) for w in self.words for n in compositions(len(w) + c, self.L)]
        return states, {x: k for k, x in enumerate(states)}

    @lru_cache(None)
    def H(self, c):
        states, idx = self.basis(c)
        H = np.zeros((len(states), len(states)), complex)
        for x, col in idx.items():
            w, n = x
            H[col, col] = sum(self.delta[v] for v in w)
            for r in range(self.L - 1):
                if n[r]:
                    q = list(n); q[r] -= 1; q[r + 1] += 1
                    row = idx[(w, tuple(q))]
                    z = -sqrt(n[r] * (n[r + 1] + 1))
                    H[row, col] += z; H[col, row] += z
            for v in range(self.M):
                if v not in w:
                    q = list(n); q[0] += 1
                    row = idx[(w + (v,), tuple(q))]
                    z = -self.lam * sqrt(n[0] + 1)
                    H[row, col] += z; H[col, row] += z
        assert np.linalg.norm(H - H.conj().T) == 0
        return H

    @lru_cache(None)
    def eig(self, c):
        return np.linalg.eigh(self.H(c))

    def evolve(self, c, x, t):
        e, v = self.eig(c)
        return v @ (np.exp(-1j * t * e) * (v.conj().T @ x))

    def ann(self, c, r):
        source, _ = self.basis(c); _, target = self.basis(c - 1)
        A = np.zeros((len(target), len(source)), complex)
        for col, (w, n) in enumerate(source):
            if n[r]:
                q = list(n); q[r] -= 1
                A[target[(w, tuple(q))], col] = sqrt(n[r])
        return A

    def append(self, c, phases=None):
        source, _ = self.basis(c); _, target = self.basis(c - 1)
        S = np.zeros((len(target), len(source)), complex)
        if phases is None:
            phases = np.ones(self.M)
        for col, (w, n) in enumerate(source):
            for v in range(self.M):
                if v not in w:
                    S[target[(w + (v,), n)], col] = phases[v]
        return S


def run():
    rows = []
    for M, L, delta, g in [(2, 4, [-.75, .75], .9), (3, 3, [.71, -.62, .13], 1.1)]:
        b = Blocks(M, L, delta, g)
        H, Hminus = b.H(0), b.H(-1)
        S, A = b.append(0), b.ann(0, 0)
        a1 = b.ann(0, 1)
        # a H_0 - H_-1 a = (h a)_0 - lambda S.
        identity = np.linalg.norm(A @ H - Hminus @ A + a1 + b.lam * S)
        assert identity < 1e-13
        seed = np.zeros(len(H), complex); seed[0] = 1
        hc = np.diag(-np.ones(L - 1), 1) + np.diag(-np.ones(L - 1), -1)
        eh, vh = np.linalg.eigh(hc)
        points, weights = np.polynomial.legendre.leggauss(80)
        # Schrodinger-space convolution, rather than the earlier Heisenberg code.
        for t in [.23, 1.7, 5.3]:
            psi = b.evolve(0, seed, t)
            integral = np.zeros(len(Hminus), complex)
            for p, weight in zip(points, weights):
                u = (p + 1) * t / 2
                C = np.dot(abs(vh[0]) ** 2, np.exp(-1j * eh * u))
                integral += weight * t / 2 * C * b.evolve(-1, S @ b.evolve(0, seed, t - u), u)
            rhs = 1j * b.lam * integral
            err = float(np.linalg.norm(A @ psi - rhs))
            assert err < 5e-12
            k = np.array([len(w) for w, n in b.basis(0)[0]])
            current = float((1j * np.vdot(psi, H @ (k * psi) - k * (H @ psi))).real)
            field_current = 2 * b.lam * np.vdot(S @ psi, A @ psi).imag
            assert abs(current - field_current) < 3e-13
            creation = b.ann(1, 0).conj().T @ psi
            ccr = abs(np.linalg.norm(creation) ** 2 - np.linalg.norm(A @ psi) ** 2 - 1)
            assert ccr < 2e-13
            rows.append(dict(N=M + 1, L=L, charge_dimensions={str(c): len(b.basis(c)[0]) for c in [-1, 0, 1]},
                t=t, intertwiner_error=float(identity), retarded_vector_error=err,
                current_error=float(abs(current - field_current)), vacuum_plus_one_error=float(ccr)))
        phases = np.exp(-1j * b.delta * .73)
        comm = b.append(0) @ b.append(1, phases) - b.append(0, phases) @ b.append(1)
        omitted = comm @ (b.ann(1, 0).conj().T @ seed)
        assert np.linalg.norm(omitted) > .01
        rows.append(dict(N=M + 1, omitted_creation_commutator_norm=float(np.linalg.norm(omitted)),
                         commuting_append_mutation_rejected=True))

    # Exact integers count closed walks, using a recurrence rather than matrix powers.
    walks = []
    for L in [1, 2, 4, 7, 12]:
        finite = {0: 1}; half = {0: 1}
        for n in range(2 * L + 1):
            if n < 2 * L:
                assert finite.get(0, 0) == half.get(0, 0)
            else:
                assert half.get(0, 0) - finite.get(0, 0) == 1
            def advance(state, end):
                out = {}
                for x, count in state.items():
                    for y in [x - 1, x + 1]:
                        if 0 <= y < end:
                            out[y] = out.get(y, 0) + count
                return out
            finite, half = advance(finite, L), advance(half, 2 * L + 2)
        walks.append(dict(L=L, first_boundary_difference=2 * L, exact_integer=True))

    grid = np.r_[np.geomspace(1e-5, 1, 80), np.geomspace(1, 2e4, 300)]
    C = abs(j1(2 * grid) / grid)
    bound = np.minimum(1, 2 * grid ** -1.5)
    assert np.all(C <= bound + 1e-14)
    rational = Q(1, 128) + Q(35, 72) + Q(146, 2**28) + Q(1, 2**100)
    assert rational < Q(1, 2)
    assert sum(Q(729, 125)**k / factorial(k) for k in range(41)) > 256
    assert Q(16, 7) * Q(3, 8)**128 < Q(1, 2**120)
    assert Q(4, 9) + Q(2, 3) + Q(4, 5) < 2
    return dict(status='PASS', charge_block_checks=rows, closed_walk_checks=walks,
                kernel_grid_max_ratio=float(np.max(C / bound)), formation_lower_exact=str(1 - rational),
                scope='fresh implementation diagnostics and exact rational constants; not external review or all-N numerical certification')
