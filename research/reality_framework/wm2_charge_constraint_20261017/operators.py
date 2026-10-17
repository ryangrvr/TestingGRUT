"""Exact physical-kernel isometry for small WM2-UC1 diagnostic models.

All rooted label words, storage masks and charge-allowed two-mode receiver
occupations are kept. Electric fields are eliminated by the proved unique
tree Gauss solution, not by an energy or numerical truncation. These small
receivers do not replace the physical S9 receiver in a formation claim.
"""
from functools import lru_cache
from itertools import permutations
from math import sqrt
import numpy as np
from scipy.sparse import coo_matrix

GMAX = 1.1
EPSILON = .25
JA, JB, KAPPA = 1 / 64, 1 / 96, .25


@lru_cache(None)
def pairs(N):
    return tuple((i, j) for i in range(N) for j in range(i + 1, N))


def bonds(w):
    return {tuple(sorted((u, v))) for u, v in zip(w, w[1:])}


@lru_cache(None)
def electric(N, w, pA, pB):
    out = {edge: 0 for edge in pairs(N)}
    for p, charge in [(pA, 1), (pB, 2)]:
        if p is None:
            continue
        stop = w.index(p)
        for u, v in zip(w[:stop], w[1:stop + 1]):
            out[tuple(sorted((u, v)))] += -charge if u < v else charge
    return tuple(out[edge] for edge in pairs(N))


def constraint_values(N, w, pA, pB, fields):
    divergence = [0] * N
    for (u, v), value in zip(pairs(N), fields):
        divergence[u] += value; divergence[v] -= value
    for p, charge in [(pA, 1), (pB, 2)]:
        if p is not None:
            divergence[p] -= charge; divergence[w[0]] += charge
    b = bonds(w)
    support = [(0 if edge in b else e * e) for edge, e in zip(pairs(N), fields)]
    return tuple(divergence), tuple(support)


def field_difference(a, b):
    return tuple(k for k, (x, y) in enumerate(zip(a, b)) if x != y)


class KernelModel:
    def __init__(self, N, delta=None, g=1., check=True):
        self.N, self.g = N, g
        self.delta = tuple(.23 + np.linspace(-.25, .25, N) if delta is None else delta)
        self.words = [(0,) + w for k in range(N) for w in permutations(range(1, N), k)]
        self.states = []
        for w in self.words:
            a_states = [None] + [(p, level) for p in w for level in range(2)]
            b_states = [None] + [(p, level) for p in w for level in range(3)]
            for q in range(1, 1 << N, 2):
                for n0 in range(q.bit_count()):
                    for a in a_states:
                        for b in b_states:
                            self.states.append((w, q, n0, a, b))
        self.index = {state: i for i, state in enumerate(self.states)}
        self.audit = {'checked_geometry_edges': 0, 'checked_charge_hops': 0,
                      'maximum_electric_magnitude': 0, 'blocked_nonzero_flux_pops': 0}
        self.H, self.Hmove = self._build(check)

    def fields(self, state):
        w, q, n0, a, b = state
        return electric(self.N, w, None if a is None else a[0], None if b is None else b[0])

    def _build(self, check):
        rows, cols, values = [], [], []
        mr, mc, mv = [], [], []
        def add(target, col, value, moving=False):
            row = self.index[target]
            rows.append(row); cols.append(col); values.append(value)
            if moving:
                mr.append(row); mc.append(col); mv.append(value)

        for col, state in enumerate(self.states):
            w, q, n0, a, b = state
            mask = sum(1 << v for v in w)
            efield = self.fields(state)
            if check:
                gv, cv = constraint_values(self.N, w, None if a is None else a[0],
                                            None if b is None else b[0], efield)
                assert not any(gv) and not any(cv)
                self.audit['maximum_electric_magnitude'] = max(self.audit['maximum_electric_magnitude'], max(map(abs, efield), default=0))
                assert max(map(abs, efield), default=0) <= 3
            diagonal = 16 * (q.bit_count() - len(w)) + 64 * (q ^ mask).bit_count()
            diagonal += sum(self.delta[v] for v in range(self.N) if q >> v & 1) - self.delta[0]
            if b is not None and b[1] == 2:
                diagonal += 1  # the third B level is physically present
            add(state, col, diagonal)
            total = q.bit_count() - 1
            for step in [-1, 1]:
                newn = n0 + step
                if 0 <= newn <= total:
                    z = -sqrt(n0 * (total - n0 + 1)) if step == -1 else -sqrt((n0 + 1) * (total - n0))
                    add((w, q, newn, a, b), col, z)

            def geom(target, z):
                if check:
                    assert self.fields(target) == efield, 'append/pop must not silently rewrite flux'
                    self.audit['checked_geometry_edges'] += 1
                add(target, col, z)

            for v in range(self.N):
                if v in w:
                    continue
                edge = tuple(sorted((w[-1], v)))
                assert efield[pairs(self.N).index(edge)] == 0
                if not q >> v & 1:
                    geom((w + (v,), q | (1 << v), n0 + 1, a, b), -self.g / self.N * sqrt(n0 + 1))
                gate = sum(p is not None and p[0] in [w[-1], v] for p in [a, b])
                if gate:
                    geom((w + (v,), q, n0, a, b), -EPSILON / sqrt(self.N) * gate)
            if len(w) > 1:
                edge = tuple(sorted((w[-2], w[-1])))
                if efield[pairs(self.N).index(edge)] == 0:
                    v = w[-1]
                    if q >> v & 1 and n0:
                        geom((w[:-1], q ^ (1 << v), n0 - 1, a, b), -self.g / self.N * sqrt(n0))
                    gate = sum(p is not None and p[0] in w[-2:] for p in [a, b])
                    if gate:
                        geom((w[:-1], q, n0, a, b), -EPSILON / sqrt(self.N) * gate)
                else:
                    self.audit['blocked_nonzero_flux_pops'] += 1
            for species, carrier, charge, hopping in [('A', a, 1, JA), ('B', b, 2, JB)]:
                if carrier is None:
                    continue
                position, level = carrier
                k = w.index(position)
                for j in [k - 1, k + 1]:
                    if not 0 <= j < len(w):
                        continue
                    new = (w[j], level)
                    target = (w, q, n0, new, b) if species == 'A' else (w, q, n0, a, new)
                    newfield = self.fields(target)
                    edge = tuple(sorted((position, w[j])))
                    link = pairs(self.N).index(edge)
                    if check:
                        assert field_difference(efield, newfield) == (link,)
                        assert newfield[link] - efield[link] == (-charge if position < w[j] else charge)
                        self.audit['checked_charge_hops'] += 1
                    add(target, col, -hopping, moving=True)
            if a is not None and b is not None and a[0] != b[0] and b[1] in [0, 1]:
                sign = 1 if a[1] == 0 else -1
                add((w, q, n0, a, (b[0], 1 - b[1])), col, KAPPA * sign)
        shape = (len(self.states),) * 2
        H = coo_matrix((values, (rows, cols)), shape=shape).tocsr()
        Hmove = coo_matrix((mv, (mr, mc)), shape=shape).tocsr()
        assert not (H - H.getH()).nnz
        assert not (Hmove - Hmove.getH()).nnz
        return H, Hmove

    def source(self, event):
        rows, cols, data = [], [], []
        for col, (w, q, n0, a, b) in enumerate(self.states):
            if a is not None and a[0] == event:
                rows.append(self.index[(w, q, n0, (event, 1 - a[1]), b)])
                cols.append(col); data.append(.5)
        return coo_matrix((data, (rows, cols)), shape=self.H.shape).tocsr()

    def receiver(self, event, prefix=None):
        rows, cols, data = [], [], []
        for col, (w, q, n0, a, b) in enumerate(self.states):
            if prefix is not None and w[:len(prefix)] != prefix:
                continue
            if b is not None and b[0] == event and b[1] in [0, 1]:
                rows.append(self.index[(w, q, n0, a, (event, 1 - b[1]))])
                cols.append(col); data.append(1j if b[1] == 0 else -1j)
        return coo_matrix((data, (rows, cols)), shape=self.H.shape).tocsr()

    def prepared(self, w, receiver_event):
        q = sum(1 << v for v in w)
        x = np.zeros(len(self.states), complex)
        for n0, fieldphase in [(0, 1), (len(w) - 1, 1j)]:
            for level, sourcephase in [(0, 1), (1, 1j)]:
                x[self.index[(w, q, n0, (w[0], level), (receiver_event, 0))]] += fieldphase * sourcephase / 2
        assert abs(np.linalg.norm(x) - 1) < 1e-14
        return x
