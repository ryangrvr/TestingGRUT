"""Full complex interaction-picture word images on a common state space.

Two receiver modes are finite diagnostic instances only. No receiver state
is truncated inside these instances. Reachable images are constructed by
exact transition rules; no amplitude threshold or history quotient is used.
"""
from collections import defaultdict
from functools import lru_cache
from math import sqrt
import numpy as np


class Comparison:
    def __init__(self, N, m, ell):
        self.N, self.m, self.ell = N, m, ell
        self.prefix = tuple(range(ell))
        self.g, self.J, self.epsilon = 1.03, .25, .25
        self.delta = np.linspace(-.73, .71, N)

    # state = (word, storage mask, 'P' or 'R', event or dummy position, contact occupation)
    def initial(self):
        w = tuple(range(self.m)); q = (1 << self.m) - 1
        return {(w, q, 'P', p, n): phase / 2
                for p in [0, 1] for n, phase in [(0, 1), (self.m - 1, 1j)]}

    @lru_cache(None)
    def field_unitary(self, total, t):
        h = np.zeros((total + 1, total + 1))
        for n in range(total):
            h[n + 1, n] = h[n, n + 1] = -sqrt((n + 1) * (total - n))
        e, v = np.linalg.eigh(h)
        return (v * np.exp(-1j * e * t)) @ v.T

    def D(self, state, t):
        groups = {}
        for (w, q, kind, p, n), z in state.items():
            key = (w, q, kind, p)
            if key not in groups:
                groups[key] = np.zeros(q.bit_count(), complex)
            groups[key][n] += z
        out = {}
        for (w, q, kind, p), x in groups.items():
            mask = sum(1 << v for v in w)
            energy = 16 * (q.bit_count() - len(w)) + 64 * (q ^ mask).bit_count()
            energy += sum(self.delta[v] for v in range(self.N) if q >> v & 1) - self.delta[w[0]]
            y = np.exp(-1j * energy * t) * (self.field_unitary(q.bit_count() - 1, t) @ x)
            for n, z in enumerate(y):
                if z != 0:
                    out[(w, q, kind, p, n)] = z
        return out

    def reference_index(self, kind, p):
        return p if kind == 'R' else self.prefix.index(p)

    def key_at(self, w, q, x, n):
        return (w, q, 'P', self.prefix[x], n) if x < self.ell else (w, q, 'R', x, n)

    def neighbours(self, state, reference=False, broken=False):
        w, q, kind, p, n = state
        if reference:
            if w[:self.ell] != self.prefix or q != sum(1 << v for v in w):
                return []
            if kind == 'P' and p not in self.prefix:
                return []
        elif kind == 'R':
            return []
        out = []
        for v in range(self.N):
            if v in w:
                continue
            neww = w + (v,)
            if not q >> v & 1:
                out.append(((neww, q | (1 << v), kind, p, n + 1), -self.g / self.N * sqrt(n + 1)))
            if not reference and p in [w[-1], v]:
                out.append(((neww, q, kind, p, n), -self.epsilon / sqrt(self.N)))
        if len(w) > 1:
            v = w[-1]
            if q >> v & 1 and n and (not reference or len(w) > self.ell):
                out.append(((w[:-1], q ^ (1 << v), kind, p, n - 1), -self.g / self.N * sqrt(n)))
            if not reference and p in [w[-1], w[-2]]:
                out.append(((w[:-1], q, kind, p, n), -self.epsilon / sqrt(self.N)))
        if reference:
            x = self.reference_index(kind, p)
            for y in [x - 1, x + 1]:
                if y < 0 or (broken and {x, y} == {0, 1}):
                    continue
                out.append((self.key_at(w, q, y, n), -self.J))
        elif p in w:
            x = w.index(p)
            for y in [x - 1, x + 1]:
                if 0 <= y < len(w):
                    out.append(((w, q, kind, w[y], n), -self.J))
        return out

    def V(self, state, reference=False, broken=False):
        out = defaultdict(complex)
        for key, z in state.items():
            for target, weight in self.neighbours(key, reference, broken):
                out[target] += z * weight
        return dict(out)

    def interaction(self, state, t, reference=False, broken=False):
        return self.D(self.V(self.D(state, t), reference, broken), -t)


def difference(a, b):
    return sqrt(sum(abs(a.get(k, 0) - b.get(k, 0))**2 for k in a.keys() | b.keys()))


def run():
    rows = []
    for N, m, ell in [(12, 6, 3), (16, 8, 4), (12, 10, 5)]:
        model = Comparison(N, m, ell)
        p = r = model.initial()
        errors, sizes = [], []
        for n in range(1, ell):
            t = [.13, .41, .82, 1.19][n - 1]
            p = model.interaction(p, t)
            r = model.interaction(r, t, reference=True)
            err = difference(p, r)
            errors.append(err); sizes.append([len(p), len(r)])
            if n < ell - 1:
                assert err < 2e-12
            else:
                assert err > 1e-6
        bad = model.interaction(model.initial(), .13, reference=True, broken=True)
        good = model.interaction(model.initial(), .13)
        assert difference(good, bad) > .05
        rows.append(dict(N=N, initial_word_length=m, prefix=ell, S9_q=ell - 2,
                         actual_S9_size_sequence=(ell == N // 4 and m >= (N + 1) // 2),
                         complex_word_image_errors=errors, image_sizes=sizes,
                         first_difference_order=ell - 1,
                         shared_edge_deletion_detected_at_order_one=True,
                         receiver_modes=2, contact_occupations='all in each charge sector'))
    return dict(status='PASS', tests=rows,
                scope='amplitude/phase diagnostic for full changing interactions and exact configuration-preserving field evolution; no formation trajectory')
