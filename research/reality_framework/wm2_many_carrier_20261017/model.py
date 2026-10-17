"""Exact finite D1 diagnostics: hard-core site slots and tree Gauss isometry.

Small complete-geometry controls use two receiver modes, not the physical
L_N=64N^2 receiver. No formation trajectories. Fixed-path Fock controls are
separately identified and do not substitute for the moving-geometry proof.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from math import sqrt
import numpy as np
from scipy.sparse import coo_matrix

JA, JB, KAPPA, EPSILON = 1 / 64, 1 / 96, 1 / 4, 1 / 4


def species_states(word, number, levels):
    return [tuple(zip(sites, colors)) for sites in combinations(sorted(word), number)
            for colors in product(range(levels), repeat=number)]


def moved(carriers, old, new):
    return tuple(sorted((new if pos == old else pos, level) for pos, level in carriers))


@lru_cache(None)
def pairs(N):
    return tuple(combinations(range(N), 2))


@lru_cache(None)
def electric(N, word, positions_A, positions_B):
    fields = dict.fromkeys(pairs(N), 0)
    charge = {x: positions_A.count(x) + 2 * positions_B.count(x) for x in word}
    remaining = sum(charge.values())
    for left, right in zip(word, word[1:]):
        remaining -= charge[left]
        fields[tuple(sorted((left, right)))] = -remaining if left < right else remaining
    return tuple(fields[edge] for edge in pairs(N))


def constraints(N, word, a, b, fields):
    gauss = [0] * N
    for (u, v), e in zip(pairs(N), fields):
        gauss[u] += e; gauss[v] -= e
    for positions, z in [(a, 1), (b, 2)]:
        for x in positions:
            gauss[x] -= z; gauss[word[0]] += z
    bonds = {tuple(sorted(edge)) for edge in zip(word, word[1:])}
    support = [0 if edge in bonds else e * e for edge, e in zip(pairs(N), fields)]
    return tuple(gauss), tuple(support)


def matter_actions(word, a, b):
    """All physical matter matrix transitions; duplicate coherent terms retained."""
    apos, bpos = {p for p, c in a}, {p for p, c in b}
    for species, slots, occupied, j in [('A', a, apos, JA), ('B', b, bpos, JB)]:
        for position, color in slots:
            place = word.index(position)
            for target_place in [place - 1, place + 1]:
                if not 0 <= target_place < len(word):
                    continue
                target = word[target_place]
                if target in occupied:
                    continue
                new = moved(slots, position, target)
                aa, bb = (new, b) if species == 'A' else (a, new)
                yield aa, bb, -j, ('hop', species, position, target)
    for pos_b, level_b in b:
        if level_b == 2:
            continue
        for pos_a, level_a in a:
            if abs(word.index(pos_a) - word.index(pos_b)) != 1:
                continue
            new_b = tuple((p, 1 - level if p == pos_b else level) for p, level in b)
            yield a, new_b, KAPPA * (1 if level_a == 0 else -1), ('internal',)


class FockChain:
    def __init__(self, length, number_A, number_B):
        self.length = length
        word = tuple(range(length))
        self.states = [(a, b) for a in species_states(word, number_A, 2)
                       for b in species_states(word, number_B, 3)]
        self.index = {state: i for i, state in enumerate(self.states)}
        rows, cols, values = [], [], []
        self.clock = np.zeros(len(self.states))
        for col, (a, b) in enumerate(self.states):
            self.clock[col] = sum(level == 2 for pos, level in b)
            for aa, bb, value, tag in matter_actions(word, a, b):
                rows.append(self.index[(aa, bb)]); cols.append(col); values.append(value)
        shape = (len(self.states),) * 2
        self.V = coo_matrix((values, (rows, cols)), shape=shape).tocsr()
        from scipy.sparse import diags
        self.H = self.V + diags(self.clock)
        assert not (self.H - self.H.getH()).nnz

    def number(self, species, event):
        which = 0 if species == 'A' else 1
        return np.array([float(any(p == event for p, level in state[which])) for state in self.states])

    def internal(self, species, event, matrix):
        which = 0 if species == 'A' else 1
        levels = 2 if species == 'A' else 3
        rows, cols, values = [], [], []
        for col, state in enumerate(self.states):
            slots = state[which]
            for p, old in slots:
                if p != event:
                    continue
                for new in range(levels):
                    value = matrix[new, old]
                    if value == 0:
                        continue
                    changed = tuple((x, new if x == event else c) for x, c in slots)
                    target = (changed, state[1]) if which == 0 else (state[0], changed)
                    rows.append(self.index[target]); cols.append(col); values.append(value)
        return coo_matrix((values, (rows, cols)), shape=self.H.shape).tocsr()


class KernelModel:
    def __init__(self, N, number_A, number_B, g=1.):
        self.N, self.g = N, g
        self.NA, self.NB = number_A, number_B
        self.delta = tuple(.23 + np.linspace(-.25, .25, N))
        self.words = [(0,) + tail for k in range(N) for tail in permutations(range(1, N), k)]
        self.states = []
        for word in self.words:
            for a in species_states(word, number_A, 2):
                for b in species_states(word, number_B, 3):
                    for q in range(1, 1 << N, 2):
                        for n0 in range(q.bit_count()):
                            self.states.append((word, q, n0, a, b))
        self.index = {state: i for i, state in enumerate(self.states)}
        self.audit = {'geometry_entries': 0, 'hopping_entries': 0,
                      'blocked_charged_pops': 0, 'largest_electric_magnitude': 0}
        self.H = self._build()

    def fields(self, state):
        word, q, n0, a, b = state
        return electric(self.N, word, tuple(p for p, c in a), tuple(p for p, c in b))

    def _build(self):
        rows, cols, values = [], [], []
        def add(target, col, value):
            rows.append(self.index[target]); cols.append(col); values.append(value)
        for col, state in enumerate(self.states):
            word, q, n0, a, b = state
            field = self.fields(state)
            self.audit['largest_electric_magnitude'] = max(self.audit['largest_electric_magnitude'], max(map(abs, field), default=0))
            gv, cv = constraints(self.N, word, tuple(p for p, c in a), tuple(p for p, c in b), field)
            assert not any(gv) and not any(cv) and max(map(abs, field), default=0) <= 3 * self.N
            mask = sum(1 << x for x in word)
            diagonal = 16 * (q.bit_count() - len(word)) + 64 * (q ^ mask).bit_count()
            diagonal += sum(self.delta[x] for x in range(self.N) if q >> x & 1) - self.delta[0]
            diagonal += sum(level == 2 for p, level in b)
            add(state, col, diagonal)
            total = q.bit_count() - 1
            for step in [-1, 1]:
                newn = n0 + step
                if 0 <= newn <= total:
                    value = -sqrt(n0 * (total - n0 + 1)) if step == -1 else -sqrt((n0 + 1) * (total - n0))
                    add((word, q, newn, a, b), col, value)
            def geometry(target, value):
                assert self.fields(target) == field, 'geometry may not silently rewrite a charge string'
                self.audit['geometry_entries'] += 1
                add(target, col, value)
            for x in range(self.N):
                if x in word:
                    continue
                if not q >> x & 1:
                    geometry((word + (x,), q | (1 << x), n0 + 1, a, b), -self.g / self.N * sqrt(n0 + 1))
                gate = sum(p == word[-1] for p, level in a + b)
                if gate:
                    geometry((word + (x,), q, n0, a, b), -EPSILON / sqrt(self.N) * gate)
            if len(word) > 1:
                edge = tuple(sorted(word[-2:]))
                if field[pairs(self.N).index(edge)] == 0:
                    x = word[-1]
                    if q >> x & 1 and n0:
                        geometry((word[:-1], q ^ (1 << x), n0 - 1, a, b), -self.g / self.N * sqrt(n0))
                    gate = sum(p in word[-2:] for p, level in a + b)
                    if gate:
                        geometry((word[:-1], q, n0, a, b), -EPSILON / sqrt(self.N) * gate)
                else:
                    self.audit['blocked_charged_pops'] += 1
            for aa, bb, value, tag in matter_actions(word, a, b):
                target = (word, q, n0, aa, bb)
                newfield = self.fields(target)
                if tag[0] == 'hop':
                    edge = tuple(sorted(tag[2:]))
                    link = pairs(self.N).index(edge)
                    changed = tuple(k for k, (x, y) in enumerate(zip(field, newfield)) if x != y)
                    z = 1 if tag[1] == 'A' else 2
                    assert changed == (link,)
                    assert newfield[link] - field[link] == (-z if tag[2] < tag[3] else z)
                    self.audit['hopping_entries'] += 1
                else:
                    assert field == newfield
                add(target, col, value)
        H = coo_matrix((values, (rows, cols)), shape=(len(self.states),) * 2).tocsr()
        assert not (H - H.getH()).nnz
        return H


def common_reference(model, prefix):
    ell = len(prefix)
    inside = np.array([word[:ell] == prefix and all(p in prefix for p, c in a + b)
                       for word, q, n0, a, b in model.states])
    rr, cc, vv = [], [], []
    for row, col, value in zip(model.H.tocoo().row, model.H.tocoo().col, model.H.tocoo().data):
        if inside[row] != inside[col]:
            continue
        if inside[col]:
            wr, qr, nr, ar, br = model.states[row]
            wc, qc, nc, ac, bc = model.states[col]
            if wr != wc and qr == qc:  # every bare geometry entry
                continue
        rr.append(row); cc.append(col); vv.append(value)
    ref = coo_matrix((vv, (rr, cc)), shape=model.H.shape).tocsr()
    assert not (ref - ref.getH()).nnz
    return ref, inside
