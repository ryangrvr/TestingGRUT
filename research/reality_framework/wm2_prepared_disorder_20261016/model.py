"""Independent sparse construction of S3 with frozen on-event detuning.

Root is event 0 in each enumerated component. Other roots are handled by
permuting the entire detuning assignment, never by choosing a better root.
No occupation cutoff is imposed: the exact charge fixes n_R=sum(q)-1.
"""
from itertools import combinations_with_replacement
from math import factorial, sqrt

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import expm_multiply

EPSILON = J = 0.25
NU = 2.1


def rooted_words(n):
    out = []

    def visit(w):
        out.append(w)
        for v in range(n):
            if v not in w:
                visit(w + (v,))

    visit((0,))
    return out


def occupation_moves(o, modes):
    """o is a sorted tuple of occupied mode labels, including repeats."""
    for r in set(o):
        count = o.count(r)
        for s in (r - 1, r + 1):
            if 0 <= s < modes:
                a = list(o)
                a.remove(r)
                a.append(s)
                yield tuple(sorted(a)), -sqrt(count * (o.count(s) + 1))


def vacuum(n, modes, delta, g=1.0, classes=False):
    """Complete probe-vacuum component, optionally quotienting equal deltas.

    Equal-type words retain their entire type ordering. The quotient is
    an invariant isometry, not a radial or incoherent approximation.
    """
    if classes:
        values = tuple(sorted(set(delta[1:])))
        counts = tuple(sum(x == a for x in delta[1:]) for a in values)
        words = []

        def visit(w):
            words.append(w)
            for a, capacity in enumerate(counts):
                if w.count(a) < capacity:
                    visit(w + (a,))

        visit(())
    else:
        words = rooted_words(n)
        values, counts = tuple(delta), None

    occupations = [list(combinations_with_replacement(range(modes), k))
                   for k in range(n)]
    states = [(w, o) for w in words
              for o in occupations[len(w) if classes else len(w) - 1]]
    index = {s: j for j, s in enumerate(states)}
    row, col, data = [], [], []
    erow, ecol, edata = [], [], []

    def add(i, j, a):
        row.append(i)
        col.append(j)
        data.append(a)

    for j, (w, o) in enumerate(states):
        add(j, j, sum(values[a] for a in (w if classes else w[1:])))
        for oo, a in occupation_moves(o, modes):
            add(index[w, oo], j, a)
        if classes:
            candidates = [(a, sqrt(capacity - w.count(a)))
                          for a, capacity in enumerate(counts)
                          if w.count(a) < capacity]
        else:
            candidates = [(v, 1.0) for v in range(n) if v not in w]
        for v, factor in candidates:
            oo = (0,) + o
            i = index[w + (v,), oo]
            a = -g / n * factor * sqrt(o.count(0) + 1)
            add(i, j, a)
            add(j, i, a)
            erow.extend((i, j))
            ecol.extend((j, i))
            edata.extend((a, a))
    h = coo_matrix((data, (row, col)), shape=(len(states), len(states))).tocsr()
    emission = coo_matrix((edata, (erow, ecol)), shape=h.shape).tocsr()
    seed = index[((), ()) if classes else ((0,), ())]
    return dict(H=h, states=states, index=index, seed=seed, n=n, modes=modes,
                classes=classes, values=values, counts=counts, delta=tuple(delta), g=g,
                H_emission=emission)


def full_probe(n, modes, delta, g=1.0):
    states = []
    for w in rooted_words(n):
        for q in range(1, 1 << n, 2):
            for o in combinations_with_replacement(range(modes), q.bit_count() - 1):
                for p in range(n):
                    states.append((w, q, o, p))
    index = {s: j for j, s in enumerate(states)}
    rows, cols, data, support = [], [], [], []
    consistent = []

    def add(j, target, a, kind):
        i = index[target]
        rows.append(i)
        cols.append(j)
        data.append(a)
        if i != j:
            support.append((j, i, 0 if kind == "receiver" else 1, kind))

    for j, (w, q, o, p) in enumerate(states):
        mask = sum(1 << v for v in w)
        if q == mask:
            consistent.append(j)
        diagonal = (16 * (q.bit_count() - len(w))
                    + 64 * (q ^ mask).bit_count()
                    + sum(delta[v] * ((q >> v & 1) - int(v == 0)) for v in range(n)))
        add(j, (w, q, o, p), diagonal, "diagonal")
        for oo, a in occupation_moves(o, modes):
            add(j, (w, q, oo, p), a, "receiver")
        u = w[-1]
        for v in range(n):
            if v in w:
                continue
            if not (q >> v & 1):
                add(j, (w + (v,), q | 1 << v, (0,) + o, p),
                    -g / n * sqrt(o.count(0) + 1), "emission_append")
            if p in (u, v):
                add(j, (w + (v,), q, o, p), -EPSILON / sqrt(n), "bare_append")
        if len(w) > 1:
            v, parent = w[-1], w[-2]
            if q >> v & 1 and 0 in o:
                a = list(o)
                a.remove(0)
                add(j, (w[:-1], q ^ 1 << v, tuple(a), p),
                    -g / n * sqrt(o.count(0)), "emission_pop")
            if p in (parent, v):
                add(j, (w[:-1], q, o, p), -EPSILON / sqrt(n), "bare_pop")
        if p in w:
            k = w.index(p)
            for l in (k - 1, k + 1):
                if 0 <= l < len(w):
                    add(j, (w, q, o, w[l]), -J, "matter")
    h = coo_matrix((data, (rows, cols)), shape=(len(states), len(states))).tocsr()
    p = np.asarray(consistent, dtype=int)
    q = np.setdiff1d(np.arange(len(states)), p)
    return dict(H=h, states=states, index=index, P=p, Q=q, support=support,
                n=n, modes=modes, delta=tuple(delta), g=g)


def evolve(h, x, time):
    a = -1j * time * h
    return expm_multiply(a, x, traceA=a.diagonal().sum())


def norm2(x):
    return float(np.vdot(x, x).real)


def unit(size, j):
    x = np.zeros(size, dtype=complex)
    x[j] = 1
    return x


def multinomial_lift(model, z):
    """Only valid for one detuning class; used as a symmetry control."""
    n, modes = model["n"], model["modes"]
    m = n - 1
    x = np.zeros(len(model["states"]), dtype=complex)
    for j, (w, o) in enumerate(model["states"]):
        k = len(o)
        multiplicities = [m - k] + [o.count(r) for r in range(modes)]
        coefficient = factorial(m)
        for a in multiplicities:
            coefficient /= factorial(a)
        value = sqrt(coefficient) * z[0] ** (m - k)
        for r in o:
            value *= z[r + 1]
        if not model["classes"]:
            value /= sqrt(factorial(m) / factorial(m - k))
        x[j] = value
    return x


def receiver_cut_error(n, modes, time, mu=1.0):
    """Norm difference between the cut and physical longer receiver.

    Applies to vacuum preparation, every allowed detuning pattern and all
    reverse processes; see MATHEMATICS.md. Not a roundoff certificate.
    """
    c = 2 * np.sinh(mu)
    return float(sqrt(n - 1) / c * np.expm1(c * time)
                 / sqrt(np.expm1(2 * mu * (modes - 1))))
