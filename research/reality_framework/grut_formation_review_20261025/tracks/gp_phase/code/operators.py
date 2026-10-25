#!/usr/bin/env python3
"""gp_phase (b): enumerate non-derivative gauge-invariant NON-HERMITIAN scalar operators of the
7737487 sector (charges X=1, Y=3, Sigma=6) up to a given mass dimension, and count the physical
(rephasing-invariant) CP phases.

A monomial X^a X*^a' Y^b Y*^b' S^c S*^c' (S = Sigma) times neutral singlets has rephasing vector
n = (a-a', b-b', c-c'); gauge invariance requires n_X + 3 n_Y + 6 n_S = 0. Its coefficient picks up
exp(-i n.alpha) under field rephasings alpha. Non-Hermitian means n != 0; each operator is paired
with its conjugate (-n). With N independent complex couplings whose vectors span a rank-r lattice,
the number of physical phases is N - r (the gauge direction (1,3,6) acts trivially on all of them).
Neutral singlets: with Candidate 1's Z2_phi x Z2_D: phi^2, D^2, H^dag H, curvature R (all dim 2).
Variant: Z2_phi broken (e.g. by the proposed wall bias): phi (dim 1) also allowed.
Derivative operators are omitted (reducible to these by integration by parts and field equations
at the order considered; they do not change the counting of independent phase classes).
"""
import itertools, json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / 'results'
Q = np.array([1, 3, 6])


def monomials(dmax):
    """Field monomials (a,a',b,b',c,c') with dim = sum, gauge invariant, n != 0, one per conjugate pair."""
    res = []
    for exps in itertools.product(range(dmax + 1), repeat=6):
        d = sum(exps)
        if d == 0 or d > dmax:
            continue
        a, ac, b, bc, c, cc = exps
        n = np.array([a - ac, b - bc, c - cc])
        if not n.any() or int(n @ Q) != 0:
            continue
        # canonical representative of the conjugate pair: first nonzero component of n positive... use
        # lexicographic comparison of exps with its conjugate
        conj = (ac, a, bc, b, cc, c)
        if conj < exps:
            continue
        res.append((d, exps, tuple(int(x) for x in n)))
    return res


def neutral_products(dmax, z2_phi=True):
    """Multisets of neutral singlets with total dimension <= dmax (no field factors)."""
    singles = ([('phi^2', 2)] if z2_phi else [('phi', 1)]) + [('D^2', 2), ('HdagH', 2), ('R', 2)]
    out = [((), 0)]
    frontier = [((), 0, 0)]
    while frontier:
        new = []
        for combo, d, start in frontier:
            for i in range(start, len(singles)):
                nm, ds = singles[i]
                if d + ds <= dmax:
                    c2 = combo + (nm,)
                    out.append((c2, d + ds)); new.append((c2, d + ds, i))
        frontier = new
    return out


def name(exps):
    a, ac, b, bc, c, cc = exps
    parts = []
    for sym, p in (('X', a), ('X*', ac), ('Y', b), ('Y*', bc), ('S', c), ('S*', cc)):
        if p:
            parts.append(sym + ('^%d' % p if p > 1 else ''))
    return ' '.join(parts)


def count(dmax, z2_phi=True):
    mons = monomials(dmax)
    neut = neutral_products(dmax, z2_phi)
    ops = []
    for d, exps, n in mons:
        for combo, dn in neut:
            if d + dn <= dmax:
                ops.append(dict(dim=d + dn, field_part=name(exps), singlets=list(combo), n=n))
    ops.sort(key=lambda o: (o['dim'], o['field_part']))
    table = []
    for D in range(1, dmax + 1):
        sel = [o for o in ops if o['dim'] <= D]
        if not sel:
            table.append(dict(dim_max=D, complex_couplings=0, rank=0, physical_phases=0)); continue
        M = np.array([o['n'] for o in sel])
        r = int(np.linalg.matrix_rank(M))
        table.append(dict(dim_max=D, complex_couplings=len(sel), rank=r, physical_phases=len(sel) - r,
                          new_at_this_dim=[o['field_part'] + (' * ' + ' '.join(o['singlets']) if o['singlets'] else '') for o in ops if o['dim'] == D]))
    return table


def lattice_class(n):
    """Write n = A u + B w, u = (3,-1,0) [Y* X^3], w = (0,2,-1) [S* Y^2]; unimodular basis."""
    nX, nY, nS = n
    B = -nS; A = -nY + 2 * B
    assert (A * 3, -A + 2 * B, -B) == (nX, nY, nS)
    return A, B


def main():
    out = dict(note='Non-derivative operators; counts are N complex couplings, rank of their rephasing vectors, physical phases N - rank.')
    for z2 in (True, False):
        tab = count(7, z2)
        out['Z2_phi_imposed' if z2 else 'Z2_phi_broken'] = tab
        print('Z2_phi imposed' if z2 else 'Z2_phi broken (phi allowed linearly)')
        for row in tab:
            new = row.get('new_at_this_dim', [])
            print('  dim<=%d: N=%d rank=%d phases=%d  new: %s' % (row['dim_max'], row['complex_couplings'], row['rank'], row['physical_phases'],
                                                              (', '.join(new[:12]) + (' ...(+%d)' % (len(new) - 12) if len(new) > 12 else ''))))
    examples = {'Y* X^3': (3, -1, 0), 'S* Y^2': (0, 2, -1), 'S* Y X^3': (3, 1, -1), 'S* X^6': (6, 0, -1), 'S* Y^3 X*^3': (-3, 3, -1),
                'S*^2 Y^4': (0, 4, -2), 'Y*^2 X^6': (6, -2, 0)}
    out['phase_classes'] = {k: dict(n=v, kappa_power=lattice_class(v)[0], mu_power=lattice_class(v)[1],
                                    invariant='arg c - %d arg kappa - %d arg mu' % lattice_class(v)) for k, v in examples.items()}
    for k, v in out['phase_classes'].items():
        print(k, v['invariant'])
    (OUT / 'operators.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
