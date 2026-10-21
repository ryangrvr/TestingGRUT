"""Independent periodic-winding reconstruction from a final solver state (q,p of 3 fields, C-order N^3 f64).

Algorithm (different from Codex's cell-BFS with unwrapped coordinates):
  1. masks: plus = (F2 < f2cut) & (S > scut), minus = (F2 < f2cut) & (S < -scut), S=q2/a, F2=(q0^2+q1^2)/a^2
  2. scipy.ndimage.label with 6-connectivity on the NON-periodic box -> open-box clusters
  3. periodic gluing: for each axis, every (index N-1, index 0) face pair with both cells in the mask is a
     graph edge between their open-box labels carrying the lattice translation +e_axis
  4. BFS over the label graph assigning each label an image offset; an edge whose endpoint offsets are
     inconsistent closes a non-contractible cycle with integer winding vector w; the component's winding
     bitmask is the OR of axes with w_axis != 0 (same bit convention as Codex: 1=x, 2=y, 4=z).
Reports components, occupied/largest fractions, any_winding and largest_winding, for the packet threshold
(|S|>0.5, F2<0.5) and stricter thresholds.
Usage: python3 -I topo.py STATE N a OUT.json [SOLVER_TOPOLOGY.jsonl]
If the solver's own topology JSONL is given, its LAST record (same state) is compared field by field.
"""
import json
import sys
from collections import deque

import numpy as np
from scipy import ndimage


def periodic_components(mask):
    N = mask.shape[0]
    lab, n = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    edges = [[] for _ in range(n + 1)]
    for ax in range(3):
        last = np.take(lab, N - 1, axis=ax)
        first = np.take(lab, 0, axis=ax)
        both = (last > 0) & (first > 0)
        pairs = np.unique(np.stack([last[both], first[both]], axis=1), axis=0)
        off = np.zeros(3, dtype=int)
        off[ax] = 1
        for u, v in pairs:
            edges[u].append((v, off))
            edges[v].append((u, -off))
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    shift = [None] * (n + 1)
    comps = []
    for s in range(1, n + 1):
        if shift[s] is not None:
            continue
        shift[s] = np.zeros(3, dtype=int)
        q = deque([s])
        members = []
        wind = 0
        while q:
            u = q.popleft()
            members.append(u)
            for v, o in edges[u]:
                target = shift[u] + o
                if shift[v] is None:
                    shift[v] = target
                    q.append(v)
                else:
                    w = target - shift[v]
                    for ax in range(3):
                        if w[ax] != 0:
                            wind |= 1 << ax
        comps.append((int(sizes[members].sum()), wind))
    return comps


def summarize(mask, v):
    comps = periodic_components(mask)
    if not comps:
        return {'components': 0, 'occupied_fraction': 0.0, 'largest_fraction': 0.0, 'any_winding': 0, 'largest_winding': 0}
    occ = sum(c[0] for c in comps)
    big = max(comps, key=lambda c: c[0])
    anyw = 0
    for c in comps:
        anyw |= c[1]
    return {'components': len(comps), 'occupied_fraction': occ / v, 'largest_fraction': big[0] / v,
            'any_winding': anyw, 'largest_winding': big[1],
            'winding_components_fraction': sum(c[0] for c in comps if c[1] == 7) / v}


def main():
    path, N, a, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    st = np.memmap(path, dtype='<f8', mode='r', shape=(2, 3, N, N, N))
    q = st[0]
    S = np.asarray(q[2]) / a
    F2 = (np.asarray(q[0]) ** 2 + np.asarray(q[1]) ** 2) / a ** 2
    v = N ** 3
    res = {'state': path, 'N': N, 'a': a, 'raw_fplus_05': float((S > .5).mean()), 'raw_fminus_05': float((S < -.5).mean())}
    for scut, f2cut in [(.5, .5), (.7, .3), (.8, .25), (.9, .1)]:
        key = f'S{scut}_F2{f2cut}'
        res[key] = {'plus': summarize((F2 < f2cut) & (S > scut), v), 'minus': summarize((F2 < f2cut) & (S < -scut), v)}
    # Charge in each region (ell = q1 p0 - q0 p1 in conformal variables)
    p = st[1]
    ell = np.asarray(q[1]) * np.asarray(p[0]) - np.asarray(q[0]) * np.asarray(p[1])
    Q = ell.sum()
    res['Q_sum_cells'] = float(Q)
    res['q_fraction_F2_lt_0.5'] = float(ell[F2 < .5].sum() / Q)
    res['q_fraction_absS_lt_0.5'] = float(ell[np.abs(S) < .5].sum() / Q)
    if len(sys.argv) > 5:
        last = [json.loads(l) for l in open(sys.argv[5]) if l.strip()][-1]
        keys = ('components', 'occupied_fraction', 'largest_fraction', 'any_winding', 'largest_winding')
        res['solver_record_match'] = all(res['S0.5_F20.5'][s][k] == last[s][k] for s in ('plus', 'minus') for k in keys)
        res['solver_record_a'] = last['a']
    open(out, 'w').write(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
