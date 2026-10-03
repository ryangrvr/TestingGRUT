"""DA0 / C3-D0 -- architecture detector qualification.  Primary detector Phi (coarse-skeleton recurrence) exactly as
chartered in C3_D0_ARCHITECTURE_DETECTOR.md sec. D0-F4 + REPAIR 01 (sec. R); control constructions P1 and the
epsilon-controls (sec. D0-F2.1 / F2.2); secondary diagnostics (sec. R, D0-O6).  Independent code path, not independent
reviewer.  Detector functions read ONLY a parent array (root has parent -1): no coordinates, no labels, no history.
"""
import math
from collections import deque
import numpy as np


# ------------------------------------------------------------------ rooted-tree primitives (parent array only)
def _structure(parent):
    parent = np.asarray(parent)
    N = len(parent); roots = np.nonzero(parent < 0)[0]
    assert len(roots) == 1, "exactly one root required"
    root = int(roots[0]); kids = [[] for _ in range(N)]
    for v in range(N):
        if parent[v] >= 0:
            kids[int(parent[v])].append(v)
    order = [root]; i = 0
    while i < len(order):
        order.extend(kids[order[i]]); i += 1
    assert len(order) == N, "parent array is not a tree"
    S = np.ones(N, dtype=np.int64)
    for v in reversed(order[1:]):
        S[parent[v]] += S[v]
    return N, root, kids, order, S


def _heavy(kids, S):
    """Heavy child = unique child of maximum S; tie -> none (canonical, label-free)."""
    h = [-1] * len(kids)
    for v, ks in enumerate(kids):
        if ks:
            m = max(S[c] for c in ks)
            arg = [c for c in ks if S[c] == m]
            if len(arg) == 1:
                h[v] = arg[0]
    return h


def _ceil_sqrt(n):
    r = math.isqrt(int(n)); return r if r * r == n else r + 1


class _Interner:
    def __init__(self):
        self.d = {}

    def __call__(self, key):
        return self.d.setdefault(key, len(self.d))


def phi(parent, return_detail=False):
    """Primary detector Phi (sec. D0-F4).  Returns Phi, or (Phi, |U|, #classes) with return_detail."""
    N, root, kids, order, S = _structure(parent); h = _heavy(kids, S)
    thr = _ceil_sqrt(N)
    par = np.asarray(parent)
    U = [v for v in order[1:] if h[int(par[v])] != v and S[v] >= thr]
    intern = _Interner(); classes = set()
    for v in U:
        mu = _ceil_sqrt(S[v])
        # kept node set: upward-closed {w in T_v : S_w >= mu}; collect in pre-order
        kept_order = []; st = [v]
        while st:
            w = st.pop(); kept_order.append(w)
            for c in kids[w]:
                if S[c] >= mu:
                    st.append(c)
        red = {}
        for w in reversed(kept_order):                 # post-order: children before parents
            kc = [c for c in kids[w] if S[c] >= mu]
            if w != v and len(kc) == 1:
                red[w] = red[kc[0]]                    # homeomorphic reduction: suppress unary non-root node
            else:
                red[w] = intern(tuple(sorted(red[c] for c in kc)))
        classes.add(red[v])
    val = 0.0 if len(U) <= 1 else 1.0 - len(classes) / len(U)
    return (val, len(U), len(classes)) if return_detail else val


# ------------------------------------------------------------------ secondary diagnostics (cannot rescue Phi)
def exact_classes(parent):
    """5a: number of distinct rooted-isomorphism classes among all subtrees (minimal-DAG size)."""
    N, root, kids, order, S = _structure(parent); intern = _Interner(); cid = np.zeros(N, dtype=np.int64)
    for v in reversed(order):
        cid[v] = intern(tuple(sorted(int(cid[c]) for c in kids[v])))
    return len(set(cid.tolist()))


def light_depth(parent):
    N, root, kids, order, S = _structure(parent); h = _heavy(kids, S); par = np.asarray(parent)
    ld = np.zeros(N)
    for v in order[1:]:
        p = int(par[v]); ld[v] = ld[p] + (0 if h[p] == v else 1)
    return float(ld.mean() / math.log2(N))


def light_fraction_top(parent):
    """q_v = 1 - S_h(v)/(S_v - 1) (q = 1 if tied maximum) for nodes with S_v in [sqrt N, N)."""
    N, root, kids, order, S = _structure(parent); h = _heavy(kids, S); lo = math.sqrt(N); q = []
    for v in range(N):
        if kids[v] and lo <= S[v] < N:
            q.append(1.0 if h[v] < 0 else 1.0 - S[h[v]] / (S[v] - 1))
    return np.array(q)


def horton_tokunaga(parent):
    """Strahler branches; Horton R_B = geometric mean N_k/N_{k+1}; Tokunaga T_1, T_2 (mean side tributaries of order
    j-1, j-2 per branch of order j), c = T_2/T_1."""
    N, root, kids, order, S = _structure(parent); o = np.ones(N, int)
    for v in reversed(order):
        if kids[v]:
            ks = [o[c] for c in kids[v]]; m = max(ks); o[v] = m + 1 if ks.count(m) >= 2 else m
    par = np.asarray(parent); K = int(o[root])
    # branch heads: root, or nodes whose parent has a different order
    nb = np.zeros(K + 2); side = {}
    for v in order:
        if v == root or o[int(par[v])] != o[v]:
            nb[o[v]] += 1
        if v != root and o[int(par[v])] != o[v]:
            p = int(par[v]); j, i = o[p], o[v]
            merging = (i == j - 1) and not any(o[c] == j for c in kids[p])           # p starts the order-j branch
            if not merging:                                                          # else: side tributary
                side[(i, j)] = side.get((i, j), 0) + 1
    ratios = [nb[k] / nb[k + 1] for k in range(1, K) if nb[k + 1] > 0]
    RB = float(np.exp(np.mean(np.log(ratios)))) if ratios else float('nan')
    def T(k):
        vals = [side.get((j - k, j), 0) / nb[j] for j in range(k + 1, K + 1) if nb[j] > 0]
        return float(np.mean(vals)) if vals else float('nan')
    T1, T2 = T(1), T(2)
    return RB, T1, T2, (T2 / T1 if T1 and T1 > 0 else float('nan')), K


# ------------------------------------------------------------------ control constructions (sec. D0-F2)
def _wilson_block(cells, root, L, rng, parent):
    """Wilson UST on the lattice subgraph induced by `cells` (a rectangle), rooted at `root`; writes into parent."""
    cs = set(cells); intree = {root}; nxt = {}
    def nbrs(u):
        x, y = divmod(u, L); out = []
        for (a, b) in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
            if 0 <= a < L and 0 <= b < L and a * L + b in cs:
                out.append(a * L + b)
        return out
    for s in sorted(cells):
        u = s
        while u not in intree:
            nb = nbrs(u); nxt[u] = nb[rng.integers(len(nb))]; u = nxt[u]
        u = s
        while u not in intree:
            parent[u] = nxt[u]; intree.add(u); u = nxt[u]


def bisection_eps(L, rng, b):
    """P0-epsilon: Stage-A bisection_tree recursion, stopped at the first block with mass <= b; UST inside it."""
    idx = lambda x, y: x * L + y; parent = -np.ones(L * L, int)
    def build(x0, x1, y0, y1, r):
        mass = (x1 - x0) * (y1 - y0)
        if mass == 1:
            return
        if mass <= b:
            _wilson_block([idx(x, y) for x in range(x0, x1) for y in range(y0, y1)], r, L, rng, parent); return
        rx, ry = divmod(r, L)
        if x1 - x0 >= y1 - y0:
            xm = (x0 + x1) // 2; A = (x0, xm, y0, y1) if rx < xm else (xm, x1, y0, y1)
            B = (xm, x1, y0, y1) if rx < xm else (x0, xm, y0, y1); ym = (y0 + y1) // 2
            bcell = idx(xm, ym) if rx < xm else idx(xm - 1, ym); acell = idx(xm - 1, ym) if rx < xm else idx(xm, ym)
        else:
            ym = (y0 + y1) // 2; A = (x0, x1, y0, ym) if ry < ym else (x0, x1, ym, y1)
            B = (x0, x1, ym, y1) if ry < ym else (x0, x1, y0, ym); xm = (x0 + x1) // 2
            bcell = idx(xm, ym) if ry < ym else idx(xm, ym - 1); acell = idx(xm, ym - 1) if ry < ym else idx(xm, ym)
        build(*A, r); build(*B, bcell); parent[bcell] = acell
    build(0, L, 0, L, 0)
    return parent


def ternary_tree(L, rng=None, b=None):
    """P1 (sec. D0-F2.1); with rng and b: P1-epsilon (recursion stopped at first block with mass <= b, UST inside)."""
    k = round(math.log(L, 3)); assert 3 ** k == L
    idx = lambda x, y: x * L + y; parent = -np.ones(L * L, int)
    def build(x0, y0, s, r):
        if s == 1:
            return
        if b is not None and s * s <= b:
            _wilson_block([idx(x, y) for x in range(x0, x0 + s) for y in range(y0, y0 + s)], r, L, rng, parent); return
        t = s // 3; rx, ry = divmod(r, L); qa, qb = (rx - x0) // t, (ry - y0) // t
        seen = {(qa, qb)}; dq = deque([(qa, qb)]); exits = {(qa, qb): r}
        while dq:
            a, bb = dq.popleft()
            for (a2, b2) in sorted([(a - 1, bb), (a + 1, bb), (a, bb - 1), (a, bb + 1)], key=lambda q: 3 * q[0] + q[1]):
                if 0 <= a2 < 3 and 0 <= b2 < 3 and (a2, b2) not in seen:
                    seen.add((a2, b2)); dq.append((a2, b2))
                    cx0, cy0 = x0 + a2 * t, y0 + b2 * t
                    if a2 == a + 1:   e = (cx0, cy0 + t // 2); pc = (cx0 - 1, cy0 + t // 2)
                    elif a2 == a - 1: e = (cx0 + t - 1, cy0 + t // 2); pc = (cx0 + t, cy0 + t // 2)
                    elif b2 == bb + 1: e = (cx0 + t // 2, cy0); pc = (cx0 + t // 2, cy0 - 1)
                    else:             e = (cx0 + t // 2, cy0 + t - 1); pc = (cx0 + t // 2, cy0 + t)
                    parent[idx(*e)] = idx(*pc); exits[(a2, b2)] = idx(*e)
        for (a, bb), ex in sorted(exits.items(), key=lambda kv: 3 * kv[0][0] + kv[0][1]):
            build(x0 + a * t, y0 + bb * t, t, ex)
    build(0, 0, L, 0)
    return parent


def icbrt_floor(n):
    c = int(round(n ** (1 / 3)))
    while c ** 3 > n: c -= 1
    while (c + 1) ** 3 <= n: c += 1
    return c


def is_lattice_spanning_tree(parent, L):
    N = L * L; parent = np.asarray(parent)
    if len(parent) != N or parent[0] != -1 or np.sum(parent < 0) != 1:
        return False
    for v in range(1, N):
        p = int(parent[v]); (x, y), (a, b) = divmod(v, L), divmod(p, L)
        if abs(x - a) + abs(y - b) != 1:
            return False
    try:
        _structure(parent); return True
    except AssertionError:
        return False
