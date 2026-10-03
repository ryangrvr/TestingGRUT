"""C3-D0 implementation-firewall unit tests (sec. R, D0-O6).  Uses SYNTHETIC trees only for invariance tests; control
trees appear only in (c) the analytic N1 test and (e) construction-validity checks, where no P0/P1/N0/epsilon Phi value
is computed or printed.  Independent code path, not independent reviewer."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d0_detector import phi, exact_classes, ternary_tree, bisection_eps, icbrt_floor, is_lattice_spanning_tree
from pm2_stageA import lattice, sp_tree, bisection_tree, wilson_ust

fails = 0
def check(name, cond):
    global fails
    print(("PASS " if cond else "FAIL ") + name, flush=True); fails += (not cond)


def relabel(parent, perm):
    """New label perm[v] for old node v."""
    parent = np.asarray(parent); new = -np.ones(len(parent), int)
    for v, p in enumerate(parent):
        new[perm[v]] = -1 if p < 0 else perm[p]
    return new


def random_recursive(n, rng):
    p = -np.ones(n, int)
    for v in range(1, n):
        p[v] = rng.integers(v)
    return p


def complete_kary(k, depth):
    p = [-1]; frontier = [0]
    for _ in range(depth):
        nf = []
        for u in frontier:
            for _ in range(k):
                p.append(u); nf.append(len(p) - 1)
        frontier = nf
    return np.array(p)


def path(n):
    return np.array([-1] + list(range(n - 1)))


rng = np.random.default_rng(12345)
# (a) relabelling invariance, (b) rooted-isomorphism invariance (children order is label order; relabelling permutes it)
syn = [("random recursive n=3000", random_recursive(3000, rng)),
       ("complete binary depth 10", complete_kary(2, 10)), ("complete ternary depth 6", complete_kary(3, 6))]
for name, t in syn:
    base = phi(t, True); base_ex = exact_classes(t)
    ok = True
    for r in range(5):
        perm = rng.permutation(len(t)); t2 = relabel(t, perm)
        ok &= phi(t2, True) == base and exact_classes(t2) == base_ex
    check(f"(a,b) relabelling / child-order invariance: {name} (Phi detail {tuple(round(x, 6) if isinstance(x, float) else x for x in base)})", ok)
# (b) isomorphic trees built differently: graft two copies of a random tree under a new root in both orders
t = random_recursive(800, rng); n = len(t)
def graft(a, b):
    p = -np.ones(1 + len(a) + len(b), int)
    for i, q in enumerate(a): p[1 + i] = 0 if q < 0 else 1 + q
    for i, q in enumerate(b): p[1 + len(a) + i] = 0 if q < 0 else 1 + len(a) + q
    return p
u = random_recursive(500, rng)
check("(b) isomorphic trees built in different orders give equal Phi", phi(graft(t, u), True) == phi(graft(u, t), True))
# (c) analytic N1: Phi = 0 with |U| = 0 on every grid L
for L in [16, 32, 64, 128, 27, 81, 243]:
    N, E = lattice(L); par, _ = sp_tree(N, E); par = np.asarray(par)
    v, nU, nc = phi(par, True)
    check(f"(c) N1 L={L}: Phi = 0 and |U| = 0 (LEMMA D0-N1)", v == 0.0 and nU == 0)
    check(f"(c) N1 L={L}: exact classes = 2L-1 (LEMMA D0-N1)", exact_classes(par) == 2 * L - 1)
# (d) analytic hand cases
check("(d) path n=500: Phi = 0", phi(path(500)) == 0.0)
v, nU, nc = phi(complete_kary(2, 6), True)
check(f"(d) complete binary depth 6: Phi = 6/7 with |U| = 14, 2 classes (got {v:.6f}, {nU}, {nc})",
      abs(v - 6 / 7) < 1e-12 and nU == 14 and nc == 2)
# (e) construction validity (no Phi)
for L in [16, 32, 64, 128]:
    check(f"(e) P0 L={L} lattice spanning tree", is_lattice_spanning_tree(bisection_tree(L)[0], L))
    b = icbrt_floor(L * L)
    check(f"(e) P0-eps L={L} (b={b}) seeds 1-3 lattice spanning trees",
          all(is_lattice_spanning_tree(bisection_eps(L, np.random.default_rng(s), b), L) for s in (1, 2, 3)))
for L in [27, 81, 243]:
    check(f"(e) P1 L={L} lattice spanning tree", is_lattice_spanning_tree(ternary_tree(L), L))
    b = icbrt_floor(L * L)
    check(f"(e) P1-eps L={L} (b={b}) seeds 1-3 lattice spanning trees",
          all(is_lattice_spanning_tree(ternary_tree(L, np.random.default_rng(s), b), L) for s in (1, 2, 3)))
# (e') P0-eps with b < 2 reproduces P0 exactly (recursion code path identical to the Stage-A construction)
check("(e') P0-eps with b=1 equals Stage-A bisection_tree at L=32",
      np.array_equal(bisection_eps(32, np.random.default_rng(0), 1), np.asarray(bisection_tree(32)[0])))
# (e'') N0 Wilson trees are valid lattice spanning trees at one size (construction only)
N, E = lattice(27)
check("(e'') N0 L=27 seed 1 lattice spanning tree", is_lattice_spanning_tree(np.asarray(wilson_ust(N, E, np.random.default_rng(1))[0]), 27))
print("TOTAL FAILURES:", fails)
sys.exit(1 if fails else 0)
