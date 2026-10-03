"""RA0 / G1 -- finite reflexive autonomy: exhaustive checks of the G1 propositions.

Independent code path, not independent reviewer.  NUMERICAL ILLUSTRATION of propositions proved in
G1_FINITE_REFLEXIVE_AUTONOMY.md; it adjudicates nothing on its own.

Order: Pi <= Pi' iff Pi refines Pi'.  BOT = discrete, TOP = indiscrete.
Version B (strict future):  x ~+ x'  iff  law(Y_{t+1}, Y_{t+2}, ...) agree, computed exactly via the span of
                            word-probability vectors  P D_{b1} P D_{b2} ... P D_{bk} 1.
Version A (present-incl.):  F_A(Pi) = Pi  meet  F_B(Pi).
"""
import itertools
import numpy as np

rng = np.random.default_rng(2026)
TOL = 1e-9


def partitions(n):
    def rec(i, lab, m):
        if i == n:
            yield tuple(lab); return
        for b in range(m + 1):
            lab.append(b); yield from rec(i + 1, lab, max(m, b + 1)); lab.pop()
    yield from rec(0, [], 0)


def canon(lab):
    mp, out = {}, []
    for x in lab:
        mp.setdefault(x, len(mp)); out.append(mp[x])
    return tuple(out)


def classes_from_rows(M):
    lab, reps = [], []
    for x in range(M.shape[0]):
        for k, r in enumerate(reps):
            if np.allclose(M[x], M[r], atol=1e-8):
                lab.append(k); break
        else:
            reps.append(x); lab.append(len(reps) - 1)
    return tuple(lab)


def F_B(P, lab):
    n = len(lab); blocks = sorted(set(lab))
    ops = [P @ np.diag([1.0 if lab[i] == b else 0.0 for i in range(n)]) for b in blocks]
    basis = np.zeros((n, 0)); frontier = [np.ones(n)]
    while frontier:
        new = []
        for v in frontier:
            for O in ops:
                w = O @ v
                M = np.column_stack([basis, w])
                if np.linalg.matrix_rank(M, TOL) > basis.shape[1]:
                    basis = M; new.append(w)
        frontier = new
    return canon(classes_from_rows(basis))


def meet(a, b):
    return canon(tuple(zip(a, b)))


def F_A(P, lab):
    return meet(canon(lab), F_B(P, lab))


def lumpable(P, lab):
    """Kemeny-Snell: P(x, B_j) equal for all x in the same block, all j."""
    lab = canon(lab); k = max(lab) + 1
    S = np.column_stack([P[:, [i for i in range(len(lab)) if lab[i] == j]].sum(1) for j in range(k)])
    for b in range(k):
        rows = S[[i for i in range(len(lab)) if lab[i] == b]]
        if not np.allclose(rows, rows[0], atol=1e-9):
            return False
    return True


def lumped(P, lab):
    lab = canon(lab); k = max(lab) + 1
    S = np.column_stack([P[:, [i for i in range(len(lab)) if lab[i] == j]].sum(1) for j in range(k)])
    return np.array([S[lab.index(b)] for b in range(k)])


def distinct_rows(Q):
    return len(set(classes_from_rows(Q))) == Q.shape[0]


def refines(a, b):
    n = len(a)
    return all((a[i] == a[j]) <= (b[i] == b[j]) for i in range(n) for j in range(n))


def fixed(P, F):
    return [l for l in partitions(P.shape[0]) if F(P, l) == canon(l)]


def dirichlet(n, zeros=0.0):
    P = rng.dirichlet(np.ones(n), size=n)
    if zeros:
        P = P * (rng.random((n, n)) > zeros); P[np.arange(n), rng.integers(0, n, n)] += 1e-3
        P /= P.sum(1, keepdims=True)
    return P


def check_characterization(P):
    """G1.1: Fix_A == lumpable; Fix_B == lumpable & distinct lumped rows. Returns mismatches."""
    bad = 0
    for l in partitions(P.shape[0]):
        a = F_A(P, l) == canon(l); b = F_B(P, l) == canon(l)
        L = lumpable(P, l)
        bad += (a != L) + (b != (L and distinct_rows(lumped(P, l))))
    return bad


if __name__ == "__main__":
    n = 5
    BOT, TOP = tuple(range(n)), (0,) * n
    print("== G1.1 exact characterization (exhaustive over all 52 partitions of 5 states) ==")
    families = {
        "random Dirichlet (20)": [dirichlet(n) for _ in range(20)],
        "sparse random (20)": [dirichlet(n, zeros=0.5) for _ in range(20)],
        "all 120 permutations": [np.eye(n)[list(p)] for p in itertools.permutations(range(n))],
        "iid chain (all rows equal)": [np.tile(rng.dirichlet(np.ones(n)), (n, 1))],
        "identity": [np.eye(n)],
    }
    for name, Ps in families.items():
        print(f"  {name:30s} mismatches vs (A: Kemeny-Snell) / (B: KS + distinct lumped rows): "
              f"{sum(check_characterization(P) for P in Ps)}")

    print("== G1.2 monotonicity (all refining pairs, both versions) ==")
    parts = list(partitions(n)); pairs = [(a, b) for a in parts for b in parts if refines(a, b)]
    for name, P in [("random", dirichlet(n)), ("sparse", dirichlet(n, 0.5)), ("permutation", np.eye(n)[[1, 2, 0, 4, 3]])]:
        vA = sum(not refines(F_A(P, a), F_A(P, b)) for a, b in pairs)
        vB = sum(not refines(F_B(P, a), F_B(P, b)) for a, b in pairs)
        print(f"  {name:12s} refining pairs {len(pairs)}  violations A {vA}  B {vB}")

    print("== G1.2 endpoints ==")
    for name, P in [("random", dirichlet(n)), ("iid", np.tile(rng.dirichlet(np.ones(n)), (n, 1))),
                    ("two equal rows", (lambda Q: (Q.__setitem__(1, Q[0]), Q)[1])(dirichlet(n)))]:
        fA, fB = fixed(P, F_A), fixed(P, F_B)
        print(f"  {name:15s} A: TOP {TOP in fA} BOT {BOT in fA} |Fix_A|={len(fA)}   "
              f"B: TOP {TOP in fB} BOT {BOT in fB} |Fix_B|={len(fB)}")

    print("== G1.3 generic no-go: 1000 random Dirichlet(1) chains, n=5 ==")
    cnt = sum(any(l not in (BOT, TOP) for l in fixed(dirichlet(n), F_A)) for _ in range(1000))
    print(f"  chains with a nontrivial Version-A fixed point: {cnt} / 1000")
    print("  codimension check: rank of KS constraints on the affine row-sum space vs (k-1)(n-k)")
    ok = True
    for l in partitions(n):
        k = max(l) + 1
        if k in (1, n):
            continue
        rows = []
        for b in range(k):
            mem = [i for i in range(n) if l[i] == b]
            for x in mem[1:]:
                for j in range(k - 1):
                    v = np.zeros((n, n)); Bj = [i for i in range(n) if l[i] == j]
                    v[x, Bj] += 1; v[mem[0], Bj] -= 1
                    rows.append(v.ravel())
        # restrict to directions with zero row sums: project out row-constant directions
        Cn = np.vstack(rows)
        Z = np.vstack([np.kron(np.eye(n)[i], np.ones(n)) for i in range(n)])     # row-sum functionals
        rk = np.linalg.matrix_rank(np.vstack([Cn, Z])) - np.linalg.matrix_rank(Z)
        ok &= (rk == (k - 1) * (n - k))
    print(f"  all nontrivial partitions of 5 states: codim == (k-1)(n-k): {ok}")

    print("== G1.4 symmetry ==")
    m = 6
    c = rng.random(m); Pc = np.array([np.roll(c, i) for i in range(m)]); Pc /= Pc.sum(1, keepdims=True)
    fA = fixed(Pc, F_A); fB = fixed(Pc, F_B)
    nontriv = [l for l in fB if 1 < max(l) + 1 < m]
    print(f"  circulant (Z6-equivariant) on 6 states: |Fix_A|={len(fA)} |Fix_B|={len(fB)}; nontrivial B: {nontriv}")
    z2, z3 = canon(tuple(i % 3 for i in range(m))), canon(tuple(i % 2 for i in range(m)))
    print(f"  Z2-orbit {z2} fixed: {z2 in fB};  Z3-orbit {z3} fixed: {z3 in fB};  comparable: "
          f"{refines(z2, z3) or refines(z3, z2)}")

    print("== G1.5 orientation ==")
    same = all(set(fixed(np.eye(n)[list(p)], F)) == set(fixed(np.eye(n)[list(p)].T, F))
               for p in itertools.permutations(range(n)) for F in (F_A, F_B))
    print(f"  all 120 bijections on 5 states: Fix(T) == Fix(T^-1) for both versions: {same}")
    W = rng.random((n, n)); W = W + W.T; Pr = W / W.sum(1, keepdims=True)
    pi = W.sum(1) / W.sum(); Pstar = (Pr.T * pi[None, :]).T
    Pstar = np.diag(1 / pi) @ Pr.T @ np.diag(pi)
    print(f"  reversible chain: P* == P: {np.allclose(Pstar, Pr)}")
    # nonreversible chain, lumpable for Pi = {0,1}{2,3,4}, with nonuniform pi inside blocks
    lab = (0, 0, 1, 1, 1)
    while True:
        P = np.zeros((n, n)); Q = rng.dirichlet(np.ones(2), size=2)
        for x in range(n):
            b = lab[x]
            for j, Bj in enumerate([[0, 1], [2, 3, 4]]):
                P[x, Bj] = Q[b, j] * rng.dirichlet(np.ones(len(Bj)))
        w, v = np.linalg.eig(P.T); pi = np.real(v[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
        if pi.min() > 1e-6:
            break
    Ps = np.diag(1 / pi) @ P.T @ np.diag(pi)
    print(f"  nonreversible chain (pi derived from P): Pi lumpable forward {lumpable(P, lab)}, "
          f"for the time reversal P* {lumpable(Ps, lab)}")
    fwd = set(fixed(P, F_A)); bwd = set(fixed(Ps, F_A))
    print(f"  Fix_A(P) = {sorted(fwd)};  Fix_A(P*) = {sorted(bwd)}")
    print(f"  covariance check: computing with P* as 'forward' simply swaps the two sets: "
          f"{set(fixed(np.diag(1/pi) @ Ps.T @ np.diag(pi), F_A)) == fwd}")

    print("== G1.6 twin control: 4 states = {0,1}^2, two candidate subsystem decompositions ==")
    states = [(a, b) for a in (0, 1) for b in (0, 1)]
    D1 = canon(tuple(a for a, b in states))            # decomposition 1, factor A = bit a
    D2 = canon(tuple(a ^ b for a, b in states))        # decomposition 2, factor A' = a xor b
    pa, pb = rng.dirichlet(np.ones(2), 2), rng.dirichlet(np.ones(2), 2)
    Pprod = np.kron(pa, pb)                            # non-interacting in decomposition 1
    Pint = rng.dirichlet(np.ones(4), 4)                # generic interacting dynamics
    for name, P in [("non-interacting in D1", Pprod), ("generic interacting", Pint)]:
        print(f"  {name:22s}: D1 factor fixed (A) {F_A(P, D1) == D1}; D2 factor fixed (A) {F_A(P, D2) == D2}")
    print("  (identical P for both decompositions; F distinguishes them only through exact autonomy of their factors)")
