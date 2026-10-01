"""SCOUT-2 S2-1: system individuation from (Hilbert space, H) with NO declared tensor-product structure.

(a) No-go baseline: for a fixed total dimension, any two TPSs are related by a unitary V; (H, T) ~ (V^dag H V, T0).
    With no criterion, H fixes only its spectrum, and every TPS is equally good.
(b) Locality selector: criterion = "H is 2-local for some qubit TPS". Map from 2-local couplings c (P of them) to the
    spectrum. Local-unitary (LU) orbits give 3n spectrum-preserving directions. Locally unique 2-local TPS (mod LU)
    iff rank J = P - 3n. Jacobian by Hellmann-Feynman: dlam_j/dc_k = <v_j|O_k|v_j>.
(c) Explicit hostile at n = 3, 4: two LU-inequivalent 2-local Hamiltonians with identical spectra (=> two different
    TPSs in which the same abstract H is 2-local). LU invariants: singular values of each pair's 3x3 coupling matrix,
    norms of the local fields (compared as multisets, to allow qubit permutations).
"""
import itertools
import numpy as np
import scipy.sparse as sps
from scipy.optimize import least_squares

rng = np.random.default_rng(2)
P1 = [sps.csr_matrix(np.array(m, dtype=complex)) for m in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
I2 = sps.identity(2, format="csr", dtype=complex)


def op(n, factors):
    out = sps.identity(1, format="csr", dtype=complex)
    for q in range(n):
        out = sps.kron(out, factors.get(q, I2), format="csr")
    return out


def basis_ops(n):
    ops, labels = [], []
    for i in range(n):
        for a in range(3):
            ops.append(op(n, {i: P1[a]})); labels.append(("1", i, a))
    for i, j in itertools.combinations(range(n), 2):
        for a in range(3):
            for b in range(3):
                ops.append(op(n, {i: P1[a], j: P1[b]})); labels.append(("2", (i, j), (a, b)))
    return ops, labels


def ham(ops, c):
    H = sum(ck * O for ck, O in zip(c, ops))
    return H.toarray()


def jac_rank(n, tol=1e-8):
    ops, _ = basis_ops(n)
    P = len(ops)
    c = rng.standard_normal(P)
    w, V = np.linalg.eigh(ham(ops, c))
    J = np.empty((2 ** n, P))
    for k, O in enumerate(ops):
        J[:, k] = np.real(np.einsum("ij,ij->j", V.conj(), O @ V))
    s = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(s > tol * s[0]))
    gap = np.min(np.diff(w))
    return P, rank, gap


def invariants(c, n):
    k = 0
    fields, pairs = [], []
    for i in range(n):
        fields.append(np.linalg.norm(c[k:k + 3])); k += 3
    for _ in itertools.combinations(range(n), 2):
        pairs.append(tuple(np.round(np.linalg.svd(c[k:k + 9].reshape(3, 3), compute_uv=False), 6))); k += 9
    return sorted(np.round(fields, 6)), sorted(pairs)


if __name__ == "__main__":
    print("=== (b) Jacobian rank of 2-local couplings -> spectrum, modulo local unitaries ===")
    print("   n  P(2-local)  P-3n  2^n-1  rank J  fibre dim = P-3n-rank   verdict        min level gap")
    for n in range(3, 11):
        P, rank, gap = jac_rank(n)
        fib = P - 3 * n - rank
        verdict = "LOCALLY UNIQUE" if fib == 0 else "NON-UNIQUE"
        print(f"  {n:2d}  {P:9d}  {P - 3 * n:5d}  {2 ** n - 1:5d}  {rank:6d}  {fib:20d}   {verdict:14s}  {gap:.2e}", flush=True)
    print("  image dimension of the 2-local family in spectrum space = rank; spectrum space dim = 2^n - 1")
    print("  -> once P-3n < 2^n-1 (here n >= 8) the image is a proper submanifold: a Lebesgue-null set of spectra admits ANY 2-local TPS")

    print("\n=== (c) explicit same-spectrum, LU-inequivalent 2-local Hamiltonians ===")
    for n in (3, 4):
        ops, _ = basis_ops(n)
        P = len(ops)
        c = rng.standard_normal(P)
        target = np.linalg.eigvalsh(ham(ops, c))
        c0 = rng.standard_normal(P)              # independent random start
        res = least_squares(lambda x: np.linalg.eigvalsh(ham(ops, x)) - target, c0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        c2 = res.x
        spec_err = np.max(np.abs(np.linalg.eigvalsh(ham(ops, c2)) - target))
        f1, p1 = invariants(c, n); f2, p2 = invariants(c2, n)
        print(f"  n = {n}: max |spectrum difference| = {spec_err:.2e}")
        print(f"     H : local-field norms {f1}")
        print(f"     H': local-field norms {f2}")
        print(f"     pair coupling singular values differ: {p1 != p2};  H: {p1[0]}  H': {p2[0]} (first pair, sorted)")
        print(f"     => LU-inequivalent (and not a qubit permutation): {f1 != f2 or p1 != p2}")
