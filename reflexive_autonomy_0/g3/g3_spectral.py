"""RA0 / G3 -- canonical asymptotic spectral hierarchy.  NUMERICAL ILLUSTRATION; independent code path, not reviewer.

Allowed data only (cross-N firewall G3.4): ordered eigenvalues, ranks, adjacent ratios, spectral counting functions,
spectral projectors.  No eigenvector is matched across sizes.  Degenerate eigenvalues are merged ("distinct levels"),
so every cut sits between distinct levels and its projector is basis-free.

For each family we report, per size, the largest adjacent ratio  r = lambda_(next distinct level) / lambda_(level)
over the NONZERO part of the spectrum (the zero eigenspace is handled separately), the rank of the cut (dimension of
the projector below it, zero modes included), and the trend in size.  A cut is "diverging" iff its ratio grows without
bound along a rank-defined sequence.
"""
import itertools
import numpy as np
import scipy.sparse as sp

np.set_printoptions(precision=4, suppress=True)


def levels(ev, tol=1e-9):
    ev = np.sort(np.real(ev)); lv, mult = [], []
    scale = max(abs(ev).max(), 1e-300)
    for x in ev:
        if lv and abs(x - lv[-1]) <= tol * abs(x) + 1e-14 * scale:     # relative merge tolerance
            mult[-1] += 1
        else:
            lv.append(x); mult.append(1)
    return np.array(lv), np.array(mult)


def cuts(ev, top=3):
    lv, mult = levels(ev)
    nz = lv > 1e-13 * max(lv.max(), 1e-300)          # zero modes: relative threshold
    z = int(mult[~nz].sum())                    # zero-eigenspace dimension
    lvn, mn = lv[nz], mult[nz]
    r = lvn[1:] / lvn[:-1]
    rank = z + np.cumsum(mn)[:-1]               # projector dimension below each cut
    order = np.argsort(r)[::-1][:top]
    return z, [(float(r[i]), int(rank[i])) for i in order], lvn


# ------------------------------------------------------------------ A independent particles
def family_A():
    print("== A  independent two-state particles (a=0.1, b=0.25) ==")
    a, b = 0.1, 0.25; p = np.array([[1 - a, a], [b, 1 - b]])
    for N in [2, 4, 6, 8, 10]:
        P = p
        for _ in range(N - 1):
            P = np.kron(P, p)
        z, top, _ = cuts(1 - np.linalg.eigvals(P))
        print(f"  N={N:2d}: zero-dim {z}; largest adjacent ratios (ratio, rank): {[(round(x, 3), k) for x, k in top]}")


# ------------------------------------------------------------------ B SSEP ring
def ssep(L, N):
    confs = list(itertools.combinations(range(L), N)); idx = {c: i for i, c in enumerate(confs)}
    rows, cols, vals = [], [], []
    for i, c in enumerate(confs):
        occ = set(c)
        for bnd in range(L):
            u, v = bnd, (bnd + 1) % L
            if (u in occ) != (v in occ):
                new = set(occ); new.symmetric_difference_update({u, v})
                rows.append(i); cols.append(idx[tuple(sorted(new))]); vals.append(1.0 / L)
            else:
                rows.append(i); cols.append(i); vals.append(1.0 / L)
    return confs, sp.csr_matrix((vals, (rows, cols)), shape=(len(confs), len(confs))).toarray()


def family_B():
    print("== B  SSEP ring, N = L/2 (full many-particle kernel) ==")
    for L in [8, 10, 12, 14]:
        confs, P = ssep(L, L // 2)
        w, V = np.linalg.eigh(P); ev = 1 - w
        z, top, lvn = cuts(ev)
        cnt = [int(np.sum(ev[ev > 1e-13] <= x * lvn[0] * (1 + 1e-9))) for x in (1, 4, 9, 16)]
        print(f"  L={L:2d} ({len(confs)} states): zero-dim {z}; top adjacent ratios {[(round(x, 3), k) for x, k in top]};"
              f"  spread lambda_max/lambda_1 = {lvn[-1] / lvn[0]:.1f};  counting N(x*lambda_1), x=1,4,9,16: {cnt}")
    print("  one-particle sector, large L (exact formula): largest adjacent distinct-level ratio and spread")
    for L in [100, 1000, 10000]:
        lam = (2.0 / L) * 2 * np.sin(np.pi * np.arange(0, L // 2 + 1) / L) ** 2      # 1-cos = 2 sin^2, no cancellation
        z, top, lvn = cuts(lam)
        cnt = [int(np.sum(lam[1:L // 2] <= x * lvn[0] * (1 + 1e-9))) * 2 for x in (1, 4, 9, 16)]
        print(f"    L={L:5d}: top ratio {top[0][0]:.3f} at rank {top[0][1]};  spread {lvn[-1] / lvn[0]:.3g};"
              f"  counting (x=1,4,9,16, both +/-q) {cnt}")


def algebra_test_ssep(L=10):
    print(f"== G3.6 algebra test, SSEP L={L}: lowest nonzero eigenspace E1 (+ constants) ==")
    confs, P = ssep(L, L // 2); w, V = np.linalg.eigh(P); ev = 1 - w
    lv, mult = levels(ev)
    lam1 = lv[lv > 1e-10][0]
    E = V[:, np.abs(ev - lam1) < 1e-9]; C = np.ones((len(confs), 1)) / np.sqrt(len(confs))
    B = np.column_stack([C, E]); Q, _ = np.linalg.qr(B)
    worst = 0.0
    for i in range(B.shape[1]):
        for j in range(i, B.shape[1]):
            u = B[:, i] * B[:, j]; worst = max(worst, np.linalg.norm(u - Q @ (Q.T @ u)) / np.linalg.norm(u))
    print(f"  dim(E1 + constants) = {B.shape[1]}; P-invariant by construction; product-closure residual = {worst:.3f}")
    # algebra generated: functions of the E1 coordinates -> level-set partition of configuration space
    coords = np.round(E, 8)
    keys = {}
    lab = np.array([keys.setdefault(tuple(r), len(keys)) for r in coords])
    k = lab.max() + 1
    S = np.column_stack([P[:, lab == j].sum(1) for j in range(k)])
    dfc = max(np.abs(S[lab == i] - S[lab == i][0]).max() for i in range(k))
    print(f"  generated algebra = functions of E1 coordinates: {k} level sets of {len(confs)} states (proper: {k < len(confs)});"
          f"  its one-step lumpability defect = {dfc:.3f} (semigroup-closed iff 0)")


# ------------------------------------------------------------------ C two-lane ladder (one particle)
def ladder(L, g):
    n = 2 * L; Gm = np.zeros((n, n))
    for lane in range(2):
        for x in range(L):
            i = lane * L + x
            for d in (1, -1):
                Gm[i, lane * L + (x + d) % L] += 1.0
            Gm[i, (1 - lane) * L + x] += g
    np.fill_diagonal(Gm, -Gm.sum(1))
    return -np.linalg.eigvalsh(Gm)              # positive relaxation rates


def family_C():
    print("== C  two-lane ladder (one particle), inter-lane rate gamma_L = L^-beta ==")
    for beta in [0, 1, 2, 3]:
        row = []
        for L in [50, 100, 200, 400]:
            z, top, lvn = cuts(ladder(L, L ** (-float(beta))))
            row.append(f"L={L}: {top[0][0]:.3g}@rank{top[0][1]}")
        print(f"  beta={beta}: largest adjacent ratio @ rank:  " + "  ".join(row))


# ------------------------------------------------------------------ D metastable positive control
def family_D():
    print("== D  metastable positive control ==")
    # D1: K=3 basins, complete graphs of size M, intra rate 1 (to each state, normalised by M), inter rate M^-1
    for M in [5, 10, 20, 40, 80]:
        K = 3; n = K * M; Gm = np.zeros((n, n))
        for i in range(n):
            bi = i // M
            for j in range(n):
                if i == j:
                    continue
                Gm[i, j] = (1.0 / M) if j // M == bi else (1.0 / M) * (1.0 / M) / ((K - 1))
        np.fill_diagonal(Gm, -Gm.sum(1))
        w, V = np.linalg.eigh(-Gm)
        z, top, _ = cuts(w)
        # slow projector: eigenvalues below the top cut
        rank = top[0][1]; E = V[:, np.argsort(w)[:rank]]
        Q, _ = np.linalg.qr(E); worst = 0.0
        for i in range(rank):
            for j in range(i, rank):
                u = E[:, i] * E[:, j]; worst = max(worst, np.linalg.norm(u - Q @ (Q.T @ u)) / np.linalg.norm(u))
        ind = np.column_stack([(np.arange(n) // M == b).astype(float) for b in range(K)])
        resid = np.linalg.norm(ind - Q @ (Q.T @ ind)) / np.linalg.norm(ind)
        nxt = f"{top[1][0]:.2f}" if len(top) > 1 else "none (only two distinct nonzero levels)"
        print(f"  D1 3 basins, M={M:3d}: top ratio {top[0][0]:8.2f} at rank {top[0][1]} (next {nxt});"
              f"  product-closure residual {worst:.2e};  basin-indicator residual (diagnostic) {resid:.2e}")
    # D2: nested: 4 basins in 2 super-basins; intra 1, within-super inter M^-1, between-super M^-2
    for M in [5, 10, 20, 40]:
        K = 4; n = K * M; Gm = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                bi, bj = i // M, j // M
                if bi == bj:
                    Gm[i, j] = 1.0 / M
                elif bi // 2 == bj // 2:
                    Gm[i, j] = (1.0 / M) / M
                else:
                    Gm[i, j] = (1.0 / M) / M ** 2
        np.fill_diagonal(Gm, -Gm.sum(1))
        z, top, _ = cuts(np.linalg.eigvalsh(-Gm))
        print(f"  D2 nested 4 basins / 2 super, M={M:3d}: top cuts (ratio, rank) {[(round(x, 1), k) for x, k in top[:2]]}")


def family_D3():
    print("== D3 non-symmetric metastable control: 3 basins, random symmetric intra-basin rates, inter rates random * M^-1 ==")
    rng = np.random.default_rng(7)
    for M in [5, 10, 20, 40, 80]:
        K = 3; n = K * M
        W = rng.random((n, n)); W = (W + W.T) / 2
        same = (np.arange(n)[:, None] // M) == (np.arange(n)[None, :] // M)
        Gm = np.where(same, W / M, W / M ** 2); np.fill_diagonal(Gm, 0.0)
        np.fill_diagonal(Gm, -Gm.sum(1))
        w, V = np.linalg.eigh(-Gm)
        z, top, _ = cuts(w)
        rank = top[0][1]; E = V[:, np.argsort(w)[:rank]]; Q, _ = np.linalg.qr(E); worst = 0.0
        for i in range(rank):
            for j in range(i, rank):
                u = E[:, i] * E[:, j]; worst = max(worst, np.linalg.norm(u - Q @ (Q.T @ u)) / np.linalg.norm(u))
        print(f"  M={M:3d}: top ratio {top[0][0]:8.2f} at rank {top[0][1]};  product-closure residual {worst:.3e}")


if __name__ == "__main__":
    family_A(); family_B(); algebra_test_ssep(10); family_C(); family_D(); family_D3()
