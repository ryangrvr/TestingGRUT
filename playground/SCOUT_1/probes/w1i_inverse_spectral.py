"""SCOUT-1 W1-I: does access (retained-site spectral measure mu_r) + topology fix the one-particle contraction K?

 1. end-site readout on a chain: mu_r -> K unique up to the sign gauge b_i -> +-b_i (Jacobi/Stieltjes/Lanczos).
 2. interior-site readout on a chain: a continuous family of different chains with identical mu_r.
 3. non-chain graphs: different local graphs (different topology) with identical mu_r at the root;
    every (graph, root) has a chain representative (Lanczos) => "chain" is a gauge slice of the U-orbit.
 4. finite access: moments m_0..m_{2k-1} fix exactly the first k diagonal and k-1(+1) off-diagonal coefficients.
 5. stability: reconstruction of deep coefficients from noisy mu_r (ill-conditioning).
"""
import numpy as np
from mpmath import mp, mpf, matrix as mpmatrix

rng = np.random.default_rng(7)
np.set_printoptions(precision=6, suppress=True, linewidth=140)


def jacobi(a, b):
    return np.diag(a) + np.diag(b, 1) + np.diag(b, -1)


def measure(K, r, tol=1e-9):
    """retained-site spectral measure; degenerate eigenvalues merged (their total weight is basis-independent)"""
    lam, V = np.linalg.eigh(K)
    w = V[r, :] ** 2
    out_l, out_w = [lam[0]], [w[0]]
    for x, y in zip(lam[1:], w[1:]):
        if abs(x - out_l[-1]) < tol:
            out_w[-1] += y
        else:
            out_l.append(x); out_w.append(y)
    return np.array(out_l), np.array(out_w)


def stieltjes(lam, w, m):
    """Lanczos/Stieltjes on a discrete measure: recover Jacobi (a_0..a_{m-1}, b_0..b_{m-2})."""
    lam = np.asarray(lam, float); w = np.asarray(w, float)
    p_prev = np.zeros_like(lam); p = np.ones_like(lam) / np.sqrt(w.sum())
    a, b = [], []
    P = [p]
    for j in range(m):
        aj = np.sum(w * lam * p * p); a.append(aj)
        q = (lam - aj) * p - (b[-1] * p_prev if b else 0)
        for pp in P:  # reorthogonalize in L2(w)
            q -= np.sum(w * q * pp) * pp
        bj = np.sqrt(np.sum(w * q * q))
        if j < m - 1:
            b.append(bj)
            p_prev, p = p, q / bj
            P.append(p)
    return np.array(a), np.array(b)


print("=== 1. end-site readout on a chain: uniqueness up to sign gauge ===")
n = 30
a = 2.0 + 0.3 + rng.uniform(-0.5, 0.5, n)
b = -rng.uniform(0.5, 1.5, n - 1)          # negative couplings (as in K = pin + Laplacian)
K = jacobi(a, b)
lam, w = measure(K, 0)
a_r, b_r = stieltjes(lam, w, n)
print("  max|a - a_rec| =", np.max(np.abs(a - a_r)), "  max| |b| - b_rec | =", np.max(np.abs(np.abs(b) - b_r)))
S = np.diag(rng.choice([-1, 1], n)); S[0, 0] = 1
lam2, w2 = measure(S @ K @ S, 0)
print("  sign-gauge S K S (S_00 = 1): same mu_r:", np.allclose(lam, lam2) and np.allclose(w, w2),
      "  (S is a hidden-sector U preserving the chain frame)")

print("\n=== 2. interior-site readout: continuous family of distinct chains, same mu_r ===")
nL, nR = 6, 7
aL = 2.3 + rng.uniform(-0.3, 0.3, nL); bL = -rng.uniform(0.6, 1.2, nL - 1)
aR = 2.3 + rng.uniform(-0.3, 0.3, nR); bR = -rng.uniform(0.6, 1.2, nR - 1)
a0, cL, cR = 2.3, -1.0, -1.0


def glue(aL, bL, aR, bR, a0, cL, cR):
    # site order: L-chain reversed (its end site adjacent to centre), centre, R-chain
    nl, nr = len(aL), len(aR)
    N = nl + 1 + nr
    K = np.zeros((N, N))
    K[:nl, :nl] = jacobi(aL[::-1], bL[::-1])
    K[nl, nl] = a0
    K[nl + 1:, nl + 1:] = jacobi(aR, bR)
    K[nl - 1, nl] = K[nl, nl - 1] = cL
    K[nl, nl + 1] = K[nl + 1, nl] = cR
    return K, nl


K1, c1 = glue(aL, bL, aR, bR, a0, cL, cR)
# combined "environment" measure sigma = cL^2 mu_L + cR^2 mu_R  (mu_X = end-site measure of each half)
lL, wL = measure(jacobi(aL, bL), 0)
lR, wR = measure(jacobi(aR, bR), 0)
atoms = np.concatenate([lL, lR]); mass = np.concatenate([cL ** 2 * wL, cR ** 2 * wR])
# re-split: move atoms between halves (keep counts nL, nR) and rescale weights
perm = rng.permutation(nL + nR)
iL, iR = np.sort(perm[:nL]), np.sort(perm[nL:])
mL, mR = mass[iL].sum(), mass[iR].sum()
aL2, bL2 = stieltjes(atoms[iL], mass[iL] / mL, nL)
aR2, bR2 = stieltjes(atoms[iR], mass[iR] / mR, nR)
K2, c2 = glue(aL2, -bL2, aR2, -bR2, a0, -np.sqrt(mL), -np.sqrt(mR))
l1, w1 = measure(K1, c1); l2, w2 = measure(K2, c2)
print("  chain A couplings (left|right):", np.round(np.abs(bL), 3), "|", np.round(np.abs(bR), 3))
print("  chain B couplings (left|right):", np.round(bL2, 3), "|", np.round(bR2, 3))
print("  centre couplings A:", (abs(cL), abs(cR)), "  B:", (round(np.sqrt(mL), 4), round(np.sqrt(mR), 4)))
print("  same mu_r at the centre:", np.allclose(l1, l2, atol=1e-10) and np.allclose(w1, w2, atol=1e-10),
      "   spectra of K equal:", np.allclose(np.linalg.eigvalsh(K1), np.linalg.eigvalsh(K2)))
print("  isospectral but not isometric: |K_A - K_B|_F =", np.linalg.norm(K1 - K2).round(4),
      "; mirror image of B equals A?", np.allclose(K1, K2[::-1, ::-1]))
from math import comb
print("  atom partitions of sigma (13 atoms) into a left half of 6 and right half of 7:", comb(13, 6),
      "-> that many distinct chains (finite case); all splits n_L'+n_R'=13:", 2 ** 13 - 2, "(ordered)")
print("  infinite-chain case: sigma = f*sigma + (1-f)*sigma for any measurable 0<=f<=1 -> continuous family")

print("\n=== 3. different topology, same mu_r at the root ===")
# graph G: root 0 - 1, site 1 branches to leaves 2,3 (all couplings -1, diag 2.3)
G = np.diag([2.3, 2.3, 2.3, 2.3]);
for i, j in ((0, 1), (1, 2), (1, 3)):
    G[i, j] = G[j, i] = -1.0
lg, wg = measure(G, 0)
mask = wg > 1e-10
ag, bg = stieltjes(lg[mask], wg[mask], mask.sum())
print("  star-tree root measure has", mask.sum(), "atoms; Lanczos chain: a =", ag.round(6), " b =", bg.round(6))
C = jacobi(ag, -bg)
lc, wc = measure(C, 0)
print("  chain with b1 = sqrt(2) reproduces mu_root:", np.allclose(lc, lg[mask]) and np.allclose(wc, wg[mask]))
# uniform-coupling example: 2D square-lattice corner/boundary site vs its Lanczos chain (non-uniform b)
L = 12
Kh = np.zeros((L * L, L * L))
for i in range(L):
    for j in range(L):
        s = i * L + j; Kh[s, s] = 4.3
        for di, dj in ((1, 0), (0, 1)):
            if i + di < L and j + dj < L:
                t = (i + di) * L + j + dj; Kh[s, t] = Kh[t, s] = -1.0
lh, wh = measure(Kh, 0)
keep = wh > 1e-14
ah, bh = stieltjes(lh[keep], wh[keep], 8)
print("  12x12 square-lattice corner: Lanczos chain b_0..b_6 =", bh.round(4),
      " (a 1D chain with these couplings has identical corner mu_r)")

print("\n=== 4. finite access depth: 2k moments <-> k Jacobi layers ===")
mp.dps = 50
lam_mp = [mpf(x) for x in lam]; w_mp = [mpf(x) for x in w]


def moments_mp(lam, w, M):
    return [sum(wi * li ** k for li, wi in zip(lam, w)) for k in range(M)]


def jacobi_from_moments(mom, k):
    """Chebyshev/Gram via Hankel Cholesky (mpmath): returns a_0..a_{k-1}, b_0..b_{k-2}"""
    Hk = mpmatrix(k + 1, k + 1)
    for i in range(k + 1):
        for j in range(k + 1):
            Hk[i, j] = mom[i + j]
    # Cholesky H = R^T R ; a_j = R[j,j+1]/R[j,j] - R[j-1,j]/R[j-1,j-1] ; b_j = R[j+1,j+1]/R[j,j]
    R = mp.cholesky(Hk).T
    a_, b_ = [], []
    for j in range(k):
        aj = R[j, j + 1] / R[j, j] - (R[j - 1, j] / R[j - 1, j - 1] if j > 0 else 0)
        a_.append(aj)
        if j < k - 1:
            b_.append(R[j + 1, j + 1] / R[j, j])
    return np.array([float(x) for x in a_]), np.array([float(x) for x in b_])


for k in (3, 6, 10):
    mom = moments_mp(lam_mp, w_mp, 2 * k + 1)
    ak, bk = jacobi_from_moments(mom, k)
    print(f"  k={k:2d}: from m_0..m_{2*k}: max|a-a_true| over first {k} = {np.max(np.abs(ak - a[:k])):.2e}, "
          f"max|b-|b_true|| over first {k-1} = {np.max(np.abs(bk - np.abs(b[:k-1]))):.2e}")
# two chains sharing the first k layers have identical low moments but different mu_r
k = 5
a_alt = a.copy(); b_alt = b.copy()
a_alt[k:] += rng.uniform(-0.4, 0.4, n - k); b_alt[k:] *= rng.uniform(0.7, 1.3, n - k - 1)
lam_alt, w_alt = measure(jacobi(a_alt, b_alt), 0)
mA = [np.sum(w * lam ** j) for j in range(14)]; mB = [np.sum(w_alt * lam_alt ** j) for j in range(14)]
print("  chains equal on layers < 5, different deeper: relative moment differences m_0..m_13:")
print("   ", " ".join(f"{(x - y) / x:.1e}" for x, y in zip(mA, mB)))

print("\n=== 5. stability: noisy mu_r -> deep coefficients ===")
for eps in (1e-12, 1e-8, 1e-4):
    errs = []
    for rep in range(20):
        wn = w * (1 + eps * rng.standard_normal(n)); wn /= wn.sum()
        ln = lam + eps * rng.standard_normal(n) * (lam.max() - lam.min())
        ar, br = stieltjes(np.sort(ln), wn[np.argsort(ln)], n)
        errs.append(np.abs(br - np.abs(b)))
    e = np.median(errs, axis=0)
    print(f"  eps={eps:.0e}: median |b_j error| at j = 0, 5, 10, 20, 28:",
          " ".join(f"{e[j]:.1e}" for j in (0, 5, 10, 20, 28)))

print("\n=== 5b. moment route conditioning (Hankel matrix of mu_r, mpmath) ===")
for k in (4, 8, 12, 16):
    mom = moments_mp(lam_mp, w_mp, 2 * k + 1)
    Hk = mpmatrix(k + 1, k + 1)
    for i in range(k + 1):
        for j in range(k + 1):
            Hk[i, j] = mom[i + j]
    ev = mp.eigsy(Hk)[0]
    print(f"  k={k:2d}: cond(H_k) = {float(max(ev) / min(ev)):.2e}")
