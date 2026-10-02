"""B1 shared helpers (extracted from b1_net_reconstruction.py). GRUT BRIDGE-1 / B1 — can the canonical GRUT Level-0 local net be reconstructed with the net DELETED from the input?
Canonical models (READ-ONLY source GRUT-RAI master-w25bu9 @ b935099; re-implemented here, not imported):
  C1-a  : K = 0.3 I + L(path, unit springs), N = 24; retained site 0 + bath   (L0_1A_CHARTER_01.md)
  L0-1b : K = 0.3 I + L(w), w_ij = c_a |i-j|^-a                               (L0_1B_CHARTER_01.md)
  L0-1c : V = 1/2 x^T K_b x + beta sum_i x_i^4 (on-site convex quartic)        (L0_1C_CHARTER_01.md)
  L0-1e : additive noise Q = 2 diag(T_i)                                      (L0_ACCESS_BRIDGE_01.md)
Net deleted = the generator is given only as an abstract operator, i.e. in an unknown orthonormal frame W (K' = W K W^T).
Firewall: no reconstruction criterion may name the canonical sites, the path graph, or the site-resolved access."""
import numpy as np
rng = np.random.default_rng(20261002)

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def lap(w):
    L = -w.copy(); np.fill_diagonal(L, 0); np.fill_diagonal(L, -L.sum(1)); return L
def path_w(N):
    w = np.zeros((N, N))
    for i in range(N - 1): w[i, i + 1] = w[i + 1, i] = 1.0
    return w
def haar_orth(N):
    q, r = np.linalg.qr(rng.normal(size=(N, N))); return q * np.sign(np.diag(r))
def geometry_predicate(w):
    """L0-1b predicate on a spring network: resistance distances R_ij from L(w); Q = max/min over i != j; line ordering
    R_{1,j} strictly increasing in j (as frozen in L0_1B_CHARTER_01.md); SURVIVES iff Q >= N/2 and ordering monotone."""
    N = len(w); Lp = np.linalg.pinv(lap(w)); d = np.diag(Lp)
    R = d[:, None] + d[None, :] - 2 * Lp; off = R[~np.eye(N, dtype=bool)]
    Q = off.max() / off.min(); mono = all(R[1, j + 1] > R[1, j] for j in range(1, N - 1))
    return Q, mono, (Q >= N / 2 and mono)
def lanczos(A, q):
    N = len(A); Qm = np.zeros((N, N)); a = np.zeros(N); b = np.zeros(N - 1); q = q / np.linalg.norm(q); Qm[:, 0] = q
    for k in range(N):
        v = A @ Qm[:, k]; a[k] = Qm[:, k] @ v
        v -= Qm[:, :k + 1] @ (Qm[:, :k + 1].T @ v); v -= Qm[:, :k + 1] @ (Qm[:, :k + 1].T @ v)   # full reorthogonalization
        if k < N - 1:
            b[k] = np.linalg.norm(v)
            if b[k] < 1e-10: return None
            Qm[:, k + 1] = v / b[k]
    return Qm, np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
def as_spring_network(T):
    """read a symmetric matrix as pin*I + L(w): flip basis signs on a tree so off-diagonals are <= 0, then w = -offdiag,
    pins = diag - row sums of w. Returns (w, pins, cooperative)."""
    S = np.ones(len(T))
    for i in range(len(T) - 1):              # tridiagonal: sequential sign fixing along the chain
        if T[i, i + 1] * S[i] * S[i + 1] > 0: S[i + 1] *= -1
    Ts = S[:, None] * T * S[None, :]; w = -Ts.copy(); np.fill_diagonal(w, 0); w = np.where(w > 1e-12, w, 0.0)
    pins = np.diag(Ts) - w.sum(1); return w, pins, bool((Ts - np.diag(np.diag(Ts)) <= 1e-12).all())
def signed_perm_equiv(A, B, tol=1e-8):
    """are A and B related by a signed permutation (graph relabeling + sign gauge)? compare sorted spectra of |A| rows / multiset invariants"""
    inv = lambda M: (np.sort(np.round(np.abs(M[~np.eye(len(M), dtype=bool)]), 6)), np.sort(np.round(np.diag(M), 6)))
    a, b = inv(A), inv(B); return len(a[0]) == len(b[0]) and np.allclose(a[0], b[0], atol=1e-5) and np.allclose(a[1], b[1], atol=1e-5)

