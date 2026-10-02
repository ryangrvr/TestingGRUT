"""B1 supplementary: (i) GLOBAL search in the L0-1b class (K = 0.3 I + weighted Laplacian) for an inequivalent isospectral net;
(ii) can the 2D-grid generator be read as a PASSIVE chain (pins >= 0) in some frame?  Same canonical models as b1_net_reconstruction.py."""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from b1_lib import lap, path_w, lanczos, as_spring_network, geometry_predicate, signed_perm_equiv
rng = np.random.default_rng(7)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def laplacian_frames(L0, restarts, scale):
    """minimize the positive part of off-diagonals of E L0 E^T over rotations E fixing the uniform vector (all generators)"""
    n = len(L0); one = np.ones(n) / np.sqrt(n)
    B = np.linalg.qr(np.c_[one, rng.normal(size=(n, n - 1))])[0][:, 1:]           # basis of the complement of the uniform vector
    m = n - 1; iu = np.triu_indices(m, 1)
    def build(a):
        A = np.zeros((m, m)); A[iu] = a; A = A - A.T; E = B @ expm(A) @ B.T + np.outer(one, one); return E @ L0 @ E.T
    def viol(a):
        off = build(a)[~np.eye(n, dtype=bool)]; return np.sum(np.clip(off, 0, None) ** 2)
    out = []
    for r_ in range(restarts):
        r = minimize(viol, rng.normal(size=len(iu[0])) * scale, method="L-BFGS-B", options={"maxiter": 3000})
        Lw = build(r.x); out.append((viol(r.x), Lw))
    return out

hdr("(i) GLOBAL search, L0-1b class: other frames where K - 0.3 I is a weighted graph Laplacian (canonical: path, N = 12 and 24)")
for N in (12, 24):
    L0 = lap(path_w(N)); K = 0.3 * np.eye(N) + L0
    res = laplacian_frames(L0, restarts=12 if N == 12 else 6, scale=1.0)
    feas = [(v, Lw) for v, Lw in res if v < 1e-14]
    ineq = [Lw for v, Lw in feas if not signed_perm_equiv(Lw + 0.3 * np.eye(N), K)]
    print(f"  N={N}: restarts {len(res)}; exactly feasible Laplacian frames {len(feas)}; of these INEQUIVALENT to the path {len(ineq)};"
          f" best residual violation {min(v for v, _ in res):.2e}")
    for Lw in ineq[:2]:
        w_ = -Lw.copy(); np.fill_diagonal(w_, 0); w_ = np.where(w_ > 1e-10, w_, 0)
        print(f"     inequivalent net: {np.count_nonzero(w_)//2} edges; geometry predicate (Q, monotone, SURVIVES) = {geometry_predicate(w_)}")

hdr("(ii) the 2D grid generator: is there ANY passive chain frame (all pins >= 0)?  maximize the minimum pin over Lanczos start vectors")
Lg = 4; n = Lg * Lg; w2 = np.zeros((n, n))
for x in range(Lg):
    for y in range(Lg):
        i = x * Lg + y
        if x + 1 < Lg: w2[i, i + Lg] = w2[i + Lg, i] = rng.uniform(0.8, 1.2)
        if y + 1 < Lg: w2[i, i + 1] = w2[i + 1, i] = rng.uniform(0.8, 1.2)
K2 = 0.3 * np.eye(n) + lap(w2)
def minpin(q):
    out = lanczos(K2, q)
    if out is None: return -1e3
    return as_spring_network(out[1])[1].min()
best = -1e9; bq = None
for r_ in range(10):
    r = minimize(lambda q: -minpin(q), rng.normal(size=n), method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-9})
    if -r.fun > best: best = -r.fun; bq = r.x
print(f"  weighted 4x4 grid (pin 0.3): best achievable minimum pin over Lanczos chain frames = {best:.4f} (passive chain exists iff >= 0)")
for s in range(n):
    e = np.zeros(n); e[s] = 1
    print(f"  start at canonical site {s:2d}: min pin of its chain frame {minpin(e):+.4f}") if s in (0, 5) else None

# --- follow-ups
hdr("(iii) inspect the near-feasible N=12 frame and the passive chain frame of the 2D grid")
N = 12; L0 = lap(path_w(N)); K = 0.3 * np.eye(N) + L0
res = sorted(laplacian_frames(L0, restarts=12, scale=1.0), key=lambda t: t[0])
v, Lw = res[0]; Lc = np.where(np.abs(Lw) < 1e-4, 0.0, Lw)
print(f"  N=12 best frame: violation {v:.1e}; after zeroing |entries| < 1e-4: signed-permutation-equivalent to the path: {signed_perm_equiv(Lc + 0.3*np.eye(N), K, tol=1e-4)};"
      f" off-diagonal nonzeros {np.count_nonzero(Lc[~np.eye(N, dtype=bool)])} (path: {2*(N-1)})")
out = lanczos(K2, bq); wch, pins, _ = as_spring_network(out[1])
print(f"  4x4 grid passive chain frame: min pin {pins.min():+.4f}; chain geometry predicate (Q, monotone, SURVIVES) = {geometry_predicate(wch)};"
      f" canonical 2D grid net predicate = {geometry_predicate(w2)}")
hdr("(iv) verify the N=12 isospectral weighted Laplacian is a genuine L0-1b-class member")
w12 = -Lw.copy(); np.fill_diagonal(w12, 0)
print(f"  min off-diagonal weight {w12[~np.eye(N, dtype=bool)].min():+.2e} (>= 0 required); max |row sum of Lw| {np.abs(Lw.sum(1)).max():.1e};"
      f" max |spectrum difference| vs path Laplacian {np.abs(np.linalg.eigvalsh(Lw) - np.linalg.eigvalsh(L0)).max():.1e}")
w12c = np.where(w12 > 1e-12, w12, 0.0)
print(f"  positive-weight edges {np.count_nonzero(w12c)//2} of {N*(N-1)//2}; weight range [{w12c[w12c>0].min():.3e}, {w12c.max():.3f}];"
      f" geometry predicate (Q, monotone, SURVIVES) = {geometry_predicate(w12c)}")
