"""SCOUT-2 S2-G3: spacetime / causal dimension.  Firewall: only causal order / propagation cones count."""
import numpy as np
from scipy.special import gamma
from scipy.optimize import brentq
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

rng = np.random.default_rng(5)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def f_mm(d): return 2 * gamma(d + 1) * gamma(d / 2) / (4 * gamma(3 * d / 2))   # expected ordering fraction, interval
def mm_dim(r): return brentq(lambda d: f_mm(d) - r, 1.01, 12)
def sprinkle_interval(N, d, T=1.0):
    pts = []
    while len(pts) < N:
        t = rng.uniform(0, T, 4 * N); x = rng.uniform(-T / 2, T / 2, (4 * N, d - 1))
        ok = np.linalg.norm(x, axis=1) < np.minimum(t, T - t)
        pts += list(np.c_[t[ok], x[ok]])
    return np.array(pts[:N])
def sprinkle_slab(N, d, T=1.0, W=4.0):
    return np.c_[rng.uniform(0, T, N), rng.uniform(-W / 2, W / 2, (N, d - 1))]
def relation(P):
    dt = P[None, :, 0] - P[:, None, 0]; dx = np.linalg.norm(P[None, :, 1:] - P[:, None, 1:], axis=2)
    return (dt > 0) & (dt ** 2 > dx ** 2)                      # R[i,j]: i precedes j
def ordering_fraction(R):
    N = R.shape[0]; return R.sum() / (N * (N - 1) / 2)
def midpoint_dim(R):
    N = R.shape[0]; past = R.sum(0); fut = R.sum(1); m = np.minimum(past, fut).max()
    return np.log2(N / m)
def longest_chain(R, P):
    order = np.argsort(P[:, 0]); Rs = R[np.ix_(order, order)]; L = np.ones(len(P), int)
    for j in range(len(P)):
        prev = np.where(Rs[:j, j])[0]
        if len(prev): L[j] = 1 + L[prev].max()
    return L.max()

hdr("G3-1/G3-2 sprinkled causal intervals in d-dim Minkowski: three order-theoretic estimators")
print("  d (true) |   N   | Myrheim-Meyer        | midpoint scaling     | longest chain L")
chains = {}
for d in (2, 3, 4):
    for N in (500, 2000):
        mmv, mid, lc = [], [], []
        for _ in range(4):
            P = sprinkle_interval(N, d); R = relation(P)
            mmv.append(mm_dim(ordering_fraction(R))); mid.append(midpoint_dim(R)); lc.append(longest_chain(R, P))
        chains[(d, N)] = np.mean(lc)
        print(f"     {d}     | {N:5d} | {np.mean(mmv):.3f} +- {np.std(mmv):.3f}     | {np.mean(mid):.3f} +- {np.std(mid):.3f}     | {np.mean(lc):.1f}")
    dlc = np.log(2000 / 500) / np.log(chains[(d, 2000)] / chains[(d, 500)])
    print(f"           chain-scaling estimate d = log(N2/N1)/log(L2/L1) = {dlc:.3f}")
print("\n  shape hostile: sprinkling into a SLAB (t in (0,1), |x_i| < 2), MM formula assumes an interval:")
for d in (2, 3, 4):
    P = sprinkle_slab(2000, d); R = relation(P); r = ordering_fraction(R)
    try: est = f"{mm_dim(r):.3f}"
    except ValueError: est = "out of range"
    print(f"     true d = {d}: ordering fraction {r:.4f} -> MM 'dimension' {est};  midpoint {midpoint_dim(R):.3f}")

hdr("G3-3 same spatial graph (Z^2 or Z), different causal structure: lattice spacetimes, interval between tips")
def lattice_interval(space_dim, T, norm, c=1, absolute_time=False):
    pts = []
    rng_ = range(-c * T, c * T + 1)
    for t in range(0, 2 * T + 1):
        rad = c * min(t, 2 * T - t)
        if space_dim == 1:
            for x in range(-rad, rad + 1): pts.append((t, x))
        else:
            for x in rng_:
                for y in rng_:
                    dist = abs(x) + abs(y) if norm == "L1" else max(abs(x), abs(y))
                    if dist <= rad: pts.append((t, x, y))
    P = np.array(pts, float); dt = P[None, :, 0] - P[:, None, 0]
    dX = np.abs(P[None, :, 1:] - P[:, None, 1:])
    dist = dX.sum(2) if norm == "L1" else dX.max(2)
    R = (dt > 0) & (dist <= c * dt) if not absolute_time else (dt > 0)
    return R
for sd, norm, T, c, ab in [(1, "L1", 30, 1, False), (2, "L1", 10, 1, False), (2, "L1", 16, 1, False), (2, "Linf", 10, 1, False),
                           (2, "Linf", 14, 1, False), (2, "L1", 6, 2, False), (2, "L1", 10, 1, True)]:
    R = lattice_interval(sd, T, norm, c, ab); r = ordering_fraction(R)
    try: est = f"{mm_dim(r):.3f}"
    except ValueError: est = "out of range"
    print(f"  spatial Z^{sd}, {norm:4s} {'cone' if not ab else 'ABSOLUTE TIME (instantaneous)'}, speed c={c}, T={T}: elements {R.shape[0]:5d}, ordering fraction {r:.4f} -> MM d = {est}; midpoint d = {midpoint_dim(R):.3f}")
print("  continuum reference (G3-1): 1+1 -> 2, 2+1 -> 3.")

hdr("G3-4 / G3-5 propagation cones from local hopping Hamiltonians (single-particle amplitude, tail threshold 1e-10)")
def lattice_H(kind, L):
    idx = lambda x, y: (x % L) * L + (y % L)
    rows, cols = [], []
    if kind == "chain":
        for x in range(L): rows += [x]; cols += [(x + 1) % L]
        n = L
    elif kind == "chain+NNN":
        for x in range(L): rows += [x, x]; cols += [(x + 1) % L, (x + 2) % L]
        n = L
    else:
        for x in range(L):
            for y in range(L):
                rows += [idx(x, y), idx(x, y)]; cols += [idx(x + 1, y), idx(x, y + 1)]
                if kind == "triangular": rows += [idx(x, y)]; cols += [idx(x + 1, y + 1)]
        n = L * L
    vals = np.ones(len(rows)); vals = vals if kind != "chain+NNN" else np.tile([1.0, 0.5], len(rows) // 2)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr(); return (A + A.T).tocsr(), n
def front(kind, L, dirs, ts):
    H, n = lattice_H(kind, L); c = (L // 2) * L + L // 2 if kind in ("square", "triangular") else L // 2
    psi0 = np.zeros(n, complex); psi0[c] = 1
    out = expm_multiply(-1j * H, psi0, start=0, stop=ts[-1], num=len(ts), endpoint=True)
    res = {}
    for dname, step in dirs.items():
        reach = []
        for k, t in enumerate(ts):
            amp = np.abs(out[k]) ** 2; r = 0
            for m in range(1, L // 2 - 1):
                if kind in ("square", "triangular"):
                    x, y = L // 2 + m * step[0], L // 2 + m * step[1]; site = (x % L) * L + (y % L); dist = m * np.hypot(*step)
                else:
                    site = (L // 2 + m) % L; dist = m
                if amp[site] > 1e-10: r = dist
            reach.append(r)
        v = np.polyfit(ts[2:], reach[2:], 1)[0]; res[dname] = v
    vol = [(np.abs(out[k]) ** 2 > 1e-10).sum() for k in range(len(ts))]
    dvol = np.polyfit(np.log(ts[3:]), np.log(np.maximum(vol[3:], 1)), 1)[0]
    return res, dvol
ts = np.linspace(2, 12, 6)
for kind, L, dirs in [("chain", 401, {"x": (1, 0)}), ("chain+NNN", 401, {"x": (1, 0)}),
                      ("square", 101, {"axis (1,0)": (1, 0), "diagonal (1,1)": (1, 1)}),
                      ("triangular", 101, {"axis (1,0)": (1, 0), "bond diag (1,1)": (1, 1), "anti-diag (1,-1)": (1, -1)})]:
    v, dv = front(kind, L, dirs, ts)
    print(f"  {kind:<11}: front speeds " + ", ".join(f"{k} {x:.2f}" for k, x in v.items()) + f";  cone-volume growth exponent {dv:.2f}"
          + (f";  anisotropy max/min {max(v.values())/min(v.values()):.2f}" if len(v) > 1 else ""))
print("  square lattice low-energy dispersion E(k) = -4 + |k|^2 + O(k^4): quadratic (non-relativistic); cone set by the")
print("  lattice maximum group velocity, direction-dependent -> preferred lattice frame: LORENTZ STRUCTURE NOT DERIVED.")
