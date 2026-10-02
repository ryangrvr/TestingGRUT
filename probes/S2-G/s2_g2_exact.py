"""SCOUT-2 S2-G2 (exact / large-size companion to s2_g2_graph.py).
Generator = combinatorial Laplacian (rate 1 per edge).  Lattice return probabilities from exact product / circulant
spectra; local estimators (fixed site) for comb and Sierpinski via sparse heat kernels."""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import numpy as np
import networkx as nx
from scipy.sparse.linalg import expm_multiply
import scipy.sparse as sp

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def ring_spec(N, w=None):
    """eigenvalues of a translation-invariant ring generator with rates w[r] for |jump| = r"""
    k = np.arange(N)
    if w is None: w = {1: 1.0}
    lam = np.zeros(N)
    for r, wr in w.items(): lam += 2 * wr * (1 - np.cos(2 * np.pi * k * r / N))
    return lam
def P1(lam, t): return np.exp(-np.outer(t, lam)).mean(1)
def ds_from_P(Pfun, ts):
    return [-2 * (np.log(Pfun(np.array([t * 1.05]))[0]) - np.log(Pfun(np.array([t / 1.05]))[0])) / (2 * np.log(1.05)) for t in ts]
ts = (1, 3, 10, 30, 100, 300, 1000)

hdr("G2-1 exact lattice controls (P_d(t) = P_1(t)^d; L = 4000)")
lam = ring_spec(4000)
for d in (1, 2, 3, 4):
    print(f"  Z^{d} torus: d_s(t={ts}) = " + " ".join(f"{x:5.3f}" for x in ds_from_P(lambda t: P1(lam, t) ** d, ts)))
G3 = nx.grid_graph(dim=[41, 41, 41], periodic=True); c = (20, 20, 20)
dist = nx.single_source_shortest_path_length(G3, c, cutoff=16); Nr = np.cumsum(np.bincount(list(dist.values())))
print("  3D torus 41^3 growth d_g(r=4,8,12,15): " + " ".join(f"{np.log(Nr[r+1]/Nr[r-1])/np.log((r+1)/(r-1)):.3f}" for r in (4, 8, 12, 15)))
idx = {v: i for i, v in enumerate(G3.nodes())}; A = nx.to_scipy_sparse_array(G3); Lc = sp.diags(np.asarray(A.sum(1)).ravel()) - A
r2 = np.zeros(A.shape[0]); dfull = nx.single_source_shortest_path_length(G3, c)
for v, dv in dfull.items(): r2[idx[v]] = dv ** 2
p0 = np.zeros(A.shape[0]); p0[idx[c]] = 1
tt = np.array([4.0, 8.0, 16.0, 32.0]); Pt = expm_multiply(-Lc.tocsr(), p0, start=0, stop=32, num=65, endpoint=True)
grid = np.linspace(0, 32, 65); m2 = [Pt[np.argmin(abs(grid - t))] @ r2 for t in tt]
print("  3D torus 41^3 walk d_w (t=4-8,8-16,16-32): " + " ".join(f"{2/(np.log(m2[i+1]/m2[i])/np.log(2)):.3f}" for i in range(3)))

hdr("G2-3 exact crossovers (product spectra; eps = 0.01)")
lam_big = ring_spec(3000)
print("  2D bundled chains (strong x, weak y): d_s(t) = " + " ".join(f"{x:5.2f}" for x in ds_from_P(lambda t: P1(lam_big, t) * P1(lam_big, 0.01 * t), ts + (3000, 10000))))
lam_pl = ring_spec(1500)
print("  3D layered (2D planes, weak z):       d_s(t) = " + " ".join(f"{x:5.2f}" for x in ds_from_P(lambda t: P1(lam_pl, t) ** 2 * P1(lam_pl, 0.01 * t), ts + (3000, 10000))))
print(f"  (t = {ts + (3000, 10000)})")

hdr("G2-4 diffusion-operator hostile, exact (same vertex set)")
lam2 = ring_spec(3000)
for a in (2.0, 1.0, 0.5):
    print(f"  Z^2, fractional generator Delta^(a/2), a={a}: d_s = " + " ".join(f"{x:5.2f}" for x in ds_from_P(lambda t: (np.exp(-np.outer(t, np.add.outer(lam2[::3], lam2[::3]).ravel() ** (a / 2))).mean(1) - 1 / 1000 ** 2), (3, 10, 30, 100))) + f"   (theory {4/a:.0f})")
N = 40000; r = np.arange(1, N // 2)
for al in (None, 1.5, 1.0, 0.5):
    w = {1: 1.0} if al is None else None
    if al is None:
        lamr = ring_spec(N)
    else:
        k = np.arange(N); wr = r ** (-(1 + al))
        lamr = 2 * (wr[None, :] * (1 - np.cos(2 * np.pi * np.outer(k[: N // 2 + 1], r) / N))).sum(1)
        lamr = np.r_[lamr, lamr[1:N - N // 2][::-1]]
    print(f"  1D ring, {'nearest neighbour' if al is None else f'Levy rates r^-(1+{al})':<22} d_s = " + " ".join(f"{x:5.2f}" for x in ds_from_P(lambda t: P1(lamr, t), (10, 30, 100, 300))) + ("" if al is None else f"   (theory 2/alpha = {2/al:.2f})"))

hdr("G2-2/G2-5 local estimators: comb backbone site; Sierpinski corner")
def local(G, x, ts_, rs_):
    nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
    A = nx.to_scipy_sparse_array(G); L = (sp.diags(np.asarray(A.sum(1)).ravel()) - A).tocsr()
    p0 = np.zeros(len(nodes)); p0[idx[x]] = 1
    grid = np.unique(np.r_[0, np.geomspace(min(ts_) / 1.05, max(ts_) * 1.05, 200)])
    out = []
    for t in ts_:
        a = expm_multiply(-L, p0, start=0, stop=t / 1.05, num=2, endpoint=True)[-1][idx[x]]
        b = expm_multiply(-L, p0, start=0, stop=t * 1.05, num=2, endpoint=True)[-1][idx[x]]
        out.append(-2 * np.log(b / a) / (2 * np.log(1.05)))
    dist = nx.single_source_shortest_path_length(G, x, cutoff=max(rs_) + 1); Nr = np.cumsum(np.bincount(list(dist.values())))
    dg = [np.log(Nr[r + 1] / Nr[r - 1]) / np.log((r + 1) / (r - 1)) for r in rs_]
    return out, dg
def comb(back, tooth):
    G = nx.Graph()
    for i in range(back):
        if i: G.add_edge((i, 0), (i - 1, 0))
        for j in range(1, tooth + 1): G.add_edge((i, j), (i, j - 1))
    return G
Gc = comb(301, 300)
ds_c, dg_c = local(Gc, (150, 0), (10, 30, 100, 300, 1000), (10, 30, 60, 100))
print("  comb 301x300, backbone centre: local d_s(t=10..1000) = " + " ".join(f"{x:.3f}" for x in ds_c) + "  (theory 3/2);  d_g(r=10,30,60,100) = " + " ".join(f"{x:.3f}" for x in dg_c) + "  (theory 2)")
ds_t, dg_t = local(Gc, (150, 150), (10, 30, 100, 300, 1000), (10, 30, 60, 100))
print("  comb, mid-tooth site:          local d_s(t=10..1000) = " + " ".join(f"{x:.3f}" for x in ds_t) + ";  d_g = " + " ".join(f"{x:.3f}" for x in dg_t))
import importlib.util, sys
spec = importlib.util.spec_from_file_location("g2", "s2_g2_graph.py")
src = open("s2_g2_graph.py").read(); ns = {}
exec(src.split("def comb(")[0].split("def sierpinski(")[0], ns)  # imports only
exec("def sierpinski" + src.split("def sierpinski")[1].split("def comb(")[0], ns)
Gs = ns["sierpinski"](8); corner = (0.0, 0.0)
dist = nx.single_source_shortest_path_length(Gs, corner); Nr = np.cumsum(np.bincount(list(dist.values())))
print("  Sierpinski gen 8 corner growth d_g(r=2^k, k=3..7): " + " ".join(f"{np.log(Nr[2**(k+1)]/Nr[2**k])/np.log(2):.3f}" for k in range(3, 7)) + "  (theory 1.585)")
ds_s, _ = local(Gs, corner, (10, 30, 100, 300, 1000), (4,))
print("  Sierpinski gen 8 corner local d_s(t=10..1000): " + " ".join(f"{x:.3f}" for x in ds_s) + "  (theory 1.365; log-periodic oscillation expected)")
