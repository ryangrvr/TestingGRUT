"""SCOUT-2 S2-G2: graph / spectral dimension.  Estimators kept separate:
  d_g  growth:    N(r) ~ r^{d_g}              (BFS balls, averaged over centres; local slope)
  d_s  spectral:  P(t) = (1/N) sum_{mu>0} e^{-t mu} ~ t^{-d_s/2}  (normalized Laplacian = standard random walk)
  d_w  walk:      <r^2(t)> ~ t^{2/d_w}        (heat-kernel second moment in graph distance)
"""
import itertools as it
import numpy as np
import networkx as nx
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

rng = np.random.default_rng(11)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def lap_sym(W):
    W = sp.csr_matrix(W, dtype=float); d = np.asarray(W.sum(1)).ravel()
    Dm = sp.diags(1 / np.sqrt(d)); return (sp.eye(W.shape[0]) - Dm @ W @ Dm).toarray(), d

def spectral_ds(mu, N, ts):
    mu = np.sort(mu)[1:]                                   # drop the zero mode (finite-size plateau)
    P = lambda t: np.exp(-t * mu).sum() / N
    return [-2 * (np.log(P(t * 1.1)) - np.log(P(t / 1.1))) / (np.log(1.1) * 2) for t in ts]

def growth(G, rs, centres=12):
    nodes = list(G.nodes()); cs = [nodes[i] for i in rng.choice(len(nodes), min(centres, len(nodes)), replace=False)]
    Nr = np.zeros(max(rs) + 2)
    for c in cs:
        dist = nx.single_source_shortest_path_length(G, c, cutoff=max(rs) + 1)
        cnt = np.bincount(list(dist.values()), minlength=max(rs) + 2)[: max(rs) + 2]
        Nr += np.cumsum(cnt)
    Nr /= len(cs)
    return [np.log(Nr[r + 1] / Nr[r - 1]) / np.log((r + 1) / (r - 1)) for r in rs]

def walk_dw(G, W, ts, starts=4):
    nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
    W = sp.csr_matrix(W, dtype=float); d = np.asarray(W.sum(1)).ravel()
    Lrw_T = (sp.eye(W.shape[0]) - sp.diags(1 / d) @ W).T.tocsr()
    m2 = np.zeros(len(ts))
    for s in rng.choice(len(nodes), starts, replace=False):
        dist = nx.single_source_shortest_path_length(G, nodes[s]); r2 = np.zeros(len(nodes))
        for v, dv in dist.items(): r2[idx[v]] = dv ** 2
        p0 = np.zeros(len(nodes)); p0[s] = 1
        P = expm_multiply(-Lrw_T, p0, start=0, stop=ts[-1], num=len(ts) * 0 + 400, endpoint=True)
        tt = np.linspace(0, ts[-1], 400)
        for i, t in enumerate(ts):
            m2[i] += P[np.argmin(abs(tt - t))] @ r2
    m2 /= starts
    sl = [np.log(m2[i + 1] / m2[i]) / np.log(ts[i + 1] / ts[i]) for i in range(len(ts) - 1)]
    return [2 / s if s > 1e-9 else np.inf for s in sl]

def sierpinski(gen):
    tri = [((0.0, 0.0), (1.0, 0.0), (0.5, np.sqrt(3) / 2))]
    for _ in range(gen):
        new = []
        for a, b, c in tri:
            ab = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2); bc = ((b[0] + c[0]) / 2, (b[1] + c[1]) / 2); ca = ((c[0] + a[0]) / 2, (c[1] + a[1]) / 2)
            new += [(a, ab, ca), (ab, b, bc), (ca, bc, c)]
        tri = new
    G = nx.Graph(); key = lambda p: (round(p[0], 9), round(p[1], 9))
    for a, b, c in tri:
        G.add_edges_from([(key(a), key(b)), (key(b), key(c)), (key(c), key(a))])
    return G

def comb(back, tooth):
    G = nx.Graph()
    for i in range(back):
        if i: G.add_edge((i, 0), (i - 1, 0))
        for j in range(1, tooth + 1): G.add_edge((i, j), (i, j - 1))
    return G

def analyse(name, G, W=None, ts=(2, 5, 10, 20, 50, 100), rs=(3, 5, 8, 12), dw_ts=(5, 10, 20, 40, 80)):
    W = nx.to_scipy_sparse_array(G, weight="weight") if W is None else W
    L, _ = lap_sym(W); mu = np.linalg.eigvalsh(L)
    ds = spectral_ds(mu, len(mu), ts); dg = growth(G, rs); dw = walk_dw(G, W, list(dw_ts))
    print(f"  {name:<34} N={G.number_of_nodes():5d} | d_g(r={rs}): " + " ".join(f"{x:5.2f}" for x in dg)
          + f" | d_s(t={ts}): " + " ".join(f"{x:5.2f}" for x in ds) + " | d_w: " + " ".join(f"{x:5.2f}" for x in dw))
    return mu

hdr("G2-1 positive controls and G2-2 same-N-order, different connectivity")
analyse("1D ring", nx.cycle_graph(2000))
analyse("2D torus 45x45", nx.grid_2d_graph(45, 45, periodic=True))
analyse("3D torus 13^3", nx.grid_graph(dim=[13, 13, 13], periodic=True), rs=(2, 3, 4, 5))
analyse("random 3-regular", nx.random_regular_graph(3, 2000, seed=1), rs=(2, 3, 4, 5))
analyse("small world WS(k=2,p=0.05)", nx.connected_watts_strogatz_graph(2000, 2, 0.05, seed=2))
analyse("binary tree depth 10", nx.balanced_tree(2, 10), rs=(2, 3, 4, 5))
analyse("Sierpinski gasket gen 7", sierpinski(7), rs=(4, 8, 16, 24), ts=(5, 10, 20, 50, 100, 200), dw_ts=(10, 20, 40, 80, 160))
analyse("comb 45 x 45", comb(45, 44), ts=(5, 10, 20, 50, 100, 200))
print("  known: Sierpinski d_f = log3/log2 = 1.585, d_s = 2 log3/log5 = 1.365, d_w = log5/log2 = 2.322;  comb: d_g = 2, d_s = 3/2")

hdr("G2-3 scale dependence: crossovers")
def aniso_torus(L, eps, dims=2):
    G = nx.grid_graph(dim=[L] * dims, periodic=True)
    for u, v in G.edges():
        diff = [a != b for a, b in zip(u, v)]
        G[u][v]["weight"] = 1.0 if (dims == 2 and diff[0]) or (dims == 3 and not diff[2]) else eps
    return G
ts_long = (2, 5, 10, 30, 100, 300, 1000, 3000)
for nm, G in [("2D bundled chains, eps = 0.01", aniso_torus(60, 0.01)), ("3D layered, interlayer eps = 0.01", aniso_torus(14, 0.01, 3))]:
    W = nx.to_scipy_sparse_array(G, weight="weight"); L, _ = lap_sym(W); mu = np.linalg.eigvalsh(L)
    print(f"  {nm:<34} d_s(t={ts_long}): " + " ".join(f"{x:5.2f}" for x in spectral_ds(mu, len(mu), ts_long)))
for p in [0.0, 0.01, 0.1]:
    G = nx.connected_watts_strogatz_graph(2000, 2, p, seed=3) if p else nx.cycle_graph(2000)
    W = nx.to_scipy_sparse_array(G); L, _ = lap_sym(W); mu = np.linalg.eigvalsh(L)
    print(f"  ring + small-world rewiring p={p:<5} d_s(t={ts_long[:6]}): " + " ".join(f"{x:5.2f}" for x in spectral_ds(mu, len(mu), ts_long[:6])))

hdr("G2-4 diffusion-operator hostile (same vertex set / same graph)")
G = nx.grid_2d_graph(45, 45, periodic=True); A = nx.to_scipy_sparse_array(G)
L, _ = lap_sym(A); mu = np.linalg.eigvalsh(L)
ts = (5, 10, 20, 50, 100)
print(f"  2D torus, standard walk          d_s: " + " ".join(f"{x:5.2f}" for x in spectral_ds(mu, len(mu), ts)))
Gw = G.copy()
for u, v in Gw.edges(): Gw[u][v]["weight"] = rng.uniform(0.5, 1.5)
Lw, _ = lap_sym(nx.to_scipy_sparse_array(Gw, weight="weight")); muw = np.linalg.eigvalsh(Lw)
print(f"  2D torus, random weights U(.5,1.5) d_s: " + " ".join(f"{x:5.2f}" for x in spectral_ds(muw, len(muw), ts)))
Ga = aniso_torus(45, 0.01); La, _ = lap_sym(nx.to_scipy_sparse_array(Ga, weight="weight")); mua = np.linalg.eigvalsh(La)
print(f"  2D torus, anisotropic eps=0.01    d_s: " + " ".join(f"{x:5.2f}" for x in spectral_ds(mua, len(mua), ts)))
for alpha in [1.0, 0.5]:
    mf = np.clip(mu, 0, None) ** (alpha / 2)
    print(f"  2D torus, fractional Delta^{alpha/2:.2f}    d_s: " + " ".join(f"{x:5.2f}" for x in spectral_ds(mf, len(mf), ts)) + f"   (theory 2d/alpha = {4/alpha:.0f})")
Nr = 1500; idx = np.arange(Nr)
for alpha in [None, 1.5, 1.0]:
    if alpha is None:
        Wr = nx.to_scipy_sparse_array(nx.cycle_graph(Nr)).toarray().astype(float); lab = "nearest neighbour"
    else:
        dd = np.abs(idx[:, None] - idx[None, :]); dd = np.minimum(dd, Nr - dd).astype(float); np.fill_diagonal(dd, np.inf)
        Wr = dd ** (-(1 + alpha)); lab = f"Levy rates 1/r^(1+{alpha})"
    d = Wr.sum(1); Lr = np.eye(Nr) - Wr / np.sqrt(d[:, None] * d[None, :]); mur = np.linalg.eigvalsh(Lr)
    print(f"  1D ring, {lab:<26} d_s: " + " ".join(f"{x:5.2f}" for x in spectral_ds(mur, Nr, (5, 10, 20, 50, 100))) +
          ("" if alpha is None else f"   (theory 2/alpha = {2/alpha:.2f})"))

hdr("G2-5 structural hostiles")
print("  (a) d_s(t) depends ONLY on the normalized-Laplacian spectrum: isospectral non-isomorphic graphs have identical")
print("      d_s(t) at every t. Search over all connected graphs with 7 nodes (networkx atlas):")
atlas = [g for g in nx.graph_atlas_g() if g.number_of_nodes() == 7 and nx.is_connected(g)]
spec = {}
for g in atlas:
    L, _ = lap_sym(nx.to_scipy_sparse_array(g)); key = tuple(np.round(np.linalg.eigvalsh(L), 8)); spec.setdefault(key, []).append(g)
pairs = [(v[0], v[1]) for v in spec.values() if len(v) > 1]
diff = []
for a, b in pairs:
    ga = sorted(np.bincount(list(nx.single_source_shortest_path_length(a, n).values())).tolist() for n in a)
    gb = sorted(np.bincount(list(nx.single_source_shortest_path_length(b, n).values())).tolist() for n in b)
    if ga != gb: diff.append((a, b))
print(f"      connected 7-node graphs: {len(atlas)}; normalized-Laplacian-cospectral non-isomorphic classes: {len(pairs)};"
      f" of which with DIFFERENT ball-growth profiles: {len(diff)}")
if diff:
    a, b = diff[0]
    print(f"      example: edges A = {sorted(a.edges())}\n               edges B = {sorted(b.edges())}")
print("  (b) same growth, different diffusion: comb (d_g ~ 2, d_s ~ 1.5) vs square torus (d_g = d_s = 2): see G2-2 rows.")
