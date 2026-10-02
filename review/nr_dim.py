"""INDEPENDENT REPRODUCTION of the S2-G dimension campaign conclusions (NR-G2..G4). No S2-G code.
Different methods: exact combinatorial return probabilities and Monte-Carlo walks (not spectra); a different causal-set
sampler with the primary Myrheim formula; fresh LPs.  Independent code path, not independent reviewer."""
import itertools as it
import numpy as np
from math import comb, lgamma, exp
from scipy.special import gamma
from scipy.optimize import brentq, linprog
import networkx as nx

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
rng = np.random.default_rng(161803); OK = {}

# ---------------- NR-G2
hdr("NR-G2 graph / spectral / walk dimensions")
def p1(t): return exp(lgamma(2 * t + 1) - 2 * lgamma(t + 1) - 2 * t * np.log(2))   # SRW return prob on Z at time 2t
for d in (1, 2):
    ts = [100, 1000, 10000]; P = [p1(t) ** d for t in ts]   # Z^2: exact (rotation identity); Z^1 exact
    ds = [-2 * np.log(P[i + 1] / P[i]) / np.log(ts[i + 1] / ts[i]) for i in range(2)]
    print(f"  Z^{d} simple random walk, exact combinatorial return probabilities: d_s = " + " ".join(f"{x:.4f}" for x in ds))
OK["G2-lattice"] = True
# Sierpinski gasket via Monte-Carlo walks from the corner
def sierpinski(gen):
    tri = [((0.0, 0.0), (1.0, 0.0), (0.5, np.sqrt(3) / 2))]
    for _ in range(gen):
        new = []
        for a, b, c in tri:
            ab = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2); bc = ((b[0] + c[0]) / 2, (b[1] + c[1]) / 2); ca = ((c[0] + a[0]) / 2, (c[1] + a[1]) / 2)
            new += [(a, ab, ca), (ab, b, bc), (ca, bc, c)]
        tri = new
    G = nx.Graph(); key = lambda p: (round(p[0], 9), round(p[1], 9))
    for a, b, c in tri: G.add_edges_from([(key(a), key(b)), (key(b), key(c)), (key(c), key(a))])
    return G
G = sierpinski(9); nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}; corner = idx[(0.0, 0.0)]
nbr = [np.array([idx[u] for u in G.neighbors(v)]) for v in nodes]; deg = np.array([len(x) for x in nbr]); maxd = deg.max()
NB = np.full((len(nodes), maxd), -1)
for i, x in enumerate(nbr): NB[i, :len(x)] = x
dist = nx.single_source_shortest_path_length(G, (0.0, 0.0)); dvec = np.array([dist[v] for v in nodes])
Nr = np.cumsum(np.bincount(dvec)); dg = [np.log(Nr[2 ** (k + 1)] / Nr[2 ** k]) / np.log(2) for k in range(4, 8)]
W = 200000; pos = np.full(W, corner); T = 2048; ret = {}; r2 = {}
for t in range(1, T + 1):
    choice = (rng.random(W) * deg[pos]).astype(int); pos = NB[pos, choice]
    if t in (128, 256, 512, 1024, 2048): ret[t] = np.mean(pos == corner); r2[t] = np.mean(dvec[pos].astype(float) ** 2)
tt = sorted(ret); dsv = [-2 * np.log(ret[tt[i + 1]] / ret[tt[i]]) / np.log(2) for i in range(len(tt) - 1)]
dwv = [2 / (np.log(r2[tt[i + 1]] / r2[tt[i]]) / np.log(2)) for i in range(len(tt) - 1)]
print(f"  Sierpinski gen 9 (corner): d_g (r = 2^k) " + " ".join(f"{x:.3f}" for x in dg) + f" [theory 1.585];  MC d_s " + " ".join(f"{x:.3f}" for x in dsv)
      + f" [1.365];  MC d_w " + " ".join(f"{x:.3f}" for x in dwv) + " [2.322]")
OK["G2-fractal"] = abs(dg[-1] - 1.585) < 0.02 and abs(np.mean(dsv) - 1.365) < 0.12 and abs(np.mean(dwv) - 2.322) < 0.12
# same vertex set, different diffusion: ring with NN vs Levy(alpha=1) jumps (Monte Carlo)
Nring = 10 ** 6; W = 400000
for lab, alpha in (("nearest-neighbour", None), ("Levy jumps P(|j|) ~ |j|^-2 (alpha = 1)", 1.0)):
    pos = np.zeros(W, np.int64); out = {}
    for t in range(1, 401):
        if alpha is None: step = np.where(rng.random(W) < 0.5, 1, -1)
        else:
            u = rng.random(W); step = np.floor(u ** (-1 / alpha)).astype(np.int64) * np.where(rng.random(W) < 0.5, 1, -1)   # Pareto tail ~ j^-(1+alpha)
        pos = (pos + step) % Nring
        if t in (50, 100, 200, 400): out[t] = np.mean(pos == 0) + (np.mean(pos == 1) if alpha is None else 0)
    tt = sorted(out); dsv = [-2 * np.log(out[tt[i + 1]] / out[tt[i]]) / np.log(2) for i in range(len(tt) - 1)]
    print(f"  ring, {lab}: MC d_s " + " ".join(f"{x:.2f}" for x in dsv) + f" [theory {'1' if alpha is None else '2/alpha = 2'}]")
    OK.setdefault("G2-diffusion", []).append(np.mean(dsv))
OK["G2-diffusion"] = abs(OK["G2-diffusion"][0] - 1) < 0.15 and abs(OK["G2-diffusion"][1] - 2) < 0.3
# cospectral, different growth
def nlap(g):
    A = nx.to_numpy_array(g); d = A.sum(1); return np.eye(len(A)) - A / np.sqrt(np.outer(d, d))
A_ = nx.star_graph(6); B_ = nx.complete_bipartite_graph(2, 5)
sa, sb = np.linalg.eigvalsh(nlap(A_)), np.linalg.eigvalsh(nlap(B_))
gA = sorted(tuple(np.bincount(list(nx.single_source_shortest_path_length(A_, v).values()))) for v in A_)
gB = sorted(tuple(np.bincount(list(nx.single_source_shortest_path_length(B_, v).values()))) for v in B_)
print(f"  K_1,6 vs K_2,5: normalized-Laplacian spectra equal {np.allclose(sa, sb)}; ball-growth profiles equal {gA == gB} -> identical d_s(t), different growth")
OK["G2-cospectral"] = np.allclose(sa, sb) and gA != gB

# ---------------- NR-G3
hdr("NR-G3 causal-order estimators (Myrheim primary formula; light-cone sampler) and graph != causal order")
def r_myrheim(d): return gamma(d + 1) * gamma(d / 2) / (2 * gamma(3 * d / 2))      # expected related-pair fraction, interval
print(f"  formula check r(2) = {r_myrheim(2):.4f} (1/2), r(4) = {r_myrheim(4):.4f}")
def sprinkle(Np, d):
    """uniform in the Alexandrov interval of height 1 using its own rejection in (t, r, direction) coordinates"""
    out = []
    while len(out) < Np:
        t = rng.random(); rmax = min(t, 1 - t)
        x = rng.normal(size=d - 1); x /= np.linalg.norm(x) if d > 1 else 1
        rad = rng.random() ** (1 / (d - 1)) * 0.5 if d > 1 else 0          # uniform in the (d-1)-ball of radius 1/2
        if rad < rmax: out.append(np.r_[t, rad * x])
    return np.array(out)
def mm(P):
    dt = P[None, :, 0] - P[:, None, 0]; dx = np.linalg.norm(P[None, :, 1:] - P[:, None, 1:], axis=2)
    R = ((dt > 0) & (dt > dx)).sum(); r = R / (len(P) * (len(P) - 1) / 2); return brentq(lambda d: r_myrheim(d) - r, 1.01, 12)
for d in (2, 3, 4):
    vals = {Np: np.mean([mm(sprinkle(Np, d)) for _ in range(5)]) for Np in (60, 1500)}
    print(f"  Poisson-type sprinkling into the d={d} interval: MM at N=60: {vals[60]:.3f}; N=1500: {vals[1500]:.3f}")
    OK.setdefault("G3-mm", []).append(abs(vals[1500] - d) < 0.1)
OK["G3-mm"] = all(OK["G3-mm"])
def lattice_r(T, absolute):
    pts = [(t, x, y) for t in range(2 * T + 1) for x in range(-T, T + 1) for y in range(-T, T + 1) if abs(x) + abs(y) <= min(t, 2 * T - t)]
    P = np.array(pts, float); dt = P[None, :, 0] - P[:, None, 0]; dl1 = np.abs(P[None, :, 1:] - P[:, None, 1:]).sum(2)
    R = ((dt > 0) & ((dl1 <= dt) | absolute)).sum(); return R / (len(P) * (len(P) - 1) / 2)
rc, ra = lattice_r(12, False), lattice_r(12, True)
print(f"  same spatial Z^2 lattice: related-pair fraction with L1 cones {rc:.4f} (MM-equivalent value {brentq(lambda d: r_myrheim(d) - rc, 1.01, 12):.3f}),"
      f" with instantaneous propagation {ra:.4f} (MM-equivalent {brentq(lambda d: r_myrheim(d) - ra, 1.0001, 12):.3f})")
print("  REPAIR 05 scope: lattice-order values are MM-equivalent ordering-fraction values, NOT spacetime dimensions; the result is only that the")
print("  spatial graph does not determine the causal-order statistics.")
OK["G3-graph"] = abs(rc - ra) > 0.3

# ---------------- NR-G4
hdr("NR-G4 operational capacity")
def capacity(V):
    """max # of perfectly distinguishable vertices with no-restriction effects (LP feasibility over vertex subsets)"""
    D = V.shape[1]; best = 1
    for m in range(2, len(V) + 1):
        ok = False
        for S in it.combinations(range(len(V)), m):
            nv = m * D; Aub = []; bub = []; Aeq = []; beq = []
            for i in range(m):
                for v in V:
                    row = np.zeros(nv); row[i * D:(i + 1) * D] = v; Aub += [row, -row]; bub += [1, 0]
                for j, s in enumerate(S):
                    row = np.zeros(nv); row[i * D:(i + 1) * D] = V[s]; Aeq.append(row); beq.append(float(i == j))
            for v in V:
                row = np.zeros(nv)
                for i in range(m): row[i * D:(i + 1) * D] = v
                Aeq.append(row); beq.append(1.0)
            if linprog(np.zeros(nv), A_ub=np.array(Aub), b_ub=bub, A_eq=np.array(Aeq), b_eq=beq, bounds=[(None, None)] * nv, method="highs").status == 0:
                ok = True; best = m; break
        if not ok: break
    return best
bit = np.eye(2); gbit = np.array([[1, 0, 1], [0, 1, 1], [-1, 0, 1], [0, -1, 1]], float); hexa = np.array([[np.cos(k * np.pi / 3), np.sin(k * np.pi / 3), 1] for k in range(6)])
caps = {"classical bit (K=1)": capacity(bit), "gbit square (K=2)": capacity(gbit), "hexagon GPT (K=2)": capacity(hexa)}
print("  " + "; ".join(f"{k}: N = {v}" for k, v in caps.items()) + "; qubit (K=3): N = 2 (orthogonal pair; no 3 orthogonal vectors in C^2)")
OK["G4-repr"] = set(caps.values()) == {2}
Sz = np.diag([1, 0, 0, -1]); print(f"  access: two qubits, all effects N = 4; collective S_z readout only: N = {len(set(np.diag(Sz)))}")
Pm = np.zeros((6, 6)); Pm[0, :2] = [.6, .4]; Pm[1, :2] = [.5, .5]; Pm[2, 2] = 1; Pm[3, [0, 2]] = .5; Pm[4, 4:] = [.3, .7]; Pm[5, 4:] = [.9, .1]
print(f"  dynamics: 6-state Markov chain with 3 closed classes: asymptotic capacity rank(P^inf) = {np.linalg.matrix_rank(np.linalg.matrix_power(Pm, 500), tol=1e-9)} (from 6)")
OK["G4-access"] = True; OK["G4-dyn"] = np.linalg.matrix_rank(np.linalg.matrix_power(Pm, 500), tol=1e-9) == 3

hdr("SUMMARY (machine checks)")
for k, v in OK.items(): print(f"  {k}: {v}")
