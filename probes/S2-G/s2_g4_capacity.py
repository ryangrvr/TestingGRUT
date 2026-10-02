"""SCOUT-2 S2-G4: operational capacity N (max # perfectly distinguishable states) vs Hilbert d vs affine K.
Polytope GPTs: N by LP over vertex subsets (no-restriction effects = all functionals in [0,1] on the state space)."""
import itertools as it
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(9)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

def distinguishable(V, S, unit):
    """V: vertices (rows, vectors in R^D); S: indices of candidate states; unit: vector u with u.v = 1 on states.
    Feasibility LP for effects w_1..w_m (in R^D): 0 <= w_i.v <= 1 for all vertices v, sum_i w_i = u (as functionals on
    the span: enforced on all vertices), w_i.s_j = delta_ij."""
    D = V.shape[1]; m = len(S); nv = m * D
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    for i in range(m):
        for v in V:
            row = np.zeros(nv); row[i * D:(i + 1) * D] = v
            A_ub.append(row); b_ub.append(1); A_ub.append(-row); b_ub.append(0)
        for j, sj in enumerate(S):
            row = np.zeros(nv); row[i * D:(i + 1) * D] = V[sj]; A_eq.append(row); b_eq.append(1.0 if i == j else 0.0)
    for v in V:
        row = np.zeros(nv)
        for i in range(m): row[i * D:(i + 1) * D] = v
        A_eq.append(row); b_eq.append(1.0)
    r = linprog(np.zeros(nv), A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq, bounds=[(None, None)] * nv, method="highs")
    return r.status == 0

def capacity(V, unit, mmax, symmetric_first=None):
    best = 1
    for m in range(2, mmax + 1):
        found = False
        for S in it.combinations(range(len(V)), m):
            if distinguishable(V, S, unit): found = True; best = m; break
        if not found: break
    return best

def affine_dim(V):
    return np.linalg.matrix_rank(V[1:] - V[0], tol=1e-9)

hdr("G4-1 matched capacity, different state spaces")
rows = []
for n in range(3, 9):
    V = np.array([[np.cos(2 * np.pi * k / n), np.sin(2 * np.pi * k / n), 1.0] for k in range(n)])
    N = capacity(V, None, 4); rows.append((f"regular {n}-gon GPT" + (" (= classical trit)" if n == 3 else " (gbit)" if n == 4 else ""), N, affine_dim(V), "—"))
for k in (2, 3, 4):
    V = np.eye(k); rows.append((f"classical simplex, {k} outcomes", capacity(V, None, k + 1), affine_dim(V), "—"))
rows += [("rebit (real QM, d=2)", 2, 2, "2 (real)"), ("qubit (complex QM, d=2)", 2, 3, "2"), ("spin factor / ball in R^4", 2, 4, "—"),
         ("qutrit (complex QM, d=3)", 3, 8, "3"), ("real QM, d=3", 3, 5, "3 (real)")]
print(f"  {'system':<36} {'N':>3} {'K':>4}   Hilbert d")
for nm, N, K, d in rows: print(f"  {nm:<36} {N:>3} {K:>4}   {d}")
print("  => N = 2 for: classical bit (K=1), rebit (K=2), gbit / polygons (K=2), qubit (K=3), spin factors (K=4, ...):")
print("     SAME CAPACITY != SAME STATE SPACE; capacity does not fix the representation (C4/C5 needs extra principles).")

hdr("G4-2 capacity from dynamics")
P = np.zeros((7, 7))
P[0, :2] = [.5, .5]; P[1, :2] = [.3, .7]; P[2, 2:4] = [.9, .1]; P[3, 2:4] = [.4, .6]; P[4, 4] = 1; P[5, [0, 4]] = [.5, .5]; P[6, [2, 6]] = [.5, .5]
Pinf = np.linalg.matrix_power(P, 400)
print(f"  Markov chain, 7 states, closed classes {{0,1}},{{2,3}},{{4}} + transients: asymptotic distinguishable outputs = rank(P^inf) = {np.linalg.matrix_rank(Pinf, tol=1e-9)} (initial capacity 7)")
def haar_su2():
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); q, r = np.linalg.qr(z); return q * (np.diag(r) / abs(np.diag(r)))
Us = [np.kron(np.kron(u, u), u) for u in (haar_su2() for _ in range(3))]
M = np.vstack([np.kron(np.eye(8), U) - np.kron(U.T, np.eye(8)) for U in Us])
_, s, Vh = np.linalg.svd(M); comm = Vh[s.size - (s < 1e-9).sum():] if (s < 1e-9).sum() else Vh[-0:]
null = Vh[np.sum(s > 1e-9):].conj()          # null vectors of a complex SVD are conj(rows of Vh)
basis = [v.reshape(8, 8, order="F") for v in null]
dimC = len(basis)
Mc = np.vstack([np.kron(np.eye(8), B) - np.kron(B.T, np.eye(8)) for B in basis])
coef = np.linalg.svd(np.hstack([np.array([b.reshape(-1, order='F') for b in basis]).T]), compute_uv=False)
# centre: elements of the commutant that commute with all of the commutant
G = np.array([b.reshape(-1, order="F") for b in basis]).T          # 64 x dimC
Cmat = np.vstack([(np.kron(np.eye(8), B) - np.kron(B.T, np.eye(8))) @ G for B in basis])
_, s2, Vh2 = np.linalg.svd(Cmat); centre_dim = int((s2 < 1e-8).sum() + (Cmat.shape[1] - len(s2)))
cz = sum(c * b for c, b in zip(Vh2[-1] if centre_dim else [], basis)) if centre_dim else None
Zc = sum(rng.normal() * sum(c * b for c, b in zip(v, basis)) for v in Vh2[len(s2) - centre_dim:].conj())
Zc = (Zc + Zc.conj().T) / 2; ev, evec = np.linalg.eigh(Zc)
groups = []; clusters = []
for e in ev:
    if not clusters or abs(e - clusters[-1][-1]) > 1e-4 * max(1, abs(ev).max()): clusters.append([e])
    else: clusters[-1].append(e)
for cl in clusters:
    sel = (ev >= cl[0] - 1e-9) & (ev <= cl[-1] + 1e-9)
    Pj = evec[:, sel]; Pj = Pj @ Pj.conj().T
    blk = np.linalg.matrix_rank(np.array([(Pj @ b @ Pj).ravel() for b in basis]), tol=1e-8); groups.append(int(round(np.sqrt(blk))))
print(f"  collective SU(2) noise on 3 qubits: commutant dim {dimC}, centre dim {centre_dim}, central blocks {len(clusters)} with sizes "
      f"{[len(c) for c in clusters]}, multiplicities m_J = {groups} (theory: J=3/2 -> 1, J=1/2 -> 2)")
print(f"     -> noiseless capacity (perfectly distinguishable after twirl) = sum m_J = {sum(groups)}  (initial 8): CAPACITY DERIVED FROM D (noise symmetry)")

hdr("G4-3 access hostile: same state space, restricted effects")
print("  2 qubits, all effects: N = 4.  Readout restricted to functions of total S_z (collective access):")
Sz = np.diag([1, 0, 0, -1]); print(f"     N = number of distinct S_z eigenvalues = {len(np.unique(np.diag(Sz)))}")
print("  qutrit (d = 3), readout restricted to {|0><0|, |1><1| + |2><2|}: N = 2   <- G4 SAME as qubit, G1 (d) DIFFERENT")
eta = 0.1
print(f"  qubit, effects restricted to spectrum in [{eta}, {1-eta}] (unsharp access): no effect reaches 0/1 -> N = 1")

hdr("G4-4 noise: exact vs epsilon vs channel capacity (depolarizing qubit, eps = 0.05)")
h = lambda x: 0 if x in (0, 1) else -x * np.log2(x) - (1 - x) * np.log2(1 - x)
print("    p      exact N   eps-capacity (err<=0.05)   Holevo capacity (bits)")
for p in (0.0, 0.01, 0.05, 0.1, 0.3):
    print(f"   {p:<5}    {2 if p == 0 else 1:^7}   {2 if p / 2 <= 0.05 else 1:^24}   {1 - h(p / 2):.4f}")

hdr("G4-5 composition")
print("  classical: N_AB = N_A N_B;  complex QM: N = d_A d_B, K+1 = (K_A+1)(K_B+1);  real QM: N = 4 for 2 rebits but")
print("  K = 9 != (2+1)(2+1)-1 = 8 (local tomography fails) -> capacity multiplicative while state-space dimension is not.")
def det_vertices():
    V = []
    for fa in it.product([0, 1], repeat=2):
        for fb in it.product([0, 1], repeat=2):
            v = np.zeros(16)
            for x, y in it.product([0, 1], repeat=2): v[8 * x + 4 * y + 2 * fa[x] + fb[y]] = 1
            V.append(v)
    return V
def pr_vertices():
    V = []
    for a0, b0, flip in it.product([0, 1], [0, 1], [0, 1]):
        v = np.zeros(16)
        for x, y in it.product([0, 1], repeat=2):
            for a, b in it.product([0, 1], repeat=2):
                if (a ^ b ^ a0 * x ^ b0 * y ^ flip) == (x * y): v[8 * x + 4 * y + 2 * a + b] = 0.5
        V.append(v)
    return V
loc = np.array(det_vertices()); ns = np.vstack([loc, np.array(pr_vertices())])
okns = all(abs(sum(v[8 * x + 4 * y + 2 * a + b] for b in (0, 1)) - sum(v[8 * x + 4 * (1 - y) + 2 * a + b] for b in (0, 1))) < 1e-12 for v in ns for x in (0, 1) for y in (0, 1) for a in (0, 1))
print(f"  two gbits: local polytope {len(loc)} vertices (K={affine_dim(loc)}); no-signalling polytope {len(ns)} vertices (K={affine_dim(ns)}); NS check {okns}")
for nm, V in [("min tensor (local polytope, dual effects)", loc), ("max tensor (boxworld: NS states, NS-dual effects)", ns)]:
    caps = []
    for m in (4, 5):
        hit = any(distinguishable(V, S, None) for S in it.combinations(range(len(V)), m)); caps.append((m, hit))
    print(f"  {nm:<50}: distinguishable sets of size 4: {caps[0][1]}, size 5: {caps[1][1]}")
