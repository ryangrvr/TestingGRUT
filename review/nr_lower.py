"""INDEPENDENT REPRODUCTION of lower-priority SCOUT-2 load-bearing conclusions (NR-L1..L8). No SCOUT-2 probe code.
Independent code path, not independent reviewer."""
import itertools as it
import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
PA = [X, Y, Z]
def op(n, d):
    out = np.eye(1, dtype=complex)
    for q in range(n): out = np.kron(out, d.get(q, I2))
    return out
GR = {}
def S_vn(r):
    e = np.linalg.eigvalsh((r + r.conj().T) / 2); e = e[e > 1e-13]; return float(-(e * np.log2(e)).sum())
def rdm(psi, n, keep):
    t = psi.reshape([2] * n); rest = [q for q in range(n) if q not in keep]
    A = np.transpose(t, list(keep) + rest).reshape(2 ** len(keep), -1); return A @ A.conj().T
rng = np.random.default_rng(31415)

# ===================================================================== NR-L1
hdr("NR-L1 S2-1 individuation: 2-local qubit Hamiltonians; local fibre of the spectrum map modulo LU (Hellmann-Feynman Jacobian)")
def basis2(n):
    B = [op(n, {i: P}) for i in range(n) for P in PA] + [op(n, {i: P, j: Q}) for i, j in it.combinations(range(n), 2) for P in PA for Q in PA]
    return B
def fibre(n, seed):
    B = basis2(n); r = np.random.default_rng(seed); th = r.normal(size=len(B)); H = sum(t * b for t, b in zip(th, B))
    ev, V = np.linalg.eigh(H)
    J = np.array([[np.real(np.vdot(V[:, k], b @ V[:, k])) for b in B] for k in range(len(ev))])
    rank = np.linalg.matrix_rank(J, tol=1e-8 * np.abs(J).max())
    return len(B), rank, len(B) - rank - 3 * n
for n in (3, 4, 5, 8):
    P, rk, fd = fibre(n, n)
    print(f"  n={n}: parameters {P}, rank of spectrum Jacobian {rk}, local fibre dimension modulo LU (P - rank - 3n) = {fd}")
    GR.setdefault("L1-fibre", []).append((n, fd))
# explicit same-spectrum, LU-inequivalent pair at n = 3
n = 3; B = basis2(n); th0 = np.random.default_rng(7).normal(size=len(B))
Hm = lambda th: sum(t * b for t, b in zip(th, B))
ev0 = np.linalg.eigvalsh(Hm(th0))
def jac(th):
    e, V = np.linalg.eigh(Hm(th)); return e, np.array([[np.real(np.vdot(V[:, k], b @ V[:, k])) for b in B] for k in range(len(e))])
G = [op(n, {i: P}) for i in range(n) for P in PA]           # LU generators
H0 = Hm(th0); lu_tan = np.array([[np.real(np.trace(b.conj().T @ (1j * (g @ H0 - H0 @ g)))) / 2 ** n for b in B] for g in G])
e, J = jac(th0); null = np.linalg.svd(J)[2][np.linalg.matrix_rank(J):]
proj = null - (null @ np.linalg.pinv(lu_tan) @ lu_tan)      # remove LU tangent components
d = proj[np.argmax(np.linalg.norm(proj, axis=1))]; d /= np.linalg.norm(d)
th = th0 + 0.8 * d
for _ in range(60):
    e, J = jac(th); th = th - np.linalg.pinv(J) @ (e - ev0)
def profile(th):
    """LU-invariant support profile: sorted single-site weights and sorted pair weights"""
    w = np.array(th) ** 2; s = [w[3 * i:3 * i + 3].sum() for i in range(n)]; pr = [w[3 * n + 9 * k: 3 * n + 9 * k + 9].sum() for k in range(3)]
    return sorted(np.round(s, 6)), sorted(np.round(pr, 6))
print(f"  explicit pair (n=3): spectrum difference {np.abs(np.linalg.eigvalsh(Hm(th)) - ev0).max():.1e}; LU/permutation-invariant support profiles"
      f" differ: {profile(th) != profile(th0)} (pair weights {profile(th0)[1]} vs {profile(th)[1]})")
GR["L1-pair"] = np.abs(np.linalg.eigvalsh(Hm(th)) - ev0).max() < 1e-9 and profile(th) != profile(th0)
print("  NOTE: zero LOCAL fibre dimension (n >= 8) is local uniqueness near the sample; it is NOT global uniqueness (CPR scope, IR-01).")

# ===================================================================== NR-L2
hdr("NR-L2 S2-1b locality criteria: local-dimension conflict, fixed-d frame disagreement, stability measure dependence")
def blocks_of(g, n): return {q: i for i, b in enumerate(g) for q in b}
def pauli_terms(H, n, tol=1e-9):
    terms = {}
    for lab in it.product(range(4), repeat=n):
        if not any(lab): continue
        P = op(n, {q: [I2, X, Y, Z][a] for q, a in enumerate(lab) if a})
        c = np.real(np.trace(P @ H)) / 2 ** n
        if abs(c) > tol: terms[lab] = c
    return terms
def crit(H, n, g, t0=0.3):
    T = pauli_terms(H, n); blk = blocks_of(g, n)
    sup = [frozenset(blk[q] for q, a in enumerate(l) if a) for l in T]; wts = [c ** 2 for c in T.values()]
    k = max(len(s) for s in sup); pairs = {frozenset(p) for s in sup for p in it.combinations(sorted(s), 2)}
    present = set(sup)
    mdl = sum(1 for l in it.product(range(4), repeat=n) if any(l) and any(frozenset(blk[q] for q, a in enumerate(l) if a) <= s for s in present))
    auto = sum(w for s, w in zip(sup, wts) if len(s) == 1) / sum(wts)
    U = expm(-1j * H * t0); lr = []
    for a, b in it.combinations(range(len(g)), 2):
        acc = []
        for _ in range(6):   # random traceless local operators (independent of the original all-Pauli average)
            def rand_local(bl):
                Mk = rng.normal(size=(2 ** len(bl), 2 ** len(bl))) + 1j * rng.normal(size=(2 ** len(bl), 2 ** len(bl))); Mk = Mk + Mk.conj().T
                Mk -= np.trace(Mk) / len(Mk) * np.eye(len(Mk)); Mk /= np.linalg.norm(Mk, 2)
                full = np.kron(Mk, np.eye(2 ** (n - len(bl)))); order = list(bl) + [q for q in range(n) if q not in bl]; perm = np.argsort(order)
                return full.reshape([2] * 2 * n).transpose(list(perm) + [p + n for p in perm]).reshape(2 ** n, 2 ** n)
            A = U.conj().T @ rand_local(g[a]) @ U; Bo = rand_local(g[b]); C = A @ Bo - Bo @ A; acc.append(np.linalg.norm(C, 2))
        lr.append(np.mean(acc))
    nb = len(g)
    return dict(k=k, edges=len(pairs) / (nb * (nb - 1) / 2), mdl=mdl, auto=auto, lr=float(np.mean(lr)))
n = 4; h12 = (lambda A: A + A.conj().T)(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))); h34 = (lambda A: A + A.conj().T)(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
HA = np.kron(h12, np.eye(4)) + np.kron(np.eye(4), h34) + 0.05 * op(4, {1: X, 2: X})
cq = crit(HA, 4, [[0], [1], [2], [3]]); cQ = crit(HA, 4, [[0, 1], [2, 3]])
print(f"  Case A (h12 + h34 + 0.05 X2X3): qubits {cq}\n                                    ququarts {cQ}")
fine_wins = [k for k in ("mdl",) if cq[k] < cQ[k]] + [k for k in ("edges",) if cq["edges"] < cQ["edges"]]
coarse_wins = [k for k in ("auto",) if cQ[k] > cq[k]] + [k for k in ("lr",) if cQ["lr"] < cq["lr"]]
print(f"  -> criteria favouring FINER factors: {fine_wins}; favouring COARSER factors: {coarse_wins}")
GR["L2-d"] = bool(fine_wins) and bool(coarse_wins)
# fixed-d frame disagreement: n=5, 3-body ZZZ chain + X fields, random H/S/CZ Clifford frames
n = 5; HB = sum(op(n, {i: Z, i + 1: Z, i + 2: Z}) * (0.7 + 0.2 * i) for i in range(n - 2)) + sum((0.8 + 0.1 * i) * op(n, {i: X}) for i in range(n))
Hd = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2); Sg = np.diag([1, 1j])
def cliff(r, depth):
    U = np.eye(2 ** n, dtype=complex)
    for _ in range(depth):
        k = r.integers(3)
        if k == 0: a, b = r.choice(n, 2, replace=False); U = (op(n, {a: np.diag([1, 0]).astype(complex)}) + op(n, {a: np.diag([0, 1]).astype(complex), b: Z})) @ U
        elif k == 1: U = op(n, {int(r.integers(n)): Hd}) @ U
        else: U = op(n, {int(r.integers(n)): Sg}) @ U
    return U
qub = [[q] for q in range(n)]; base = None; found = None; r2 = np.random.default_rng(99)
def exact_crit(H):
    T = pauli_terms(H, n); sup = [frozenset(q for q, a in enumerate(l) if a) for l in T]
    k = max(len(s) for s in sup); pairs = {frozenset(p) for s in sup for p in it.combinations(sorted(s), 2)}
    mdl = sum(1 for l in it.product(range(4), repeat=n) if any(l) and any(frozenset(q for q, a in enumerate(l) if a) <= s for s in set(sup)))
    return dict(k=k, edges=len(pairs), mdl=mdl)
base = exact_crit(HB)
for trial in range(400):
    U = cliff(r2, int(r2.integers(1, 9))); c2 = exact_crit(U.conj().T @ HB @ U)
    better = {k for k in c2 if c2[k] < base[k]}; worse = {k for k in c2 if c2[k] > base[k]}
    if better and worse: found = (c2, better, worse); break
print(f"  fixed d=2, n=5 (ZZZ chain + X): computational {base}; " + (f"Clifford frame after {trial+1} draws {found[0]}: better on {sorted(found[1])}, worse on {sorted(found[2])}" if found else "no strict disagreement in 400 draws"))
GR["L2-frame"] = found is not None
# stability measure dependence (Case A)
G = (lambda A: A + A.conj().T)(rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16))); G *= np.linalg.norm(HA) / np.linalg.norm(G)
Vl = sum(rng.normal() * op(4, {i: PA[a], j: PA[b]}) for i, j in it.combinations(range(4), 2) for a in range(3) for b in range(3)); Vl *= np.linalg.norm(HA) / np.linalg.norm(Vl)
for nm, Vp in (("GUE (TPS-neutral)", G), ("2-local in the qubit TPS", Vl)):
    Hp = HA + 1e-3 * Vp; kq = crit(Hp, 4, [[0], [1], [2], [3]])["k"]; kQ = crit(Hp, 4, [[0, 1], [2, 3]])["k"]
    print(f"  stability: perturbation {nm:<24} eps=1e-3: exact interaction order qubits k={kq}, ququarts k={kQ}")
    GR.setdefault("L2-stab", []).append(kq)
GR["L2-stab"] = GR["L2-stab"][0] == 4 and GR["L2-stab"][1] == 2

# ===================================================================== NR-L3
hdr("NR-L3 S2-2 composition discriminator: discrete polytope symmetries vs continuous entangling quantum dynamics")
sq = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float); sq[:, 2] = 1   # gbit square in (x, y, 1)
cnt = 0
for perm in it.permutations(range(4)):
    M = np.linalg.lstsq(sq[:3], sq[list(perm)][:3], rcond=None)[0]
    if np.allclose(sq @ M, sq[list(perm)]): cnt += 1
print(f"  gbit (square) state space: vertex permutations extending to linear maps = {cnt} (finite dihedral group); a reversible map of a")
print("  polytope permutes its finitely many vertices, so the reversible group is finite and contains no continuous one-parameter")
print("  family: no continuous interaction. (GMCD primary result: boxworld reversible dynamics are trivial; scoped, IR overlay.)")
def concurrence_pure(psi):
    yy = np.kron(Y, Y); return abs(np.vdot(psi, yy @ np.conj(psi)))
U_real = expm(-1j * 0.5 * np.kron(Y, X)); print(f"  real QM: exp(-i t Y(x)X) is real orthogonal: {np.allclose(U_real.imag, 0)}; concurrence of U|00>: {concurrence_pure(U_real @ np.array([1, 0, 0, 0], complex)):.3f}")
plus = np.ones(4, complex) / 2; U_c = expm(-1j * 0.5 * np.kron(Z, Z))
print(f"  complex QM: exp(-i t Z(x)Z) on |++>: concurrence {concurrence_pure(U_c @ plus):.3f}")
GR["L3"] = cnt == 8 and concurrence_pure(U_real @ np.array([1, 0, 0, 0], complex)) > 0.5 and concurrence_pure(U_c @ plus) > 0.5

# ===================================================================== NR-L4
hdr("NR-L4 S2-3/3b: measure selection vs preparation forgetting")
nb = 400; M = 20000
bits = (rng.random((M, nb)) < 0.5).astype(np.int8)
x = np.sqrt(rng.random(M))     # density 2x on [0,1): first 50 bits from x, rest random
for j in range(50): x *= 2; bits[:, j] = (x >= 1); x -= bits[:, j]
freq = [bits[:, t].mean() for t in (0, 2, 4, 8, 16)]
bern = (rng.random((M, nb)) < 0.3).astype(np.int8); fbern = [bern[:, t].mean() for t in (0, 8, 16)]
print(f"  doubling map (exact bit shift): a.c. density-2x preparation, P(x_t >= 1/2) at t = 0,2,4,8,16: " + " ".join(f"{f:.3f}" for f in freq)
      + f";  Bernoulli(0.3) invariant preparation: " + " ".join(f"{f:.3f}" for f in fbern) + ";  fixed point 0: stays 0")
phi = (np.sqrt(5) - 1) / 2; Tn = 200000
avgs = [np.mean(((x0 + phi * np.arange(Tn)) % 1) < 0.5) for x0 in (0.0, 0.123456, 1 / 3, 0.5)]
ens = rng.random(5000) * 0.1; ensf = [np.mean(((ens + phi * t) % 1) < 0.5) for t in (0, 1, 2, 3, 10, 100, 1000)]
print(f"  golden rotation (uniquely ergodic): time averages from 4 starts: " + " ".join(f"{a:.4f}" for a in avgs)
      + f";  ensemble on [0, 0.1): P(<1/2) at t = 0,1,2,3,10,100,1000: " + " ".join(f"{e:.2f}" for e in ensf))
def lagrange(u, v):
    while True:
        if u @ u > v @ v: u, v = v, u
        m = np.round((u @ v) / (u @ u))
        if m == 0: return u, v
        v = v - m * u
def horo_avg(u, v, T=2000.0, dt=0.01):
    hits = 0; cnt = 0
    for t in np.arange(0, T, dt):
        a = np.array([u[0] + t * u[1], u[1]]); b = np.array([v[0] + t * v[1], v[1]])
        s, _ = lagrange(a, b); hits += (s @ s) < 0.5; cnt += 1
    return hits / cnt
al = np.sqrt(2) - 1; ca, sa = np.cos(al), np.sin(al)
print(f"  horocycle on SL(2,R)/SL(2,Z) (own Lagrange reduction): time fraction |v1|^2 < 0.5 from Z^2 rotated by sqrt2-1 = "
      f"{horo_avg(np.array([ca, sa]), np.array([-sa, ca])):.4f} (Haar 3*0.5/pi = {1.5/np.pi:.4f}); from Z^2 (periodic) = {horo_avg(np.array([1.0, 0]), np.array([0, 1.0]), T=50):.4f}")
print("  compact-quotient unique ergodicity is IMPORTED (Furstenberg), not simulated.")
GR["L4"] = abs(freq[-1] - 0.5) < 0.02 and abs(fbern[-1] - 0.3) < 0.02 and max(abs(a - 0.5) for a in avgs) < 1e-3 and (ensf[-1] in (0.0, 1.0) or abs(ensf[-1] - 0.5) > 0.3)

# ===================================================================== NR-L5
hdr("NR-L5 S2-4 convexity: deterministic hidden coin, independent vs dependent tester")
N5 = 10 ** 6; k = np.arange(N5); lam, mu = 0.3, 0.4
c = ((k * phi) % 1) < lam; s_ind = ((k * np.sqrt(2)) % 1) < mu; s_dep = ((2 * k * phi) % 1) < mu
print(f"  P(c=1) = {c.mean():.4f} (= lambda {lam}: the weight is the supplied threshold);  independent tester: P(c,s) = {np.mean(c & s_ind):.4f} vs product {c.mean()*s_ind.mean():.4f};"
      f"  P(c=1 | s=0) = {c[~s_ind].mean():.4f}")
print(f"  rationally dependent tester (2n phi): P(c,s) = {np.mean(c & s_dep):.4f} vs product {c.mean()*s_dep.mean():.4f};  P(c=1 | s=0) = {c[~s_dep].mean():.4f} != lambda")
GR["L5"] = abs(c[~s_ind].mean() - lam) < 2e-3 and abs(c[~s_dep].mean() - lam) > 0.05

# ===================================================================== NR-L6
hdr("NR-L6 S2-5/S2-6 shared reference: J(x)J witness; Ising reference-field coherence (Wolff updates)")
Jm = np.array([[0, -1], [1, 0]], float)
def random_real_sym(d): A = rng.normal(size=(d, d)); return A + A.T
rho = lambda sgn: (np.eye(4) + sgn * np.kron(Jm, Jm)) / 4
w, v = np.linalg.eigh(np.kron(Jm, Jm)); ref_shared = np.outer(v[:, -1], v[:, -1])
ref_prod = np.diag([1.0, 0, 0, 0]); ref_mixed = np.eye(4) / 4
def maxdiff(ref, trials=2000):
    """max |<O_A (x) O_B>_+ - <...>_-| over random local real observables on (A R_A) and (B R_B); plus the J-witness"""
    # ordering: A, B, R_A, R_B -> regroup to (A, R_A, B, R_B)
    def full(sgn): T = np.kron(rho(sgn), ref).reshape([2] * 8); return T.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
    Dp = full(1) - full(-1); best = abs(np.trace(np.kron(np.kron(Jm, Jm), np.kron(Jm, Jm)).T @ Dp))   # witness (J_A J_RA)(x)(J_B J_RB) is real symmetric locally
    for _ in range(trials):
        OA, OB = random_real_sym(4), random_real_sym(4); OA /= np.linalg.norm(OA, 2); OB /= np.linalg.norm(OB, 2)
        best = max(best, abs(np.trace(np.kron(OA, OB) @ Dp)))
    return best
for nm, ref in (("shared correlated reference", ref_shared), ("product reference |00>", ref_prod), ("maximally mixed reference", ref_mixed)):
    print(f"  {nm:<28}: max local distinguishability of rho+ vs rho- = {maxdiff(ref):.3f}")
GR["L6-J"] = maxdiff(ref_shared) > 1.9 and maxdiff(ref_prod) < 1e-12 and maxdiff(ref_mixed) < 1e-12
def wolff(L, T, sweeps):
    s = np.where(rng.random((L, L)) < 0.5, 1, -1); p = 1 - np.exp(-2 / T); out = []
    for sw in range(sweeps):
        i, j = rng.integers(L, size=2); seed = s[i, j]; stack = [(i, j)]; s[i, j] = -seed
        while stack:
            a, b = stack.pop()
            for c_, d_ in (((a + 1) % L, b), ((a - 1) % L, b), (a, (b + 1) % L), (a, (b - 1) % L)):
                if s[c_, d_] == seed and rng.random() < p: s[c_, d_] = -seed; stack.append((c_, d_))
        if sw > sweeps // 2: out.append(np.mean(s[:, 0] * s[:, L // 2]))
    return 2 * np.mean(out)
for T in (1.8, 3.2):
    print(f"  Ising reference field L=32, T={T} (Wolff): 2<f_A f_B> at distance 16 = {wolff(32, T, 3000):.3f}")
print("  scoped reading (Y-06): asymptotic shared coherence in the ordered phase; NOT 'LT iff ordered phase'; full LT not derived (Y-05).")

# ===================================================================== NR-L7
hdr("NR-L7 S2-7 consistent histories: multiplicity, final-basis freedom, exact global consistent set")
n = 5; g7 = rng.uniform(0.8, 1.2, 4); H7 = op(n, {0: X}) + sum(g7[k] * op(n, {0: Z, k + 1: Y}) for k in range(4))
U1, U2 = expm(-1j * H7 * 1.3), expm(-1j * H7 * 2.9); psi0 = np.zeros(2 ** n, complex); psi0[0] = 1
def projs(theta, q):
    v0 = np.array([np.cos(theta / 2), np.sin(theta / 2)]); v1 = np.array([-np.sin(theta / 2), np.cos(theta / 2)])
    return [op(n, {q: np.outer(v, v).astype(complex)}) for v in (v0, v1)]
def decoh(P1, P2):
    C = {(a, b): (U2.conj().T @ P2[b] @ U2) @ (U1.conj().T @ P1[a] @ U1) @ psi0 for a in range(2) for b in range(2)}
    p = {k: np.vdot(v, v).real for k, v in C.items()}; worst = 0; same_b_zero = 0
    for k1, k2 in it.combinations(C, 2):
        val = abs(np.vdot(C[k2], C[k1]))
        if k1[1] != k2[1]: same_b_zero = max(same_b_zero, val)
        elif p[k1] > 1e-10 and p[k2] > 1e-10: worst = max(worst, val / np.sqrt(p[k1] * p[k2]))
    return worst, same_b_zero
grid = np.linspace(0, np.pi, 19); cons = 0; maxbb = 0
for a in grid:
    for b in grid:
        wv, zb = decoh(projs(a, 0), projs(b, 0)); cons += wv < 0.05; maxbb = max(maxbb, zb)
print(f"  S + 4 env, two-time system histories on a 19x19 basis grid: consistent sets (eps 0.05): {cons}; max |D| between different FINAL outcomes: {maxbb:.1e} (identically 0)")
chi = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n); chi /= np.linalg.norm(chi); vt = U1 @ psi0 + 0.5 * chi; vt /= np.linalg.norm(vt)
P1 = [np.outer(vt, vt.conj()), np.eye(2 ** n) - np.outer(vt, vt.conj())]; wv2 = U2 @ U1.conj().T @ (P1[0] @ (U1 @ psi0)); wv2 /= np.linalg.norm(wv2)
P2 = [np.outer(wv2, wv2.conj()), np.eye(2 ** n) - np.outer(wv2, wv2.conj())]
Cg = {(a, b): U2.conj().T @ P2[b] @ U2 @ U1.conj().T @ P1[a] @ U1 @ psi0 for a in range(2) for b in range(2)}
gmax = max(abs(np.vdot(Cg[k2], Cg[k1])) for k1, k2 in it.combinations(Cg, 2))
print(f"  global (non-local, non-record) projector set: max |off-diagonal D| = {gmax:.1e} -> exactly consistent")
GR["L7"] = cons > 1 and maxbb < 1e-12 and gmax < 1e-12

# ===================================================================== NR-L8
hdr("NR-L8 S2-8 Darwinism: pointer given split; split not selected; frame dependence; fragment dependence")
n = 7; H8 = sum(op(n, {0: Z, k: Y}) for k in range(1, n)); p0 = op(n, {0: np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)}) @ np.eye(2 ** n, dtype=complex)[:, 0]
psi = expm(-1j * H8 * np.pi / 4) @ p0
def MI(psi, A, Bs): return S_vn(rdm(psi, n, A)) + S_vn(rdm(psi, n, Bs)) - S_vn(rdm(psi, n, sorted(A + Bs)))
def holevo(psi, basis, frag):
    cond = []; ps = []
    for v in basis:
        P = op(n, {0: np.outer(v, v.conj()).astype(complex)}); phi_ = P @ psi; p = np.vdot(phi_, phi_).real
        cond.append(rdm(phi_ / np.sqrt(p), n, frag)); ps.append(p)
    avg = sum(p * r for p, r in zip(ps, cond)); return S_vn(avg) - sum(p * S_vn(r) for p, r in zip(ps, cond))
chiZ = holevo(psi, [np.array([1, 0]), np.array([0, 1])], [1]); chiX = holevo(psi, [np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)], [1])
print(f"  given S|E: Holevo chi(S basis : E1): Z {chiZ:.3f}, X {chiX:.3f} -> pointer Z")
def R_count(psi, sys_, frags, delta=0.1):
    HS = S_vn(rdm(psi, n, sys_)); need = None
    for j in range(1, len(frags) + 1):
        if np.mean([MI(psi, sys_, sum(c, [])) for c in it.combinations(frags, j)]) >= (1 - delta) * HS: need = j; break
    return len(frags) / need if need else 0.0
Rs = [R_count(psi, [s], [[q] for q in range(n) if q != s]) for s in range(n)]
print(f"  redundancy with slot j as 'system', j = 0..6: {Rs}  (thresholded R_delta: numerical, not covered by D0 — IR-07)")
Ucl = np.eye(2 ** n, dtype=complex)
for kk in range(1, n): Ucl = (op(n, {0: np.diag([1, 0]).astype(complex)}) + op(n, {0: np.diag([0, 1]).astype(complex), kk: X})) @ Ucl
ghz = np.zeros(2 ** n, complex); ghz[0] = ghz[-1] = 1 / np.sqrt(2)
R_ghz = R_count(ghz, [0], [[q] for q in range(1, n)]); psi_f = Ucl @ ghz
HSf = S_vn(rdm(psi_f, n, [0]))
print(f"  same global GHZ state: frame 1 R = {R_ghz}; CNOT-ladder frame: H(S) = {HSf:.3f} -> R = 0 (no information to record)")
env = list(range(1, n))
for size in (1, 2, 3):
    frags = [env[i:i + size] for i in range(0, len(env), size)]
    print(f"  fragments of {size} qubit(s): R = {R_count(psi, [0], frags)}")
GR["L8"] = chiZ > 0.9 and chiX < 1e-9 and len(set(Rs)) == 1 and R_ghz > 0 and HSf < 1e-9

hdr("SUMMARY (machine checks)")
for k, v in GR.items(): print(f"  {k}: {v}")
