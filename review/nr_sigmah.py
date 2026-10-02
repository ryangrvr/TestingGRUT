"""INDEPENDENT REPRODUCTION of S2-SigmaH key claims (NR-SigmaH-1..5).  Shares NO code with probes/S2-Sigma or S2-SigmaH.
Different routines:
  locality  : exact-support weights by Moebius inversion of partial-trace projections (not Pauli transforms)
  entropies : eigvalsh of explicitly traced reduced density matrices (not SVD)
  W_psi     : QR completion of psi (not Householder)
  Cliffords : fresh random H/S/CZ circuits, new seeds (original: H/S/CNOT)
  LU test   : invariance of each block's operator algebra under W (not operator-Schmidt rank)
Run for 3 independent frame seeds."""
import itertools as it
import numpy as np
from scipy.linalg import expm

n = 6; N = 2 ** n
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Z = np.diag([1.0, -1.0]).astype(complex)
Hg = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2); Sg = np.diag([1, 1j])

def kron_list(ms):
    out = np.eye(1, dtype=complex)
    for m in ms: out = np.kron(out, m)
    return out
def op(d): return kron_list([d.get(q, I2) for q in range(n)])
def ptrace_keep(M, keep):
    """partial trace of an operator on n qubits, keeping qubits in `keep` (sorted)."""
    t = M.reshape([2] * (2 * n)); cur = n
    for q in sorted(set(range(n)) - set(keep), reverse=True):
        t = np.trace(t, axis1=q, axis2=q + cur); cur -= 1
    k = len(keep); return t.reshape(2 ** k, 2 ** k)
def embed(Mk, keep):
    """Mk on qubits `keep` (sorted) tensored with identity elsewhere."""
    rest = [q for q in range(n) if q not in keep]
    full = np.kron(Mk, np.eye(2 ** len(rest)))
    order = list(keep) + rest; perm = np.argsort(order)
    t = full.reshape([2] * (2 * n)); t = np.transpose(t, list(perm) + [p + n for p in perm])
    return t.reshape(N, N)
SUBSETS = [tuple(s) for r in range(0, n + 1) for s in it.combinations(range(n), r)]

def support_weights(H):
    """w(S) = squared HS norm of the component of H with exact qubit-support S (S nonempty), via Moebius inversion."""
    wle = {}
    for S in SUBSETS:
        if len(S) == 0:
            P = np.trace(H) / N * np.eye(N)
        else:
            P = embed(ptrace_keep(H, list(S)) / 2 ** (n - len(S)), list(S))
        wle[S] = np.vdot(P, P).real
    w = {}
    for S in SUBSETS:
        if not S: continue
        val = 0.0
        for r in range(len(S) + 1):
            for T in it.combinations(S, r): val += (-1) ** (len(S) - r) * wle[T]
        w[S] = max(val, 0.0)
    return w

def set_partitions(s):
    if not s: yield []; return
    f, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)): yield p[:i] + [[f] + p[i]] + p[i + 1:]
        yield [[f]] + p
GROUPS = [sorted(sorted(b) for b in p) for p in set_partitions(list(range(n))) if len(p) >= 2]
def dtype(g): return "x".join(str(2 ** len(b)) for b in sorted(g, key=len, reverse=True))

def keff(w, g):
    blk = {q: i for i, b in enumerate(g) for q in b}; num = den = 0.0
    for S, v in w.items():
        if v < 1e-12: continue
        k = len({blk[q] for q in S}); num += k * v; den += v
    return num / den

def entropy_bits(r):
    ev = np.linalg.eigvalsh((r + r.conj().T) / 2); ev = ev[ev > 1e-13]; return float(-(ev * np.log2(ev)).sum())
def block_entropies(psi, blocks):
    rho = np.outer(psi, psi.conj()); return [entropy_bits(ptrace_keep(rho, sorted(b))) for b in blocks]

rng0 = np.random.default_rng(777)
def rand_clifford(rng, depth=14):
    U = np.eye(N, dtype=complex)
    for _ in range(depth):
        r = rng.integers(3)
        if r == 0:
            a, b = rng.choice(n, 2, replace=False); U = (op({a: np.diag([1, 0]).astype(complex)}) + op({a: np.diag([0, 1]).astype(complex), b: Z})) @ U
        elif r == 1: U = op({int(rng.integers(n)): Hg}) @ U
        else: U = op({int(rng.integers(n)): Sg}) @ U
    return U
def haar(rng):
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); q, r = np.linalg.qr(A); return q * (np.diag(r) / abs(np.diag(r)))
def product_frame(psi, rng):
    """unitary W with W psi = |0..0> (up to phase), via QR completion: columns of Q start with psi."""
    A = np.c_[psi, rng.normal(size=(N, N - 1)) + 1j * rng.normal(size=(N, N - 1))]
    Q, R = np.linalg.qr(A); Q[:, 0] *= R[0, 0] / abs(R[0, 0]); return Q.conj().T   # W = Q^dag; frame U = W^dag = Q
def lu_equivalent(U1, U2, g):
    W = U1.conj().T @ U2
    for b in g:
        for q in b:
            for P in (X, Z):
                A = W.conj().T @ op({q: P}) @ W
                proj = embed(ptrace_keep(A, sorted(b)) / 2 ** (n - len(b)), sorted(b))
                if np.linalg.norm(A - proj) > 1e-7: return False
    return True

def mfi(): return sum(op({i: Z, i + 1: Z}) for i in range(n - 1)) + 0.9045 * sum(op({i: X}) for i in range(n)) + 0.809 * sum(op({i: Z}) for i in range(n))
def prod_state(a):
    return kron_list([np.array([[np.cos(a / 2)], [np.sin(a / 2)]], complex) for _ in range(n)]).ravel()
TG = np.arange(-10, 10.0001, 0.25); T0 = int(np.argmin(abs(TG)))

def evaluate(H, psi, frames):
    ev, V = np.linalg.eigh(H); rows = []
    for fname, U in frames:                      # frame U: coordinates of a state in this TPS are U^dag psi
        w = support_weights(U.conj().T @ H @ U)
        traj = [U.conj().T @ (V @ (np.exp(-1j * ev * t) * (V.conj().T @ psi))) for t in TG]
        ent_cache = {}
        for g in GROUPS:
            Ct = []
            for k, st in enumerate(traj):
                tot = 0.0
                for b in g:
                    key = (k, tuple(b))
                    if key not in ent_cache: ent_cache[key] = block_entropies(st, [b])[0]
                    tot += ent_cache[key]
                Ct.append(tot)
            Ct = np.array(Ct)
            rows.append(dict(frame=fname, g=g, d=dtype(g), L=-keff(w, g), C0=Ct[T0], Cmin=Ct.min(), tstar=TG[int(np.argmin(Ct))]))
    return rows

def pareto(rows, Ck):
    V = np.round(np.array([[r["L"], -r[Ck]] for r in rows]), 7); fr = []
    for i in range(len(V)):
        if not any((V[j] >= V[i]).all() and (V[j] > V[i]).any() for j in range(len(V))): fr.append(i)
    mx = V.max(0); dom = [i for i in range(len(V)) if (V[i] == mx).all()]
    return fr, dom
def lex(rows, first, Ck):
    key = (lambda r: (round(r["L"], 7), -round(r[Ck], 7))) if first == "L" else (lambda r: (-round(r[Ck], 7), round(r["L"], 7)))
    best = max(key(r) for r in rows); return [i for i, r in enumerate(rows) if key(r) == best]
def lu_classes(idx, rows, frames):
    F = dict(frames); classes = []
    for i in idx:
        for cl in classes:
            j = cl[0]
            if rows[j]["g"] == rows[i]["g"] and lu_equivalent(F[rows[j]["frame"]], F[rows[i]["frame"]], rows[i]["g"]): cl.append(i); break
        else: classes.append([i])
    return classes
def lab(r): return f"{r['frame']}:{r['d']}:{''.join('(' + ''.join(map(str, b)) + ')' for b in r['g'])}"

if __name__ == "__main__":
    H = mfi(); p = prod_state(1.4)
    print("=" * 78); print("INDEPENDENT REPRODUCTION NR-SigmaH-1..5 (3 fresh frame seeds)"); print("=" * 78)
    # sanity: Moebius support weights reproduce ||H - tr part||^2
    w = support_weights(H); tl = H - np.trace(H) / N * np.eye(N)
    print(f"sanity: sum of exact-support weights {sum(w.values()):.6f} vs ||H - tr||^2 {np.vdot(tl, tl).real:.6f}; "
          f"max support size with weight: {max(len(S) for S, v in w.items() if v > 1e-10)} (expect 2)")
    summary = {k: [] for k in ("NR1", "NR2", "NR3", "NR4", "NR5")}
    for seed in (11, 22, 33):
        rng = np.random.default_rng(seed)
        C = [rand_clifford(rng) for _ in range(3)]; Hu = haar(rng)
        haar_psi = (lambda v: v / np.linalg.norm(v))(rng.normal(size=N) + 1j * rng.normal(size=N))
        def frames_for(Hm, psi):
            return [("id", np.eye(N, dtype=complex)), ("cl1", C[0]), ("cl2", C[1]), ("cl3", C[2]), ("haar", Hu),
                    ("comm", expm(-1j * Hm * 0.7)), ("Wpsi", product_frame(psi, rng).conj().T)]
        print(f"\n--- frame seed {seed} ---")
        # NR1 / NR4 positive compatible case
        fr1 = frames_for(H, p); R1 = evaluate(H, p, fr1); front, dom = pareto(R1, "Cmin"); cls = lu_classes(front, R1, fr1)
        frames_on = sorted({R1[c[0]]["frame"] for c in cls}); dts = sorted({R1[c[0]]["d"] for c in cls})
        la, lb = lex(R1, "L", "Cmin"), lex(R1, "C", "Cmin")
        ok1 = bool(dom) and frames_on == ["id"] and {R1[i]["frame"] for i in la} <= {R1[c[0]]["frame"] for c in cls} | {R1[i]["frame"] for i in front}
        print(f"NR1 positive: front {len(front)} -> {len(cls)} LU classes; frames {frames_on}; dominant {'yes' if dom else 'NO'}; "
              f"rule A frames {sorted({R1[i]['frame'] for i in la})}, rule B frames {sorted({R1[i]['frame'] for i in lb})}")
        summary["NR1"].append(bool(dom) and frames_on == ["id"])
        tie = len(dts) > 1 and len({round(R1[c[0]]['L'], 7) for c in cls}) == 1 and len({round(R1[c[0]]['Cmin'], 7) for c in cls}) == 1
        print(f"NR4 d-tie: dominant LU classes span dtypes {dts}; identical (L, Cmin) vectors: {tie}")
        summary["NR4"].append(tie)
        # NR2 conflict
        psi_c = C[0] @ p; fr2 = frames_for(H, psi_c); R2 = evaluate(H, psi_c, fr2); front2, dom2 = pareto(R2, "Cmin")
        cls2 = lu_classes(front2, R2, fr2); f2 = sorted({R2[c[0]]['frame'] for c in cls2})
        la2, lb2 = lex(R2, "L", "Cmin"), lex(R2, "C", "Cmin")
        A_fr = sorted({R2[i]['frame'] for i in la2}); B_fr = sorted({R2[i]['frame'] for i in lb2})
        ok2 = (not dom2) and "id" in f2 and "cl1" in f2 and A_fr != B_fr
        print(f"NR2 conflict: front {len(front2)} -> {len(cls2)} LU classes; frames {f2}; dominant {'yes' if dom2 else 'NONE'}; rule A {A_fr} vs rule B {B_fr}")
        summary["NR2"].append(ok2)
        # NR3 Haar
        fr3 = frames_for(H, haar_psi); R3 = evaluate(H, haar_psi, fr3); front3, dom3 = pareto(R3, "Cmin")
        lb3 = lex(R3, "C", "Cmin"); B3 = sorted({R3[i]['frame'] for i in lb3}); la3 = lex(R3, "L", "Cmin"); A3 = sorted({R3[i]['frame'] for i in la3})
        ok3 = (not dom3) and B3 == ["Wpsi"]
        print(f"NR3 Haar: front {len(front3)}; dominant {'yes' if dom3 else 'NONE'}; rule B (correlation first) frames {B3}; rule A frames {A3}")
        summary["NR3"].append(ok3)
        # NR5 epoch covariance (positive case), s = 2.25
        ev, V = np.linalg.eigh(H); ps = V @ (np.exp(-1j * ev * 2.25) * (V.conj().T @ p))
        fr5 = frames_for(H, ps); fr5 = [f for f in fr5 if f[0] != "Wpsi"] + [("Wpsi", fr1[-1][1])]   # keep W frame identical to NR1 for a like-for-like comparison
        R5 = evaluate(H, ps, fr5)
        _, dom0_before = pareto(R1, "C0"); _, dom0_after = pareto(R5, "C0"); _, domm_after = pareto(R5, "Cmin")
        same = {lab(R1[i]) for i in dom} == {lab(R5[i]) for i in domm_after}
        tb = sorted({R1[i]['tstar'] for i in dom}); ta = sorted({R5[i]['tstar'] for i in domm_after})
        print(f"NR5 epoch: C0 dominance before {'yes' if dom0_before else 'NONE'} / after shift {'yes' if dom0_after else 'NONE'};"
              f" Cmin dominant set identical after shift: {same}; t* before {tb} after {ta}")
        summary["NR5"].append((not dom0_after) and same and ta == [t - 2.25 for t in tb])
    print("\n" + "=" * 78)
    claims = {"NR1": "positive compatible case: dominant set exists and lies in the native (id) frame after LU quotient",
              "NR2": "conflict: no dominance; both the H-local and the psi-product frames on the front; rules A and B disagree",
              "NR3": "Haar: no dominance; correlation-first rule returns the trivial psi-product frame",
              "NR4": "local-dimension tie: dominant classes span several dtypes with identical (L, Cmin)",
              "NR5": "epoch: C0 loses dominance under time shift; Cmin dominant set unchanged with t* shifted by -2.25"}
    for k, v in summary.items():
        print(f"{k}: {'REPRODUCED' if all(v) else 'NOT REPRODUCED' if not any(v) else 'PARTIAL'} ({sum(v)}/3 seeds) — {claims[k]}")
