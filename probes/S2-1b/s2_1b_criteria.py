"""SCOUT-2 S2-1b: can locality be selected without supplying graph distance, local dimension or k?

The SAME abstract dynamics H (an N x N Hermitian matrix) is scored in several candidate tensor-product
structures (TPS = a unitary frame U plus a grouping of qubit slots into factors). Seven criteria:

  1 minimal k        : max support size (in factors) of the exact-support decomposition of U^dag H U
  2 sparse graph     : # factor pairs that interact (contained in some support)
  3 Lieb-Robinson    : mean over factor pairs of normalized ||[P_a(t0), Q_b]||^2 (smaller = sharper cone)
  4 MDL              : dim of the support-closed operator space containing H (parameters needed)
  5 stability        : does the winner survive a small perturbation? (perturbation measure must be named)
  6 locality-max     : fraction of ||H||^2 on supports of size <= 2 factors
  7 autonomy         : fraction of ||H||^2 on single-factor supports (factors evolve ~ autonomously)

The exact-support decomposition (component of H on strings of exact factor-support S) is LU-invariant, so
every criterion is a TPS invariant (not a basis artefact).
"""
import itertools as it
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(31)
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
PA = [I2, X, Y, Z]


def kron(*ms):
    out = np.eye(1)
    for m in ms:
        out = np.kron(out, m)
    return out


def pauli_basis(n):
    labels = list(it.product(range(4), repeat=n))
    mats = np.array([kron(*[PA[a] for a in l]) for l in labels])
    return labels, mats


BASES = {}


def decompose(H, n):
    if n not in BASES:
        BASES[n] = pauli_basis(n)
    labels, mats = BASES[n]
    c = np.einsum('kij,ji->k', mats, H).real / 2 ** n
    return labels, c


def op(n, site_ops):
    return kron(*[site_ops.get(q, I2) for q in range(n)])


class TPS:
    def __init__(self, name, n, groups, U=None):
        self.name, self.n, self.groups = name, n, groups
        self.U = np.eye(2 ** n) if U is None else U
        self.g_of = {q: gi for gi, g in enumerate(groups) for q in g}

    def fsupp(self, label):
        return frozenset(self.g_of[q] for q, a in enumerate(label) if a)


def scores(H, T, tol=1e-9, t0=0.3, lr=True):
    n = T.n
    Hl = T.U.conj().T @ H @ T.U
    labels, c = decompose(Hl, n)
    w = {}
    for l, ci in zip(labels, c):
        s = T.fsupp(l)
        if s and abs(ci) > tol:
            w[s] = w.get(s, 0) + ci ** 2
    tot = sum(w.values())
    supports = list(w)
    k = max(len(s) for s in supports)
    pairs = {frozenset(p) for s in supports for p in it.combinations(sorted(s), 2)}
    # MDL: number of Pauli strings whose factor-support is non-empty and inside some support
    mdl = sum(1 for l in labels if (s := T.fsupp(l)) and any(s <= S for S in supports))
    loc2 = sum(v for s, v in w.items() if len(s) <= 2) / tot
    auto = sum(v for s, v in w.items() if len(s) == 1) / tot
    out = dict(k=k, edges=len(pairs), npairs=len(T.groups) * (len(T.groups) - 1) // 2, mdl=mdl, loc2=loc2, auto=auto)
    if lr:
        Ut = expm(-1j * Hl * t0)
        vals = []
        for a, b in it.combinations(range(len(T.groups)), 2):
            Pa = [l for l in it.product(range(4), repeat=len(T.groups[a])) if any(l)]
            Qb = [l for l in it.product(range(4), repeat=len(T.groups[b])) if any(l)]
            acc = []
            for la in Pa:
                P = op(n, {q: PA[x] for q, x in zip(T.groups[a], la)})
                Pt = Ut.conj().T @ P @ Ut
                for lb in Qb:
                    Q = op(n, {q: PA[x] for q, x in zip(T.groups[b], lb)})
                    C = Pt @ Q - Q @ Pt
                    acc.append(np.vdot(C, C).real / (4 * 2 ** n))
            vals.append(np.mean(acc))
        out['lr'] = float(np.mean(vals))
    return out


def rnd_herm(d, s=1.0):
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return s * (A + A.conj().T) / 2


def table(H, tpss, **kw):
    res = {T.name: scores(H, T, **kw) for T in tpss}
    keys = ['k', 'edges', 'mdl', 'lr', 'loc2', 'auto']
    print(f"  {'TPS':<26}" + "".join(f"{k:>10}" for k in keys))
    for nm, r in res.items():
        print(f"  {nm:<26}" + "".join(f"{(r[k] if not isinstance(r.get(k), float) else round(r[k], 4))!s:>10}" for k in keys))
    best = {
        'minimal k': min(res, key=lambda t: res[t]['k']),
        'sparse graph (edge fraction)': min(res, key=lambda t: res[t]['edges'] / max(res[t]['npairs'], 1)),
        'Lieb-Robinson': min(res, key=lambda t: res[t]['lr']),
        'MDL': min(res, key=lambda t: res[t]['mdl']),
        'locality-max': max(res, key=lambda t: res[t]['loc2']),
        'autonomy': max(res, key=lambda t: res[t]['auto']),
    }
    for c, t in best.items():
        ties = [u for u in res if u != t and {
            'minimal k': lambda u: res[u]['k'] == res[t]['k'],
            'sparse graph (edge fraction)': lambda u: abs(res[u]['edges'] / max(res[u]['npairs'], 1) - res[t]['edges'] / max(res[t]['npairs'], 1)) < 1e-12,
            'Lieb-Robinson': lambda u: abs(res[u]['lr'] - res[t]['lr']) < 1e-9,
            'MDL': lambda u: res[u]['mdl'] == res[t]['mdl'],
            'locality-max': lambda u: abs(res[u]['loc2'] - res[t]['loc2']) < 1e-9,
            'autonomy': lambda u: abs(res[u]['auto'] - res[t]['auto']) < 1e-9}[c](u)]
        print(f"    {c:<30} -> {t}" + (f"   (TIE with {', '.join(ties)})" if ties else ""))
    return res, best


print("=" * 78)
print("S2-1b  competing locality criteria on the same abstract dynamics")
print("=" * 78)

# ---------------- 0. no local dimension, no number of factors: degeneracy ----------------
print("\n--- 0. Minimal k without a supplied number of factors (N = 16) ---")
H0 = rnd_herm(16)
for name, groups in [("1 factor (16)", [[0, 1, 2, 3]]), ("4 (x) 4", [[0, 1], [2, 3]]), ("2 (x) 8", [[0], [1, 2, 3]]),
                     ("2 (x) 2 (x) 4", [[0], [1], [2, 3]]), ("2^4 qubits", [[0], [1], [2], [3]])]:
    r = scores(H0, TPS(name, 4, groups), lr=False)
    print(f"  random H, {name:<16}: k = {r['k']}  (#factors {len(groups)})")
print("  => over all factorizations minimal k is won by the TRIVIAL factorization (k = 1). The criterion is")
print("     ill-posed unless the number of factors or the local dimension is supplied (finest = prime: 2^4).")

# ---------------- Case A: qubits vs ququarts (local dimension conflict) ----------------
print("\n--- Case A. N = 16: strong generic pair terms + weak bridge.  H = h12 + h34 + g X2 X3 ---")
g = 0.05
h12 = rnd_herm(4); h34 = rnd_herm(4)
HA = np.kron(h12, np.eye(4)) + np.kron(np.eye(4), h34) + g * op(4, {1: X, 2: X})
tq = TPS("qubits 2(x)2(x)2(x)2", 4, [[0], [1], [2], [3]])
tQ = TPS("ququarts 4(x)4", 4, [[0, 1], [2, 3]])
t28 = TPS("2 (x) 8", 4, [[0], [1, 2, 3]])
resA, bestA = table(HA, [tq, tQ, t28])

# ---------------- Case B: two qubit TPSs related by a Clifford frame ----------------
print("\n--- Case B. n = 5 qubits, same local dimension, same #factors; TPS2 = Clifford-rotated frame ---")
n = 5
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2); Sg = np.diag([1, 1j])
def cnot(c, t):
    P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
    return op(n, {c: P0}) + op(n, {c: P1, t: X})
J = rng.uniform(0.5, 1.5, n - 1); hx = rng.uniform(0.5, 1.5, n)
HB = sum(J[i] * op(n, {i: Z, i + 1: Z}) for i in range(n - 1)) + sum(hx[i] * op(n, {i: X}) for i in range(n))
t1 = TPS("computational", n, [[q] for q in range(n)])
HB3 = sum(op(n, {i: Z, i + 1: Z, i + 2: Z}) * rng.uniform(0.5, 1.5) for i in range(n - 2)) + sum(hx[i] * op(n, {i: X}) for i in range(n))
for hname, Hb in [("transverse-field Ising chain (2-local NN)", HB), ("3-body ZZZ chain + X fields", HB3)]:
    found = {}
    for trial in range(1500):
        U = np.eye(2 ** n, dtype=complex)
        for _ in range(rng.integers(1, 9)):
            r_ = rng.integers(3)
            if r_ == 0:
                c, t = rng.choice(n, 2, replace=False); U = cnot(c, t) @ U
            elif r_ == 1:
                U = op(n, {int(rng.integers(n)): Hd}) @ U
            else:
                U = op(n, {int(rng.integers(n)): Sg}) @ U
        s1 = scores(Hb, t1, lr=False); s2 = scores(Hb, TPS("rot", n, [[q] for q in range(n)], U), lr=False)
        better = {c for c in ['k', 'edges', 'mdl'] if s2[c] < s1[c]}
        worse = {c for c in ['k', 'edges', 'mdl'] if s2[c] > s1[c]}
        if better and worse:
            key = (tuple(sorted(better)), tuple(sorted(worse)))
            if key not in found:
                found[key] = U
    print(f"\n  H = {hname}: 1500 random Clifford frames; strict exact-criterion disagreement patterns: {len(found)}")
    for (b, w_), U in list(found.items())[:2]:
        print(f"\n  pattern: rotated frame BETTER on {b}, WORSE on {w_}")
        table(Hb, [t1, TPS("Clifford-rotated", n, [[q] for q in range(n)], U)])

# ---------------- Stability: the perturbation measure decides ----------------
print("\n--- Stability (Case A). Perturb H -> H + eps V. Which measure for V? ---")
for label, Vgen in [("GUE (TPS-neutral)", lambda: rnd_herm(16)),
                    ("2-local in QUBIT TPS", lambda: sum(rng.normal() * op(4, {i: PA[a], j: PA[b]}) for i, j in it.combinations(range(4), 2) for a in range(1, 4) for b in range(1, 4))),
                    ("local in QUQUART TPS", lambda: np.kron(rnd_herm(4), np.eye(4)) + np.kron(np.eye(4), rnd_herm(4)))]:
    for eps in [1e-3, 1e-1]:
        V = Vgen(); V = V / np.linalg.norm(V) * np.linalg.norm(HA)
        Hp = HA + eps * V
        sq = scores(Hp, tq, lr=False); sQ = scores(Hp, tQ, lr=False)
        print(f"  V ~ {label:<22} eps={eps:<6}: k q/Q = {sq['k']}/{sQ['k']}, MDL q/Q = {sq['mdl']}/{sQ['mdl']}, "
              f"loc2 q/Q = {sq['loc2']:.3f}/{sQ['loc2']:.3f}, auto q/Q = {sq['auto']:.3f}/{sQ['auto']:.3f}")
print("  => exact criteria (k, edges, MDL) are destroyed by any TPS-neutral perturbation; they 'survive' only under a")
print("     perturbation measure that is itself local in a chosen TPS (circular). Weighted criteria survive but still")
print("     disagree (loc2 vs autonomy). STABILITY IS MEASURE-PRICED.")
