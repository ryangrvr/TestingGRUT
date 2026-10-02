"""SCOUT-2 S2-Sigma (+ G1): factorization AND local dimension selection without hidden weights.

Pre-registered (probes/PROBE_CHARTERS.md, Wave 2).  n = 6 qubit slots, N = 64.
Candidates = (set partition of the 6 slots into >= 2 blocks: 202 groupings, local dims 2..32, mixed)
           x (frames: identity, 4 random Clifford circuits, 1 Haar unitary, commutant frame W = exp(-iHs)).
Criteria (larger = better, each LU-invariant):
  L = -k_eff      (||H||^2-weighted mean factor-support size)
  Q = -mean_f [S_f(tau) - S_f(0)]           (Carroll-Singh-type entanglement-growth penalty)
  R = max_f #{g != f : I(f:g)(tau) >= (1-delta) S_f(tau)}, only if S_f(tau) > 0.1 bit (record redundancy)
  P = ||H||^2 fraction on single-factor supports (predictive autonomy)
  M = -dim(support-closed operator space containing H) / N^2 (description length)
"""
import itertools as it
import sys
import numpy as np
from scipy.linalg import expm

n = 6; N = 2 ** n
rng = np.random.default_rng(2026)
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
PA = np.array([I2, X, Y, Z])
NAMES = ['L', 'Q', 'R', 'P', 'M']


def op(d):
    out = np.eye(1)
    for q in range(n):
        out = np.kron(out, d.get(q, I2))
    return out


def pauli_coeffs(H):
    """c[a1..an] = tr(P_a H)/N via sequential contraction."""
    t = H.reshape([2] * (2 * n))
    # order axes as (i1, j1, i2, j2, ...)
    t = np.transpose(t, [x for q in range(n) for x in (q, q + n)])
    for q in range(n):
        # contract pair (i, j) at positions q, q+1 with P[a, j, i]
        t = np.tensordot(t, PA, axes=([q, q + 1], [2, 1]))   # appends axis a at end
        t = np.moveaxis(t, -1, q)
    return (t.reshape(-1).real / N)


LABELS = np.array(list(it.product(range(4), repeat=n)))
NONID = LABELS != 0                     # (4096, 6)


def set_partitions(s):
    if not s:
        yield []
        return
    first, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


GROUPINGS = [sorted(sorted(b) for b in p) for p in set_partitions(list(range(n))) if len(p) >= 2]
assert len(GROUPINGS) == 202


def dtype_of(g):
    return "x".join(str(2 ** len(b)) for b in sorted(g, key=len, reverse=True))


def h_criteria(c, g):
    w = c ** 2; w[0] = 0.0
    T = np.stack([NONID[:, b].any(1) for b in g], 1)
    size = T.sum(1)
    code = T.astype(np.int64) @ (1 << np.arange(len(g)))
    tot = w.sum()
    L = -(w * size).sum() / tot
    P = w[size == 1].sum() / tot
    present = np.unique(code[(np.abs(c) > 1e-9) & (code != 0)])
    covered = np.zeros(len(code), bool)
    for S in present:
        covered |= (code & ~S) == 0
    covered &= code != 0
    M = -covered.sum() / N ** 2
    return L, P, M


SUBSETS = [tuple(s) for r in range(1, n + 1) for s in it.combinations(range(n), r)]


def entropies(state):
    """S(subset) for all 63 subsets; state = pure vector or density matrix."""
    out = {}
    if state.ndim == 1:
        t = state.reshape([2] * n)
        for s in SUBSETS:
            rest = [q for q in range(n) if q not in s]
            m = np.transpose(t, list(s) + rest).reshape(2 ** len(s), -1)
            p = np.linalg.svd(m, compute_uv=False) ** 2
            p = p[p > 1e-14]; out[s] = float(-(p * np.log2(p)).sum())
    else:
        t = state.reshape([2] * (2 * n))
        for s in SUBSETS:
            rest = [q for q in range(n) if q not in s]
            m = np.transpose(t, list(s) + rest + [x + n for x in s] + [x + n for x in rest])
            m = m.reshape(2 ** len(s), 2 ** len(rest), 2 ** len(s), 2 ** len(rest))
            r = np.einsum('ajbj->ab', m)
            p = np.linalg.eigvalsh(r); p = p[p > 1e-14]
            out[s] = float(-(p * np.log2(p)).sum())
    return out


def state_criteria(E0, Et, g, delta=0.1, frag='single'):
    blocks = [tuple(b) for b in g]
    Q = -np.mean([Et[b] - E0[b] for b in blocks])
    R = 0
    for f in blocks:
        Sf = Et[f]
        if Sf <= 0.1:
            continue
        others = [b for b in blocks if b != f]
        frags = others if frag == 'single' else [tuple(sorted(a + b)) for a, b in it.combinations(others, 2)]
        cnt = sum(1 for h in frags if Sf + Et[h] - Et[tuple(sorted(f + h))] >= (1 - delta) * Sf - 1e-12)
        R = max(R, cnt)
    return Q, R


# ---------------- frames ----------------
Hd = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2); Sg = np.diag([1, 1j])
def cnot(c, t):
    return op({c: np.diag([1, 0]).astype(complex)}) + op({c: np.diag([0, 1]).astype(complex), t: X})
def rand_clifford(depth):
    U = np.eye(N, dtype=complex)
    for _ in range(depth):
        r = rng.integers(3)
        if r == 0:
            a, b = rng.choice(n, 2, replace=False); U = cnot(a, b) @ U
        elif r == 1:
            U = op({int(rng.integers(n)): Hd}) @ U
        else:
            U = op({int(rng.integers(n)): Sg}) @ U
    return U
def haar_u():
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    q, r = np.linalg.qr(A); return q * (np.diag(r) / np.abs(np.diag(r)))
CLIFFS = [rand_clifford(12) for _ in range(4)]
HAAR = haar_u()


def frames_for(H, s=0.7):
    return [("id", np.eye(N)), ("cliff1", CLIFFS[0]), ("cliff2", CLIFFS[1]), ("cliff3", CLIFFS[2]),
            ("cliff4", CLIFFS[3]), ("haar", HAAR), ("commutant e^-iHs", expm(-1j * H * s))]


# ---------------- models ----------------
def rh(d, s=1.0):
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); return s * (A + A.conj().T) / 2
def embed_pair(h4, a):
    # h4 acts on qubits (a, a+1)
    return np.kron(np.kron(np.eye(2 ** a), h4), np.eye(2 ** (n - a - 2)))
MODELS = {}
MODELS['M1 clustered pairs + weak bridges'] = sum(embed_pair(rh(4), a) for a in (0, 2, 4)) + 0.05 * (op({1: X, 2: X}) + op({3: X, 4: X}))
MODELS['M2 TFIM ring (translation-symmetric)'] = sum(op({i: Z, (i + 1) % n: Z}) for i in range(n)) + sum(op({i: X}) for i in range(n))
gk = rng.uniform(0.8, 1.2, n - 1)
MODELS['M3 record-forming star'] = sum(gk[k - 1] * op({0: Z, k: Y}) for k in range(1, n)) + sum(0.1 * rng.uniform(.5, 1.5) * op({q: Z}) for q in range(1, n))
M4 = sum(rng.normal() * op({i: PA[a], i + 1: PA[b]}) for i in range(n - 1) for a in range(1, 4) for b in range(1, 4))
MODELS['M4 generic random 2-local chain'] = M4 + sum(rng.normal() * op({i: PA[a]}) for i in range(n) for a in range(1, 4))


def states_for(H):
    ev, V = np.linalg.eigh(H)
    prod = np.zeros(N, complex); prod[0] = 1
    lp = np.array([1.0 + 0j])
    for q in range(n):
        v = rng.normal(size=2) + 1j * rng.normal(size=2); lp = np.kron(lp, v / np.linalg.norm(v))
    hv = rng.normal(size=N) + 1j * rng.normal(size=N); hv /= np.linalg.norm(hv)
    rec = op({0: Hd}) @ prod
    gibbs = (V * np.exp(-1.0 * (ev - ev.min()))) @ V.conj().T; gibbs /= np.trace(gibbs)
    return {"product |0..0>": prod, "random local product": lp, "ground": V[:, 0], "mid-spectrum eigenstate": V[:, N // 2],
            "Haar": hv, "record-forming |+>|0..>": rec, "Gibbs beta=1": gibbs}


def evolve(state, H, tau):
    U = expm(-1j * H * tau)
    return U @ state if state.ndim == 1 else U @ state @ U.conj().T


def in_frame(state, U):
    return U.conj().T @ state if state.ndim == 1 else U.conj().T @ state @ U


def vectors(H, state, tau=1.0, delta=0.1, frag='single', frames=None):
    frames = frames or frames_for(H)
    st_t = evolve(state, H, tau)
    rows, meta = [], []
    for fname, U in frames:
        c = pauli_coeffs(U.conj().T @ H @ U)
        E0 = entropies(in_frame(state, U)); Et = entropies(in_frame(st_t, U))
        for g in GROUPINGS:
            L, P, M = h_criteria(c, g)
            Q, R = state_criteria(E0, Et, g, delta, frag)
            rows.append([L, Q, R, P, M]); meta.append((fname, dtype_of(g), g))
    return np.round(np.array(rows), 9), meta


def pareto(V):
    ge = (V[None, :, :] >= V[:, None, :]).all(2)      # ge[i,j]: j >= i everywhere
    gt = (V[None, :, :] > V[:, None, :]).any(2)
    dominated = (ge & gt).any(1)
    return np.where(~dominated)[0]


def dominant(V):
    mx = V.max(0)
    return np.where((V == mx).all(1))[0]


def normalize(V, how):
    if how == 'minmax':
        span = V.max(0) - V.min(0); span[span == 0] = 1; return (V - V.min(0)) / span
    if how == 'zscore':
        sd = V.std(0); sd[sd == 0] = 1; return (V - V.mean(0)) / sd
    if how == 'rank':
        return np.argsort(np.argsort(V, 0), 0).astype(float) / len(V)


TRANSFORMS = {'identity': lambda x: x, 'log1p': lambda x: np.log1p(x - x.min(0)), 'square': lambda x: (x - x.min(0)) ** 2,
              'sqrt': lambda x: np.sqrt(x - x.min(0))}


def scalar_winners(V):
    wins = set()
    for tn, T in TRANSFORMS.items():
        TV = T(V)
        for nm in ['minmax', 'zscore', 'rank']:
            wins.add(int(np.argmax(normalize(TV, nm).sum(1))))
    return wins


def simplex(V, how='minmax', samples=20000):
    W = rng.dirichlet(np.ones(V.shape[1]), samples)
    win = np.argmax(normalize(V, how) @ W.T, 0)
    u, cnt = np.unique(win, return_counts=True)
    return sorted(zip(cnt / samples, u), reverse=True)


def label(meta, i):
    f, d, g = meta[i]
    return f"{f}:{d}:{''.join('(' + ''.join(map(str, b)) + ')' for b in g)}"


def report(V, meta, tag, short=False):
    front = pareto(V); dom = dominant(V)
    dtypes = sorted({meta[i][1] for i in front}); fr = sorted({meta[i][0] for i in front})
    print(f"\n[{tag}]  candidates {len(V)};  Pareto front size {len(front)};  weakly dominant candidate: "
          f"{'YES ' + label(meta, dom[0]) if len(dom) else 'NONE'}")
    print(f"   front local-dim types ({len(dtypes)}): {', '.join(dtypes[:12])}{' ...' if len(dtypes) > 12 else ''}")
    print(f"   front frames: {fr}")
    for k, nm in enumerate(NAMES):
        best = np.where(V[:, k] == V[:, k].max())[0]
        print(f"   best {nm}: {V[best[0], k]:+.4f}  by {len(best)} candidates, e.g. {label(meta, best[0])}")
    sw = scalar_winners(V)
    print(f"   Sigma-2 equal-weight scalar winners over 4 transforms x 3 normalizations: {len(sw)} distinct: "
          + "; ".join(label(meta, i) for i in list(sw)[:4]))
    for how in ['minmax', 'rank']:
        sp = simplex(V, how)
        print(f"   Sigma-3 Dirichlet(1) simplex [{how}]: {len(sp)} distinct winners; top shares: "
              + "; ".join(f"{s:.2f} {label(meta, i)}" for s, i in sp[:3]))
    return front, dom


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    print("=" * 78); print("S2-Sigma + G1: factorization and local-dimension selection without hidden weights"); print("=" * 78)
    print(f"groupings {len(GROUPINGS)}; frames 7; candidates per run {len(GROUPINGS) * 7}")
    print("NOTE (Sigma-2, exact): Pareto dominance is invariant under every strictly increasing componentwise transform;")
    print("      the front below is therefore objective-scale-invariant by construction. Scalar winners are not.")

    if part in ('all', 'main'):
        print("\n" + "#" * 78 + "\n# Sigma-1/2/3/4/6: models x states (tau = 1, delta = 0.1, eps = 0)\n" + "#" * 78)
        for mname, H in MODELS.items():
            print("\n" + "=" * 30 + f" {mname} " + "=" * 30)
            for sname, st in states_for(H).items():
                V, meta = vectors(H, st)
                report(V, meta, f"{mname} | {sname}")

    if part in ('all', 'scale'):
        print("\n" + "#" * 78 + "\n# Sigma-5: scale hostile (tau, delta, fragment, perturbation eps)\n" + "#" * 78)
        for mname in ['M1 clustered pairs + weak bridges', 'M4 generic random 2-local chain']:
            H0 = MODELS[mname]
            G = rh(N); G = G / np.linalg.norm(G) * np.linalg.norm(H0)
            st = states_for(H0)["record-forming |+>|0..>"]
            for eps in [0.0, 0.05, 0.3]:
                H = H0 + eps * G
                for tau, delta, frag in [(0.2, 0.1, 'single'), (1.0, 0.1, 'single'), (5.0, 0.1, 'single'),
                                         (1.0, 0.3, 'single'), (1.0, 0.1, 'pairs')]:
                    V, meta = vectors(H, st, tau, delta, frag)
                    front = pareto(V); dom = dominant(V)
                    sp = simplex(V, 'minmax', 8000)
                    sw = scalar_winners(V)
                    print(f"  {mname[:3]} eps={eps:<4} tau={tau:<3} delta={delta} frag={frag:<6}: front {len(front):4d}, "
                          f"dominant {'YES' if len(dom) else 'no '}, front dtypes {len({meta[i][1] for i in front}):2d}, "
                          f"scalar winners {len(sw)}, simplex top {sp[0][0]:.2f} {label(meta, sp[0][1])}")

    if part in ('all', 'sym'):
        print("\n" + "#" * 78 + "\n# Sigma-7: symmetry / degeneracy\n" + "#" * 78)
        H = MODELS['M2 TFIM ring (translation-symmetric)']
        st = states_for(H)["product |0..0>"]
        V, meta = vectors(H, st)
        idx = {(m[0], str(m[2])): i for i, m in enumerate(meta)}
        a = idx[("id", str([[0, 1], [2, 3], [4, 5]]))]; b = idx[("id", str([[0, 5], [1, 2], [3, 4]]))]
        print(f"  (a) ring, translation-related pairings (01)(23)(45) vs (12)(34)(05): V_a = {V[a]}, V_b = {V[b]}, "
              f"identical: {np.allclose(V[a], V[b])}")
        groups = {}
        for i in range(len(V)):
            groups.setdefault(tuple(V[i]), []).append(i)
        front = pareto(V)
        degen = [len(groups[tuple(V[i])]) for i in front]
        print(f"      front size {len(front)}; front members sharing an identical vector with another candidate: "
              f"{sum(1 for d in degen if d > 1)}")
        # (b) commutant frame
        for mname, H in MODELS.items():
            st = states_for(H)["random local product"]
            s = 0.7; W = expm(-1j * H * s)
            c_id = pauli_coeffs(H); c_W = pauli_coeffs(W.conj().T @ H @ W)
            dh = max(np.abs(np.array(h_criteria(c_id, g)) - np.array(h_criteria(c_W, g))).max() for g in GROUPINGS)
            # psi-criteria in W frame with psi  ==  psi-criteria in id frame with psi(-s)
            Vw, _ = vectors(H, st, frames=[("W", W)])
            Vs, _ = vectors(H, expm(1j * H * s) @ st, frames=[("id", np.eye(N))])
            Vi, _ = vectors(H, st, frames=[("id", np.eye(N))])
            print(f"  (b) {mname}: max |H-only criteria (L,P,M): id frame - commutant frame| = {dh:.1e};  "
                  f"max |V(W frame, psi) - V(id frame, psi(-s))| = {np.abs(Vw - Vs).max():.1e};  "
                  f"max |V(W frame, psi) - V(id frame, psi)| = {np.abs(Vw - Vi).max():.3f}")
