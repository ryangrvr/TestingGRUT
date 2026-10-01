"""SCOUT-2 S2-8: quantum Darwinism / records.  Does redundancy select the pointer observable AND the
system/environment split, or does it presuppose the split?  (Pre-registered firewall: if redundancy is defined
only after the split S | E1 (x) ... (x) Em is specified -> A-PRICED.)

Model: qubit S (slot 0) + m environment qubits.  H = sum_k g_k Z_S Y_k  (S's Z is copied into each E_k).
S starts in |+>, E in |0...0>.
"""
import itertools as it
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(8)
m = 8; n = m + 1
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
Hd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


def op(site_ops):
    out = np.eye(1)
    for q in range(n):
        out = np.kron(out, site_ops.get(q, I2))
    return out


def rho_of(psi, keep):
    t = psi.reshape([2] * n)
    rest = [q for q in range(n) if q not in keep]
    t = np.transpose(t, list(keep) + rest).reshape(2 ** len(keep), -1)
    return t @ t.conj().T


def S_vn(r):
    w = np.linalg.eigvalsh(r); w = w[w > 1e-12]
    return float(-(w * np.log2(w)).sum())


def MI(psi, A, B):
    return S_vn(rho_of(psi, A)) + S_vn(rho_of(psi, B)) - S_vn(rho_of(psi, sorted(A + B)))


def plateau(psi, sys, env, samples=12):
    """average I(sys : F) over random fragments F of each size"""
    HS = S_vn(rho_of(psi, sys))
    out = []
    for f in range(0, len(env) + 1):
        combs = list(it.combinations(env, f))
        if len(combs) > samples:
            combs = [combs[i] for i in rng.choice(len(combs), samples, replace=False)]
        out.append(np.mean([MI(psi, sys, list(c)) if f else 0.0 for c in combs]))
    return HS, out


def redundancy(HS, curve, nfrag, delta=0.1):
    for f, v in enumerate(curve):
        if f and v >= (1 - delta) * HS:
            return nfrag / f
    return 0.0


def evolve(g, t, extra=None):
    H = sum(g[k] * op({0: Z, k + 1: Y}) for k in range(m))
    if extra is not None:
        H = H + extra
    psi0 = np.zeros(2 ** n, complex); psi0[0] = 1
    psi0 = op({0: Hd}) @ psi0
    return expm(-1j * H * t) @ psi0


def holevo(psi, basis, frag):
    """chi(basis on S : fragment)"""
    cond = []; ps = []
    for v in basis:
        P = op({0: np.outer(v, v.conj())})
        phi = P @ psi; p = np.vdot(phi, phi).real
        if p < 1e-12:
            continue
        phi = phi / np.sqrt(p)
        cond.append(rho_of(phi, frag)); ps.append(p)
    avg = sum(p * r for p, r in zip(ps, cond))
    return S_vn(avg) - sum(p * S_vn(r) for p, r in zip(ps, cond))


print("=" * 78); print(f"S2-8  quantum Darwinism: S + {m} environment qubits, H = sum g_k Z_S Y_k"); print("=" * 78)
g = rng.uniform(0.8, 1.2, m)
t = np.pi / 4 / g.mean() * 0.8          # partial records (overlap cos(2 g t) != 0)
psi = evolve(g, t)
sysA = [0]; envA = list(range(1, n))

print("\n--- 1. GIVEN the split S | E1..E8 (qubits): mutual-information plateau ---")
HS, cur = plateau(psi, sysA, envA)
print("  |F| : " + " ".join(f"{f:5d}" for f in range(m + 1)))
print("  I   : " + " ".join(f"{v:5.3f}" for v in cur) + f"    H(S) = {HS:.3f}")
print(f"  redundancy R_0.1 (fragments = single qubits) = {redundancy(HS, cur, m):.2f}")

print("\n--- 2. Pointer selection GIVEN the split: Holevo chi(basis on S : one env qubit) ---")
for name, b in [("Z basis", [np.array([1, 0]), np.array([0, 1])]),
                ("X basis", [np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)]),
                ("(Z+X)/sqrt2", [np.array([np.cos(np.pi / 8), np.sin(np.pi / 8)]), np.array([-np.sin(np.pi / 8), np.cos(np.pi / 8)])])]:
    print(f"  {name:<12}: chi(S:E1) = {holevo(psi, b, [1]):.4f},  chi(S:E1..E3) = {holevo(psi, b, [1, 2, 3]):.4f}")
print("  => pointer = Z (the observable the interaction couples to): SELECTOR (pointer, given split).")

print("\n--- 3. Can redundancy select WHICH slot is the system? (same state, same qubit TPS) ---")
gfull = np.full(m, 1.0); psiG = evolve(gfull, np.pi / 4)   # perfect records -> GHZ
for label, ps in [("partial records (t as above)", psi), ("perfect records (GHZ)", psiG)]:
    Rs = []
    for s in range(n):
        env = [q for q in range(n) if q != s]
        HS_, cur_ = plateau(ps, [s], env, samples=6)
        Rs.append(redundancy(HS_, cur_, m))
    print(f"  {label:<30}: R_0.1 with slot j as 'system', j = 0..{m}: " + " ".join(f"{r:4.1f}" for r in Rs))
print("  => perfect records: every slot is equally a 'system' with maximal redundancy (GHZ symmetry) - split NOT selected.")

print("\n--- 4. Redundancy across TPSs of the same global state ---")
# any two pure states of the same dimension are related by a global unitary, so redundancy can be anything
# depending on the frame. Explicit: GHZ in frame 1; product state in frame 2 (disentangling CNOT ladder + H).
def cnot(c, tq):
    P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
    return op({c: P0}) + op({c: P1, tq: X})
Ucl = np.eye(2 ** n, dtype=complex)
for k in range(1, n):
    Ucl = cnot(0, k) @ Ucl
ghz = np.zeros(2 ** n, complex); ghz[0] = ghz[-1] = 1 / np.sqrt(2)
for label, st in [("frame 1 (GHZ)", ghz), ("frame 2 = CNOT-ladder frame", Ucl @ ghz)]:
    HS_, cur_ = plateau(st, [0], list(range(1, n)), samples=6)
    print(f"  {label:<30}: H(S) = {HS_:.3f}, I(S:F) by |F| = " + " ".join(f"{v:.2f}" for v in cur_) +
          f", R = {redundancy(HS_, cur_, m) if HS_ > 1e-9 else 0:.1f}")
print("  => same global state: maximal redundancy in one TPS, NONE in another. Redundancy is a TPS-relative quantity.")

print("\n--- 5. Fragment individuation (access): group the 8 env qubits into fragments ---")
HS, cur = plateau(psi, sysA, envA)
for size in [1, 2, 4]:
    frags = [envA[i:i + size] for i in range(0, m, size)]
    # smallest number j of fragments whose union gives <I> >= 0.9 H(S) (averaged over choices, as in section 1)
    need = None
    for j in range(1, len(frags) + 1):
        avg = np.mean([MI(psi, sysA, sum(c, [])) for c in it.combinations(frags, j)])
        if avg >= 0.9 * HS:
            need = j; break
    print(f"  fragments of {size} qubit(s) ({len(frags)} fragments): R_0.1 = {len(frags) / need if need else 0:.2f}")
print("  => R counts fragments; the number depends on how the observer's accessible fragments are individuated (A).")

print("\n--- 6. Environment self-interaction (D): scrambling destroys local records ---")
Henv = sum(rng.normal() * op({a: P, b: Q}) for a, b in it.combinations(range(1, n), 2)
           for P in (X, Y, Z) for Q in (X, Y, Z)) * 0.3
for tag, extra, tt in [("no env-env coupling", None, t), ("env all-to-all scrambling, same t", Henv, t),
                       ("env scrambling, record then scramble 3t", None, None)]:
    if tt is None:
        ps = expm(-1j * Henv * 3 * t) @ psi
    else:
        ps = evolve(g, tt, extra)
    HS_, cur_ = plateau(ps, sysA, envA, samples=6)
    print(f"  {tag:<38}: I(S:F) = " + " ".join(f"{v:.2f}" for v in cur_) + f"  H(S)={HS_:.2f}  R={redundancy(HS_, cur_, m):.2f}")
print("  => redundancy needs a non-scrambling environment of separately accessible fragments (a D + A condition).")
