"""SCOUT-2 S2-7: consistent histories.  Does consistency select a unique quasi-classical history set?

Two-time histories: projectors P_a(theta1) at t1 and P_b(theta2) at t2 onto a qubit basis at angle theta in the
x-z plane, acting on a chosen slot.  Decoherence functional D(ab, a'b') = <psi0| C_a'b'^dag C_ab |psi0>,
C_ab = P_b(t2) P_a(t1) (Heisenberg).  Consistency measure: max_{(ab)!=(a'b')} |D| / sqrt(p p').
Model: S + 4 environment qubits, H = omega X_S + sum_k g_k Z_S Y_k (system has its own dynamics).
"""
import itertools as it
import numpy as np
from scipy.linalg import expm

rng = np.random.default_rng(7)
m = 4; n = m + 1
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)


def op(d):
    out = np.eye(1)
    for q in range(n):
        out = np.kron(out, d.get(q, I2))
    return out


def proj(theta):
    v0 = np.array([np.cos(theta / 2), np.sin(theta / 2)]); v1 = np.array([-np.sin(theta / 2), np.cos(theta / 2)])
    return [np.outer(v0, v0).astype(complex), np.outer(v1, v1).astype(complex)]


def heis(P, U):
    return U.conj().T @ P @ U


def dmax(psi0, P1s, P2s, U1, U2):
    C = {}
    for a, P1 in enumerate(P1s):
        for b, P2 in enumerate(P2s):
            C[(a, b)] = heis(P2, U2) @ heis(P1, U1) @ psi0
    keys = list(C)
    p = {k: np.vdot(C[k], C[k]).real for k in keys}
    worst = 0.0
    for k1, k2 in it.combinations(keys, 2):
        if p[k1] > 1e-10 and p[k2] > 1e-10:
            worst = max(worst, abs(np.vdot(C[k2], C[k1])) / np.sqrt(p[k1] * p[k2]))
    return worst


omega = 1.0
g = rng.uniform(0.8, 1.2, m)
t1, t2 = 1.3, 2.9
psi0 = np.zeros(2 ** n, complex); psi0[0] = 1
psi0 = op({0: np.array([[np.cos(.4), -np.sin(.4)], [np.sin(.4), np.cos(.4)]], complex)}) @ psi0
grid = np.linspace(0, np.pi, 37)
EPS = 0.05

print("=" * 78); print("S2-7  consistent histories: S + 4 env qubits, two-time histories"); print("=" * 78)
for label, gs in [("closed system (g = 0)", 0 * g), ("monitored system (g ~ 1)", g), ("strongly monitored (g ~ 4)", 4 * g)]:
    H = omega * op({0: X}) + sum(gs[k] * op({0: Z, k + 1: Y}) for k in range(m))
    U1, U2 = expm(-1j * H * t1), expm(-1j * H * t2)
    for slot in [0, 1]:
        M = np.array([[dmax(psi0, [op({slot: P}) for P in proj(a)], [op({slot: P}) for P in proj(b)], U1, U2)
                       for b in grid] for a in grid])
        frac = (M < EPS).mean()
        zz = M[0, 0]
        cons = [(round(grid[i] / np.pi, 3), round(grid[j] / np.pi, 3)) for i, j in zip(*np.where(M < EPS))]
        # record quality: overlap of environment branch states at t1 (product over env qubits)
        print(f"  {label:<26} slot {slot} ({'S' if slot == 0 else 'E1'}): fraction of (theta1,theta2) consistent at eps={EPS}: "
              f"{frac:.3f};  Z-Z histories: {zz:.3f};  #consistent sets {len(cons)}")
        if slot == 0 and gs.any():
            far = [c for c in cons if min(c[0], 1 - c[0]) > 0.15 or min(c[1], 1 - c[1]) > 0.15]
            print(f"     consistent sets far (>0.15 pi) from the pointer Z: {len(far)}  e.g. {far[:4]}")

print("\n--- Dowker-Kent-type sets: global projectors onto the evolving state itself ---")
H = omega * op({0: X}) + sum(g[k] * op({0: Z, k + 1: Y}) for k in range(m))
U1, U2 = expm(-1j * H * t1), expm(-1j * H * t2)
# projector onto a random superposition |chi> containing psi(t) component: P = |chi><chi|, 1-P
for trial in range(3):
    chi1 = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n); chi1 /= np.linalg.norm(chi1)
    psi_t1 = U1 @ psi0
    # build P1 = projector onto span{psi(t1)} rotated slightly into random direction; P2 onto psi(t2)-based
    v = psi_t1 + 0.6 * chi1; v /= np.linalg.norm(v)
    P1 = np.outer(v, v.conj()); P1s = [P1, np.eye(2 ** n) - P1]
    # Schroedinger-picture projectors -> Heisenberg via heis inside dmax expects Schroedinger ops at t1, t2
    w = U2 @ U1.conj().T @ (P1 @ psi_t1); w /= np.linalg.norm(w)
    P2 = np.outer(w, w.conj()); P2s = [P2, np.eye(2 ** n) - P2]
    print(f"  trial {trial}: max |D| = {dmax(psi0, P1s, P2s, U1, U2):.2e}  (P1 = projector onto a random-tilted"
          f" global vector; P2 = its evolved image) -> exactly consistent, no TPS, no record structure")
print("  => consistency admits sets built from arbitrary global (non-local, non-record) projectors.")
print("     Selecting the quasi-classical set needs a supplied projector family: local slot (split), a basis class")
print("     (coarse-graining), and the time sequence - all access data.")

print("\n--- Persistent records: Z-Z consistency vs environment size (omega = 0.3, g_k ~ U(2,6)) ---")
grid13 = np.linspace(0, np.pi, 13)
for mm in [1, 2, 3, 5, 7]:
    m = mm; n = mm + 1
    gk = rng.uniform(2, 6, mm)
    H = 0.3 * op({0: X}) + sum(gk[k] * op({0: Z, k + 1: Y}) for k in range(mm))
    U1, U2 = expm(-1j * H * t1), expm(-1j * H * t2)
    p0 = np.zeros(2 ** n, complex); p0[0] = 1
    p0 = op({0: np.array([[np.cos(.4), -np.sin(.4)], [np.sin(.4), np.cos(.4)]], complex)}) @ p0
    zz = dmax(p0, [op({0: P}) for P in proj(0)], [op({0: P}) for P in proj(0)], U1, U2)
    M = np.array([[dmax(p0, [op({0: P}) for P in proj(a)], [op({0: P}) for P in proj(b)], U1, U2) for b in grid13] for a in grid13])
    cons = [(round(grid13[i] / np.pi, 3), round(grid13[j] / np.pi, 3)) for i, j in zip(*np.where(M < EPS))]
    nonZ = [c for c in cons if c not in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]]
    print(f"  m = {mm}: Z-Z max|D| = {zz:.3f};  consistent (13x13 grid, eps={EPS}): {len(cons)}; of which not Z-Z: {len(nonZ)};"
          f" distinct theta1/pi among consistent: {sorted(set(c[0] for c in cons))}")
print("  note: off-diagonal terms with b != b' vanish identically (orthogonal final projectors), so the FINAL-time")
print("  basis is never constrained by consistency: it is a free (access) choice in every consistent set.")
