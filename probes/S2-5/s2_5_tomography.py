"""SCOUT-2 S2-5: why does local tomography fail in real QM, and what restores it?

Real QM: rho_pm = (I +- J (x) J)/4 on two rebits, J = [[0,-1],[1,0]] (so J (x) J = -Y (x) Y, real symmetric).
They agree on all local real statistics but differ on the global observable J (x) J (SCOUT-1 W1-L).

Mechanism tested: J is the 'imaginary unit'. Complex QM = real QM with a SHARED reference rebit carrying J
(ABW-type encoding). Give each party a local reference rebit R_A, R_B:
  party A measures real symmetric observables on A (x) R_A, party B on B (x) R_B (local real tomography of A R_A | B R_B).
  (1) references in a shared (correlated) real state  -> can rho_+ and rho_- be told apart by LOCAL statistics?
  (2) references in a product state                    -> ?
  (3) the ABW encoding: complex statistics reproduced by real QM + one global reference rebit.
"""
import itertools
import numpy as np

I2 = np.eye(2); X = np.array([[0., 1], [1, 0]]); Z = np.diag([1., -1]); J = np.array([[0., -1], [1, 0]])
JJ = np.kron(J, J)
rho_p = (np.eye(4) + JJ) / 4; rho_m = (np.eye(4) - JJ) / 4


def real_sym_basis(k):
    """basis of real symmetric 2^k x 2^k matrices built from {I, X, Z} and products with an even number of J's"""
    singles = [I2, X, Z, J]
    out = []
    for combo in itertools.product(range(4), repeat=k):
        if sum(1 for c in combo if c == 3) % 2 == 0:
            M = np.array([[1.]])
            for c in combo:
                M = np.kron(M, singles[c])
            out.append(M)
    return out


def local_stats_difference(ref_state):
    # system order A, B, R_A, R_B -> reorder to (A, R_A) | (B, R_B)
    full_p = np.kron(rho_p, ref_state); full_m = np.kron(rho_m, ref_state)
    perm = [0, 2, 1, 3]

    def reorder(M):
        T = M.reshape([2] * 8)
        T = T.transpose(perm + [p + 4 for p in perm])
        return T.reshape(16, 16)
    fp, fm = reorder(full_p), reorder(full_m)
    locA = real_sym_basis(2); locB = real_sym_basis(2)
    diffs = [abs(np.trace((fp - fm) @ np.kron(a, b))) for a in locA for b in locB]
    return max(diffs)


print("=== real QM: are rho_+ and rho_- distinguishable by LOCAL statistics? ===")
print("  without references (A | B only):",
      max(abs(np.trace((rho_p - rho_m) @ np.kron(a, b))) for a in real_sym_basis(1) for b in real_sym_basis(1)))
shared = np.zeros((4, 4)); v = np.array([0, 1, -1, 0]) / np.sqrt(2); shared = np.outer(v, v)   # real 'singlet' reference
print(f"  references in the shared real state (|01>-|10>)/sqrt2: <J (x) J>_ref = {np.trace(shared @ JJ):+.3f};"
      f" max local-statistics difference = {local_stats_difference(shared):.3f}")
prod = np.kron(np.diag([1., 0]), np.diag([1., 0]))
print(f"  references in a product state |00>: <J (x) J>_ref = {np.trace(prod @ JJ):+.3f};"
      f" max local-statistics difference = {local_stats_difference(prod):.3f}")
mixed = np.eye(4) / 4
print(f"  references maximally mixed: max local-statistics difference = {local_stats_difference(mixed):.3f}")

print("\n=== ABW-type encoding: complex QM = real QM + one GLOBAL reference rebit ===")
rng = np.random.default_rng(4)
psi = rng.standard_normal(4) + 1j * rng.standard_normal(4); psi /= np.linalg.norm(psi)
A = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4)); A = A + A.conj().T
psi_R = np.concatenate([np.kron(psi.real, [1, 0]) + np.kron(psi.imag, [0, 1])])
A_R = np.kron(A.real, I2) + np.kron(A.imag, J)
print(f"  <psi|A|psi> complex = {np.real(psi.conj() @ A @ psi):+.6f};  real encoding = {psi_R @ A_R @ psi_R:+.6f};"
      f"  A_R symmetric: {np.allclose(A_R, A_R.T)}")
print("  -> the imaginary unit is a physical, globally shared reference rebit; local tomography holds exactly when every")
print("     party has (correlated) access to that frame. LT is therefore a property of a SHARED REFERENCE STATE (C5-H),")
print("     not of the single-system kinematics.")
