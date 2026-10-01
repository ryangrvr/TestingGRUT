"""SCOUT-1 W1-L: does local tomography (with composition) select the lift's number field?
 1. parameter counting K(d): real d(d+1)/2, complex d^2, quaternionic d(2d-1); test K_AB = K_A K_B.
 2. explicit real-QM failure: two rebit-pair states with identical statistics for ALL local real observables,
    distinguishable by a global real observable (sigma_y x sigma_y).
 3. complex QM: local Pauli products are a complete operator basis (rank check).
"""
import itertools
import numpy as np

K = {"real": lambda d: d * (d + 1) // 2, "complex": lambda d: d * d, "quaternionic": lambda d: d * (2 * d - 1)}
print("=== 1. state-space dimension K(d) (number of real parameters of unnormalized states) ===")
for f, Kf in K.items():
    rows = []
    for dA, dB in ((2, 2), (2, 3), (3, 3), (2, 4)):
        KAB = Kf(dA * dB) if f != "quaternionic" else None
        rows.append(f"(dA,dB)=({dA},{dB}): K_A*K_B={Kf(dA) * Kf(dB)}" + (f", K(dA*dB)={KAB}" if KAB is not None else ""))
    print(f"  {f:12s}: " + "; ".join(rows))
print("  complex: K_AB = K_A K_B exactly (d^2 is multiplicative). real: K(dA dB) > K_A K_B (global-only parameters).")
print("  quaternionic: no tensor product of quaternionic modules over H (non-commutative); composite proposals")
print("  (literature: quaternionic composites are not locally tomographic) -- NOT checked in this script.")
for d in (2, 3, 4, 5):
    print(f"   d={d}: real deficit K(d^2) - K(d)^2 = {K['real'](d * d) - K['real'](d) ** 2}")

print("\n=== 2. real QM fails local tomography (rebits) ===")
I = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1., -1]); Y = np.array([[0, -1j], [1j, 0]])
YY = np.real(np.kron(Y, Y))      # real symmetric!
rho_p = (np.eye(4) + YY) / 4
rho_m = (np.eye(4) - YY) / 4
for name, r in (("rho+", rho_p), ("rho-", rho_m)):
    print(f"  {name}: real={np.allclose(r.imag if np.iscomplexobj(r) else 0, 0)}, eigenvalues {np.round(np.linalg.eigvalsh(r), 4)}")
local_real = [I, X, Z]           # real symmetric 2x2 basis (sigma_y is not real)
diffs = [abs(np.trace((rho_p - rho_m) @ np.kron(a, b))) for a in local_real for b in local_real]
print("  max |<A x B>_+ - <A x B>_-| over all local real observables:", max(diffs))
print("  global real observable sigma_y x sigma_y:", np.trace(rho_p @ YY).real, "vs", np.trace(rho_m @ YY).real)
print("  => identical local statistics, globally distinguishable: local tomography FAILS for real QM")

print("\n=== 3. complex QM: local products span all Hermitian operators ===")
paulis = [I, X, Y, Z]
for n in (2, 3):
    ops = [np.kron(np.kron(*c[:2]), c[2]) if n == 3 else np.kron(*c) for c in itertools.product(paulis, repeat=n)]
    M = np.array([o.ravel() for o in ops])
    print(f"  {n} qubits: rank of local Pauli products = {np.linalg.matrix_rank(M)} = 4^{n} = {4 ** n}")
