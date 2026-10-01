"""SCOUT-1 D3 computational ablation checks (literature statements are SECONDARY; these are the parts verifiable here).

 1. Purification: does it hold in classical / real / complex theory? (existence + uniqueness up to reversible maps on
    the purifying system).  -> tests whether purification is the FIELD selector or the quantum-vs-classical selector.
 2. Local tomography in fermionic quantum theory (parity superselection): two even states with identical local
    (parity-even) statistics but different global coherence.
 3. Bloch-ball bits (spin factors V_d): K = d + 1. A locally tomographic composite of two such bits has K_AB = (d+1)^2.
    Compare with the composite dimension of the Jordan family on 4 levels: real 10, complex 16, quaternionic 28.
 4. Born accounting: in a GPT the outcome probability is the affine pairing e(omega). With a self-dual Jordan cone,
    e(omega) = Tr(E rho). The W1-G non-Born family h_eps is not affine in the state, so it is excluded by the
    FRAMEWORK premise (linearity), not by any reconstruction axiom.
"""
import numpy as np

rng = np.random.default_rng(8)

print("=== 1. purification ===")
# classical bit: a pure joint state is a point mass, so its marginals are point masses -> a mixed marginal has no purification
p = np.array([0.3, 0.7])
print("  classical: mixed marginal", p, "has a pure (point-mass) joint extension:", False,
      "(every marginal of a point mass is a point mass)")
for field in ("real", "complex"):
    d = 3
    if field == "real":
        A = rng.standard_normal((d, d))
    else:
        A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    rho = A @ A.conj().T; rho /= np.trace(rho).real
    w, V = np.linalg.eigh(rho)
    psi = sum(np.sqrt(max(w[i], 0)) * np.kron(V[:, i], np.eye(d)[i]) for i in range(d))
    M = psi.reshape(d, d)
    marg = M @ M.conj().T
    # second purification: rotate the purifier by a random orthogonal / unitary
    Q = np.linalg.qr(rng.standard_normal((d, d)) + (0 if field == "real" else 1j * rng.standard_normal((d, d))))[0]
    psi2 = (M @ Q.T).reshape(-1)
    M2 = psi2.reshape(d, d)
    # recover the connecting reversible map: M2 = M R  =>  R = pinv(M) M2 ; check R is orthogonal/unitary
    R = np.linalg.pinv(M) @ M2
    print(f"  {field:7s}: purification in the {field} theory reproduces rho: {np.allclose(marg, rho)};  "
          f"purification entries real: {np.allclose(psi.imag, 0) if field == 'real' else 'n/a'};  "
          f"two purifications related by a {('orthogonal' if field == 'real' else 'unitary')} map on the purifier: "
          f"{np.allclose(R @ R.conj().T, np.eye(d), atol=1e-8)}")
print("  => purification holds in real AND complex QM: it separates quantum from classical, NOT complex from real")

print("\n=== 2. fermionic local tomography (2 modes, parity superselection) ===")
# basis |00>, |01>, |10>, |11>; even sector {00, 11}
rho_p = np.zeros((4, 4)); rho_m = np.zeros((4, 4))
for r, s in ((rho_p, +1), (rho_m, -1)):
    r[0, 0] = r[3, 3] = 0.5; r[0, 3] = r[3, 0] = 0.5 * s
n = np.diag([0., 1.]); I2 = np.eye(2)
local_even = [I2, n]                     # parity-even local operators of one mode (odd ones a, a+ are not observable)
diffs = [abs(np.trace((rho_p - rho_m) @ np.kron(A, B))) for A in local_even for B in local_even]
# global even observable built from two odd local operators: a1+ a2+ + h.c. (with JW sign irrelevant for |00>-|11>)
G = np.zeros((4, 4)); G[0, 3] = G[3, 0] = 1.0
print("  max local-even statistic difference:", max(diffs), "; global even observable:",
      np.trace(rho_p @ G), "vs", np.trace(rho_m @ G))
print("  => fermionic QT fails local tomography (like real QM) while satisfying purification (SECONDARY: D'Ariano et al.)")

print("\n=== 3. Bloch-ball bits V_d and locally tomographic composites ===")
for d, name in ((2, "real (rebit)"), (3, "complex (qubit)"), (5, "quaternionic")):
    K1 = d + 1
    print(f"  d = {d} {name:16s}: K_bit = {K1}, locally tomographic composite needs K_AB = {K1 ** 2}")
for name, K2, K4, KAB in (("real", 3, 10, 9), ("complex", 4, 16, 16), ("quaternionic", 6, 28, 36)):
    print(f"  {name:12s}: bit = 2x2 Hermitian, dim {K2}; natural 4-level composite = 4x4 Hermitian, dim {K4}; "
          f"locally tomographic needs {KAB} -> {'MATCH' if K4 == KAB else 'mismatch'}")
print("  real composite is too big (10 > 9: a global-only parameter, W1-L); quaternionic 4x4 is too small (28 < 36);")
print("  only complex matches -> within Jordan/HSD + ball bits, local tomography singles out complex QM")
print("  (consistent with Barnum-Wilce, SECONDARY). The dimension count is checked here; the theorem itself is not.")

print("\n=== 4. Born accounting ===")
h = lambda p, eps: p + eps * np.sin(2 * np.pi * p) / (2 * np.pi)
E = np.diag([1., 0.])
r1 = np.diag([1., 0.]); r2 = np.array([[0.5, 0.5], [0.5, 0.5]])
mix = 0.5 * (r1 + r2)
for eps in (0.0, 0.3):
    lhs = h(np.trace(mix @ E), eps); rhs = 0.5 * (h(np.trace(r1 @ E), eps) + h(np.trace(r2 @ E), eps))
    print(f"  eps = {eps}: h(p(mixture)) = {lhs:.4f}, mixture of h(p) = {rhs:.4f}  -> affine: {np.isclose(lhs, rhs)}")
print("  => the GPT framework (outcome probability affine in the state = operational mixing) already excludes every")
print("     non-Born h; given self-duality, p(E|rho) = Tr(E rho). Born is bought by the FRAMEWORK premise, not by an axiom")
