"""GRAVITY-SCOUT-1 / G2 -- finite control for gravitational subsystem structure.

FINITE-DIMENSIONAL ILLUSTRATION ONLY (type I throughout). It does not model gravity; it isolates the algebraic
mechanisms that the G2 sources attribute to gravity, so that each priced ingredient can be switched on/off.
Independent code path, not independent reviewer.

Inside  H_in  = two qubits, "charge" H_in = diag(0,1,1,2)  (sectors E = 0, 1, 2 of dims 1, 2, 1)
Outside H_out = five levels, nondegenerate, no accidental degeneracy with H_in.

T0  G = 0             : A_out = 1 (x) B(K).                     (flat control: ordinary split / tensor factor)
T1  charge dressing   : A_out = alg{1 (x) B(K), H_in (x) 1}.   (Gauss law: the outside measures the inside charge;
                        caricature of Donnelly-Giddings first-order dressing with total charges)
T2  exact boundary H  : A_out = alg{1 (x) B(K), P_0},  P_0 = ground-state projector of
                        H_tot = H_in + H_out + lam V.          (caricature of Raju's holography-of-information
                        ingredients: H a boundary term, unique vacuum, outside-cyclic vacuum)
Algebras are computed as bicommutants  alg(S) = S''  (finite dims).
"""
import numpy as np

rng = np.random.default_rng(7)
n_in, n_out = 4, 5
n = n_in*n_out
I_in, I_out = np.eye(n_in), np.eye(n_out)
H_in = np.diag([0., 1., 1., 2.])
H_out = np.diag([0., 1.3, 2.1, 2.9, 3.7])


def herm(k):
    a = rng.normal(size=(k, k)) + 1j*rng.normal(size=(k, k)); return (a + a.conj().T)/2


def nullspace(M, tol=1e-10):
    """right null space; singular values below tol * s_max count as zero.  tol is the RESOLUTION epsilon."""
    if M.shape[0] > M.shape[1]:
        M = np.linalg.qr(M, mode="r")                    # same right null space, square k^2 x k^2
    u, s, vh = np.linalg.svd(M)
    return vh[int(np.sum(s > max(tol*s[0], 1e-12))):].conj().T      # absolute floor: M = 0 -> everything


def commutant(gens, tol=1e-10):
    """basis (as list of n x n matrices) of {X : [X, g] = 0 for all g}; row-major vec: vec(gX - Xg)."""
    k = gens[0].shape[0]; I = np.eye(k)
    # unit-norm generators with the identity part removed (it commutes with everything); a generator whose
    # traceless part is below eps is resolution-indistinguishable from a scalar and is dropped
    gens = [g/np.linalg.norm(g) for g in gens]
    gens = [g - np.trace(g)/k*I for g in gens]
    gens = [g for g in gens if np.linalg.norm(g) > tol]
    if not gens:
        return [x.reshape(k, k) for x in np.eye(k*k)]
    M = np.vstack([np.kron(g, I) - np.kron(I, g.T) for g in gens])
    N = nullspace(M, tol)
    return [N[:, j].reshape(k, k) for j in range(N.shape[1])]


def algebra(gens, tol=1e-10):
    return commutant(commutant(gens, tol), tol)


def center_dim(A):
    Ac = commutant(A)
    # dim(A cap A') via dim(A) + dim(A') - dim(A + A')
    V = np.array([x.ravel() for x in A + Ac])
    return len(A) + len(Ac) - np.linalg.matrix_rank(V, tol=1e-8)


out_gens = [np.kron(I_in, herm(n_out)), np.kron(I_in, herm(n_out))]     # generate 1 (x) B(K)
kin = lambda x: np.kron(x, I_out)

print("== T0  G = 0 (flat control) ==")
A_out = algebra(out_gens); A_in = commutant(A_out)
print(f"dim A_out = {len(A_out)} (= n_out^2 = {n_out**2}),  dim A_in = {len(A_in)} (= n_in^2 = {n_in**2}),"
      f"  center(A_in) = {center_dim(A_in)}  -> type-I factor, tensor split; any marginals extend to a product")

print("== T1  charge dressing: outside also holds H_in ==")
A_out1 = algebra(out_gens + [kin(H_in)]); A_in1 = commutant(A_out1)
print(f"dim A_out = {len(A_out1)},  dim A_in = {len(A_in1)} (= 1+4+1 = sum_E n_E^2),  center(A_in) = {center_dim(A_in1)}"
      f"  (= number of charge sectors)")
H_in_in_both = all(np.allclose(kin(H_in) @ x, x @ kin(H_in)) for x in A_in1 + A_out1)
print(f"H_in commutes with A_in and A_out (shared central element): {H_in_in_both}")
rho = herm(n_in); rho = rho @ rho.conj().T; rho /= np.trace(rho)            # charge-mixed inside marginal
sig = herm(n_out); sig = sig @ sig.conj().T; sig /= np.trace(sig)
w = np.kron(rho, sig)
ev = lambda X: np.trace(w @ X).real
print(f"naive tensor state rho (x) sigma: w(H_in H_in) - w(H_in)w(H_in) = {ev(kin(H_in@H_in)) - ev(kin(H_in))**2:.4f}"
      f"  = Var_rho(H_in)  (nonzero => NOT a product state for the dressed pair)")
P1 = np.diag([0., 1., 1., 0.]); rho1 = P1 @ rho @ P1; rho1 /= np.trace(rho1)
w = np.kron(rho1, sig)
print(f"charge-sharp marginal (sector E = 1):  Var(H_in) = {ev(kin(H_in@H_in)) - ev(kin(H_in))**2:.1e}"
      f"  -> product across the dressed pair; A_in restricted to sector E=1 is B(C^2), a type-I factor")

print("== T2  exact boundary Hamiltonian: outside holds the vacuum projector P_0 of H_tot ==")
V = herm(n); V /= np.linalg.norm(V, 2)
H0 = kin(H_in) + np.kron(I_in, H_out)
for lam in [0.0, 1e-3, 1e-2, 1e-1, 3e-1]:
    Ht = H0 + lam*V
    e, U = np.linalg.eigh(Ht)
    Om = U[:, 0]; P0 = np.outer(Om, Om.conj())
    s = np.linalg.svd(Om.reshape(n_in, n_out), compute_uv=False)              # Schmidt coefficients in|out
    dims = []
    for eps in [1e-2, 1e-4, 1e-6, 1e-10]:
        A2 = algebra(out_gens + [P0], eps); dims.append((eps, len(A2), len(commutant(A2, eps))))
    # splitting of the E_in = 1 doublets (degenerate at lam = 0): needed energy resolution
    doublets = [e[i+1]-e[i] for i in range(n-1) if abs(e[i+1]-e[i]) < 0.15]
    print(f"lam={lam:<6g} Schmidt(Omega) = " + " ".join(f"{x:.2e}" for x in s) +
          "   [eps: dim alg{1(x)B(K),P_0} / dim A_in] " + " ".join(f"{e_:.0e}:{a_}/{c_}" for e_, a_, c_ in dims) +
          f"   min doublet splitting = {min(doublets) if doublets else float('nan'):.2e}")
print(f"(n^2 = {n*n}; eps = relative singular-value resolution used to decide linear (in)dependence)")
print("Reading: lam = 0 -> P_0 is a product projector; the outside gains only |0><0|_in (x) B(K) (charge-like data).")
print("lam != 0 -> Omega has full Schmidt rank n_in (outside-cyclic); b1 P_0 b2 spans all matrix units, A_out = B(H),")
print("A_in = C: no interior subsystem survives. The smallest Schmidt coefficient ~ lam^k sets the precision needed")
print("to exploit it (operator-norm amplification ~ 1/s_min): reconstruction is exact algebraically but resolution-priced.")
print("At coarse eps the computed algebra collapses to the charge-like algebra: eps is a numerical proxy for observational")
print("resolution (an ANALOGY, not a derivation of any gravitational resolution limit).")
