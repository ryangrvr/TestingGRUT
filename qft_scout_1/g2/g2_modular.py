"""QFT-SCOUT-1 G2 -- modular flow, H_epoch and orientation: finite/type-I controls + lattice Bisognano-Wichmann illustration.
Conventions: Tomita-Takesaki sigma_t(A) = Delta^{it} A Delta^{-it}; for (B(C^n), rho) this is rho^{it} A rho^{-it}.
Lattice parts are ILLUSTRATIONS of continuum theorems (BW), never adjudicating."""
import numpy as np
from scipy.linalg import expm, logm, sqrtm
rng = np.random.default_rng(20261006)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def rand_rho(n):
    G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); R = G @ G.conj().T; return R / np.trace(R).real
def rand_herm(n):
    G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return (G + G.conj().T) / 2
def mpow(R, z):
    w, U = np.linalg.eigh(R); return (U * w.astype(complex) ** z) @ U.conj().T
n = 4; I = np.eye(n)

hdr("G2-1 finite / type-I control: standard form of (B(C^n), rho) on C^n (x) C^n")
rho = rand_rho(n); p, U = np.linalg.eigh(rho)
Omega = sum(np.sqrt(p[k]) * np.kron(U[:, k], U[:, k].conj()) for k in range(n))          # purification, |Omega> in H (x) H
Delta = np.kron(rho, np.linalg.inv(rho).T)                                                 # modular operator for M = B (x) 1
def J(v):                                                                                    # modular conjugation: a (x) b -> conj-swap in rho's basis
    Vm = v.reshape(n, n); return (Vm.conj().T).reshape(-1)
A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); B = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
MA = np.kron(A, I)
S_lhs = np.kron(A.conj().T, I) @ Omega                                                      # S A Omega = A* Omega
S_rhs = J(mpow(Delta, 0.5) @ (MA @ Omega))
print(f"  Tomita relation S = J Delta^(1/2): |J Delta^(1/2) A Omega - A* Omega| = {np.linalg.norm(S_lhs - S_rhs):.1e}; Delta Omega = Omega: {np.linalg.norm(Delta @ Omega - Omega):.1e}")
t = 0.37; Dt = mpow(Delta, 1j * t)
flowM = Dt @ MA @ mpow(Delta, -1j * t); target = np.kron(mpow(rho, 1j * t) @ A @ mpow(rho, -1j * t), I)
print(f"  sigma_t on M is rho^(it) A rho^(-it): error {np.linalg.norm(flowM - target):.1e}  (an INNER flow: no spacetime, no physical-time content)")
w = lambda X: (Omega.conj() @ X @ Omega)
kms = abs(w(MA @ np.kron(rho @ B @ np.linalg.inv(rho), I)) - w(np.kron(B, I) @ MA))
print(f"  KMS (Tomita convention): omega(A sigma_(-i)(B)) = omega(B A): error {kms:.1e}  -> omega is KMS at beta = -1 for sigma_t, i.e. beta = +1 for tau_t = sigma_(-t)")
hdr("G2-2 state dependence (same algebra, different faithful states)")
rho2 = rand_rho(n)
K1, K2 = -logm(rho), -logm(rho2)
print(f"  modular Hamiltonians K = -log rho differ: |K1 - K2| = {np.linalg.norm(K1 - K2):.3f}")
Kany = rand_herm(n); rhoK = expm(-Kany); rhoK /= np.trace(rhoK)
print(f"  ANY Hermitian K is the modular Hamiltonian of rho_K = e^-K/Z: |(-log rho_K) - K| mod scalar = "
      f"{np.linalg.norm((-logm(rhoK) - Kany) - np.trace(-logm(rhoK) - Kany) / n * I):.1e} -> in type I every inner flow is modular for some state")
ut = mpow(rho2, 1j * t) @ mpow(rho, -1j * t)
s2 = mpow(rho2, 1j * t) @ A @ mpow(rho2, -1j * t); s1 = mpow(rho, 1j * t) @ A @ mpow(rho, -1j * t)
print(f"  Connes cocycle (finite): sigma2_t = Ad(u_t) sigma1_t with u_t = rho2^(it) rho1^(-it) in M: error {np.linalg.norm(s2 - ut @ s1 @ ut.conj().T):.1e}; u_t unitary: {np.linalg.norm(ut @ ut.conj().T - I):.1e}")
print("  -> type I: modular flows are inner and STATE-PRICED (the flow carries exactly the state's information).")
hdr("G2-3 sign inversion")
MBp = np.kron(I, B)                                                                          # element of the commutant M'
flowMp = Dt @ MBp @ mpow(Delta, -1j * t); rT = rho.T
print(f"  commutant M': Delta^(it) (1 (x) B) Delta^(-it) = 1 (x) rho^T^(-it) B rho^T^(it): error "
      f"{np.linalg.norm(flowMp - np.kron(I, mpow(rT, -1j * t) @ B @ mpow(rT, 1j * t))):.1e} -> the commutant's own modular flow (Delta' = Delta^-1) runs the opposite way")
Jop = np.array([J(e) for e in np.eye(n * n)]).T                                              # matrix of J composed with conj
v = rng.normal(size=n * n) + 1j * rng.normal(size=n * n)
lhs = J(Dt @ v); rhs = Dt @ J(v)
print(f"  J commutes with Delta^(it) (antilinear: J Delta^(it) J = Delta^(it)): error {np.linalg.norm(lhs - rhs):.1e}; J maps M to M' (J (A(x)1) J = 1(x)conj-transpose form)")
# anti-automorphism alpha(A) = A^T of B(C^n): sigma^(phi o alpha)_t = alpha^-1 o sigma^phi_(-t) o alpha
lhs_a = mpow(rho.T, 1j * t) @ A @ mpow(rho.T, -1j * t)                     # sigma_t of the state phi o alpha = tr(rho^T .)
rhs_a = (mpow(rho, -1j * t) @ A.T @ mpow(rho, 1j * t)).T                    # alpha^-1 sigma_-t^phi alpha (A)
print(f"  anti-automorphism alpha(A) = A^T: sigma^(phi o alpha)_t = alpha^-1 sigma^phi_(-t) alpha: error {np.linalg.norm(lhs_a - rhs_a):.1e} -> anti-isomorphisms reverse modular time")
print("  -> the parametrization sign of modular time flips under algebra <-> commutant (which side is 'the system': a Sigma / A_partition choice)")
print("  passivity test (Pusz-Woronowicz form, finite): work extractable by cyclic unitaries U from rho under generator H")
def max_work(rho, H, trials=4000):
    best = 0.0; E0 = np.trace(rho @ H).real
    for _ in range(trials):
        Uq = expm(1j * rand_herm(n) * rng.uniform(0, 3)); best = max(best, E0 - np.trace(Uq @ rho @ Uq.conj().T @ H).real)
    w_, V_ = np.linalg.eigh(rho); e_, W_ = np.linalg.eigh(H)                                 # exact optimum: pair largest p with lowest e
    exact = E0 - np.sum(np.sort(w_)[::-1] * np.sort(e_))
    return best, exact
for nm, H in [("H = +K (modular Hamiltonian, tau_t = sigma_-t)", K1), ("H = -K (reversed orientation)", -K1)]:
    b, ex = max_work(rho, H)
    print(f"   {nm}: max extractable work (random search {b:.4f}; exact {ex:.4f})")
print("  -> passivity selects H = +K, i.e. the orientation tau_t = sigma_(-t) (positive temperature); the reversed orientation is maximally active.")
print("     But passivity (no work from cycles) is itself a Kelvin-type arrow postulate: orientation is RELOCATED into it, not derived.")

hdr("G2-4 lattice illustration (80-digit precision): modular Hamiltonian of an interval vs boost-like / conformal profiles")
import mpmath as mp
mp.mp.dps = 80
def interval_Hp(Nn, Rn, mu, beta=None):
    K = mp.matrix(Nn, Nn)
    for i in range(Nn):
        K[i, i] = mu ** 2 + 2
        if i + 1 < Nn: K[i, i + 1] = K[i + 1, i] = -1
    lam, V = mp.eigsy(K)
    f = lambda l: (1 if beta is None else mp.coth(beta * mp.sqrt(l) / 2))
    R = list(range(Nn // 2, Nn // 2 + Rn)); X = mp.matrix(Rn, Rn); P = mp.matrix(Rn, Rn)
    for a, i in enumerate(R):
        for b, j in enumerate(R):
            X[a, b] = sum(V[i, k] * V[j, k] * lam[k] ** -0.5 * f(lam[k]) for k in range(Nn)) / 2
            P[a, b] = sum(V[i, k] * V[j, k] * lam[k] ** 0.5 * f(lam[k]) for k in range(Nn)) / 2
    w, O = mp.eigsy(X); Xh = O * mp.diag([mp.sqrt(x) for x in w]) * O.T
    nu2, Q = mp.eigsy(Xh * P * Xh); nu = [mp.sqrt(x) for x in nu2]
    eps = [mp.log((2 * v + 1) / (2 * v - 1)) for v in nu]
    return Xh * Q * mp.diag([e / v for e, v in zip(eps, nu)]) * Q.T * Xh
Nn, Rn, mu = 120, 30, mp.mpf("0.02"); js = [0, 1, 2, 5, 10, 15, 20, 25]
for label, beta in [("vacuum", None), ("thermal beta=4 (same algebra, different state)", mp.mpf(4))]:
    Hp = interval_Hp(Nn, Rn, mu, beta)
    lin = [float(Hp[j, j] / (2 * mp.pi * (j + 0.5))) for j in js]
    par = [float(Hp[j, j] / (2 * mp.pi * (j + 0.5) * (Rn - j - 0.5) / Rn)) for j in js]
    print(f"  {label}: H_p[j,j] / 2pi x            at j = {js}: {np.round(lin, 3)}")
    print(f"  {' ' * len(label)}  H_p[j,j] / 2pi x(l-x)/l   at j = {js}: {np.round(par, 3)}")
print("  -> vacuum: near the cut the modular Hamiltonian is 2 pi x (boost-like, BW); across the interval it tracks the conformal")
print("     two-ended profile 2 pi x(l-x)/l (geometric only for conformal theories; the small mass mu = 0.02 gives deviations);")
print("     a thermal state on the SAME algebra is not boost-like. The geometric identification is a property of the vacuum +")
print("     covariance + spectrum condition (supplied), not of the algebra alone.")
