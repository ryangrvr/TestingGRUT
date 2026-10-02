"""QFT-SCOUT-1 G3 -- half-sided modular inclusions (HSMI): finite/type-I obstruction + light-ray sign/generator toy.
Convention (Wiesbrock CMP 157, 1993): N subset M is (-)hsm iff Delta_M^{-it} N Delta_M^{it} subset N for t >= 0;
(+)hsm iff the same holds for t <= 0.  Code convention sigma_t = Ad Delta^{it}."""
import numpy as np
from scipy.linalg import expm, logm
rng = np.random.default_rng(20261007)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def rand_rho(n):
    G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); R = G @ G.conj().T; return R / np.trace(R).real
def rand_unitary(n):
    Q, R = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))); return Q * (np.diag(R) / abs(np.diag(R)))
def mpow(R, z):
    w, U = np.linalg.eigh(R); return (U * w.astype(complex) ** z) @ U.conj().T
def basis_of(mats):
    V = np.array([m.reshape(-1) for m in mats]).T; Q, _ = np.linalg.qr(V); return Q
def defect(rho, Nbasis_mats, t):
    """size of sigma_t(N) outside N: max over a basis of || (1 - P_N) vec(sigma_t(x)) ||"""
    Q = basis_of(Nbasis_mats); Ut = mpow(rho, 1j * t)
    worst = 0.0
    for x in Nbasis_mats:
        v = (Ut @ x @ Ut.conj().T).reshape(-1); worst = max(worst, np.linalg.norm(v - Q @ (Q.conj().T @ v)) / np.linalg.norm(v))
    return worst

hdr("G3-2 FINITE / TYPE-I CONTROL: can a non-trivial half-sided modular inclusion exist?")
print("  Lemma F1 (proved): if log Delta is bounded and sigma_t(N) subset N for t in [0, eps), the derivation d = i[log Delta, .]")
print("  maps N into N (derivative at 0+ of a curve staying in the closed subspace N); hence sigma_t(N) = exp(t d)(N) subset N for")
print("  ALL real t. In finite dimensions log Delta is bounded: HALF-SIDED INVARIANCE => TWO-SIDED INVARIANCE.")
print("  Lemma F2 (proved): finite dimensions, Omega separating for M and cyclic for N subset M => dim N >= dim H >= dim M, so N = M.")
n = 4; E = lambda i, j, d: np.outer(np.eye(d)[i], np.eye(d)[j])
I2 = np.eye(2)
rho = rand_rho(n)
print("\n  (a) random subalgebras N = U (B(C^2) (x) 1) U^+ inside M = B(C^4), random faithful rho: defect of sigma_t(N) from N")
asym = []
for trial in range(5):
    U = rand_unitary(n); Nm = [U @ np.kron(E(i, j, 2), I2) @ U.conj().T for i in range(2) for j in range(2)]
    ds = {t: defect(rho, Nm, t) for t in (-1.0, -0.1, -0.01, 0.01, 0.1, 1.0)}
    asym.append(abs(ds[0.01] - ds[-0.01]) / ds[0.01])
    print("   trial", trial, " ".join(f"t={t:+.2f}:{v:.2e}" for t, v in ds.items()))
print(f"  -> never half-sided: the defect is present for BOTH signs, symmetric to first order (relative asymmetry at |t|=0.01: max {max(asym):.2e})")
print("\n  (b) a modular-invariant subalgebra (rho = rho_a (x) rho_b, N = B(C^2) (x) 1): defect for both signs")
ra, rb = rand_rho(2), rand_rho(2); rho_p = np.kron(ra, rb)
Nm = [np.kron(E(i, j, 2), I2) for i in range(2) for j in range(2)]
print("   " + " ".join(f"t={t:+.2f}:{defect(rho_p, Nm, t):.1e}" for t in (-1.0, -0.1, 0.1, 1.0)) + "  -> invariance is two-sided (consistent with F1)")
print("\n  (c) standardness: in the standard form H = C^4 (x) C^4 (dim 16), the cyclic subspace of N = B(C^2) (x) 1 (x) 1 has dimension")
p, V = np.linalg.eigh(rho_p)
Om = sum(np.sqrt(p[k]) * np.kron(V[:, k], V[:, k].conj()) for k in range(n))
NOm = np.array([np.kron(x, np.eye(n)) @ Om for x in Nm]).T
print(f"   dim(N Omega) = {np.linalg.matrix_rank(NOm)} < 16 = dim H: Omega is NOT cyclic for N (F2) -> N cannot be standard unless N = M")
print("  -> FINITE / TYPE-I OBSTRUCTION: a half-sided modular inclusion needs an UNBOUNDED modular generator; with standard N it")
print("     is impossible in finite dimensions (it collapses to N = M). Import (Wiesbrock 1993): a non-trivial HSMI forces type III_1.")

hdr("G3-1 / G3-6 light-ray toy: chiral modular flows act geometrically (BW / Hislop-Longo, supplied)")
print("  modular flow of the right ray (a, inf): x -> a + e^{-2 pi t}(x - a);  of the left ray (-inf, b) [commutant]: x -> b + e^{+2 pi t}(x - b)")
def flow(kind, c, t, x): return c + np.exp((-1 if kind == "R" else 1) * 2 * np.pi * t) * (x - c)
def inside(kind, c, x): return (x > c) if kind == "R" else (x < c)
def hsm_sign(M, N, ts=np.linspace(0.01, 2, 50), xs=None):
    """check Ad Delta_M^{-it} N subset N for t >= 0 (minus) or t <= 0 (plus), geometrically on sample points of N"""
    (kM, cM), (kN, cN) = M, N
    xs = (cN + np.linspace(0.001, 5, 400)) if kN == "R" else (cN - np.linspace(0.001, 5, 400))
    minus = all(inside(kN, cN, flow(kM, cM, -t, xs)).all() for t in ts)          # Delta^{-it}: parameter -t
    plus = all(inside(kN, cN, flow(kM, cM, +t, xs)).all() for t in ts)
    return ("-hsm" if minus else "") + ("+hsm" if plus else "") or "not hsm"
def gen(kind, c):   # vector field of t -> flow (coefficient of d/dx as a function of x): d/dt flow at t=0
    return (lambda x: -2 * np.pi * (x - c)) if kind == "R" else (lambda x: 2 * np.pi * (x - c))
xs = np.linspace(-3, 3, 7)
for nm, M, N in [("right rays: N=(1,inf) subset M=(0,inf)", ("R", 0.0), ("R", 1.0)),
                 ("left rays (mirror): N=(-inf,-1) subset M=(-inf,0)", ("L", 0.0), ("L", -1.0)),
                 ("right rays, larger gap: N=(2.5,inf) subset M=(0,inf)", ("R", 0.0), ("R", 2.5))]:
    d = (gen(*N)(xs) - gen(*M)(xs)) / (2 * np.pi)
    print(f"  {nm}: {hsm_sign(M, N)};  (G_N - G_M)/2pi as a vector field = {np.round(d, 6)} (constant -> a TRANSLATION, coefficient {d[0]:+.2f})")
print("  -> one inclusion yields a translation generator (scale = the inclusion gap): the generator of x -> x + s, i.e. the same")
print("     rightward null translation for BOTH mirrored configurations; the hsm SIGN records whether N = U(+gap) M U(-gap) (-hsm)")
print("     or N = U(-gap) M U(+gap) (+hsm), i.e. the relation between the inclusion order and the positive-energy direction.")
print("  -> the family M(s) = U(s) M U(-s), s in R, is a totally ordered continuum generated from ONE inclusion (a light-ray net).")
