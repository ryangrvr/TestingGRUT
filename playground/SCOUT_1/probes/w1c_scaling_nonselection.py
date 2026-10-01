"""SCOUT-1 W1-C: deformation-group invariance of the earned predicates, and what that forbids.

Deformation generators acting on a declared model (K, retained site r):
  S_lam : K -> lam*K, t -> t/lam                (time/rate rescaling, lam > 0)
  T_s   : K -> K + s*I                          (pin shift inside the gapped class, s > -lambda_0)
  U     : K -> U K U^T, U orthogonal, U e_r = e_r (hidden-sector relabeling)

Checks:
  A. computable earned predicates (E-1, E-2, E-3, E-4-type, CM) are invariant under S_lam and T_s;
     E-2 locality is invariant under U only when the declared local frame is transported with U.
  B. mu_r is the complete U-invariant: moments agree to machine precision after a random hidden U.
  C. equivariance weights of quotient components (lambda_0, moments, bandwidth, kappa, edge exponent,
     standardized cumulants, lambda_0/W).
  D. weighted homogeneity of the record's dimensionful identities (E-14 kappa, E-15 Delta c_3).
  E. fixed-point lemma: an invariant predicate set can only select values fixed by rho(G).
"""
import numpy as np
import sympy as sp

rng = np.random.default_rng(1)
np.set_printoptions(precision=6, suppress=True)


def chain(n=400, pin=0.3, g=1.0):
    K = np.diag(np.full(n, 2 * g + pin)) - g * (np.eye(n, k=1) + np.eye(n, k=-1))
    return K, 0


def halfplane(L=30, pin=0.3):
    n = L * L
    K = np.zeros((n, n))
    idx = lambda i, j: i * L + j
    for i in range(L):
        for j in range(L):
            a = idx(i, j)
            K[a, a] = 4 + pin
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < L and 0 <= jj < L:
                    K[a, idx(ii, jj)] = -1.0
    return K, idx(0, L // 2)  # boundary site, middle of the edge


def ring_affinity(n=12, pin=0.3, d=0.4):
    """asymmetric ring (cycle affinity, E-4 type): forward 1+d, backward 1-d"""
    K = np.diag(np.full(n, 2 + pin))
    for i in range(n):
        K[i, (i + 1) % n] -= 1 + d
        K[i, (i - 1) % n] -= 1 - d
    return K, 0


def spectral_measure(K, r):
    w, V = np.linalg.eigh(K)
    return w, V[r, :] ** 2


def bandwidth(K, tol=1e-12):
    i, j = np.nonzero(np.abs(K) > tol)
    return int(np.max(np.abs(i - j)))


def predicates(K, r):
    sym = np.allclose(K, K.T, atol=1e-12)
    ev = np.linalg.eigvals(K)
    out = {
        "symmetric(passive/reciprocal)": sym,
        "gap>0 (E-1)": bool(np.min(ev.real) > 1e-12),
        "real spectrum": bool(np.max(np.abs(ev.imag)) < 1e-9),
        "local, declared frame (E-2)": bandwidth(K) <= 1 if K.shape[0] < 500 else None,
    }
    if sym:
        lam, wt = spectral_measure(K, r)
        out["CM kernel: weights>=0 & lam>=0"] = bool(np.all(wt > -1e-14) and np.all(lam > -1e-12))
    return out


def moments(lam, wt, kmax=4):
    return np.array([np.sum(wt * lam ** k) for k in range(1, kmax + 1)])


def standardized(lam, wt):
    m = np.sum(wt * lam)
    c2 = np.sum(wt * (lam - m) ** 2)
    c3 = np.sum(wt * (lam - m) ** 3)
    c4 = np.sum(wt * (lam - m) ** 4)
    return np.array([c3 / c2 ** 1.5, c4 / c2 ** 2 - 3])


def edge_exponent(lam, wt, frac=(0.002, 0.03)):
    """fit N(lambda_0 + e) ~ e^{gamma+1} for the cumulative retained weight near the bottom"""
    lam0 = lam[0]
    W = lam[-1] - lam[0]
    # merge (near-)degenerate eigenvalues: the retained weight of a degenerate eigenspace is basis-independent,
    # its split among eigh's eigenvectors is not (half-plane has exact degeneracies)
    key = np.round((lam - lam0) / W, 9)
    uk, inv = np.unique(key, return_inverse=True)
    wt = np.bincount(inv, weights=wt)
    e = uk * W
    cum = np.cumsum(wt)
    sel = (e > frac[0] * W) & (e < frac[1] * W)
    p = np.polyfit(np.log(e[sel]), np.log(cum[sel]), 1)
    return p[0] - 1


def components(K, r):
    lam, wt = spectral_measure(K, r)
    lam0, W = lam[0], lam[-1] - lam[0]
    m = moments(lam, wt)
    return {
        "lambda_0": lam0,
        "W (support width)": W,
        "m1=K_rr": m[0],
        "m2": m[1],
        "kappa=1/(2 sqrt m1)": 1 / (2 * np.sqrt(m[0])),
        "kappa^2*m1": (1 / (2 * np.sqrt(m[0]))) ** 2 * m[0],
        "lambda_0/W": lam0 / W,
        "edge exponent gamma": edge_exponent(lam, wt),
        "skew": standardized(lam, wt)[0],
        "excess kurtosis": standardized(lam, wt)[1],
    }


def report(name, K, r):
    print(f"\n=== {name} ===")
    base_p = predicates(K, r)
    base_c = components(K, r)
    ops = {
        "S_lam (lam=3)": (3.0 * K, r),
        "S_lam (lam=0.2)": (0.2 * K, r),
        "T_s (s=+0.5)": (K + 0.5 * np.eye(len(K)), r),
        "T_s (s=-0.2)": (K - 0.2 * np.eye(len(K)), r),
    }
    print("predicate invariance:")
    for op, (K2, r2) in ops.items():
        p2 = predicates(K2, r2)
        same = all(p2[k] == base_p[k] for k in base_p)
        print(f"  {op:18s} all predicates unchanged: {same}")
    print("component equivariance (ratio to base):")
    keys = list(base_c)
    print("  " + " | ".join(f"{k}" for k in keys))
    print("  base: " + " ".join(f"{base_c[k]:.6g}" for k in keys))
    for op, (K2, r2) in ops.items():
        c2 = components(K2, r2)
        print(f"  {op:18s}" + " ".join(f"{c2[k]:.6g}" for k in keys))
    return base_p, base_c


# ---------------- A + C: chain and half-plane ----------------
Kc, rc = chain()
report("pinned chain n=400, pin=0.3 (S5-1 / C1-a parent)", Kc, rc)
Kh, rh = halfplane()
report("pinned half-plane 30x30, boundary site (EDA W-A)", Kh, rh)

# E-4-type asymmetric ring: the predicate "not CM / complex modes" is S_lam and T_s invariant
Kr, rr = ring_affinity()
ev = np.linalg.eigvals(Kr)
print("\n=== asymmetric ring (E-4 type) ===")
for op, K2 in {"base": Kr, "S_lam 3": 3 * Kr, "T_s +0.5": Kr + 0.5 * np.eye(12)}.items():
    e2 = np.linalg.eigvals(K2)
    print(f"  {op:10s} complex modes present: {np.max(np.abs(e2.imag)) > 1e-9}  "
          f"max|Im|/max|Re-lam0| = {np.max(np.abs(e2.imag)) / np.max(np.abs(e2.real - e2.real.min())):.6f}")

# ---------------- B: hidden-sector unitary relabeling ----------------
n = 120
Kb, rb = chain(n=n)
Q, _ = np.linalg.qr(rng.standard_normal((n - 1, n - 1)))
U = np.eye(n)
U[1:, 1:] = Q
Ku = U @ Kb @ U.T
lam1, w1 = spectral_measure(Kb, rb)
lam2, w2 = spectral_measure(Ku, rb)
print("\n=== hidden-sector relabeling U (U e_r = e_r), chain n=120 ===")
print("  max |moment_k difference| k=1..8:",
      np.max(np.abs(moments(lam1, w1, 8) - moments(lam2, w2, 8)) / moments(lam1, w1, 8)))
print("  spectra equal:", np.allclose(lam1, lam2, atol=1e-10),
      " retained weights equal (sorted):", np.allclose(np.sort(w1), np.sort(w2), atol=1e-10))
print("  bandwidth in the ORIGINAL frame: before", bandwidth(Kb), " after", bandwidth(Ku))
print("  bandwidth in the TRANSPORTED frame (U^T Ku U):", bandwidth(U.T @ Ku @ U))
print("  => E-2 locality is U-invariant only with the frame transported; mu_r is U-invariant outright.")
# Lanczos from e_r recovers the same Jacobi matrix in both cases (mu_r is the complete invariant of (K, e_r))


def lanczos(K, r, m):
    q = np.zeros(len(K)); q[r] = 1.0
    Qs = [q]; a = []; b = []
    for j in range(m):
        v = K @ Qs[-1]
        aj = Qs[-1] @ v; a.append(aj)
        v = v - aj * Qs[-1] - (b[-1] * Qs[-2] if j > 0 else 0)
        for qq in Qs:  # full reorthogonalization
            v -= (qq @ v) * qq
        bj = np.linalg.norm(v); b.append(bj)
        Qs.append(v / bj)
    return np.array(a), np.array(b[:-1])


a1, b1 = lanczos(Kb, rb, 30)
a2, b2 = lanczos(Ku, rb, 30)
print("  Lanczos(K,e_r) vs Lanczos(UKU^T,e_r), first 30 coefficients: max diff",
      max(np.max(np.abs(a1 - a2)), np.max(np.abs(b1 - b2))))

# ---------------- D: weighted homogeneity of record identities ----------------
print("\n=== weighted homogeneity of record identities (sympy) ===")
lam_, t = sp.symbols("lam t", positive=True)
wb, wa, wT, wx = sp.symbols("w_beta w_a w_T1 w_x")
# rate weight 1: K11 -> lam*K11, t -> t/lam
# E-15: Delta c3 = 24 beta T1 a (44 beta a^2 + 5 K11): needs w(beta a^2) = w(K11) = 1
# E-15: m1^S - m1^D = -12 beta T1 a t^2 + O(t^3): w(x) = w_beta + w_T1 + w_a - 2
sol = sp.solve([sp.Eq(wb + 2 * wa, 1)], [wb], dict=True)
print("  E-15 Delta c3 homogeneous iff w_beta + 2 w_a = 1:", sol)
print("  E-15 drift identity sets w_x = w_beta + w_T1 + w_a - 2 (consistent for any w_T1); pure numbers 12, 24, 44, 5 are weight-0")
K11 = sp.symbols("K11", positive=True)
kappa = 1 / (2 * sp.sqrt(K11))
print("  E-14 kappa(lam K11)/kappa(K11) =", sp.simplify(kappa.subs(K11, lam_ * K11) / kappa), " (weight -1/2; kappa^2 K11 = 1/4 is weight 0)")

# ---------------- E: fixed-point lemma, explicit ----------------
print("\n=== fixed-point lemma ===")
print("  If P is G-invariant and P => Q = q* with Q(g m) = rho(g) Q(m), then rho(g) q* = q* for all admissible g.")
for wgt in (1, -sp.Rational(1, 2), 2):
    q = sp.symbols("q", nonnegative=True)
    fix = sp.solve(sp.Eq(lam_ ** wgt * q, q), q)
    print(f"  weight {wgt}: solutions of lam^w q = q for generic lam:", fix, "(plus q = infinity)")
print("  pin shift: lambda_0 + s = lambda_0 has no solution for s != 0 -> no finite lambda_0 selectable by T_s-invariant P")
print("  => invariant P can only impose G-invariant sets: for a weight-w quantity {0}, (0,inf), {inf} and unions")
print("     (sign/class statements), or exact values of weight-0 (G-invariant) combinations.")
