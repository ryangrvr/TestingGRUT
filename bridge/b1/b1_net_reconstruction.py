"""GRUT BRIDGE-1 / B1 — can the canonical GRUT Level-0 local net be reconstructed with the net DELETED from the input?
Canonical models (READ-ONLY source GRUT-RAI master-w25bu9 @ b935099; re-implemented here, not imported):
  C1-a  : K = 0.3 I + L(path, unit springs), N = 24; retained site 0 + bath   (L0_1A_CHARTER_01.md)
  L0-1b : K = 0.3 I + L(w), w_ij = c_a |i-j|^-a                               (L0_1B_CHARTER_01.md)
  L0-1c : V = 1/2 x^T K_b x + beta sum_i x_i^4 (on-site convex quartic)        (L0_1C_CHARTER_01.md)
  L0-1e : additive noise Q = 2 diag(T_i)                                      (L0_ACCESS_BRIDGE_01.md)
Net deleted = the generator is given only as an abstract operator, i.e. in an unknown orthonormal frame W (K' = W K W^T).
Firewall: no reconstruction criterion may name the canonical sites, the path graph, or the site-resolved access."""
import numpy as np
from b1_lib import *
N = 24; K = 0.3 * np.eye(N) + lap(path_w(N)); W = haar_orth(N); Kp = W @ K @ W.T     # Kp: generator with the net deleted
ev = np.linalg.eigvalsh(K)
hdr("B1-0 canonical net, for reference only (NOT an input to any reconstruction below)")
print(f"  C1-a: N={N}, spectrum [{ev.min():.4f}, {ev.max():.4f}], nondegenerate: {np.min(np.diff(ev)) > 1e-9}; nnz(K) = {np.count_nonzero(np.abs(K) > 1e-12)}")
print(f"  canonical geometry predicate (L0-1b): Q, monotone, SURVIVES = {geometry_predicate(path_w(N))}")

hdr("B1-1 interaction support / sparsity (generator-only criteria on Kp)")
evals, U = np.linalg.eigh(Kp)
print(f"  sparsest frame: the normal-mode frame makes Kp DIAGONAL (nnz = {N}) < canonical chain (nnz = {3*N-2}).")
print("  -> pure sparsity / autonomy / response-factorization optimum = decoupled normal modes: NO interaction edges, the")
print("     L0-1b geometry predicate is undefined (empty graph). The canonical net is NOT the sparsity optimum.")
print("  Requiring a CONNECTED sparsest net (bandwidth-1 chain): Lanczos from ANY unit start vector q gives a frame where Kp is")
print("  tridiagonal (a nearest-neighbour chain), so the constraint selects a CONTINUUM (q on S^{N-1} modulo sign).")
stats = dict(n=0, coop=0, pins_ok=0, geom=0, inequiv=0, uniform_pin=0)
examples = []
for trial in range(400):
    q = rng.normal(size=N); out = lanczos(Kp, q)
    if out is None: continue
    Qm, T = out; stats["n"] += 1
    w, pins, coop = as_spring_network(T); stats["coop"] += coop
    ok_pins = bool((pins >= -1e-9).all()); stats["pins_ok"] += ok_pins
    Qg, mono, surv = geometry_predicate(w); stats["geom"] += surv and ok_pins
    ineq = not signed_perm_equiv(T, K); stats["inequiv"] += ineq and surv and ok_pins
    stats["uniform_pin"] += bool(np.allclose(pins, 0.3, atol=1e-6))
    if surv and ok_pins and ineq and len(examples) < 3: examples.append((np.round(np.diag(T)[:5], 3), np.round(w[np.arange(4), np.arange(1, 5)], 3), np.round(pins[:5], 3), round(Qg, 1)))
print(f"  400 random start vectors: {stats['n']} chains; sign-gauge cooperative (Kamke) {stats['coop']}; all pins >= 0 (passive spring reading) {stats['pins_ok']};"
      f" L0-1b geometry predicate SURVIVES with pins >= 0: {stats['geom']}; of these NOT signed-permutation-equivalent to the canonical chain: {stats['inequiv']};"
      f" uniform pin 0.3 (canonical C1-a form): {stats['uniform_pin']}")
for e in examples: print(f"    e.g. diag {e[0]}  springs {e[1]}  pins {e[2]}  geometry Q = {e[3]}")
# existence of passive chains NEAR the canonical one (pins strictly positive there, so an open set of start vectors works)
from scipy.linalg import expm
from scipy.optimize import minimize
e0 = np.zeros(N); e0[0] = 1
for eps in (0.02, 0.05, 0.1, 0.2):
    good = 0; ineq = 0; geo = 0; tries = 50
    for _ in range(tries):
        qq = W @ (e0 + eps * rng.normal(size=N)); out = lanczos(Kp, qq)
        if out is None: continue
        w, pins, coop = as_spring_network(out[1]); okp = bool((pins >= -1e-9).all()); Qg, mono, surv = geometry_predicate(w)
        good += okp; geo += okp and surv; ineq += okp and surv and not signed_perm_equiv(out[1], K)
    print(f"  start q = canonical end-site + eps*noise, eps={eps}: passive (pins >= 0) chains {good}/{tries}; + geometry SURVIVES {geo}; of these inequivalent to canonical {ineq}")
print("  -> GENERAL PASSIVE COOPERATIVE CLASS (weighted springs, non-uniform pins >= 0): a CONTINUUM of inequivalent chain nets for")
print("     the same generator, all passing the L0-1b geometry predicate. Random q fails (pins < 0) only because the admissible set is")
print("     a neighbourhood, not the whole sphere.")
# L0-1b class proper: uniform pin 0.3 + arbitrary non-negative weights (K - 0.3 I must be a weighted graph Laplacian)
L0 = K - 0.3 * np.eye(N); one = np.ones(N) / np.sqrt(N); P1 = np.eye(N) - np.outer(one, one)
cnt = 0
for _ in range(3000):
    A = rng.normal(size=(N, N)); A = A - A.T; A = P1 @ A @ P1
    Wl = expm(0.05 * A); Lw = Wl @ L0 @ Wl.T
    if (Lw[~np.eye(N, dtype=bool)] <= 1e-12).all(): cnt += 1
print(f"  L0-1b class (uniform pin 0.3, weighted springs >= 0): random small rotations fixing the uniform mode that keep K - 0.3I a weighted Laplacian: {cnt}/3000")
pairs = [(i, j) for i in range(N) for j in range(i + 1, min(N, i + 3))]
def build(a):
    A = np.zeros((N, N))
    for k, (i, j) in enumerate(pairs): A[i, j] = a[k]; A[j, i] = -a[k]
    A = P1 @ A @ P1; E = expm(A); return E @ L0 @ E.T
def viol(a):
    Lw = build(a); off = Lw[~np.eye(N, dtype=bool)]
    return np.sum(np.clip(off, 0, None) ** 2) + 1e-2 * max(0.0, 0.5 - np.linalg.norm(a)) ** 2
best = None
for trial in range(8):
    r = minimize(viol, rng.normal(size=len(pairs)) * 0.3, method="L-BFGS-B")
    Lw = build(r.x)
    if viol(r.x) < 1e-14 and not signed_perm_equiv(Lw + 0.3 * np.eye(N), K):
        w_ = -Lw.copy(); np.fill_diagonal(w_, 0); w_ = np.where(w_ > 1e-12, w_, 0)
        best = (geometry_predicate(w_), np.count_nonzero(w_) // 2, np.abs(np.linalg.eigvalsh(Lw) - np.linalg.eigvalsh(L0)).max()); break
if best:
    print(f"  optimization: an INEQUIVALENT weighted Laplacian, isospectral to the canonical one (max eigenvalue diff {best[2]:.1e}), exists in the"
          f" L0-1b class: edges {best[1]}, geometry predicate Q={best[0][0]:.1f}, monotone={best[0][1]}, SURVIVES={best[0][2]}")
else:
    print("  optimization: no inequivalent isospectral weighted Laplacian found in 8 restarts -> reported as NOT FOUND (the L0-1b class may pin the net)")
print("  NARROW C1-a class (uniform pin 0.3, UNIT springs): the path is determined by its Laplacian spectrum (DLS, known result), so the")
print("  net is unique up to graph automorphism and the K-commutant (H-relative gauge) -> CONDITIONAL on the supplied class (CRITERION-PRICED).")

hdr("B1-2 the same generator spectrum, a 2D canonical-style net vs a 1D chain net (geometry / dimension not reconstructable)")
Lg = 5; w2 = np.zeros((Lg * Lg, Lg * Lg)); w2u = np.zeros_like(w2)
for x in range(Lg):
    for y in range(Lg):
        i = x * Lg + y
        if x + 1 < Lg: w2u[i, i + Lg] = w2u[i + Lg, i] = 1; w2[i, i + Lg] = w2[i + Lg, i] = rng.uniform(0.8, 1.2)
        if y + 1 < Lg: w2u[i, i + 1] = w2u[i + 1, i] = 1; w2[i, i + 1] = w2[i + 1, i] = rng.uniform(0.8, 1.2)
ev_u = np.linalg.eigvalsh(lap(w2u)); deg = int(np.sum(np.diff(np.round(ev_u, 9)) == 0))
print(f"  UNIFORM 5x5 grid: {deg} degenerate eigenvalue pairs -> it CANNOT be a connected chain in any frame (a Jacobi matrix has a simple")
print("  spectrum): spectral multiplicity constrains which nets are possible (a partial, not complete, constraint).")
K2 = 0.3 * np.eye(Lg * Lg) + lap(w2); out = lanczos(K2, rng.normal(size=Lg * Lg)); T2 = out[1]
wch, pch, _ = as_spring_network(T2)
print(f"  5x5 grid net: geometry predicate {geometry_predicate(w2)[:2]} (Q, line-ordering) ; Lanczos chain frame of the SAME operator:"
      f" {geometry_predicate(wch)[:2]}; chain pins >= 0: {bool((pch >= -1e-9).all())}")
print("  -> the generator alone does not fix even the DIMENSION of the net-derived geometry (2D lattice vs 1D chain, same spectrum).")

hdr("B1-3 invariant-sector structure and response kernels")
print(f"  Kp has nondegenerate spectrum: its invariant subspaces are spans of eigenvectors (global modes), none localized -> no net.")
u = W[:, 0]; tau = np.array([1, 5, 10, 20, 40.0])
kern = lambda v: [float(v @ (np.linalg.eigh(Kp)[1] @ np.diag(np.exp(-np.linalg.eigh(Kp)[0] * t)) @ np.linalg.eigh(Kp)[1].T) @ v) for t in tau]
q = rng.normal(size=N); q /= np.linalg.norm(q)
print(f"  memory kernel k(tau) of the canonical retained site:  " + " ".join(f"{x:.2e}" for x in kern(u)))
print(f"  memory kernel of a Lanczos-chain 'site 0' (random q):  " + " ".join(f"{x:.2e}" for x in kern(q)))
print("  -> every unit vector is the 'retained site' of some chain net; the retained/bath partition is not fixed by K.")

hdr("B1-4 noise layer (L0-1e: Q = 2 diag(T_i)) as a net carrier")
for lab, T_ in (("uniform T_i = 1", np.ones(N)), ("non-uniform T_i", rng.uniform(0.5, 2.0, N))):
    Qn = W @ (2 * np.diag(T_)) @ W.T
    ev_q = np.linalg.eigvalsh(Qn); degenerate = np.max(ev_q) - np.min(ev_q) < 1e-9
    if degenerate:
        print(f"  {lab}: Q' is proportional to I in every frame -> carries NO net information")
    else:
        _, Vq = np.linalg.eigh(Qn); overlap = np.abs(Vq.T @ W).max(0)
        print(f"  {lab}: the eigenframe of Q' recovers the canonical site axes (min |overlap| = {overlap.min():.6f}) -> net carried by the NOISE layer")
print("  -> with non-uniform site temperatures the net is RECONSTRUCTABLE from Q alone: RELOCATION (net information encoded in the")
print("     supplied environment layer), not derivation.")

hdr("B1-5 nonlinear drift (L0-1c on-site convex quartic) as a net carrier")
beta = 0.7
T4 = np.zeros((N, N, N, N))
for i in range(N): v = W[:, i]; T4 += 24 * beta * np.einsum('a,b,c,d->abcd', v, v, v, v)   # 4th derivative tensor of beta sum (W^T y)_i^4
def power_iter(T, iters=200):
    x = rng.normal(size=N); x /= np.linalg.norm(x)
    for _ in range(iters): x = np.einsum('abcd,b,c,d->a', T, x, x, x); x /= np.linalg.norm(x)
    lam = np.einsum('abcd,a,b,c,d->', T, x, x, x, x); return x, lam
Td = T4.copy(); found = []
for _ in range(N):
    x, lam = power_iter(Td); found.append(x); Td = Td - lam * np.einsum('a,b,c,d->abcd', x, x, x, x)
F = np.array(found).T; ov = np.abs(F.T @ W); match = ov.max(1)
print(f"  quartic tensor of V in the deleted-net frame -> tensor power iteration + deflation recovers {N} axes;"
      f" min over recovered axes of max |overlap with a canonical site| = {match.min():.8f}; residual ||T|| after deflation = {np.linalg.norm(Td):.1e}")
print("  orthogonally decomposable symmetric tensors have a unique decomposition up to sign/permutation (known result) -> the")
print("  on-site quartic drift DETERMINES the canonical site axes. The drift was supplied in site-local form, so this is")
print("  RELOCATION / REDUNDANT SUPPLY (net recoverable from the drift layer), not derivation from un-netted structure.")
Wr = haar_orth(N); Tg = np.zeros((N, N, N, N))
for i in range(N): v = Wr[:, i]; Tg += 24 * beta * np.einsum('a,b,c,d->abcd', v, v, v, v)
print("  control: a quartic that is on-site in a DIFFERENT frame recovers THAT frame instead (overlap with canonical sites "
      f"{np.abs(Wr.T @ W).max():.3f}) -> the drift selects whatever net it was written in.")
