"""RA0 / G4 -- asymptotic metastable algebra.  NUMERICAL ILLUSTRATION; independent code path, not independent reviewer.

Firewall by code structure (G4.8):
  * `blind(lam, U, pi, rank)` sees ONLY the spectral data of P_N / G_N (eigenvalues, pi-orthonormal eigenfunctions, the
    stationary law pi derived from P_N) and a rank.  It never sees labels, coordinates, K, or a threshold.
  * The rank is a cut rank whose divergence is CERTIFIED analytically (Weyl lower bound, Prop G4-T) or REFUSED
    analytically (Props G3-A / G3-B; ladder beta <= 2).  Finite-size ratio trends are printed but are not the certificate.
  * `audit(...)` and `thm_T(...)` use the hidden construction labels, and are called only AFTER `blind` has returned.
Norm: L^2(pi_N).  pi_N is the unique stationary law of the irreducible reversible P_N (priced: reversibility,
irreducibility).  Spectral projectors are pi-orthogonal because P_N is pi-self-adjoint.  Nothing here is measure-free.
Computational tolerances (numerical-equality only, not modelling choices): zero-mode detection 1e-10 relative,
power-iteration convergence 1e-13, de-duplication of converged fixed points 1e-6.
"""
import itertools
import numpy as np

np.set_printoptions(precision=4, suppress=True)
RNG_COMPUTE = 12345           # seeds only the starting points of numerical searches, never a family


# ------------------------------------------------------------------ generic reversible machinery
def from_conductances(C, nu):
    """Continuous-time reversible chain: rate x->y = C[x,y]/nu[x];  pi ∝ nu.  Returns (Lsym, pi)."""
    C = np.array(C, float); np.fill_diagonal(C, 0.0)
    s = np.sqrt(nu)
    Lsym = -C / np.outer(s, s)
    Lsym[np.diag_indices_from(Lsym)] = C.sum(1) / nu
    return Lsym, nu / nu.sum()


def from_kernel(P, pi):
    """Discrete-time reversible kernel: relaxation operator I - P, symmetrised in L^2(pi)."""
    s = np.sqrt(pi)
    return np.eye(len(pi)) - (s[:, None] * P) / s[None, :], pi


def spectrum(Lsym, pi):
    lam, phi = np.linalg.eigh((Lsym + Lsym.T) / 2)
    U = phi / np.sqrt(pi)[:, None]                 # pi-orthonormal eigenfunctions:  sum_x pi u_i u_j = delta_ij
    U[:, 0] = 1.0                                  # irreducible: zero mode is the constant (fix sign/rounding)
    return np.maximum(lam, 0.0), U


def top_cuts(lam, top=3):
    z = int(np.sum(lam <= 1e-10 * lam.max()))
    r = lam[z + 1:] / lam[z:-1]                    # ratio at rank k = lambda_k / lambda_{k-1}, k = z+1 ...
    order = np.argsort(r)[::-1][:top]
    return [(float(r[i]), int(z + 1 + i)) for i in order]


def ratio_at(lam, k):
    return lam[k] / lam[k - 1]


# ------------------------------------------------------------------ the blind procedure (no labels, no K, no epsilon)
def pin(pi, f, g):
    return float(np.sum(pi * f * g))


def defects(V, pi):
    """V: pi-orthonormal basis (n x k).  Returns (HS defect, sup defect estimate).  Both basis-independent.
    HS^2 = sum_{i,j} ||(I-E)(u_i u_j)||^2  (Hilbert-Schmidt norm of the product-residual map, an upper bound for sup).
    sup  = max over unit f,g in V of ||(I-E)(fg)||  (alternating maximisation; a lower estimate of the true sup)."""
    n, k = V.shape
    w = np.sqrt(pi)
    R = np.zeros((k, k, n))
    for i in range(k):
        for j in range(k):
            h = V[:, i] * V[:, j]
            R[i, j] = (h - V @ (V.T @ (pi * h))) * w          # residual, scaled so Euclidean norm = pi-norm
    hs = float(np.sqrt(np.sum(R ** 2)))
    rng = np.random.default_rng(RNG_COMPUTE); best = 0.0
    for _ in range(30):
        a = rng.standard_normal(k); a /= np.linalg.norm(a)
        for _ in range(200):
            M = np.einsum('i,ijn->nj', a, R); _, s_, vt = np.linalg.svd(M, full_matrices=False); b = vt[0]
            M = np.einsum('j,ijn->ni', b, R); _, s_, vt = np.linalg.svd(M, full_matrices=False); a = vt[0]
        best = max(best, float(s_[0]))
    return hs, best


def christoffel(V):
    """C_V = sup_{f in V, ||f||_pi = 1} ||f||_inf = sqrt(max_x sum_i u_i(x)^2)   (basis-free)."""
    return float(np.sqrt(np.max(np.sum(V ** 2, 1))))


def robust_idempotents(V, pi, starts=300):
    """Attracting fixed points of the tensor power map v -> tau(v,v)/|tau(v,v)|, tau(a,b,c) = pi(f_a f_b f_c) on V.
    A fixed point tau(v,v) = lam v gives the idempotent e = (V v)/lam of the compressed product  f*g = E(fg):
    e*e = e.  Attracting fixed points <-> primitive idempotents for an exact partition algebra (odeco tensor)."""
    n, k = V.shape
    rng = np.random.default_rng(RNG_COMPUTE)
    found, nonconv = [], 0
    for _ in range(starts):
        v = rng.standard_normal(k); v /= np.linalg.norm(v)
        ok = False
        for _ in range(5000):
            f = V @ v
            t = V.T @ (pi * f * f)
            nt = np.linalg.norm(t)
            if nt == 0:
                break
            v2 = t / nt
            if np.linalg.norm(v2 - v) < 1e-13:
                v = v2; ok = True; break
            v = v2
        if not ok:
            nonconv += 1; continue
        f = V @ v; lam = float(np.sum(pi * f ** 3))
        if lam <= 0:
            continue
        T = V.T @ ((pi * f)[:, None] * V)
        Pp = np.eye(k) - np.outer(v, v)
        rho = float(np.max(np.abs(np.linalg.eigvals(Pp @ (2 * T / lam) @ Pp))))
        if rho < 1 and not any(np.linalg.norm(v - u) < 1e-6 for u, _ in found):
            found.append((v, lam))
    E = np.array([(V @ v) / lam for v, lam in found]).T if found else np.zeros((n, 0))
    return E, nonconv


def partition_basis(lab, pi):
    blocks = np.unique(lab)
    A = np.column_stack([(lab == b).astype(float) / np.sqrt(pi[lab == b].sum()) for b in blocks])
    return A                                             # pi-orthonormal basis of the partition algebra


def gap_metric(V, A, pi):
    """||E_V - E_A|| in L^2(pi) (basis-free).  Equal dims: sin of the largest principal angle; else 1."""
    if V.shape[1] != A.shape[1]:
        return 1.0
    R = (V - A @ (A.T @ (pi[:, None] * V))) * np.sqrt(pi)[:, None]     # (I - E_A) on a pi-orthonormal basis of V
    return float(np.linalg.norm(R, 2))                               # = ||(I-E_A)E_V|| = ||E_V - E_A|| (equal dims)


def blind(lam, U, pi, rank):
    """G4.5.  Input: spectral data of P_N only, and a certified cut rank.  Output: recovered partition + diagnostics."""
    V = U[:, :rank]
    hs, sup = defects(V, pi)
    E, nonconv = robust_idempotents(V, pi)
    m = E.shape[1]
    out = dict(rank=rank, hs=hs, sup=sup, CV=christoffel(V), m=m, nonconv=nonconv)
    if m == 0:
        out.update(lab=None, gap=1.0); return out
    lab = np.argmax(E, axis=1)                                   # canonical rounding: argmax over idempotents
    A = partition_basis(lab, pi)
    out.update(lab=lab, nblocks=A.shape[1], gap=gap_metric(V, A, pi),
               unit=float(np.sqrt(pin(pi, E.sum(1) - 1, E.sum(1) - 1))),
               idem=float(max(np.sqrt(pin(pi, E[:, i] ** 2 - E[:, i], E[:, i] ** 2 - E[:, i])) for i in range(m))))
    if rank == 2:                                                 # Prop G4-P (rank 2): canonical two-valued rounding
        f = V[:, 1]; s3 = float(np.sum(pi * f ** 3))
        labP = (f > s3 / 2).astype(int)
        out['pearson_bound'] = float(np.sqrt(max(0.0, np.sum(pi * f ** 4) - 1 - s3 ** 2)))
        out['pearson_gap'] = gap_metric(V, partition_basis(labP, pi), pi) if len(np.unique(labP)) == 2 else 1.0
        out['pearson_agrees_with_argmax'] = bool(len(np.unique(lab)) == 2 and
                                                 (np.all(labP == lab) or np.all(labP == 1 - lab)))
    return out


# ------------------------------------------------------------------ post-hoc diagnostics (hidden labels; AFTER blind)
def audit(lab_rec, lab_true, pi):
    if lab_rec is None:
        return dict(mis=1.0, exact=False)
    mis = 0.0
    for b in np.unique(lab_rec):
        sel = lab_rec == b
        vals, counts = np.unique(lab_true[sel], return_counts=True)
        maj = vals[np.argmax([pi[sel & (lab_true == v)].sum() for v in vals])]
        mis += pi[sel & (lab_true != maj)].sum()
    exact = len(np.unique(lab_rec)) == len(np.unique(lab_true)) and mis == 0.0
    return dict(mis=float(mis), exact=bool(exact))


def thm_T(C, nu, lab_true, pi, V):
    """Prop G4-T quantities (proof side; uses hidden labels).  g = min intra-block spectral gap of G0, eta = ||G1||_pi."""
    C = np.array(C, float); same = lab_true[:, None] == lab_true[None, :]
    L0, _ = from_conductances(np.where(same, C, 0.0), nu)
    L1 = from_conductances(C, nu)[0] - L0
    ev0 = np.linalg.eigvalsh((L0 + L0.T) / 2); K = len(np.unique(lab_true))
    g = float(ev0[K]); eta = float(np.max(np.abs(np.linalg.eigvalsh((L1 + L1.T) / 2))))
    A = partition_basis(lab_true, pi)
    pmin = min(pi[lab_true == b].sum() for b in np.unique(lab_true))
    s = gap_metric(V, A, pi)
    return dict(g=g, eta=eta, weyl=(g - eta) / eta if eta > 0 else np.inf,
                dk=eta / (g - eta) if g > eta else np.inf, s_true=s, pmin=pmin,
                dbound=s * (4 / np.sqrt(pmin) + christoffel(V)))


# ------------------------------------------------------------------ families (deterministic rules across N)
def h(x, y):                                       # deterministic positive symmetric modulation, no symmetry group
    return 1.0 + 0.5 * np.cos(0.37 * (x + y)) + 0.3 * np.cos(1.1 * np.abs(x - y))


def nu_rule(n):
    return 1.0 + 0.5 * np.cos(1.7 * np.arange(n))


def F1(M):
    """G4.4: K=3 sectors of non-identical sizes M, floor(1.5M), 2M; intra conductance h/n_s; inter h*kappa_st*M^-2."""
    sizes = [M, int(1.5 * M), 2 * M]; lab = np.repeat(np.arange(3), sizes); n = lab.size
    kap = np.array([[0, 1.0, 0.5], [1.0, 0, 2.0], [0.5, 2.0, 0]])
    x = np.arange(n); H = h(x[:, None], x[None, :])
    ns = np.array(sizes)[lab]
    C = np.where(lab[:, None] == lab[None, :], H / ns[:, None], H * kap[lab][:, lab] / M ** 2)
    return C, nu_rule(n), lab


def F2(M):
    """Nested: 4 basins (sizes M, 1.25M, 1.5M, 1.75M) in super-basins {0,1},{2,3}; within-super inter ~M^-2, across ~M^-3."""
    sizes = [M, int(1.25 * M), int(1.5 * M), int(1.75 * M)]; lab = np.repeat(np.arange(4), sizes); n = lab.size
    sup = lab // 2; x = np.arange(n); H = h(x[:, None], x[None, :]); ns = np.array(sizes)[lab]
    C = np.where(lab[:, None] == lab[None, :], H / ns[:, None],
                 np.where(sup[:, None] == sup[None, :], H / M ** 2, H / M ** 3))
    return C, nu_rule(n), lab, sup


def ladder(L, beta, hop=(1.0, 1.0), mod=False):
    """mod=False: G3 ladder (translation-invariant: lane indicators are EXACTLY lumpable for every beta).
    mod=True: site-dependent hops hop_lane*(1+0.3cos(1.1x)) and rungs gamma*(1+0.5cos(0.37x)); no exact lumpability."""
    n = 2 * L; C = np.zeros((n, n)); g = L ** (-float(beta))
    for lane in range(2):
        for xx in range(L):
            i, j = lane * L + xx, lane * L + (xx + 1) % L
            C[i, j] = C[j, i] = hop[lane] * ((1 + 0.3 * np.cos(1.1 * xx)) if mod else 1.0)
            C[i, (1 - lane) * L + xx] = C[(1 - lane) * L + xx, i] = g * ((1 + 0.5 * np.cos(0.37 * xx)) if mod else 1.0)
    return C, np.ones(n), np.repeat([0, 1], L)


def ssep(L):
    N = L // 2; confs = list(itertools.combinations(range(L), N)); idx = {c: i for i, c in enumerate(confs)}
    P = np.zeros((len(confs), len(confs)))
    for i, c in enumerate(confs):
        occ = set(c)
        for b in range(L):
            u, v = b, (b + 1) % L
            if (u in occ) != (v in occ):
                new = set(occ); new.symmetric_difference_update({u, v}); P[i, idx[tuple(sorted(new))]] += 1.0 / L
            else:
                P[i, i] += 1.0 / L
    return P, np.full(len(confs), 1.0 / len(confs))


def indep(N, a=0.1, b=0.25):
    p = np.array([[1 - a, a], [b, 1 - b]]); q = np.array([b, a]) / (a + b); P, pi = p, q
    for _ in range(N - 1):
        P = np.kron(P, p); pi = np.kron(pi, q)
    return P, pi


# ------------------------------------------------------------------ runs
def fmt(o):
    s = (f"rank {o['rank']}: Delta_HS {o['hs']:.3e}  Delta_sup~{o['sup']:.3e}  C_V {o['CV']:.2f}  "
         f"robust idempotents found {o['m']}")
    if o['lab'] is not None:
        s += f" -> {o['nblocks']} blocks;  gap(V, A_blind) {o['gap']:.3e};  |sum e - 1| {o['unit']:.1e};  max|e^2-e| {o['idem']:.1e}"
    return s


def run_metastable(name, sizes, build, certified_ranks, labels_for_rank):
    print(f"== {name} ==")
    for M in sizes:
        C, nu, *labs = build(M)
        Lsym, pi = from_conductances(C, nu); lam, U = spectrum(Lsym, pi)
        print(f"  M={M:3d} (n={len(pi)}): top cuts (ratio, rank) {[(round(r, 2), k) for r, k in top_cuts(lam)]}")
        recs = {}
        for k in certified_ranks:
            o = blind(lam, U, pi, k); recs[k] = o
            print(f"    BLIND  {fmt(o)}")
            if 'pearson_bound' in o:
                print(f"    rank-2 Prop G4-P: Pearson bound sqrt(kurt-1-skew^2) {o['pearson_bound']:.3e} >= gap(V, A_f) "
                      f"{o['pearson_gap']:.3e};  agrees with argmax partition: {o['pearson_agrees_with_argmax']}")
        for k in certified_ranks:                                   # ---- post-hoc: labels enter only here
            lt = labels_for_rank(k, labs); a = audit(recs[k]['lab'], lt, pi); t = thm_T(C, nu, lt, pi, U[:, :k])
            print(f"    AUDIT  rank {k}: misassigned pi-mass {a['mis']:.2e}, exact match {a['exact']};  "
                  f"Thm T: eta/g {t['eta'] / t['g']:.2e}, Weyl ratio bound {t['weyl']:.2f} (actual {ratio_at(lam, k):.2f}), "
                  f"DK bound {t['dk']:.2e} >= gap(V, A_true) {t['s_true']:.2e};  Delta bound {t['dbound']:.2e} >= Delta_sup")
        if len(certified_ranks) == 2:
            k1, k2 = certified_ranks; l1, l2 = recs[k1]['lab'], recs[k2]['lab']
            viol = sum(pi[(l2 == b) & (l1 != np.bincount(l1[l2 == b]).argmax())].sum() for b in np.unique(l2))
            print(f"    NESTING (blind vs blind, no labels): pi-mass of fine blocks outside their coarse block {viol:.2e};"
                  f" exact refinement {viol == 0}")


def run_refusal(name, items):
    print(f"== {name} ==")
    for tag, (Lsym, pi) in items:
        lam, U = spectrum(Lsym, pi); tc = top_cuts(lam)
        print(f"  {tag} (n={len(pi)}): top cuts {[(round(r, 3), k) for r, k in tc]}  -> no certified diverging cut: REFUSE")
        o = blind(lam, U, pi, tc[0][1])                              # STRESS ONLY: force the largest finite-N cut
        print(f"    stress (forced rank {tc[0][1]}, not a G4 output): {fmt(o)}")


if __name__ == "__main__":
    run_metastable("F1  deterministic non-symmetric metastable family (3 sectors)", [5, 10, 20, 40, 80], F1, [3],
                   lambda k, labs: labs[0])
    run_metastable("F2  nested deterministic family (4 basins / 2 super-basins)", [5, 10, 20, 40, 80],
                   F2, [2, 4], lambda k, labs: labs[1] if k == 2 else labs[0])
    run_metastable("LADDER beta=3, identical lanes", [25, 50, 100, 200], lambda L: ladder(L, 3), [2],
                   lambda k, labs: labs[0])
    run_metastable("LADDER beta=3, non-identical, site-modulated lanes (hop 1.0 / 1.6; not lumpable)", [25, 50, 100, 200],
                   lambda L: ladder(L, 3, (1.0, 1.6), True), [2], lambda k, labs: labs[0])
    run_refusal("FALSE-POSITIVE A: SSEP ring, N=L/2", [(f"L={L}", from_kernel(*ssep(L))) for L in (8, 10, 12)])
    run_refusal("FALSE-POSITIVE B: independent two-state particles",
                [(f"N={N}", from_kernel(*indep(N))) for N in (2, 4, 6, 8)])
    for beta in (0, 1, 2):
        run_refusal(f"FALSE-POSITIVE C: ladder beta={beta} (non-identical, site-modulated lanes)",
                    [(f"L={L}", from_conductances(*ladder(L, beta, (1.0, 1.6), True)[:2])) for L in (25, 50, 100, 200)])
