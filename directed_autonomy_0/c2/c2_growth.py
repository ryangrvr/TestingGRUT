"""DA0 / C2 -- growing endogenous architecture.  NUMERICAL ILLUSTRATION; independent code path, not independent reviewer.

Firewall (C2-F1..F7): nothing is optimised; outputs are structural (K_N, cut ranks, filtration, nesting defect, macro
currents).  BLIND side: `blind_level(G or Lsym, pi, rank)` -> slow space V (Riesz / spectral), canonical PRIMITIVE
idempotents of the compressed product (C2-A3: attracting fixed points of v -> T(v,v)/|T(v,v)|, kept iff exactly one
eigenvalue of T(c,.,.) exceeds 1/2), argmax partition.  Random starts are a numerical solver only; completeness is
evidenced (count == dim V, sum e ~ 1), not proved, unless Lemma G5-3' applies.  PROOF / AUDIT side (hidden labels):
s = gap(V, A_true), Lambda, kappa = K p_max, C_V, epsilon bounds, m_tot, m_rel (Hungarian matching -> an upper bound on
the min-max m_rel), nesting defect.
Computational tolerances (numerical equality only): power-iteration convergence 1e-12, de-duplication 1e-6 (relative).
"""
import os, sys, itertools
import numpy as np
from scipy.optimize import linear_sum_assignment
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'reflexive_autonomy_0', 'g4'))
from g4_algebra import h, nu_rule, gap_metric, partition_basis, christoffel, defects

np.set_printoptions(precision=4, suppress=True, linewidth=170)
SEED = 20261003


# ------------------------------------------------------------------ reversible machinery
def from_conductances(C, nu):
    C = np.array(C, float); np.fill_diagonal(C, 0.0); s = np.sqrt(nu)
    L = -C / np.outer(s, s); L[np.diag_indices_from(L)] = C.sum(1) / nu
    return L, nu / nu.sum()


def spectrum(Lsym, pi):
    lam, phi = np.linalg.eigh((Lsym + Lsym.T) / 2)
    U = phi / np.sqrt(pi)[:, None]; U[:, 0] *= np.sign(U[0, 0]) or 1.0
    return np.maximum(lam, 0.0), U


def cut_ratios(lam):
    z = int(np.sum(lam <= 1e-12 * lam.max())); r = lam[z + 1:] / np.maximum(lam[z:-1], 1e-300)
    return {k: float(r[k - z - 1]) for k in range(z + 1, len(lam))}


# ------------------------------------------------------------------ BLIND: canonical primitive idempotents (C2-A3)
def primitive_idempotents(V, pi, starts_per_dim=24, iters=300):
    n, k = V.shape; rng = np.random.default_rng(SEED); found = []
    for _ in range(starts_per_dim * k):
        v = rng.standard_normal(k); v /= np.linalg.norm(v)
        for _ in range(iters):
            f = V @ v; t = V.T @ (pi * f * f); nt = np.linalg.norm(t)
            if nt == 0:
                break
            v2 = t / nt
            if np.linalg.norm(v2 - v) < 1e-12:
                v = v2; break
            v = v2
        f = V @ v; lam = float(np.sum(pi * f ** 3))
        if lam <= 0:
            continue
        c = v / lam; Tc = V.T @ ((pi * (V @ c))[:, None] * V)
        if np.sum(np.linalg.eigvalsh((Tc + Tc.T) / 2) > 0.5) != 1:
            continue                                   # not primitive by the canonical spectral count
        if np.linalg.norm(V.T @ (pi * (V @ c) ** 2) - c) > 1e-8 * max(1, np.linalg.norm(c)):
            continue                                   # not an idempotent to numerical precision
        if not any(np.linalg.norm(c - c0) < 1e-6 * np.linalg.norm(c0) for c0 in found):
            found.append(c)
    E = np.column_stack([V @ c for c in found]) if found else np.zeros((n, 0))
    return E


def blind_level(V, pi):
    E = primitive_idempotents(V, pi); k = V.shape[1]
    out = dict(k=k, n_prim=E.shape[1], E=E)
    if E.shape[1] != k:
        out['lab'] = None; return out                  # refuse
    out['lab'] = np.argmax(E, 1); out['sum_e'] = float(np.sqrt(np.sum(pi * (E.sum(1) - 1) ** 2)))
    out['gap'] = gap_metric(V, partition_basis(out['lab'], pi), pi); return out


# ------------------------------------------------------------------ PROOF / AUDIT side
def audit(lab_rec, lab, pi, V):
    K = len(np.unique(lab)); A = partition_basis(lab, pi); s = gap_metric(V, A, pi)
    p = np.array([pi[lab == b].sum() for b in np.unique(lab)]); Lam = p.min() ** -0.5; kappa = K * p.max(); CV = christoffel(V)
    eps_bd = s ** 2 * (11 * Lam + 2 * CV)
    out = dict(s=s, Lam=Lam, kappa=kappa, CV=CV, eps_bd=eps_bd,
               eps_glob_old=1 / (36 * K ** 1.5 * Lam ** 2), eps_glob_new=1 / (36 * kappa ** 1.5 * Lam ** 2))
    if lab_rec is None:
        out.update(m_tot=np.nan, m_rel=np.nan); return out
    tb, rb = np.unique(lab), np.unique(lab_rec)
    O = np.array([[pi[(lab == a) & (lab_rec == b)].sum() for b in rb] for a in tb])
    r, c = linear_sum_assignment(-O)
    mrel = 0.0; mtot = 0.0
    for a, b in zip(r, c):
        sym = pi[(lab == tb[a]) ^ (lab_rec == rb[b])].sum(); mrel = max(mrel, sym / p[a])
    mtot = 1.0 - O[r, c].sum()
    out.update(m_tot=float(mtot), m_rel=float(mrel),
               mrel_bd=((256 / 9) * K * eps_bd ** 2 + 8 * s ** 2) / p.min())
    return out


def nesting_defect(coarse, fine, pi):
    """nu = minimum pi-mass to reassign so that every fine block lies inside some coarse block (blind vs blind)."""
    nu = 0.0
    for f in np.unique(fine):
        sel = fine == f; nu += pi[sel].sum() - max(pi[sel & (coarse == c)].sum() for c in np.unique(coarse[sel]))
    return float(nu)


# ------------------------------------------------------------------ families
def FK(K, m, zeta):
    """C2-A2 rate family (ARCHITECTURE SUPPLIED BY FAMILY): K balanced basins of m states, intra conductance h/m,
    all-to-all inter conductance zeta*h/(m^2 (K-1)) (exit rate ~ zeta), weights nu."""
    lab = np.repeat(np.arange(K), m); n = K * m; x = np.arange(n); H = h(x[:, None], x[None, :])
    C = np.where(lab[:, None] == lab[None, :], H / m, zeta * H / (m * m * (K - 1))); return C, nu_rule(n), lab


def H1_bits(nbits, m=3, sw=0.01):
    """n independent two-basin systems; each: two complete graphs of m states (rate h/m), switching sw*h/m^2 per pair."""
    lab1 = np.repeat([0, 1], m); x = np.arange(2 * m); H = h(x[:, None], x[None, :])
    C1 = np.where(lab1[:, None] == lab1[None, :], H / m, sw * H / m); np.fill_diagonal(C1, 0)
    nu1 = 1.0 + 0.3 * np.cos(1.3 * x); L1, pi1 = from_conductances(C1, nu1)
    L = np.zeros((1, 1)); pi = np.ones(1); labs = []
    for b in range(nbits):
        L = np.kron(L, np.eye(2 * m)) + np.kron(np.eye(L.shape[0]), L1); pi = np.kron(pi, pi1)
    for b in range(nbits):                                        # hidden bit label of each micro state
        idx = np.array(list(itertools.product(range(2 * m), repeat=nbits)))[:, b]; labs.append(lab1[idx])
    lab = sum(labs[b] * 2 ** (nbits - 1 - b) for b in range(nbits))
    return L, pi, lab, labs, L1, pi1


def H2_ring(K, m=6, fwd=1.0, bwd=0.2, zeta=None):
    """K basins on a driven ring (non-reversible; ARCHITECTURE SUPPLIED BY FAMILY)."""
    zeta = zeta if zeta is not None else 1.0 / (10 * m)
    lab = np.repeat(np.arange(K), m); n = K * m; x = np.arange(n); H = h(x[:, None], x[None, :]); nu = nu_rule(n)
    G = np.where(lab[:, None] == lab[None, :], H / m / nu[:, None], 0.0)
    nxt = (lab[:, None] + 1) % K == lab[None, :]; prv = (lab[:, None] - 1) % K == lab[None, :]
    G = G + np.where(nxt, fwd * zeta * H / m, 0.0) + np.where(prv, bwd * zeta * H / m, 0.0)
    np.fill_diagonal(G, 0.0); np.fill_diagonal(G, -G.sum(1)); return G, lab


def P1_tree(d, m=4, zeta=0.02, b=2):
    """Nested ultrametric hierarchy (ARCHITECTURE SUPPLIED BY FAMILY): b^d leaf basins of m states; conductance between
    leaves whose lowest common ancestor is l levels up: zeta^l * h /(m^2 * (#leaves at that distance))."""
    nl = b ** d; lab = np.repeat(np.arange(nl), m); n = nl * m; x = np.arange(n); H = h(x[:, None], x[None, :])
    def lca(a, c):
        l = 0
        while a != c:
            a //= b; c //= b; l += 1
        return l
    LC = np.array([[lca(i, j) for j in range(nl)] for i in range(nl)])
    cnt = {l: (b ** l - b ** (l - 1)) for l in range(1, d + 1)}
    C = np.zeros((n, n))
    for i in range(nl):
        for j in range(nl):
            l = LC[i, j]; blk = np.ix_(lab == i, lab == j)
            C[blk] = H[blk] / m if l == 0 else (zeta ** l) * H[blk] / (m * m * cnt[l])
    labs = {j: lab // (b ** (d - j)) for j in range(0, d + 1)}     # hidden level-j partition (b^j blocks)
    return C, nu_rule(n), labs


def east(L, q):
    """E1 East model on sites 1..L, facilitating boundary eta_0 = 0; site x updates (->0 rate q, ->1 rate 1-q) iff
    eta_{x-1} = 0.  Reversible w.r.t. product Bernoulli(P[0]=q).  Returns symmetrised relaxation operator and pi."""
    n = 2 ** L; G = np.zeros((n, n)); conf = ((np.arange(n)[:, None] >> np.arange(L)[None, :]) & 1)   # bit x = eta_{x+1}
    for s in range(n):
        for x in range(L):
            left = 0 if x == 0 else conf[s, x - 1]
            if left == 0:
                t = s ^ (1 << x); G[s, t] = q if conf[s, x] == 1 else 1 - q
    np.fill_diagonal(G, -G.sum(1))
    pi = np.prod(np.where(conf == 0, q, 1 - q), axis=1); sq = np.sqrt(pi)
    Lsym = -(sq[:, None] * G / sq[None, :]); return (Lsym + Lsym.T) / 2, pi, conf


# ------------------------------------------------------------------ runs
def run_A2():
    print("== C2-A2 rate experiment (family FK; ARCHITECTURE SUPPLIED BY FAMILY): zeta_K = 0.3 K^-alpha, m = 6 ==")
    print("   columns: K, s, s*K^(3/4), s*K^(3/2), #primitive found (blind), m_tot, m_rel (audit), eps_bd, eps_glob old / new (G5-3 vs G5-3')")
    for alpha in [0.0, 0.5, 1.0]:
        for K in [4, 8, 16, 32, 64]:
            C, nu, lab = FK(K, 6, 0.3 * K ** (-alpha)); Ls, pi = from_conductances(C, nu); lam, U = spectrum(Ls, pi)
            V = U[:, :K]; bl = blind_level(V, pi); a = audit(bl['lab'], lab, pi, V); rr = cut_ratios(lam)
            print(f"  alpha={alpha:.1f} K={K:3d}: cut ratio@K {rr[K]:7.2f}; s {a['s']:.3e}; sK^.75 {a['s'] * K ** .75:.3e}; sK^1.5 {a['s'] * K ** 1.5:.3e};"
                  f" prim {bl['n_prim']:3d}/{K}; m_tot {a['m_tot']:.2e}; m_rel {a['m_rel']:.2e} (bd {a.get('mrel_bd', np.nan):.1e});"
                  f" eps_bd {a['eps_bd']:.2e}; glob old {a['eps_glob_old']:.1e} new {a['eps_glob_new']:.1e}; kappa {a['kappa']:.2f}")


def run_A2b():
    print("== C2-A2b breakdown sweep (FK, m = 6): where does blind recovery fail as the coupling zeta grows? ==")
    for K in [8, 32, 64]:
        for zeta in [0.3, 1.0, 3.0, 10.0]:
            C, nu, lab = FK(K, 6, zeta); Ls, pi = from_conductances(C, nu); lam, U = spectrum(Ls, pi)
            V = U[:, :K]; bl = blind_level(V, pi); a = audit(bl['lab'], lab, pi, V); rr = cut_ratios(lam)
            print(f"  K={K:3d} zeta={zeta:5.1f}: cut ratio@K {rr[K]:6.2f}; s {a['s']:.3e}; primitive {bl['n_prim']:3d}/{K};"
                  f" m_tot {a['m_tot']:.2e}; m_rel {a['m_rel']:.2e}")


def run_H1():
    print("== H1 independent metastable bits (K = 2^n; dual grade: LARGE K TRIVIAL + SUPPLIED) ==")
    for nb in [1, 2, 3, 4]:
        Ls, pi, lab, labs, L1, pi1 = H1_bits(nb); lam, U = spectrum(Ls, pi); K = 2 ** nb
        l1 = np.linalg.eigvalsh(L1); g1, ls = l1[2], l1[1]
        rr = cut_ratios(lam); slow = lam[:K]
        inner = max(rr[k] for k in range(2, K)) if K > 2 else float('nan')
        V = U[:, :K]; bl = blind_level(V, pi); a = audit(bl['lab'], lab, pi, V)
        print(f"  n={nb} (states {len(pi)}, K={K}): single bit lambda_s {ls:.3e}, g {g1:.3f}; n*lambda_s/g {nb * ls / g1:.3e};"
              f" cut ratio@K {rr[K]:.2f}; max ratio inside slow sector {inner:.2f} (bounded => depth 1)")
        print(f"      slow eigenvalues / lambda_s: {np.round(slow / ls, 3)}")
        print(f"      BLIND: primitive {bl['n_prim']}/{K}" + (f"; gap(V,A_blind) {bl['gap']:.2e}" if bl['lab'] is not None else " REFUSE")
              + f"; AUDIT m_rel {a['m_rel']:.2e}, s {a['s']:.2e}")
        if bl['lab'] is not None and nb >= 2:                 # post-hoc product test (hidden bit labels; audit only)
            mp = {bb: int(np.bincount(lab[bl['lab'] == bb], weights=pi[bl['lab'] == bb]).argmax()) for bb in np.unique(bl['lab'])}
            lab_h = np.array([mp[v] for v in bl['lab']]); M = np.zeros((len(pi), K)); M[np.arange(len(pi)), lab_h] = 1
            sq = np.sqrt(pi); G = -(Ls / sq[None, :] * sq[:, None]); F = M.T @ (pi[:, None] * G) @ M; pa = M.T @ pi; Q = F / pa[:, None]
            # Kronecker-sum reference from the single-bit lumped 2x2 generator
            M1 = np.zeros((len(pi1), 2)); M1[np.arange(len(pi1)), np.repeat([0, 1], len(pi1) // 2)] = 1
            sq1 = np.sqrt(pi1); G1 = -(L1 / sq1[None, :] * sq1[:, None]); q1 = (M1.T @ (pi1[:, None] * G1) @ M1) / (M1.T @ pi1)[:, None]
            Qref = np.zeros((1, 1))
            for _ in range(nb):
                Qref = np.kron(Qref, np.eye(2)) + np.kron(np.eye(Qref.shape[0]), q1)
            print(f"      post-hoc: max|Q - Kronecker-sum of single-bit macro generators| {np.abs(Q - Qref).max():.1e}"
                  f" (|Q| max {np.abs(Q).max():.1e}) -> independent flips")


def run_H2_H3():
    print("== H2 driven ring of K basins (dual grade: CLOCK + SUPPLIED) and H3 Markov causal-state count ==")
    for K in [4, 8, 16, 32]:
        G, lab = H2_ring(K); n = len(lab)
        w, Vl = np.linalg.eig(G.T); pi = np.real(Vl[:, np.argmin(np.abs(w))]); pi /= pi.sum()
        lam, Vr = np.linalg.eig(G); o = np.argsort(-lam.real); lam, Vr = lam[o], Vr[:, o]
        r = -lam.real; ratio = r[K] / r[K - 1]
        E = np.real(Vr[:, :K] @ np.linalg.inv(Vr)[:K, :]); sq = np.sqrt(pi)
        Q_, _ = np.linalg.qr(sq[:, None] * (E @ np.random.default_rng(0).standard_normal((n, K)))); V = Q_ / sq[:, None]
        bl = blind_level(V, pi); a = audit(bl['lab'], lab, pi, V)
        line = f"  K={K:3d} n={n}: decay-rate ratio@K {ratio:.2f}; primitive {bl['n_prim']}/{K}; AUDIT m_rel {a['m_rel']:.2e}"
        if bl['lab'] is not None:
            M = np.zeros((n, K)); M[np.arange(n), bl['lab']] = 1; F = M.T @ (pi[:, None] * G) @ M; pa = M.T @ pi; Q = F / pa[:, None]
            J = pa[:, None] * Q; np.fill_diagonal(J, 0); J = J - J.T
            sup = np.abs(pa[:, None] * Q) > 1e-14 * np.abs(pa[:, None] * Q).max(); np.fill_diagonal(sup, False)
            E_ = int(np.sum(np.triu(sup | sup.T, 1))); cyc = E_ - K + 1
            P = np.eye(K) + Q / np.max(-np.diag(Q)); rows = {tuple(np.round(P[i], 12)) for i in range(K)}
            line += f"; macro cycle-space dim {cyc}; |J| on ring edges min {np.min(np.abs(J[np.abs(J) > 1e-14])):.2e} max {np.abs(J).max():.2e}"
            line += f"; H3: distinct rows of uniformised macro kernel (= causal states of the 1st-order Markov macro-process) {len(rows)}"
        print(line)


def run_P1():
    print("== P1 nested ultrametric hierarchy (b=2; ARCHITECTURE SUPPLIED BY FAMILY): depth d, zeta = 0.02 ==")
    for d in [1, 2, 3, 4, 5]:
        C, nu, labs = P1_tree(d); Ls, pi = from_conductances(C, nu); lam, U = spectrum(Ls, pi); rr = cut_ratios(lam)
        ranks = [2 ** j for j in range(1, d + 1)]
        print(f"  d={d} (states {len(pi)}): cut ratios at ranks {ranks}: {[round(rr[k], 1) for k in ranks]};"
              f" largest ratio elsewhere {max(v for k, v in rr.items() if k not in ranks):.2f}")
        recs = {}
        for j, k in enumerate(ranks, start=1):
            V = U[:, :k]; bl = blind_level(V, pi); a = audit(bl['lab'], labs[j], pi, V); recs[k] = bl
            print(f"      level {j} (rank {k}): primitive {bl['n_prim']}/{k}; AUDIT s {a['s']:.2e}, m_tot {a['m_tot']:.1e}, m_rel {a['m_rel']:.1e}")
        for k1, k2 in zip(ranks[:-1], ranks[1:]):
            if recs[k1]['lab'] is not None and recs[k2]['lab'] is not None:
                print(f"      nesting defect nu(rank {k1} ⊃ rank {k2}) {nesting_defect(recs[k1]['lab'], recs[k2]['lab'], pi):.1e} (blind vs blind)")


def run_P1_scaling():
    print("== P1 cut ratios vs zeta (depth d = 3): the level cuts diverge only as zeta -> 0 (supplied scaling) ==")
    for zeta in [0.1, 0.05, 0.02, 0.01]:
        C, nu, labs = P1_tree(3, zeta=zeta); Ls, pi = from_conductances(C, nu); lam, U = spectrum(Ls, pi); rr = cut_ratios(lam)
        print(f"  zeta={zeta}: ratios at ranks 2,4,8: {[round(rr[k], 1) for k in (2, 4, 8)]}")


def run_E1():
    print("== E1 East model (bounded-description, no explicit hierarchy) ==")
    print("  (a) fixed q, growing L (no supplied scaling): top decay-rate cut ratios")
    for q in [0.1, 0.3]:
        for L in [4, 6, 8, 10, 12]:
            Ls, pi, conf = east(L, q); lam = np.linalg.eigvalsh(Ls); lam = np.maximum(lam, 0); rr = cut_ratios(lam)
            top = sorted(rr.items(), key=lambda t: -t[1])[:3]
            print(f"    q={q} L={L:2d}: top (ratio, rank) {[(round(v, 2), k) for k, v in top]}; gap {lam[1]:.3e}")
    print("  (b) fixed L, q -> 0 (supplied scaling): cut ratios at fixed ranks and their growth exponent in 1/q")
    qs = [0.1, 0.05, 0.02, 0.01]
    for L in [3, 4, 5, 6, 8, 10]:
        R = {}
        for q in qs:
            Ls, pi, conf = east(L, q); lam = np.maximum(np.linalg.eigvalsh(Ls), 0); R[q] = cut_ratios(lam)
        ks = sorted(R[qs[0]].keys())
        slope = {k: np.polyfit(np.log(1 / np.array(qs)), np.log([R[q][k] for q in qs]), 1)[0] for k in ks}
        div = [k for k in ks if slope[k] > 0.5]
        print(f"    L={L}: ranks whose ratio grows like q^-a with a > 0.5 (numerical classification): "
              f"{[(k, round(slope[k], 2), round(R[qs[-1]][k], 1)) for k in div]} -> depth d(L) = {len(div)}")
        if L in (4, 6, 8, 10):
            Ls, pi, conf = east(L, qs[-1]); lam, U = spectrum(Ls, pi)
            for k in div:
                V = U[:, :k]; bl = blind_level(V, pi)
                if bl['lab'] is None:
                    print(f"      q={qs[-1]} rank {k}: primitive {bl['n_prim']}/{k} -> REFUSE; Delta_HS {defects(V, pi)[0]:.2e};"
                          f" smallest pi-mass of a primitive idempotent's support {min([pi[bl['E'][:, i] > 0.5].sum() for i in range(bl['E'].shape[1])] or [0]):.2e}"); continue
                pm = np.array([pi[bl['lab'] == b].sum() for b in np.unique(bl['lab'])])
                print(f"      (Delta_HS of this slow space {defects(V, pi)[0]:.2e})")
                tops = []
                for bb in np.unique(bl['lab']):                     # post-hoc description only: modal configuration
                    sel = np.nonzero(bl['lab'] == bb)[0]; c0 = conf[sel[np.argmax(pi[sel])]]
                    tops.append(''.join(str(v) for v in c0))
                print(f"      post-hoc modal configuration (eta_1..eta_L) of each recovered block: {sorted(tops)}")
                print(f"      q={qs[-1]} rank {k}: primitive {bl['n_prim']}/{k}; recovered block masses {np.round(np.sort(pm), 5)};"
                      f" gap(V,A_blind) {bl['gap']:.2e}; |sum e - 1| {bl['sum_e']:.1e}")


if __name__ == "__main__":
    run_A2(); run_A2b(); run_H1(); run_H2_H3(); run_P1(); run_P1_scaling(); run_E1()
