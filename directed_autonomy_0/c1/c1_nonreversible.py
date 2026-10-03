"""DA0 / C1 -- non-reversible endogenous macrostructure.  NUMERICAL ILLUSTRATION; independent code path, not reviewer.

Firewall by code structure (C0, DC1-01..04):
  * BLIND side (`blind_pipeline`) sees only the generator G and its derived objects: pi, G* (pi-adjoint), J, EP,
    spectrum, Riesz projectors of conjugation-closed spectral sets, and a CERTIFIED cut rank.  It calls the frozen RA0/G5
    recovery map R(V, pi), forms the lumped macro-generator Q, the macro current, macro EP and the macro cycle space.
  * The certified rank comes from Theorem C1-T (proof side, hidden labels allowed in the PROOF), never from labels in
    the pipeline.  `proof_side` (hidden labels: G0/G1 split, E0, eta, g, a, bounds, audit) runs after the blind side.
  * No cycle basis is used anywhere: only J, the cycle-space dimension, and (when that dimension is 1) the unique cycle.
Norm: L^2(pi).  pi is the unique stationary law (irreducibility priced).  All edges are mutually supported, so the
finite EP formula applies (DC1-01).
Computational tolerances (numerical equality only): zero-mode 1e-10 relative; support of macro edges 1e-14 relative.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'reflexive_autonomy_0', 'g4'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'reflexive_autonomy_0', 'g5'))
from g4_algebra import F1, h, nu_rule, gap_metric, partition_basis, christoffel
from g5_recovery import recover

np.set_printoptions(precision=4, suppress=True)
ZERO_J = 1e-13      # computational equality tolerance for "macro current is zero" (round-off level; declared, not a model input)


# ------------------------------------------------------------------ canonical objects of (G, pi)
def stationary(G):
    w, V = np.linalg.eig(G.T); i = np.argmin(np.abs(w)); p = np.real(V[:, i]); return p / p.sum()


def pinorm(X, pi):
    s = np.sqrt(pi); return float(np.linalg.norm(s[:, None] * X / s[None, :], 2))


def adjoint(G, pi):
    return (G.T * pi[None, :]) / pi[:, None]                        # G*(x,y) = pi_y G(y,x) / pi_x


def current(G, pi):
    F = pi[:, None] * G; np.fill_diagonal(F, 0.0); return F - F.T


def entropy_production(G, pi):
    F = pi[:, None] * G; np.fill_diagonal(F, 0.0); m = (F > 0) & (F.T > 0)
    assert np.all(m == (F > 0)), "one-way edge present: finite EP formula not applicable (DC1-01)"
    J = F - F.T; return 0.5 * float(np.sum(J[m] * np.log(F[m] / F.T[m])))


def cycle_space_dim(Adj):
    """dim of cycle space of the undirected support graph = E - V + (#components)."""
    n = Adj.shape[0]; E = int(np.sum(np.triu(Adj | Adj.T, 1)))
    seen, comp = np.zeros(n, bool), 0
    for s0 in range(n):
        if seen[s0]:
            continue
        comp += 1; stack = [s0]; seen[s0] = True
        while stack:
            u = stack.pop()
            for v in np.nonzero(Adj[u] | Adj[:, u])[0]:
                if not seen[v]:
                    seen[v] = True; stack.append(v)
    return E - n + comp


def spectrum_sorted(G):
    lam, Vr = np.linalg.eig(G); r = -lam.real; o = np.argsort(r, kind='stable')
    return lam[o], Vr[:, o]


def cuts(lam, top=3):
    r = -lam.real; z = int(np.sum(r <= 1e-10 * r.max())); out = []
    for k in range(z + 1, len(r)):
        if abs(lam[k] - np.conj(lam[k - 1])) < 1e-9 * abs(lam[k]) and abs(lam[k].imag) > 0:
            continue                                                  # would split a conjugate pair: not admissible
        out.append((float(r[k] / r[k - 1]), k))
    out.sort(reverse=True); return out[:top]


def riesz(lam, Vr, k):
    """Riesz projector of the conjugation-closed set {lam_0..lam_{k-1}} (eigen-expansion; real by closure)."""
    Vinv = np.linalg.inv(Vr); E = Vr[:, :k] @ Vinv[:k, :]
    return np.real(E), float(np.max(np.abs(np.imag(E))))


def contour_K(G, pi, lam, k, npts=240):
    """K_N(Gamma) on the circle |z| = rho, rho = sqrt(max|slow| * min|fast|) (spectrally defined contour)."""
    a, b = np.max(np.abs(lam[:k])), np.min(np.abs(lam[k:])); rho = np.sqrt(max(a, 1e-300) * b)
    s = np.sqrt(pi); Gs = s[:, None] * G / s[None, :]; n = len(pi); worst = 0.0
    for t in np.linspace(0, 2 * np.pi, npts, endpoint=False):
        z = rho * np.exp(1j * t); worst = max(worst, 1.0 / np.linalg.svd(z * np.eye(n) - Gs, compute_uv=False)[-1])
    return rho * worst, rho, a < rho < b


def pi_orthonormal(B, pi):
    s = np.sqrt(pi); Q, _ = np.linalg.qr(s[:, None] * B); return Q / s[:, None]


def lump(G, pi, lab):
    """pi-weighted lumped macro-generator Q(a,b) = sum_{x in a, y in b} pi_x G_xy / pi(a)  (canonical given Pi, pi)."""
    blocks = np.unique(lab); K = len(blocks); M = np.zeros((len(pi), K)); M[np.arange(len(pi)), np.searchsorted(blocks, lab)] = 1
    F = M.T @ (pi[:, None] * G) @ M; pa = M.T @ pi; Q = F / pa[:, None]; return Q, pa


# ------------------------------------------------------------------ BLIND pipeline
def blind_pipeline(G, rank):
    pi = stationary(G); out = dict(pi=pi, EP=entropy_production(G, pi))
    lam, Vr = spectrum_sorted(G); out['cuts'] = cuts(lam); out['lam'] = lam
    if rank is None:
        return out
    E, imres = riesz(lam, Vr, rank)
    out.update(E=E, normE=pinorm(E, pi), imres=imres, idem_res=float(np.max(np.abs(E @ E - E))),
               comm_res=float(np.max(np.abs(G @ E - E @ G))), condV=float(np.linalg.cond(np.sqrt(pi)[:, None] * Vr)))
    out['K'], out['rho'], out['sep'] = contour_K(G, pi, lam, rank)
    V = pi_orthonormal(E @ np.random.default_rng(0).standard_normal((len(pi), rank)), pi)   # a basis of Ran E (R is basis-free)
    out['V'] = V; r = recover(V, pi); out['rec'] = r
    if r['lab'] is None:
        return out
    Q, pa = lump(G, pi, r['lab']); JQ = pa[:, None] * Q; np.fill_diagonal(JQ, 0); JQ = JQ - JQ.T
    sup = (np.abs(pa[:, None] * Q) > 1e-14 * np.abs(pa[:, None] * Q).max()); np.fill_diagonal(sup, False)
    out.update(Q=Q, pa=pa, JQ=JQ, EPQ=entropy_production(Q, pa), cyc=cycle_space_dim(sup))
    return out


def macro_cycle(JQ):
    """If the macro support graph has a 1-dim cycle space, J_Q is a multiple of the unique cycle: report its oriented
    vertex sequence (by following positive current) and the circulation value.  No basis is chosen."""
    K = JQ.shape[0]; start = 0; seq = [start]; cur = start
    if np.abs(JQ).max() <= ZERO_J:                                   # zero macro current (computational equality)
        return None, 0.0
    for _ in range(K):
        nxt = int(np.argmax(JQ[cur])); 
        if JQ[cur, nxt] <= 0:
            return None, 0.0
        if nxt == start:
            break
        seq.append(nxt); cur = nxt
    if len(seq) < 3:
        return None, 0.0
    return seq, float(JQ[seq[0], seq[1]])


def canon_perm(lab_a, lab_b):
    """Do two partitions coincide up to block relabelling?  (blind vs blind comparison; no hidden labels)."""
    if lab_a is None or lab_b is None:
        return False, None
    m = {}
    for x, y in zip(lab_a, lab_b):
        if m.setdefault(x, y) != y:
            return False, None
    return len(set(m.values())) == len(m), m


# ------------------------------------------------------------------ PROOF side (hidden labels; after blind side)
def proof_side(G, pi, lab, rank, E):
    K = len(np.unique(lab)); same = lab[:, None] == lab[None, :]
    G0 = np.where(same, G, 0.0); np.fill_diagonal(G0, 0.0); np.fill_diagonal(G0, -G0.sum(1)); G1 = G - G0
    eta = pinorm(G1, pi)
    lam0, V0 = spectrum_sorted(G0); E0, _ = riesz(lam0, V0, K); e = pinorm(E0, pi)
    g = 0.5 * float(np.min(-lam0[K:].real))                                   # half the fast abscissa of G0
    s = np.sqrt(pi); G0s = s[:, None] * G0 / s[None, :]; P = s[:, None] * (np.eye(len(pi)) - E0) / s[None, :]; n = len(pi)
    Y = 3 * float(np.linalg.norm(G0s, 2)); a = 0.0
    for y in np.concatenate([np.linspace(-Y, Y, 161), np.linspace(-3 * g, 3 * g, 81)]):
        z = -g + 1j * y; a = max(a, float(np.linalg.norm(np.linalg.solve(z * np.eye(n) - G0s, P), 2)))
    ok = eta * a < 1
    ell = eta * e / (1 - eta * a) if ok else np.inf
    rho = np.sqrt(ell * g) if ok else np.nan; R0 = e / rho + a if ok else np.inf; d = eta * R0
    Ebd = 2 * e * a * eta + (e + rho * a) * d ** 2 / (1 - d) if ok and d < 1 else np.inf   # first-order residue + remainder
    A = partition_basis(lab, pi); V = pi_orthonormal(E @ np.random.default_rng(0).standard_normal((n, rank)), pi)
    return dict(eta=eta, e=e, g=g, a=a, ag=a * g, ell=ell, ratio_bd=g / ell if ok else 0.0, Ebd=Ebd,
                EmE0=pinorm(E - E0, pi), gapVA=gap_metric(V, A, pi), slow_ok=bool(np.all(np.abs(np.linalg.eigvals(G)[np.argsort(-np.linalg.eigvals(G).real)][:rank]) <= ell * (1 + 1e-9))) if ok else False)


def audit(lab_rec, lab, pi):
    mis = 0.0
    for b in np.unique(lab_rec):
        sel = lab_rec == b; vals = np.unique(lab[sel])
        maj = vals[np.argmax([pi[sel & (lab == v)].sum() for v in vals])]; mis += pi[sel & (lab != maj)].sum()
    return float(mis)


# ------------------------------------------------------------------ families (deterministic rules)
def F1_flux(M):
    C, nu, lab = F1(M); pi = nu / nu.sum(); Fs = C / nu.sum(); np.fill_diagonal(Fs, 0.0); return Fs, pi, lab


def gen_from_flux(F, pi):
    G = F / pi[:, None]; np.fill_diagonal(G, 0.0); np.fill_diagonal(G, -G.sum(1)); return G


def R0(M):
    Fs, pi, lab = F1_flux(M); return gen_from_flux(Fs, pi), lab


def circ_macro(Fs, lab, theta=0.9):
    """R1a: divergence-free inter-block circulation 0->1->2->0, F_a = theta*c*sigma/(n_b n_b'); |F_a| <= F_s."""
    n_b = np.bincount(lab); K = len(n_b); sig = np.zeros((K, K))
    for b in range(K):
        sig[b, (b + 1) % K] = 1; sig[(b + 1) % K, b] = -1
    inter = lab[:, None] != lab[None, :]; nn = n_b[lab][:, None] * n_b[lab][None, :]
    c = float(np.min((Fs * nn)[inter])); return theta * c * sig[lab][:, lab] / nn


def circ_intra(Fs, lab, theta=0.9):
    """R1b: per block, Hamiltonian cycle in index order with flux theta*min F_s on its edges; no inter-block current."""
    Fa = np.zeros_like(Fs)
    for b in np.unique(lab):
        idx = np.nonzero(lab == b)[0]; m = len(idx)
        f = theta * min(Fs[idx[i], idx[(i + 1) % m]] for i in range(m))
        for i in range(m):
            x, y = idx[i], idx[(i + 1) % m]; Fa[x, y] += f; Fa[y, x] -= f
    return Fa


def R1a(M):
    Fs, pi, lab = F1_flux(M); return gen_from_flux(Fs + circ_macro(Fs, lab), pi), lab


def R1b(M):
    Fs, pi, lab = F1_flux(M); return gen_from_flux(Fs + circ_intra(Fs, lab), pi), lab


def R2(M, fwd=1.0, bwd=0.2):
    """4 basins on a driven ring (sizes M, 1.25M, 1.5M, 1.75M); intra reversible conductances h/n_b, weights nu;
    inter only between ring neighbours: b->b+1 rate fwd*h/M^2, b->b-1 rate bwd*h/M^2 (non-reversible)."""
    sizes = [M, int(1.25 * M), int(1.5 * M), int(1.75 * M)]; lab = np.repeat(np.arange(4), sizes); n = lab.size
    x = np.arange(n); H = h(x[:, None], x[None, :]); nu = nu_rule(n); ns = np.array(sizes)[lab]
    G = np.where(lab[:, None] == lab[None, :], H / ns[:, None] / nu[:, None], 0.0)
    nxt = (lab[:, None] + 1) % 4 == lab[None, :]; prv = (lab[:, None] - 1) % 4 == lab[None, :]
    G = G + np.where(nxt, fwd * H / M ** 2, 0.0) + np.where(prv, bwd * H / M ** 2, 0.0)
    np.fill_diagonal(G, 0.0); np.fill_diagonal(G, -G.sum(1)); return G, lab


def drift_ring(L, off=0):
    xx = np.arange(L) + off; G = np.zeros((L, L))
    for i in range(L):
        G[i, (i + 1) % L] = 1.0 * (1 + 0.3 * np.cos(1.1 * xx[i])); G[i, (i - 1) % L] = 0.4 * (1 + 0.2 * np.cos(0.37 * xx[i]))
    return G


def R3(L):
    G = drift_ring(L); np.fill_diagonal(G, -G.sum(1)); return G, None


def R4(M):
    """Purpose: exercise DC1-02 (non-normal blocks).  3 drifting site-modulated rings of sizes M, 1.5M, 2M;
    all-to-all inter-block rates h/M^4."""
    sizes = [M, int(1.5 * M), 2 * M]; lab = np.repeat(np.arange(3), sizes); n = lab.size; G = np.zeros((n, n)); o = 0
    for b, m in enumerate(sizes):
        G[o:o + m, o:o + m] = drift_ring(m, off=o); o += m
    x = np.arange(n); H = h(x[:, None], x[None, :])
    G = G + np.where(lab[:, None] != lab[None, :], H / M ** 4, 0.0)
    np.fill_diagonal(G, 0.0); np.fill_diagonal(G, -G.sum(1)); return G, lab


def K1_controls():
    print("== K1 controls ==")
    a, b = 0.3, 0.7; G = np.array([[-a, a], [b, -b]]); pi = stationary(G)
    print(f"  two-state bit (any rates): EP = {entropy_production(G, pi):.2e}  (every 2-state chain satisfies detailed balance)")
    f, r = 1.0, 0.1; G = np.array([[-(f + r), f, r], [r, -(f + r), f], [f, r, -(f + r)]]); pi = stationary(G); J = current(G, pi)
    lam = np.sort_complex(np.linalg.eigvals(G))
    print(f"  driven 3-state loop (thermostat/clock): EP = {entropy_production(G, pi):.3f}; cycle-space dim "
          f"{cycle_space_dim(np.abs(G) > 0)}; J(0->1) = {J[0, 1]:.3f} = J(1->2) = {J[1, 2]:.3f} = J(2->0) = {J[2, 0]:.3f};"
          f" eigenvalues {np.round(lam, 3)} -> canonical directed cycle with NO metastable cut")


# ------------------------------------------------------------------ runs
def run_partition_family(name, build, sizes, K):
    print(f"== {name} ==")
    for M in sizes:
        G, lab = build(M); n = G.shape[0]
        b = blind_pipeline(G, K); pi = b['pi']
        print(f"  M={M:3d} n={n}: EP {b['EP']:.3e}; top decay-rate cuts (ratio, rank) {[(round(x, 2), k) for x, k in b['cuts']]}")
        r = b['rec']
        print(f"    BLIND rank {K}: ||E||_pi {b['normE']:.3f}, K(Gamma) {b['K']:.3f} (circle separates: {b['sep']}), cond(eigvecs) {b['condV']:.1e},"
              f" |Im E| {b['imres']:.1e}, |E^2-E| {b['idem_res']:.1e}, |GE-EG| {b['comm_res']:.1e}")
        print(f"          R(V): idempotents {r['n_idem']}, primitive {r['n_prim']}" + (f", gap(V,A_blind) {r['gap']:.3e}" if r['lab'] is not None else " -> REFUSE"))
        if r['lab'] is not None:
            seq, j = macro_cycle(b['JQ']) if b['cyc'] == 1 else (None, 0.0)
            print(f"          macro: EP_Q {b['EPQ']:.3e}; macro cycle-space dim {b['cyc']}; |J_Q| max {np.abs(b['JQ']).max():.3e};"
                  + (f" canonical macro cycle {seq} with circulation {j:.3e}" if seq else " no canonical macro cycle / zero macro current"))
            Gs = adjoint(G, pi); bs = blind_pipeline(Gs, K); same, mp = canon_perm(r['lab'], bs['rec']['lab'])
            if same:
                Js = np.zeros_like(b['JQ']); inv = {v: k for k, v in mp.items()}
                for i in range(K):
                    for jj in range(K):
                        Js[i, jj] = bs['JQ'][mp[i], mp[jj]]
                cov = float(np.max(np.abs(Js + b['JQ'])))
            else:
                cov = np.nan
            print(f"          COVARIANCE (G -> G*, blind): same partition {same}; max|J_Q* + J_Q| {cov:.1e}; EP(G*) {bs['EP']:.3e};"
                  f" gap(V, V*) {gap_metric(b['V'], bs['V'], pi):.3e}")
        if lab is not None:                                               # ---- proof side, after blind side
            JB = np.zeros((K, K)); Jm = current(G, pi)
            for i in range(K):
                for jj in range(K):
                    JB[i, jj] = Jm[np.ix_(lab == i, lab == jj)].sum()
            p = proof_side(G, pi, lab, K, b['E'])
            mis = audit(r['lab'], lab, pi) if r['lab'] is not None else 1.0
            qmax = float(np.max(-np.diag(G)))
            mp_h = {bl: int(np.bincount(lab[r['lab'] == bl], weights=pi[r['lab'] == bl]).argmax()) for bl in np.unique(r['lab'])} if r['lab'] is not None else {}
            if r['lab'] is not None and b['cyc'] == 1:
                seq, _ = macro_cycle(b['JQ']); seqB, _ = macro_cycle(JB)
                if seq:
                    print(f"    AUDIT blind macro cycle in hidden labels {[mp_h[i] for i in seq]}; hidden J_B cycle {seqB}")
            print(f"    PROOF eta {p['eta']:.2e}, ||E0|| {p['e']:.3f}, g {p['g']:.3e}, a*g {p['ag']:.2f}, eta/g {p['eta'] / p['g']:.2e};"
                  f" ell {p['ell']:.2e} (slow |lam|<=ell: {p['slow_ok']}); ratio bound g/ell {p['ratio_bd']:.2f}")
            print(f"          ||E-E0|| {p['EmE0']:.3e} <= bound {p['Ebd']:.3e};  gap(V,A_true) {p['gapVA']:.3e};"
                  f"  AUDIT misclassified {mis:.2e};  hidden macro current |J_B| max {np.abs(JB).max():.3e},"
                  f" 4*qmax*mis {4 * qmax * mis:.1e}")


def run_flux_family(name, build, sizes):
    print(f"== {name} ==")
    for L in sizes:
        G, _ = build(L); b = blind_pipeline(G, None); pi = b['pi']; J = current(G, pi)
        lam = b['lam']; r = -lam.real
        print(f"  L={L:3d}: EP {b['EP']:.3e}; top decay-rate cuts {[(round(x, 3), k) for x, k in b['cuts']]};"
              f" cycle-space dim {cycle_space_dim(np.abs(G) > 0)}; J on ring edges min {min(J[i, (i + 1) % L] for i in range(L)):.4e}"
              f" max {max(J[i, (i + 1) % L] for i in range(L)):.4e}; slowest pair {np.round(lam[1:3], 4)}")


def matched_pair_check():
    print("== R1 matched-pair check (DC1-03): same pi, same symmetric part S, different antisymmetric part A ==")
    for M in [20, 80]:
        Ga, _ = R1a(M); Gb, lab = R1b(M); pa, pb = stationary(Ga), stationary(Gb)
        Sa = (Ga + adjoint(Ga, pa)) / 2; Sb = (Gb + adjoint(Gb, pb)) / 2; Aa = (Ga - adjoint(Ga, pa)) / 2; Ab = (Gb - adjoint(Gb, pb)) / 2
        cos = float(np.sum(pa[:, None] * Aa * Ab) / np.sqrt(np.sum(pa[:, None] * Aa ** 2) * np.sum(pa[:, None] * Ab ** 2)))
        print(f"  M={M}: max|pi_a - pi_b| {np.abs(pa - pb).max():.1e}; max|S_a - S_b| {np.abs(Sa - Sb).max():.1e}; "
              f"cosine(A_a, A_b) in pi-weighted HS {cos:.2e} (non-collinear)")


def probe_defect_rings():
    print("== DC1-02 probe (diagnostic): driven rings with one slow defect bond -- does non-normal conditioning explode? ==")
    for (f, bk, sl) in [(1, 0.05, 0.05), (1, 0.02, 0.01)]:
        for L in [10, 20, 40]:
            G = np.zeros((L, L))
            for i in range(L):
                G[i, (i + 1) % L] = f; G[i, (i - 1) % L] = bk
            G[L - 1, 0] = sl; G[0, L - 1] = bk * sl; np.fill_diagonal(G, -G.sum(1))
            pi = stationary(G); lam, V = spectrum_sorted(G); E, _ = riesz(lam, V, 3)
            print(f"  fwd {f} bwd {bk} slow-bond {sl} L={L}: cond_pi(eigvecs) {np.linalg.cond(np.sqrt(pi)[:, None] * V):.1e};"
                  f" ||E(rank 3)||_pi {pinorm(E, pi):.2f}; pi max/min {pi.max() / pi.min():.1e}")


if __name__ == "__main__":
    K1_controls(); matched_pair_check(); probe_defect_rings()
    run_partition_family("R0  F1 reversible (RA0 baseline, K4)", R0, [20, 40, 80], 3)
    run_partition_family("R1a F1 skeleton + inter-block circulation (same pi, same S)", R1a, [20, 40, 80], 3)
    run_partition_family("R1b F1 skeleton + intra-block circulation (same pi, same S)", R1b, [20, 40, 80], 3)
    run_partition_family("R2  4 metastable basins on a driven ring (fwd 1, bwd 0.2)", R2, [10, 20, 40, 80], 4)
    run_partition_family("R4  non-normal stress: 3 drifting rings, inter h/M^4", R4, [10, 20, 40, 80], 3)
    run_flux_family("R3  biased site-modulated ring (no metastability)", R3, [25, 50, 100, 200])
