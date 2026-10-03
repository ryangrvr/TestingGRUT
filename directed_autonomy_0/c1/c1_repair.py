"""DA0 / C1 REPAIR 01 -- (i) N4a/N4b audit of R1a, R2, R4; (ii) R5 multi-cycle metastable macro-current control.
NUMERICAL ILLUSTRATION; independent code path, not independent reviewer.  The original C1 script and log are preserved.

R5 purpose (C1R-07): close the DC1-01 multi-cycle case.  K=4 basins, complete inter-basin graph (macro cycle space
dimension 3), non-reversible current = superposition of two independent circulations.  NO spanning tree, NO cycle basis,
NO decomposition of J into cycles is used anywhere; the macro current matrix J_Q itself is reported.
"""
import numpy as np
from c1_nonreversible import (blind_pipeline, adjoint, current, cycle_space_dim, canon_perm, gen_from_flux, F1_flux,
                              R1a, R2, R4, h, nu_rule, spectrum_sorted, pinorm, audit, ZERO_J)
from g4_algebra import gap_metric

np.set_printoptions(precision=4, suppress=True, linewidth=160)


def n4_audit(name, build, sizes, K):
    """Proof side (hidden labels): eta = ||G1||_pi, g = half fast abscissa of G0, q_max, j_min = min nonzero |J_B|."""
    print(f"== N4a/N4b audit: {name} ==")
    for M in sizes:
        G, lab = build(M); pi = blind_pipeline(G, None)['pi']
        same = lab[:, None] == lab[None, :]
        G0 = np.where(same, G, 0.0); np.fill_diagonal(G0, 0.0); np.fill_diagonal(G0, -G0.sum(1)); eta = pinorm(G - G0, pi)
        lam0, _ = spectrum_sorted(G0); g = 0.5 * float(np.min(-lam0[K:].real)); q = float(np.max(-np.diag(G)))
        Jm = current(G, pi); JB = np.array([[Jm[np.ix_(lab == i, lab == j)].sum() for j in range(K)] for i in range(K)])
        nz = np.abs(JB[np.triu_indices(K, 1)]); jmin = float(nz[nz > ZERO_J].min()) if np.any(nz > ZERO_J) else 0.0
        print(f"  M={M:3d}: eta {eta:.3e}, g {g:.3e}, q_max {q:.2f}, j_min {jmin:.3e}; j_min/eta {jmin / eta:.3e} (N4a: bounded below?);"
              f" q_max*eta/g^2 {q * eta / g ** 2:.3e} (N4b specialisation: -> 0?)")


def circ3(Fs, lab, cyc, theta):
    """Divergence-free circulation around basins cyc=(a,b,c): F_a(x,y) = theta*c0*sigma/(n_bx n_by), sigma=+1 along
    a->b->c->a.  Product weights 1/n_b make each node's in/out circulation equal (divergence-free).  This is a
    CONSTRUCTION of the supplied D (two chosen circulations), not a decomposition used by the analysis."""
    n_b = np.bincount(lab); K = len(n_b); sig = np.zeros((K, K))
    for i in range(3):
        u, v = cyc[i], cyc[(i + 1) % 3]; sig[u, v] = 1; sig[v, u] = -1
    inter = lab[:, None] != lab[None, :]; nn = n_b[lab][:, None] * n_b[lab][None, :]
    c0 = float(np.min((Fs * nn)[inter])); return theta * c0 * sig[lab][:, lab] / nn


def R5(M):
    """K=4 basins (sizes M, 1.25M, 1.5M, 1.75M), reversible skeleton: intra h/n_b, inter h*kappa/M^2 for ALL basin pairs
    (complete macro graph), stationary weights nu; plus two independent circulations 0->1->2->0 (theta 0.5) and
    0->2->3->0 (theta 0.3)."""
    sizes = [M, int(1.25 * M), int(1.5 * M), int(1.75 * M)]; lab = np.repeat(np.arange(4), sizes); n = lab.size
    x = np.arange(n); H = h(x[:, None], x[None, :]); nu = nu_rule(n); ns = np.array(sizes)[lab]
    kap = np.array([[0, 1.0, 0.7, 1.3], [1.0, 0, 0.8, 0.6], [0.7, 0.8, 0, 1.1], [1.3, 0.6, 1.1, 0]])
    C = np.where(lab[:, None] == lab[None, :], H / ns[:, None], H * kap[lab][:, lab] / M ** 2); np.fill_diagonal(C, 0)
    pi = nu / nu.sum(); Fs = C / nu.sum()
    Fa = circ3(Fs, lab, (0, 1, 2), 0.5) + circ3(Fs, lab, (0, 2, 3), 0.3)
    assert np.all(np.abs(Fa) <= Fs + 1e-300), "circulation exceeds symmetric flux"
    return gen_from_flux(Fs + Fa, pi), lab


def single_cycle_supported(JQ):
    """Basis-free test: is J_Q a multiple of ONE simple cycle?  True iff its support edges form a single cycle and every
    vertex of the support has degree 2 (then J is constant in magnitude along it)."""
    K = JQ.shape[0]; E = [(i, j) for i in range(K) for j in range(i + 1, K) if abs(JQ[i, j]) > ZERO_J]
    deg = np.zeros(K, int)
    for i, j in E:
        deg[i] += 1; deg[j] += 1
    V = np.nonzero(deg)[0]
    return len(E) == len(V) and np.all(deg[V] == 2), len(E)


def run_R5():
    print("== R5 multi-cycle metastable macro-current (C1R-07) ==")
    K = 4
    for M in [10, 20, 40, 80]:
        G, lab = R5(M); b = blind_pipeline(G, K); r = b['rec']; pi = b['pi']
        print(f"  M={M:3d} n={G.shape[0]}: micro EP {b['EP']:.3e}; top decay-rate cuts {[(round(x, 2), k) for x, k in b['cuts']]};"
              f" ||E||_pi {b['normE']:.3f}, K(Gamma) {b['K']:.3f}")
        print(f"    BLIND R(V): idempotents {r['n_idem']}, primitive {r['n_prim']}"
              + (f", gap(V,A_blind) {r['gap']:.3e}" if r['lab'] is not None else " -> REFUSE"))
        if r['lab'] is None:
            continue
        one, ne = single_cycle_supported(b['JQ'])
        print(f"    macro cycle-space dim {b['cyc']}; macro EP_Q {b['EPQ']:.3e}; J_Q supported on {ne} macro edges;"
              f" J_Q is a multiple of a single cycle: {one}")
        print("    canonical macro current J_Q (blind labels):\n" + "\n".join("      " + str(row) for row in b['JQ']))
        div = np.abs(b['JQ'].sum(1)).max(); print(f"    divergence of J_Q (row sums) {div:.1e}")
        Gs = adjoint(G, pi); bs = blind_pipeline(Gs, K); same, mp = canon_perm(r['lab'], bs['rec']['lab'])
        if same:
            Js = np.array([[bs['JQ'][mp[i], mp[j]] for j in range(K)] for i in range(K)])
            print(f"    COVARIANCE (G -> G*): same partition True; max|J_Q* + J_Q| {np.abs(Js + b['JQ']).max():.1e}")
        else:
            print("    COVARIANCE (G -> G*): partitions differ")
        mp_h = {bl: int(np.bincount(lab[r['lab'] == bl], weights=pi[r['lab'] == bl]).argmax()) for bl in np.unique(r['lab'])}
        Jm = current(G, pi); JB = np.array([[Jm[np.ix_(lab == i, lab == j)].sum() for j in range(K)] for i in range(K)])
        JQh = np.zeros((K, K))
        for i in range(K):
            for j in range(K):
                JQh[mp_h[i], mp_h[j]] = b['JQ'][i, j]
        print(f"    AUDIT (hidden labels, after blind): misclassified {audit(r['lab'], lab, pi):.2e}; max|J_Q - J_B| {np.abs(JQh - JB).max():.1e};"
              f" hidden J_B(0,1) {JB[0, 1]:.3e}, J_B(0,2) {JB[0, 2]:.3e}, J_B(0,3) {JB[0, 3]:.3e}, J_B(1,2) {JB[1, 2]:.3e},"
              f" J_B(1,3) {JB[1, 3]:.3e}, J_B(2,3) {JB[2, 3]:.3e}")


if __name__ == "__main__":
    n4_audit("R1a", R1a, [20, 40, 80], 3)
    n4_audit("R2", R2, [10, 20, 40, 80], 4)
    n4_audit("R4", R4, [10, 20, 40, 80], 3)
    n4_audit("R5", R5, [10, 20, 40, 80], 4)
    run_R5()
