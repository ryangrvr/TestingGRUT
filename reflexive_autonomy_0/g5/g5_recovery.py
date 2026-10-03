"""RA0 / G5 -- metastable blind-recovery theorem: NUMERICAL AUDIT of the bounds.  Independent code path, not reviewer.

Firewall by code structure:
  * `recover(V, pi)` is the G5 recovery map.  It sees only the slow space V (pi-orthonormal basis) and pi.  It computes
    the cubic moment tensor T_V, ALL idempotents of the compressed product a*b = E_V(ab) (i.e. all solutions of
    T_V(c,c) = c), classifies primitive ones by the spectral count "#eig(T_V(c,.,.)) > 1/2 == 1", and rounds by argmax.
    Random Newton starts are a numerical solver for a polynomial system whose solution set is fixed by (V, pi); the
    theory (G5_METASTABLE_RECOVERY_THEOREM.md) does not use them.
  * `proof_side(...)` uses the hidden partition to evaluate the theorem's quantities (s, W, epsilon, bounds).  It is
    called only after `recover` has returned.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'g4'))
from g4_algebra import F1, F2, from_conductances, spectrum, gap_metric, partition_basis, christoffel, thm_T, RNG_COMPUTE

np.set_printoptions(precision=4, suppress=True)


def cubic(V, pi):
    return np.einsum('x,xa,xb,xc->abc', pi, V, V, V)


def all_idempotents(T, starts=4000, radius=2.5):
    """All solutions of T(c,c) = c by Newton from many starts (numerical solver; tolerances: residual 1e-12,
    de-duplication 1e-8)."""
    k = T.shape[0]; rng = np.random.default_rng(RNG_COMPUTE); sols = [np.zeros(k)]
    for _ in range(starts):
        c = rng.standard_normal(k); c *= radius * rng.random() / np.linalg.norm(c)
        for _ in range(100):
            F = np.einsum('abc,b,c->a', T, c, c) - c
            if np.linalg.norm(F) < 1e-13:
                break
            J = 2 * np.einsum('abc,b->ac', T, c) - np.eye(k)
            try:
                c = c - np.linalg.solve(J, F)
            except np.linalg.LinAlgError:
                break
            if np.linalg.norm(c) > 1e6:
                break
        if np.linalg.norm(np.einsum('abc,b,c->a', T, c, c) - c) < 1e-12 and not any(np.linalg.norm(c - s) < 1e-8 for s in sols):
            sols.append(c)
    return sols


def recover(V, pi):
    """G5 recovery map: V -> T_V -> idempotents -> primitive idempotents -> argmax partition.  No labels, no K."""
    k = V.shape[1]; T = cubic(V, pi); sols = all_idempotents(T)
    prim = [c for c in sols if np.sum(np.linalg.eigvalsh(np.einsum('abc,a->bc', T, c)) > 0.5) == 1]
    E = np.column_stack([V @ c for c in prim]) if prim else np.zeros((len(pi), 0))
    out = dict(k=k, n_idem=len(sols), n_prim=len(prim), E=E, C=prim)
    if len(prim) != k:
        out['lab'] = None; return out                           # recovery map undefined: refuse
    srt = np.sort(E, 1)
    out['ties'] = float(pi[srt[:, -1] == srt[:, -2]].sum())
    out['lab'] = np.argmax(E, 1)
    out['gap'] = gap_metric(V, partition_basis(out['lab'], pi), pi)
    return out


def inj_norm_sym(E, starts=400):
    """Injective norm of a symmetric 3-tensor = max_{|a|=1} |E(a,a,a)| (Banach).  Projected-gradient ascent."""
    k = E.shape[0]; rng = np.random.default_rng(RNG_COMPUTE); best = 0.0
    for _ in range(starts):
        a = rng.standard_normal(k); a /= np.linalg.norm(a)
        for _ in range(300):
            g = np.einsum('abc,b,c->a', E, a, a); val = a @ g
            a2 = np.sign(val) * g if np.linalg.norm(g) > 0 else a
            a2 = a + 0.5 * (a2 / max(np.linalg.norm(a2), 1e-300) - a); a2 /= np.linalg.norm(a2)
            if np.linalg.norm(a2 - a) < 1e-14:
                break
            a = a2
        best = max(best, abs(np.einsum('abc,a,b,c->', E, a, a, a)))
    return best


def proof_side(V, pi, lab, rec, C=None, nu=None):
    """Theorem quantities.  Uses the hidden partition `lab` (proof only)."""
    K = len(np.unique(lab)); Ab = partition_basis(lab, pi)
    p = np.array([pi[lab == b].sum() for b in np.unique(lab)]); Lam = p.min() ** -0.5
    s = gap_metric(V, Ab, pi)
    M = V.T @ (pi[:, None] * Ab); U, _, Wt = np.linalg.svd(M); Q = U @ Wt          # polar factor of E_V E_A
    WA = V @ Q                                                                      # W applied to the A-basis
    WmI = float(np.linalg.norm((WA - Ab) * np.sqrt(pi)[:, None], 2))
    Tt = cubic(WA, pi); TA = np.zeros((K, K, K)); TA[np.arange(K), np.arange(K), np.arange(K)] = p ** -0.5
    eps = inj_norm_sym(Tt - TA); CV = christoffel(V)
    eps_bd = s ** 2 * (11 * Lam + 2 * CV)        # LEMMA G5-2: first-order terms cancel because A is an algebra
    out = dict(K=K, Lam=Lam, s=s, WmI=WmI, eps=eps, eps_bd=eps_bd, eps_loc=1 / (8 * Lam),
               eps_glob=1 / (36 * K ** 1.5 * Lam ** 2), CV=CV, s2CV=s ** 2 * CV)
    if rec['lab'] is not None:
        E = rec['E']; err2, sup = [], []
        for i in range(K):                                                          # match by proximity (proof side)
            ind = (lab == np.unique(lab)[i]).astype(float)
            d = [np.sqrt(np.sum(pi * (E[:, j] - ind) ** 2)) for j in range(E.shape[1])]
            j = int(np.argmin(d)); err2.append(d[j]); sup.append(np.max(np.abs(E[:, j] - ind)))
        out['idem_err'] = max(err2); out['idem_err_bd'] = (8 / 3) * eps + np.sqrt(2) * s * np.sqrt(p.max())
        out['sup_err'] = max(sup)
        mis = 0.0
        for b in np.unique(rec['lab']):
            sel = rec['lab'] == b; vals = np.unique(lab[sel])
            maj = vals[np.argmax([pi[sel & (lab == v)].sum() for v in vals])]; mis += pi[sel & (lab != maj)].sum()
        out['mis'] = float(mis); out['mis_bd'] = (256 / 9) * K * eps ** 2 + 8 * s ** 2
        out['mis_bd_from_s'] = (256 / 9) * K * eps_bd ** 2 + 8 * s ** 2
    if C is not None:
        t = thm_T(C, nu, lab, pi, V); out['dk'] = t['dk']; out['eta_g'] = t['eta'] / t['g']
    return out


def run_F1():
    print("== G5.9 main control F1 (K=3 hidden; rank 3 certified by G4-T) ==")
    print("   blind:  #idempotents (theory: 2^k = 8 when eps < eps_glob), #primitive (theory: k), tie mass, gap(V,A_blind)")
    print("   proof:  s=||E_V-E_A|| vs DK bound; ||W-I|| vs sqrt2*s; eps=||T~-T_A||_inj vs bound(s); local/global thresholds;")
    print("           idempotent L2 error vs bound; sup error (margin diagnostic); misclassified mass vs bound (eps) and (eps_bd)")
    for M in [5, 10, 20, 40, 80, 160]:
        C, nu, lab = F1(M); Lsym, pi = from_conductances(C, nu); lam, U = spectrum(Lsym, pi); V = U[:, :3]
        r = recover(V, pi)
        g = f"ties {r['ties']:.1e}, gap {r['gap']:.3e}" if r['lab'] is not None else "REFUSE (primitive count != k)"
        print(f"  M={M:3d} n={len(pi)}: BLIND idempotents {r['n_idem']}, primitive {r['n_prim']}; {g}")
        q = proof_side(V, pi, lab, r, C, nu)
        print(f"     PROOF s {q['s']:.3e} (DK {q['dk']:.2e}, eta/g {q['eta_g']:.2e}); ||W-I|| {q['WmI']:.3e} <= {np.sqrt(2) * q['s']:.3e};"
              f" Lam {q['Lam']:.2f}, C_V {q['CV']:.2f}, s^2 C_V {q['s2CV']:.2e}")
        print(f"     PROOF eps {q['eps']:.3e} <= eps_bd {q['eps_bd']:.3e};  eps_loc {q['eps_loc']:.3e}, eps_glob {q['eps_glob']:.2e}"
              f"  -> theorem hypotheses met at this M: {q['eps'] < q['eps_glob']} (with true eps), {q['eps_bd'] < q['eps_glob']} (with bound)")
        if 'mis' in q:
            print(f"     PROOF idempotent L2 err {q['idem_err']:.3e} <= {q['idem_err_bd']:.3e};  sup err {q['sup_err']:.3f} (exact iff < 1/2 suffices);"
                  f"  misclassified {q['mis']:.2e} <= {q['mis_bd']:.2e} (eps) / {q['mis_bd_from_s']:.2e} (eps_bd)")


def run_rank2():
    print("== G5.8 rank-2 reduction: G5 map vs PROP G4-P rounding {f > s/2} on F2 coarse level (any V of dim 2) ==")
    for M in [5, 20, 80]:
        C, nu, lab, sup = F2(M); Lsym, pi = from_conductances(C, nu); lam, U = spectrum(Lsym, pi); V = U[:, :2]
        r = recover(V, pi); f = V[:, 1]; s3 = float(np.sum(pi * f ** 3))
        labP = (f > s3 / 2).astype(int)
        same = r['lab'] is not None and (np.all(labP == r['lab']) or np.all(labP == 1 - r['lab']))
        # closed form idempotents e_+- = (1 -+ s/sqrt(s^2+4))/2 +- f/sqrt(s^2+4)
        d = np.sqrt(s3 ** 2 + 4); ep = (1 - s3 / d) / 2 + f / d
        cf = min(np.max(np.abs(r['E'][:, j] - ep)) for j in range(r['E'].shape[1]))
        print(f"  M={M:3d}: idempotents {r['n_idem']} (theory 4 for every 2-dim V), primitive {r['n_prim']};"
              f" partition == G4-P rounding: {same};  max|e_+ - closed form| {cf:.1e}")


def run_refusal():
    print("== refusal check (stress, not a G5 output): G5 map on forced non-metastable slow spaces ==")
    from g4_algebra import ssep, from_kernel, ladder
    for tag, (Lsym, pi) in [("SSEP L=10", from_kernel(*ssep(10))),
                            ("ladder beta=0 L=50 (mod)", from_conductances(*ladder(50, 0, (1.0, 1.6), True)[:2]))]:
        lam, U = spectrum(Lsym, pi); V = U[:, :3]; r = recover(V, pi)
        print(f"  {tag}: forced rank 3: idempotents found {r['n_idem']}, primitive {r['n_prim']} -> "
              f"{'REFUSE' if r['lab'] is None else 'partition with gap %.3f' % r['gap']}")


if __name__ == "__main__":
    run_F1(); run_rank2(); run_refusal()
