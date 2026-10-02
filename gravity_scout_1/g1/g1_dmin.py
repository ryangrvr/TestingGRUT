"""GRAVITY-SCOUT-1 / G1 -- minimum admissible splitting distance d_min.

Independent code path, not independent reviewer. Lattice runs are ILLUSTRATIONS of continuum statements.

Part A  (flat control, 1+1): exact minimum energy of a PRODUCT state across a collar of d sites,
        E_prod(d) = min { <H>_w - <H>_vac : w restricted to A(left) v A(right) is a product state }.
        Gaussian free scalar H = 1/2 sum p^2 + 1/2 phi.K.phi.  The minimum over ALL states (mixed,
        non-Gaussian) equals an SDP over second moments: Gaussianising a product state keeps its energy
        and keeps the A|B cross blocks zero; time-reversal averaging removes x-p blocks; first moments
        only add energy.  So  min 1/2 tr P + 1/2 tr K X  s.t. [[X, I/2],[I/2, P]] >= 0, X_AB = P_AB = 0.
        Continuum check: fixed physical collar D and mass M, lattice spacing a halved.
Part B  (3+1 planar collar, massless scalar): transverse Fourier modes decouple (transverse-translation
        averaging preserves feasibility and energy), each mode is a 1+1 field of mass k.  Per unit area
            E/A = (1/2pi) int k e(d,k) dk ,   alpha_3 := d^3 E/A   (units hbar c),
        and the vacuum mutual information across the collar per unit area, kappa_I := d^2 I/A.
        Longitudinal direction is a lattice (spacing 1); d in sites.  Convergence in d = 2,3,4,6 is the
        continuum check.
The gravity read-off (Part C) is analytic and lives in G1_DMIN_RESULT.md.
"""
import numpy as np, cvxpy as cp, sys, time


def kmat(N, mu):
    return (mu**2 + 2)*np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)


def vac(K):
    lam, V = np.linalg.eigh(K)
    return 0.5*(V*lam**-0.5) @ V.T, 0.5*(V*lam**0.5) @ V.T


def emin(args):
    L, d, mu = args
    N = 2*L + d; K = kmat(N, mu); A = np.arange(L); B = np.arange(L + d, N)
    X = cp.Variable((N, N), symmetric=True); P = cp.Variable((N, N), symmetric=True)
    I = np.eye(N)/2
    cons = [cp.bmat([[X, I], [I, P]]) >> 0, X[np.ix_(A, B)] == 0, P[np.ix_(A, B)] == 0]
    pr = cp.Problem(cp.Minimize(0.5*cp.trace(P) + 0.5*cp.trace(K @ X)), cons)
    pr.solve(solver="CLARABEL")
    Xv, Pv = vac(K); Ev = 0.5*np.trace(Pv) + 0.5*np.trace(K @ Xv)
    # excess energy-density profile of the minimiser (site energies)
    h = 0.5*np.diag(P.value - Pv) + 0.5*np.diag(K @ (X.value - Xv))
    c = L + (d - 1)/2; w = 3*max(d, 1)
    frac = h[np.abs(np.arange(N) - c) <= w].sum()/h.sum()
    return pr.value - Ev, frac, pr.status


def S(X, P, idx):
    nu = np.sqrt(np.clip(np.linalg.eigvals(X[np.ix_(idx, idx)] @ P[np.ix_(idx, idx)]).real, 0.25, None))
    nu = nu[nu > 0.5 + 1e-12]
    return float(np.sum((nu + .5)*np.log(nu + .5) - (nu - .5)*np.log(nu - .5)))


def mi(L, d, mu):
    N = 2*L + d; X, P = vac(kmat(N, mu)); A = np.arange(L); B = np.arange(L + d, N)
    return S(X, P, A) + S(X, P, B) - S(X, P, np.concatenate([A, B]))


if __name__ == "__main__":
    # usage: python3 g1_dmin.py A        (flat 1+1 control)
    #        python3 g1_dmin.py B <d>    (3+1 planar collar of d sites)
    # Run serially; g1_dmin.log is the concatenation of the runs A, B 2, B 3, B 4, B 6.
    t0 = time.time()
    if sys.argv[1] == "A":
        print("== Part A: 1+1 flat control; physical collar D=1, mass M=0.5, region L_phys=6 per side ==")
        jobs = [(int(6/a), int(round(1/a)), 0.5*a) for a in [1/2, 1/3, 1/4]]
        jobs += [(int(6/a), 0, 0.5*a) for a in [1/2, 1/3, 1/4]]   # sharp cut d = 0 (non-split)
        for (L, d, mu) in jobs:
            e, f, st = emin((L, d, mu)); a = mu/0.5
            print(f"a={a:.4f}  d_sites={d:2d}  D_phys={d*a:.3f}  E_prod_phys={e/a:.5f}  "
                  f"(lattice {e:.6f})  excess_within_3d={f:.3f}  [{st}]", flush=True)
    else:
        d = int(sys.argv[2]); L = 20
        if d == 2:
            print(f"== Part B: 3+1 planar collar via transverse modes, massless scalar, L={L} sites per side ==")
        us = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0])   # u = k d
        res = [emin((L, d, u/d)) for u in us]
        f = np.array([d*r[0] for r in res])                        # f(u) = d * e(d, k=u/d)
        # solver noise floor ~1e-6: modes with f below FLOOR are dropped and replaced by an exponential tail
        FLOOR = 1e-5
        ok = f > FLOOR; uo, g = us[ok], (us*f)[ok]
        head = 0.5*uo[0]*g[0]                                      # int_0^u0 u f du  (u f ~ linear)
        rate = np.log(g[-2]/g[-1])/(uo[-1] - uo[-2])
        tail = g[-1]/rate                                          # exponential tail beyond the last kept mode
        alpha3 = (np.trapezoid(g, uo) + head + tail)/(2*np.pi)
        fr = [r[1] for r, x in zip(res, f) if x > 1e-3]            # localisation diagnostic, modes above 1e-3 only
        m = np.array([mi(L, d, u/d) for u in us])                  # per-mode MI; I/A = (1/2pi d^2) int u m du
        kap = (np.trapezoid(us*m, us) + 0.5*us[0]*us[0]*m[0])/(2*np.pi)
        print(f"d={d}  f(u)=d*e(d,u/d): " + " ".join(f"{x:.4g}" for x in f))
        print(f"      mode-MI m(u):     " + " ".join(f"{x:.4g}" for x in m))
        print(f"      alpha_3 = d^3 E/A = {alpha3:.5f}  (head {head/(2*np.pi):.1e}, tail {tail/(2*np.pi):.1e})"
              f"   kappa_I = d^2 I/A = {kap:.5f}   min excess-energy fraction within 3d (modes f>1e-3): "
              f"{min(fr):.3f}   statuses: {sorted(set(r[2] for r in res))}")
    print(f"[runtime {time.time()-t0:.0f}s]")
