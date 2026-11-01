"""m1_vacuum SKEPTIC, part 3: fixed-charge radial Hessian of the Q = 200000 M1 object, built differently from both
Codex (3-point Laplacian on u = rho dF, w = rho dS, Dirichlet nodes, dense eigh of L+ + K) and the track (same
operator, own profile).  Independent code (Claude skeptic, 1 Nov 2026).

Here the second variation is taken DIRECTLY of the discrete cell-centred energy functional of sk_qball.py,
    E_Q[F, S] = 4 pi { sum_faces (1/2) rho_f^2 (dF/h)^2 h + ... + sum_cells rho_i^2 h [U(F_i,S_i)] } + Qbar^2 / (2 I[F]),
    I[F] = 4 pi sum_i rho_i^2 h F_i^2,
as a generalised symmetric eigenproblem  A v = lambda M v,  M = diag(rho_i^2 h)  (the L^2(rho^2 d rho) norm, so lambda
is the same quantity as the continuum L+ eigenvalue; u = rho dF is an isometry).  Checks:
  (1) the discrete Hessian A + K is verified against second differences of E_Q itself in random directions
      (this tests the rank-one charge term K = 4 Omega^2 |MF><MF| / <F,MF> without relying on any derivation);
  (2) L+ alone: count of negative eigenvalues (expect exactly 1);
  (3) the secular-equation identity f(0) = 1 + c w^T A^{-1} w = (1/I) dQbar/dOmega: the lowest eigenvalue of L+ + K is
      positive iff dQbar/dOmega < 0 (Vakhitov-Kolokolov) -- computed from a sparse solve and compared with the branch;
  (4) the lowest L+ + K eigenvalue vs box size L and vs the lowest eigenvalue of the FREE operator -d^2 + 1 - Omega^2
      on the same grid (object removed): the box-mode interpretation of the track;
  (5) L- (phase) and ell = 1 (translation) zero modes.
"""
import json, math, sys, time
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import splu
from scipy.linalg import eigh

sys.path.insert(0, sys.argv[1])
import sk_qball as q
a, b, c = q.a, q.b, q.c
QBAR = 2000.0


def build(sv, F, S, Om):
    """returns M (diag), A_plus (2N x 2N: F,S blocks, no K), K (rank-one vector w and coefficient), A_minus (N x N),
    gradient block G (N x N) and the energy functional E_Q for the second-difference test"""
    N, h, r2, Af = sv.N, sv.h, sv.r2, sv.Af
    M = r2 * h
    # gradient quadratic form (1/2) sum_j Af_j (X_j - X_{j-1})^2 / h, faces j = 1..N, ghost X_N = -X_{N-1}
    diag = (Af[:-1] + Af[1:]) / h
    diag[-1] = (Af[-2] + 4 * Af[-1]) / h
    off = -Af[1:-1] / h
    G = sps.diags([off, diag, off], [-1, 0, 1], format='csr')
    UFF = S * S + 3 * a * F * F
    UFS = 2 * S * F
    USS = q.Uphi_pp(S) + F * F
    Ap = sps.bmat([[G + sps.diags(M * (UFF - Om * Om)), sps.diags(M * UFS)],
                   [sps.diags(M * UFS), G + sps.diags(M * USS)]], format='csr')
    Am = G + sps.diags(M * (S * S + a * F * F - Om * Om))
    I = np.sum(M * F * F)
    w = np.r_[M * F, np.zeros(N)]
    cK = 4 * Om * Om / I
    def E_Q(Fv, Sv):
        Fg = np.r_[Fv, -Fv[-1]]; Sg = np.r_[Sv, 2 - Sv[-1]]
        Eg = 0.5 * np.sum(Af[1:] * ((Fg[1:] - Fg[:-1])**2 + (Sg[1:] - Sg[:-1])**2)) / h
        EU = np.sum(M * (q.Uphi(Sv) + 0.5 * Sv * Sv * Fv * Fv + a / 4 * Fv**4))
        Iv = np.sum(M * Fv * Fv)
        return Eg + EU + (QBAR / (4 * math.pi))**2 / (2 * Iv)      # energy / 4 pi at fixed Qbar (Omega eliminated)
    return M, Ap, Am, G, w, cK, I, E_Q


def main():
    out = {'boxes': []}
    t0 = time.time()
    # reference object at several box sizes (h = 0.1) and h = 0.05 at L = 80; the object is re-solved on each grid
    prev = None
    for (L, h) in [(40.0, 0.1), (80.0, 0.1), (120.0, 0.1), (200.0, 0.1), (80.0, 0.05)]:
        sv = q.CC(L, h)
        if prev is None:
            F, S = sv.guess(); Om = 0.8817
        else:
            F, S = sv.interp(prev[0], prev[1], prev[2]); Om = prev[3]
        F, S, Om, ob = q.fixed_Q(sv, QBAR, Om, F, S)
        prev = (sv, F, S, Om)
        N = sv.N
        M, Ap, Am, G, w, cK, I, E_Q = build(sv, F, S, Om)
        Md = np.diag(np.r_[M, M])
        ApK = Ap.toarray() + cK * np.outer(w, w)
        ev_p = eigh(Ap.toarray(), Md, eigvals_only=True, subset_by_index=[0, 3])
        ev_pk, vec_pk = eigh(ApK, Md, subset_by_index=[0, 3])
        ev_m = eigh(Am.toarray(), np.diag(M), eigvals_only=True, subset_by_index=[0, 2])
        # ell = 1: + 2/rho^2 on both blocks (translation zero mode)
        A1 = Ap.toarray() + np.diag(np.r_[M * 2 / sv.r2, M * 2 / sv.r2])
        ev_1 = eigh(A1, Md, eigvals_only=True, subset_by_index=[0, 2])
        # free operator (object removed): -d^2 + (1 - Omega^2), same grid, F block only
        Afree = (G + sps.diags(M * (1 - Om * Om))).toarray()
        ev_free = eigh(Afree, np.diag(M), eigvals_only=True, subset_by_index=[0, 2])
        n_neg_p = int(np.sum(eigh(Ap.toarray(), Md, eigvals_only=True) < 0))
        edge = 1 - Om * Om
        n_below_edge_pk = int(np.sum(eigh(ApK, Md, eigvals_only=True) < edge - 1e-9))
        # secular identity: f(0) = 1 + cK w^T Ap^{-1} w
        y = splu(Ap.tocsc()).solve(w)
        f0 = 1 + cK * float(w @ y)
        # dQbar/dOmega from the branch (two fixed-Omega solves at +-dOm)
        dOm = 2e-4
        Qs = []
        for s in (-1, 1):
            Fs_, Ss_, _, _ = sv.solve(F, S, Om + s * dOm); Qs.append(sv.observables(Fs_, Ss_, Om + s * dOm)['Qbar'])
        dQ_dOm = (Qs[1] - Qs[0]) / (2 * dOm)
        f0_from_branch = dQ_dOm / (4 * math.pi * I)        # (1/I) dQbar/dOmega with I = 4 pi int rho^2 F^2
        # second-difference test of the Hessian (random directions, fixed Qbar)
        rng = np.random.default_rng(3)
        errs = []
        for _ in range(4):
            eta = rng.standard_normal(2 * N) * np.exp(-sv.r / 8.0).repeat(1)[np.r_[np.arange(N), np.arange(N)]]
            t = 1e-4
            d2 = (E_Q(F + t * eta[:N], S + t * eta[N:]) - 2 * E_Q(F, S) + E_Q(F - t * eta[:N], S - t * eta[N:])) / t**2
            quad = float(eta @ (ApK @ eta))
            errs.append(abs(d2 - quad) / abs(quad))
        # localisation of the lowest L+ + K mode
        v0 = vec_pk[:, 0]
        inside = float(np.sum(Md.diagonal() * v0**2 * (np.r_[sv.r, sv.r] < 15)) / np.sum(Md.diagonal() * v0**2))
        row = dict(L=L, h=h, Omega=Om, edge_1_minus_Om2=edge,
                   Lplus_lowest=ev_p.tolist(), Lplus_negative_count=n_neg_p,
                   LplusK_lowest=ev_pk.tolist(), LplusK_count_below_edge=n_below_edge_pk,
                   free_box_lowest=ev_free.tolist(), LplusK_minus_free=float(ev_pk[0] - ev_free[0]),
                   Lminus_lowest=ev_m.tolist(), ell1_lowest=ev_1.tolist(),
                   secular_f0=f0, secular_f0_from_dQdOmega=f0_from_branch, dQbar_dOmega=dQ_dOm,
                   hessian_second_difference_relative_errors=errs,
                   lowest_LplusK_mode_fraction_inside_rho15=inside)
        out['boxes'].append(row)
        print('L=%4.0f h=%.2f | L+ min %.5f (neg count %d) | L+ +K lowest %.5f %.5f (below edge: %d) | free box %.5f | edge %.5f | L- %.1e | ell1 %.1e | f(0) %.4f vs %.4f | hess test %.1e | inside15 %.3f  t=%.0fs' % (
            L, h, ev_p[0], n_neg_p, ev_pk[0], ev_pk[1], n_below_edge_pk, ev_free[0], edge, ev_m[0], ev_1[0], f0, f0_from_branch, max(errs), inside, time.time() - t0), flush=True)
    out['codex'] = dict(L80_h0p1=0.2243399467102812, L120_h0p1=0.2233198384209463, L80_h0p2=0.2243396857099204,
                        phase_zero_h0p1=-2.7709672563027838e-05, ell1_zero_h0p1=-0.00012476784885012988)
    out['track'] = dict(Lplus_negative=-0.09515, lowest_LpK={'40': 0.23078, '80': 0.22434, '120': 0.22332, '200': 0.22283}, edge=0.22257)
    out['runtime_s'] = time.time() - t0
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
