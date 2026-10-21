"""Independent 2D (phi, |X|) stochastic check of the X 'entropic' correction (verifier v_inflation).

Written from scratch; no module code imported or copied.  D is set to zero (its factor R_D = 1.005 is checked
in entropic.py).  Fields: y = phi/H (odd sector: Dirichlet at y = 0 by oddness), r = rho/H with
rho^2 = q1^2 + q2^2 = 2|X|^2 (two real X components, radial measure r dr).
w = V/H^4, V = (1/2)(m2_phi + 12 xi_phi H^2) phi^2 + lam_phi phi^4/4
              + (1/2)(m2_X + g^2 phi^2 + 12 xi_X H^2) rho^2 + lam_X rho^4/4,  couplings at mu = H (rg_xi.py).
Generator L f = -(1/3) grad w . grad f + Dc Laplacian_2D(q1,q2) f, reversible w.r.t. pi = r exp(-8 pi^2 w/3).
Discretization: cell-centred finite volume with face conductances c_face = pi(face) (Scharfetter-free,
second order), so pi_i k_ij = (Dc/h^2) c_face is symmetric.  Lambda_1 = lowest eigenvalue of the symmetrized
-L (shift-invert eigsh).  Fresh start (phi = v, X = 0): backward implicit Euler with splu, evaluated at the
cell (y0, r = h_r/2) and extrapolated to r = 0 using the r = 3h_r/2 cell (f even in r).
"""
import json, math, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, splu
import rg_xi

DC = 1.0 / (8 * math.pi ** 2)
V = 1000.0


def build(H, hy=0.04, hr=0.05, smax=45.0, g2scale=1.0):
    r_ = rg_xi.run_to(np.array([max(H, rg_xi.mu0)]))
    c = rg_xi.couplings(*r_[float(max(H, rg_xi.mu0))])
    a2 = c['m2_phi'] / H ** 2 + 12 * c['xi_phi']
    b2 = c['m2_X'] / H ** 2 + 12 * c['xi_X']
    lp, lx, g2 = c['lam_phi'], c['lam_X'], c['gX2'] * g2scale
    S = 8 * math.pi ** 2 / 3

    def s(y, r):
        return S * (0.5 * a2 * y ** 2 + 0.25 * lp * y ** 4 + 0.5 * (b2 + g2 * y ** 2) * r ** 2 + 0.25 * lx * r ** 4)
    smin = s(math.sqrt(max(-a2 / lp, 0.0)), 0.0)
    y0 = V / H
    # y grid: cell centres (m + 1/2) hy, choose hy so that y0 is a cell centre
    k = max(3, int(round(y0 / hy - 0.5)))
    hy = y0 / (k + 0.5)
    Ymax = y0 + 1.0
    while s(Ymax, 0.0) - smin < smax:
        Ymax += 0.25
    Rmax = 1.0
    while s(0.0, Rmax) - smin < smax:
        Rmax += 0.25
    ny = int(math.ceil(Ymax / hy)); nr = int(math.ceil(Rmax / hr))
    yc = (np.arange(ny) + 0.5) * hy; rc = (np.arange(nr) + 0.5) * hr
    Y, R = np.meshgrid(yc, rc, indexing='ij')
    lpi = np.log(R) - (s(Y, R) - smin)
    idx = np.arange(ny * nr).reshape(ny, nr)
    rows, cols, vals = [], [], []
    diag = np.zeros(ny * nr)
    # y-faces
    yf = (np.arange(ny - 1) + 1.0) * hy
    YF, RF = np.meshgrid(yf, rc, indexing='ij')
    lcf = np.log(RF) - (s(YF, RF) - smin)
    a = idx[:-1, :].ravel(); b = idx[1:, :].ravel(); lc = lcf.ravel()
    rows.append(a); cols.append(b); vals.append(lc)
    # r-faces
    rf = (np.arange(nr - 1) + 1.0) * hr
    YG, RG = np.meshgrid(yc, rf, indexing='ij')
    lcg = np.log(RG) - (s(YG, RG) - smin)
    a2_ = idx[:, :-1].ravel(); b2_ = idx[:, 1:].ravel()
    rows.append(a2_); cols.append(b2_); vals.append(lcg.ravel())
    A = np.concatenate(rows); Bc = np.concatenate(cols); LC = np.concatenate(vals)
    hA = np.concatenate([np.full(len(a), hy), np.full(len(a2_), hr)])
    LPI = lpi.ravel()
    # symmetrized off-diagonal S_ij = (DC/h^2) c_ij / sqrt(pi_i pi_j)
    off = (DC / hA ** 2) * np.exp(LC - 0.5 * (LPI[A] + LPI[Bc]))
    # rates k_ij = (DC/h^2) c_ij / pi_i
    kab = (DC / hA ** 2) * np.exp(LC - LPI[A]); kba = (DC / hA ** 2) * np.exp(LC - LPI[Bc])
    np.add.at(diag, A, kab); np.add.at(diag, Bc, kba)
    # killing at y = 0 face (odd sector, ghost f = -f): rate 2 * (DC/hy^2) c(0,r)/pi_0
    lc0 = np.log(rc) - (s(0.0, rc) - smin)
    kill = 2 * (DC / hy ** 2) * np.exp(lc0 - lpi[0, :])
    diag[idx[0, :]] += kill
    n = ny * nr
    Ssym = sp.coo_matrix((np.concatenate([off, off]), (np.concatenate([A, Bc]), np.concatenate([Bc, A]))), shape=(n, n))
    Ssym = (sp.diags(diag) - Ssym).tocsc()
    # backward generator L (non-symmetric): L f_i = sum_j k_ij (f_j - f_i) - kill_i f_i
    Lm = sp.coo_matrix((np.concatenate([kab, kba]), (np.concatenate([A, Bc]), np.concatenate([Bc, A]))), shape=(n, n))
    Lm = (Lm - sp.diags(diag)).tocsc()
    kill_vec = np.zeros(n); kill_vec[idx[0, :]] = kill
    return dict(S=Ssym, L=Lm, kill=kill_vec, idx=idx, k=k, hy=hy, hr=hr, ny=ny, nr=nr, B=None)


def lambda1_2d(H, **kw):
    m = build(H, **kw)
    ev = eigsh(m['S'], k=2, sigma=0.0, which='LM', return_eigenvectors=False)
    return float(np.sort(ev)[0]), m


def fresh_2d(H, Ns, dt=0.02, **kw):
    """f = 1 - E[sgn phi_N]; df/dN = L f + kill*(1) [the killed mass maps to f_ghost = 1 - (-1)... ]."""
    m = build(H, **kw)
    n = m['L'].shape[0]
    # u = E[sgn]: du/dN = L u (kill term encodes u_ghost = -u). f = 1 - u: df/dN = L f + kill * 1 * 2? check:
    # L u = sum k (u_j - u_i) - kill u_i; with u = 1 - f: L u = -sum k (f_j - f_i) - kill (1 - f_i) = -(L f) - kill
    # so df/dN = -du/dN = L f + kill.
    M = (sp.identity(n, format='csc') - dt * m['L']).tocsc()
    lu = splu(M)
    f = np.zeros(n)
    out = {}
    done = 0
    for N in sorted(Ns):
        for _ in range(int(round(N / dt)) - done):
            f = lu.solve(f + dt * m['kill'])
        done = int(round(N / dt))
        f1 = f[m['idx'][m['k'], 0]]; f2 = f[m['idx'][m['k'], 1]]
        # f even in r: f(r) ~ f0 + c r^2; r1 = hr/2, r2 = 3hr/2 -> f0 = (9 f1 - f2)/8
        out[N] = dict(p_r0=0.5 * (9 * f1 - f2) / 8, p_cell=0.5 * f1)
    return out


if __name__ == '__main__':
    t0 = time.time()
    res = {}
    import fp1d_verify as F1
    for H in (1000.0, 330.0):
        L2, m = lambda1_2d(H, g2scale=0.0)
        res['decoupled_test_%g' % H] = dict(Lambda1_2D_g0=L2, Lambda1_1D=F1.lambda1(H, hgoal=0.002)[0])
        print(res['decoupled_test_%g' % H], flush=True)
    for H in (1000.0, 500.0, 400.0, 330.0, 280.0, 250.0):
        L2, m = lambda1_2d(H)
        L2f, m2 = lambda1_2d(H, hy=0.03, hr=0.0375)
        res['%g' % H] = dict(Lambda1_2D=L2, Lambda1_2D_fine=L2f, grid=[m['ny'], m['nr']])
        print(H, res['%g' % H], time.time() - t0, flush=True)
    fr = {}
    for dt in (0.02, 0.01):
        o = fresh_2d(330.0, [50.0, 60.0], dt=dt)
        fr['dt=%g' % dt] = o
        print('fresh 330', dt, o, time.time() - t0, flush=True)
    res['fresh_330'] = fr
    json.dump(res, open('results_fp2d.json', 'w'), indent=1)
