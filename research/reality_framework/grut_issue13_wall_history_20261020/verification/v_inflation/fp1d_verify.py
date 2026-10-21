"""Independent 1D Starobinsky-Yokoyama Fokker-Planck check for the phi spectator (verifier v_inflation).

Written from scratch; no module code imported or copied.

Convention (stated): metric (-,+,+,+), V_R = (1/2) xi_phi R phi^2, R = 12 H^2 in de Sitter, so
  V(phi) = (1/2) m2(H) phi^2 + (lam(H)/4) phi^4 + 6 xi(H) H^2 phi^2  (+ const),
with m2, lam, xi the one-loop MSbar couplings run to mu = H by rg_xi.py (m2(mu0) = -lam v^2).
Units y = phi/H, w = V/H^4, time N = e-folds.  Langevin dy/dN = -w'(y)/3 + eta/(2 pi), <eta eta> = delta(N-N').
Backward Kolmogorov generator  L f = A f' + Dc f'',  A = -w'/3,  Dc = 1/(8 pi^2).

Fresh start at +v: f(y,N) = 1 - E_y[sgn phi_N] obeys df/dN = L f on y > 0 with f(0,N) = 1 (oddness of
E[sgn]), f'(Y) = 0 (reflecting), f(y,0) = 0.  p_- = f(v/H, N)/2.

Discretization: second-order central differences on a uniform grid with y0 = v/H on a node; Peclet number
checked < 1 so that I - dt L is an M-matrix; implicit Euler in N (every solve sums positive terms only, so tiny
p_- keep relative precision); Richardson extrapolation in dt and in h.

Odd-sector eigenvalue Lambda_1: smallest eigenvalue of -L with Dirichlet at y = 0 (killed half chain),
by power iteration on the explicit birth-death Green function (cumulative positive sums; valid at any barrier),
cross-checked by a dense symmetric tridiagonal eigensolver where Lambda_1 > 1e-11.
"""
import json, math, sys, time
import numpy as np
from scipy.linalg import solve_banded, eigh_tridiagonal
from scipy.optimize import brentq

sys.path.insert(0, '.')
import rg_xi

DC = 1.0 / (8.0 * math.pi ** 2)
V_GEV = 1000.0
_cache = {}


def params(H):
    key = float(H)
    if key not in _cache:
        r = rg_xi.run_to(np.array([max(H, rg_xi.mu0)]))
        c = rg_xi.couplings(*r[float(max(H, rg_xi.mu0))])
        _cache[key] = c
    return _cache[key]


def wfun(H, use_running=True, xi_on=True):
    c = params(H) if use_running else dict(lam_phi=0.01, m2_phi=-0.01 * V_GEV ** 2, xi_phi=0.0)
    a2 = c['m2_phi'] / H ** 2 + (12.0 * c['xi_phi'] if xi_on else 0.0)  # coefficient of y in w'
    lam = c['lam_phi']
    w = lambda y: 0.5 * a2 * y ** 2 + 0.25 * lam * y ** 4
    wp = lambda y: a2 * y + lam * y ** 3
    return w, wp, a2, lam


def grid(H, hgoal, **kw):
    w, wp, a2, lam = wfun(H, **kw)
    y0 = V_GEV / H
    k = max(4, int(round(y0 / hgoal)))
    h = y0 / k
    wmin = min(w(0.0), w(math.sqrt(max(-a2, 0) / lam)) if a2 < 0 else 0.0)
    # upper edge where Boltzmann weight relative to the minimum is < e^-80
    Y = max(1.05 * y0, 1.0)
    while (8 * math.pi ** 2 / 3) * (w(Y) - wmin) < 80.0:
        Y *= 1.05
    n = int(math.ceil(Y / h))
    y = h * np.arange(n + 1)
    A = -wp(y) / 3.0
    pe = np.max(np.abs(A[1:]) * h / (2 * DC))
    if pe >= 1.0:
        raise ValueError('Peclet %.3f >= 1, refine h' % pe)
    return y, h, k, A, pe


def p_minus_fresh(H, Nlist, hgoal=0.005, dt=0.01, **kw):
    """returns p_-(N) for N in Nlist (implicit Euler, one grid, one dt)."""
    y, h, k, A, pe = grid(H, hgoal, **kw)
    n = len(y) - 1
    up = DC / h ** 2 + A / (2 * h)
    lo = DC / h ** 2 - A / (2 * h)
    # unknowns nodes 1..n
    diag = np.full(n, 1.0 + dt * 2 * DC / h ** 2)
    sup = -dt * up[1:n]          # (i, i+1) for i = 1..n-1
    sub = -dt * lo[2:n + 1]      # (i, i-1) for i = 2..n
    sub = sub.copy(); sub[-1] = -dt * 2 * DC / h ** 2  # reflecting ghost at node n: lower coefficient doubles
    ab = np.zeros((3, n))
    ab[0, 1:] = sup
    ab[1, :] = diag
    ab[2, :-1] = sub
    f = np.zeros(n)
    bc = dt * lo[1] * 1.0  # f_0 = 1
    Ns = sorted(Nlist)
    out = {}
    steps_done = 0
    for N in Ns:
        nst = int(round(N / dt))
        for _ in range(nst - steps_done):
            r = f.copy(); r[0] += bc
            f = solve_banded((1, 1), ab, r, overwrite_ab=False, overwrite_b=True, check_finite=False)
        steps_done = nst
        out[N] = 0.5 * f[k - 1]  # node k is y0 (unknown index k-1)
    return out, dict(h=h, n=n, peclet=pe, Y=y[-1])


def p_minus_extrap(H, Nlist, hgoal=0.005, dt=0.02, **kw):
    """Richardson in dt (implicit Euler, first order) and in h (second order)."""
    res = {}
    for hg in (hgoal, hgoal / 2):
        a, _ = p_minus_fresh(H, Nlist, hg, dt, **kw)
        b, info = p_minus_fresh(H, Nlist, hg, dt / 2, **kw)
        res[hg] = {N: 2 * b[N] - a[N] for N in Nlist}
        res[('raw', hg)] = (a, b)
    out = {}
    for N in Nlist:
        ph, ph2 = res[hgoal][N], res[hgoal / 2][N]
        out[N] = dict(p=(4 * ph2 - ph) / 3, p_h=ph, p_h2=ph2,
                      p_raw_dt=res[('raw', hgoal / 2)][0][N], p_raw_dt2=res[('raw', hgoal / 2)][1][N])
    return out


def lambda1(H, hgoal=0.005, iters=8, **kw):
    """Odd-sector smallest eigenvalue by power iteration on the explicit Green function."""
    y, h, k, A, pe = grid(H, hgoal, **kw)
    n = len(y) - 1
    up = DC / h ** 2 + A / (2 * h)
    lo = DC / h ** 2 - A / (2 * h)
    up = up.copy(); lo = lo.copy()
    up[n] = 0.0; lo[n] = 2 * DC / h ** 2
    # detailed balance pi_{i+1} lo_{i+1} = pi_i up_i, i = 1..n-1 (log space)
    lpi = np.zeros(n + 1)
    for i in range(1, n):
        lpi[i + 1] = lpi[i] + math.log(up[i]) - math.log(lo[i + 1])
    lpi[1:] -= lpi[1:].max()
    pi = np.exp(lpi[1:])            # nodes 1..n
    lo_ = lo[1:]                     # killing/leftward rates at nodes 1..n
    def G(v):
        F = np.cumsum((pi * v)[::-1])[::-1]   # sum_{j>=k} pi_j v_j
        return np.cumsum(F / (pi * lo_))
    x = np.ones(n)
    for _ in range(iters):
        gx = G(x)
        ratio = x / gx
        lo_b, hi_b = ratio.min(), ratio.max()
        x = gx / gx.max()
    return 0.5 * (lo_b + hi_b), (lo_b, hi_b), dict(h=h, n=n)


def lambda1_eig(H, hgoal=0.005, **kw):
    """Dense symmetric tridiagonal check (valid only where Lambda_1 >> 1e-13)."""
    y, h, k, A, pe = grid(H, hgoal, **kw)
    n = len(y) - 1
    up = DC / h ** 2 + A / (2 * h)
    lo = DC / h ** 2 - A / (2 * h)
    up[n] = 0.0; lo[n] = 2 * DC / h ** 2
    d = up[1:] + lo[1:]
    e = -np.sqrt(up[1:n] * lo[2:n + 1])  # symmetrized off-diagonal
    ev = eigh_tridiagonal(d, e, eigvals_only=True, select='i', select_range=(0, 2))
    return ev


def B_HM(H, **kw):
    w, wp, a2, lam = wfun(H, **kw)
    if a2 >= 0:
        return 0.0
    ym2 = -a2 / lam
    return (8 * math.pi ** 2 / 3) * (w(0.0) - w(math.sqrt(ym2)))


def kramers_lambda1(H, **kw):
    w, wp, a2, lam = wfun(H, **kw)
    a = -a2
    return math.sqrt(2.0) * a / (3 * math.pi) * math.exp(-B_HM(H, **kw))


def pure_quartic_test(lam=0.01):
    """w = lam y^4/4: Lambda_1/sqrt(lam) (SY: 0.08892) and the even-sector Lambda_2."""
    # emulate with a2 = 0 by a custom grid
    h = 0.004
    Y = 1.0
    while (8 * math.pi ** 2 / 3) * 0.25 * lam * Y ** 4 < 80:
        Y *= 1.05
    n = int(math.ceil(Y / h)); y = h * np.arange(n + 1)
    A = -(lam * y ** 3) / 3
    up = DC / h ** 2 + A / (2 * h); lo = DC / h ** 2 - A / (2 * h)
    up[n] = 0.0; lo[n] = 2 * DC / h ** 2
    # odd sector: Dirichlet at 0
    d = up[1:] + lo[1:]; e = -np.sqrt(up[1:n] * lo[2:n + 1])
    odd = eigh_tridiagonal(d, e, eigvals_only=True, select='i', select_range=(0, 1))
    # even sector: Neumann at 0 (ghost f_{-1} = f_1)
    upE = up.copy(); upE[0] = 2 * DC / h ** 2
    dE = upE[:] + np.concatenate([[0.0], lo[1:]])
    eE = -np.sqrt(upE[0:n] * lo[1:n + 1])
    even = eigh_tridiagonal(dE, eE, eigvals_only=True, select='i', select_range=(0, 1))
    return dict(Lambda1_over_sqrtlam=odd[0] / math.sqrt(lam), Lambda3_over_sqrtlam=odd[1] / math.sqrt(lam),
                even0=even[0], Lambda2_over_sqrtlam=even[1] / math.sqrt(lam),
                two_Lambda2_over_sqrtlam=2 * even[1] / math.sqrt(lam))


def spectral_check(H, N, hgoal=0.005):
    """p_-(fresh) from the eigen-expansion of the odd sector (independent of time stepping)."""
    y, h, k, A, pe = grid(H, hgoal)
    n = len(y) - 1
    up = DC / h ** 2 + A / (2 * h); lo = DC / h ** 2 - A / (2 * h)
    up[n] = 0.0; lo[n] = 2 * DC / h ** 2
    lpi = np.zeros(n + 1)
    for i in range(1, n):
        lpi[i + 1] = lpi[i] + math.log(up[i]) - math.log(lo[i + 1])
    lpi[1:] -= lpi[1:].max()
    s = np.exp(0.5 * lpi[1:])     # sqrt(pi) nodes 1..n
    d = up[1:] + lo[1:]; e = -np.sqrt(up[1:n] * lo[2:n + 1])
    lam, V = eigh_tridiagonal(d, e)
    # u(y,N) = E[sgn] with u(y,0)=1 (y>0): in symmetric variables g = s*u ; g = sum c_k e^{-lam N} V_k
    c = V.T @ (s * 1.0)
    g0 = (V[k - 1, :] * c)
    twop = np.sum(g0 * (-np.expm1(-lam * N))) / s[k - 1]
    return 0.5 * twop, lam[:3]


if __name__ == '__main__':
    t0 = time.time()
    out = {}
    out['pure_quartic'] = pure_quartic_test()
    # p_- at requested H values and module table values
    Hs = [250.0, 251.19, 298.5, 335.0, 350.0, 355.0, 398.1, 501.2, 1000.0, 1607.0, 2000.0, 2037.0]
    tab = {}
    for H in Hs:
        r = p_minus_extrap(H, [50.0, 60.0], hgoal=0.005 if H < 600 else min(0.005, 0.005 * 600 / H))
        tab['%g' % H] = {('N%d' % int(N)): r[N] for N in r}
        tab['%g' % H]['B_HM'] = B_HM(H)
        print(H, {N: r[N]['p'] for N in r}, time.time() - t0, flush=True)
    out['p_minus_fresh'] = tab
    # spectral cross-check of time stepping
    sc = {}
    for H in (250.0, 350.0, 2000.0):
        p, l3 = spectral_check(H, 60.0, hgoal=0.0025 if H < 600 else 0.0015)
        sc['%g' % H] = dict(p_minus_N60_spectral=p, lowest_odd_eigs=list(l3))
    out['spectral_check'] = sc
    print(sc, flush=True)

    def solve_thr(N, target, lo_H, hi_H):
        f = lambda H: math.log(p_minus_extrap(H, [N], hgoal=0.005 if H < 600 else 0.005 * 600 / H)[N]['p']) - math.log(target)
        return brentq(f, lo_H, hi_H, xtol=0.05)
    thr = {}
    thr['N60_p_1e-3'] = solve_thr(60.0, 1e-3, 300.0, 420.0)
    thr['N50_p_1e-3'] = solve_thr(50.0, 1e-3, 300.0, 420.0)
    thr['N60_p_pc'] = solve_thr(60.0, 0.3116077, 1300.0, 2000.0)
    thr['N50_p_pc'] = solve_thr(50.0, 0.3116077, 1400.0, 2200.0)
    thr['N60_p_1e-10'] = solve_thr(60.0, 1e-10, 200.0, 280.0)
    out['thresholds_GeV'] = thr
    print(thr, time.time() - t0, flush=True)
    json.dump(out, open('results_fp1d.json', 'w'), indent=1, default=str)
    print('done', time.time() - t0)
