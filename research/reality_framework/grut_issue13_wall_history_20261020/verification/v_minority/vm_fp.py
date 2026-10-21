"""v_minority: independent 1D stochastic (Starobinsky-Yokoyama) production rate dp_-/dn, relic integral and thresholds.

Written from scratch by the verifier; it does NOT import inflation_spectator or minority_relic code.
Physical inputs (same convention as the packet, X = D = 0 path):
  V(phi) = 1/2 (m_phi^2(mu) + 12 xi_phi(mu) H^2) phi^2 + lambda_phi(mu)/4 phi^4,  mu = H,
  one-loop running of all couplings from mu_0 = 100 GeV via the generic real-scalar tensor formulas
  16pi^2 beta(L_abcd) = L_abef L_efcd + 2 perms, 16pi^2 beta(m2_ab) = L_abef m2_ef, 16pi^2 beta(xi_ab) = L_abef (xi_ef - d_ef/6).
FP in y = phi/H: dP/dn = d/dy [ D e^{-s} d/dy (e^{s} P) ], D = 1/(8 pi^2), s = (8 pi^2/3) V/H^4.
Discretization: nodes y_i = i h (i = 1..M), absorbing node y_0 = 0, reflecting far end; EXACT 1D steady conductance
  c_{i+1/2} = D / ∫_{y_i}^{y_{i+1}} e^{s} dy  (Gauss-Legendre), flux J = -c (rho_{i+1} - rho_i), rho = e^{s} P.
  Time evolution by eigendecomposition of the symmetrized generator (exact in n).
dp_-/dn = (absorbed flux at y = 0)/2  (reflection principle: P(phi_n < 0) = P(hit 0 by n)/2).
Lambda_1 for large barriers: MFPT formula in log space.
"""
import json, math, os, sys, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import logsumexp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vm_cosmo as C

OUT = os.path.dirname(os.path.abspath(__file__))
Dc = 1 / (8 * math.pi ** 2)
v, mu0 = 1000.0, 100.0
P0 = dict(lp=0.01, g2=0.01, lX=0.001, h2=1e-4, lD=0.001, al=0.0, m2p=-0.01 * v ** 2, m2X=0.0, m2D=100.0 ** 2,
          xp=0.0, xX=0.0, xD=0.0)

# ---------------- running ----------------
def tensors(p):
    L = np.zeros((4, 4, 4, 4))
    def setsym(idx, val):
        import itertools
        for perm in set(itertools.permutations(idx)):
            L[perm] = val
    setsym((0, 0, 0, 0), 6 * p['lp'])
    setsym((1, 1, 1, 1), 6 * p['lX']); setsym((2, 2, 2, 2), 6 * p['lX'])
    setsym((1, 1, 2, 2), 2 * p['lX'])
    setsym((3, 3, 3, 3), 6 * p['lD'])
    setsym((0, 0, 1, 1), 2 * p['g2']); setsym((0, 0, 2, 2), 2 * p['g2'])
    setsym((0, 0, 3, 3), 2 * p['h2'])
    setsym((1, 1, 3, 3), 2 * p['al']); setsym((2, 2, 3, 3), 2 * p['al'])
    M2 = np.diag([p['m2p'], p['m2X'], p['m2X'], p['m2D']])
    XI = np.diag([p['xp'], p['xX'], p['xX'], p['xD']])
    return L, M2, XI

KEYS = list(P0.keys())

def betas(p):
    L, M2, XI = tensors(p)
    k = 1 / (16 * math.pi ** 2)
    LL = np.einsum('abef,efcd->abcd', L, L)
    bL = k * (LL + np.einsum('acef,efbd->abcd', L, L) + np.einsum('adef,efbc->abcd', L, L))
    bM = k * np.einsum('abef,ef->ab', L, M2)
    bX = k * np.einsum('abef,ef->ab', L, XI - np.eye(4) / 6)
    return dict(lp=bL[0, 0, 0, 0] / 6, g2=bL[0, 0, 1, 1] / 2, lX=bL[1, 1, 1, 1] / 6, h2=bL[0, 0, 3, 3] / 2,
                lD=bL[3, 3, 3, 3] / 6, al=bL[1, 1, 3, 3] / 2, m2p=bM[0, 0], m2X=bM[1, 1], m2D=bM[3, 3],
                xp=bX[0, 0], xX=bX[1, 1], xD=bX[3, 3]), bL

def run_to(mu):
    t1 = math.log(mu / mu0)
    f = lambda t, z: [betas(dict(zip(KEYS, z)))[0][k] for k in KEYS]
    if t1 == 0:
        return dict(P0)
    s = solve_ivp(f, [0, t1], [P0[k] for k in KEYS], rtol=1e-10, atol=1e-14)
    return dict(zip(KEYS, s.y[:, -1]))

_RUN_CACHE = {}
def couplings(H):
    key = round(H, 6)
    if key not in _RUN_CACHE:
        p = run_to(H)
        _RUN_CACHE[key] = (p['m2p'] + 12 * p['xp'] * H ** 2, p['lp'], p)
    return _RUN_CACHE[key]

def s_fun(H, tree=False):
    if tree:
        m2e, lp = -0.01 * v ** 2, 0.01
    else:
        m2e, lp, _ = couplings(H)
    a2 = m2e / H ** 2
    return (lambda y: (8 * math.pi ** 2 / 3) * (0.5 * a2 * y ** 2 + 0.25 * lp * y ** 4)), a2, lp

# ---------------- FP chain ----------------
_gl_x, _gl_w = np.polynomial.legendre.leggauss(12)

def build(H, M=None, hy=None, tree=False, fresh_node=True):
    s, a2, lp = s_fun(H, tree)
    ym = math.sqrt(-a2 / lp)                     # minimum
    smin = s(ym)
    # far end where s - smin = 90
    Y = brentq(lambda y: s(y) - smin - 90, ym, 50 * ym + 50)
    y_start = v / H
    if hy is None:
        hy = 0.004
    if fresh_node:
        k = max(4, int(round(y_start / hy)))
        hy = y_start / k
    Mn = int(math.ceil(Y / hy))
    y = hy * np.arange(1, Mn + 1)
    sv = s(y)
    # conductances between node i-1 and i, i = 1..Mn (node 0 is y=0 absorbing)
    yl = np.concatenate([[0.0], y[:-1]])
    mid = 0.5 * (yl + y); half = 0.5 * hy
    pts = mid[:, None] + half * _gl_x[None, :]
    sp = s(pts)
    ref = sp.max(axis=1)
    integ = np.log(half * (np.exp(sp - ref[:, None]) * _gl_w[None, :]).sum(axis=1)) + ref   # log ∫ e^{s}
    logc = math.log(Dc) - integ                                                              # log c_{i-1/2}
    return dict(y=y, s=sv, logc=logc, hy=hy, ym=ym, B=s(0) - smin, a2=a2, lp=lp, y_start=y_start)


def generator(G):
    """symmetrized generator S = E^{1/2} A E^{-1/2}, A P = dP/dn, E = diag(e^{s})."""
    y, sv, logc, hy = G['y'], G['s'], G['logc'], G['hy']
    Mn = len(y)
    c = np.exp(logc)                      # c[i] = conductance between node i-1 and i (i=0: absorbing face)
    # dP_i/dn = (1/h)[ c_i (rho_{i-1} - rho_i) + c_{i+1} (rho_{i+1} - rho_i) ], rho = e^{s} P
    # symmetrized: q = e^{s/2} P;  S_ii = -(c_i + c_{i+1}) e^{s_i}/h ; S_{i,i+1} = c_{i+1} e^{(s_i + s_{i+1})/2}/h
    cn = np.concatenate([c[1:], [0.0]])   # right face of each node, reflecting at end
    diag = -(c + cn) * np.exp(sv) / hy
    off = c[1:] * np.exp(0.5 * (sv[:-1] + sv[1:])) / hy
    S = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    return S


def solve_rates(G, n_grid):
    """returns dict with Lambda_1 (lowest eigenvalue), dp/dn on n_grid for relaxed and fresh starts."""
    S = generator(G)
    lam, U = np.linalg.eigh(S)
    lam = -lam                             # decay rates (positive)
    order = np.argsort(lam); lam = lam[order]; U = U[:, order]
    y, sv, hy = G['y'], G['s'], G['hy']
    c0 = math.exp(G['logc'][0])
    # flux out at y=0: J = c0 rho_1 = c0 e^{s_1} P_1 = c0 e^{s_1/2} q_1
    # relaxed start: P_i ∝ e^{-s_i}, normalized sum P h = 1  -> q = e^{s/2} P = e^{-s/2}/Z
    lz = logsumexp(-sv) + math.log(hy)
    q0 = np.exp(-0.5 * sv - lz)
    coef_rel = U.T @ q0
    # fresh start: delta at node y_start: P = 1/h at that node
    j = int(round(G['y_start'] / hy)) - 1
    q0f = np.zeros_like(sv); q0f[j] = math.exp(0.5 * sv[j]) / hy
    coef_fr = U.T @ q0f
    pref = c0 * math.exp(0.5 * sv[0]) * U[0, :]
    E = np.exp(-np.outer(n_grid, lam))
    J_rel = E @ (pref * coef_rel)
    J_fr = E @ (pref * coef_fr)
    # cumulative p_- = F/2 = (1 - survival)/2
    surv_rel = E @ ((U.T @ np.exp(-0.5 * sv)).T * coef_rel) * hy
    surv_fr = E @ ((U.T @ np.exp(-0.5 * sv)).T * coef_fr) * hy
    return dict(Lambda1=lam[0], Lambda2=lam[1], dpdn_rel=0.5 * J_rel, dpdn_fr=0.5 * J_fr,
                pminus_rel=0.5 * (1 - surv_rel), pminus_fr=0.5 * (1 - surv_fr), y_start_node=y[j])


def lambda1_mfpt(H, tree=False):
    """Lambda_1 ≈ 1/<tau>, <tau> = (1/D) ∫_0^∞ e^{s(u)} [∫_u^∞ e^{-s}]^2 du / ∫_0^∞ e^{-s}  (log space)."""
    s, a2, lp = s_fun(H, tree)
    ym = math.sqrt(-a2 / lp)
    Y = brentq(lambda y: s(y) - s(ym) - 120, ym, 50 * ym + 50)
    u = np.linspace(0, Y, 400001)
    su = s(u)
    du = u[1] - u[0]
    # tail integral I(u) = ∫_u^Y e^{-s}: cumulative from the right in log space
    w = -su
    logI = np.empty_like(u)
    # trapezoid in log space via cumulative logaddexp
    seg = np.logaddexp(w[:-1], w[1:]) + math.log(du / 2)
    logI[-1] = -np.inf
    logI[:-1] = np.logaddexp.accumulate(seg[::-1])[::-1]
    integrand = su + 2 * logI
    lnum = logsumexp(integrand[:-1]) + math.log(du) - math.log(Dc)
    lZ = logI[0]
    return math.exp(-(lnum - lZ)), s(0) - s(ym)


# ---------------- relic integral ----------------
def nprime_table():
    lT = np.linspace(math.log(1e-6), math.log(4e10), 4000)
    npv = np.array([C.nprime(math.exp(x)) for x in lT])
    Hs = np.array([C.H_of_T(math.exp(x)) for x in lT])
    ss = np.array([C.s_of_T(math.exp(x)) for x in lT])
    return lT, npv, Hs / ss

NPT = None

def relic_rho_over_s(H, Nstar, rates_fn, f=1.0, CE=None, Tcut=None, dnp=0.005):
    """rho/s = f C_E sigma ∫_{n'(Tcut)}^{N_0} dn' (dp/dn)(n'+DeltaN) H/s(T(n'))."""
    global NPT
    if NPT is None:
        NPT = nprime_table()
    lT, npv, Hs = NPT
    if CE is None:
        CE = 3 * 7.7866 / (4 * math.pi)
    if Tcut is None:
        Tcut = C.T_dom_iso()
    N0 = C.N0_instant(H)
    dN = Nstar - N0
    npc = C.nprime(Tcut)
    grid = np.arange(npc, N0, dnp)
    Hs_g = np.exp(np.interp(grid, npv, np.log(Hs)))
    n = grid + dN
    dp = rates_fn(n)
    val = np.trapezoid(dp * Hs_g, grid)
    return f * CE * C.sigma * val, dict(N0=N0, DeltaN=dN, npcut=npc)


def rates_for(H, M_h=0.004, tree=False):
    G = build(H, hy=M_h, tree=tree)
    ngrid = np.arange(0, 70.0001, 0.01)
    R = solve_rates(G, ngrid)
    rel = lambda n: np.interp(n, ngrid, R['dpdn_rel'])
    fr = lambda n: np.interp(n, ngrid, R['dpdn_fr'])
    return G, R, rel, fr, ngrid


if __name__ == '__main__':
    t0 = time.time()
    out = {}
    # running check
    p1 = run_to(1000.0)
    out['running_mu_1TeV'] = {k: p1[k] for k in ['lp', 'xp', 'xX', 'xD', 'm2p']}
    p3 = run_to(362.0)
    out['running_mu_362'] = {k: p3[k] for k in ['lp', 'xp', 'm2p']}
    out['beta_lp_check_mu0'] = dict(code=betas(P0)[0]['lp'],
                                    hand=(18 * 0.01 ** 2 + 4 * 0.01 ** 2 + 2 * 1e-8) / (16 * math.pi ** 2))
    # validation: p_- fresh N50/N60 at 250, 335, 1000, 2000 GeV and Lambda_1 at 1000, 501
    val = {}
    for H in [250.0, 335.0, 1000.0, 2000.0]:
        G = build(H, hy=0.004)
        R = solve_rates(G, np.array([50.0, 60.0]))
        val[f'{H:g}'] = dict(pm_fresh_N50=R['pminus_fr'][0], pm_fresh_N60=R['pminus_fr'][1], Lambda1=R['Lambda1'],
                             B=G['B'])
    for H in [501.0, 362.0]:
        G = build(H, hy=0.004)
        R = solve_rates(G, np.array([50.0]))
        val[f'{H:g}'] = dict(Lambda1=R['Lambda1'], B=G['B'])
    # grid convergence at 250 GeV fresh
    conv = {}
    for hy in [0.008, 0.004, 0.002]:
        G = build(250.0, hy=hy)
        R = solve_rates(G, np.array([50.0, 60.0]))
        conv[f'h={hy}'] = [R['pminus_fr'][0], R['pminus_fr'][1], R['Lambda1']]
    val['grid_conv_250'] = conv
    out['validation'] = val
    print('validation done', time.time() - t0, flush=True)

    # relic threshold: relaxed start, N* = 50, f = 1
    target = C.R_bound if hasattr(C, 'R_bound') else 0.12 * C.rho_c_over_h2 / (2 * math.pi ** 2 / 45 * 3.909 * C.T0 ** 3)
    out['relic_bound_used'] = target
    scan = []
    Hs_scan = [300, 320, 340, 350, 355, 360, 365, 370, 380, 400, 420, 440, 460, 480, 500, 520, 540]
    store = {}
    for H in Hs_scan:
        G, R, rel, fr, ng = rates_for(float(H))
        r50, meta = relic_rho_over_s(H, 50, rel)
        r60, _ = relic_rho_over_s(H, 60, rel)
        f50, _ = relic_rho_over_s(H, 50, fr)
        f60, _ = relic_rho_over_s(H, 60, fr)
        # closed form, quasi-stationary
        cf = 3 * 7.7866 / (4 * math.pi) * R['Lambda1'] / 2 * C.sigma * C.H_of_T(C.T_dom_iso()) / C.s_of_T(C.T_dom_iso())
        rec = dict(H=H, Lambda1=R['Lambda1'], B=G['B'], rel50=r50, rel60=r60, fr50=f50, fr60=f60, closed=cf,
                   DeltaN50=meta['DeltaN'],
                   dpdn_rel_at_n=[float(rel(np.array([x]))[0]) for x in [0, 4, 10, 14.3, 20, 30, 50]],
                   dpdn_fr_at_n=[float(fr(np.array([x]))[0]) for x in [4, 10, 14.3, 20, 30, 50]])
        scan.append(rec)
        print(H, rec['Lambda1'], r50, r60, f50, f60, cf, time.time() - t0, flush=True)
    out['scan'] = scan

    def thr(key):
        Hs = np.array([r['H'] for r in scan], float); val_ = np.array([r[key] for r in scan])
        lv = np.log(val_ / target)
        for i in range(len(Hs) - 1):
            if lv[i] < 0 <= lv[i + 1]:
                x = np.log(Hs[i]) + (0 - lv[i]) * (np.log(Hs[i + 1]) - np.log(Hs[i])) / (lv[i + 1] - lv[i])
                return float(np.exp(x))
        return None
    out['thresholds_f1'] = {k: thr(k) for k in ['rel50', 'rel60', 'fr50', 'fr60', 'closed']}
    # refine relaxed N50 by root-finding
    def g(H):
        G, R, rel, fr, ng = rates_for(H)
        return math.log(relic_rho_over_s(H, 50, rel)[0] / target)
    Hr = brentq(g, 340, 400, xtol=0.05)
    out['threshold_relaxed_N50_f1_root'] = Hr
    # grid-convergence of the threshold quantity at the root
    convr = {}
    for hy in [0.008, 0.004, 0.002]:
        G, R, rel, fr, ng = rates_for(Hr, M_h=hy)
        convr[f'h={hy}'] = relic_rho_over_s(Hr, 50, rel)[0] / target
    out['threshold_grid_conv'] = convr
    # tree-level (no running, no xi) variant
    def gt(H):
        G, R, rel, fr, ng = rates_for(H, tree=True)
        return math.log(relic_rho_over_s(H, 50, rel)[0] / target)
    out['threshold_relaxed_N50_f1_tree'] = brentq(gt, 330, 400, xtol=0.05)
    # f = 0.1 relaxed N50
    def g01(H):
        G, R, rel, fr, ng = rates_for(H)
        return math.log(0.1 * relic_rho_over_s(H, 50, rel)[0] / target)
    out['threshold_relaxed_N50_f0.1_root'] = brentq(g01, 400, 600, xtol=0.1)
    # Lambda1 for critic comparison: Lambda1/2 = 1e-11
    def gl(H):
        G = build(H, hy=0.004)
        return math.log(lambda1_mfpt(H)[0] / 2e-11)
    out['H_Lambda1_half_1e-11'] = brentq(gl, 200, 300, xtol=0.01)
    # Lambda_1 via MFPT vs eigen at 362 and 300
    out['mfpt_vs_eig'] = {f'{H}': dict(mfpt=lambda1_mfpt(H)[0], eig=solve_rates(build(H, hy=0.004), np.array([0.0]))['Lambda1'])
                          for H in [300.0, 362.0, 250.0]}
    print('relic done', time.time() - t0, flush=True)

    # regime A: N(all sizes) = ∫ dp/dn e^{3n'} dn' ≈ (Lambda1/2) e^{3 N0}/3 in today's Hubble volume, < 1e-3
    def NA(H, N0=None):
        L1, B = lambda1_mfpt(H)
        if N0 is None:
            N0 = C.N0_instant(H)
        return (L1 / 2) * math.exp(3 * N0) / 3
    out['regimeA_instant'] = brentq(lambda H: math.log(NA(H) / 1e-3), 120, 200, xtol=0.01)
    out['regimeA_exit_N60'] = brentq(lambda H: math.log(NA(H, 60.0) / 1e-3), 110, 200, xtol=0.01)
    out['regimeA_exit_N50'] = brentq(lambda H: math.log(NA(H, 50.0) / 1e-3), 110, 200, xtol=0.01)
    out['regimeA_kination_max'] = brentq(lambda H: math.log(NA(H, C.N0_instant(H) + math.log(3 * C.MP ** 2 * H ** 2 / C.rho_rad(4e-3)) / 12 + math.log(106.75 / 10.76) / 12) / 1e-3), 110, 200, xtol=0.01)
    L145, B145 = lambda1_mfpt(145.6)
    out['Lambda1_at_145.6'] = dict(Lambda1=L145, B=B145, N0=C.N0_instant(145.6), N=NA(145.6))
    out['regimeA_tree'] = brentq(lambda H: math.log((lambda1_mfpt(H, tree=True)[0] / 2) * math.exp(3 * C.N0_instant(H)) / 3 / 1e-3), 120, 200, xtol=0.01)
    out['runtime_s'] = time.time() - t0
    with open(os.path.join(OUT, 'results_fp.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({k: out[k] for k in out if k != 'scan'}, indent=1, default=float))
