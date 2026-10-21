#!/usr/bin/env python3
"""v_heating check 5: de Sitter stochastic spectators (phi, q1, q2, D) of the unchanged potential.

Independent verifier code.  Starobinsky-Yokoyama: per e-fold, d z_a = -dU/dz_a/(3H^2) dN + (H/2pi) dW_a for each
canonical real field; P_eq ~ exp(-8 pi^2 U/(3 H^4)).  Curvature terms (1/2) xi_a R z_a^2 with R = 12 H^2 and xi run
at one loop from xi(mu0 = 100 GeV) = 0 to mu = H (own derivation of the beta matrix; couplings frozen otherwise).

Parts:
 A. marginal s(phi) = a U_phi - ln I_X(phi) - ln I_D(phi) (closed forms), local/global restoration thresholds;
 B. odd eigenvalue Lambda_1 of the 1D FP on the marginal and on the X = D = 0 slice; validation on pure quartic;
 C. domain-count thresholds for both definitions (hidden_heating N_dom(Delta), inflation_spectator N_reg);
 D. Langevin: <r^2>(N) of X at phi = 0 from X = 0 (finite-duration restoration);
 E. Langevin: fresh start phi = +v, X = D = 0, p_- after N = 50, 60 e-folds with X, D fluctuating vs X, D frozen.
"""
import json, os, sys, time
import numpy as np
from scipy import special, integrate, linalg, optimize

OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
v = 1000.0
lp, g2, lX, h2, lD, mD0 = 0.01, 0.01, 0.001, 1e-4, 0.001, 100.0
mu0 = 100.0
pc = 0.3116077


def xi_at(H):
    """one-loop xi running, 16pi^2 dxi_a/dt = sum_b Bm_ab (xi_b - 1/6); couplings frozen at mu0 (alpha = 0)."""
    if H <= mu0:
        return np.zeros(3)
    Bm = np.array([[6 * lp, 4 * g2, 2 * h2],
                   [2 * g2, 8 * lX, 0.0],
                   [2 * h2, 0.0, 6 * lD]]) / (16 * np.pi**2)
    t = np.log(H / mu0)
    return (np.eye(3) - linalg.expm(Bm * t)) @ np.full(3, 1 / 6)


# ---------------- A. marginal -----------------
def s_marginal(phi, H, withX=True, withD=True, xi=None):
    if xi is None:
        xi = xi_at(H)
    a = 8 * np.pi**2 / (3 * H**4)
    s = a * (lp / 4 * (phi**2 - v**2)**2 + 6 * xi[0] * H**2 * phi**2)
    if withX:
        A = g2 * phi**2 / 2 + 6 * xi[1] * H**2          # coefficient of r^2
        x = A * np.sqrt(a / lX)
        s = s - np.log(special.erfcx(x))                # + phi-independent constant dropped
    if withD:
        M2 = mD0**2 + h2 * phi**2 + 12 * xi[2] * H**2
        al = a * M2 / 2; be = a * lD / 4
        z = al**2 / (8 * be)
        ID = 0.5 * np.sqrt(al / (2 * be)) * special.kve(0.25, z)
        s = s - np.log(ID)
    return s


def restoration_thresholds():
    def d2s0(H):
        e = 1e-3 * v
        return (s_marginal(e, H) - s_marginal(0.0, H))       # sign of curvature at 0 (s even)
    def glob(H):
        phis = np.linspace(0, 2.5 * v, 50001)
        s = s_marginal(phis, H)
        return s[0] - np.min(s[1:])
    Hloc = optimize.brentq(d2s0, 300, 450, xtol=1e-6)
    Hglob = optimize.brentq(glob, 380, 500, xtol=1e-6)
    # leading estimate
    rr = np.sqrt(3 / (2 * np.pi**3 * lX))
    lead = np.sqrt(lp / (g2 * rr))
    # same without xi and without D
    Hloc_noxi = optimize.brentq(lambda H: s_marginal(1.0, H, xi=np.zeros(3)) - s_marginal(0.0, H, xi=np.zeros(3)), 300, 450)
    Hglob_noxi = optimize.brentq(lambda H: (lambda s: s[0] - np.min(s[1:]))(s_marginal(np.linspace(0, 2.5 * v, 50001), H, xi=np.zeros(3))), 380, 500)
    return dict(H_local_GeV=Hloc, H_global_GeV=Hglob, leading_estimate_H_over_v=lead, r2_eq_over_H2=rr,
                H_local_noxi=Hloc_noxi, H_global_noxi=Hglob_noxi)


# ---------------- B. Lambda_1 -----------------
D_ = 1 / (8 * np.pi**2)


def lambda1(sfun, ymax, n=6000):
    """odd eigenvalue: Dirichlet at y = 0 on [0, ymax]; s(y) dimensionless potential in y = phi/H.
    returns (eig-based, mfpt-based)."""
    h = ymax / n
    y = (np.arange(n) + 0.5) * h
    s = sfun(y)
    k = D_ / h**2
    up = np.exp((s[:-1] - s[1:]) / 2)      # rate factor i -> i+1
    dn = np.exp((s[1:] - s[:-1]) / 2)      # rate factor i+1 -> i
    diag = np.zeros(n)
    diag[:-1] -= k * up
    diag[1:] -= k * dn
    diag[0] -= 2 * k                        # killing into the mirror cell (odd sector)
    off = np.full(n - 1, k)
    w = linalg.eigh_tridiagonal(-diag, -off, select='i', select_range=(0, 0), eigvals_only=True)
    lam_eig = float(w[0])
    # MFPT from the well minimum to y = 0 (continuous formula, log space)
    yy = np.linspace(0, ymax, 200001)
    ss = sfun(yy)
    m = np.argmin(ss)
    # inner(y) = int_y^ymax e^{-s(z)} dz
    smin = ss.min()
    ez = np.exp(-(ss - smin))
    inner = integrate.cumulative_trapezoid(ez[::-1], yy[::-1], initial=0)[::-1] * -1
    inner = np.maximum(inner, 1e-300)
    logint = ss[:m + 1] - smin + np.log(inner[:m + 1])     # log of e^{s} * inner, shifted
    mx = logint.max()
    val = np.trapezoid(np.exp(logint - mx), yy[:m + 1]) if m > 0 else 0.0
    if val <= 0:
        return lam_eig, float('nan'), float('nan'), 0.0
    logtau = mx + np.log(val) - np.log(D_)
    lam_mfpt = float(np.exp(-logtau)) if logtau < 700 else 0.0
    return lam_eig, lam_mfpt, float(logtau), float(ss[0] - smin)


def lambda1_marg(H, withX=True, withD=True):
    ymax = max(3.0 * v / H, 7.0)
    f = lambda y: s_marginal(y * H, H, withX, withD)
    return lambda1(f, ymax)


def best(le, lm, logtau):
    # eigenvalue reliable above ~1e-9; else 1/tau (natural log returned for extreme values)
    if le > 1e-8 or not np.isfinite(logtau):
        return np.log(le)
    return -logtau


def validation():
    out = {}
    lam = 0.01
    f = lambda y: 8 * np.pi**2 / 3 * lam / 4 * y**4
    le, lm, lt, S = lambda1(f, 8.0)
    out['pure_quartic_L1_over_sqrt_lam'] = le / np.sqrt(lam)
    out['SY_value'] = 1.36859 * np.sqrt(1 / (24 * np.pi**2))
    # double well without X, D at H = 250 (eigen vs MFPT vs Kramers SY eq. 58 form)
    for H in (250.0, 300.0):
        ff = lambda y: s_marginal(y * H, H, False, False, xi=np.zeros(3))
        le, lm, lt, S = lambda1(ff, max(3 * v / H, 7.0))
        kram = np.sqrt(2) * lp * v**2 / (3 * np.pi * H**2) * np.exp(-2 * np.pi**2 * lp * v**4 / (3 * H**4))
        out[f'dw_H{int(H)}'] = dict(eig=le, mfpt=lm, kramers_SY58=kram, S=S)
    return out


def scan_lambda1():
    rows = []
    for Hv in [0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40,
               0.42, 0.44, 0.46, 0.50]:
        H = Hv * v
        le, lm, lt, S = lambda1_marg(H)
        le1, lm1, lt1, S1 = lambda1_marg(H, False, False)
        rows.append(dict(H_over_v=Hv, S_marg=S, L1_marg_eig=le, L1_marg_mfpt=lm, ln_L1_marg=best(le, lm, lt),
                         S_1d=S1, L1_1d_eig=le1, L1_1d_mfpt=lm1, ln_L1_1d=best(le1, lm1, lt1)))
    return rows


def thresholds(rows):
    Hs = np.array([r['H_over_v'] for r in rows]) * v
    out = {}
    for key in ('ln_L1_marg', 'ln_L1_1d'):
        lnL = np.array([r[key] for r in rows])
        def H_for(target_ln):
            # monotone increasing in H in barrier regime: interpolate in H vs lnL
            idx = np.where(lnL >= target_ln)[0]
            if idx.size == 0 or idx[0] == 0:
                return None
            i = idx[0]
            return float(np.interp(target_ln, lnL[i - 1:i + 1], Hs[i - 1:i + 1]))
        d = {}
        # hidden_heating definition: N_dom(Delta) = (L1/2)(e^{3 Delta} - 1)/3 = 1
        for Dl in (5, 10, 20):
            d[f'Ndom1_Delta{Dl}'] = H_for(np.log(6.0) - np.log(np.expm1(3 * Dl)))
        # inflation_spectator definitions (relaxed-start count): N_reg = (L1/2)(e^{3n}-1)/3 < 1e-3
        for nn, lab in ((50, 'Nreg_all_N50'), (60, 'Nreg_all_N60'), (17.77, 'Nreg_after_1MeV')):
            d[lab + '_lt_1e-3'] = H_for(np.log(6e-3) - np.log(np.expm1(3 * nn)))
        d['L1_N60_eq_0.1'] = H_for(np.log(0.1 / 60))
        d['L1_N60_eq_1'] = H_for(np.log(1 / 60))
        out[key] = d
    return out


# ---------------- D/E. Langevin -----------------
def langevin(H, N_list, nsamp, fluct=True, phi0=None, seed=1, dN=0.01, track_r2=False):
    rng = np.random.default_rng(seed)
    xi = xi_at(H)
    phi = np.full(nsamp, v if phi0 is None else phi0)
    q1 = np.zeros(nsamp); q2 = np.zeros(nsamp); Dd = np.zeros(nsamp)
    sq = H / (2 * np.pi) * np.sqrt(dN)
    c = dN / (3 * H**2)
    out = {}
    Nmax = max(N_list)
    nsteps = int(round(Nmax / dN))
    marks = {int(round(N / dN)): N for N in N_list}
    r2hist = []
    for it in range(1, nsteps + 1):
        r2 = q1 * q1 + q2 * q2
        dphi = lp * (phi * phi - v * v) * phi + 12 * xi[0] * H**2 * phi
        if fluct:
            dphi += g2 * phi * r2 + h2 * phi * Dd * Dd
            mq = g2 * phi * phi + lX * r2 + 12 * xi[1] * H**2
            dq1 = mq * q1; dq2 = mq * q2
            dD = (mD0**2 + h2 * phi * phi + 12 * xi[2] * H**2) * Dd + lD * Dd**3
            q1 = q1 - c * dq1 + sq * rng.standard_normal(nsamp)
            q2 = q2 - c * dq2 + sq * rng.standard_normal(nsamp)
            Dd = Dd - c * dD + sq * rng.standard_normal(nsamp)
        phi = phi - c * dphi + sq * rng.standard_normal(nsamp)
        if it in marks:
            out[marks[it]] = float(np.mean(phi < 0))
        if track_r2 and it % 100 == 0:
            r2hist.append((it * dN, float(np.mean(q1 * q1 + q2 * q2)) / H**2))
    return out, r2hist


def x_growth():
    """<r^2>/H^2 vs N for X at phi = 0 (only X fluctuating, quartic lX), from X = 0."""
    rng = np.random.default_rng(7)
    n = 200000; dN = 0.01; H = 1.0
    q1 = np.zeros(n); q2 = np.zeros(n)
    res = {}
    lamq = lX
    for it in range(1, 30001):
        r2 = q1 * q1 + q2 * q2
        q1 = q1 - dN / 3 * lamq * r2 * q1 + np.sqrt(dN) / (2 * np.pi) * rng.standard_normal(n)
        q2 = q2 - dN / 3 * lamq * r2 * q2 + np.sqrt(dN) / (2 * np.pi) * rng.standard_normal(n)
        if it in (3000, 6000, 10000, 20000, 30000):
            res[it * dN] = float(np.mean(q1 * q1 + q2 * q2))
    return res


def main():
    t0 = time.time()
    part = sys.argv[1] if len(sys.argv) > 1 else 'ABC'
    res = {}
    fn = os.path.join(OUT, 'vh4_desitter.json')
    if os.path.exists(fn):
        res = json.load(open(fn))
    if 'A' in part:
        res['xi_at_H'] = {str(H): xi_at(H).tolist() for H in (379.0, 436.0, 1e3, 1e13)}
        res['restoration'] = restoration_thresholds()
    if 'B' in part:
        res['validation'] = validation()
        res['scan'] = scan_lambda1()
    if 'C' in part:
        res['thresholds'] = thresholds(res['scan'])
    if 'D' in part:
        g = x_growth()
        res['X_growth_r2_over_H2'] = g
        # finite-duration restoration H/v from g2 <r^2> = lp v^2 (xi neglected), free vs quartic growth
        res['finite_N_restoration_H_over_v'] = {
            str(N): dict(free=float(np.pi * np.sqrt(2 / N)), quartic=float(np.sqrt(lp / (g2 * g[N]))))
            for N in (30.0, 60.0, 100.0)}
    if 'E' in part:
        rows = []
        for H in [400.0, 500.0, 700.0, 1000.0, 1400.0]:
            pf, _ = langevin(H, [50, 60], 60000, fluct=True, seed=11)
            p0, _ = langevin(H, [50, 60], 60000, fluct=False, seed=12)
            rows.append(dict(H=H, pminus_fluct_N50=pf[50], pminus_fluct_N60=pf[60],
                             pminus_frozen_N50=p0[50], pminus_frozen_N60=p0[60]))
            print(rows[-1], time.time() - t0, flush=True)
        res['langevin_fresh_start'] = rows
    if 'F' in part:
        rows = res.get('langevin_fresh_start', [])
        for H in [2000.0, 2500.0]:
            pf, _ = langevin(H, [50, 60], 60000, fluct=True, seed=13)
            p0, _ = langevin(H, [50, 60], 60000, fluct=False, seed=14)
            rows.append(dict(H=H, pminus_fluct_N50=pf[50], pminus_fluct_N60=pf[60],
                             pminus_frozen_N50=p0[50], pminus_frozen_N60=p0[60]))
            print(rows[-1], time.time() - t0, flush=True)
        res['langevin_fresh_start'] = sorted(rows, key=lambda r: r['H'])
    res.setdefault('runtime_s', {})[part] = time.time() - t0
    os.makedirs(OUT, exist_ok=True)
    with open(fn, 'w') as f:
        json.dump(res, f, indent=1, default=float)
    print(json.dumps({k: res[k] for k in res if k not in ('scan',)}, indent=1, default=float))


if __name__ == '__main__':
    main()
