"""sp_point track (Claude, own code): model, inflationary minimum, run wrapper for sprel.c.

Physical potential in H_i = 1 units (complex fields, kinetic |dz|^2):
  V = lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(k Y* X^3) + 2Re(mu S* Y^2) - cH|X|^2 (inflation only)
Couplings fixed at the stated point (unchanged): M_Y=1e12, h=1e11, mu=5e11 GeV, k=1e-3, lx=1e-3, ly=ls=1.
"""
import os
import subprocess
import numpy as np
from scipy.optimize import fsolve, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/sp_point'
BIN = os.path.join(SCR, 'bin', 'sprel')
LX, LY, LS, KAP = 1e-3, 1.0, 1.0, 1e-3
M_GEV, H_GEV, MU_GEV = 1e12, 1e11, 5e11
MRED, GSTAR, YQ = 2.435e18, 106.75, 4.8294009076988175e-9
H_STAR = 1.1211167304108865e13
QN = ["NX", "NY", "NS", "Neff", "EX", "epsX", "NcircX", "JY", "JYeig", "EYexc", "EYtot", "Gauss", "deps", "dNX",
      "lnT", "absY"]


def build():
    src = os.path.join(HERE, 'sprel.c')
    os.makedirs(os.path.dirname(BIN), exist_ok=True)
    if not os.path.exists(BIN) or os.path.getmtime(BIN) < os.path.getmtime(src):
        subprocess.run(['gcc', '-O3', '-o', BIN, src, '-lm'], check=True)


SCALE = dict(M=1.0, h=1.0, mu=1.0)  # diagnostic multipliers only (mechanism tests); 1 = stated couplings


def pars(H):
    return SCALE['M'] * M_GEV / H, SCALE['h'] * H_GEV / H, SCALE['mu'] * MU_GEV / H


def V(z, H, cH=0.0):
    M, h, mu = pars(H)
    x, y, s = z
    return float(LX * abs(x)**4 + LY * abs(y)**4 + LS * (abs(s)**2 - h * h / 2)**2 + M * M * abs(y)**2
                 + 2 * np.real(KAP * np.conj(y) * x**3 + mu * np.conj(s) * y * y) - cH * abs(x)**2)


def radial_grad(r, H, cH, cos_mu=-1.0):
    """Gradient of V on X=x, Y=-y, S=-s e^{i delta} (kappa phase aligned); cos_mu = cos(theta_mu) = -cos(delta)."""
    M, h, mu = pars(H)
    x, y, s = r
    return np.array([4 * LX * x**3 - 6 * KAP * y * x * x - 2 * cH * x,
                     4 * LY * y**3 + 2 * M * M * y - 2 * KAP * x**3 + 4 * mu * s * y * cos_mu,
                     4 * LS * s * (s * s - h * h / 2) + 2 * mu * y * y * cos_mu])


def radial_V(r, H, cH, cos_mu=-1.0):
    M, h, mu = pars(H)
    x, y, s = r
    return (LX * x**4 + LY * y**4 + LS * (s * s - h * h / 2)**2 + M * M * y * y - 2 * KAP * y * x**3
            + 2 * mu * s * y * y * cos_mu - cH * x * x)


def minimum(H, cH=1.0, cos_mu=-1.0):
    """Global radial minimum on X=x, Y=-y, S=-s e^{i delta} (kappa phase aligned), cos_mu = -cos(delta).

    Own method: eliminate x(y) and s(y) by their stationarity conditions (unique positive roots; s may sit at the
    s=0 boundary when cos_mu>0), scan y on a log grid for every sign change of the remaining condition, refine each
    root with brentq, and keep the root of lowest V. Independent multi-start 6D complex minimisation is a check
    (min_check below).
    """
    from scipy.optimize import brentq
    M, h, mu = pars(H)
    def x_of(y):
        return (3 * KAP * y + np.sqrt(9 * KAP**2 * y * y + 8 * LX * cH)) / (4 * LX)
    def s_of(y):
        b = -mu * y * y * cos_mu / (2 * LS)   # s(s^2 - h^2/2) = b
        f = lambda s: s * (s * s - h * h / 2) - b
        if b >= 0:
            lo, hi = h / np.sqrt(2), h / np.sqrt(2) + abs(b)**(1 / 3) + 1e-14
            return brentq(f, lo, hi, xtol=1e-16, rtol=1e-15) if b > 0 else h / np.sqrt(2)
        # b<0: candidate roots in (0, h/sqrt2); minimum of the conditional quartic over s>=0
        cands = [0.0]
        grid = np.linspace(0, h / np.sqrt(2), 2001)
        fv = f(grid)
        for k in range(len(grid) - 1):
            if fv[k] * fv[k + 1] < 0:
                cands.append(brentq(f, grid[k], grid[k + 1], xtol=1e-16, rtol=1e-15))
        Vs = [LS * (c * c - h * h / 2)**2 + 2 * mu * y * y * cos_mu * c for c in cands]
        return cands[int(np.argmin(Vs))]
    def g(y):
        return 4 * LY * y**3 + 2 * M * M * y - 2 * KAP * x_of(y)**3 + 4 * mu * s_of(y) * y * cos_mu
    ys = np.geomspace(1e-8, 1e3, 4000)
    gv = np.array([g(y) for y in ys])
    roots = []
    for k in range(len(ys) - 1):
        if gv[k] * gv[k + 1] < 0:
            roots.append(brentq(g, ys[k], ys[k + 1], xtol=1e-15, rtol=1e-15))
    best = None
    for y in roots:
        r = np.array([x_of(y), y, s_of(y)])
        v = radial_V(r, H, cH, cos_mu)
        if best is None or v < best[0]:
            best = (v, r)
    return best[1]


def min_check(H, cH=1.0, n=24, seed=1):
    """Multi-start unconstrained 6D minimisation of the full complex inflationary potential (check only)."""
    rng = np.random.default_rng(seed)
    ref = minimum(H, cH)
    vref = radial_V(ref, H, cH)
    best = np.inf
    for _ in range(n):
        u0 = rng.normal(size=6) * np.repeat(np.abs(ref) * 2 + 0.1, 2)
        f = lambda u: V(u[0::2] + 1j * u[1::2], H, cH)
        r = minimize(f, u0, method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
        best = min(best, r.fun)
    return dict(V_reference=float(vref), V_best_multistart=float(best), rel_diff=float((best - vref) / abs(vref)))


def masses(H, cH=1.0):
    """Radial and angular mass^2 (canonical, in H^2) at the aligned inflationary minimum."""
    x, y, s = minimum(H, cH)
    M, h, mu = pars(H)
    Hs = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            e = 1e-5
            di = np.zeros(3); dj = np.zeros(3); di[i] = e; dj[j] = e
            r = np.array([x, y, s])
            Hs[i, j] = (radial_V(r + di + dj, H, cH) - radial_V(r + di - dj, H, cH) - radial_V(r - di + dj, H, cH)
                        + radial_V(r - di - dj, H, cH)) / (4 * e * e)
    rad = np.linalg.eigvalsh(Hs / 2)  # canonical radial field sqrt2 r
    K = np.diag([x * x, y * y, s * s])
    v1 = np.array([3., -1., 0.]); v2 = np.array([0., 2., -1.])
    W = KAP * y * x**3 * np.outer(v1, v1) + mu * s * y * y * np.outer(v2, v2)
    ev, evec = np.linalg.eig(np.linalg.solve(K, W))
    order = np.argsort(ev.real)
    return dict(amplitudes=[x, y, s], radial_mass2=rad.tolist(), angular_mass2=ev.real[order].tolist(),
                angular_vectors=evec.real[:, order].T.tolist())


def initial_z(H, delta, cH=1.0, conditional=False):
    if conditional:
        x, y, s = minimum(H, cH, cos_mu=-np.cos(delta))
    else:
        x, y, s = minimum(H, cH)
    return np.array([x, -y, -s * np.exp(1j * delta)], complex)


def yield_factor(H):
    T = (90 * MRED**2 * H * H / (np.pi**2 * GSTAR))**0.25
    return H**3 / (2 * np.pi**2 / 45 * GSTAR * T**3)


def run(H, z0, a_end=1000., method=8, cpar=0.05, windows=(1, 100, 300, 600, 1000), ns=200, tangent=None,
        v0=None):
    """tangent: 12-vector (dchi_re/im x3, dpi_re/im x3) or None."""
    build()
    M, h, mu = pars(H)
    args = [BIN, str(method), repr(float(cpar)), repr(M), repr(h), repr(mu), repr(float(a_end)), str(ns),
            '1' if tangent is not None else '0', str(len(windows))] + [repr(float(w)) for w in windows]
    for zz in z0:
        args += [repr(float(zz.real)), repr(float(zz.imag))]
    if v0 is not None or tangent is not None:
        vv = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
        for zz in vv:
            args += [repr(float(zz.real)), repr(float(zz.imag))]
    if tangent is not None:
        args += [repr(float(t)) for t in tangent]
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(p.stderr)
    hist, win, final, end = [], {}, None, None
    for line in p.stdout.splitlines():
        f = line.split()
        v = list(map(float, f[1:]))
        if f[0] == 'S':
            hist.append(dict(a=v[0], **dict(zip(QN, v[1:]))))
        elif f[0] == 'W':
            d = dict(zip(QN, v[3:3 + len(QN)]))
            sqs = v[3 + len(QN):]
            d['NX_std'] = float(np.sqrt(max(0, sqs[0] - d['NX']**2)))
            d['Neff_std'] = float(np.sqrt(max(0, sqs[1] - d['Neff']**2)))
            d['epsX_std'] = float(np.sqrt(max(0, sqs[2] - d['epsX']**2)))
            d['YX_over_YQ'] = d['NX'] * yield_factor(H) / YQ
            d['eps_of_means'] = d['NX'] / d['NcircX']
            d['pairs_per_charge'] = 1.5 * d['JY'] / abs(d['NX']) if d['NX'] else None
            win['%g-%g' % (v[0], v[1])] = d
        elif f[0] == 'F':
            final = v
        elif f[0] == 'END':
            end = dict(steps=int(v[0]), rejects=int(v[1]), max_gauss_rel=v[2], tangent_rescaled=bool(v[3]) if len(v) > 3 else False)
    return dict(H=H, windows=win, history=hist, final=final, controls=end)


def delta_tangent(H, delta, cH=1.0):
    """Tangent of the initial state w.r.t. the supplied source phase (aligned radii held)."""
    x, y, s = minimum(H, cH)
    d = -1j * s * np.exp(1j * delta)
    t = np.zeros(12)
    t[4], t[5] = d.real, d.imag      # dchi_S
    t[10], t[11] = d.real, d.imag    # dpi_S (zero velocity: pi = chi)
    return t
