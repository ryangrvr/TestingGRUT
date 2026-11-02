"""Inflationary stage of the mf_hybrid models [C given the model choices].

- calibrate(): tree clock coefficients (c1, m2) so that the full one-loop potential gives A_s and the target n_s at the
  pivot, N_back e-folds before sigma_c (slow roll, second order in the spectral index).
- ms(): exact background in e-folds plus Mukhanov-Sasaki modes for A_s, n_s, alpha_s, r (check of the slow-roll fit).
- the tuning diagnostics: the loop share of V' and V'' at the pivot and at sigma_c.
Usage: python3 -B inflation.py scan   (CW-flatness scan -> ../results/cw_scan.json)
"""
import math, json, sys, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, brentq
from scipy.integrate import quad
import mfh_common as mc

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
C2 = 4 * (math.log(2) + 0.5772156649) - 5      # Stewart-Lyth constant


def N_of_sigma(inf, s):
    """e-folds from sigma to sigma_c in slow roll: N = int V/(mP^2 V') d sigma (eps ~ 1e-8, exact to O(eps))."""
    p = inf.p
    f = lambda z: inf.V(z) / (p['mP'] ** 2 * inf.dV(z))
    lo = p['sigma_c'] * (1 + 1e-12)
    return quad(f, lo, s, limit=400, epsrel=1e-10)[0]


def sigma_of_N(inf, N):
    p = inf.p
    hi = p['sigma_c'] * 1.01 + 1e-3 * p['mP']
    while N_of_sigma(inf, hi) < N:
        hi = p['sigma_c'] + 2 * (hi - p['sigma_c'])
        if hi > 50 * p['mP']:
            raise ValueError('no pivot')
    return brentq(lambda z: N_of_sigma(inf, z) - N, p['sigma_c'] * (1 + 1e-12), hi, xtol=1e-14 * hi, rtol=1e-13)


class _Path:
    def __init__(self, inf, Nback, npts):
        Ns = np.linspace(0, Nback, npts)
        self.t = Ns
        self.y = np.array([[sigma_of_N(inf, n) if n > 0 else inf.p['sigma_c'] for n in Ns]])


def sr_path(inf, Nback, npts=400):
    return _Path(inf, Nback, npts)


def sr_obs(inf, s):
    p = inf.p
    mP = p['mP']
    V = inf.V(s)
    V1, V2, V3 = inf.derivs(s)
    eps = 0.5 * (mP * V1 / V) ** 2
    eta = mP ** 2 * V2 / V
    xi2 = mP ** 4 * V1 * V3 / V ** 2
    ns = 1 - 6 * eps + 2 * eta + (2.0 / 3) * eta ** 2 + (2.0 / 3) * (3 * C2 - 1) * xi2 + (4 * C2 - 14) * eps * eta
    As = V / (24 * math.pi ** 2 * eps * mP ** 4)
    alpha = 16 * eps * eta - 24 * eps ** 2 - 2 * xi2
    return dict(eps=eps, eta=eta, xi2=xi2, ns=ns, As=As, alpha_s=alpha, r=16 * eps, H=math.sqrt(V / 3) / mP)


def calibrate(p, ns_t, Nback, mu2=None, c_bias=0.0, guess=None, loop=True, As_t=None):
    """Solve (c1, m2) for A_s and n_s at the pivot Nback e-folds before sigma_c.

    Parametrised by the pivot value sigma_*: the total V' and V'' there are fixed by A_s and n_s (eta re-targeted for the
    second-order terms), and sigma_* is found by N(sigma_*) = Nback (1D root)."""
    mP, V0 = p['mP'], p['V0']
    eps_t = V0 / (24 * math.pi ** 2 * (As_t or mc.AS) * mP ** 4)
    eta_t = 0.5 * (ns_t - 1)

    def build(ss, eta_t):
        base = mc.Inflation(p, 0.0, 0.0, c_bias=c_bias, mu2=mu2, loop=loop)
        l1, l2, _ = base.derivs(ss)
        V = base.V(ss)
        V1 = math.sqrt(2 * eps_t) * V / mP
        V2 = eta_t * V / mP ** 2
        m2 = l2 - V2
        c1 = V1 - l1 + m2 * (ss - p['sigma_c'])
        return mc.Inflation(p, c1, m2, c_bias=c_bias, mu2=mu2, loop=loop)

    def nres(ss, eta_t):
        inf = build(ss, eta_t)
        xs = np.linspace(p['sigma_c'] * (1 + 1e-9), ss, 200)
        if min(inf.dV(z) for z in xs) <= 0:
            return -1e3
        return N_of_sigma(inf, ss) - Nback

    for it in range(4):
        lo = p['sigma_c'] * (1 + 1e-6)
        hi = p['sigma_c'] + 0.003 * mP
        while nres(hi, eta_t) < 0:
            hi = p['sigma_c'] + 1.6 * (hi - p['sigma_c'])
            if hi > p['sigma_c'] + 5 * mP:
                raise ValueError('no pivot')
        # bracket from below: smallest sigma_* with a monotonic path and N > 0 is near sigma_c
        ss = brentq(lambda z: nres(z, eta_t), p['sigma_c'] + 1e-6 * (hi - p['sigma_c']), hi, xtol=1e-13 * hi)
        inf = build(ss, eta_t)
        o = sr_obs(inf, ss)
        eta_t = eta_t - 0.5 * (o['ns'] - ns_t)
    inf = build(ss, eta_t)
    sol = sr_path(inf, Nback, npts=40)
    s_star = ss
    o = sr_obs(inf, s_star)
    o.update(c1=inf.c1, m2=inf.m2, sigma_star=s_star, sigma_star_MP=s_star / p['mP'],
             Delta_sigma_MP=(s_star - p['sigma_c']) / p['mP'], residual=[math.log(o['As'] / mc.AS), o['ns'] - ns_t],
             Nback=Nback)
    x = [0, 0]
    # check the slope is positive on the whole path (sigma rolls monotonically to sigma_c)
    ss = np.linspace(p['sigma_c'] * (1 + 1e-6), s_star + 0.05 * (s_star - p['sigma_c']), 400)
    dv = np.array([inf.dV(z) for z in ss])
    o.update(monotonic=bool(np.all(dv > 0)))
    pa = inf.parts(s_star)
    pc = inf.parts(p['sigma_c'] * 1.0001)
    o['tuning'] = dict(loop_over_tree_dV_star=pa['loop_dV'] / pa['tree_dV'], loop_over_tree_d2V_star=pa['loop_d2V'] / pa['tree_d2V'],
                       loop_d2V_star_H2=pa['loop_d2V'], tree_d2V_star_H2=pa['tree_d2V'],
                       loop_over_tree_dV_c=pc['loop_dV'] / pc['tree_dV'], g2_sigma_star=p['g'] ** 2 * s_star, X=p['X'])
    return inf, o, sol


def ms(inf, Nback_pivot, Nextra_before=12.0, dlnk=(-1.0, 0.0, 1.0)):
    """Exact background (fields in mP-free H units) + Mukhanov-Sasaki zeta_k; returns A_s, n_s, alpha_s at the pivot."""
    p = inf.p
    mP = p['mP']
    s0 = None
    for ne in (Nextra_before, 9.0, 7.5, 6.8):
        try:
            s0 = sigma_of_N(inf, Nback_pivot + ne)
            break
        except ValueError:
            continue
    # exact background in N (forward), start on the slow-roll attractor
    V0 = inf.V(s0)
    pi0 = -mP ** 2 * inf.dV(s0) / V0

    def bg(N, y):
        s, ps = y[0], y[1]
        eps = 0.5 * ps * ps / mP ** 2
        V = inf.V(s)
        H2 = V / (mP ** 2 * (3 - eps))
        return [ps, -(3 - eps) * ps - inf.dV(s) / H2]

    ev = lambda N, y: y[0] - p['sigma_c'] * (1 + 1e-7)
    ev.terminal = True
    b = solve_ivp(bg, (0, 200), [s0, pi0], rtol=1e-11, atol=[1e-8 * p['sigma_c'], 1e-14 * abs(pi0)], dense_output=True,
                  events=ev, method='DOP853')
    Nend = b.t_events[0][0]
    Npiv = Nend - Nback_pivot

    def HH(N):
        s, ps = b.sol(N)
        eps = 0.5 * ps * ps / mP ** 2
        return math.sqrt(inf.V(s) / (mP ** 2 * (3 - eps))), eps

    Hp, ep = HH(Npiv)
    kp = math.exp(Npiv) * Hp     # a = e^N
    out = []
    for dk in dlnk:
        k = kp * math.exp(dk)
        # start at k/aH = 200
        Ns = Npiv + dk - math.log(200.0)
        Ne = Npiv + dk + math.log(1e3)

        def mode(N, y):
            s, ps = b.sol(N)
            eps = 0.5 * ps * ps / mP ** 2
            H = math.sqrt(inf.V(s) / (mP ** 2 * (3 - eps)))
            # d eps/dN
            dps = -(3 - eps) * ps - inf.dV(s) * mP ** 2 * (3 - eps) / inf.V(s)
            etaH = 2 * dps / ps
            kaH = k / (math.exp(N) * H)
            zr, zi, dr, di = y
            return [dr, di, -(3 - eps + etaH) * dr - kaH ** 2 * zr, -(3 - eps + etaH) * di - kaH ** 2 * zi]

        s, ps = b.sol(Ns)
        eps = 0.5 * ps * ps / mP ** 2
        H = math.sqrt(inf.V(s) / (mP ** 2 * (3 - eps)))
        a = math.exp(Ns)
        z = a * math.sqrt(2 * eps) * mP
        # v = e^{-ik tau}/sqrt(2k), tau = -1/(aH) (leading); zeta = v/z, d zeta/dN = (dv/dtau)/(aH z) - zeta (1 + eta_H/2)
        tau = -1.0 / (a * H) * (1 + eps)
        v = np.exp(-1j * k * tau) / math.sqrt(2 * k)
        dv = -1j * k * v
        dps = -(3 - eps) * ps - inf.dV(s) * mP ** 2 * (3 - eps) / inf.V(s)
        etaH = 2 * dps / ps
        zeta = v / z
        dzeta = dv / (a * H * z) - zeta * (1 + 0.5 * etaH)
        m = solve_ivp(mode, (Ns, Ne), [zeta.real, zeta.imag, dzeta.real, dzeta.imag], rtol=1e-9, atol=1e-30,
                      method='DOP853')
        zf = m.y[0][-1] + 1j * m.y[1][-1]
        P = k ** 3 * abs(zf) ** 2 / (2 * math.pi ** 2)
        out.append((dk, P))
    lk = np.array([o[0] for o in out])
    lP = np.log([o[1] for o in out])
    c = np.polyfit(lk, lP, 2)
    rt = 2 * Hp ** 2 / (math.pi ** 2 * mP ** 2) / out[1][1]
    return dict(As=out[1][1], ns=1 + c[1], alpha_s=2 * c[0], r=rt, H_star=Hp, eps_star=ep, N_total_bg=Nend), b


def scan():
    """CW-flatness scan: for each (m_chi, m_sig) solve the clock and record the loop share and alpha_s."""
    rows = []
    H = 1.30e11
    for m_chi in (300.0, 1000.0, 4000.0):
        for m_sig in (3.0, 10.0, 30.0, 100.0, 300.0, 750.0, 1500.0, 3000.0):
            p = mc.derived(H, m_chi, m_sig)
            try:
                inf, o, _ = calibrate(p, mc.NS_PLANCK, 45.0, mu2=None)
                ok = abs(o['residual'][0]) < 1e-6 and abs(o['residual'][1]) < 1e-6
            except Exception as e:  # noqa
                rows.append(dict(m_chi=m_chi, m_sig=m_sig, error=str(e)))
                continue
            rows.append(dict(m_chi=m_chi, m_sig=m_sig, lam=p['lam'], g=p['g'], sigma_c_MP=p['sigma_c_MP'],
                             v_MP=p['v_MP'], solved=ok, monotonic=o['monotonic'], ns=o['ns'], alpha_s=o['alpha_s'],
                             sigma_star_MP=o['sigma_star_MP'], **o['tuning']))
            r = rows[-1]
            print('m_chi %6.0f m_sig %6.0f X %.3g g2s* %.3g solved %s mono %s alpha %.2e loop/tree dV %.3g d2V %.3g' % (
                m_chi, m_sig, p['X'], r['g2_sigma_star'], ok, o['monotonic'], o['alpha_s'],
                r['loop_over_tree_dV_star'], r['loop_over_tree_d2V_star']))
    json.dump(dict(label='C (given the model choices)', note='full chi one-loop (MS-bar, mu = M_chi(2 sigma_c)); clock c1, m2 '
                   'solved for A_s and n_s = 0.9649 at 45 e-folds before sigma_c; m_chi and m_sig are the vacuum masses in H units',
                   rows=rows), open(os.path.join(RES, 'cw_scan.json'), 'w'), indent=1)


if __name__ == '__main__':
    if 'scan' in sys.argv:
        scan()
