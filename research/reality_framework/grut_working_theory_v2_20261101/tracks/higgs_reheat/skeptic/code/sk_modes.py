"""Skeptic of track higgs_reheat: independent re-computation of the post-inflation long-wave Higgs variance and the
danger exponent D, in a different formulation (conformal time, rescaled field u = a h, own background, scipy Hankel
Bunch-Davies start) and with an explicit choice of field normalisation.  Own code (Claude, skeptic pass, 2 Nov 2026).

Only the SM input (2-loop lambda_eff(h) and the Higgs anomalous dimension Gamma(mu)) is taken from the track's smpot.py
(cf_vacuum's pipeline; the brief says not to redo it).  Everything else is written here.

Field normalisation (the point under test).  smpot's lambda_eff_flat(h) includes e^{4 Gamma(h)}: V = lambda_eff h0^4/4
with h0 the field renormalised at M_t (Buttazzo et al. convention).  The canonically normalised fluctuation field at the
scale of the problem (mu ~ H_i) is phi = e^{Gamma(H_i)} h0; its mass in de Sitter is xi R (tree), the mass the mode
equations use.  Hence the barrier for the mode variance sits at phi_max = e^{Gamma} h_max and the tachyonic curvature
is -d^2V/dphi^2 = -e^{-2 Gamma} d^2V/dh0^2.  The track's modes.py compares the canonical variance with h_max in h0 units
and uses d^2V/dh0^2 ('track' normalisation below); 'canon' is the consistent one.

Equations (H_i = 1, M_P = 1 for the background, a(0) = 1, eta = conformal time with eta(0) = 0):
    a' = a^2 H, I' = a J, J' = a[-(3H + Gam) J - b^2 I], r' = a[Gam J^2 - 4 H r], tau' = a
    u_k'' = -[k^2 + a^2 ((xi - 1/6) R + m_p^2)] u_k,      h_k = u_k / a
    R = 2 b^2 I^2 - J^2,  H^2 = (J^2/2 + b^2 I^2/2 + r)/3
BD start at tau = 0 (de Sitter, eta_dS = -1): u = (sqrt(pi)/2) H_nu(k), u' = (sqrt(pi)/2)(-H_nu(k)/2 - k H_nu'(k)),
nu = sqrt(9/4 - 12 xi) (real here).
Danger (same criterion as the track, re-implemented): for a smoothed field value fbar >= barrier, the first time tau_c
with m_c^2 = -V''(fbar) - xi <R> > (cH H)^2; k_c = a m_c; sigma_c^2 = int_{k<k_c} dlnk k^3/(2 pi^2) |h_k|^2;
D = log10[e^{3N} k_c^3 e^{-x}(1+x)], x = fbar^2/(2 sigma_c^2); D_max over fbar.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import hankel1

HERE = os.path.dirname(os.path.abspath(__file__))
TRACK_CODE = os.path.abspath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, TRACK_CODE)
import smpot as S  # noqa: E402  (SM input only)

LN10 = np.log(10.0)


class SMIn:
    """lambda_eff(h0) (incl. e^{4Gamma}), Gamma(mu), barrier h_max; V and V'' in h0 and in canonical phi units."""

    def __init__(self, Mt):
        self.sm = S.SM(Mt)
        self.Mt = Mt
        self.flat = self.sm.flat_scales()
        lh = np.linspace(np.log(1e9), np.log(1e15), 1500)
        le = self.sm.lam_eff_flat(np.exp(lh))
        self.V = CubicSpline(lh, le * np.exp(4 * lh) / 4)     # V(h0) in GeV^4 as a function of ln h0

    def Gam(self, mu):
        return float(self.sm.cpl(mu)['Gam'])

    def d2V_dh0(self, h0):
        lh = np.log(h0)
        V1 = self.V(lh, 1)
        V2 = self.V(lh, 2)
        return (V2 - V1) / h0 ** 2                           # d^2V/dh^2 = (V_ll - V_l)/h^2

    def xi_H(self, H_GeV, xi_EW):
        return float(self.sm.xi(H_GeV, xi_EW))

    def xi_EW(self, H_GeV, xi_H):
        c, c0 = self.sm.cpl(H_GeV), self.sm.cpl(self.Mt)
        fac = float(np.exp(c['Ixi'] - c0['Ixi']))
        return 1 / 6 + (xi_H - 1 / 6) / fac


def bd(ks, m2):
    if m2 >= 9 / 4 - 1e-9:                                   # heavy: complex order, mpmath fallback
        import mpmath as mp
        nu = mp.sqrt(mp.mpf(9) / 4 - mp.mpf(m2))
        pref = mp.sqrt(mp.pi) / 2 * mp.e ** (-mp.pi * mp.im(nu) / 2)
        u, up = [], []
        for k in ks:
            k = mp.mpf(float(k))
            Hn = mp.hankel1(nu, k)
            Hd = mp.hankel1(nu - 1, k) - nu / k * Hn
            u.append(complex(pref * Hn))
            up.append(complex(pref * (-Hn / 2 - k * Hd)))
        return np.array(u), np.array(up)
    nu = np.sqrt(9 / 4 - m2)
    Hn = hankel1(nu, ks)
    Hd = 0.5 * (hankel1(nu - 1, ks) - hankel1(nu + 1, ks))
    c = np.sqrt(np.pi) / 2
    return c * Hn, c * (-0.5 * Hn - ks * Hd)


def evolve(xi, b=3.0, gam=0.5, ks=None, tau_end=45.0, n_eta=9000, rtol=1e-9, atol=1e-12, m2_init=None):
    ks = np.exp(np.linspace(np.log(1e-5), np.log(2.0), 240)) if ks is None else ks
    nk = len(ks)
    u0, up0 = bd(ks, 12 * xi if m2_init is None else m2_init)
    wr = np.max(np.abs(u0 * np.conj(up0) - np.conj(u0) * up0 - 1j))
    y0 = np.concatenate([[1.0, np.sqrt(6.0) / b, 0.0, 0.0, 0.0], u0.real, u0.imag, up0.real, up0.imag])

    def rhs(eta, y):
        a, I, J, r = y[0], y[1], y[2], y[3]
        H = np.sqrt((0.5 * J * J + 0.5 * b * b * I * I + r) / 3.0)
        R = 2 * b * b * I * I - J * J
        d = np.empty_like(y)
        d[0] = a * a * H
        d[1] = a * J
        d[2] = a * (-(3 * H + gam) * J - b * b * I)
        d[3] = a * (gam * J * J - 4 * H * r)
        d[4] = a
        w2 = ks ** 2 + a * a * (xi - 1 / 6) * R
        o = 5
        d[o:o + nk] = y[o + 2 * nk:o + 3 * nk]
        d[o + nk:o + 2 * nk] = y[o + 3 * nk:o + 4 * nk]
        d[o + 2 * nk:o + 3 * nk] = -w2 * y[o:o + nk]
        d[o + 3 * nk:o + 4 * nk] = -w2 * y[o + nk:o + 2 * nk]
        return d

    def stop(eta, y):
        return y[4] - tau_end
    stop.terminal = True
    # conformal time to tau_end: radiation era a ~ 0.7 tau^0.5... eta(45) < 25
    etas = np.linspace(0, 30.0, n_eta)
    sol = solve_ivp(rhs, (0, 30.0), y0, method='DOP853', t_eval=etas, rtol=rtol, atol=atol, events=stop,
                    max_step=0.01)
    Y = sol.y
    a, I, J, r, tau = Y[0], Y[1], Y[2], Y[3], Y[4]
    H = np.sqrt((0.5 * J * J + 0.5 * b * b * I * I + r) / 3.0)
    R = 2 * b * b * I * I - J * J
    o = 5
    u2 = Y[o:o + nk] ** 2 + Y[o + nk:o + 2 * nk] ** 2
    spec = u2 / a ** 2                                        # |h_k|^2 (nk x nt)
    return dict(ks=ks, tau=tau, a=a, H=H, R=R, r=r, spec=spec, wronskian_err=float(wr), xi=xi, b=b, gam=gam)


def running_mean(tau, R, P=np.pi / 3):
    cs = np.concatenate([[0.0], np.cumsum(0.5 * (R[1:] + R[:-1]) * np.diff(tau))])
    out = np.empty_like(R)
    for i, t in enumerate(tau):
        lo, hi = max(t - P / 2, tau[0]), min(t + P / 2, tau[-1])
        out[i] = (np.interp(hi, tau, cs) - np.interp(lo, tau, cs)) / max(hi - lo, 1e-12)
    return out


def sigma2_below(ev, kc, i):
    ks = ev['ks']
    lk = np.log(ks)
    f = ks ** 3 / (2 * np.pi ** 2) * ev['spec'][:, i]
    m = ks < kc
    if m.sum() < 2:
        return 0.0
    s = np.trapezoid(f[m], lk[m])
    # partial last interval
    j = np.where(m)[0][-1]
    if j + 1 < len(ks):
        fk = np.interp(np.log(kc), lk[j:j + 2], f[j:j + 2])
        s += 0.5 * (f[j] + fk) * (np.log(kc) - lk[j])
    return float(s)


def sigma2_gauss(ev, kc, i):
    ks = ev['ks']
    f = ks ** 3 / (2 * np.pi ** 2) * ev['spec'][:, i] * np.exp(-(ks / kc) ** 2)
    return float(np.trapezoid(f, np.log(ks)))


def danger(ev, smi, H_GeV, norm='canon', N=60.0, cH=1.0, nf=120, region_factor=1.0, kc_factor=1.0,
           window='sharp', tau_th=np.inf):
    tau, a, H, R = ev['tau'], ev['a'], ev['H'], ev['R']
    xi = ev['xi']
    Rbar = running_mean(tau, R)
    G = smi.Gam(H_GeV)
    hmax = smi.flat['h_max']
    h0s = hmax * np.exp(np.linspace(0, np.log(20.0), nf))
    V2h0 = smi.d2V_dh0(h0s) / H_GeV ** 2                    # H_i^2 units
    if norm == 'canon':
        fbar = np.exp(G) * h0s / H_GeV                       # canonical field values
        V2 = np.exp(-2 * G) * V2h0
    else:                                                    # the track's convention
        fbar = h0s / H_GeV
        V2 = V2h0
    best = dict(D=-np.inf)
    for j in range(nf):
        mc2 = -V2[j] - xi * Rbar
        ok = np.where((mc2 > (cH * H) ** 2) & (tau < tau_th))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = kc_factor * a[i] * np.sqrt(mc2[i])
        s2 = sigma2_below(ev, kc, i) if window == 'sharp' else sigma2_gauss(ev, kc, i)
        if s2 <= 0:
            continue
        x = fbar[j] ** 2 / (2 * s2)
        D = (3 * N - x + np.log1p(x)) / LN10 + 3 * np.log10(kc) + np.log10(region_factor)
        if D > best['D']:
            best = dict(D=float(D), fbar_over_barrier=float(h0s[j] / hmax), tau_c=float(tau[i]), k_c=float(kc),
                        sigma2_c=s2, x=float(x), m_c2=float(mc2[i]))
    best['Gamma_HI'] = G
    best['V2_at_hmax_H2'] = float(V2[0])
    return best


def corners():
    return [(172.60, 1.29e11), (172.60, 1.51e11), (172.87, 1.29e11), (172.87, 1.51e11), (172.33, 1.51e11)]


def main(which='bench'):
    out_path = os.path.join(HERE, '..', 'results', 'sk_modes.json')
    out = json.load(open(out_path)) if os.path.exists(out_path) else {}
    t0 = time.time()
    if which == 'bench':
        rows = {}
        for Mt, Hg in corners():
            smi = SMIn(Mt)
            xh = smi.xi_H(Hg, 0.1)
            ev = evolve(xh)
            r = dict(xi_H=xh, wronskian_err=ev['wronskian_err'])
            # sigma^2 at tau = 0 for k < 1
            r['sigma2_L_tau0_canon'] = sigma2_below(ev, 1.0, 0)
            for norm in ('track', 'canon'):
                for cH in (0.5, 1.0, 2.0):
                    r[f'{norm} cH={cH}'] = danger(ev, smi, Hg, norm=norm, cH=cH)
            r['canon cH=1 region (2pi/kc)^3'] = danger(ev, smi, Hg, norm='canon', region_factor=1 / (2 * np.pi) ** 3)
            r['canon cH=1 N=61 + past-light-cone volume x36'] = danger(ev, smi, Hg, norm='canon', N=61.0,
                                                                     region_factor=36.0)
            rows[f'Mt={Mt} H={Hg:.3g}'] = r
            print(Mt, Hg, 'xiH %.5f' % xh, 'track D %.1f' % r['track cH=1.0']['D'], 'canon D %.1f' % r['canon cH=1.0']['D'],
                  'tau_c %.2f/%.2f' % (r['track cH=1.0']['tau_c'], r['canon cH=1.0']['tau_c']),
                  's2c %.3e/%.3e' % (r['track cH=1.0']['sigma2_c'], r['canon cH=1.0']['sigma2_c']),
                  '%.0fs' % (time.time() - t0), flush=True)
        out = json.load(open(out_path)) if os.path.exists(out_path) else {}
        out['benchmark_xiEW=0.1'] = rows
    elif which == 'edges':
        rows = out.get('lower_edges_N=60', {})
        sel = [int(v) for v in sys.argv[2:]] if len(sys.argv) > 2 else range(len(corners()))
        for Mt, Hg in [corners()[i] for i in sel]:
            smi = SMIn(Mt)
            cache = {}

            def D(xh, norm):
                if xh not in cache:
                    cache[xh] = evolve(xh)
                return danger(cache[xh], smi, Hg, norm=norm)['D']
            res = {}
            for norm in ('track', 'canon'):
                lo, hi = 0.015, 0.16
                dlo, dhi = D(lo, norm), D(hi, norm)
                if not (dlo > 0 > dhi):
                    res[norm] = dict(note='no sign change', D_lo=dlo, D_hi=dhi)
                    continue
                for _ in range(12):
                    mid = 0.5 * (lo + hi)
                    dm = D(mid, norm)
                    if dm > 0:
                        lo = mid
                    else:
                        hi = mid
                xe = 0.5 * (lo + hi)
                res[norm] = dict(xi_H_lower=xe, xi_EW_lower=float(smi.xi_EW(Hg, xe)))
            rows[f'Mt={Mt} H={Hg:.3g}'] = res
            print(Mt, Hg, res, '%.0fs' % (time.time() - t0), flush=True)
        out_path = os.path.join(HERE, '..', 'results', 'sk_edges_%s.json' % '_'.join(str(i) for i in sel))
        out = {'lower_edges_N=60': rows}
    elif which == 'criteria':
        rows = {}
        for Mt, Hg in corners():
            smi = SMIn(Mt)
            xh = smi.xi_H(Hg, 0.1)
            ev = evolve(xh)
            r = {}
            for norm in ('track', 'canon'):
                for f in (0.5, 1.0, 2.0):
                    r[f'{norm} k_c={f} a m_c sharp'] = danger(ev, smi, Hg, norm=norm, kc_factor=f)['D']
                r[f'{norm} gauss window R=1/k_c'] = danger(ev, smi, Hg, norm=norm, window='gauss')['D']
                r[f'{norm} gauss window R=1/(2k_c)'] = danger(ev, smi, Hg, norm=norm, window='gauss', kc_factor=2.0)['D']
                # plasma deadline: D with only tau_c < tau_th counted
                r[f'{norm} D(tau_th)'] = {f'{t}': danger(ev, smi, Hg, norm=norm, tau_th=t)['D'] for t in
                                          (5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0)}
            rows[f'Mt={Mt} H={Hg:.3g}'] = r
            print(Mt, Hg, json.dumps(r, default=float), '%.0fs' % (time.time() - t0), flush=True)
        out = json.load(open(out_path)) if os.path.exists(out_path) else {}
        out['criterion_sensitivity_xiEW=0.1'] = rows
    elif which == 'history':
        # sensitivity of the benchmark D to the supplied inflaton history (m_I, Gamma_I) = (3, Gamma)
        rows = {}
        for Mt, Hg in ((172.87, 1.29e11), (172.87, 1.51e11), (172.60, 1.51e11)):
            smi = SMIn(Mt)
            xh = smi.xi_H(Hg, 0.1)
            for gam in (0.1, 0.25, 0.5, 1.0, 2.0):
                ev = evolve(xh, gam=gam, tau_end=60.0)
                d = danger(ev, smi, Hg, norm='canon')
                dt = danger(ev, smi, Hg, norm='track')
                rows[f'Mt={Mt} H={Hg:.3g} Gamma={gam}'] = dict(canon=d, track=dt)
                print(Mt, Hg, gam, 'canon D %.1f tau_c %.1f' % (d['D'], d.get('tau_c', np.nan)),
                      'track D %.1f' % dt['D'], '%.0fs' % (time.time() - t0), flush=True)
        out = json.load(open(out_path)) if os.path.exists(out_path) else {}
        out['history_sensitivity'] = rows
    json.dump(out, open(out_path, 'w'), indent=1, default=float)
    print('total %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'bench')
