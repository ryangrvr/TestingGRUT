#!/usr/bin/env python3
"""mf_higgs SKEPTIC: Higgs modes after inflation on the hybrid histories and the danger exponent D, own code.

Different from the track (mh_modes.py: h_k, hdot_k, fixed-step RK4, track's Background):
  * variable u_k = a^{3/2} h_k, u'' + Omega_k^2 u = 0 with Omega_k^2 = k^2/a^2 + (xi - 1/4) R + (3/4) H^2 + M_x^2
    (exact identity -(9/4)H^2 - (3/2)Hdot = (3/4)H^2 - R/4, using R = 6(2H^2 + Hdot));
  * adaptive DOP853 (scipy) on the stacked real system, piecewise per stage (homogeneous raw-R stage, fluid stage),
    cubic-spline background from the history JSON (no rebuilt R lines: the track showed they move D by <= 0.04);
  * Bunch-Davies from scipy Hankel functions at t0 = -3 (exact de Sitter there);
  * my own implementation of the argued scale-resolved criterion D (higgs_reheat / inf_higgs, with the inf_higgs
    skeptic's abundance-limit grid x <= 3N + 40):
      a region of canonical value phibar >= phi_max rolls at the first save tau_c >= tau(eps_H = 1) with
      m_c^2 = -V''(phibar) - xi Rbar - m_p^2 > H^2;  k_c = a m_c;  sigma_c^2 = variance of the modes k < k_c;
      x = phibar^2/(2 sigma_c^2);  D = log10[e^{3N} k_c^3 e^{-x}(1 + x)], maximised over phibar.
  * Hartree in the inf_higgs convention (6 lambda_eff e^{-4 Gamma} sigma^2, adiabatic subtraction per mode).
Histories (model choices): 'EXP' = mf_hybrid bench_HFa as exported (K = 2.05); 'MC_P0' = mf_chain's rebuilt Pauli
history bench_itHFa_P0_it2 (Gamma_th = 0, tails free, K = 7.84, H_i = 2.40e11); 'MC_TR' = mf_chain's rebuilt unblocked
history bench_itHFa_track_it1 (K = 2.04, H_i = 1.284e11).
Usage: python3 -B sk_modes.py runs   -> the four scratch runs (EXP_base, MC_TR_Hi1.284e+11, MC_P0_Hi1.3e+11, MC_P0_Hi2.4e+11)
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import hankel1, h1vp

import sk_common as C

HPATH = {'EXP': C.TRACKS / 'mf_hybrid' / 'results' / 'bench_HFa.json',
         'EXPb': C.TRACKS / 'mf_hybrid' / 'results' / 'bench_HFb.json',
         'MC_P0': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json',
         'MC_TR': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_track_it1.json'}
NSTAR_BENCH, K_BENCH, H_BENCH = 52.644453, 2.0481855, 1.30002426e11


class Hist:
    def __init__(self, tag):
        d, h = C.load_hist(HPATH[tag])
        self.tag, self.d, self.h = tag, d, h
        self.Hi = float(d['H_i_GeV'])
        rh = d['reheating']
        self.K = float(rh['K'])
        self.t_sw = float(rh.get('t_sw', 0.6))
        self.te = float(d['waterfall']['t_eps1'])
        t = h['t']
        self.N_obs = NSTAR_BENCH - 0.5 * math.log(self.K / K_BENCH) + 0.5 * math.log(self.Hi / H_BENCH) + C.LN_K_H0
        # stage splines: homogeneous (raw R) and fluid (cycle-mean R); junction at t_sw
        m1 = t <= self.t_sw + 1e-9
        m2 = t >= self.t_sw - 1e-9
        self.sp = []
        for m in (m1, m2):
            tt = t[m]
            keep = np.concatenate([[True], np.diff(tt) > 1e-12])
            tt = tt[keep]
            self.sp.append(dict(t0=tt[0], t1=tt[-1], lna=CubicSpline(tt, h['N'][m][keep]),
                                H=CubicSpline(tt, h['H'][m][keep]), R=CubicSpline(tt, h['R'][m][keep])))
        # criterion mean R: running mean over +-0.015 H^-1 in the raw stage, the table (cycle mean) after
        self.tab_t = t
        R = h['R'].copy()
        Rbar = R.copy()
        raw = m1 & (t > self.te - 0.05)
        tr, Rr = t[m1], R[m1]
        cR = np.concatenate([[0.0], np.cumsum(0.5 * (Rr[1:] + Rr[:-1]) * np.diff(tr))])
        for i in np.where(raw)[0]:
            lo, hi = max(t[i] - 0.015, tr[0]), min(t[i] + 0.015, tr[-1])
            Rbar[i] = (np.interp(hi, tr, cR) - np.interp(lo, tr, cR)) / (hi - lo)
        self.Rbar_tab = Rbar

    def bg(self, t, stage):
        s = self.sp[stage]
        return math.exp(float(s['lna'](t))), float(s['H'](t)), float(s['R'](t))

    def Rbar(self, t):
        return np.interp(t, self.tab_t, self.Rbar_tab)


def bd(k, xi):
    """Bunch-Davies u = a^{3/2} h at a = e^{t0}, H = 1, m^2 = 12 xi (exact de Sitter)."""
    nu = math.sqrt(9 / 4 - 12 * xi)
    return lambda x: (math.sqrt(math.pi) / 2 * hankel1(nu, x), -math.sqrt(math.pi) / 2 * x * h1vp(nu, x))


def run(tag, xi, Hi=None, t0=-3.0, t_end=80.0, n=240, kmin=1e-4, Kmax=10.0, hartree=True, Mt_hartree=172.57,
        m2x=None, rtol=1e-8, label=''):
    tw = time.time()
    H = Hist(tag)
    Hi = H.Hi if Hi is None else Hi
    fl = C.Flat(Mt_hartree)
    eG = math.exp(fl.Gam(Hi))
    ks = np.exp(np.linspace(math.log(kmin), math.log(Kmax), n))
    w = np.gradient(np.log(ks)) * ks ** 3 / (2 * math.pi ** 2)
    a0, H0, R0 = H.bg(t0, 0)
    x0 = ks / (a0 * H0)
    if 12 * xi < 9 / 4 - 1e-6:
        u0, ud0 = bd(ks, xi)(x0)
    else:   # heavy field (m^2 = 12 xi H^2 > 9/4 H^2): adiabatic vacuum, accurate for m >> H
        Om2 = ks ** 2 / a0 ** 2 + 12 * xi * H0 ** 2 - 9 / 4 * H0 ** 2
        Om = np.sqrt(Om2)
        dOm = -(ks ** 2 / a0 ** 2) * H0 / Om
        u0 = (1 / np.sqrt(2 * Om)).astype(complex)
        ud0 = (-1j * Om - dOm / (2 * Om)) * u0
    y0 = np.concatenate([u0.real, u0.imag, ud0.real, ud0.imag])
    saves = np.unique(np.concatenate([np.arange(-1.0, H.te - 0.2, 0.1), np.arange(H.te - 0.2, H.te + 0.6, 0.005),
                                      np.arange(H.te + 0.6, t_end, 0.04), [t_end]]))
    saves = saves[(saves > t0) & (saves <= t_end)]

    def rhs_factory(stage):
        def rhs(t, y):
            a, Hh, R = H.bg(t, stage)
            ur, ui, vr, vi = y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:]
            M2 = xi * R + (m2x(t) if m2x else 0.0)
            mh = 0.0
            if hartree:
                om = np.sqrt(ks * ks / (a * a) + max(M2, 0.0) + Hh * Hh)
                hh2 = (ur * ur + ui * ui) / a ** 3
                s2 = float(np.dot(w, np.maximum(hh2 - 1 / (2 * a ** 3 * om), 0.0)))
                if s2 > 0:
                    mh = 6 * float(fl.lam(max(2 * math.sqrt(s2) * Hi / eG, 1e8))) * eG ** -4 * s2
            Om2 = ks * ks / (a * a) + (xi - 0.25) * R + 0.75 * Hh * Hh + M2 - xi * R + mh
            # note: M2 already contains xi R; Om2 = k^2/a^2 + (xi - 1/4) R + 3/4 H^2 + m2x + mh
            return np.concatenate([vr, vi, -Om2 * ur, -Om2 * ui])
        return rhs

    rec_t, rec_u2, rec_bg = [], [], []
    y = y0
    segs = [(t0, min(H.t_sw, t_end), 0), (H.t_sw, t_end, 1)]
    nfev = 0
    for (ta, tb, st) in segs:
        if tb <= ta:
            continue
        te_ = saves[(saves > ta) & (saves <= tb)]
        sol = solve_ivp(rhs_factory(st), (ta, tb), y, method='DOP853', rtol=rtol, atol=1e-12, t_eval=te_,
                        dense_output=False, max_step=0.05)
        nfev += sol.nfev
        for j, tt in enumerate(sol.t):
            yy = sol.y[:, j]
            a, Hh, R = H.bg(tt, st)
            rec_t.append(tt)
            rec_u2.append(((yy[:n] ** 2 + yy[n:2 * n] ** 2) / a ** 3).astype(np.float64))
            rec_bg.append((a, Hh, R, float(H.Rbar(tt))))
        y = sol.y[:, -1]
    a, Hh, R = H.bg(t_end, 1)
    ur, ui, vr, vi = y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:]
    wr = (ur * vi - ui * vr) * 2   # u u'* - u* u' = -2i Im(u* u')... check |.| = 1
    hk = (ur + 1j * ui) / a ** 1.5
    hd = ((vr + 1j * vi) - 1.5 * Hh * (ur + 1j * ui)) / a ** 1.5
    Rb_end = float(H.Rbar(t_end))
    om_end = np.sqrt(ks ** 2 / a ** 2 + max(xi * Rb_end, 0.0) + Hh ** 2)
    nk = a ** 3 / (2 * om_end) * (np.abs(hd) ** 2 + om_end ** 2 * np.abs(hk) ** 2) - 0.5
    out = dict(tag=tag, xi=xi, Hi=Hi, N_obs=H.N_obs, K=H.K, te=H.te, ks=ks, w=w, t=np.array(rec_t), nk_end=nk,
               h2=np.array(rec_u2).T, bg=np.array(rec_bg), nfev=nfev, seconds=time.time() - tw,
               wronskian_err=float(np.max(np.abs(np.abs(wr) - 1))))
    return out


def save(o, name):
    p = C.SP / 'runs' / (name + '.npz')
    p.parent.mkdir(exist_ok=True)
    np.savez_compressed(p, **{k: v for k, v in o.items() if k not in ('tag',)}, tag=o['tag'])
    return p


def load(name):
    z = np.load(C.SP / 'runs' / (name + '.npz'), allow_pickle=True)
    o = {k: z[k] for k in z.files}
    for k in ('xi', 'Hi', 'N_obs', 'K', 'te'):
        o[k] = float(o[k])
    o['tag'] = str(o['tag'])
    return o


def danger(o, Mt, N=None, Mh=C.PDG['Mh'], als=C.PDG['als'], mp2=None, cH=1.0, kcf=1.0, nh=300, x_extra=40.0,
           tau_th=np.inf, Hi=None):
    """own implementation of the scale-resolved criterion.  mp2: array of plasma mass^2 at the saves (or None)."""
    N = o['N_obs'] if N is None else N
    Hi = o['Hi'] if Hi is None else Hi
    fl = C.Flat(Mt, Mh, als)
    eG = math.exp(fl.Gam(Hi))
    t = o['t']
    a, H, Rb = o['bg'][:, 0], o['bg'][:, 1], o['bg'][:, 3]
    ks, w = o['ks'], o['w']
    cum = np.cumsum(w[:, None] * o['h2'], axis=0)
    post = t >= o['te']
    s2max = float(np.max(cum[-1, post]))
    xlim = 3 * N + x_extra
    phimax = eG * fl.h_max / Hi
    philim = math.sqrt(2 * xlim * s2max)
    hm = max(philim / phimax, 1.5)
    phis = phimax * np.exp(np.linspace(0, math.log(hm), nh))
    V2 = np.array([fl.V2_canon(p * Hi / eG, Hi) for p in phis])
    mp2 = np.zeros(len(t)) if mp2 is None else mp2
    best = dict(D=-np.inf)
    nroll = 0
    for j, ph in enumerate(phis):
        mc2 = -V2[j] - o['xi'] * Rb - mp2
        ok = np.where((mc2 > (cH * H) ** 2) & post & (t < tau_th))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = kcf * a[i] * math.sqrt(mc2[i])
        jj = int(np.searchsorted(ks, kc)) - 1
        if jj < 0:
            continue
        nroll += 1
        s2 = float(cum[jj, i])
        x = ph ** 2 / (2 * s2)
        Dv = (math.log1p(x) - x + 3 * N) / math.log(10) + 3 * math.log10(kc)
        if Dv > best['D']:
            best = dict(D=Dv, phi_over_phimax=ph / phimax, tau_c=float(t[i]), k_c=kc, sigma2_c=s2, x=x, m_c2=float(mc2[i]))
    best['n_rolled'] = nroll
    best['phimax_canon'] = phimax
    return best


def crit(Ms, Ds):
    for i in range(len(Ms) - 1):
        if np.isfinite(Ds[i]) and np.isfinite(Ds[i + 1]) and Ds[i] < 0 <= Ds[i + 1]:
            return float(Ms[i] + (0 - Ds[i]) * (Ms[i + 1] - Ms[i]) / (Ds[i + 1] - Ds[i]))
    return None


def make_runs():
    """the four runs used by sk_postD / sk_medium / sk_infl (gaussian tail), with fixed-dimensionless N_obs at a rescaled H_i"""
    for tag, Hi, name in (('EXP', None, 'EXP_base'), ('MC_TR', None, 'MC_TR_Hi1.284e+11'),
                          ('MC_P0', 1.30002426e11, 'MC_P0_Hi1.3e+11'), ('MC_P0', None, 'MC_P0_Hi2.4e+11')):
        Hh = Hist(tag)
        Hi_ = Hh.Hi if Hi is None else Hi
        xi = C.Flat(172.57).xi(Hi_, 0.1)
        o = run(tag, xi, Hi=Hi_)
        o['N_obs'] = NSTAR_BENCH - 0.5 * math.log(Hh.K / K_BENCH) + 0.5 * math.log(Hi_ / H_BENCH) + C.LN_K_H0
        save(o, name)
        print(name, 'xi %.5f N_obs %.3f K %.2f nfev %d %.1fs wronskian %.1e' % (xi, o['N_obs'], o['K'], o['nfev'], o['seconds'],
                                                                                o['wronskian_err']), flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'runs':
        make_runs()
