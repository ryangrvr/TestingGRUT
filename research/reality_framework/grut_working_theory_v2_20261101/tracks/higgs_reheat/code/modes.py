"""Item 3b of track higgs_reheat: k-resolved Higgs modes through the end of inflation, the oscillating-R era and
reheating, for (m_I, Gamma_I) = (3, 0.5) H_i.  Own code.  Writes ../results/modes.json.

Equations [D] (each of the 4 real Higgs components, H_i = 1, comoving k with a(0) = 1):
    h_k'' + 3 H h_k' + [k^2/a^2 + xi R(tau) + m_p^2(tau) + m_Hartree^2] h_k = 0,
    m_Hartree^2 = 6 lam_eff(mu) sigma^2   (V = lam (H^dag H)^2, per-component variance sigma^2, <h^2> = 4 sigma^2,
                                            mu = 2 sigma, lam_eff = flat 2-loop of smpot/sm_rge),
background (lna, I, I', rho_r) integrated in the same system (identical to background.py, which agrees with our
rpreheat.c engine to 1e-13).  Initial state: Bunch-Davies de Sitter modes at tau = 0 with m^2 = 12 xi H^2 (exact
Hankel functions, heavy fields included; Wronskian checked).
Variance definitions (stated, not unique):
    sigma_L^2 = int_{k<1} d^3k/(2pi)^3 |h_k|^2        long-wavelength part (modes outside the horizon at the end of
                                                      inflation: they carry the inflationary distribution);
    sigma_P^2 = int_{k>=1} d^3k/(2pi)^3 max(|h_k|^2 - 1/(2 a^3 w_k), 0),  w_k^2 = k^2/a^2 + max(M^2, 0) + H^2
                                                      (produced quanta, zeroth-order adiabatic subtraction).
The Hartree mass uses sigma^2 = sigma_L^2 + sigma_P^2.
Plasma mass m_p^2 = kappa_T T(tau)^2 (kappa_T from the couplings at mu = pi T; T from the actual rho_r(tau) with
g* = 106.75): tau_th = 0 means it follows rho_r from the start (prompt, instant thermalisation); tau_th > 0 multiplies it
by a smooth step (1 + tanh((tau - tau_th)/0.2))/2 (delayed thermalisation, a model choice).  tau_th = inf is the
purely gravitational (literature-standard, conservative) case.  With a plasma the integration stops at tau_th + 0.5;
after that m_p^2 >= 1e4 H_i^2 >> |xi R|, |V''| and every long mode only oscillates with decaying amplitude [A].
With a plasma, the zeroth-order subtraction in sigma_P is not meaningful (a fast mass rise squeezes sub-horizon
modes: n_k ~ m_p/4k without any growth of the field amplitude), so plasma runs report the unsubtracted long-wave
sigma_L (max over time) as their variance diagnostic.
Decay criteria:
    mean field: 2 sigma > h_max (flat 2-loop barrier) at any time, or Hartree runaway (event: 2 sigma > 30 h_max);
    tail at the end (crude, Hubble-patch): log10[e^{3N} e^{-x}(1+x)], x = h_max^2/(2 sigma_end^2);
    scale-resolved danger D (our criterion, [E]): a region with smoothed field hbar > h_max can roll only once
        m_c^2(hbar, tau) = -V''_flat(hbar) - xi <R>(tau) - m_p^2(tau) exceeds H(tau)^2, and only if it is larger than
        1/m_c.  At the first such time tau_c: k_c = a(tau_c) m_c, sigma_c^2 = int_{k<k_c} d^3k/(2pi)^3 |h_k(tau_c)|^2
        (all those modes are super-horizon then), N_reg = e^{3N} k_c^3 regions in our past light cone,
        D(hbar) = log10[N_reg e^{-x}(1+x)], x = hbar^2/(2 sigma_c^2).  Survival needs max_hbar D < 0.
        <R> is the running one-period mean of R.  A plasma mass (>= 1e4 H_i^2 >> |V''|) removes every hbar whose
        tau_c is later than tau_th.
"""
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

import smpot as S
from background import B_I, GAM, T_GeV

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


class Lam:
    def __init__(self, Mt):
        self.sm = S.SM(Mt)
        hs = np.exp(np.linspace(np.log(1e8), np.log(1e18), 400))
        self.f = interp1d(np.log(hs), self.sm.lam_eff_flat(hs), bounds_error=False,
                          fill_value=(float(self.sm.lam_eff_flat(np.array([1e8]))[0]), -0.01))
        self.flat = self.sm.flat_scales()
        # thermal mass coefficient m_T^2/T^2 = 3g^2/16 + g'^2/16 + yt^2/4 + lam/2  (h << T)
        Ts = np.exp(np.linspace(np.log(1e10), np.log(1e16), 61))
        k = []
        for T in Ts:
            c = self.sm.cpl(np.pi * T)
            k.append(3 * c['G2'] / 16 + 0.6 * c['G1'] / 16 + c['T'] / 4 + c['lam'] / 2)
        self.kT = interp1d(np.log(Ts), np.array(k), bounds_error=False, fill_value=(k[0], k[-1]))

    def __call__(self, h):
        return self.f(np.log(np.maximum(h, 1.0)))


def bd_modes(ks, m2):
    """Bunch-Davies de Sitter modes at tau = 0 (a = 1, H = 1): returns h, h' (cosmic time) arrays."""
    nu = mp.sqrt(mp.mpf(9) / 4 - mp.mpf(m2))
    pref = mp.sqrt(mp.pi) / 2 * mp.e ** (-mp.pi * mp.im(nu) / 2)
    h, hp = [], []
    for k in ks:
        k = mp.mpf(k)
        Hn = mp.hankel1(nu, k)
        Hd = mp.hankel1(nu - 1, k) - nu / k * Hn
        h.append(complex(pref * Hn))
        hp.append(complex(-pref * (mp.mpf(3) / 2 * Hn + k * Hd)))
    return np.array(h), np.array(hp)


def kgrid(Kmax, n=260):
    return np.exp(np.linspace(np.log(1e-4), np.log(Kmax), n))


def run(xi, Mt=172.60, H_GeV=1.29e11, tau_th=np.inf, tau_end=40.0, hartree=True, n=260, Kmax=None, lam=None,
        rtol=1e-8, atol=1e-12, nsave=400, cHs=(), plasma_extra=0.5):
    lam = Lam(Mt) if lam is None else lam
    Kmax = (6.0 + 4.0 * np.sqrt(12 * max(xi, 0))) if Kmax is None else Kmax
    ks = kgrid(Kmax, n)
    lk = np.log(ks)
    wq = np.gradient(lk) * ks ** 3 / (2 * np.pi ** 2)          # d^3k/(2pi)^3 = dlnk k^3/(2 pi^2)
    IR = ks < 1.0
    h0, hp0 = bd_modes(ks, 12 * xi)
    wr = (h0 * np.conj(hp0) - np.conj(h0) * hp0)
    wr_err = float(np.max(np.abs(wr - 1j)))
    hmax = lam.flat['h_max'] / H_GeV
    nk = len(ks)
    y0 = np.concatenate([[0.0, np.sqrt(6.0) / B_I, 0.0, 0.0], h0.real, h0.imag, hp0.real, hp0.imag])
    tend = min(tau_end, tau_th + plasma_extra)
    state = {'runaway': False}

    def parts(t, y):
        lna, I, J, r = y[:4]
        a = np.exp(lna)
        H = np.sqrt((0.5 * J * J + 0.5 * B_I ** 2 * I * I + r) / 3.0)
        R = 2 * B_I ** 2 * I * I - J * J
        hr, hi, pr, pi = y[4:4 + nk], y[4 + nk:4 + 2 * nk], y[4 + 2 * nk:4 + 3 * nk], y[4 + 3 * nk:]
        mod2 = hr * hr + hi * hi
        mp2 = 0.0
        if np.isfinite(tau_th) and r > 0:
            T = T_GeV(r, H_GeV) / H_GeV
            ramp = 1.0 if tau_th <= 0 else 0.5 * (1 + np.tanh((t - tau_th) / 0.2))
            mp2 = float(lam.kT(np.log(max(T * H_GeV, 1e10)))) * T * T * ramp
        M2c = xi * R + mp2
        w2 = ks ** 2 / a ** 2 + max(M2c, 0.0) + H * H
        sub = 1.0 / (2 * a ** 3 * np.sqrt(w2))
        sL = float(np.sum(wq[IR] * mod2[IR]))
        sP = float(np.sum(wq[~IR] * np.maximum(mod2[~IR] - sub[~IR], 0.0)))
        s2 = sL + sP
        mH = 6 * float(lam(2 * np.sqrt(s2) * H_GeV)) * s2 if hartree else 0.0
        return a, H, R, J, r, hr, hi, pr, pi, M2c + mH, sL, sP, mp2, mH

    def rhs(t, y):
        a, H, R, J, r, hr, hi, pr, pi, M2, sL, sP, mp2, mH = parts(t, y)
        w2 = ks ** 2 / a ** 2 + M2
        d = np.empty_like(y)
        d[0] = H
        d[1] = J
        d[2] = -(3 * H + GAM) * J - B_I ** 2 * y[1]
        d[3] = GAM * J * J - 4 * H * r
        d[4:4 + nk] = pr
        d[4 + nk:4 + 2 * nk] = pi
        d[4 + 2 * nk:4 + 3 * nk] = -3 * H * pr - w2 * hr
        d[4 + 3 * nk:] = -3 * H * pi - w2 * hi
        return d

    def blow(t, y):
        _, _, _, _, _, _, _, _, _, _, sL, sP, _, _ = parts(t, y)
        return 2 * np.sqrt(sL + sP) - 30 * hmax          # far beyond the barrier: runaway
    blow.terminal = True
    ts = np.linspace(0, tend, nsave)
    t0 = time.time()
    sol = solve_ivp(rhs, (0, tend), y0, method='DOP853', t_eval=ts, rtol=rtol, atol=atol, events=blow,
                    max_step=0.05)
    rows = []
    for j, t in enumerate(sol.t):
        a, H, R, J, r, hr, hi, pr, pi, M2, sL, sP, mp2, mH = parts(t, sol.y[:, j])
        rows.append((float(t), float(a), float(H), float(R), float(M2), float(sL), float(sP), float(mp2), float(mH)))
    rows = np.array(rows)
    runaway = len(sol.t_events[0]) > 0
    spec = sol.y[4:4 + nk, :] ** 2 + sol.y[4 + nk:4 + 2 * nk, :] ** 2        # |h_k|^2 at saved times
    dang = danger(ks, wq, rows, spec, xi, lam, H_GeV, tau_th)
    dang_c = {f'cH={c}': danger(ks, wq, rows, spec, xi, lam, H_GeV, tau_th, cH=c)['best'] for c in cHs}
    # final occupation numbers for sub-horizon modes
    yN = sol.y[:, -1]
    a, H, R, J, r, hr, hi, pr, pi, M2, sL, sP, mp2, mH = parts(sol.t[-1], yN)
    w = np.sqrt(np.maximum(ks ** 2 / a ** 2 + max(M2, 0.0), 1e-30))
    nkf = a ** 3 / (2 * w) * ((pr ** 2 + pi ** 2) + w ** 2 * (hr ** 2 + hi ** 2)) - 0.5
    sub_h = ks / a > 3 * H
    s2f = sL + sP
    sig0 = float(rows[0, 5] + rows[0, 6])
    x = hmax ** 2 / (2 * s2f)
    lp = float(mp.log10(mp.e ** (-x) * (1 + x)))
    res = dict(xi=xi, Mt=Mt, H_GeV=H_GeV, tau_th=(None if not np.isfinite(tau_th) else tau_th), tau_end=float(sol.t[-1]),
               hartree=hartree, nk=nk, Kmax=Kmax, wronskian_err=wr_err, runaway=bool(runaway),
               sigma2_0=sig0, sigma2_L_end=float(sL), sigma2_P_end=float(sP), sigma2_end=float(s2f),
               amplification_sigma2=float(s2f / sig0), max_sigma2=float(np.max(rows[:, 5] + rows[:, 6])),
               max_sigma2_L=float(np.max(rows[:, 5])), sigma2_L_0=float(rows[0, 5]),
               tau_at_max=float(rows[np.argmax(rows[:, 5] + rows[:, 6]), 0]),
               h_max_over_H=float(hmax), x_tail=float(x), log10_P_patch=lp,
               log10_expected_patches={N: lp + 3 * N / np.log(10) for N in (50, 55, 60)},
               min_M2=float(np.min(rows[:, 4])), max_nk_subhorizon=(float(np.max(nkf[sub_h])) if sub_h.any() else None),
               k_at_max_nk=(float(ks[sub_h][np.argmax(nkf[sub_h])]) if sub_h.any() else None),
               danger=dang, danger_cH=dang_c, seconds=time.time() - t0, nfev=int(sol.nfev))
    return res, rows


def danger(ks, wq, rows, spec, xi, lam, H_GeV, tau_th=np.inf, Ns=(50, 55, 60), nh=160, cH=1.0):
    """scale-resolved danger function (see module docstring).  rows: tau a H R M2 sL sP mp2 mH."""
    tau, a, H, R = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
    # running one-period mean of R (R oscillates at 2 m_I = 6 H_i: period pi/3)
    P = np.pi / 3
    Rbar = np.array([np.mean(R[(tau >= t - P / 2) & (tau <= t + P / 2)]) for t in tau])
    hmax = lam.flat['h_max']
    lh = np.linspace(np.log(hmax), np.log(hmax) + np.log(20.0), nh)
    hh = np.exp(lh)
    d = 1e-3
    Vf = lambda x: lam(x) * x ** 4 / 4
    V2 = (Vf(hh * (1 + d)) - 2 * Vf(hh) + Vf(hh * (1 - d))) / (hh * d) ** 2     # GeV^2
    out = {N: [] for N in Ns}
    best = {N: dict(D=-np.inf) for N in Ns}
    for j, hb in enumerate(hh):
        mc2 = -V2[j] / H_GeV ** 2 - xi * Rbar - (rows[:, 7] if True else 0.0)
        ok = np.where((mc2 > (cH * H) ** 2) & (tau < tau_th))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = a[i] * np.sqrt(mc2[i])
        m = ks < kc
        s2 = float(np.sum(wq[m] * spec[m, i]))
        if s2 <= 0:
            continue
        x = (hb / H_GeV) ** 2 / (2 * s2)
        lp = float(mp.log10(mp.e ** (-x) * (1 + x)))
        for N in Ns:
            D = lp + 3 * N / np.log(10) + 3 * np.log10(kc)
            if D > best[N]['D']:
                best[N] = dict(D=D, hbar_GeV=float(hb), hbar_over_hmax=float(hb / hmax), tau_c=float(tau[i]),
                               k_c=float(kc), sigma2_c=s2, x=float(x), m_c2=float(mc2[i]))
    first_tau_c = None
    for j, hb in enumerate(hh):
        mc2 = -V2[j] / H_GeV ** 2 - xi * Rbar - rows[:, 7]
        ok = np.where(mc2 > H ** 2)[0]
        if len(ok):
            first_tau_c = float(tau[ok[0]]) if first_tau_c is None else min(first_tau_c, float(tau[ok[0]]))
    return dict(best={N: best[N] for N in Ns}, earliest_tau_c=first_tau_c,
                V2_at_hmax_GeV2=float(V2[0]), sqrt_negV2_at_hmax_over_H=float(np.sqrt(max(-V2[0], 0)) / H_GeV))


def main(which='scan'):
    out = {}
    path = os.path.join(RES, 'modes.json')
    if os.path.exists(path):
        out = json.load(open(path))
    t0 = time.time()
    if which == 'scan':
        for Mt in (172.60, 172.87, 172.33):
            lam = Lam(Mt)
            for H_GeV in (1.29e11, 1.51e11):
                key = f'grav_only Mt={Mt} H={H_GeV:.3g}'
                rows = out.get(key, {})
                for xi in (0.05, 0.088, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0):
                    if f'{xi}' in rows:
                        continue
                    r, _ = run(xi, Mt=Mt, H_GeV=H_GeV, lam=lam)
                    rows[f'{xi}'] = r
                    print(key, xi, "amp %.3g" % r["amplification_sigma2"], "runaway", r["runaway"], "D60 %.1f" % r["danger"]["best"][60]["D"],
                          'log10N60 %.1f' % r['log10_expected_patches'][60], 'maxn %s' % r['max_nk_subhorizon'],
                          '%.0fs' % r['seconds'], flush=True)
                    out[key] = rows
                    json.dump(out, open(path, 'w'), indent=1, default=float)
    elif which == 'series':
        # time series for the benchmark xi(H_i) = 0.0878 (xi_EW = 0.1), gravitational only
        lam = Lam(172.60)
        r, rows = run(0.0878, lam=lam, nsave=800)
        out['series_xi=0.0878_Mt=172.6_H=1.29e11'] = dict(summary=r, cols='tau a H R M2 sigma2_L sigma2_P mp2 mHartree',
                                                          rows=rows[::2].tolist())
        r2, rows2 = run(10.0, lam=lam, nsave=800)
        out['series_xi=10_Mt=172.6_H=1.29e11'] = dict(summary=r2, cols='tau a H R M2 sigma2_L sigma2_P mp2 mHartree',
                                                      rows=rows2[::2].tolist())
        json.dump(out, open(path, 'w'), indent=1, default=float)
    elif which == 'plasma':
      out.pop('plasma Mt=172.6 H=1.29e11', None)        # superseded first version (abrupt ramp)
      for Mt, H_GeV in ((172.60, 1.29e11), (172.87, 1.51e11)):
        lam = Lam(Mt)
        pk = f'plasma Mt={Mt} H={H_GeV:.3g}'
        rows = out.get(pk, {})
        for tau_th in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
            for xi in (0.05, 0.088, 1.0, 10.0, 30.0, 50.0, 100.0):
                key = f'xi={xi},tau_th={tau_th}'
                if key in rows:
                    continue
                r, _ = run(xi, Mt=Mt, H_GeV=H_GeV, lam=lam, tau_th=tau_th)
                rows[key] = r
                print(pk, key, "runaway", r["runaway"], "maxsigL2/sigL2(0) %.3g" % (r["max_sigma2_L"] / r["sigma2_L_0"]),
                      "2sigLmax/hmax %.3g" % (2 * r["max_sigma2_L"] ** 0.5 / r["h_max_over_H"]),
                      "D60 %.1f" % r["danger"]["best"][60]["D"], '%.0fs' % r['seconds'], flush=True)
                out[pk] = rows
                json.dump(out, open(path, 'w'), indent=1, default=float)
    elif which == 'prompt_large':
        lam = Lam(172.87)
        rows = {}
        for xi in (100.0, 300.0, 1000.0):
            r, rr = run(xi, Mt=172.87, H_GeV=1.51e11, lam=lam, tau_th=0.0, plasma_extra=2.5, n=200, nsave=300)
            mH = rr[:, 4]
            rows[f'{xi}'] = dict(tau_end=r['tau_end'], max_sigma2_L_over_initial=r['max_sigma2_L'] / r['sigma2_L_0'],
                                 two_sigma_L_max_over_hmax=2 * r['max_sigma2_L'] ** 0.5 / r['h_max_over_H'],
                                 min_total_mass2=float(mH.min()), runaway=r['runaway'], seconds=r['seconds'])
            print('prompt', xi, rows[f'{xi}'], flush=True)
        out['prompt plasma, large xi, Mt=172.87 H=1.51e11, to tau=2.5'] = rows
        json.dump(out, open(path, 'w'), indent=1, default=float)
    print('total %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'scan')
