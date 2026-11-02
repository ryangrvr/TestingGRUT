#!/usr/bin/env python3
"""Track mf_higgs: k-resolved Higgs modes on the hybrid backgrounds and the scale-resolved danger exponent D (own code).

Equations [D] (each real Higgs component, canonical field phi = e^{Gamma(H_i)} h0, H_i = 1, a = 1 at sigma = sigma_c):
    h_k'' + 3 H h_k' + [k^2/a^2 + xi R(t) + m_x^2(t) + 6 lam_eff e^{-4 Gamma} sigma^2] h_k = 0.
  * Start at t0 < 0 in the exact de Sitter phase of the hybrid (H = H_i to 1e-7, R = 12 H^2 there): exact
    Bunch-Davies modes for every k (including the super-horizon ones, which is exact for constant H and m^2).
  * R(t): raw (exact homogeneous solution of mf_hybrid) through the waterfall, then the fluid mean plus the rebuilt
    2 m_chi and 2 m_sigma lines (mh_common.Background).  Fixed-step vectorised RK4 with a step cap that resolves
    k/a and the R lines.
  * Hartree as in inf_higgs (6 lam_eff(2 sigma) e^{-4 Gamma} sigma^2, adiabatic subtraction for every mode).
  * m_x^2(t): optional extra Higgs mass in the modes (portal variants; plasma-in-modes variant).
Danger D [E] (higgs_reheat / inf_higgs scale-resolved criterion, canonical normalisation): a region with smoothed
h0 = hbar >= h_max rolls once m_c^2 = -V''(hbar) - xi <R> - <m_x^2> - m_p^2 > (cH H)^2; at that first time tau_c >= tau_min
(tau_min = eps_H = 1), k_c = a m_c, sigma_c^2 = sum_{k < k_c} |h_k|^2 d^3k/(2 pi)^3, x = (e^Gamma hbar)^2 / (2 sigma_c^2),
D(hbar) = log10[e^{3N} k_c^3 e^{-x} (1 + x)], maximised over hbar.
N1 fix (inf_higgs skeptic): the hbar grid runs from h_max up to the abundance limit x = 3N + 40 (not 20 h_max), and a
run in which no hbar inside that limit rolls returns D = -inf flagged 'no roll within the abundance limit'.
"""
import sys
sys.dont_write_bytecode = True
import math
import time

import mpmath as mp
import numpy as np
from scipy.special import hankel1

import mh_common as C


def bd(ks_over_aH, m2_over_H2):
    """Bunch-Davies de Sitter modes at a = 1, H = 1: (h, hdot) for x = k/(aH)."""
    nu2 = 9 / 4 - m2_over_H2
    x = np.asarray(ks_over_aH, float)
    if nu2 > 1e-12:
        nu = math.sqrt(nu2)
        Hn = hankel1(nu, x)
        Hd = hankel1(nu - 1, x) - nu / x * Hn
        c = math.sqrt(math.pi) / 2
        return c * Hn, -c * (1.5 * Hn + x * Hd)
    nu = mp.sqrt(mp.mpf(nu2))
    pref = mp.sqrt(mp.pi) / 2 * mp.e ** (-mp.pi * mp.im(nu) / 2)
    h, hp = [], []
    for k in x:
        k = mp.mpf(float(k))
        Hn = mp.hankel1(nu, k)
        Hd = mp.hankel1(nu - 1, k) - nu / k * Hn
        h.append(complex(pref * Hn))
        hp.append(complex(-pref * (mp.mpf(3) / 2 * Hn + k * Hd)))
    return np.array(h), np.array(hp)


def run(bg, xi, smt, t0=-3.0, t_end=80.0, n=260, Kmax=None, kmin=1e-4, hartree=True, m2x=None, cK=0.05,
        dt_max=0.02, save_dt=0.04, keep_spec=True):
    """Evolve the Higgs modes on Background bg.  smt: ih_common.SMTab (Hartree lambda, Gamma).  m2x: callable t -> extra
    Higgs mass^2 (H_i^2) in the mode equations.  Returns dict(res, rows, spec, ks, wq)."""
    tw = time.time()
    H_GeV = bg.Hi_GeV
    Kmax = (6.0 + 4.0 * math.sqrt(12 * max(xi, 0))) if Kmax is None else Kmax
    ks = np.exp(np.linspace(math.log(kmin), math.log(Kmax), n))
    wq = np.gradient(np.log(ks)) * ks ** 3 / (2 * math.pi ** 2)
    eG = math.exp(smt.G(H_GeV))
    hmax_c = eG * smt.flat['h_max'] / H_GeV
    a0, H0, R0 = bg.a_H_R(t0)
    mx0 = m2x(t0) if m2x else 0.0
    h, p = bd(ks / (a0 * H0), (xi * R0 + mx0) / H0 ** 2)
    h = h * a0 ** -1.5 * H0 ** -0.5
    p = p * a0 ** -1.5 * H0 ** 0.5
    wr = h * np.conj(p) - np.conj(h) * p
    wr_err = float(np.max(np.abs(wr * a0 ** 3 - 1j)))
    lamf = smt.lam

    def accel(t, h, p):
        a, H, R = bg.a_H_R(t)
        M2 = xi * R + (m2x(t) if m2x else 0.0)
        s2 = 0.0
        mH = 0.0
        if hartree:
            w2 = ks * ks / (a * a) + max(M2, 0.0) + H * H
            mod2 = h.real ** 2 + h.imag ** 2
            s2 = float(np.dot(wq, np.maximum(mod2 - 1.0 / (2 * a ** 3 * np.sqrt(w2)), 0.0)))
            if s2 > 0:
                mH = 6 * float(lamf(2 * math.sqrt(s2) * H_GeV)) * eG ** -4 * s2
        w2 = ks * ks / (a * a) + M2 + mH
        return -3 * H * p - w2 * h, (a, H, R, M2, s2, mH)

    t = t0
    rows, spec = [], []
    next_save = t0
    runaway = False
    nstep = 0
    while t < t_end - 1e-12:
        a, H, _ = bg.a_H_R(t)
        dt = min(cK * a / Kmax, bg.dt_limit(t), dt_max, t_end - t)
        k1p, info = accel(t, h, p)
        if t >= next_save - 1e-12:
            aa, HH, RR, M2, s2, mH = info
            rows.append((t, aa, HH, RR, float(bg.Rbar(t)), M2, s2, mH, float(bg.plasma_m2(t))))
            if keep_spec:
                spec.append((h.real ** 2 + h.imag ** 2).astype(np.float64))
            near = (bg.tau_endinf - 0.2 < t < bg.tau_endinf + 0.6)
            next_save = t + (0.005 if near else (save_dt if t > 0 else 0.1))
            if 2 * math.sqrt(max(s2, 0)) > 30 * hmax_c:
                runaway = True
                break
        k1h = p
        k2h = p + 0.5 * dt * k1p
        k2p, _ = accel(t + 0.5 * dt, h + 0.5 * dt * k1h, k2h)
        k3h = p + 0.5 * dt * k2p
        k3p, _ = accel(t + 0.5 * dt, h + 0.5 * dt * k2h, k3h)
        k4h = p + dt * k3p
        k4p, _ = accel(t + dt, h + dt * k3h, k4h)
        h = h + dt / 6 * (k1h + 2 * k2h + 2 * k3h + k4h)
        p = p + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)
        t += dt
        nstep += 1
    rows = np.array(rows)
    a, H, R = bg.a_H_R(t)
    wr = h * np.conj(p) - np.conj(h) * p
    # occupation numbers at the end (adiabatic, mean curvature mass)
    w = np.sqrt(ks ** 2 / a ** 2 + max(xi * float(bg.Rbar(t)), 0.0) + H * H)
    nk = a ** 3 / (2 * w) * (np.abs(p) ** 2 + w ** 2 * np.abs(h) ** 2) - 0.5
    res = dict(xi=xi, H_GeV=H_GeV, n=n, Kmax=Kmax, t0=t0, t_end=float(t), nstep=nstep, runaway=runaway,
               wronskian_err0=wr_err, wronskian_err_end=float(np.max(np.abs(wr * a ** 3 - 1j))),
               tau_endinf=bg.tau_endinf, max_sigma2=float(np.max(rows[:, 6])), seconds=time.time() - tw,
               nk_end=nk.tolist(), hist=bg.hist, name=bg.name, N_obs=bg.N_obs, K=bg.K)
    return dict(res=res, rows=rows, spec=(np.array(spec).T if keep_spec else None), ks=ks, wq=wq)


def danger(R, smt, xi, N, mp2=None, m2xbar=None, cH=1.0, kc_f=1.0, tau_min=None, tau_th=np.inf, nh=300,
           x_extra=40.0, hmult_min=20.0, H_GeV=None):
    """scale-resolved danger exponent (see module docstring).  rows: t a H R Rbar M2 s2 mH plasma_m2.
    mp2: plasma mass^2 at the saves (None = gravitational only); 'bg' = the stored hard-gas plasma column."""
    rows, spec, ks, wq = R['rows'], R['spec'], R['ks'], R['wq']
    H_GeV = R['res']['H_GeV'] if H_GeV is None else H_GeV
    tau, a, H, Rb = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 4]
    tau_min = R['res']['tau_endinf'] if tau_min is None else tau_min
    if mp2 is None:
        mp2 = np.zeros(len(tau))
    elif isinstance(mp2, str) and mp2 == 'bg':
        mp2 = rows[:, 8]
    m2xbar = np.zeros(len(tau)) if m2xbar is None else m2xbar
    eG = math.exp(smt.G(H_GeV))
    hmax = smt.flat['h_max']
    cum = np.cumsum(wq[:, None] * spec, axis=0)
    post = tau >= tau_min
    s2max = float(np.max(cum[-1, post])) if post.any() else float(np.max(cum[-1]))
    x_lim = 3 * N + x_extra
    h_lim = math.sqrt(2 * x_lim * s2max) * H_GeV / eG
    hm = max(hmult_min, h_lim / hmax)
    hh = hmax * np.exp(np.linspace(0.0, math.log(hm), nh))
    d = 1e-3
    Vf = lambda x: smt.lam(x) * x ** 4 / 4
    V2 = (Vf(hh * (1 + d)) - 2 * Vf(hh) + Vf(hh * (1 - d))) / (hh * d) ** 2 / eG ** 2 / H_GeV ** 2   # canonical, H_i^2
    best = dict(D=-np.inf)
    rolled = 0
    for j, hb in enumerate(hh):
        mc2 = -V2[j] - xi * Rb - m2xbar - mp2
        ok = np.where((mc2 > (cH * H) ** 2) & (tau >= tau_min) & (tau < tau_th))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = kc_f * a[i] * math.sqrt(mc2[i])
        jj = int(np.searchsorted(ks, kc)) - 1
        if jj < 0:
            continue
        rolled += 1
        s2 = float(cum[jj, i])
        x = (eG * hb / H_GeV) ** 2 / (2 * s2)
        Dv = (math.log1p(x) - x) / math.log(10) + 3 * math.log10(kc) + 3 * N / math.log(10)
        if Dv > best['D']:
            best = dict(D=Dv, hbar_over_hmax=hb / hmax, tau_c=float(tau[i]), k_c=kc, sigma2_c=s2, x=x, m_c2=float(mc2[i]))
    best['n_rolled'] = rolled
    best['hbar_grid_max_over_hmax'] = hm
    best['note'] = 'no roll within the abundance limit' if rolled == 0 else ''
    # bound for -inf: the largest D a roll at the grid edge could have had
    best['D_bound_if_none'] = (math.log1p(x_lim) - x_lim) / math.log(10) + 3 * N / math.log(10) + 3 * math.log10(max(Kmax_of(R), 1.0))
    return best


def Kmax_of(R):
    return float(R['ks'][-1])


def crit(Mts, Ds):
    """first zero crossing of D(M_t), increasing M_t"""
    for i in range(len(Mts) - 1):
        a, b = Ds[i], Ds[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a < 0 <= b:
            return float(Mts[i] + (0 - a) * (Mts[i + 1] - Mts[i]) / (b - a))
        if (not np.isfinite(a)) and np.isfinite(b) and b >= 0:
            return float(Mts[i + 1])
    return None
