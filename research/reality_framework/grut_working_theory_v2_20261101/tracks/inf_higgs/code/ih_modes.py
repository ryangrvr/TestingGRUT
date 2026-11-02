#!/usr/bin/env python3
"""Track inf_higgs: k-resolved Higgs modes on the inflaton benchmarks (own code; generalises
tracks/higgs_reheat/skeptic/code/modes_canon.py, i.e. higgs_reheat's modes.py with the skeptic's canonical-normalisation
fix, to an explicit inflaton potential and a start deep in the de Sitter phase).

Equations [D] (each real Higgs component, canonical field phi = e^{Gamma(H_i)} h0, H_i = 1, a(S1) = 1):
    h_k'' + 3 H h_k' + [k^2/a^2 + xi R(tau) + m_x^2(tau) + 6 lam_eff e^{-4 Gamma} sigma^2] h_k = 0,
background (lna, phi, phidot, rho_r, rho_m) integrated in the same system (ih_common.Background equations).
  * Start: tau0 = Nb0 e-folds before eps_H = 1 (Nb0 = 5 default; H(tau0) = 1.0003 H_i), exact Bunch-Davies de Sitter
    modes with the LOCAL H(tau0), a(tau0) and m^2 = xi R(tau0).  From there the modes see the benchmark's own H(N), R(N)
    through the last e-folds (R/H^2 falls from 12 to 6 at eps_H = 1), the end of inflation and the anharmonic oscillation.
    For the M1r control: tau0 = 0, a = H = 1, m^2 = 12 xi (exactly higgs_reheat's setup).
  * m_x^2: optional extra Higgs masses: R-even portal (lam_IH/2) phi_I^2 |H|^2 -> m^2 = lam_IH phi_I^2/2 (canonical
    doublet, inf_reheat skeptic's convention), and an R-odd trilinear A phi_I |H|^2 -> m^2 = A phi_I (both in H_i units).
  * Hartree: sigma^2 = sigma_L^2 (k < k_split) + sigma_P^2 (k >= k_split), k_split = a H at the end of inflation; both with
    the zeroth-order adiabatic subtraction for the benchmarks (the start is deep in inflation, where k < k_split modes are
    still sub-horizon; the subtraction is negligible once they are super-horizon); the M1r control keeps modes_canon's
    convention (k < 1 unsubtracted) so that it reproduces the skeptic's numbers.  lam_eff is the flat 2-loop value at mu = 2 sigma (h0 units)
    times e^{-4 Gamma(H_i)} (canonical; skeptic's fix).
Danger D [E] (higgs_reheat's scale-resolved criterion, canonical normalisation): a region with smoothed h0 = hbar >= h_max
rolls once m_c^2 = -d2V/dphi2(hbar) - xi <R> - <m_x^2> - m_p^2 > (cH H)^2; at that first time tau_c, k_c = a m_c,
sigma_c^2 = sum_{k<k_c} |h_k|^2 d^3k/(2pi)^3, D(hbar) = log10[e^{3N} k_c^3 e^{-x}(1+x)], x = (e^Gamma hbar)^2/(2 sigma_c^2).
N counts e-folds from the exit of the observable horizon to the a = 1 slice.  Only tau_c >= tau_min counts (benchmarks:
tau_min = time of eps_H = 1; the inflationary phase is treated by ih_inflation.py).  <R> is the running mean over one local
oscillation period of the benchmark's inflaton (fixed window pi/3 for the M1r control, as in modes_canon).
A plasma onset at tau_th removes every hbar with tau_c >= tau_th (sharp; the skeptic's convention); a supplied smooth
m_p^2(tau) enters m_c^2 instead.  The deadline curve D(tau_th) is returned for every run.
"""
import sys
sys.dont_write_bytecode = True
import math
import time

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

import ih_common as C


def bd_modes(x, m2_over_H2):
    """Bunch-Davies de Sitter modes at a = 1, H = 1 for k = x: (h, hdot)."""
    nu = mp.sqrt(mp.mpf(9) / 4 - mp.mpf(m2_over_H2))
    pref = mp.sqrt(mp.pi) / 2 * mp.e ** (-mp.pi * mp.im(nu) / 2)
    h, hp = [], []
    for k in x:
        k = mp.mpf(float(k))
        Hn = mp.hankel1(nu, k)
        Hd = mp.hankel1(nu - 1, k) - nu / k * Hn
        h.append(complex(pref * Hn))
        hp.append(complex(-pref * (mp.mpf(3) / 2 * Hn + k * Hd)))
    return np.array(h), np.array(hp)


def kgrid(Kmax, n=260, kmin=1e-4):
    return np.exp(np.linspace(math.log(kmin), math.log(Kmax), n))


def run(bg, xi, sm, H_GeV, n=260, Kmax=None, rtol=1e-8, atol=1e-13, nsave=900, hartree=True, lamIH=0.0, Aterm=0.0,
        tau_end=None, max_step=None, keep_spec=True):
    """Evolve the Higgs modes on Background bg.  Returns (res dict, rows array, spec array, ks, wq)."""
    t0w = time.time()
    tau_end = bg.tau_end if tau_end is None else tau_end
    Kmax = (6.0 + 4.0 * math.sqrt(12 * max(xi, 0))) if Kmax is None else Kmax
    ks = kgrid(Kmax, n)
    wq = np.gradient(np.log(ks)) * ks ** 3 / (2 * math.pi ** 2)
    nk = len(ks)
    y0b = bg.y0.copy()
    if bg.kind == 'M1r':
        a0, H0, R0 = 1.0, 1.0, 12.0
        ksplit = 1.0
    else:
        a0 = math.exp(y0b[0]); H0 = bg.H_of(y0b); R0 = bg.R_of(y0b)
        Fe = bg.fields(np.array([bg.tau_endinf]))
        ksplit = float(Fe['a'][0] * Fe['H'][0])
    m2x0 = lamIH * (y0b[1] - bg.phimin) ** 2 / (2 * C_eta2(bg)) + Aterm * (y0b[1] - bg.phimin) / math.sqrt(C_eta2(bg))
    h0, hp0 = bd_modes(ks / (a0 * H0), (xi * R0 + m2x0) / H0 ** 2)
    h0 = h0 * a0 ** -1.5 * H0 ** -0.5
    hp0 = hp0 * a0 ** -1.5 * H0 ** 0.5
    wr = h0 * np.conj(hp0) - np.conj(h0) * hp0
    wr_err = float(np.max(np.abs(wr * a0 ** 3 - 1j)))
    IR = ks < ksplit
    subIR = (bg.kind == 'bench')     # deep-inflation start: IR modes are sub-horizon at tau0, subtract them too
    eG = math.exp(sm.G(H_GeV))
    eta2 = C_eta2(bg)
    eta = math.sqrt(eta2)
    hmax_c = eG * sm.flat['h_max'] / H_GeV
    V, dV, gam, phimin = bg.V, bg.dV, bg.gam, bg.phimin

    def parts(t, y):
        lna, f, fd, rr, rm = y[:5]
        a = math.exp(lna)
        Vf = V(f)
        H = math.sqrt(max((0.5 * fd * fd + Vf + rr + rm) / 3.0, 1e-300))
        R = 4 * Vf - fd * fd + rm
        hr, hi = y[5:5 + nk], y[5 + nk:5 + 2 * nk]
        mod2 = hr * hr + hi * hi
        fI = f - phimin
        m2x = (lamIH * fI * fI / (2 * eta2) if lamIH else 0.0) + (Aterm * fI / eta if Aterm else 0.0)
        M2c = xi * R + m2x
        if hartree:
            w2 = ks ** 2 / a ** 2 + max(M2c, 0.0) + H * H
            sub = 1.0 / (2 * a ** 3 * np.sqrt(w2))
            sL = float(np.dot(wq[IR], np.maximum(mod2[IR] - sub[IR], 0.0) if subIR else mod2[IR]))
            sP = float(np.dot(wq[~IR], np.maximum(mod2[~IR] - sub[~IR], 0.0)))
            s2 = sL + sP
            mH = 6 * float(sm.lam(2 * math.sqrt(max(s2, 1e-300)) * H_GeV)) * eG ** -4 * s2
        else:
            sL = sP = mH = 0.0
        return a, H, R, rr, rm, M2c + mH, sL, sP, mH, m2x, fI

    def rhs(t, y):
        a, H, R, rr, rm, M2, sL, sP, mH, m2x, fI = parts(t, y)
        lna, f, fd = y[0], y[1], y[2]
        w2 = ks ** 2 / (a * a) + M2
        d = np.empty_like(y)
        d[0] = H
        d[1] = fd
        d[2] = -(3 * H + gam) * fd - dV(f)
        d[3] = gam * fd * fd + gam * rm - 4 * H * rr
        d[4] = -(3 * H + gam) * rm
        d[5:5 + 2 * nk] = y[5 + 2 * nk:5 + 4 * nk]
        d[5 + 2 * nk:5 + 3 * nk] = -3 * H * y[5 + 2 * nk:5 + 3 * nk] - w2 * y[5:5 + nk]
        d[5 + 3 * nk:] = -3 * H * y[5 + 3 * nk:] - w2 * y[5 + nk:5 + 2 * nk]
        return d

    def blow(t, y):
        s = parts(t, y)
        return 2 * math.sqrt(s[6] + s[7]) - 30 * hmax_c
    blow.terminal = True
    y0 = np.concatenate([y0b, h0.real, h0.imag, hp0.real, hp0.imag])
    te = bg.tau_endinf
    tsave = np.unique(np.concatenate([np.linspace(bg.tau0, min(te + 0.5, tau_end), 120),
                                      np.linspace(min(te + 0.5, tau_end), tau_end, nsave)]))
    ms = max_step if max_step else (0.05 if bg.kind == 'M1r' else 0.02)
    segs_t = [(bg.tau0, tau_end)]
    if bg.frag is not None and bg.frag < tau_end:
        segs_t = [(bg.tau0, bg.frag), (bg.frag, tau_end)]
    Ys, Ts = [], []
    y = y0
    runaway = False
    nfev = 0
    for j, (ta, tb) in enumerate(segs_t):
        if j == 1:   # fragmentation: homogeneous inflaton -> w = 0 fluid with the same energy
            y = y.copy()
            y[4] += 0.5 * y[2] ** 2 + V(y[1]); y[1] = phimin; y[2] = 0.0
            ms = 0.05
        te_ = tsave[(tsave >= ta) & (tsave <= tb)]
        sol = solve_ivp(rhs, (ta, tb), y, method='DOP853', t_eval=te_, rtol=rtol, atol=atol, events=blow,
                        max_step=ms)
        nfev += sol.nfev
        Ys.append(sol.y); Ts.append(sol.t)
        y = sol.y[:, -1]
        if len(sol.t_events[0]):
            runaway = True
            break
    Yall = np.concatenate(Ys, axis=1)
    Tall = np.concatenate(Ts)
    _, iu = np.unique(Tall, return_index=True)
    Yall, Tall = Yall[:, iu], Tall[iu]
    rows = []
    for j, t in enumerate(Tall):
        a, H, R, rr, rm, M2, sL, sP, mH, m2x, fI = parts(t, Yall[:, j])
        rows.append((t, a, H, R, M2, sL, sP, 0.0, mH, m2x, fI, rr, rm))
    rows = np.array(rows)
    spec = Yall[5:5 + nk, :] ** 2 + Yall[5 + nk:5 + 2 * nk, :] ** 2
    # occupation numbers at the final time (adiabatic there: xi <R> and its oscillation are small against k^2/a^2 for the
    # produced modes), with the period-mean curvature mass
    yN = Yall[:, -1]
    aN = math.exp(yN[0]); HN = bg.H_of(yN[:5])
    RbN = float(bg.Rbar(np.array([Tall[-1]]))[0]) if bg.kind == 'bench' else float(rows[-1, 3])
    wN = np.sqrt(ks ** 2 / aN ** 2 + max(xi * RbN, 0.0) + HN * HN)
    hrN, hiN = yN[5:5 + nk], yN[5 + nk:5 + 2 * nk]
    prN, piN = yN[5 + 2 * nk:5 + 3 * nk], yN[5 + 3 * nk:]
    nkN = aN ** 3 / (2 * wN) * ((prN ** 2 + piN ** 2) + wN ** 2 * (hrN ** 2 + hiN ** 2)) - 0.5
    res_nk = dict(nk_end=nkN.tolist(), Rbar_end=RbN)
    res = dict(xi=xi, H_GeV=H_GeV, nk=nk, Kmax=Kmax, k_split=ksplit, wronskian_err=wr_err, runaway=runaway,
               tau0=bg.tau0, tau_endinf=bg.tau_endinf, tau_end=float(Tall[-1]), sigma2_tau0=float(rows[0, 5] + rows[0, 6]),
               sigma2_L_endinf=float(np.interp(te, rows[:, 0], rows[:, 5])),
               sigma2_end=float(rows[-1, 5] + rows[-1, 6]), max_sigma2=float(np.max(rows[:, 5] + rows[:, 6])),
               min_M2=float(np.min(rows[:, 4])), lamIH=lamIH, Aterm=Aterm, seconds=time.time() - t0w, nfev=int(nfev), **res_nk)
    return res, rows, (spec if keep_spec else None), ks, wq


def sigma2_produced(res, rows, Rb, ks, wq, xi):
    """variance carried by the quanta counted at the final time, mapped back to earlier (adiabatic) times:
    sigma^2(tau) = sum wq n_k / (a^3 w_k), w_k^2 = k^2/a^2 + max(xi <R>, 0) + H^2 [E: adiabatic-era estimate]."""
    nk = np.maximum(np.array(res['nk_end']), 0.0)
    a, H = rows[:, 1], rows[:, 2]
    out = np.empty(len(rows))
    for i in range(len(rows)):
        w = np.sqrt(ks ** 2 / a[i] ** 2 + max(xi * Rb[i], 0.0) + H[i] ** 2)
        out[i] = float(np.sum(wq * nk / (a[i] ** 3 * w)))
    return out


def C_eta2(bg):
    return (bg.bench['eta'] ** 2) if bg.kind == 'bench' else (1.29e11 / C.MPL) ** 2


def rbar_rows(bg, rows):
    if bg.kind == 'M1r':
        return bg.Rbar(rows[:, 0], mode='fixed', P=math.pi / 3)
    return bg.Rbar(rows[:, 0], mode='period')


def mean_extra(bg, rows, lamIH=0.0, Aterm=0.0):
    """period-mean of the extra Higgs mass (portal and A-term) at the save times"""
    if not lamIH and not Aterm:
        return np.zeros(len(rows))
    F = bg.fine()
    eta2 = C_eta2(bg)
    m2 = lamIH * F['phiI'] ** 2 / (2 * eta2) + Aterm * F['phiI'] / math.sqrt(eta2)
    # reuse Rbar machinery: temporarily average m2 with the same windows
    t = F['tau']
    c = np.concatenate([[0.0], np.cumsum(0.5 * (m2[1:] + m2[:-1]) * np.diff(t))])
    tt = F['turn_tau']
    out = np.empty(len(rows))
    per = (tt[2:] - tt[:-2]) if len(tt) >= 3 else None
    for i, x in enumerate(rows[:, 0]):
        if per is None or x < tt[1]:
            out[i] = np.interp(x, t, m2)
        else:
            w = np.interp(x, tt[1:-1], per)
            lo, hi = max(x - w / 2, t[0]), min(x + w / 2, t[-1])
            out[i] = (np.interp(hi, t, c) - np.interp(lo, t, c)) / (hi - lo)
    return out


def danger(rows, spec, ks, wq, Rbar, xi, sm, H_GeV, Ns=(50, 55, 60), nh=160, cH=1.0, kc_f=1.0, mp2=None,
           m2xbar=None, tau_th=np.inf, tau_min=-np.inf,
           tau_ths=(1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 40, 60, 80, 120)):
    """Danger exponent (module docstring).  rows: tau a H R M2 sL sP mp2 mH m2x phiI rr rm."""
    tau, a, H = rows[:, 0], rows[:, 1], rows[:, 2]
    mp2 = np.zeros(len(tau)) if mp2 is None else mp2
    m2xbar = np.zeros(len(tau)) if m2xbar is None else m2xbar
    hmax = sm.flat['h_max']
    hh = np.exp(np.linspace(math.log(hmax), math.log(hmax) + math.log(20.0), nh))
    d = 1e-3
    Vf = lambda x: sm.lam(x) * x ** 4 / 4
    eG = math.exp(sm.G(H_GeV))
    V2 = (Vf(hh * (1 + d)) - 2 * Vf(hh) + Vf(hh * (1 - d))) / (hh * d) ** 2 / eG ** 2      # d2V/dphi2, GeV^2
    per_h = []
    for j, hb in enumerate(hh):
        mc2 = -V2[j] / H_GeV ** 2 - xi * Rbar - m2xbar - mp2
        ok = np.where((mc2 > (cH * H) ** 2) & (tau < tau_th) & (tau >= tau_min))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = kc_f * a[i] * math.sqrt(mc2[i])
        m = ks < kc
        s2 = float(np.sum(wq[m] * spec[m, i]))
        if s2 <= 0:
            continue
        x = (eG * hb / H_GeV) ** 2 / (2 * s2)
        lp = float(mp.log10(mp.e ** (-x) * (1 + x)))
        per_h.append((float(hb), float(tau[i]), float(kc), s2, float(x), float(mc2[i]), lp + 3 * math.log10(kc)))
    best = {}
    for N in Ns:
        if not per_h:
            best[N] = dict(D=-np.inf)
            continue
        Ds = [p[6] + 3 * N / math.log(10) for p in per_h]
        j = int(np.argmax(Ds))
        p = per_h[j]
        best[N] = dict(D=Ds[j], hbar_over_hmax=p[0] / hmax, tau_c=p[1], k_c=p[2], sigma2_c=p[3], x=p[4], m_c2=p[5])
    curve = {}
    for tt in tau_ths:
        sel = [p for p in per_h if p[1] < tt]
        curve[str(tt)] = {N: (max(p[6] for p in sel) + 3 * N / math.log(10)) if sel else -np.inf for N in Ns}
    first = min((p[1] for p in per_h), default=None)
    return dict(best=best, deadline_curve=curve, earliest_tau_c=first,
                sqrt_negV2_at_hmax_over_H=float(math.sqrt(max(-V2[0], 0)) / H_GeV), h_max_GeV=hmax,
                h_max_canon_over_H=eG * hmax / H_GeV)


def deadline(curve, N, target=0.0):
    """tau_th at which D(tau_th) crosses target (linear interpolation on the curve); None if never."""
    items = sorted(((float(k), v[N]) for k, v in curve.items()), key=lambda z: z[0])
    ts = [z[0] for z in items]
    Ds = [z[1] for z in items]
    prev = None
    for t, D in zip(ts, Ds):
        if D >= target:
            if prev is None or not np.isfinite(prev[1]):
                return t
            return prev[0] + (target - prev[1]) * (t - prev[0]) / (D - prev[1])
        prev = (t, D)
    return None
