#!/usr/bin/env python3
"""Skeptic of track inf_higgs (key inf_higgs_skeptic; Claude, adversarial pass, 3 November 2026).  Own code.

Independent re-implementation of the Higgs-mode / danger-exponent computation on the inflaton benchmarks.  What is
shared with the track (on purpose, it is the SM input the brief says not to redo): higgs_reheat/code/smpot.py
(3-loop running, flat 2-loop lambda_eff, Gamma(mu), the SMDR matching), imported read-only.  Everything else is
re-written with a different formulation:
  * potentials written out from the formulas of inf_potential NOTES section 1 (parameters read from bench_*.json);
  * background pre-run in cosmic time with LSODA from slow roll ~12 e-folds before the end (track: e-fold time + DOP853);
  * modes in the friction-free variable u_k = a^{3/2} h_k:  u'' + [k^2/a^2 + M^2 - 9H^2/4 - 3 Hdot/2] u = 0,
    integrated together with the background by a vectorised fixed-step classical RK4 (track: h_k with 3H h', DOP853);
  * Bunch-Davies start from scipy Hankel functions of real order (track: mpmath);
  * own danger function (same argued criterion, written from higgs_reheat SKEPTIC's definition) and own period mean.
Units: fields in reduced Planck units, time in 1/H_i, a(S1) = 1 (S1 = one e-fold before eps_H = 1).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import hankel1

HERE = Path(__file__).resolve().parent
SKD = HERE.parent
RES = SKD / 'results'
TRACK = SKD.parent
TRACKS = TRACK.parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_higgs_skeptic')
SP.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(TRACKS / 'higgs_reheat' / 'code'))
import smpot as S          # noqa: E402  read-only SM input

MPL = 2.435e18
LN10 = math.log(10.0)


# --------------------------------------------------------------------------------------------- SM input
class SMIn:
    """lambda_eff(h0) incl. e^{4 Gamma} (flat 2 loop), V(h0) spline in ln h0, Gamma(mu), h_max, xi running."""

    def __init__(self, Mt, Mh=125.20, als=0.1180):
        self.Mt, self.Mh, self.als = Mt, Mh, als
        key = SP / 'sm' / f'sm_{Mt:.4f}_{Mh:.3f}_{als:.5f}.json'
        key.parent.mkdir(exist_ok=True)
        if key.exists():
            d = json.load(open(key))
        else:
            sm = S.SM(Mt, Mh=Mh, als=als)
            lh = np.linspace(math.log(1e8), math.log(1e17), 1800)
            le = sm.lam_eff_flat(np.exp(lh))
            lm = np.linspace(math.log(1e2), math.log(1e18), 400)
            cp = [sm.cpl(math.exp(x)) for x in lm]
            d = dict(lh=lh.tolist(), le=le.tolist(), lm=lm.tolist(), Gam=[float(c['Gam']) for c in cp],
                     Ixi=[float(c['Ixi']) for c in cp], lam=[float(c['lam']) for c in cp],
                     G1=[float(c['G1']) for c in cp], G2=[float(c['G2']) for c in cp], T=[float(c['T']) for c in cp],
                     Ixi_Mt=float(sm.cpl(Mt)['Ixi']))
            json.dump(d, open(key, 'w'))
        self.d = d
        lh = np.array(d['lh'])
        le = np.array(d['le'])
        self.le = CubicSpline(lh, le)
        self.Vs = CubicSpline(lh, le * np.exp(4 * lh) / 4)
        lm = np.array(d['lm'])
        self.sp = {k: CubicSpline(lm, np.array(d[k])) for k in ('Gam', 'Ixi', 'lam', 'G1', 'G2', 'T')}
        # own h_max: maximum of V(h0) (spline root of V_l)
        hs = np.exp(np.linspace(math.log(1e9), math.log(1e15), 6001))
        V = self.Vs(np.log(hs))
        k = int(np.argmax(V))
        # refine with the derivative
        from scipy.optimize import brentq
        f = lambda x: float(self.Vs(x, 1))
        self.h_max = math.exp(brentq(f, math.log(hs[max(k - 2, 0)]), math.log(hs[min(k + 2, len(hs) - 1)])))
        j = np.where(le < 0)[0][0]
        self.Lambda_I = float(np.exp(np.interp(0.0, [le[j], le[j - 1]], [lh[j], lh[j - 1]])))

    def Gam(self, mu):
        return float(self.sp['Gam'](math.log(mu)))

    def lam_run(self, mu):
        return float(self.sp['lam'](math.log(mu)))

    def d2V(self, h0):
        """d^2V/dh0^2 (GeV^2) from the spline in ln h0"""
        lh = np.log(h0)
        return (self.Vs(lh, 2) - self.Vs(lh, 1)) / h0 ** 2

    def xi_at(self, mu, xi_EW):
        return 1 / 6 + (xi_EW - 1 / 6) * math.exp(float(self.sp['Ixi'](math.log(mu))) - self.d['Ixi_Mt'])

    def xi_EW_of(self, mu, xi_mu):
        return 1 / 6 + (xi_mu - 1 / 6) * math.exp(-(float(self.sp['Ixi'](math.log(mu))) - self.d['Ixi_Mt']))


_SMC = {}


def sm_in(Mt, Mh=125.20, als=0.1180):
    k = (round(Mt, 4), round(Mh, 3), round(als, 5))
    if k not in _SMC:
        _SMC[k] = SMIn(Mt, Mh, als)
    return _SMC[k]


# --------------------------------------------------------------------------------------------- potentials (own)
class Pot:
    """v(phi) shape with V0 = 1 from the formulas (inf_potential NOTES section 1); phi in M_P."""

    def __init__(self, name, iH=0):
        d = json.load(open(TRACKS / 'inf_potential' / 'results' / f'bench_{name}.json'))
        r = d['by_Hi'][iH]
        self.name, self.rec = name, r
        p = r['params']
        self.H_GeV = r['H_i_GeV']
        self.eta = self.H_GeV / MPL
        self.amp = r['V0_MP4'] / self.eta ** 2          # V in units H_i^2 M_P^2
        self.m_I = r['m_I_over_Hi']
        if name.startswith('HT'):
            self.kind = 'HT'
            self.vev, self.pw, self.beta = p['vev'], p['p'], p['beta']
            b, pw = self.beta, self.pw
            self.xtop = -(b / pw) ** (1.0 / (pw - 1))
            self.ftop = 1 - b * self.xtop - self.xtop ** pw
            from scipy.optimize import brentq
            self.xmin = brentq(lambda x: 1 - b * x - x ** pw, 0.5, 1.0)
            self.phimin = self.xmin * self.vev
        elif name.startswith('T'):
            self.kind = 'T'; self.M = math.sqrt(6 * p['alpha']); self.phimin = 0.0
        elif name.startswith('E'):
            self.kind = 'E'; self.M = math.sqrt(1.5 * p['alpha']); self.phimin = 0.0
        else:
            raise ValueError(name)

    def V(self, f):
        if self.kind == 'HT':
            x = f / self.vev
            return self.amp * (1 - self.beta * x - x ** self.pw) ** 2 / self.ftop ** 2
        if self.kind == 'T':
            return self.amp * np.tanh(f / self.M) ** 2
        return self.amp * (-np.expm1(-f / self.M)) ** 2

    def dV(self, f):
        if self.kind == 'HT':
            x = f / self.vev
            g = 1 - self.beta * x - x ** self.pw
            return self.amp * 2 * g * (-self.beta - self.pw * x ** (self.pw - 1)) / self.vev / self.ftop ** 2
        if self.kind == 'T':
            y = f / self.M
            return self.amp * 2 * np.tanh(y) / np.cosh(y) ** 2 / self.M
        e = np.exp(-f / self.M)
        return self.amp * 2 * (1 - e) * e / self.M


class Quad:
    """M1r surrogate: V = b^2 phi^2/2 released at rest at a = H = 1"""

    def __init__(self, b=3.0, H_GeV=1.29e11):
        self.b, self.kind, self.name, self.m_I = b, 'Q', 'M1r', b
        self.phimin = 0.0
        self.H_GeV = H_GeV
        self.eta = H_GeV / MPL

    def V(self, f):
        return 0.5 * self.b ** 2 * f * f

    def dV(self, f):
        return self.b ** 2 * f


# --------------------------------------------------------------------------------------------- background pre-run
def initial_state(pot, Nb0=5.0):
    """slow-roll start ~12 e-folds before the end (from the bench history table), LSODA in cosmic time with ln a;
    returns (tau0, y0) at Nb0 e-folds before eps_H = 1, normalised so that a(S1 = end - 1 e-fold) = 1, tau(S1) = 0,
    plus the end time."""
    hist = pot.rec['history']
    row = max(hist, key=lambda r: r['N_before_end'] if r['N_before_end'] <= 14 else -1)
    f0 = row['phi_MP']
    H0 = math.sqrt(pot.V(f0) / 3)
    fd0 = -pot.dV(f0) / (3 * H0)

    def rhs(t, y):
        la, f, fd = y
        H = math.sqrt((0.5 * fd * fd + pot.V(f)) / 3)
        return [H, fd, -3 * H * fd - pot.dV(f)]

    def eps1(t, y):
        H2 = (0.5 * y[2] ** 2 + pot.V(y[1])) / 3
        return 0.5 * y[2] ** 2 / H2 - 1.0
    eps1.terminal = True; eps1.direction = 1
    s = solve_ivp(rhs, (0, 1e4), [0.0, f0, fd0], method='LSODA', rtol=1e-11, atol=1e-14, dense_output=True,
                  events=eps1)
    te = float(s.t_events[0][0]); lae = float(s.y_events[0][0][0])
    # times of lna = lae - 1 (S1) and lae - Nb0 (start)
    from scipy.optimize import brentq
    tS1 = brentq(lambda t: s.sol(t)[0] - (lae - 1.0), 0, te)
    t0 = brentq(lambda t: s.sol(t)[0] - (lae - Nb0), 0, te)
    y0 = s.sol(t0).copy()
    y0[0] -= (lae - 1.0)
    return dict(tau0=t0 - tS1, y0=y0, tau_end_inf=te - tS1, N_pre=lae)


# --------------------------------------------------------------------------------------------- modes (RK4, u = a^{3/2} h)
def bd_u(ks_phys_over_H, m2_over_H2, H):
    """Bunch-Davies de Sitter: u = (sqrt(pi)/2) H^{-1/2} H_nu(x), du/dt = -(sqrt(pi)/2) H^{1/2} x H_nu'(x), x = k/(aH)."""
    nu = math.sqrt(9 / 4 - m2_over_H2)
    x = ks_phys_over_H
    Hn = hankel1(nu, x)
    Hd = hankel1(nu - 1, x) - nu / x * Hn
    c = math.sqrt(math.pi) / 2
    return c * H ** -0.5 * Hn, -c * H ** 0.5 * x * Hd


def evolve(pot, xi, smi=None, H_GeV=None, gam=0.0, Nb0=5.0, tau_end=80.0, n=240, kmin=1e-4, Kmax=None, nper=160,
           hartree=True, Aterm=0.0, lamIH=0.0, frag_after_end=None, frag_w=0.0, dsave=0.04, cK=0.05):
    """Higgs modes on the benchmark background (bench) or the M1r surrogate (Quad).
    Aterm, lamIH in H_i units (A phi_I in H_i^2 with phi_I in M_P: A_Hi = A_GeV * M_P / H_i^2 ... see A_to_code).
    frag_after_end: switch the homogeneous inflaton into a fluid with w = frag_w (0 or 1/3) at tau_end_inf + value."""
    H_GeV = pot.H_GeV if H_GeV is None else H_GeV
    Kmax = (6.0 + 4.0 * math.sqrt(12 * max(xi, 0))) if Kmax is None else Kmax
    ks = np.exp(np.linspace(math.log(kmin), math.log(Kmax), n))
    dlk = np.gradient(np.log(ks))
    wq = dlk * ks ** 3 / (2 * math.pi ** 2)
    if pot.kind == 'Q':
        y = np.array([0.0, math.sqrt(6.0) / pot.b, 0.0, 0.0, 0.0])     # lna, phi, phidot, rho_r, rho_f
        tau0, te_inf = 0.0, 0.0
    else:
        ini = initial_state(pot, Nb0)
        y = np.array([ini['y0'][0], ini['y0'][1], ini['y0'][2], 0.0, 0.0])
        tau0, te_inf = ini['tau0'], ini['tau_end_inf']
    fmin = pot.phimin
    eta2 = pot.eta ** 2
    eG = math.exp(smi.Gam(H_GeV)) if smi is not None else 1.0
    wf = [frag_w]
    fragged = [False]

    def bgq(y):
        la, f, fd, rr, rf = y
        Vf = pot.V(f)
        rho = 0.5 * fd * fd + Vf + rr + rf
        H = math.sqrt(rho / 3)
        R = 4 * Vf - fd * fd + (1 - 3 * wf[0]) * rf
        Hd = -0.5 * (fd * fd + 4 * rr / 3 + (1 + wf[0]) * rf)
        return H, R, Hd

    def extra(f):
        fI = f - fmin
        return (lamIH * fI * fI / (2 * eta2) if lamIH else 0.0) + (Aterm * fI / pot.eta if Aterm else 0.0)

    def deriv(y, U, P, sig2_holder):
        la, f, fd, rr, rf = y
        a = math.exp(la)
        H, R, Hd = bgq(y)
        M2 = xi * R + extra(f)
        if hartree and smi is not None:
            mod2 = (U.real ** 2 + U.imag ** 2) / a ** 3
            w2 = ks ** 2 / a ** 2 + max(M2, 0.0) + H * H
            s2 = float(np.dot(wq, np.maximum(mod2 - 1 / (2 * a ** 3 * np.sqrt(w2)), 0.0)))
            mu = max(2 * math.sqrt(max(s2, 1e-300)) * H_GeV, 1e3)
            M2 += 6 * float(smi.le(math.log(mu))) * eG ** -4 * s2
            sig2_holder[0] = s2
        Om2 = ks ** 2 / a ** 2 + M2 - 2.25 * H * H - 1.5 * Hd
        dy = np.array([H, fd, -(3 * H + gam) * fd - pot.dV(f), gam * fd * fd + gam * rf - 4 * H * rr,
                       -(3 * (1 + wf[0]) * H + gam) * rf])
        return dy, P, -Om2 * U, M2

    # initial modes
    a0 = math.exp(y[0])
    H0, R0, _ = bgq(y)
    m20 = (xi * R0 + extra(y[1])) / H0 ** 2
    U, P = bd_u(ks / (a0 * H0), m20, H0)
    W = U * np.conj(P) - np.conj(U) * P
    wr_err = float(np.max(np.abs(W - 1j)))
    t = tau0
    dt_osc = 2 * math.pi / (max(pot.m_I, 1.0) * nper)
    saves, spec = [], []
    fine = []
    next_save = tau0
    sig = [0.0]
    runaway = False
    hmax_c = (eG * smi.h_max / H_GeV) if smi is not None else 1e9
    nstep = 0
    while t < tau_end - 1e-12:
        a = math.exp(y[0])
        dt = min(cK * a / Kmax, dt_osc, tau_end - t)
        if frag_after_end is not None and (not fragged[0]) and t >= te_inf + frag_after_end:
            rho_I = 0.5 * y[2] ** 2 + pot.V(y[1])
            y = y.copy(); y[4] += rho_I; y[1] = fmin; y[2] = 0.0
            fragged[0] = True
        # RK4
        k1y, k1U, k1P, M2 = deriv(y, U, P, sig)
        if t >= next_save - 1e-12:
            H, R, Hd = bgq(y)
            saves.append((t, a, H, R, M2, sig[0], y[1] - fmin, y[2], y[3], y[4]))
            spec.append(((U.real ** 2 + U.imag ** 2) / a ** 3).astype(np.float64))
            next_save = t + (dsave if t > te_inf + 0.3 else min(dsave, 0.02))
            if 2 * math.sqrt(max(sig[0], 0)) > 30 * hmax_c:
                runaway = True
                break
        k2y, k2U, k2P, _ = deriv(y + 0.5 * dt * k1y, U + 0.5 * dt * k1U, P + 0.5 * dt * k1P, sig)
        k3y, k3U, k3P, _ = deriv(y + 0.5 * dt * k2y, U + 0.5 * dt * k2U, P + 0.5 * dt * k2P, sig)
        k4y, k4U, k4P, _ = deriv(y + dt * k3y, U + dt * k3U, P + dt * k3P, sig)
        y = y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
        U = U + dt / 6 * (k1U + 2 * k2U + 2 * k3U + k4U)
        P = P + dt / 6 * (k1P + 2 * k2P + 2 * k3P + k4P)
        t += dt
        nstep += 1
        if nstep % 2 == 0:
            H, R, Hd = bgq(y)
            fine.append((t, R, y[1] - fmin, y[2], H, extra(y[1])))
    rows = np.array(saves)
    fine = np.array(fine)
    W = U * np.conj(P) - np.conj(U) * P
    return dict(ks=ks, wq=wq, rows=rows, spec=np.array(spec).T, fine=fine, tau0=tau0, tau_end_inf=te_inf,
                wronskian_err0=wr_err, wronskian_err_end=float(np.max(np.abs(W - 1j))), runaway=runaway,
                nstep=nstep, xi=xi, H_GeV=H_GeV, gam=gam, kind=pot.kind, name=pot.name,
                frag=(te_inf + frag_after_end) if frag_after_end is not None else None)


def period_mean(ev, col=1, fixed=None):
    """running mean of fine[:, col] (R or extra mass) over one local inflaton period, centred.  The local period comes
    from successive maxima of phi_I after the end of inflation; before that (inflation) the instantaneous value.
    fixed: a fixed window (M1r control: pi/3)."""
    F = ev['fine']
    t, x = F[:, 0], F[:, col]
    c = np.concatenate([[0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t))])
    tq = ev['rows'][:, 0]
    out = np.interp(tq, t, x)
    if fixed is not None:
        w = np.full(len(tq), fixed)
    else:
        fI = F[:, 2]
        ii = np.where((fI[1:-1] > fI[:-2]) & (fI[1:-1] >= fI[2:]) & (t[1:-1] > ev['tau_end_inf']))[0] + 1
        if ev.get('frag') is not None:
            ii = ii[t[ii] < ev['frag']]
        if len(ii) < 2:
            return out
        tm = 0.5 * (t[ii[1:]] + t[ii[:-1]])
        per = np.diff(t[ii])
        w = np.interp(tq, tm, per)
        w = np.where(tq < ev['tau_end_inf'], 0.0, w)
        if ev.get('frag') is not None:
            w = np.where(tq > ev['frag'], 0.0, w)
    for i, (x0, ww) in enumerate(zip(tq, w)):
        if ww > 0:
            lo, hi = max(x0 - ww / 2, t[0]), min(x0 + ww / 2, t[-1])
            if hi > lo:
                out[i] = (np.interp(hi, t, c) - np.interp(lo, t, c)) / (hi - lo)
    return out


def danger(ev, smi, Rbar, Ns=(60.0,), cH=1.0, kc_f=1.0, mp2=None, m2xbar=None, nh=200, tau_min=None, tau_th=np.inf,
           H_GeV=None, hmult=20.0):
    """argued scale-resolved criterion (higgs_reheat skeptic's definition, canonical normalisation):
    for each hbar >= h_max: first tau >= tau_min with m_c^2 = -V''_phi(hbar) - xi <R> - <m_x^2> - m_p^2 > (cH H)^2;
    k_c = kc_f a m_c; sigma_c^2 = sum_{k<k_c} |h_k|^2 d^3k/(2 pi)^3; x = (e^Gamma hbar)^2/(2 sigma_c^2);
    D = log10[e^{3N} k_c^3 e^{-x}(1+x)], maximised over hbar."""
    H_GeV = ev['H_GeV'] if H_GeV is None else H_GeV
    rows, spec, ks, wq = ev['rows'], ev['spec'], ev['ks'], ev['wq']
    tau, a, H = rows[:, 0], rows[:, 1], rows[:, 2]
    tau_min = ev['tau_end_inf'] if tau_min is None else tau_min
    mp2 = np.zeros(len(tau)) if mp2 is None else mp2
    m2xbar = np.zeros(len(tau)) if m2xbar is None else m2xbar
    eG = math.exp(smi.Gam(H_GeV))
    hh = smi.h_max * np.exp(np.linspace(0.0, math.log(hmult), nh))
    V2c = smi.d2V(hh) / eG ** 2 / H_GeV ** 2                 # canonical V'' in H_i^2
    cum = np.cumsum(wq[:, None] * spec, axis=0)               # sigma^2(k <= ks[j]) at every save
    best = {N: dict(D=-np.inf) for N in Ns}
    for j, hb in enumerate(hh):
        mc2 = -V2c[j] - ev['xi'] * Rbar - m2xbar - mp2
        ok = np.where((mc2 > (cH * H) ** 2) & (tau >= tau_min) & (tau < tau_th))[0]
        if not len(ok):
            continue
        i = ok[0]
        kc = kc_f * a[i] * math.sqrt(mc2[i])
        jj = np.searchsorted(ks, kc) - 1
        if jj < 0:
            continue
        s2 = float(cum[jj, i])
        x = (eG * hb / H_GeV) ** 2 / (2 * s2)
        base = (math.log1p(x) - x) / LN10 + 3 * math.log10(kc)
        for N in Ns:
            D = base + 3 * N / LN10
            if D > best[N]['D']:
                best[N] = dict(D=D, hbar_over_hmax=hb / smi.h_max, tau_c=float(tau[i]), k_c=kc, sigma2_c=s2, x=x)
    return best


def crit(Mts, Ds):
    for i in range(len(Mts) - 1):
        if Ds[i] < 0 <= Ds[i + 1]:
            return Mts[i] + (0 - Ds[i]) * (Mts[i + 1] - Mts[i]) / (Ds[i + 1] - Ds[i])
    return None


def jdump(o, p):
    json.dump(o, open(p, 'w'), indent=1, default=lambda z: float(z) if isinstance(z, (np.floating, np.integer)) else str(z))
