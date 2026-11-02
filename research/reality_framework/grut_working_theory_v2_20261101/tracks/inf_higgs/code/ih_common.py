#!/usr/bin/env python3
"""Track inf_higgs (Claude, own code; inflaton-sector pass, issue #13 context, 3 November 2026).

Common module: the inflaton benchmarks of track inf_potential (HT6t, T1, E1) as backgrounds for the Higgs, in H_i units,
from deep in the de Sitter phase (with the benchmark's own H(N), R(N)) through the end of inflation and the
post-inflation oscillation of the benchmark's potential (not a quadratic released at rest); plus the M1r surrogate
(quadratic m_I = 3 H_i released at rest at a = 1, H = 1; the higgs_reheat background) as the regression control.

Units: field in reduced Planck units (M_P = 1); time in 1/H_i; energy densities in M_P^2 H_i^2; R in H_i^2.
H_i := H at the slice S1 (one e-fold before eps_H = 1, inf_potential's primary convention); a(S1) = 1, tau(S1) = 0.
  rho = phidot^2/2 + A v(phi) + rho_r + rho_m,   A = V0/eta^2, eta = H_i/M_P,
  R = 4 A v - phidot^2 + rho_m   (radiation is traceless; Jordan-frame hidden fields neglected: eta^2 V ~ 1e-10 [D]).
Decay: constant width Gamma (H_i units) on the homogeneous inflaton (friction Gamma phidot, source Gamma phidot^2), the
engine's convention; the realistic channels of inf_reheat give Gamma/m_I ~ 1e-6 or less, i.e. Gamma ~ 0 over tau <~ 1e2.
Fragmentation variant [model choice]: at tau_frag the homogeneous inflaton is replaced by a w = 0 fluid with the same
energy (oscillons / fragmented n = 1 attractors, Lozanov-Amin 1608.01213 [R]); used for T1/E1 only.

Reused read-only (imported, bytecode writing disabled): tracks/inf_potential/code/inf_common.py (potentials,
slow-roll background), tracks/higgs_reheat/code/smpot.py (+ sm_rge.py, c1_fit_grid.json; SM running and lambda_eff).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
RES = TRACK / 'results'
TRACKS = TRACK.parent
SP = Path(os.environ.get('IH_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                         'scratchpad/wt2/inf_higgs'))
SP.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(TRACKS / 'inf_potential' / 'code'))
sys.path.insert(0, str(TRACKS / 'higgs_reheat' / 'code'))
import inf_common as IC          # noqa: E402  (read-only reuse)
import smpot as S                # noqa: E402  (read-only reuse)

MPL = 2.435e18
GSTAR = 106.75
PDG24 = dict(Mt=172.57, sMt=0.29, Mh=125.20, sMh=0.11, als=0.1180)   # PDG 2024 [R]
MH_HR = 125.13                                                       # higgs_reheat's input (PDG 2026 there)
BENCH = ('HT6t', 'T1', 'E1')
KIND = {'HT6t': 'HT', 'T1': 'T', 'E1': 'E'}
# fragmentation times (tau after eps_H = 1) from inf_potential NOTES section 3 (linear backreaction) [Q]
TFRAG_AFTER_END = {'T1': 0.43, 'E1': 0.35}


def T_GeV(rho_r, H_GeV, gstar=GSTAR):
    """temperature of radiation rho_r (M_P^2 H_i^2 units)"""
    rho = np.maximum(rho_r, 0.0) * MPL ** 2 * H_GeV ** 2
    return (30 * rho / (math.pi ** 2 * gstar)) ** 0.25


# ------------------------------------------------------------------------------------------- benchmarks
def load_bench(name, iH=0):
    d = json.load(open(TRACKS / 'inf_potential' / 'results' / f'bench_{name}.json'))
    r = d['by_Hi'][iH]
    m = IC.Model(KIND[name], **r['params'])
    eta = r['H_i_GeV'] / MPL
    A = r['V0_MP4'] / eta ** 2
    out = dict(name=name, model=m, A=A, H_GeV=r['H_i_GeV'], rec=r, phimin=m.phi_min(), m_I=r['m_I_over_Hi'],
               Nstar=r['Nstar'], eta=eta)
    return out


class Background:
    """Homogeneous background in H_i units.  kind = 'bench' (needs bench dict) or 'M1r' (quadratic, b, at rest)."""

    def __init__(self, bench=None, gam=0.0, Nb0=5.0, tau_end=80.0, frag=None, b=3.0, rtol=1e-11, atol=1e-14,
                 max_step=None):
        self.bench, self.gam, self.Nb0, self.frag = bench, gam, Nb0, frag
        self.kind = 'M1r' if bench is None else 'bench'
        self.b = b
        if self.kind == 'M1r':
            self.A = 1.0
            self.V = lambda f: 0.5 * b * b * f * f
            self.dV = lambda f: b * b * f
            self.phimin = 0.0
            self.y0 = np.array([0.0, math.sqrt(6.0) / b, 0.0, 0.0, 0.0])
            self.tau0 = 0.0
            self.m_I = b
        else:
            m, A = bench['model'], bench['A']
            self.A = A
            self.V = lambda f, m=m, A=A: A * m.v(f)
            self.dV = lambda f, m=m, A=A: A * m.dv(f)
            self.phimin = bench['phimin']
            self.m_I = bench['m_I']
            bg = IC.background(m, Nwant=max(75.0, Nb0 + 10))
            self.bg_slow = bg
            f, fN = IC.bg_at(bg, Nb0)
            eps = 0.5 * fN * fN
            H = math.sqrt(A * m.v(f) / (3 - eps))
            y0 = np.array([-(Nb0 - 1.0), f, H * fN, 0.0, 0.0])
            # time from Nb0 to S1 (lna = 0) with the background alone
            ev = lambda t, y: y[0]
            ev.terminal = True; ev.direction = 1
            s = solve_ivp(self._rhs, (0, 50), y0, method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
            self.tau0 = -float(s.t_events[0][0])
            self.y0 = y0
        self.tau_end = tau_end
        self.max_step = max_step if max_step else (0.02 if self.kind == 'bench' else 0.02)
        self._integrate(rtol, atol)

    # ------------------------------------------------------------------ equations
    def H_of(self, y):
        lna, f, fd, rr, rm = y[:5]
        return math.sqrt(max((0.5 * fd * fd + self.V(f) + rr + rm) / 3.0, 1e-300))

    def R_of(self, y):
        f, fd, rm = y[1], y[2], y[4]
        return 4 * self.V(f) - fd * fd + rm

    def _rhs(self, t, y):
        lna, f, fd, rr, rm = y[:5]
        H = math.sqrt(max((0.5 * fd * fd + self.V(f) + rr + rm) / 3.0, 1e-300))
        g = self.gam
        return [H, fd, -(3 * H + g) * fd - self.dV(f), g * fd * fd + g * rm - 4 * H * rr, -(3 * H + g) * rm]

    def _integrate(self, rtol, atol):
        """dense background from tau0 to tau_end; optional fragmentation switch at frag (tau)."""
        segs = []
        t0 = self.tau0
        y = self.y0.copy()
        # end-of-inflation event eps_H = 1 (bench) to define tau_endinf
        def epsH(t, y):
            H = self.H_of(y)
            return (y[2] ** 2 + 4 * (y[3]) / 3 + y[4]) / (2 * H * H) - 1.0
        epsH.direction = 1
        tf = self.frag if (self.frag is not None and self.frag < self.tau_end) else None
        t1 = tf if tf is not None else self.tau_end
        s = solve_ivp(self._rhs, (t0, t1), y, method='DOP853', rtol=rtol, atol=atol, dense_output=True,
                      max_step=self.max_step, events=epsH)
        segs.append(s)
        self.tau_endinf = float(s.t_events[0][0]) if (self.kind == 'bench' and len(s.t_events[0])) else 0.0
        if tf is not None:
            y = s.y[:, -1].copy()
            rho_I = 0.5 * y[2] ** 2 + self.V(y[1])
            y[4] += rho_I; y[1] = self.phimin; y[2] = 0.0
            s2 = solve_ivp(self._rhs, (tf, self.tau_end), y, method='DOP853', rtol=rtol, atol=atol,
                           dense_output=True, max_step=0.5)
            segs.append(s2)
        self.segs = segs

    def y_at(self, t):
        t = np.atleast_1d(t)
        out = np.empty((5, len(t)))
        for i, tt in enumerate(t):
            s = self.segs[0] if (len(self.segs) == 1 or tt <= self.segs[0].t[-1]) else self.segs[1]
            out[:, i] = s.sol(min(max(tt, s.t[0]), s.t[-1]))
        return out

    def fields(self, t):
        Y = self.y_at(t)
        lna, f, fd, rr, rm = Y
        V = np.array([self.V(x) for x in f])
        rho = 0.5 * fd ** 2 + V + rr + rm
        H = np.sqrt(rho / 3)
        R = 4 * V - fd ** 2 + rm
        return dict(tau=np.atleast_1d(t), a=np.exp(lna), H=H, R=R, phi=f, phiI=f - self.phimin, phid=fd, rho_r=rr,
                    rho_m=rm, rho=rho, eps=(fd ** 2 + 4 * rr / 3 + rm) / (2 * H * H))

    # ------------------------------------------------------------------ oscillation diagnostics
    def fine(self, dt=None):
        """fine table (for averages and turning points)"""
        if hasattr(self, '_fine'):
            return self._fine
        dt = dt if dt else min(0.002, 0.05 / max(self.m_I, 1.0))
        n = int((self.tau_end - self.tau0) / dt) + 2
        t = np.linspace(self.tau0, self.tau_end, n)
        F = self.fields(t)
        # turning points of phi (sign changes of phidot) after the end of inflation
        fd = F['phid']
        idx = np.where((fd[1:] * fd[:-1] < 0) & (t[1:] > self.tau_endinf - 0.05))[0]
        if self.frag is not None:
            idx = idx[t[idx] < self.frag]
        F['turn_idx'] = idx
        F['turn_tau'] = t[idx]
        F['turn_phiI'] = F['phiI'][idx]
        self._fine = F
        return F

    def Rbar(self, taus, mode='period', P=None):
        """running mean of R over one oscillation period (mode 'period': local phi period from turning points;
        mode 'fixed': window P).  Before the first turning point (inflation) Rbar = R."""
        F = self.fine()
        t, R = F['tau'], F['R']
        cR = np.concatenate([[0.0], np.cumsum(0.5 * (R[1:] + R[:-1]) * np.diff(t))])
        ci = interp1d(t, cR, bounds_error=False, fill_value=(cR[0], cR[-1]))
        taus = np.atleast_1d(taus)
        if mode == 'fixed':
            w = np.full_like(taus, P)
        else:
            tt = F['turn_tau']
            if len(tt) >= 3:
                per = tt[2:] - tt[:-2]               # full period between alternate turning points
                tm = tt[1:-1]
                w = np.interp(taus, tm, per, left=per[0], right=per[-1])
                # inflation (before eps_H = 1): no average; from eps_H = 1 on, average over the (first) local period
                w = np.where(taus < self.tau_endinf, 0.0, w)
            else:
                w = np.zeros_like(taus)
            if self.frag is not None:
                w = np.where(taus > self.frag, 0.0, w)  # fluid: R already smooth
        out = np.empty_like(taus, dtype=float)
        Rt = interp1d(t, R, bounds_error=False, fill_value=(R[0], R[-1]))
        for i, (x, ww) in enumerate(zip(taus, w)):
            if ww <= 0:
                out[i] = float(Rt(x))
            else:
                lo, hi = max(x - ww / 2, t[0]), min(x + ww / 2, t[-1])
                out[i] = float((ci(hi) - ci(lo)) / (hi - lo))
        return out

    def windows(self, nmax=12, xi=0.0877):
        """R < 0 windows after the end of inflation: start, end, R_min/H^2, duration in H^-1, tachyonic exponent
        int sqrt(max(-xi R, 0)) dt per window."""
        F = self.fine()
        t, R, H = F['tau'], F['R'], F['H']
        neg = R < 0
        out = []
        i = 0
        n = len(t)
        start = np.searchsorted(t, self.tau_endinf - 0.5)
        i = start
        while i < n and len(out) < nmax:
            if neg[i]:
                j = i
                while j < n and neg[j]:
                    j += 1
                seg = slice(i, j)
                k = i + int(np.argmin(R[seg]))
                expo = float(np.trapezoid(np.sqrt(np.maximum(-xi * R[seg], 0)), t[seg]))
                out.append(dict(tau_start=float(t[i]), tau_stop=float(t[min(j, n - 1)]), tau_Rmin=float(t[k]),
                                R_min=float(R[k]), H=float(H[k]), R_min_over_H2=float(R[k] / H[k] ** 2),
                                duration_H=float((t[min(j, n - 1)] - t[i]) * H[k]), exponent_xi=expo))
                i = j
            else:
                i += 1
        return out


# ------------------------------------------------------------------------------------------- SM tables (cached)
class SMTab:
    """Flat 2-loop lambda_eff(h0) (incl. e^{4 Gamma}) on a log grid, Gamma(mu), h_max, kappa_T, xi running; cached."""

    def __init__(self, Mt, Mh=PDG24['Mh'], als=PDG24['als']):
        self.Mt, self.Mh, self.als = Mt, Mh, als
        key = f'sm_{Mt:.4f}_{Mh:.3f}_{als:.4f}.json'
        path = SP / 'smcache' / key
        path.parent.mkdir(exist_ok=True)
        if path.exists():
            d = json.load(open(path))
        else:
            sm = S.SM(Mt, Mh=Mh, als=als)
            hs = np.exp(np.linspace(math.log(1e8), math.log(1e18), 400))
            le = sm.lam_eff_flat(hs)
            mus = np.exp(np.linspace(math.log(1e9), math.log(1e17), 161))
            G = [float(sm.cpl(m)['Gam']) for m in mus]
            Ts = np.exp(np.linspace(math.log(1e10), math.log(1e16), 61))
            kT = []
            for T in Ts:
                c = sm.cpl(math.pi * T)
                kT.append(float(3 * c['G2'] / 16 + 0.6 * c['G1'] / 16 + c['T'] / 4 + c['lam'] / 2))
            fl = sm.flat_scales()
            Ixi = [float(sm.cpl(m)['Ixi']) for m in mus]
            Ixi_t = float(sm.cpl(Mt)['Ixi'])
            d = dict(hs=hs.tolist(), le=le.tolist(), mus=mus.tolist(), G=G, Ts=Ts.tolist(), kT=kT, flat=fl,
                     Ixi=Ixi, Ixi_Mt=Ixi_t)
            json.dump(d, open(path, 'w'))
        self.d = d
        self.flat = d['flat']
        self._le = interp1d(np.log(d['hs']), d['le'], bounds_error=False, fill_value=(d['le'][0], -0.01))
        self._G = interp1d(np.log(d['mus']), d['G'], bounds_error=False, fill_value=(d['G'][0], d['G'][-1]))
        self._kT = interp1d(np.log(d['Ts']), d['kT'], bounds_error=False, fill_value=(d['kT'][0], d['kT'][-1]))
        self._Ixi = interp1d(np.log(d['mus']), d['Ixi'], bounds_error=False, fill_value=(d['Ixi'][0], d['Ixi'][-1]))

    def lam(self, h):
        return self._le(np.log(np.maximum(h, 1.0)))

    def G(self, mu):
        return float(self._G(math.log(mu)))

    def kT(self, T):
        return float(self._kT(math.log(max(T, 1e10))))

    def xi(self, mu, xi_EW):
        """xi(mu) from xi(M_t) = xi_EW (one-loop running along the 3-loop trajectory)"""
        return 1 / 6 + (xi_EW - 1 / 6) * math.exp(float(self._Ixi(math.log(mu))) - self.d['Ixi_Mt'])

    def xi_EW(self, mu, xi_mu):
        return 1 / 6 + (xi_mu - 1 / 6) * math.exp(-(float(self._Ixi(math.log(mu))) - self.d['Ixi_Mt']))


def jdump(obj, path):
    json.dump(obj, open(path, 'w'), indent=1, default=float)
