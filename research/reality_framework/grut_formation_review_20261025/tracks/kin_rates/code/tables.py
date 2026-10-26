#!/usr/bin/env python3
"""kin_rates: interpolation of the Monte Carlo kernel tables and their MB thermal averages (own code).

n->2 kernel Gt(s) (dimensionless) and 2->n cross section sigma(s) (GeV^-2) are tabulated on sqrt s = sqrt s_th (1+delta).
Interpolation: the regulated 'full' column is positive and smooth -> log-log in delta; the subtracted kernel is
carried as the ratio r = G_sub/G_full (linear in ln delta, bounded).  Below the first point: power law of the
full column between the first two points (threshold behaviour) and the first ratio; above the last point: last
value (the tables are flat at high energy) and last ratio.
Thermal averages (MB, mu = 0):
    n->2 : gamma(T) = (1/(S_in S_out)) (T/8 pi^3) int ds sqrt s K1(sqrt s/T) Gt(s)
    2->n : gamma(T) = (1/S_in) (T/32 pi^4) int ds lambda(s,m1^2,m2^2) sigma(s) K1(sqrt s/T)/sqrt s
(both derived in NOTES; identical for the same |M|^2, which is the table-level detailed-balance check).
"""
import numpy as np
from scipy.special import kve
import model as M


def kallen(a, b, c):
    return a * a + b * b + c * c - 2 * a * b - 2 * a * c - 2 * b * c


class Table:
    def __init__(self, chan, kind, freg=0.02):
        """chan: channel dict from kernels_n2.json / xsec_2n.json; kind 'n2' or '2n'."""
        self.kind = kind
        self.Sin, self.Sout = chan['Sin'], chan['Sout']
        self.m_in, self.m_out = chan['m_in'], chan['m_out']
        self.sq_th = chan['sqrt_s_th']
        rows = [r for r in chan['rows'] if abs(r['freg'] - freg) < 1e-12]
        rows = sorted(rows, key=lambda r: r['sqrt_s_over_mX'])
        self.sq = np.array([r['sqrt_s_over_mX'] * M.MX for r in rows])
        key, keyf, keye = ('G', 'G_full', 'err') if kind == 'n2' else ('sigma', 'sigma_full', 'err')
        self.val = np.array([r[key] for r in rows]); self.full = np.array([r[keyf] for r in rows])
        self.err = np.array([r[keye] for r in rows])
        self.delta = self.sq / self.sq_th - 1
        self.ld = np.log(self.delta)
        self.lfull = np.log(np.maximum(self.full, 1e-300))
        self.ratio = self.val / np.where(self.full > 0, self.full, 1.0)
        self.rerr = self.err / np.where(self.full > 0, self.full, 1.0)
        self.k_th = (self.lfull[1] - self.lfull[0]) / (self.ld[1] - self.ld[0])     # threshold power law

    def full_of(self, sq):
        d = np.asarray(sq) / self.sq_th - 1
        out = np.zeros_like(d, dtype=float)
        ok = d > 0
        ld = np.log(np.where(ok, d, 1.0))
        lf = np.interp(ld, self.ld, self.lfull)
        lo = ld < self.ld[0]
        lf = np.where(lo, self.lfull[0] + self.k_th * (ld - self.ld[0]), lf)
        out[ok] = np.exp(lf[ok])
        return out

    def ratio_of(self, sq, which='ratio'):
        d = np.asarray(sq) / self.sq_th - 1
        ld = np.log(np.where(d > 0, d, 1e-300))
        arr = self.ratio if which == 'ratio' else self.rerr
        return np.interp(ld, self.ld, arr)

    def value(self, sq, mode='sub', shift=0.0):
        """Kernel at sqrt s = sq (array).  mode 'sub' (subtracted), 'full' (regulated full); shift = +-1 moves the
        subtracted kernel by its MC error (error propagation).  Where the two bracketing table points are both
        positive the subtracted kernel is interpolated log-log directly (smooth power-law regions, e.g. below the
        chain threshold); elsewhere as ratio x full."""
        sq = np.asarray(sq, dtype=float)
        f = self.full_of(sq)
        if mode == 'full':
            return f
        r = self.ratio_of(sq) + shift * self.ratio_of(sq, 'rerr')
        out = r * f
        d = sq / self.sq_th - 1
        ok = d > 0
        ld = np.log(np.where(ok, d, 1.0))
        idx = np.clip(np.searchsorted(self.ld, ld) - 1, 0, len(self.ld) - 2)
        v0 = self.val[idx] + shift * self.err[idx]; v1 = self.val[idx + 1] + shift * self.err[idx + 1]
        inside = (ld >= self.ld[0]) & (ld <= self.ld[-1]) & (v0 > 0) & (v1 > 0)
        if inside.any():
            w = (ld - self.ld[idx]) / (self.ld[idx + 1] - self.ld[idx])
            lv = np.log(np.where(inside, v0, 1.0)) * (1 - w) + np.log(np.where(inside, v1, 1.0)) * w
            out = np.where(inside, np.exp(lv), out)
        return np.where(ok, out, 0.0)

    def thermal(self, T, mode='sub', shift=0.0, npts=3000, top=None):
        sq_max = self.sq_th + 60 * T if top is None else top
        sq = np.geomspace(self.sq_th * (1 + 1e-7), max(sq_max, self.sq_th * 1.001), npts)
        s = sq * sq
        val = self.value(sq, mode, shift)
        w = kve(1, sq / T) * np.exp(-sq / T)
        if self.kind == 'n2':
            integrand = 2 * sq * sq * w * val                 # ds = 2 sq d(sq)
            pref = T / (8 * np.pi**3) / (self.Sin * self.Sout)
        else:
            lam = kallen(s, self.m_in[0]**2, self.m_in[1]**2)
            integrand = 2 * sq * lam * val * w / sq           # ds lambda sigma K1/sqrt s
            pref = T / (32 * np.pi**4) / self.Sin
        return pref * np.trapezoid(integrand, sq)

    def chain_threshold(self, chan, freg=0.02):
        """First tabulated sqrt s/m_X where any line was flagged near shell (sequential chain open)."""
        rows = sorted([r for r in chan['rows'] if abs(r['freg'] - freg) < 1e-12], key=lambda r: r['sqrt_s_over_mX'])
        for r in rows:
            if r.get('lines_flagged', 0) > 0 and (r.get('sub', 0) != 0 or r.get('add', 0) != 0):
                return r['sqrt_s_over_mX']
        return None


def n_eq(m, T):
    return m * m * T * kve(2, m / T) * np.exp(-m / T) / (2 * np.pi**2)
