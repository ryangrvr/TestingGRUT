#!/usr/bin/env python3
"""kin_rates (Claude, own code): 2->2 cross sections, 1<->2 rates and their Maxwell-Boltzmann thermal averages,
with the real-intermediate-state (RIS) subtraction for the poles that lie in the physical region.

Amplitudes (tree, M = -prod V prod 1/(p^2-m^2); amp.py checks them against the recursion to 1e-15):
  X Xb -> X Xb : -4 lX - lam3^2 [1/(s-ms^2) + 1/(t-ms^2)]            s-channel pole at s = ms^2 (physical)
  X X  -> X X  : -4 lX - lam3^2 [1/(t-ms^2) + 1/(u-ms^2)]            no pole
  X s  -> X s  : -2 g^2 - lam3^2 [1/(s-mX^2) + 1/(u-mX^2)] - lam3 lam_s3/(t-ms^2)
                                                                      u-channel X pole inside the angular range for
                                                                      sqrt s between ~3.6 and 5 m_X (Coleman-Norton
                                                                      sequence s -> X Xb, Xb X -> s)
  s s  -> s s  : -lam_s4 - lam_s3^2 [1/(s-ms^2) + 1/(t-ms^2) + 1/(u-ms^2)]   no pole (s >= 4 ms^2)
  s s  -> X Xb : -2 g^2 - lam_s3 lam3/(s-ms^2) - lam3^2 [1/(t-mX^2) + 1/(u-mX^2)]   no pole (shown in NOTES)
  X Xb -> s s  : same |M|^2 with crossed kinematics
  s   -> X Xb  : M = -lam3,  Gamma = lam3^2 beta/(16 pi ms)

RIS subtraction.  Write M = A + B D with D = 1/(x + i eps), x = p^2 - m^2, eps = m Gamma, for the pole line.  Then
|M|^2 = A^2 + 2AB Re D + B^2 |D|^2 and the on-shell part B^2 (pi/eps) delta(x) is the sequential process already in
the network (inverse decay x decay, or decay x inverse decay).  The subtracted kernel replaces |D|^2 by
Re(D^2) = (x^2-eps^2)/(x^2+eps^2)^2, whose integral against a smooth h is the Hadamard finite part of int h/x^2 and
is regulator independent as eps -> 0; Re D gives the principal value.  Both are integrated with an explicit
subtraction of h(0) (+ h'(0) x) and analytic remainders (fp_quad, pv_quad), so the tiny physical width
(Gamma_s/ms = 7.7e-5) is handled exactly.

Cross sections: sigma = (1/S_out) (k_f/k_i)/(32 pi s) int_{-1}^{1} dcos |M|^2;  sigma_T uses the extra (1-cos) weight.
Thermal averages (MB, one degree per species, g = 1):
  gamma(T) = (1/S_in) (T/(32 pi^4)) int ds lambda(s,m1^2,m2^2) sigma(s) K_1(sqrt s/T)/sqrt s     [reactions/volume/time]
  <sigma v> = gamma/(n1 n2),  n_i = m_i^2 T K_2(m_i/T)/(2 pi^2).
Parameters can be replaced by in-medium values (masses, lam3 = 2 g^2 v_T, lam_s3 = 6 lP v_T, Gamma_s).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import kv, kve
import model as M
from amp import kallen

_GL = {}


def gl(n):
    if n not in _GL:
        _GL[n] = leggauss(n)
    return _GL[n]


class Params:
    def __init__(self, mX=M.MX, ms=M.MS, lam3=M.LAM3, lam_s3=M.LAM_S3, lam_s4=M.LAM_S4, g2=M.GX**2, lX=M.LX,
                 Gamma_s=None, Gamma_X_reg=1e-7):
        self.mX, self.ms, self.lam3, self.lam_s3, self.lam_s4, self.g2, self.lX = mX, ms, lam3, lam_s3, lam_s4, g2, lX
        self.mX2, self.ms2 = mX * mX, ms * ms
        self.beta = np.sqrt(max(1 - 4 * self.mX2 / self.ms2, 0))
        self.Gamma_s = Gamma_s if Gamma_s is not None else lam3**2 * self.beta / (16 * np.pi * ms)
        self.eps_s = ms * self.Gamma_s           # regulator of the s propagator (physical width)
        self.eps_X = mX * Gamma_X_reg * mX      # regulator of the X propagator: Gamma_X_reg in units of mX
        self.n = 96

    def medium(self, T, thermal_module=None):
        """Parameters at temperature T from thermal.thermal_masses (in-medium masses and v_T couplings; the width
        is the vacuum formula with in-medium masses; Bose enhancement is applied separately)."""
        import thermal as TH
        d = TH.thermal_masses(T)
        return Params(mX=d['mX'], ms=d['ms'], lam3=d['lam3'], lam_s3=d['lam_s3'])


CHANNELS = {
    # name: (m1, m2, m3, m4 keys), S_in, S_out
    'XXb_XXb': (('mX', 'mX', 'mX', 'mX'), 1, 1),
    'XX_XX': (('mX', 'mX', 'mX', 'mX'), 2, 2),
    'Xs_Xs': (('mX', 'ms', 'mX', 'ms'), 1, 1),
    'ss_ss': (('ms', 'ms', 'ms', 'ms'), 2, 2),
    'ss_XXb': (('ms', 'ms', 'mX', 'mX'), 2, 1),
    'XXb_ss': (('mX', 'mX', 'ms', 'ms'), 1, 2),
}


def masses_of(ch, P):
    keys, Sin, Sout = CHANNELS[ch]
    return [getattr(P, k) for k in keys], Sin, Sout


# ---------------------------------------------------------------- finite-part / principal-value quadratures
def fp_quad(h, a, b, eps, n=96, hprime=None):
    """int_a^b h(x) (x^2-eps^2)/(x^2+eps^2)^2 dx  (= Re D^2 weight), h vectorised.  Subtraction of h0 + h1 x when the
    pole is inside or near the interval; analytic remainders.  Exact for finite eps up to quadrature error."""
    x0, w0 = gl(n)
    L = b - a
    inside = (a - 0.25 * L) < 0 < (b + 0.25 * L)
    if not inside:
        x = a + 0.5 * L * (x0 + 1); w = 0.5 * L * w0
        return np.sum(w * h(x) * (x * x - eps * eps) / (x * x + eps * eps)**2)
    h0 = float(h(np.array([0.0]))[0])
    if hprime is None:
        d = 1e-4 * L
        h1 = float((h(np.array([d]))[0] - h(np.array([-d]))[0]) / (2 * d))
    else:
        h1 = hprime
    # smooth remainder, resolved on two sub-panels on each side of the pole
    tot = 0.0
    for (lo, hi) in ((a, 0.0), (0.0, b)):
        if hi <= lo:
            continue
        # geometric sub-panels toward the pole to resolve the eps-scale structure of the remainder
        edges = [lo, hi]
        tot += _panel(lambda x: (h(x) - h0 - h1 * x) * (x * x - eps * eps) / (x * x + eps * eps)**2, lo, hi, n)
    # analytic: int (x^2-eps^2)/(x^2+eps^2)^2 = -x/(x^2+eps^2);  int x (x^2-eps^2)/(x^2+eps^2)^2 = 1/2 ln(x^2+eps^2) + eps^2/(x^2+eps^2)
    F0 = lambda x: -x / (x * x + eps * eps)
    F1 = lambda x: 0.5 * np.log(x * x + eps * eps) + eps * eps / (x * x + eps * eps)
    return tot + h0 * (F0(b) - F0(a)) + h1 * (F1(b) - F1(a))


def pv_quad(h, a, b, eps, n=96):
    """int_a^b h(x) x/(x^2+eps^2) dx (= Re D weight)."""
    x0, w0 = gl(n)
    L = b - a
    inside = (a - 0.25 * L) < 0 < (b + 0.25 * L)
    if not inside:
        x = a + 0.5 * L * (x0 + 1); w = 0.5 * L * w0
        return np.sum(w * h(x) * x / (x * x + eps * eps))
    h0 = float(h(np.array([0.0]))[0])
    tot = 0.0
    for (lo, hi) in ((a, 0.0), (0.0, b)):
        if hi <= lo:
            continue
        tot += _panel(lambda x: (h(x) - h0) * x / (x * x + eps * eps), lo, hi, n)
    return tot + h0 * 0.5 * np.log((b * b + eps * eps) / (a * a + eps * eps))


def _panel(f, lo, hi, n):
    """GL on [lo,hi] split geometrically toward 0 (the pole) when 0 is an endpoint."""
    x0, w0 = gl(n)
    if lo == 0.0 or hi == 0.0:
        # panels [r^k, r^{k-1}] scaled to the interval
        L = hi - lo
        s = 1.0 if lo == 0.0 else -1.0
        edges = [L * 0.5**k for k in range(0, 40)] + [0.0]
        tot = 0.0
        for k in range(len(edges) - 1):
            e1, e2 = edges[k + 1], edges[k]
            if e2 - e1 < 1e-14 * L:
                break
            x = e1 + 0.5 * (e2 - e1) * (x0 + 1); w = 0.5 * (e2 - e1) * w0
            xx = (lo + x) if s > 0 else (hi - x)
            tot += np.sum(w * f(xx))
        return tot
    x = lo + 0.5 * (hi - lo) * (x0 + 1); w = 0.5 * (hi - lo) * w0
    return np.sum(w * f(x))


# ---------------------------------------------------------------- amplitudes split into (A, B, pole variable)
def kin(ch, s, P):
    m, Sin, Sout = masses_of(ch, P)
    m1, m2, m3, m4 = m
    sq = np.sqrt(s)
    ki = np.sqrt(max(kallen(s, m1 * m1, m2 * m2), 0)) / (2 * sq)
    kf = np.sqrt(max(kallen(s, m3 * m3, m4 * m4), 0)) / (2 * sq)
    E1 = (s + m1 * m1 - m2 * m2) / (2 * sq); E3 = (s + m3 * m3 - m4 * m4) / (2 * sq)
    E2 = sq - E1; E4 = sq - E3
    return dict(m=m, Sin=Sin, Sout=Sout, ki=ki, kf=kf, E1=E1, E2=E2, E3=E3, E4=E4, s=s)


def tu(K, c):
    """t = (p1-k1)^2, u = (p1-k2)^2 for cos theta = c."""
    m1, m2, m3, m4 = K['m']
    t = m1 * m1 + m3 * m3 - 2 * (K['E1'] * K['E3'] - K['ki'] * K['kf'] * c)
    u = m1 * m1 + m4 * m4 - 2 * (K['E1'] * K['E4'] + K['ki'] * K['kf'] * c)
    return t, u


def amp_parts(ch, K, c, P):
    """Return (A(c), B, x(c), eps) with M = A + B/(x + i eps); B = 0 when the channel has no physical pole."""
    s = K['s']; t, u = tu(K, c)
    if ch == 'XXb_XXb':
        A = -4 * P.lX - P.lam3**2 / (t - P.ms2)
        return A, -P.lam3**2, (s - P.ms2) * np.ones_like(c), P.eps_s
    if ch == 'XX_XX':
        return -4 * P.lX - P.lam3**2 * (1 / (t - P.ms2) + 1 / (u - P.ms2)), 0.0, None, 0.0
    if ch == 'Xs_Xs':
        A = -2 * P.g2 - P.lam3**2 / (s - P.mX2) - P.lam3 * P.lam_s3 / (t - P.ms2)
        return A, -P.lam3**2, u - P.mX2, P.eps_X
    if ch == 'ss_ss':
        return -P.lam_s4 - P.lam_s3**2 * (1 / (s - P.ms2) + 1 / (t - P.ms2) + 1 / (u - P.ms2)), 0.0, None, 0.0
    if ch in ('ss_XXb', 'XXb_ss'):
        # for XXb_ss, t = (pX - ks)^2 plays the role of the X-exchange invariant: same expression by crossing
        return -2 * P.g2 - P.lam_s3 * P.lam3 / (s - P.ms2) - P.lam3**2 * (1 / (t - P.mX2) + 1 / (u - P.mX2)), 0.0, None, 0.0
    raise KeyError(ch)


def angular(ch, s, P, weight=None, mode='sub', n=None):
    """int_{-1}^{1} dc |M|^2 w(c) with the pole prescription.  mode: 'sub' (RIS-subtracted), 'nores' (A^2 only),
    'int' (2AB ReD only), 'res' (B^2 Re D^2 only), 'full' (finite-width |M|^2, no subtraction).
    Returns the integral (float)."""
    n = n or P.n
    K = kin(ch, s, P)
    if K['ki'] <= 0 or K['kf'] <= 0:
        return 0.0
    w = (lambda c: np.ones_like(c)) if weight is None else weight
    x0, w0 = gl(n)
    A, B, x, eps = amp_parts(ch, K, x0, P)
    if B == 0.0 or ch == 'XXb_XXb':
        # no angular pole: the s-pole of XXb_XXb is handled in the s-integration (components returned by sigma_parts)
        if mode == 'nores' or B == 0.0:
            return float(np.sum(w0 * A * A * w(x0)))
        D = 1 / (x + 1j * eps)
        if mode == 'full':
            return float(np.sum(w0 * (A * A + 2 * A * B * D.real + B * B * np.abs(D)**2) * w(x0)))
        if mode == 'sub':
            return float(np.sum(w0 * (A * A + 2 * A * B * D.real + B * B * (D * D).real) * w(x0)))
        if mode == 'int':
            return float(np.sum(w0 * 2 * A * B * D.real * w(x0)))
        if mode == 'res':
            return float(np.sum(w0 * B * B * (D * D).real * w(x0)))
    # Xs_Xs: u-pole, x linear in c with slope -2 ki kf: c = (x - x(0)) / slope
    xa = float(amp_parts(ch, K, np.array([-1.0]), P)[2][0]); xb = float(amp_parts(ch, K, np.array([1.0]), P)[2][0])
    lo, hi = min(xa, xb), max(xa, xb)
    slope = (xb - xa) / 2.0
    cof = lambda xx: (xx - xa) / slope - 1.0
    hA2 = lambda xx: (lambda c: amp_parts(ch, K, c, P)[0]**2 * w(c))(cof(xx))
    hAB = lambda xx: (lambda c: 2 * amp_parts(ch, K, c, P)[0] * B * w(c))(cof(xx))
    hB2 = lambda xx: B * B * w(cof(xx))
    jac = 1.0 / abs(slope)
    if mode == 'nores':
        return float(np.sum(w0 * A * A * w(x0)))
    if mode == 'full':
        D = 1 / (x + 1j * eps)
        # the full integrand is finite for eps>0; resolve the peak by the same panel machinery via fp + delta term
        res = fp_quad(hB2, lo, hi, eps, n) + _delta_part(hB2, lo, hi, eps)
        return jac * (_panel_plain(hA2, lo, hi, n) + pv_quad(hAB, lo, hi, eps, n) + res)
    parts = 0.0
    if mode in ('sub', 'nores_int_res', 'int'):
        parts += pv_quad(hAB, lo, hi, eps, n) if mode != 'nores' else 0.0
    if mode == 'int':
        return jac * parts
    if mode == 'res':
        return jac * fp_quad(hB2, lo, hi, eps, n)
    return jac * (_panel_plain(hA2, lo, hi, n) + parts + fp_quad(hB2, lo, hi, eps, n))


def _panel_plain(h, lo, hi, n):
    x0, w0 = gl(n)
    x = lo + 0.5 * (hi - lo) * (x0 + 1); w = 0.5 * (hi - lo) * w0
    return float(np.sum(w * h(x)))


def _delta_part(h, lo, hi, eps):
    """int_lo^hi h(x) [ |D|^2 - Re D^2 ] dx = int h 2 eps^2/(x^2+eps^2)^2 dx  (the on-shell piece, -> pi/eps h(0))."""
    if not (lo < 0 < hi):
        x0, w0 = gl(96); x = lo + 0.5 * (hi - lo) * (x0 + 1); w = 0.5 * (hi - lo) * w0
        return float(np.sum(w * h(x) * 2 * eps * eps / (x * x + eps * eps)**2))
    h0 = float(h(np.array([0.0]))[0])
    F = lambda x: x / (x * x + eps * eps) + np.arctan(x / eps) / eps   # int 2 eps^2/(x^2+eps^2)^2 dx
    rem = _panel(lambda x: (h(x) - h0) * 2 * eps * eps / (x * x + eps * eps)**2, lo, 0.0, 96) + \
        _panel(lambda x: (h(x) - h0) * 2 * eps * eps / (x * x + eps * eps)**2, 0.0, hi, 96)
    return rem + h0 * (F(hi) - F(lo))


def sigma(ch, s, P, mode='sub', transport=False):
    """Cross section (GeV^-2) at Mandelstam s; for XXb_XXb the s-pole is kept as components (see sigma_parts) and
    'sub'/'full' here use the finite physical width Gamma_s pointwise."""
    K = kin(ch, s, P)
    if K['ki'] <= 0 or K['kf'] <= 0:
        return 0.0
    w = (lambda c: 1 - c) if transport else None
    I = angular(ch, s, P, weight=w, mode=mode)
    return I * (K['kf'] / K['ki']) / (32 * np.pi * s) / K['Sout']


def sigma_parts(s, P, transport=False):
    """XXb_XXb: (sigma0, sigma1, sigma2) with sigma = sigma0 + sigma1 Re D(s) + sigma2 Re D(s)^2 (sub) or |D|^2."""
    ch = 'XXb_XXb'
    K = kin(ch, s, P)
    if K['ki'] <= 0:
        return 0.0, 0.0, 0.0
    w = (lambda c: 1 - c) if transport else (lambda c: np.ones_like(c))
    x0, w0 = gl(P.n)
    A, B, x, eps = amp_parts(ch, K, x0, P)
    pref = (K['kf'] / K['ki']) / (32 * np.pi * s) / K['Sout']
    return (pref * float(np.sum(w0 * A * A * w(x0))), pref * float(np.sum(w0 * 2 * A * B * w(x0))),
            pref * float(np.sum(w0 * B * B * w(x0))))


# ---------------------------------------------------------------- MB thermodynamics and thermal averages
def n_eq(m, T, mu=0.0):
    return m * m * T * kve(2, m / T) * np.exp(-(m - mu) / T) / (2 * np.pi**2)


def thermal_rate(ch, T, P, mode='sub', transport=False, n_s=160, smax_T=40.0):
    """gamma(T) = (1/S_in)(T/32 pi^4) int ds lambda sigma K1(sqrt s/T)/sqrt s, reactions per volume per time, MB.
    The XXb_XXb s-pole is integrated with fp_quad/pv_quad; other channels on a log grid in sqrt s."""
    m, Sin, Sout = masses_of(ch, P)
    smin = max(m[0] + m[1], m[2] + m[3])**2
    sqmax = np.sqrt(smin) + smax_T * T
    pref = T / (32 * np.pi**4) / Sin

    def weight(s):
        return kallen(s, m[0] * m[0], m[1] * m[1]) * kve(1, np.sqrt(s) / T) * np.exp(-np.sqrt(s) / T) / np.sqrt(s)
    if ch != 'XXb_XXb' or mode == 'nores':
        # log grid in sqrt s with GL panels
        edges = np.geomspace(np.sqrt(smin) * (1 + 1e-12), sqmax, 25)
        x0, w0 = gl(24)
        tot = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            # map to log
            la, lb = np.log(a), np.log(b)
            l = la + 0.5 * (lb - la) * (x0 + 1); w = 0.5 * (lb - la) * w0
            sq = np.exp(l)
            vals = np.array([sigma(ch, q * q, P, mode=mode, transport=transport) * weight(q * q) for q in sq])
            tot += np.sum(w * vals * 2 * sq * sq)      # ds = 2 sq d(sq) = 2 sq^2 d(ln sq)
        return pref * tot
    # XXb_XXb with the s-pole: components sigma0 + sigma1 ReD + sigma2 ReD^2 (or |D|^2 in 'full')
    eps = P.eps_s
    xa, xb = smin - P.ms2, sqmax**2 - P.ms2
    h0 = lambda x: np.array([sigma_parts(xx + P.ms2, P, transport)[0] * weight(xx + P.ms2) for xx in np.atleast_1d(x)])
    h1 = lambda x: np.array([sigma_parts(xx + P.ms2, P, transport)[1] * weight(xx + P.ms2) for xx in np.atleast_1d(x)])
    h2 = lambda x: np.array([sigma_parts(xx + P.ms2, P, transport)[2] * weight(xx + P.ms2) for xx in np.atleast_1d(x)])
    # split the x range into panels: [xa, -W], [-W, W] (pole), [W, xb] with W = 0.5 ms^2, log-spaced outside
    W = min(0.5 * P.ms2, 0.5 * (P.ms2 - smin))
    tot = 0.0
    if mode == 'nores':
        return pref * (_log_panels(h0, xa, xb))
    # non-resonant
    tot += _log_panels(h0, xa, xb)
    # interference and resonant, with the pole window treated by pv/fp
    tot += pv_quad(h1, -W, W, eps, 96) + _log_panels(lambda x: h1(x) * x / (x * x + eps * eps), xa, -W) \
        + _log_panels(lambda x: h1(x) * x / (x * x + eps * eps), W, xb)
    if mode == 'sub':
        tot += fp_quad(h2, -W, W, eps, 96)
        tot += _log_panels(lambda x: h2(x) * (x * x - eps * eps) / (x * x + eps * eps)**2, xa, -W) \
            + _log_panels(lambda x: h2(x) * (x * x - eps * eps) / (x * x + eps * eps)**2, W, xb)
    elif mode == 'full':
        tot += fp_quad(h2, -W, W, eps, 96) + _delta_part(h2, -W, W, eps)
        tot += _log_panels(lambda x: h2(x) / (x * x + eps * eps), xa, -W) + _log_panels(lambda x: h2(x) / (x * x + eps * eps), W, xb)
    elif mode == 'ris':
        return pref * _delta_part(h2, -W, W, eps)
    return pref * tot


def _log_panels(h, a, b, npan=14, n=16):
    """int_a^b h(x) dx with panels geometric in |x| (a, b of the same sign or a<0<b handled by caller)."""
    if b <= a:
        return 0.0
    x0, w0 = gl(n)
    if a > 0:
        edges = np.geomspace(a, b, npan + 1)
    elif b < 0:
        edges = -np.geomspace(-b, -a, npan + 1)[::-1]
    else:
        edges = np.linspace(a, b, npan + 1)
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = lo + 0.5 * (hi - lo) * (x0 + 1); w = 0.5 * (hi - lo) * w0
        tot += np.sum(w * h(x))
    return tot


def sigma_v(ch, T, P, **kw):
    m, Sin, Sout = masses_of(ch, P)
    return thermal_rate(ch, T, P, **kw) * Sin / (n_eq(m[0], T) * n_eq(m[1], T))


def decay_rate_density(T, P):
    """gamma_{s -> X Xb}(T) = n_s^eq Gamma_s <ms/E> = ms^2 T Gamma_s K1(ms/T)/(2 pi^2)  (MB, mu = 0)."""
    return P.ms**2 * T * P.Gamma_s * kve(1, P.ms / T) * np.exp(-P.ms / T) / (2 * np.pi**2)


if __name__ == '__main__':
    P = Params()
    print('Gamma_s =', P.Gamma_s)
    # pole window of Xs_Xs
    for sq in np.array([3.5, 3.6, 4.0, 4.5, 5.0, 5.1]) * M.MX:
        K = kin('Xs_Xs', sq * sq, P)
        xa = amp_parts('Xs_Xs', K, np.array([-1.0]), P)[2][0]; xb = amp_parts('Xs_Xs', K, np.array([1.0]), P)[2][0]
        print(f'sqrt s/mX={sq/M.MX:.2f}  u-mX^2 range/mX^2: {xa/M.MX**2:.3f} .. {xb/M.MX**2:.3f}  pole inside: {xa<0<xb or xb<0<xa}')
    # fp_quad test against analytic: h = 1 + x + x^2 on [-1, 2], eps -> 0: FP int (1+x+x^2)/x^2 = [-1/x + ln|x| + x]
    h = lambda x: 1 + x + x * x
    exact = (-1 / 2 + np.log(2) + 2) - (1 + 0 - 1)   # at b=2 minus at a=-1 : (-1/x + ln|x| + x)
    for eps in (1e-2, 1e-4, 1e-6):
        print('fp_quad', eps, fp_quad(h, -1, 2, eps), 'exact FP', exact)
    print('pv_quad', pv_quad(h, -1, 2, 1e-6), 'exact PV', (np.log(2) + 2 + 2) - (0 - 1 + 0.5))
    # sigma values
    for sq in (2.1, 2.449, 3.0, 5.0, 20.0):
        s = (sq * M.MX)**2
        print(f'sqrt s={sq} mX: sigma XXb sub={sigma("XXb_XXb", s, P):.4e} nores={sigma("XXb_XXb", s, P, "nores"):.4e}  XX={sigma("XX_XX", s, P):.4e}')
    for sq in (3.6, 4.0, 4.5, 5.0, 8.0):
        s = (sq * M.MX)**2
        print(f'sqrt s={sq} mX: Xs sub={sigma("Xs_Xs", s, P):.4e} nores={sigma("Xs_Xs", s, P, "nores"):.4e} res={sigma("Xs_Xs", s, P, "res"):.4e} int={sigma("Xs_Xs", s, P, "int"):.4e}')
    for sq in (4.9, 5.0, 6.0, 10.0):
        s = (sq * M.MX)**2
        print(f'sqrt s={sq} mX: ss_ss={sigma("ss_ss", s, P):.4e} ss_XXb={sigma("ss_XXb", s, P):.4e} XXb_ss={sigma("XXb_ss", s, P):.4e}')
    # detailed balance / RIS normalisation
    for t in (0.3, 1.0, 3.0, 12.0):
        T = t * M.MX
        ris = thermal_rate('XXb_XXb', T, P, mode='ris')
        dec = decay_rate_density(T, P)
        sub = thermal_rate('XXb_XXb', T, P, mode='sub'); nores = thermal_rate('XXb_XXb', T, P, mode='nores')
        print(f'T/m={t}: RIS piece/decay rate = {ris/dec:.6f}; gamma_sub={sub:.4e} gamma_nores={nores:.4e} ratio={sub/nores:.4f}')
        g1 = thermal_rate('ss_XXb', T, P); g2 = thermal_rate('XXb_ss', T, P)
        print(f'      ss->XXb gamma={g1:.4e}  XXb->ss gamma={g2:.4e}  ratio={g1/g2:.6f}')
