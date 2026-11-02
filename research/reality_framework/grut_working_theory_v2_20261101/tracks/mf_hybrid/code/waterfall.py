"""Homogeneous two-field waterfall with the spectator s and linear chi modes [C given the model choices].

From ~6 e-folds before sigma_c (slow-roll attractor) through the biased waterfall and the first coherent oscillations.
Cosmic time in H_ref units. Returns a dense record: t, ln a, H, eps_H, R, sigma, chi, s/s_*, energy pieces and the
growth of linear delta chi_k (k in units of a m_chi at the slice), which marks the onset of non-linearity (fragmentation).

The chi one-loop term (MS-bar) is kept as a function of sigma for sigma > sigma_c (the inflationary valley) and dropped
below sigma_c, where it vanishes continuously (M^4 ln M^2 -> 0) and is negligible against the tree terms.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import mfh_common as mc


class Model:
    def __init__(self, p, inf, c_bias, smech='a', g_s=None, xi_s=-0.75, m_s=3.0):
        self.p, self.inf, self.cb = p, inf, c_bias
        self.smech, self.xi_s, self.m_s = smech, xi_s, m_s
        self.g_s = g_s if g_s is not None else m_s / p['v']
        # vacuum: chi ~ v, sigma ~ sigma_0; constant fixes V(vacuum) = 0 (cosmological-constant tuning, a model choice)
        self.Vc = 0.0
        self.Vc = -self.V(self.vac_sigma(), self.p['v'] + c_bias / (2 * p['lam'] * p['v'] ** 2))

    def vac_sigma(self):
        p = self.p
        return -(self.inf.c1 + self.inf.m2 * p['sigma_c']) / (p['g'] ** 2 * p['v'] ** 2)

    def V(self, s, x):
        p, inf = self.p, self.inf
        d = s - p['sigma_c']
        out = (0.25 * p['lam'] * (x * x - p['v'] ** 2) ** 2 + 0.5 * p['g'] ** 2 * s * s * x * x - self.cb * x
               + inf.c1 * d - 0.5 * inf.m2 * d * d + self.Vc)
        if s > p['sigma_c'] and inf.loop:
            out += float(mc.cw(p['g'] ** 2 * s * s - p['lam'] * p['v'] ** 2, inf.mu2))
        return out

    def grad(self, s, x):
        p, inf = self.p, self.inf
        d = s - p['sigma_c']
        Vs = p['g'] ** 2 * s * x * x + inf.c1 - inf.m2 * d
        if s > p['sigma_c'] and inf.loop:
            Vs += inf._cwd(s)[0]
        Vx = p['lam'] * (x * x - p['v'] ** 2) * x + p['g'] ** 2 * s * s * x - self.cb
        return Vs, Vx

    def mxx(self, s, x):
        p = self.p
        return p['lam'] * (3 * x * x - p['v'] ** 2) + p['g'] ** 2 * s * s

    def ms2(self, x, R):
        if self.smech == 'a':
            return (self.g_s * x) ** 2
        return self.m_s ** 2 + self.xi_s * R


def run(model, s_start, t_end_after_c, ks=(0.0, 0.3, 0.6, 1.0, 1.5), rtol=1e-9, max_step=None, nrec=20000, t_modes=0.2):
    """Integrate from sigma = s_start (slow-roll attractor) to t_c + t_end_after_c. Returns dict of arrays."""
    p, inf = model.p, model.inf
    mP = p['mP']
    V0 = inf.V(s_start)
    sd0 = -inf.dV(s_start) / math.sqrt(3 * V0) * mP      # slow roll: 3 H sdot = -V'
    x0 = model.cb / max(inf.M2(s_start), 1e-30)
    # mode normalisation: k in units of m_chi (comoving, a = 1 at the slice: fixed after the run, here a = 1 at start)
    nk = len(ks)

    def rhs(t, y):
        s, sd, x, xd, la, ss, ssd, Ur = y[:8]
        V = model.V(s, x)
        Vs, Vx = model.grad(s, x)
        K = 0.5 * (sd * sd + xd * xd)
        rr = max(Ur, 0.0) * math.exp(-4 * la)
        H = math.sqrt((K + V + rr) / (3 * mP ** 2))
        eps = (K + 2 * rr / 3) / (mP ** 2 * H * H)
        R = 6 * H * H * (2 - eps)
        ms2 = model.ms2(x, R)
        Gs = Gc = 0.0
        if model._decay_on:
            Es = 0.5 * sd * sd + 0.5 * p['g'] ** 2 * x * x * (s - model._s0) ** 2
            Ec = max(K + V - Es, 0.0)
            Gs = sum(c.rate(Es) for c in model.chans_sig)
            Gc = sum(c.rate(Ec) for c in model.chans_chi)
        out = [sd, -(3 * H + Gs) * sd - Vs, xd, -(3 * H + Gc) * xd - Vx, H, ssd, -3 * H * ssd - ms2 * ss,
               (Gs * sd * sd + Gc * xd * xd) * math.exp(4 * la)]
        if not model._modes_on:
            return out + [0.0] * (2 * nk)
        a = math.exp(la - model._lac)
        m2 = model.mxx(s, x)
        for i, k in enumerate(ks):
            u, ud = y[8 + 2 * i], y[9 + 2 * i]
            kk = (k * p['m_chi'] / a) ** 2
            out += [ud, -3 * H * ud - (kk + m2) * u]
        return out

    model._modes_on, model._lac = False, 0.0
    model._decay_on = False
    model._s0 = model.vac_sigma()
    if not hasattr(model, 'chans_sig'):
        model.chans_sig, model.chans_chi = [], []
    y0 = [s_start, sd0, x0, 0.0, 0.0, 1.0, 0.0, 0.0]
    for k in ks:
        y0 += [1.0, 0.0]
    evc = lambda t, y: y[0] - p['sigma_c']
    evc.terminal = True
    evc.direction = -1
    # phase 1: to sigma_c (chi heavy, adiabatic at chi_min)
    tc_guess = 1e3
    sol1 = solve_ivp(rhs, (0, tc_guess), y0, method='DOP853', rtol=rtol, atol=1e-12, events=evc,
                     dense_output=True, max_step=0.05)
    tc = sol1.t[-1] if sol1.status == 1 else None
    if sol1.t_events[0].size:
        tc = sol1.t_events[0][0]
    yc = sol1.sol(tc)
    # reset the linear modes at the slice to unit amplitude (growth is measured from the slice)
    yc = list(yc)
    for i in range(nk):
        yc[8 + 2 * i], yc[9 + 2 * i] = 1.0, 0.0
    model._modes_on, model._lac = True, float(yc[4])
    model._decay_on = True
    ms_ = max_step or 0.05 / p['m_chi']

    def big(t, y):
        return max(abs(y[8 + 2 * i]) for i in range(nk)) - 1e40
    big.terminal = True
    t_modes = min(t_modes, t_end_after_c)
    sol2a = solve_ivp(rhs, (tc, tc + t_modes), yc, method='DOP853', rtol=rtol, atol=1e-10, dense_output=True,
                      max_step=ms_, events=big)
    ya = list(sol2a.y[:, -1])
    t_a = sol2a.t[-1]
    model._modes_on = False
    sol2b = solve_ivp(rhs, (t_a, tc + t_end_after_c), ya, method='DOP853', rtol=rtol, atol=1e-10, dense_output=True,
                      max_step=ms_)

    class _S:
        def sol(self, tt):
            tt = np.atleast_1d(tt)
            out = np.empty((len(ya), len(tt)))
            m1 = tt <= t_a
            if m1.any():
                out[:, m1] = sol2a.sol(tt[m1])
            if (~m1).any():
                out[:, ~m1] = sol2b.sol(tt[~m1])
                out[8:, ~m1] = np.nan
            return out
    sol2 = _S()
    # records
    t1 = np.linspace(0, tc, 600)
    t2 = np.linspace(tc, tc + t_end_after_c, nrec)
    Y1 = sol1.sol(t1)
    Y2 = sol2.sol(t2)
    rec = {}
    T = np.concatenate([t1, t2[1:]])
    Y = np.concatenate([Y1, Y2[:, 1:]], axis=1)
    la_c = float(yc[4])
    rec['t'] = T - tc
    rec['lna'] = Y[4] - la_c
    s, sd, x, xd = Y[0], Y[1], Y[2], Y[3]
    Vv = np.array([model.V(a_, b_) for a_, b_ in zip(s, x)])
    Kin = 0.5 * (sd ** 2 + xd ** 2)
    rr = np.maximum(Y[7], 0) * np.exp(-4 * Y[4])
    H = np.sqrt((Kin + Vv + rr) / (3 * mP ** 2))
    eps = (Kin + 2 * rr / 3) / (mP ** 2 * H ** 2)
    rec.update(sigma=s, sigmadot=sd, chi=x, chidot=xd, H=H, eps=eps, R=6 * H ** 2 * (2 - eps), V=Vv,
               Ksig=0.5 * sd ** 2, Kchi=0.5 * xd ** 2, s_over_sstar=Y[5], sdot_over_sstar=Y[6])
    # chi-only and sigma-only energy split (potential split: V(sigma, v) vs V(0, chi) parts are not separable; use the
    # oscillation energies about the vacuum: E_sig = Ksig + (1/2) m_sig^2 (sigma - sigma0)^2 ; E_chi = rest)
    s0 = model.vac_sigma()
    Esig = rec['Ksig'] + 0.5 * p['m_sig'] ** 2 * (s - s0) ** 2 * (x / p['v']) ** 2
    rec['E_sig'] = Esig
    rec['E_tot'] = Kin + Vv
    rec['rho_r'] = rr
    rec['modes'] = {str(k): Y[8 + 2 * i] for i, k in enumerate(ks)}
    rec['t_c'] = tc
    rec['t_modes_end'] = t_a - tc
    rec['ks'] = list(ks)
    return rec, sol1, sol2
