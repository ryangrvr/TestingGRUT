"""Post-waterfall history: condensates sigma, chi -> singlet fermions N (axial couplings) -> L H, spectator s [C/E].

Fluids (cosmic time, H_i units, a = 1 at the slice): rho_sig, rho_chi (matter-like, oscillation amplitude
Sigma = sqrt(2 rho)/m), rho_N (relativistic N), rho_SM (N -> L H), rho_s (exact s field to t_wkb, then WKB).
Decay of each condensate through each axial species j (M_j, f_j):
  Gamma_j = (m/pi) eps_j(A)/(S^2/2)          non-perturbative, A = Sigma/(2 f_j), S = Sigma/m, eps_j from axial.py tables
          + (M_j/f_j)^2 m beta/(16 pi)        perturbative pseudoscalar width once Sigma < f_j (and 2 M_j < m)
The N decay to L H at Gamma_N = mt M^2/(8 pi v_EW^2) (seesaw, rest frame) times M/E_N, E_N redshifting from the
production energy ~ (3/4) mu5_max.  These rates are [E] built on the [C] half-period production of axial.py.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import axial

V_EW = 174.1
MNU = 0.05e-9   # GeV


def eps_table(Mh, As=None):
    As = As if As is not None else np.geomspace(0.3, 80.0, 28)
    eps = [axial.half_period(float(A), Mh, nk=900)['eps'] for A in As]
    return dict(Mh=Mh, A=list(map(float, As)), eps=list(map(float, eps)))


class Channel:
    def __init__(self, m, M, f, tab):
        self.m, self.M, self.f = m, M, f
        A = np.array(tab['A'])
        e = np.maximum(np.array(tab['eps']), 1e-300)
        self.lA, self.le = np.log(A), np.log(e)
        self._f = interp1d(self.lA, self.le, kind='linear', bounds_error=False, fill_value=(self.le[0], self.le[-1]))
        # small-A extrapolation: adiabatic suppression; use the tabulated slope below the grid
        self.slope_lo = (self.le[1] - self.le[0]) / (self.lA[1] - self.lA[0])
        self.slope_hi = (self.le[-1] - self.le[-2]) / (self.lA[-1] - self.lA[-2])
        self.Gp = axial.gamma_pert(M, f, m)

    def rate(self, rho):
        if rho <= 0:
            return 0.0
        Sig = math.sqrt(2 * rho) / self.m          # field amplitude
        S = Sig / self.m                            # amplitude in units of m
        A = Sig / (2 * self.f)
        lA = math.log(A)
        if lA < self.lA[0]:
            return self.Gp           # small amplitude: the perturbative pseudoscalar width (zero if 2M > m)
        if lA > self.lA[-1]:
            le = self.le[-1] + self.slope_hi * (lA - self.lA[-1])
        else:
            le = float(self._f(lA))
        return (self.m / math.pi) * math.exp(le) / (0.5 * S * S)


def run_fluid(state0, mP, chans_sig, chans_chi, H_GeV, s_model, t_end=3e4, t_wkb=60.0, rtol=1e-8):
    """state0: dict(t, lna, rho_sig, rho_chi, rho_r, s, sd).  s_model(t, rho_chi, Rbar) -> m_s^2.
    Returns record arrays."""
    MN = [c.M for c in chans_sig + chans_chi]
    GN_rest = [MNU * (M * H_GeV) ** 2 / (8 * math.pi * V_EW ** 2) / H_GeV for M in MN]
    Mavg = np.mean(MN)
    GNr = float(np.mean(GN_rest))

    def rhs(t, y, wkb):
        la, us, uc, ur, usm, s, sd, lEN = y
        a = math.exp(la)
        rs, rc = max(us, 0) / a ** 3, max(uc, 0) / a ** 3
        rr, rsm = max(ur, 0) / a ** 4, max(usm, 0) / a ** 4
        rho = rs + rc + rr + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        Gs = sum(c.rate(rs) for c in chans_sig)
        Gc = sum(c.rate(rc) for c in chans_chi)
        EN = math.exp(lEN)
        GNlab = GNr * min(1.0, Mavg / EN)
        Rbar = (rs + rc) / mP ** 2
        out = [H, -Gs * us, -Gc * uc, (Gs * rs + Gc * rc) * a ** 4 - GNlab * ur, GNlab * ur]
        if not wkb:
            ms2 = s_model(t, rc, Rbar)
            out += [sd, -3 * H * sd - ms2 * s]
        else:
            out += [0.0, 0.0]
        out += [-H]   # E_N redshifts as 1/a (new N are injected at higher energy; this is the old-population clock)
        return out

    y0 = [state0['lna'], state0['rho_sig'] * math.exp(3 * state0['lna']), state0['rho_chi'] * math.exp(3 * state0['lna']),
          state0['rho_r'] * math.exp(4 * state0['lna']), 0.0, state0['s'], state0['sd'], math.log(state0['E_N'])]
    t0 = state0['t']
    tw = max(t_wkb, t0 + 1.0)
    sol1 = solve_ivp(lambda t, y: rhs(t, y, False), (t0, tw), y0, method='DOP853', rtol=rtol, atol=1e-30,
                     dense_output=True, max_step=0.02)
    y1 = sol1.y[:, -1].copy()
    # WKB amplitude of s at tw
    ms2 = s_model(tw, max(y1[2], 0) / math.exp(3 * y1[0]), max(y1[1] + y1[2], 0) / math.exp(3 * y1[0]) / mP ** 2)
    ms = math.sqrt(max(ms2, 1e-30))
    amp_w = math.sqrt(y1[5] ** 2 + y1[6] ** 2 / ms ** 2)
    sol2 = solve_ivp(lambda t, y: rhs(t, y, True), (tw, t_end), y1, method='LSODA', rtol=rtol, atol=1e-30,
                     dense_output=True)
    T1 = np.linspace(t0, tw, 3000)
    T2 = np.geomspace(tw, t_end, 1500)[1:]
    Y1, Y2 = sol1.sol(T1), sol2.sol(T2)
    T = np.concatenate([T1, T2])
    Y = np.concatenate([Y1, Y2], axis=1)
    a = np.exp(Y[0])
    rs, rc = np.maximum(Y[1], 0) / a ** 3, np.maximum(Y[2], 0) / a ** 3
    rN, rsm = np.maximum(Y[3], 0) / a ** 4, np.maximum(Y[4], 0) / a ** 4
    rho = rs + rc + rN + rsm
    H = np.sqrt(rho / (3 * mP ** 2))
    rr = rN + rsm
    eps = 1.5 + 0.5 * rr / rho
    s = Y[5].copy()
    # WKB: amplitude ~ amp_w (a_w/a)^{3/2} (m_s constant after the turn-on)
    aw = math.exp(y1[0])
    mask = T > tw
    s_env = np.where(mask, amp_w * (aw / a) ** 1.5, np.nan)
    s[mask] = np.nan
    ms2 = np.array([s_model(t, c, (x + c) / mP ** 2) for t, c, x in zip(T, rc, rs)])
    rec = dict(t=T, lna=Y[0], a=a, H=H, eps=eps, rho_sig=rs, rho_chi=rc, rho_N=rN, rho_SM=rsm,
               R_mean=(rs + rc) / mP ** 2, R_osc_sig=3 * rs / mP ** 2, R_osc_chi=3 * rc / mP ** 2,
               s_over_sstar=s, s_env=s_env, m_s=np.sqrt(np.maximum(ms2, 0)), K=a ** 2 * H)
    rec['Gamma_sig'] = np.array([sum(c.rate(x) for c in chans_sig) for x in rs])
    rec['Gamma_chi'] = np.array([sum(c.rate(x) for c in chans_chi) for x in rc])
    rec['GammaN_rest_avg'] = GNr
    return rec
