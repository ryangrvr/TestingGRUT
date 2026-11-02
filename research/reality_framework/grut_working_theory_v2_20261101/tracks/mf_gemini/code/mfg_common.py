"""mf_gemini common module: the external two-field waterfall proposal, as written.

V(sigma, chi) = lam/4 (chi^2 - v^2)^2 + 1/2 g^2 sigma^2 chi^2 + V_SR(sigma)

Units (code): reduced M_P = 1 and H0 = 1, where H0 = sqrt(V0/3) is the Hubble rate on the
inflationary valley (V0 = lam v^4/4 = 3 M_P^2 H0^2, the proposal's normalisation).  Fields in M_P,
time tau = H0 t, energies in M_P^2 H0^2.  Physical H0 (GeV) only enters through the quantum noise
H/2pi relative to M_P (hub = H0/M_P) and the A_s-fixed slope.

Proposal relations [D]: m_chi = sqrt(2 lam) v = 3 H0 and lam v^4/4 = 3 H0^2  ->  v^2 = 8/3,
lam = 27/16 (H0/M_P)^2 (code units: lam = 27/16).  The waterfall point: g^2 sigma_c^2 = lam v^2.

Every choice here (V_SR family, sigma_c, slope, end-of-inflation surface) is a model choice.
Imports: numpy, scipy only.  Run with python3 -B.
"""
import numpy as np
from scipy.integrate import solve_ivp

MP_GEV = 2.435e18          # reduced Planck mass [R]
AS = 2.100e-9              # Planck 2018 A_s (ln 1e10 A_s = 3.044) [R]
NS_PLANCK = (0.9649, 0.0042)   # Planck 2018 TT,TE,EE+lowE+lensing [R]
NS_ACT = (0.9743, 0.0034)      # ACT DR6 P-ACT-LB, arXiv:2503.14452 [R]


def eps_star(H_GeV):
    """slow-roll epsilon at the pivot for A_s with H_* = H (V0-dominated)."""
    hub = H_GeV / MP_GEV
    return hub**2 / (8 * np.pi**2 * AS)


def r_of_H(H_GeV):
    return 16 * eps_star(H_GeV)


class Model:
    """Two-field hybrid in code units (M_P = H0 = 1).

    vsr: dict with 'kind' in {'lin','hill','log','quad','none'} and parameters:
      lin : V_SR = 3 p (sigma - sigma_c)                        (slope only, p = sqrt(2 eps_c))
      hill: V_SR = 3 [-eta/2 (sigma - sm)^2 + eta/2 sm^2]       (max at sm > sigma_c; red tilt)
      log : V_SR = 3 (kap2/16 pi^2) ln((sigma^2 + s0^2)/sigma_c^2)  (CW-like; red tilt)
      quad: V_SR = 3 eta/2 sigma^2                              (convex; blue tilt)
    """

    def __init__(self, sigc, vsr, v2=8.0 / 3.0, lam=None, hub=1.3e11 / MP_GEV):
        self.v2 = v2
        self.v = np.sqrt(v2)
        # lam fixed by lam v^4/4 = 3 (V0 = 3 in code units) unless given
        self.lam = 12.0 / v2**2 if lam is None else lam
        self.sigc = sigc
        self.g2 = self.lam * v2 / sigc**2
        self.g = np.sqrt(self.g2)
        self.vsr = dict(vsr)
        self.hub = hub   # H0/M_P (sets the quantum noise H/2pi in code units)

    # --- V_SR and derivatives ---
    def VSR(self, s):
        k = self.vsr['kind']
        if k == 'lin':
            return 3 * self.vsr['p'] * (s - self.sigc)
        if k == 'hill':
            e, sm = self.vsr['eta'], self.vsr['sm']
            return 3 * (-0.5 * e * (s - sm)**2 + 0.5 * e * sm**2)
        if k == 'log':
            return 3 * self.vsr['kap2'] / (16 * np.pi**2) * np.log((s**2 + self.vsr['s0']**2) / self.sigc**2)
        if k == 'quad':
            return 3 * 0.5 * self.vsr['eta'] * s**2
        return 0.0 * s

    def dVSR(self, s):
        k = self.vsr['kind']
        if k == 'lin':
            return 3 * self.vsr['p'] + 0.0 * s
        if k == 'hill':
            return -3 * self.vsr['eta'] * (s - self.vsr['sm'])
        if k == 'log':
            return 3 * self.vsr['kap2'] / (16 * np.pi**2) * 2 * s / (s**2 + self.vsr['s0']**2)
        if k == 'quad':
            return 3 * self.vsr['eta'] * s
        return 0.0 * s

    # --- full potential ---
    def V(self, s, c):
        return 0.25 * self.lam * (c**2 - self.v2)**2 + 0.5 * self.g2 * s**2 * c**2 + self.VSR(s)

    def Vs(self, s, c):
        return self.g2 * s * c**2 + self.dVSR(s)

    def Vc(self, s, c):
        return self.lam * (c**2 - self.v2) * c + self.g2 * s**2 * c

    def mchi2(self, s, c=0.0):
        return self.lam * (3 * c**2 - self.v2) + self.g2 * s**2

    # --- background in cosmic time: y = (N, s, ds, c, dc) ---
    def rhs(self, t, y):
        N, s, ds, c, dc = y
        h = np.sqrt((0.5 * ds**2 + 0.5 * dc**2 + self.V(s, c)) / 3.0)
        return [h, ds, -3 * h * ds - self.Vs(s, c), dc, -3 * h * dc - self.Vc(s, c)]

    def hubble(self, y):
        N, s, ds, c, dc = y
        return np.sqrt((0.5 * ds**2 + 0.5 * dc**2 + self.V(s, c)) / 3.0)

    def epsH(self, y):
        N, s, ds, c, dc = y
        h = self.hubble(y)
        return 0.5 * (ds**2 + dc**2) / h**2

    def run(self, y0, t1, events=None, max_step=np.inf, rtol=1e-10, atol=1e-13, dense=False):
        return solve_ivp(self.rhs, (0.0, t1), y0, method='DOP853', rtol=rtol, atol=atol,
                         events=events, max_step=max_step, dense_output=dense)


def slow_roll_ic(m, s0, c0=0.0, dc0=0.0):
    """initial state on the valley chi ~ 0 with slow-roll sigma velocity."""
    h = np.sqrt(m.V(s0, c0) / 3.0)
    ds0 = -m.Vs(s0, c0) / (3 * h)
    return np.array([0.0, s0, ds0, c0, dc0])


def end_of_inflation(m, y0, t1=400.0):
    """integrate to eps_H = 1 (first crossing) and return (sol, state at end)."""
    def ev(t, y):
        return m.epsH(y) - 1.0
    ev.terminal = True
    ev.direction = 1
    sol = m.run(y0, t1, events=ev)
    if sol.t_events[0].size:
        return sol, sol.y_events[0][0], sol.t_events[0][0]
    return sol, None, None


# ---------------------------------------------------------------------------------------------
# Two-stage background: full (sigma, chi) equations until sigma is heavy (g|chi| = KSW h), then the
# sigma oscillation is replaced by its adiabatic invariant J = a^3 E_osc/m_sigma (WKB average):
#   <sigma^2> = sbar^2 + J/(m_sigma a^3),  rho_osc = J m_sigma/a^3 (w = 0),
#   sbar = -V_SR'(0)/(g^2 chi^2 + V_SR''(0)) (minimum of the sigma potential at fixed chi).
# Valid while m_sigma >> H and |d m_sigma/dt| << m_sigma^2 (checked a posteriori: adiab).
# ---------------------------------------------------------------------------------------------
KSW = 30.0


def d2VSR(m, s):
    k = m.vsr['kind']
    if k == 'hill':
        return -3 * m.vsr['eta'] + 0.0 * s
    if k == 'quad':
        return 3 * m.vsr['eta'] + 0.0 * s
    if k == 'log':
        s0 = m.vsr['s0']
        return 3 * m.vsr['kap2'] / (16 * np.pi**2) * 2 * (s0**2 - s**2) / (s**2 + s0**2)**2
    return 0.0 * s


def _avg_terms(m, N, c, J):
    msig2 = m.g2 * c**2 + d2VSR(m, 0.0)
    msig2 = max(msig2, 1e-300)
    sbar = -m.dVSR(0.0) / msig2
    msig = np.sqrt(msig2)
    a3 = np.exp(3 * N)
    s2 = sbar**2 + J / (msig * a3)
    rho_osc = J * msig / a3
    return sbar, msig, s2, rho_osc


def avg_rhs(m, J):
    def f(t, y):
        N, c, dc = y
        sbar, msig, s2, rho_osc = _avg_terms(m, N, c, J)
        V = 0.25 * m.lam * (c**2 - m.v2)**2 + 0.5 * m.g2 * sbar**2 * c**2 + m.VSR(sbar)
        h = np.sqrt((0.5 * dc**2 + V + rho_osc) / 3.0)
        force = m.lam * (c**2 - m.v2) * c + m.g2 * c * s2
        return [h, dc, -3 * h * dc - force]
    return f


def avg_state(m, J, y):
    N, c, dc = y
    sbar, msig, s2, rho_osc = _avg_terms(m, N, c, J)
    V = 0.25 * m.lam * (c**2 - m.v2)**2 + 0.5 * m.g2 * sbar**2 * c**2 + m.VSR(sbar)
    h = np.sqrt((0.5 * dc**2 + V + rho_osc) / 3.0)
    eps = (dc**2 + rho_osc) / (2 * h**2)
    return dict(h=h, eps=eps, rho_chi=0.5 * dc**2 + 0.25 * m.lam * (c**2 - m.v2)**2,
                rho_sig=rho_osc + 0.5 * m.g2 * sbar**2 * c**2 + m.VSR(sbar), msig=msig)


def two_stage(m, y0, t_after_end=0.0, t_max=400.0, record=False, rtol=1e-10):
    """returns dict with N_end, state at end, and (optionally) the post-end chi(t) record.
    Stops with N_end = nan if sigma leaves the tachyonic band on the far side (sigma < -1.02 sigma_c)
    before sigma becomes heavy (the waterfall then does not complete in this crossing)."""
    def ev_exit(t, y):
        return y[1] + 1.02 * m.sigc
    ev_exit.terminal = True
    ev_exit.direction = -1

    def ev_heavy(t, y):
        return m.g * abs(y[3]) - KSW * m.hubble(y)
    ev_heavy.terminal = True
    ev_heavy.direction = 1

    def ev_end(t, y):
        return m.epsH(y) - 1.0
    ev_end.terminal = True
    ev_end.direction = 1
    s1 = m.run(y0, t_max, events=[ev_heavy, ev_end, ev_exit], rtol=rtol, atol=1e-18)
    out = dict(stage1_t=float(s1.t[-1]))
    if s1.t_events[2].size:
        out.update(N_end=np.nan, switched=False, band_exit=True)
        return out
    if s1.t_events[1].size:            # inflation ended before sigma became heavy
        ye = s1.y_events[1][0]
        out.update(N_end=float(ye[0]), switched=False, end_state=dict(sigma=float(ye[1]), chi=float(ye[3])),
                   h_end=float(m.hubble(ye)))
        return out
    if not s1.t_events[0].size:
        out.update(N_end=np.nan, switched=False)
        return out
    y1 = s1.y_events[0][0]
    t1 = s1.t_events[0][0]
    N, s, ds, c, dc = y1
    msig2 = m.g2 * c**2 + d2VSR(m, 0.0)
    sbar = -m.dVSR(0.0) / msig2
    Eosc = 0.5 * ds**2 + 0.5 * msig2 * (s - sbar)**2
    J = Eosc * np.exp(3 * N) / np.sqrt(msig2)
    f = avg_rhs(m, J)

    def ev_end2(t, y):
        return avg_state(m, J, y)['eps'] - 1.0
    ev_end2.terminal = True
    ev_end2.direction = 1
    s2 = solve_ivp(f, (t1, t1 + t_max), [N, c, dc], method='DOP853', rtol=rtol, atol=1e-18,
                   events=[ev_end2])
    if not s2.t_events[0].size:
        out.update(N_end=np.nan, switched=True)
        return out
    ye = s2.y_events[0][0]
    te = s2.t_events[0][0]
    st = avg_state(m, J, ye)
    out.update(N_end=float(ye[0]), switched=True, N_switch=float(N), J=float(J),
               end_state=dict(chi=float(ye[1]), dchi=float(ye[2]), chi_over_v=float(ye[1] / m.v),
                              rho_chi_frac=float(st['rho_chi'] / (3 * st['h']**2)),
                              rho_sig_frac=float(st['rho_sig'] / (3 * st['h']**2)),
                              msig_over_h=float(st['msig'] / st['h'])),
               h_end=float(st['h']))
    if record and t_after_end > 0:
        s3 = solve_ivp(f, (te, te + t_after_end), ye, method='DOP853', rtol=rtol, atol=1e-18,
                       dense_output=True, max_step=0.02)
        out['post'] = s3
        out['post_J'] = J
    return out
