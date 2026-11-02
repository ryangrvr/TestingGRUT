#!/usr/bin/env python3
"""inf_potential track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026).

Common module: explicit single-field inflaton potentials for M1r, background (slow roll and full ODE), exact
Mukhanov-Sasaki power spectra, the pivot e-fold number N_* from the post-inflation history, and the post-inflation
oscillation with a phenomenological perturbative width (the engine's convention).

Units: reduced Planck mass M_P = 1 (2.435e18 GeV, the repository's MPL).  V = V0 * v(phi).  In e-fold time N (code
variable N_c, N_c = 0 at the end of inflation eps_H = 1, negative before), the background shape phi(N) does not depend on
V0, and H = sqrt(V0) * h(N).  P_zeta = V0 * Phat_zeta (Phat computed at V0 = 1).

Potential families (every one is a model choice; none follows from closure):
  'T'  : alpha-attractor T-model, n = 1 (Kallosh-Linde 1306.5220):  v = tanh^2(phi / M),  M = sqrt(6 alpha).
  'E'  : alpha-attractor E-model, n = 1 (Kallosh-Linde 1306.5220; Starobinsky at alpha = 1):  v = (1 - exp(-phi/M))^2,
         M = sqrt(3 alpha / 2).
  'HT' : hilltop with a quadratic minimum (Boubekeur-Lyth hep-ph/0502047 squared form; Kohri-Lin-Lyth 0707.3826):
         v = f(x)^2 / f(x_top)^2 with f = 1 - beta x - x^p, x = phi / vev (so v = 1 at the top of the hill).
         beta = 0 is the pure BL form (top at x = 0); beta > 0 is a small linear tilt (top at x_top < 0), our deformation
         used only to move n_s.  The minimum is at f = 0 (x_min slightly below 1 for beta > 0).
  'PL' : polynomial plateau / KKLTI-type (Kachru et al. hep-th/0308055 brane form; Kallosh-Linde-Roest polynomial
         attractors): v = phi^q / (phi^q + mu^q).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MPL_GEV = 2.435e18                  # reduced Planck mass [Q repository convention]
GEV_TO_INV_MPC = 1.56374e38         # 1 GeV in Mpc^-1 (hbar c = 1.97327e-16 GeV m, 1 Mpc = 3.08568e22 m) [R]
T0_GEV = 2.7255 * 8.617333e-14      # CMB temperature today [R]
GS0 = 3.909                         # entropy dof today (3 nu at T_nu = (4/11)^{1/3} T) [R]
G_RH = 106.75                       # SM dof at reheating [R]; the hidden sector is neglected (cold, tiny fraction)
K_PIVOT = 0.05                      # Mpc^-1
AS_TARGET = 2.100e-9                # Planck 2018 [R], brief value
NS_TARGET, NS_SIG = 0.9649, 0.0042  # Planck 2018 TT,TE,EE+lowE+lensing [R]


class Model:
    def __init__(self, kind, **p):
        self.kind = kind
        self.p = dict(p)
        if kind == 'T':
            self.M = math.sqrt(6 * p['alpha'])
        elif kind == 'E':
            self.M = math.sqrt(1.5 * p['alpha'])
        elif kind == 'HT':
            self.vev = p['vev']; self.pw = p['p']; self.beta = p.get('beta', 0.0)
            b, pw = self.beta, self.pw
            # top of the hill: f'(x) = -b - p x^{p-1} = 0  ->  x_top = -(b/p)^{1/(p-1)} (p even); minimum: f(x_min) = 0
            self.xtop = -(b / pw) ** (1.0 / (pw - 1)) if b > 0 else 0.0
            self.ftop = 1 - b * self.xtop - self.xtop ** pw
            self.xmin = brentq(lambda x: 1 - b * x - x ** pw, 0.5, 1.0)
        elif kind == 'PL':
            self.mu = p['mu']; self.q = p['q']
        else:
            raise ValueError(kind)

    # ---- potential shape v(phi) (V0 = 1), its first and second derivatives ----
    def v(self, f):
        if self.kind == 'T':
            return math.tanh(f / self.M) ** 2
        if self.kind == 'E':
            return (-math.expm1(-f / self.M)) ** 2
        if self.kind == 'HT':
            x = f / self.vev
            return (1 - self.beta * x - x ** self.pw) ** 2 / self.ftop ** 2
        if self.kind == 'PL':
            r = (abs(f) / self.mu) ** self.q
            return r / (1 + r)

    def dv(self, f):
        if self.kind == 'T':
            y = f / self.M
            if abs(y) > 300:
                return 0.0
            c = math.cosh(y)
            return 2 * math.tanh(y) / (c * c) / self.M
        if self.kind == 'E':
            e = math.exp(-f / self.M)
            return 2 * (1 - e) * e / self.M
        if self.kind == 'HT':
            x = f / self.vev; pw = self.pw
            fx = 1 - self.beta * x - x ** pw
            return 2 * fx * (-self.beta - pw * x ** (pw - 1)) / self.vev / self.ftop ** 2
        if self.kind == 'PL':
            q, mu = self.q, self.mu
            r = (abs(f) / mu) ** q
            return math.copysign(q * r / (abs(f) * (1 + r) ** 2), f) if f != 0 else 0.0

    def d2v(self, f):
        if self.kind == 'T':
            y = f / self.M
            if abs(y) > 300:
                return 0.0
            c = math.cosh(y); t = math.tanh(y)
            s2 = 1 / (c * c)
            return 2 * (s2 * s2 - 2 * t * t * s2) / self.M ** 2
        if self.kind == 'E':
            e = math.exp(-f / self.M)
            return 2 * e * (2 * e - 1) / self.M ** 2
        if self.kind == 'HT':
            x = f / self.vev; pw = self.pw
            fx = 1 - self.beta * x - x ** pw
            fp = -self.beta - pw * x ** (pw - 1)
            fpp = -pw * (pw - 1) * x ** (pw - 2)
            return 2 * (fp * fp + fx * fpp) / self.vev ** 2 / self.ftop ** 2
        if self.kind == 'PL':
            q, mu = self.q, self.mu
            af = abs(f)
            if af == 0:
                return 2 / mu ** 2 if q == 2 else 0.0
            r = (af / mu) ** q
            # d/dphi [q r / (phi (1+r)^2)]
            return q * r * ((q - 1) * (1 + r) - 2 * q * r) / (af * af * (1 + r) ** 3)

    # ---- geometry of the model ----
    def phi_min(self):
        """position of the minimum (where the holder field phi_I = phi - phi_min vanishes)"""
        if self.kind == 'HT':
            return self.xmin * self.vev
        return 0.0

    def m2_min(self):
        """V''(phi_min)/V0: oscillation mass^2 at the minimum in V0 units"""
        if self.kind == 'PL' and self.q != 2:
            return 0.0
        return self.d2v(self.phi_min())

    def phi_start_guess(self, Ntot=80.0):
        """field value roughly Ntot e-folds before the end (slow roll), on the inflationary side"""
        if self.kind == 'T':
            a = self.p['alpha']
            return 0.5 * self.M * math.acosh(1 + 4 * Ntot / (3 * a))
        if self.kind == 'E':
            a = self.p['alpha']
            return self.M * math.log(1 + 4 * Ntot / (3 * a))
        if self.kind == 'PL':
            q, mu = self.q, self.mu
            # phi^{q+2} ~ q (q+2) mu^q Ntot  (v ~ 1 - (mu/phi)^q)
            return (q * (q + 2) * mu ** q * Ntot) ** (1.0 / (q + 2)) * 1.5
        if self.kind == 'HT':
            # slow-roll e-fold integral N(phi) = int v/(-v') dphi from phi to the eps_V = 1 point, on a grid in
            # d = phi - phi_top (log spaced); pick the point with N = Ntot
            ft = self.xtop * self.vev
            ds = np.logspace(math.log10(self.vev) - 12, math.log10(0.999 * (self.xmin - self.xtop) * self.vev), 4000)
            fs = ft + ds
            g = np.array([self.v(f) / max(-self.dv(f), 1e-300) for f in fs])
            ev = np.array([0.5 * (self.dv(f) / self.v(f)) ** 2 for f in fs])
            iend = int(np.argmax(ev >= 1.0)) if np.any(ev >= 1.0) else len(fs) - 1
            # cumulative from the end backwards
            seg = 0.5 * (g[1:] + g[:-1]) * np.diff(fs)
            Ncum = np.concatenate([np.cumsum(seg[:iend][::-1])[::-1], [0.0]])
            j = int(np.argmin(np.abs(Ncum - Ntot)))
            return fs[j]
        raise ValueError


def bg_rhs(N, y, m):
    f, fN = y
    eps = 0.5 * fN * fN
    vv = m.v(f)
    return [fN, -(3 - eps) * (fN + m.dv(f) / vv)]


def background(m, Nwant=75.0, phi0=None, rtol=1e-11, atol=1e-13):
    """Full background ODE in N from deep slow roll to eps_H = 1.  Returns dict with arrays on N_c (0 at the end)."""
    side = 1.0
    if phi0 is None:
        # find a start with at least Nwant e-folds: integrate from a guess, enlarge until enough
        f0 = m.phi_start_guess(Nwant)
        for it in range(25):
            r = _bg_integrate(m, f0, rtol, atol)
            if r is not None and r['Ntot'] >= Nwant:
                break
            f0 = _enlarge(m, f0)
        else:
            raise RuntimeError('no start with enough e-folds')
    else:
        r = _bg_integrate(m, phi0, rtol, atol)
    if r['Ntot'] > Nwant + 15:
        # restart on the attractor Nwant + 5 e-folds before the end (keeps the dense solution short)
        f1, fN1 = r['sol'].sol(r['Nend'] - (Nwant + 5))
        r = _bg_integrate(m, f1, rtol, atol, fN0=fN1)
    return r


def _enlarge(m, f0):
    if m.kind in ('T', 'E'):
        return f0 + 0.5 * m.M
    if m.kind == 'PL':
        return f0 * 1.3
    if m.kind == 'HT':
        # move toward the top
        xt = m.xtop * m.vev
        return xt + 0.5 * (f0 - xt)


def _bg_integrate(m, f0, rtol, atol, fN0=None):
    """Integrate from phi = f0 at slow-roll velocity (or the given fN0) to eps_H = 1."""
    if fN0 is None:
        fN0 = -m.dv(f0) / m.v(f0)
    # inflationary side: for HT the field rolls to larger phi (dv < 0 there), for T/E/PL to smaller phi
    ev = lambda N, y, m=m: 0.5 * y[1] ** 2 - 1.0
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(bg_rhs, (0.0, 400.0), [f0, fN0], args=(m,), method='DOP853', rtol=rtol, atol=atol,
                    events=ev, dense_output=True, max_step=0.05)
    if sol.status != 1:
        return None
    Nend = sol.t_events[0][0]
    out = dict(sol=sol, Nend=Nend, Ntot=Nend, phi_end=sol.y_events[0][0][0], phiN_end=sol.y_events[0][0][1])
    return out


def bg_at(bg, Nbefore):
    """background state Nbefore e-folds before the end: phi, phi_N, eps, h = H/sqrt(V0) (V0 = 1), eta_V"""
    N = bg['Nend'] - Nbefore
    f, fN = bg['sol'].sol(N)
    return f, fN


def h_of(m, f, fN):
    eps = 0.5 * fN * fN
    return math.sqrt(m.v(f) / (3 - eps)), eps


def ms_rhs(N, y, m, kap):
    f, fN, ur, ui, uNr, uNi, tr, ti, tNr, tNi = y
    eps = 0.5 * fN * fN
    vv = m.v(f)
    H2 = vv / (3 - eps)
    fNN = -(3 - eps) * (fN + m.dv(f) / vv)
    k2 = kap * kap * math.exp(-2 * N) / H2
    meff = k2 + m.d2v(f) / H2 - ((3 - eps) * fN * fN + 2 * fN * fNN)
    fr = 3 - eps
    return [fN, fNN, uNr, uNi, -fr * uNr - meff * ur, -fr * uNi - meff * ui,
            tNr, tNi, -fr * tNr - k2 * tr, -fr * tNi - k2 * ti]


def ms_power(m, bg, kap, kin=200.0, kout=1e-3):
    """Exact scalar and tensor power at comoving kap = k / sqrt(V0) (a = 1 at the end, N_c = N - Nend), V0 = 1.
    Returns (Phat_zeta, Phat_t, N_exit_before_end)."""
    sol = bg['sol']; Nend = bg['Nend']
    # start where kap e^{-Nc}/h = kin
    def ratio(N):
        f, fN = sol.sol(N)
        h, _ = h_of(m, f, fN)
        return kap * math.exp(-(N - Nend)) / h
    # bracket: the mode exits at Nend - ln(kap/h) roughly; search within 20 e-folds before that
    f_, fN_ = sol.sol(Nend); h_end_, _ = h_of(m, f_, fN_)
    Nx_guess = Nend + math.log(kap / h_end_)
    lo = max(0.0, Nx_guess - 25.0)
    Ns = brentq(lambda N: math.log(ratio(N) / kin), lo, Nend)
    Nx = brentq(lambda N: math.log(ratio(N)), Ns, Nend + 1e-9) if ratio(Nend) < 1 else Nend
    # end: where ratio = kout or at the end of inflation (super-horizon freeze), whichever first
    Ne = brentq(lambda N: math.log(ratio(N) / kout), Ns, Nend) if ratio(Nend) < kout else Nend
    f, fN = sol.sol(Ns)
    h, eps = h_of(m, f, fN)
    a = math.exp(Ns - Nend)
    k = kap
    u0 = 1.0 / (a * math.sqrt(2 * k))
    x = k / (a * math.sqrt(m.v(f) / (3 - eps)))
    y0 = [f, fN, u0, 0.0, -u0, -x * u0, u0, 0.0, -u0, -x * u0]
    # shift N so that exp(-N) in ms_rhs is a^{-1}: use Nc = N - Nend
    def rhs(Nc, y):
        return ms_rhs(Nc, y, m, kap)
    s = solve_ivp(rhs, (Ns - Nend, Ne - Nend), y0, method='DOP853', rtol=1e-10, atol=1e-30)
    yf = s.y[:, -1]
    fN = yf[1]
    u2 = yf[2] ** 2 + yf[3] ** 2
    t2 = yf[6] ** 2 + yf[7] ** 2
    Pz = k ** 3 * u2 / (2 * math.pi ** 2 * fN * fN)
    Pt = 8 * k ** 3 * t2 / (2 * math.pi ** 2)
    return Pz, Pt, Nend - Nx


def kap_of_Nexit(m, bg, Nbefore):
    """comoving kap = a H / sqrt(V0) of the mode that exits Nbefore e-folds before the end"""
    f, fN = bg_at(bg, Nbefore)
    h, _ = h_of(m, f, fN)
    return math.exp(-Nbefore) * h


def spectrum_at_pivot(m, bg, Nstar, dlnk=0.25):
    """P_zeta, n_s, alpha_s, r at the mode exiting Nstar e-folds before the end (V0 = 1 normalisation for P)."""
    k0 = kap_of_Nexit(m, bg, Nstar)
    ks = [k0 * math.exp(j * dlnk) for j in (-2, -1, 0, 1, 2)]
    P = [ms_power(m, bg, k) for k in ks]
    lz = np.log([p[0] for p in P]); lt = np.log([p[1] for p in P])
    x = np.array([-2, -1, 0, 1, 2]) * dlnk
    cz = np.polyfit(x, lz, 3); ct = np.polyfit(x, lt, 2)
    ns = 1 + cz[2]; als = 2 * cz[1]
    nt = ct[1]
    Pz, Pt = P[2][0], P[2][1]
    return dict(Pz=Pz, Pt=Pt, r=Pt / Pz, ns=ns, alpha_s=als, n_t=nt, kap=k0)


def slow_roll_at(m, bg, Nbefore):
    f, fN = bg_at(bg, Nbefore)
    vv = m.v(f); e = 0.5 * (m.dv(f) / vv) ** 2; eta = m.d2v(f) / vv
    return dict(phi=f, eps_V=e, eta_V=eta, ns_SR=1 - 6 * e + 2 * eta, r_SR=16 * e, P_SR=vv / (24 * math.pi ** 2 * e))


# ---------------- post-inflation (cosmic time, units 1/H_end... we use units 1/sqrt(V0)) ----------------
def post_rhs(t, y, m, gam):
    lna, f, fd, rr = y
    rho = 0.5 * fd * fd + m.v(f) + rr
    H = math.sqrt(rho / 3)
    return [H, fd, -(3 * H + gam) * fd - m.dv(f), gam * fd * fd - 4 * H * rr]


def post_inflation(m, bg, gam, t_end, start_Nbefore=0.0, rtol=1e-10, atol=1e-14, max_step=None, dense=False):
    """Homogeneous post-inflation evolution with a constant perturbative width gam (units sqrt(V0)), from Nbefore
    e-folds before the end (a = 1 at the end of inflation).  Returns the solution object."""
    f, fN = bg_at(bg, start_Nbefore)
    h, eps = h_of(m, f, fN)
    fd = h * fN
    y0 = [-start_Nbefore, f, fd, 0.0]
    kw = {}
    if max_step:
        kw['max_step'] = max_step
    return solve_ivp(post_rhs, (0.0, t_end), y0, args=(m, gam), method='DOP853', rtol=rtol, atol=atol,
                     dense_output=dense, **kw)


def Nstar_from_history(Hstar_GeV, a_end_over_a0):
    return math.log(Hstar_GeV * GEV_TO_INV_MPC / K_PIVOT) + math.log(a_end_over_a0)


def a_end_over_a0_wT(rho_end_GeV4, T_rh, w, g=G_RH):
    """reheating after the end at constant w, then radiation with entropy conservation"""
    rho_rh = math.pi ** 2 / 30 * g * T_rh ** 4
    if rho_rh > rho_end_GeV4:
        rho_rh = rho_end_GeV4
        T_rh = (30 * rho_rh / (math.pi ** 2 * g)) ** 0.25
    Nre = math.log(rho_end_GeV4 / rho_rh) / (3 * (1 + w))
    return math.exp(-Nre) * (GS0 / g) ** (1 / 3) * T0_GEV / T_rh, Nre, T_rh


def T_instant(rho_end_GeV4, g=G_RH):
    return (30 * rho_end_GeV4 / (math.pi ** 2 * g)) ** 0.25
