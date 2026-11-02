"""COPY of fx_formation/code/thermo.py (qm_micro, 2 Nov 2026).  Change: RES points to fx_formation/results (read-only).

fx_formation item 1b (Claude, 1 Nov 2026): equilibrium of the hidden X/Xbar/phi plasma with M1 objects.

Units m_X = 1.  Labels in NOTES.md.  Pieces:
  * vapour: ideal relativistic Bose gases with the VACUUM masses outside objects (phi = v): X (mass 1, +mu), Xbar
    (mass 1, -mu), phi (mass sqrt(6) = sqrt(2 lambda_phi) v / m_X, mu = 0).  Thermal masses are omitted [E]
    (delta m_X^2 / T^2 ~ (lambda_X/3 + g^2/12) ~ 1e-3).
  * liquid (bulk Q-matter, the Q -> infinity end of the family): zero-T grand potential -P0(mu) =
    100 min_{S,F} [U(S,F) - mu^2 F^2 / 2] (m_X^4), plus a one-loop thermal pressure of its three excitation branches
    (rotating frame; the cubic of tracks/m1_vacuum/code/coldfluid_m1.py at the bulk point; one phonon, two gapped) [E].
  * coexistence mu_c(T): P0(mu) + p_liq,th(T) = p_vap(mu, T).
  * finite-Q objects (Saha): n(Q) = g_int (E(Q) T / 2 pi)^{3/2} exp[-(E(Q) - mu Q)/T], E(Q) from results/family.json
    (spline in ln Q; the three-term fit beyond 1e8), Q in [Q_min, Q_max]; g_int varied 1e-10 .. 1e10.
"""
import json, math, os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', 'fx_formation', 'results')
a_, b_, c_ = 0.1, 3.0, 0.2
MPHI = math.sqrt(6.0)
OM_B, S_B, F_B, N_B = 0.7758647074595397, 0.151172107845, 2.406476756778, 449.31336788454803


# ------------------------------------------------------------------ ideal Bose gas, one dof
def bose_dof(m, mu, T):
    """n, p, e, s for one bosonic dof of mass m and chemical potential mu (|mu| < m), temperature T."""
    if T <= 0:
        return 0.0, 0.0, 0.0, 0.0
    x0 = (m - mu) / T
    if x0 > 700:
        return 0.0, 0.0, 0.0, 0.0
    def occ(p):
        E = math.sqrt(p * p + m * m)
        y = (E - mu) / T
        return 0.0 if y > 700 else 1.0 / math.expm1(y)
    # scale of the momentum integrand: thermal p ~ sqrt(2 m T + T^2) ... put breakpoints
    pt = math.sqrt(T * T + 2 * m * T)
    pmax = pt * 40 + 40 * T
    pts = [pt * 0.05, pt * 0.3, pt, 3 * pt, 10 * pt]
    def I(fun):
        tot = 0.0
        edges = [0.0] + [p for p in pts if p < pmax] + [pmax]
        for lo, hi in zip(edges[:-1], edges[1:]):
            tot += quad(fun, lo, hi, epsabs=0, epsrel=1e-11, limit=200)[0]
        return tot / (2 * math.pi ** 2)
    n = I(lambda p: p * p * occ(p))
    e = I(lambda p: p * p * math.sqrt(p * p + m * m) * occ(p))
    P = I(lambda p: p ** 4 / (3 * math.sqrt(p * p + m * m)) * occ(p))
    s = (e + P - mu * n) / T
    return n, P, e, s


def vapour(mu, T):
    nx, px, ex, sx = bose_dof(1.0, mu, T)
    nb, pb, eb, sb = bose_dof(1.0, -mu, T)
    nf, pf, ef, sf = bose_dof(MPHI, 0.0, T)
    return dict(n_net=nx - nb, n_X=nx, n_Xbar=nb, n_phi=nf, P=px + pb + pf, e=ex + eb + ef + mu * 0, s=sx + sb + sf)


# ------------------------------------------------------------------ bulk Q-matter
def Uphi(S):
    return b_ / 4 * (S * S - 1) ** 2 - c_ * (S - 1) ** 3


def Uphi_p(S):
    return b_ * (S * S - 1) * S - 3 * c_ * (S - 1) ** 2


def Uphi_pp(S):
    return b_ * (3 * S * S - 1) - 6 * c_ * (S - 1)


def P0(mu):
    """Zero-T pressure of homogeneous Q-matter at chemical potential mu (m_X^4), maximised over (S, F); 0 if the
    maximum is the vacuum (F = 0)."""
    f = lambda x: -(0.5 * mu * mu * x[1] ** 2 - Uphi(x[0]) - 0.5 * x[0] ** 2 * x[1] ** 2 - a_ / 4 * x[1] ** 4)
    r = minimize(f, x0=[S_B, F_B], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
    val = -r.fun * 100
    return max(val, 0.0), r.x


def liquid_branches(k):
    """Rotating-frame excitation energies (three branches) of bulk Q-matter at wavenumber k (m_X units)."""
    S, F, Om = S_B, F_B, OM_B
    D = Uphi_pp(S) + F * F
    k2 = k * k
    p1 = np.polymul([1, k2 + 2 * a_ * F * F], [1, k2])
    p = np.polymul(p1, [1, k2 + D])
    p = np.polyadd(p, 4 * Om * Om * np.polymul([1, 0], [1, k2 + D]))
    p = np.polysub(p, 4 * F * F * S * S * np.array([1, k2]))
    z = np.roots(p)
    w2 = np.sort(-z.real)
    return np.sqrt(np.maximum(w2, 0.0))


_KG = None


def _liq_table():
    global _KG
    if _KG is None:
        ks = np.concatenate([np.geomspace(1e-5, 1.0, 400), np.linspace(1.0, 400.0, 4000)[1:]])
        W = np.array([liquid_branches(k) for k in ks])
        _KG = (ks, W)
    return _KG


def liquid_thermal(T):
    """One-loop thermal pressure and entropy density of the three liquid branches (m_X^4, m_X^3)."""
    ks, W = _liq_table()
    if T <= 0:
        return 0.0, 0.0
    x = W / T
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        lnf = np.where(x > 700, 0.0, -np.log1p(-np.exp(-np.minimum(x, 700))))
        occ = np.where(x > 700, 0.0, 1.0 / np.expm1(np.maximum(np.minimum(x, 700), 1e-300)))
    p = T * np.trapezoid(ks[:, None] ** 2 * lnf, ks, axis=0).sum() / (2 * math.pi ** 2)
    e = np.trapezoid(ks[:, None] ** 2 * W * occ, ks, axis=0).sum() / (2 * math.pi ** 2)
    s = (e + p) / T
    return p, s


def mu_coex(T, thermal_liquid=True):
    """Coexistence chemical potential of bulk Q-matter with the vapour at temperature T."""
    pl = liquid_thermal(T)[0] if thermal_liquid else 0.0
    # P0(mu) ~ n_B (mu - Om_B) + 0.5 chi (mu - Om_B)^2 near Om_B; solve P0(mu) + pl = p_vap(mu) by bracketing
    def g(mu):
        return P0(mu)[0] + pl - vapour(mu, T)['P'] if mu >= OM_B else N_B * (mu - OM_B) + pl - vapour(mu, T)['P']
    lo, hi = OM_B - 0.3, 0.9999
    if g(lo) * g(hi) > 0:
        return None
    return brentq(g, lo, hi, xtol=1e-12)


# ------------------------------------------------------------------ object family
class Family:
    def __init__(self, path=None):
        d = json.load(open(path or os.path.join(RES, 'family.json')))
        rows = d['rows']
        Q = [r['Q'] for r in rows]; EQ = [r['E_over_Qm'] for r in rows]; R90 = [r['r90'] for r in rows]
        # stable branch below 6e4 (fixed Omega, Omega < Omega(Q_min))
        br = [b for b in d['branch_near_Qmin'] if b['Omega'] < d['Omega_at_Qmin'] and b['Q'] < 59999]
        for b in br:
            Q.append(b['Q']); EQ.append(b['E_over_Qm']); R90.append(b['r90'])
        Q.append(d['Q_min']); EQ.append(d['E_over_Qm_at_Qmin']); R90.append(np.interp(d['Q_min'], [b['Q'] for b in br][::-1], [b['r90'] for b in br][::-1]) if br else 7.4)
        o = np.argsort(Q)
        self.Q = np.array(Q)[o]; self.EQ = np.array(EQ)[o]; self.R90 = np.array(R90)[o]
        keep = np.r_[True, np.diff(self.Q) > 1.0]
        self.Q, self.EQ, self.R90 = self.Q[keep], self.EQ[keep], self.R90[keep]
        self.Qmin = d['Q_min']
        self.fit = d['fit3']
        self.sp = CubicSpline(np.log(self.Q), self.EQ)
        self.spr = CubicSpline(np.log(self.Q), np.log(self.R90))

    def e_over_q(self, Q):
        Q = np.asarray(Q, float)
        f = self.fit
        big = OM_B + f['A'] * Q ** (-1 / 3) + f['B'] * Q ** (-2 / 3) + f['C'] / Q
        return np.where(Q <= self.Q[-1], self.sp(np.log(np.clip(Q, self.Q[0], self.Q[-1]))), big)

    def r90(self, Q):
        Q = np.asarray(Q, float)
        bulkR = (3 * Q / (4 * math.pi * N_B)) ** (1 / 3)
        return np.where(Q <= self.Q[-1], np.exp(self.spr(np.log(np.clip(Q, self.Q[0], self.Q[-1])))), bulkR * self.R90[-1] / ((3 * self.Q[-1] / (4 * math.pi * N_B)) ** (1 / 3)))


def saha_objects(fam, mu, T, Qmax, gint=1.0, nq=20000):
    """Charge density and number density in objects Q in [Q_min, Q_max] (continuous-Q Saha sum), as logs."""
    lq = np.linspace(math.log(fam.Qmin * 1.0000001), math.log(Qmax), nq)
    Q = np.exp(lq)
    E = fam.e_over_q(Q) * Q
    f = math.log(gint) + 1.5 * np.log(E * T / (2 * math.pi)) - (E - mu * Q) / T + np.log(Q)   # log of Q n(Q)
    # integrate exp(f) dQ piecewise with linear interpolation of f in Q (exact exponential integral per cell)
    dQ = np.diff(Q); df = np.diff(f)
    fm = np.max(f)
    with np.errstate(over='ignore', invalid='ignore'):
        cell = np.where(np.abs(df) < 1e-8, dQ * np.exp(f[:-1] - fm), dQ * (np.exp(f[1:] - fm) - np.exp(f[:-1] - fm)) / df)
    lcharge = fm + math.log(np.sum(cell))
    g = f - np.log(Q)
    gm = np.max(g)
    dg = np.diff(g)
    with np.errstate(over='ignore', invalid='ignore'):
        cell2 = np.where(np.abs(dg) < 1e-8, dQ * np.exp(g[:-1] - gm), dQ * (np.exp(g[1:] - gm) - np.exp(g[:-1] - gm)) / dg)
    lnum = gm + math.log(np.sum(cell2))
    return lcharge, lnum


def saha_solve(fam, nQ, T, Qmax, gint=1.0):
    """Solve n_Q = n_vap,net(mu) + Sum_Q Q n(Q) for mu; return eta and mu."""
    def resid(mu):
        v = vapour(mu, T)['n_net']
        lc, _ = saha_objects(fam, mu, T, Qmax, gint)
        tot = v + math.exp(min(lc, 700))
        return math.log(max(tot, 1e-300)) - math.log(nQ)
    lo, hi = -0.9999, 0.9999
    if resid(hi) < 0:
        return dict(mu=None, eta=None, note='vapour+objects cannot hold the charge below mu = 0.9999 (BEC of free X)')
    mu = brentq(resid, lo, hi, xtol=1e-13)
    lc, ln = saha_objects(fam, mu, T, Qmax, gint)
    eta = math.exp(lc) / nQ
    return dict(mu=mu, eta=eta, mean_Q=math.exp(lc - ln), n_obj=math.exp(ln))
