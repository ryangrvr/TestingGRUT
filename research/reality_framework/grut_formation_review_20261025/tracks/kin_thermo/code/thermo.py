#!/usr/bin/env python3
"""kin_thermo (Claude, own code): Bose-Einstein state functions of the hidden gas WITH the light-scalar reservoir.

Species: X (charge +1, chemical potential mu_+), Xbar (charge -1, mu_-), s (neutral real scalar, mu_s), each one
bosonic degree, Bose statistics.  Units m_X = 1 (t = T/m_X, u = mu/m_X, r = m/m_X; densities in m_X^3,
energy densities in m_X^4).  For one species of mass r at (t, u), u < r (Kolb & Turner; the same modified-Bessel
series our ls_threshold/lsthermo.py uses, re-implemented here with the derivatives and the tail control):

  n_k   = (r^2 t/2pi^2) e^{k u/t} K2(k r/t)/k
  P_k   = (r^2 t^2/2pi^2) e^{k u/t} K2(k r/t)/k^2
  rho_k = (r^2 t/2pi^2) e^{k u/t} [3 t K2/k^2 + r K1/k]
  E2_k  = (r^2/2pi^2) e^{k u/t} [12 t^3 K2/k^3 + 5 r t^2 K1/k^2 + r^2 t K0/k]      (int p^2 E^2 f_k /2pi^2)
  n = sum n_k, ...; f(1+f) = sum_k k e^{-k(E-u)/t}  gives the exact susceptibilities
  dn/du|_t = (1/t) sum k n_k,     dn/dt|_u = (1/t^2) sum k (rho_k - u n_k),
  drho/du|_t = (1/t) sum k rho_k, drho/dt|_u = (1/t^2) sum k (E2_k - u rho_k),   dP/du = n, dP/dt = s.
Entropy density s = (rho + P - sum_i u_i n_i)/t.  Scaled Bessel functions kve are used, so every term carries
e^{-k(r-u)/t}; the series is cut when that factor is below e^{-40} (K up to 1.5e6 terms, cached per (t, r)).
Within |r - u|/t < 2.7e-5 (K would exceed the cap) the direct momentum quadrature `quad_species` is used instead.

Closures offered (chemical equilibrium relations):
  pairs:    mu_+ + mu_- = 0       (4<->2 and ss<->X Xbar fast: pair number equilibrated)
  fast s:   mu_s = mu_+ + mu_-    (s <-> X Xbar fast: s in equilibrium with the pair gas; = 0 in full equilibrium)
  caps:     mu_+ <= u_cap (objects at omega), mu_- <= u_cap (anti-objects), mu_s <= r_s (s condensation).
Nothing here chooses the closure; the solver does.  Codex's THERMODYNAMIC_CONTROL state uses pairs + fast s.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import kve

TWO_PI2 = 2 * np.pi**2
RS = np.sqrt(6.0)              # m_s/m_X
KMAX = 1_500_000
_cache = {}


def _bessel(t, r, umax):
    """k, K0e, K1e, K2e arrays long enough for every u <= umax at temperature t (scaled by e^{k r/t})."""
    d = max((r - umax) / t, 1e-12)
    K = int(min(40.0 / d + 60, KMAX))
    key = (round(float(t), 15), round(float(r), 15), K)
    if key not in _cache:
        if len(_cache) > 48:
            _cache.clear()
        k = np.arange(1, K + 1, dtype=float)
        z = k * r / t
        _cache[key] = (k, kve(0, z), kve(1, z), kve(2, z))
    return _cache[key]


def _needs_quad(t, r, u):
    return (r - u) / t < 40.0 / KMAX * 1.0 + 1e-15


def quad_species(t, r, u, derivs=False):
    """Direct momentum quadrature (independent of the series): n, P, rho[, dn/du, drho/du, dn/dt, drho/dt]."""
    def f_occ(x):               # x = p/t; E - u computed without cancellation for u near r
        E = np.sqrt(x * x + r * r / (t * t)) * t
        w = (x * x * t * t / (E + r) + (r - u)) / t
        if w > 700:
            return 0.0, 0.0, E
        f = 1.0 / np.expm1(w) if w > 1e-12 else 1.0 / w
        return f, f * (1 + f), E

    def integ(kind):
        def g(x):
            f, ff, E = f_occ(x)
            if kind == 'n':
                return x * x * f
            if kind == 'P':
                return x * x * f * (x * x * t * t) / (3 * E)
            if kind == 'rho':
                return x * x * f * E
            if kind == 'dn_du':
                return x * x * ff / t
            if kind == 'drho_du':
                return x * x * ff * E / t
            if kind == 'dn_dt':
                return x * x * ff * (E - u) / (t * t)
            if kind == 'drho_dt':
                return x * x * ff * E * (E - u) / (t * t)
        pts = [0.0, 0.03, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 400.0]
        return sum(quad(g, a, b, epsabs=0, epsrel=1e-11, limit=400)[0] for a, b in zip(pts[:-1], pts[1:]))

    n = t**3 / TWO_PI2 * integ('n')
    P = t**3 / TWO_PI2 * integ('P')
    rho = t**3 / TWO_PI2 * integ('rho')
    out = dict(n=n, P=P, rho=rho)
    if derivs:
        out.update(dn_du=t**3 / TWO_PI2 * integ('dn_du'), drho_du=t**3 / TWO_PI2 * integ('drho_du'),
                   dn_dt=t**3 / TWO_PI2 * integ('dn_dt'), drho_dt=t**3 / TWO_PI2 * integ('drho_dt'))
    return out


def species(t, r, u, derivs=True):
    """n, P, rho, s and the four susceptibilities for one bosonic degree of mass r at (t, u); dimensionless."""
    if u >= r:
        raise ValueError('chemical potential at or above the mass: u=%r r=%r' % (u, r))
    if _needs_quad(t, r, u):
        q = quad_species(t, r, u, derivs=True)
        q['s'] = (q['rho'] + q['P'] - u * q['n']) / t
        q['method'] = 'quad'
        return q
    k, k0, k1, k2 = _bessel(t, r, u)
    e = np.exp(-k * (r - u) / t)
    a2 = k2 * e
    a1 = k1 * e
    a0 = k0 * e
    pref = r * r / TWO_PI2
    n_k = pref * t * a2 / k
    P_k = pref * t * t * a2 / (k * k)
    rho_k = pref * t * (3 * t * a2 / (k * k) + r * a1 / k)
    n, P, rho = n_k.sum(), P_k.sum(), rho_k.sum()
    out = dict(n=n, P=P, rho=rho, s=(rho + P - u * n) / t, method='series', terms=len(k))
    if derivs:
        E2_k = pref * (12 * t**3 * a2 / k**3 + 5 * r * t * t * a1 / (k * k) + r * r * t * a0 / k)
        out.update(dn_du=np.dot(k, n_k) / t, drho_du=np.dot(k, rho_k) / t,
                   dn_dt=np.dot(k, rho_k - u * n_k) / (t * t), drho_dt=np.dot(k, E2_k - u * rho_k) / (t * t))
    return out


def log_n(t, r, u):
    """ln n(t, r, u), stable for very negative u (e^{-(r-u)/t} factored out)."""
    k, _, _, k2 = _bessel(t, r, min(u, r - 1e-9))
    e = np.exp(-(k - 1) * (r - u) / t)
    return -(r - u) / t + np.log(r * r * t / TWO_PI2 * np.dot(k2 * e, 1 / k))


def net_charge(t, u):
    """n(t,1,u) - n(t,1,-u) for the X/Xbar pair at +-u without cancellation (pairs closure)."""
    k, _, _, k2 = _bessel(t, 1.0, abs(u))
    d = np.exp(-k * (1.0 - abs(u)) / t) * (-np.expm1(-2 * k * abs(u) / t))
    return np.sign(u) * t / TWO_PI2 * np.dot(k2 * d, 1 / k)


# ---------------------------------------------------------------- gas state at given (t, u_+, u_-, u_s)
def gas(t, up, um, us=None, fast_s=True, with_s=True, ms_over_m=RS, mx_over_m=1.0):
    """Full gas state.  fast_s: u_s = u_+ + u_-.  Returns per-species dicts and totals (dimensionless)."""
    if us is None:
        us = (up + um) if fast_s else 0.0
    X = species(t, mx_over_m, up)
    Xb = species(t, mx_over_m, um)
    S = species(t, ms_over_m, us) if with_s else dict(n=0., P=0., rho=0., s=0., dn_du=0., drho_du=0., dn_dt=0., drho_dt=0.)
    rho = X['rho'] + Xb['rho'] + S['rho']
    P = X['P'] + Xb['P'] + S['P']
    nq = X['n'] - Xb['n']
    if up == -um:
        nq = net_charge(t, up) * (1.0 if mx_over_m == 1.0 else 1.0)
    ent = (rho + P - up * X['n'] - um * Xb['n'] - us * S['n']) / t
    return dict(t=t, up=up, um=um, us=us, X=X, Xb=Xb, S=S, n_plus=X['n'], n_minus=Xb['n'], n_s=S['n'],
                n_Q=nq, pair_inventory=Xb['n'] + S['n'], rho=rho, P=P, s=ent,
                scalar_energy_fraction=S['rho'] / rho if rho > 0 else 0.0)


def susceptibility(t, up, um, us=None, fast_s=True, with_s=True, ms_over_m=RS):
    """Jacobian d(n_Q, n_minus, n_s, rho, P_inv) / d(t, u_+, u_-, u_s) with P_inv = n_- + n_s (pair inventory).
    With fast_s the u_s column is folded into u_+ and u_- (u_s = u_+ + u_-) and a 3x3 block is also returned."""
    g = gas(t, up, um, us, fast_s, with_s, ms_over_m)
    X, Xb, S = g['X'], g['Xb'], g['S']
    # rows: n_Q, n_-, n_s, rho ; columns: t, u_+, u_-, u_s
    J = np.array([[X['dn_dt'] - Xb['dn_dt'], X['dn_du'], -Xb['dn_du'], 0.0],
                  [Xb['dn_dt'], 0.0, Xb['dn_du'], 0.0],
                  [S['dn_dt'], 0.0, 0.0, S['dn_du']],
                  [X['drho_dt'] + Xb['drho_dt'] + S['drho_dt'], X['drho_du'], Xb['drho_du'], S['drho_du']]])
    out = dict(J4=J, rows=['n_Q', 'n_minus', 'n_s', 'rho'], cols=['t', 'u_plus', 'u_minus', 'u_s'], state=g)
    if fast_s:
        M = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1]], float)   # (t,u+,u-,u_s) <- (t,u+,u-)
        J3 = J @ M
        out['J3'] = J3
        out['cols3'] = ['t', 'u_plus', 'u_minus']
        # pair-inventory susceptibility at fixed t and n_Q (compare MATHEMATICS.md: 1 + K(n_+ + n_-) for MB)
        # dP_inv = (row n_- + row n_s) . d(t,u+,u-);  constrain dt = 0, dn_Q = 0 -> du_+ = -(J3[0,2]/J3[0,1]) du_-
        a = J3[0, 1]; b = J3[0, 2]
        dpinv = (J3[1] + J3[2])
        dnm = J3[1]
        ratio_u = -b / a
        dP_dnm_fixed_t_nQ = (dpinv[1] * ratio_u + dpinv[2]) / (dnm[1] * ratio_u + dnm[2])
        out['dPinv_dnminus_fixed_t_nQ'] = dP_dnm_fixed_t_nQ
    return out


def K_bose(t, ms_over_m=RS):
    """n_s^eq/(n_+^eq n_-^eq) at zero chemical potentials with Bose statistics (m_X^-3)."""
    S = species(t, ms_over_m, 0.0, derivs=False)
    X = species(t, 1.0, 0.0, derivs=False)
    return S['n'] / (X['n'] ** 2)


def K_mb(t, ms_over_m=RS):
    """Maxwell-Boltzmann K(t) as used by lsrelic.K_s (relativistic MB densities)."""
    r = ms_over_m
    return TWO_PI2 * r * r * kve(2, r / t) / (t * kve(2, 1 / t)**2) * np.exp(-(r - 2) / t)


# ---------------------------------------------------------------- inverse maps (state solvers)
def _fill(t, r, target, cap):
    """u such that n(t, r, u) = target for u <= cap; returns (u, condensed excess).  Newton on ln n (convex)."""
    if target <= 0:
        return -np.inf, 0.0
    lt = np.log(target)
    lcap = log_n(t, r, cap)
    if lcap <= lt:
        return cap, target - np.exp(lcap)
    # bracket and brentq on ln n (monotone in u)
    lo = cap - 1.0
    while log_n(t, r, lo) > lt:
        lo -= max(1.0, 2 * t)
        if lo < -1e6:
            break
    u = brentq(lambda x: log_n(t, r, x) - lt, lo, cap, xtol=1e-15, rtol=1e-14, maxiter=500)
    return u, 0.0


def state_equilibrium(nq, rho, u_cap, e_cap=None, with_s=True, include_cond_energy=True, t_bounds=(1e-4, 3e3)):
    """FULL chemical equilibrium (pairs: u_- = -u_+; fast s: u_s = 0): given net charge nq and total energy rho
    (both dimensionless, condensate energy e_cap per charge included), return t, u, condensed charge n_c, state.
    This is the closure of Codex's THERMODYNAMIC_CONTROL / the old proxies (linear object energy e_cap*n_c)."""
    e_cap = u_cap if e_cap is None else e_cap

    def at(t):
        if net_charge(t, u_cap) <= nq:
            u, nc = u_cap, nq - net_charge(t, u_cap)
        else:
            u = brentq(lambda x: net_charge(t, x) - nq, 0.0, u_cap, xtol=1e-15, rtol=1e-14)
            nc = 0.0
        g = gas(t, u, -u, 0.0, with_s=with_s)
        return u, nc, g

    def resid(lt):
        u, nc, g = at(np.exp(lt))
        return g['rho'] + (e_cap * nc if include_cond_energy else 0.0) - rho

    lt = brentq(resid, np.log(t_bounds[0]), np.log(t_bounds[1]), xtol=1e-13, rtol=1e-13)
    t = np.exp(lt)
    u, nc, g = at(t)
    g.update(n_cond=nc, u_cap=u_cap, e_cap=e_cap, rho_total=g['rho'] + e_cap * nc, closure='pairs+fast_s')
    return g


def state_two_potential(nq, pinv, rho, u_cap, e_cap=None, anti_cap=None, with_s=True, fast_s=True,
                        t_bounds=(1e-4, 3e3), tguess=None):
    """TWO chemical potentials (kinetic equilibrium, fast s): given net charge nq, pair inventory
    pinv = n_- + n_s (+ anti-condensate) and total energy rho, return t, u_+, u_-, condensates and the state.
    Objects at u_+ = u_cap take the X excess (energy e_cap per charge); anti-objects at u_- = anti_cap take the
    Xbar excess.  With fast_s, u_s = u_+ + u_- and n_s follows; the inventory equation is solved for u_-."""
    e_cap = u_cap if e_cap is None else e_cap
    anti_cap = u_cap if anti_cap is None else anti_cap

    def at(t):
        # unknown u_-: n_-(u_-) + n_s(u_+ + u_-) [+ anti-cond] = pinv, with u_+ from n_+ = nq + n_- (+cond)
        def inv(um):
            nm = species(t, 1.0, um, derivs=False)['n'] if um > -np.inf else 0.0
            up, cp = _fill(t, 1.0, nq + nm, u_cap)
            ns = species(t, RS, up + um, derivs=False)['n'] if (with_s and fast_s) else 0.0
            return nm + ns, up, cp, nm, ns
        # anti-cap check: inventory at u_- = anti_cap
        top, up_top, cp_top, nm_top, ns_top = inv(anti_cap)
        if top <= pinv:
            um, cm = anti_cap, pinv - top
            up, cp, nm, ns = up_top, cp_top, nm_top, ns_top
        else:
            lo = anti_cap - 1.0
            while inv(lo)[0] > pinv:
                lo -= max(1.0, 2 * t)
                if lo < -1e6:
                    break
            um = brentq(lambda x: inv(x)[0] - pinv, lo, anti_cap, xtol=1e-15, rtol=1e-14, maxiter=500)
            _, up, cp, nm, ns = inv(um)
            cm = 0.0
        g = gas(t, up, um, (up + um) if (with_s and fast_s) else 0.0, fast_s=fast_s, with_s=with_s)
        return g, cp, cm

    def resid(lt):
        g, cp, cm = at(np.exp(lt))
        return np.log((g['rho'] + e_cap * cp + anti_cap * cm) / rho)

    lo, hi = np.log(t_bounds[0]), np.log(t_bounds[1])
    if tguess is not None:
        a_, b_ = np.log(tguess) - 0.1, np.log(tguess) + 0.1
        ra, rb = resid(a_), resid(b_)
        while ra > 0 and a_ > lo:
            b_, rb = a_, ra; a_ -= 0.4; ra = resid(a_)
        while rb < 0 and b_ < hi:
            a_, ra = b_, rb; b_ += 0.4; rb = resid(b_)
        if ra <= 0 <= rb:
            lo, hi = a_, b_
    lt = brentq(resid, lo, hi, xtol=1e-13, rtol=1e-13)
    t = np.exp(lt)
    g, cp, cm = at(t)
    g.update(n_cond=cp, n_anticond=cm, u_cap=u_cap, e_cap=e_cap, anti_cap=anti_cap,
             rho_total=g['rho'] + e_cap * cp + anti_cap * cm, closure='two-potential' + ('+fast_s' if fast_s else ''))
    return g


def state_three_potential(nq, nm, ns, rho, t0=5.0, up0=0.0, um0=0.0, us0=0.0, tol=1e-12, itmax=80):
    """THREE chemical potentials (s chemically decoupled from the pairs): given n_Q, n_-, n_s and rho, Newton
    iteration in (ln t, u_+, u_-, u_s) with the analytic Jacobian.  No caps (raises if a potential reaches
    its mass); the caller handles condensation."""
    x = np.array([np.log(t0), up0, um0, us0], float)
    target = np.array([nq, nm, ns, rho])
    for it in range(itmax):
        t, up, um, us = np.exp(x[0]), x[1], x[2], x[3]
        if up >= 1 or um >= 1 or us >= RS:
            raise ValueError('potential reached a mass during Newton; condensation needed')
        S = susceptibility(t, up, um, us, fast_s=False)
        g = S['state']
        F = np.array([g['n_Q'], g['n_minus'], g['n_s'], g['rho']]) - target
        J = S['J4'].copy()
        J[:, 0] *= t                     # d/d ln t
        dx = np.linalg.solve(J, -F)
        lam = 1.0
        for _ in range(30):               # damped step keeping potentials below masses
            xn = x + lam * dx
            if xn[1] < 1 and xn[2] < 1 and xn[3] < RS:
                break
            lam *= 0.5
        x = xn
        if np.max(np.abs(lam * dx)) < tol:
            break
    t, up, um, us = np.exp(x[0]), x[1], x[2], x[3]
    g = gas(t, up, um, us, fast_s=False)
    g.update(closure='three-potential', newton_iterations=it + 1,
             residual=float(np.max(np.abs(np.array([g['n_Q'], g['n_minus'], g['n_s'], g['rho']]) - target) / np.abs(target))))
    return g


if __name__ == '__main__':
    import json
    MX = 0.09240896373753621
    # Codex's reservoir control rows (mu = 0)
    for t in (0.3, 1.0, 3.0, 12.0, 30.0):
        g = gas(t, 0.0, 0.0, 0.0)
        print('t=%5.1f  n_s/n_- = %.12f  scalar fraction = %.12f' % (t, g['n_s'] / g['n_minus'], g['scalar_energy_fraction']))
    # matched-energy states at R = 40.5626
    import kt_common as C
    rho = C.rho_onset(40.56262288525937) / MX**4
    nq = C.N_Q0 / MX**3
    for with_s in (False, True):
        g = state_equilibrium(nq, rho, C.OMEGA_REF, with_s=with_s)
        print('with_s', with_s, 't=%.9f u=%.6f n_cond=%.9f scalar=%.9f pinv=%.9f' % (
            g['t'], g['up'], g['n_cond'] * MX**3, g['scalar_energy_fraction'], g['pair_inventory'] * MX**3))
