"""Track `bound`: full flat-FRW integrations of the q-theory (action of BRIEF.md, fixed mu) with dust + radiation,
for explicit (eps, K) families, in the three regimes of NOTES section A. Writes results/frw_regimes.json.

Units: H0 ~ 1 (time in 1/H0), K(q0) = 1 at the Minkowski point of mu0. Compressed hierarchy m_q/H0 = 0.1 ... 100
(physically m_q/H0 >= 2.5e30 for an unscreened q-mode, NOTES A.4); the runs test the RELATION, not the physical numbers.

Equations (equations track 3.2, V_m = 0):  y = (q, v = qdot, H, rho_d, rho_r)
  R_q   = 2 (eps'(q) - mu)/K'(q)
  vdot  = { [rho_d + 4 rho_V - K R_q]/3 - 3 H K' v - K'' v^2 } / K'          rho_V = eps - mu q
  Hdot  = R_q/6 - 2 H^2 ;  rho_d' = -3 H rho_d ;  rho_r' = -4 H rho_r
  constraint C = 3 H^2 K + 3 H K' v - rho_d - rho_r - rho_V = 0  (propagates as Cdot = -4 H C)

Part I  (regime i, adiabatic m >> H): mu = mu0 + dmu calibrated so that the static de Sitter root gives Omega_Lam ~ 0.685;
        quasi-static initial data at H = H_start; compare
            Delta(t) = [3 K_lab H^2 - rho_d - rho_r] - Lam_cal
        with the derived  rho_dyn = [(9/4) rho_d^2 + 2 rho_d rho_r - 4 Lam rho_r]/(3 K_lab m^2)   (adiabatic_relation.py S1).
Part II (regime ii, Hubble-frozen m <~ H0): mu = mu0 (no calibration offset), q displaced at early times by dq_i, shot so that
        H = 1 when rho_d = 0.945; w_eff(z), CPL (w0, wa), K drift (G variation), Gdot/G today.
Part III (regime iii, oscillation relic): m = 10, 30 H0 with a displacement frozen until H ~ m; dust-like scaling of the relic.
Run: python3 frw_regimes.py   (~1-3 min)
"""
import json
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OMEGA_M, OMEGA_R = 0.315, 9.1e-5           # Planck 2018 Omega_m [Q via ledger row 9]; Omega_r h^2-based value [R], ~9e-5
OMEGA_L = 1.0 - OMEGA_M - OMEGA_R
H0_PER_YR = 67.4 / 3.0856775814913673e19 * 3.15576e7   # H0 = 67.4 km/s/Mpc in 1/yr


# ----------------------------------------------------------------------------------------------- families
class Family:
    """eps(q), K(q) with derivatives; Minkowski point q0 = 1, K(q0) = 1, m^2(q0) = (2/3) K eps''/K'^2 = m_target^2."""

    def __init__(self, name, m):
        self.name, self.m = name, m
        if name == 'F1_quadratic_linearK':          # eps = (e2/2)(1 + q^2): q0 = 1, mu0 = e2 ; K = 1 + k1 (q - 1)
            self.k1 = 0.4
            self.e2 = 1.5 * m**2 * self.k1**2
            self.mu0 = self.e2
        elif name == 'F2_KV0806_quartic_Kpropq':    # KV 0806.2805 (4.7),(4.14): eps = (E/2)(-q^2 + q^4/3), K = q
            self.E = 1.5 * m**2
            self.mu0 = -self.E / 3.0
        elif name == 'F3_cubic_quartic_expK':       # eps = e2 x^2/2 + e3 x^3/6 + e4 x^4/24 + mu0 q (x = q-1), K = exp(b x)
            self.b = 0.7
            self.e2 = 1.5 * m**2 * self.b**2
            self.e3 = 3.0 * self.e2
            self.e4 = 10.0 * self.e2
            self.mu0 = 0.3 * self.e2                # any mu0: eps(q0) - mu0 q0 = 0 holds by construction
        else:
            raise ValueError(name)

    def eps(self, q):
        if self.name.startswith('F1'):
            return 0.5 * self.e2 * (1.0 + q * q)
        if self.name.startswith('F2'):
            return 0.5 * self.E * (-q * q + q**4 / 3.0)
        x = q - 1.0
        return self.e2 * x * x / 2 + self.e3 * x**3 / 6 + self.e4 * x**4 / 24 + self.mu0 * q

    def eps1(self, q):
        if self.name.startswith('F1'):
            return self.e2 * q
        if self.name.startswith('F2'):
            return self.E * (-q + 2.0 * q**3 / 3.0)
        x = q - 1.0
        return self.e2 * x + self.e3 * x * x / 2 + self.e4 * x**3 / 6 + self.mu0

    def eps2(self, q):
        if self.name.startswith('F1'):
            return self.e2
        if self.name.startswith('F2'):
            return self.E * (-1.0 + 2.0 * q * q)
        x = q - 1.0
        return self.e2 + self.e3 * x + self.e4 * x * x / 2

    def K(self, q):
        if self.name.startswith('F1'):
            return 1.0 + self.k1 * (q - 1.0)
        if self.name.startswith('F2'):
            return q
        return np.exp(self.b * (q - 1.0))

    def K1(self, q):
        if self.name.startswith('F1'):
            return self.k1
        if self.name.startswith('F2'):
            return 1.0
        return self.b * np.exp(self.b * (q - 1.0))

    def K2(self, q):
        if self.name.startswith('F1') or self.name.startswith('F2'):
            return 0.0
        return self.b**2 * np.exp(self.b * (q - 1.0))


def m2_at(fam, q, R=0.0):
    """scalaron mass^2 at the point (q, R): (f' - R f'')/(3 f''), f' = K, f'' = K'^2/(2 (eps'' - K'' R/2))."""
    D = fam.eps2(q) - fam.K2(q) * R / 2
    fpp = fam.K1(q)**2 / (2 * D)
    return (fam.K(q) - R * fpp) / (3 * fpp)


def bracket_root(f, guess, step=0.05, maxw=0.95):
    """root of f nearest to guess (expanding symmetric bracket, kept inside (guess - maxw, guess + maxw))."""
    w = step
    while w <= maxw + 1e-12:
        lo, hi = guess - w, guess + w
        if np.sign(f(lo)) != np.sign(f(hi)):
            return brentq(f, lo, hi, xtol=1e-15)
        w *= 1.5
    raise ValueError('no root bracketed near %g' % guess)


def q_R0(fam, mu, guess=1.0):
    """R = 0 point of this mu: eps'(q) = mu."""
    return bracket_root(lambda q: fam.eps1(q) - mu, guess)


def static_root(fam, mu, guess=1.0):
    """de Sitter / AdS static point: eps' - mu = 2 K' rho_V / K (root with K > 0 nearest to guess)."""
    f = lambda q: fam.eps1(q) - mu - 2 * fam.K1(q) * (fam.eps(q) - mu * q) / fam.K(q)
    try:
        qs = bracket_root(f, guess, step=0.02, maxw=0.95)
    except ValueError:
        grid = np.linspace(guess - 60, guess + 60, 240001)
        grid = grid[np.array([fam.K(g) for g in grid]) > 1e-6]
        vals = np.array([f(g) for g in grid])
        idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
        if len(idx) == 0:
            raise
        j = idx[np.argmin(np.abs(grid[idx] - guess))]
        qs = brentq(f, grid[j], grid[j + 1], xtol=1e-15)
    rhoV = fam.eps(qs) - mu * qs
    return qs, rhoV, rhoV / (3 * fam.K(qs))


def rhs_factory(fam, mu):
    def rhs(t, y):
        q, v, H, rd, rr = y
        K, k1, k2 = fam.K(q), fam.K1(q), fam.K2(q)
        rhoV = fam.eps(q) - mu * q
        Rq = 2.0 * (fam.eps1(q) - mu) / k1
        Kdd = (rd + 4.0 * rhoV - K * Rq) / 3.0 - 3.0 * H * k1 * v
        vdot = (Kdd - k2 * v * v) / k1
        Hdot = Rq / 6.0 - 2.0 * H * H
        return [v, vdot, Hdot, -3.0 * H * rd, -4.0 * H * rr]
    return rhs


def constraint_rel(fam, mu, y):
    q, v, H, rd, rr = y
    rhoV = fam.eps(q) - mu * q
    C = 3 * H * H * fam.K(q) + 3 * H * fam.K1(q) * v - rd - rr - rhoV
    return C / (rd + rr + abs(rhoV) + 3 * H * H * fam.K(q))


def q_quasistatic(fam, mu, rd, guess):
    f = lambda q: rd + 4 * (fam.eps(q) - mu * q) - fam.K(q) * 2 * (fam.eps1(q) - mu) / fam.K1(q)
    return bracket_root(f, guess, step=0.02)


def H_from_constraint(fam, mu, q, v, rd, rr):
    K, k1 = fam.K(q), fam.K1(q)
    rho = rd + rr + fam.eps(q) - mu * q
    return (-3 * k1 * v + np.sqrt(9 * k1**2 * v**2 + 12 * K * rho)) / (6 * K)


# ----------------------------------------------------------------------------------------------- Part I
def calibrate_mu(fam, target_frac=OMEGA_L):
    """dmu such that Lam_cal = 3 K_lab H_s^2 = 3 K_lab * target_frac (H0 = 1 units, K_lab = K at R = 0 of this mu)."""
    def g(dmu):
        mu = fam.mu0 + dmu
        qlab = q_R0(fam, mu)
        qs, rhoV, H2 = static_root(fam, mu, qlab)
        return 3 * fam.K(qlab) * H2 - 3 * fam.K(qlab) * target_frac
    lo = -0.1
    while np.sign(g(lo)) == np.sign(g(-1e-9)):
        lo *= 1.5
    dmu = brentq(g, lo, -1e-9, xtol=1e-14)
    return dmu


def adiabatic_state(fam, mu, rd, rd0, rr0, guess):
    """Second-order adiabatic data (q, v) on the slow manifold at dust density rd.
    q_ad = q_qs + dq2,  dq2 = (Kdd + 3 H Kd)_qs / (dL/dq),  L(q) = [rd + 4 rho_V - K R_q]/3 (trace equation);
    v_ad = (d q_ad / d ln rd) * (-3 H). Mismatch with the true slow solution is O((H/m)^4) in q."""
    def qs_of_s(s):
        return q_quasistatic(fam, mu, np.exp(s), guess)

    def H_of_s(s, q, v=0.0):
        r = np.exp(s)
        rr = rr0 * (r / rd0)**(4.0 / 3.0)
        return H_from_constraint(fam, mu, q, v, r, rr)

    def Lq(q, r):
        f = lambda qq: (r + 4 * (fam.eps(qq) - mu * qq) - fam.K(qq) * 2 * (fam.eps1(qq) - mu) / fam.K1(qq)) / 3.0
        h = 1e-6
        return (f(q + h) - f(q - h)) / (2 * h)

    def q_ad(s):
        hs = 1e-3
        qm, q0_, qp = qs_of_s(s - hs), qs_of_s(s), qs_of_s(s + hs)
        Km, K0_, Kp = fam.K(qm), fam.K(q0_), fam.K(qp)
        Ks = (Kp - Km) / (2 * hs)
        Kss = (Kp - 2 * K0_ + Km) / hs**2
        H = H_of_s(s, q0_)
        Rq = 2 * (fam.eps1(q0_) - mu) / fam.K1(q0_)
        Hdot = Rq / 6 - 2 * H * H
        sdot, sddot = -3 * H, -3 * Hdot
        Kd = Ks * sdot
        Kdd = Kss * sdot**2 + Ks * sddot
        return q0_ + (Kdd + 3 * H * Kd) / Lq(q0_, np.exp(s)), H

    s = np.log(rd)
    q, H = q_ad(s)
    hs = 1e-3
    dq_ds = (q_ad(s + hs)[0] - q_ad(s - hs)[0]) / (2 * hs)
    v = dq_ds * (-3 * H)
    for _ in range(6):
        v = dq_ds * (-3 * H)
        H = H_of_s(s, q, v)          # H last: the constraint holds exactly for the returned (q, v, H)
    return q, v, H


def run_adiabatic(fam_name, m, a_evals=(0.5, 0.7, 1.0, 1.3), start_factor=1.5):
    """For each evaluation epoch a_e: start at a_e/start_factor on the slow manifold (second-order adiabatic data),
    integrate the full equations, and compare the one-period average of Delta with that of rho_dyn at a_e.
    (A residual free oscillation enters Delta linearly; its relative size grows like (a/a_start)^3, hence the restarts.)"""
    fam = Family(fam_name, m)
    dmu = calibrate_mu(fam)
    mu = fam.mu0 + dmu
    qlab = q_R0(fam, mu)
    Klab = fam.K(qlab)
    qs, rhoV_s, H2_s = static_root(fam, mu, qlab)
    Lam_cal = 3 * Klab * H2_s
    m0sq = m2_at(fam, qlab, 0.0)
    rd0, rr0 = 3 * OMEGA_M * Klab, 3 * OMEGA_R * Klab        # values at a = 1 (defined by rho_d)
    period = 2 * np.pi / np.sqrt(m0sq)
    res = {'family': fam_name, 'm_target': m, 'm0 (R=0 of this mu)': float(np.sqrt(m0sq)), 'dmu (calibrated)': dmu,
           'K_lab': Klab, 'Lam_cal': Lam_cal, 'q_static': qs, 'H_static': float(np.sqrt(H2_s)),
           'initial_data': 'second-order adiabatic (slow-manifold) data at a_e/%.2f, one run per a_e' % start_factor}
    comp = []
    for a_e in a_evals:
        ai = a_e / start_factor
        rdi, rri = rd0 * ai**-3, rr0 * ai**-4
        qi, vi, Hi = adiabatic_state(fam, mu, rdi, rd0, rr0, qlab)
        rd_e = rd0 * a_e**-3

        def ev_end(t, y):
            return y[3] - rd_e * 0.6
        ev_end.terminal = True
        sol = solve_ivp(rhs_factory(fam, mu), (0, 200.0), [qi, vi, Hi, rdi, rri], method='DOP853', rtol=1e-12,
                        atol=1e-15, max_step=period / 40, events=[ev_end], dense_output=True)
        t = sol.t
        rd = sol.y[3]
        tc = float(np.interp(-rd_e, -rd, t))
        tt = np.linspace(tc - period / 2, tc + period / 2, 401)
        yy = sol.sol(tt)
        D = 3 * Klab * yy[2]**2 - yy[3] - yy[4] - Lam_cal
        P = (2.25 * yy[3]**2 + 2 * yy[3] * yy[4] - 4 * Lam_cal * yy[4]) / (3 * Klab * m0sq)
        Hc = float(np.mean(yy[2]))
        cons = max(abs(constraint_rel(fam, mu, sol.y[:, i])) for i in range(0, len(t), max(1, len(t) // 200)))
        comp.append({'a': a_e, 'a_start': ai, 'H/m0 at a': Hc / np.sqrt(m0sq), 'H/m0 at start': float(Hi / np.sqrt(m0sq)),
                     'ratio <Delta>/<rho_dyn>': float(np.mean(D) / np.mean(P)),
                     'free-oscillation amplitude in Delta / <rho_dyn>': float(np.std(D - P) * np.sqrt(2) / np.mean(P)),
                     'Omega_dyn_measured': float(np.mean(D) / (3 * Klab * Hc**2)),
                     'Omega_dyn_formula': float(np.mean(P) / (3 * Klab * Hc**2)),
                     'max_rel_constraint': float(cons)})
    res['comparison'] = comp
    c1 = [c for c in comp if abs(c['a'] - 1.0) < 1e-9]
    if c1:
        Om_run = rd0 / (3 * Klab * (c1[0]['H/m0 at a'] * np.sqrt(m0sq))**2)
        res['Omega_dyn_today_measured'] = c1[0]['Omega_dyn_measured']
        res['Omega_dyn_today_formula'] = c1[0]['Omega_dyn_formula']
        res['(9/4) Om_run^2 (H_run/m0)^2'] = 2.25 * Om_run**2 * c1[0]['H/m0 at a']**2
        res['Om_run'] = Om_run
    return res


# ----------------------------------------------------------------------------------------------- Part II / III
def run_light_scan(m, fam_name='F1_quadratic_linearK', a_i=1e-4, a_end=1.0):
    """Regime (ii)/(i) crossover: mu calibrated (static root gives Omega_Lam ~ 0.685), q starts at its vacuum value
    (R = 0 point of this mu, K = K_lab) with v = 0 at a_i (radiation era), dust + radiation, integrate to a = 1.
    Reports the K (i.e. 1/G) drift, Gdot/G today, the CPL projection of w_eff and Omega_m in the K(today) normalisation."""
    fam = Family(fam_name, m)
    # calibration in the de Sitter normalisation: H_s^2 = rho_V,s/(3 K(q_s)) = Omega_L (for light q the R = 0 point
    # of mu is not the local vacuum; the local q follows the cosmological one when m r << 1)
    def g(dmu):
        qs_, rV_, H2_ = static_root(fam, fam.mu0 + dmu, 1.0)
        return H2_ - OMEGA_L
    grid = -np.logspace(-4, 1.3, 400)
    vals = []
    for d in grid:
        try:
            vals.append(g(d))
        except ValueError:
            vals.append(np.nan)
    vals = np.array(vals)
    ok = np.where(np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (np.sign(vals[:-1]) != np.sign(vals[1:])))[0]
    dmu = brentq(g, grid[ok[0]], grid[ok[0] + 1], xtol=1e-14)
    mu = fam.mu0 + dmu
    qlab, _, _ = static_root(fam, mu, 1.0)          # start at the dS point (late-time attractor)
    Klab = fam.K(qlab)
    rd0, rr0 = 3 * OMEGA_M * Klab, 3 * OMEGA_R * Klab
    rdi, rri = rd0 * a_i**-3, rr0 * a_i**-4
    Hi = H_from_constraint(fam, mu, qlab, 0.0, rdi, rri)

    def ev(t, y):
        return y[3] - rd0 * a_end**-3
    ev.terminal = True
    per = 2 * np.pi / m
    sol = solve_ivp(rhs_factory(fam, mu), (0, 200.0), [qlab, 0.0, Hi, rdi, rri], method='DOP853', rtol=1e-11,
                    atol=1e-14, max_step=min(per / 30, 0.02), events=[ev], dense_output=True)
    q, v, H, rd, rr = sol.y
    a = (rd0 / rd)**(1.0 / 3.0)
    N = np.log(a)
    Kq = fam.K(q)
    Know = Kq[-1]
    rho_DE = 3 * Know * H**2 - rd - rr
    Ng = np.linspace(np.log(0.25), 0.0, 600)
    rDE = np.interp(Ng, N, rho_DE)
    w = -1 - np.gradient(np.log(np.abs(rDE)), Ng) / 3
    ag = np.exp(Ng)
    sel = ag >= 0.4
    A = np.vstack([np.ones(sel.sum()), 1 - ag[sel]]).T
    w0, wa = np.linalg.lstsq(A, w[sel], rcond=None)[0]
    Kz = lambda z: float(np.interp(-np.log1p(z), N, Kq))
    gdot = float(-fam.K1(q[-1]) * v[-1] / Know)
    # effective expansion law p = d ln a / d ln t at z = 100 and z = 10 (GR matter era: 2/3; omega_BD = 0 phi-MDE: 1/2)
    tt = sol.t
    lnt = np.log(np.maximum(tt - tt[0] + 1e-300, 1e-300) + 1.0 / (2 * Hi))   # t measured from a radiation-era origin
    p_of = lambda z: float(np.interp(-np.log1p(z), N, np.gradient(N, lnt)))
    cons = [abs(constraint_rel(fam, mu, sol.y[:, i])) for i in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 300))]
    return {'family': fam_name, 'm/H0': m, 'dmu (calibrated)': dmu, 'K(q_s) (start, normalisation)': Klab, 'a_i': a_i,
            'H(a=1) [units of the calibrated H_LCDM(1) = 1]': float(H[-1]),
            'Omega_m(a=1), K_now normalisation': float(rd[-1] / (3 * Know * H[-1]**2)),
            'K(z=1100)/K_now = G_now/G_rec': Kz(1100.0), 'K(z=1e3...BBN proxy a_i)/K_now': float(Kq[0] / Know),
            'K(z=1)/K_now': Kz(1.0),
            'Gdot/G today [H0]': gdot, 'Gdot/G today [1/yr] (H0 = 67.4)': gdot * H0_PER_YR,
            'w_eff(a=1)': float(w[-1]), 'CPL_fit_a>=0.4 (w0, wa)': [float(w0), float(wa)],
            'p = dlna/dlnt at z=100': p_of(100.0), 'p at z=10': p_of(10.0),
            'max_rel_constraint': float(max(cons))}


def rhs_fluid_factory(fam, mu, kappa):
    """q-theory FRW with one fluid of trace fraction kappa: w_f = (1 - kappa)/3 (T = -kappa rho_f)."""
    wf = (1.0 - kappa) / 3.0

    def rhs(t, y):
        q, v, H, rf = y
        K, k1, k2 = fam.K(q), fam.K1(q), fam.K2(q)
        rhoV = fam.eps(q) - mu * q
        Rq = 2.0 * (fam.eps1(q) - mu) / k1
        Kdd = (kappa * rf + 4.0 * rhoV - K * Rq) / 3.0 - 3.0 * H * k1 * v
        return [v, (Kdd - k2 * v * v) / k1, Rq / 6.0 - 2.0 * H * H, -3.0 * H * (1 + wf) * rf]
    return rhs


def run_relic_kappa(m, kappa, H_i_over_m=1e3, n_after=4.0):
    """Regime (iii): q starts at its vacuum value (mu = mu0, Minkowski point) in an early era dominated by a fluid with
    trace fraction kappa (radiation-like plasma with small trace); the trace drives K while H >> m (log growth), the field
    is left misaligned from the slow manifold at H ~ m and oscillates. Measures the oscillation (relic) energy
    E_osc = (3/4)(K'^2/K)[(v - v_slow)^2 + m^2 (q - q_slow)^2] relative to rho_f, n_after e-folds after H = m."""
    fam = Family('F1_quadratic_linearK', m)
    mu = fam.mu0
    rf_i = 3 * (H_i_over_m * m)**2
    Hi = H_from_constraint(fam, mu, 1.0, 0.0, rf_i, 0.0)
    rhs = rhs_fluid_factory(fam, mu, kappa)
    wf = (1 - kappa) / 3

    def ev(t, y):
        return y[2] - m * np.exp(-2.0 * n_after)       # H ~ a^-2: n_after e-folds after H = m
    ev.terminal = True
    per = 2 * np.pi / m
    sol = solve_ivp(rhs, (0, 1e9), [1.0, 0.0, Hi, rf_i], method='LSODA', rtol=1e-10, atol=1e-16,
                    max_step=per / 30, events=[ev], dense_output=True)
    q, v, H, rf = sol.y[:, -1]
    # slow manifold: quasi-static q for the current trace source kappa*rf (H << m here)
    f = lambda qq: kappa * rf + 4 * (fam.eps(qq) - mu * qq) - fam.K(qq) * 2 * (fam.eps1(qq) - mu) / fam.K1(qq)
    q_sl = bracket_root(f, 1.0, step=1e-4)
    k1 = fam.K1(q)
    E_osc = 0.75 * (k1**2 / fam.K(q)) * (v**2 + m**2 * (q - q_sl)**2)
    # log-growth estimate of the misalignment at onset: dK/K ~ (kappa/2) ln(t_osc/t_i), t ~ 1/(2H)
    est = 0.75 * (0.5 * kappa * np.log(H_i_over_m))**2
    return {'m': m, 'kappa': kappa, 'H_i/m': H_i_over_m, 'e-folds after H = m': n_after,
            'E_osc/rho_f at end': float(E_osc / rf),
            'E_osc/rho_f scaled back to H = m (x e^{-n(1-3wf)})': float(E_osc / rf * np.exp(-n_after * (1 - 3 * wf) - 0.0)),
            'relic fraction at onset / kappa^2': float(E_osc / rf / np.exp(n_after * (1 - 3 * wf) + 0.0) / kappa**2) if False else
            float(E_osc / rf / np.exp(n_after) / kappa**2),
            'log-growth estimate (3/4)[(kappa/2) ln(H_i/m)]^2 / kappa^2': est / kappa**2,
            'status': int(sol.status)}


def main():
    t0 = time.time()
    out = {'inputs': {'Omega_m': OMEGA_M, 'Omega_r': OMEGA_R, 'Omega_L': OMEGA_L,
                      'provenance': 'Omega_m: Planck 2018 (ledger row 9) [Q]; Omega_r ~ 9.1e-5 [R]; units H0 = 1, K(q0) = 1; '
                                    'families supplied [P], illustrative; F2 is the KV 0806.2805 eq (4.7)/(4.14) form'}}
    partI = []
    for fam_name in ['F1_quadratic_linearK', 'F2_KV0806_quartic_Kpropq', 'F3_cubic_quartic_expK']:
        for m in [10.0, 30.0, 100.0, 300.0]:
            r = run_adiabatic(fam_name, m)
            partI.append(r)
            print('I', fam_name, m, [(c['a'], round(c['H/m0 at a'], 4), round(c['ratio <Delta>/<rho_dyn>'], 5))
                                     for c in r['comparison']], flush=True)
    # error coefficient: (1 - ratio)/(H/m)^2 at a = 1, per family (should be ~constant in m if the error is O(H^2/m^2))
    coef = {}
    for r in partI:
        c1 = [c for c in r['comparison'] if abs(c['a'] - 1.0) < 1e-9][0]
        coef.setdefault(r['family'], []).append([r['m_target'], (1 - c1['ratio <Delta>/<rho_dyn>']) / c1['H/m0 at a']**2])
    out['partI_adiabatic'] = partI
    out['partI_error_coefficient (1-ratio)/(H/m)^2 at a=1'] = coef
    partII = []
    for m in [0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1000.0]:
        r = run_light_scan(m)
        partII.append(r)
        print('II', m, {k: r[k] for k in ['H(a=1) [units of the calibrated H_LCDM(1) = 1]', 'K(z=1100)/K_now = G_now/G_rec',
                                          'K(z=1e3...BBN proxy a_i)/K_now', 'p = dlna/dlnt at z=100']}, flush=True)
    out['partII_light_scan'] = partII
    partIII = []
    for kap in [1e-3, 3e-3, 1e-2]:
        for hi in [1e2, 1e3]:
            r = run_relic_kappa(10.0, kap, H_i_over_m=hi)
            partIII.append(r)
            print('III', kap, hi, r['relic fraction at onset / kappa^2'], flush=True)
    out['partIII_relic_kappa'] = partIII
    out['partIII_fit'] = {'relic fraction at onset / [kappa^2 ln^2(H_i/m)]':
                          [r['relic fraction at onset / kappa^2'] / np.log(r['H_i/m'])**2 for r in partIII]}
    out['runtime_s'] = time.time() - t0
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'frw_regimes.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
