"""Track `bound`, part B: disappearing-Lambda / chameleon f(R) models as q-theories. Writes results/hu_sawicki.json.

B2  Vacuum structure of Hu-Sawicki n=1 (arXiv:0705.1158), Starobinsky n=2 (0706.2041) and Appleby-Battye (0705.3199):
    f(0), f'(0), f''(0) (P_q's Minkowski equilibrium: healthy?), de Sitter root R_1 with its stability, and (HS) the
    point R_* where f' = 0 between Minkowski and de Sitter.
B3  Hu-Sawicki n=1 background (HS eqs. (reducedfriedmann a,b), Radau), calibrated c1/c2 = 6 Om_L/Om_m to Omega_m = 0.315
    (the single calibrated number: the effective Lambda), w_eff(z), CPL projection, future de Sitter attractor, the
    first-order analytic formula of adiabatic_relation.py S1 as a cross-check, and quasi-static linear growth.
B4  Threshold re-selection under P_q for a q-realization (linear K vs affine gauge), numerically (Appleby-Battye form,
    compressed units).
Conventions: HS action (R + f_HS)/(2 kappa^2), units m_HS^2 = kappa^2 rho_m0/3 = 1 (so kappa^2 rho_m = 3 a^-3).
Run: python3 hu_sawicki.py  (~1-2 min)
"""
import json
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OM_T = 0.315        # tilde Omega_m used in the calibration c1/c2 = 6 (1 - Om)/Om (Planck 2018, ledger row 9) [Q]


class HS:
    """f_HS(R) = -c1 x^n/(c2 x^n + 1), x = R/m^2, m^2 = 1."""

    def __init__(self, fR0, n=1, Om=OM_T):
        self.n, self.Om = n, Om
        self.c1_over_c2 = 6.0 * (1 - Om) / Om                       # calibrated: effective Lambda
        x0 = 12.0 / Om - 9.0
        self.c1_over_c22 = abs(fR0) * x0**(n + 1) / n               # HS approximation f_R0 = -n c1/c2^2 x0^{-n-1}
        self.c2 = self.c1_over_c2 / self.c1_over_c22
        self.c1 = self.c1_over_c2 * self.c2
        self.fR0_target = fR0

    def f(self, R):
        x = R ** self.n
        return -self.c1 * x / (self.c2 * x + 1)

    def fR(self, R):
        n = self.n
        return -self.c1 * n * R**(n - 1) / (self.c2 * R**n + 1)**2

    def fRR(self, R):
        n, c1, c2 = self.n, self.c1, self.c2
        x = R**n
        return -c1 * n * ((n - 1) * R**(n - 2) * (c2 * x + 1) - 2 * c2 * n * R**(2 * n - 2)) / (c2 * x + 1)**3


def hs_rhs(model):
    def rhs(N, y):
        yH, yR = y
        a3 = np.exp(-3 * N)
        R = yR + 3 * a3
        fR, fRR, f = model.fR(R), model.fRR(R), model.f(R)
        dyH = yR / 3 - 4 * yH
        dyR = 9 * a3 - (1 / (yH + a3)) * (1 / fRR) * (yH - fR * (yR / 6 - yH - a3 / 2) + f / 6)
        return [dyH, dyR]
    return rhs


def hs_jac(model):
    def jac(N, y):
        yH, yR = y
        a3 = np.exp(-3 * N)
        R = yR + 3 * a3
        fR, fRR, f = model.fR(R), model.fRR(R), model.f(R)
        h = 1e-6 * R
        fRRR = (model.fRR(R + h) - model.fRR(R - h)) / (2 * h)
        Q = 1 / ((yH + a3) * fRR)
        S = yH - fR * (yR / 6 - yH - a3 / 2) + f / 6
        dQ_dyH = -Q / (yH + a3)
        dS_dyH = 1 + fR
        dQ_dyR = -Q * fRRR / fRR
        dS_dyR = -fRR * (yR / 6 - yH - a3 / 2)
        return np.array([[-4.0, 1.0 / 3.0], [-(dQ_dyH * S + Q * dS_dyH), -(dQ_dyR * S + Q * dS_dyR)]])
    return jac


def hs_initial(model, N_i):
    """quasi-static (minimum of V_eff) initial data: R (1 - f_R) + 2 f = 3 a^-3 (HS units), y_H = -f/6 + O(f_R)."""
    a3 = np.exp(-3 * N_i)
    g = lambda R: R - model.fR(R) * R + 2 * model.f(R) - 3 * a3
    R = brentq(g, 3 * a3 * 0.5, 3 * a3 * 2 + 10 * model.c1_over_c2)
    yR = R - 3 * a3
    yH = -model.f(R) / 6 - model.fR(R) * a3 / 2          # Friedmann at first order in f_R with H^2 ~ a^-3
    return [yH, yR]


def solve_hs(model, N_i=np.log(1 / 11.), N_f=0.0, rtol=1e-7):
    sol = solve_ivp(hs_rhs(model), (N_i, N_f), hs_initial(model, N_i), method='Radau', rtol=rtol, atol=1e-13,
                    jac=hs_jac(model), dense_output=True)
    return sol


def analytic_first_order(model, N):
    """adiabatic_relation.py S1 with F = f_HS + 2 Lam (F -> 0 at high R), GR+Lam background; returns y_H = rho_DE/3."""
    a3 = np.exp(-3 * N)
    Lam = model.c1_over_c2 / 2
    rho = 3 * a3
    R = rho + 4 * Lam
    H2 = (rho + Lam) / 3
    F = model.f(R) + 2 * Lam
    Fp, Fpp = model.fR(R), model.fRR(R)
    HRdot = -3 * H2 * rho                                   # H * dR/dt with dR/dt = -3 H rho
    rhoDE = Lam + (R * Fp - F) / 2 - 3 * Fp * H2 - 3 * Fpp * HRdot
    return rhoDE / 3


def w_and_cpl(Ngrid, yH):
    w = -1 - np.gradient(np.log(yH), Ngrid) / 3
    a = np.exp(Ngrid)
    sel = a >= 1 / 3.0
    A = np.vstack([np.ones(sel.sum()), 1 - a[sel]]).T
    w0, wa = np.linalg.lstsq(A, w[sel], rcond=None)[0]
    # local Taylor projection at a = 1
    wa_loc = -np.gradient(w, a)[-1]
    return w, float(w0), float(wa), float(w[-1]), float(wa_loc)


def growth(model, sol, k_over_m, N_i, N_f=0.0, lcdm=False):
    """quasi-static sub-horizon growth: delta'' + (2 + dlnH/dN) delta' - 1.5 (a^-3/(H^2)) (G_eff/G) delta = 0.
    Background: GR + Lam_eff (HS units), R = R_GR = 3 a^-3 + 4 Lam_eff; the O(f_R0) background shift is neglected
    (it is ~1e-4 relative at |f_R0| = 1e-4, B3), so only G_eff(k, a) differs between the runs. `sol` is unused."""
    Lam = model.c1_over_c2 / 2.0

    def H2_of(N):
        return np.exp(-3 * N) + Lam / 3.0

    def rhs(N, y):
        d, dp = y
        h = 1e-4
        dlnH = 0.5 * (np.log(H2_of(N + h)) - np.log(H2_of(N - h))) / (2 * h)
        a = np.exp(N)
        if lcdm:
            geff = 1.0
        else:
            R = 3 * np.exp(-3 * N) + 4 * Lam
            fR, fRR = model.fR(R), model.fRR(R)
            u = (k_over_m / a)**2 * fRR / (1 + fR)
            geff = (1 + 4 * u) / (1 + 3 * u) / (1 + fR)
        return [dp, -(2 + dlnH) * dp + 1.5 * np.exp(-3 * N) / H2_of(N) * geff * d]
    a_i = np.exp(N_i)
    s = solve_ivp(rhs, (N_i, N_f), [a_i, a_i], method='DOP853', rtol=1e-9, atol=1e-14, dense_output=True)
    return s


def vacua_HS(model):
    """total f = R + f_HS (HS units). Returns Minkowski data, the dS root, the f' = 0 point (n = 1)."""
    ft = lambda R: R + model.f(R)
    fp = lambda R: 1 + model.fR(R)
    fpp = model.fRR
    out = {'f(0)': float(ft(0.0)), "f'(0) = 1 + f_R(0)": float(fp(1e-300) if model.n == 1 else 1.0),
           "f''(0)": float(fpp(1e-12)),
           'Minkowski status': ('ghost: f\'(0) = 1 - c1 < 0' if (model.n == 1 and model.c1 > 1) else
                                ('tachyonic: f\'\'(0) < 0' if fpp(1e-12) < 0 else 'healthy'))}
    # dS roots of R f' - 2 f = 0 for R > 0
    g = lambda R: R * fp(R) - 2 * ft(R)
    grid = np.logspace(-6, 3, 20000)
    v = g(grid)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    roots = [brentq(g, grid[i], grid[i + 1], xtol=1e-14) for i in idx]
    out['dS roots R/m^2'] = roots
    if roots:
        R1 = max(roots)
        out['R_1 (largest dS root)'] = R1
        out['R_1 / (4 Lam_eff), Lam_eff = c1/(2 c2)'] = R1 / (2 * model.c1_over_c2)
        out["f'(R_1)"] = float(fp(R1))
        out["stability f'(R_1)/(R_1 f''(R_1)) > 1"] = float(fp(R1) / (R1 * fpp(R1)))
        out['scalaron m^2(R_1)/H_1^2 = (f\'/f\'\' - R)/3 / (R_1/12)'] = float(((fp(R1) / fpp(R1)) - R1) / 3 / (R1 / 12))
    if model.n == 1 and model.c1 > 1:
        Rstar = (np.sqrt(model.c1) - 1) / model.c2
        out["R_* where f' = 0 (n=1, c1 > 1): R_*/m^2"] = float(Rstar)
        out['R_*/R_1'] = float(Rstar / max(roots)) if roots else None
    out['c1'], out['c2'] = model.c1, model.c2
    return out


def starobinsky_vacua(lam, n=2):
    """f = R + lam R0 [(1 + R^2/R0^2)^-n - 1], units R0 = 1 (0706.2041 eq (f))."""
    ft = lambda R: R + lam * ((1 + R * R)**(-n) - 1)
    fp = lambda R: 1 - 2 * n * lam * R * (1 + R * R)**(-n - 1)
    fpp = lambda R: -2 * n * lam * ((1 + R * R)**(-n - 1) - 2 * (n + 1) * R * R * (1 + R * R)**(-n - 2))
    g = lambda R: R * fp(R) - 2 * ft(R)
    grid = np.linspace(1e-4, 4 * lam + 5, 200001)
    v = g(grid)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    roots = [brentq(g, grid[i], grid[i + 1], xtol=1e-14) for i in idx]
    R1 = max(roots)
    return {'lambda': lam, 'n': n, 'f(0)': ft(0.0), "f'(0)": fp(0.0), "f''(0)": fpp(0.0),
            "Minkowski scalaron m^2 = f'/(3 f'') at 0 [units R0]": fp(0.0) / (3 * fpp(0.0)),
            'dS roots R/R0': roots, 'x_1 = R_1/R0': R1, "f'(R_1)": fp(R1),
            "stability f'(R_1) > R_1 f''(R_1)": bool(fp(R1) > R1 * fpp(R1)),
            "f''(R_1) > 0": bool(fpp(R1) > 0), 'R_2 = R0/sqrt(2n+1) (f\'\' = 0)': 1 / np.sqrt(2 * n + 1)}


def appleby_battye(a_, b_):
    """F = R/2 + (1/2a) log[cosh(aR) - tanh(b) sinh(aR)]  (0705.3199 eq (ab3)), units of their M^2 = 1."""
    F = lambda R: R / 2 + np.log(np.cosh(a_ * R) - np.tanh(b_) * np.sinh(a_ * R)) / (2 * a_)
    Fp = lambda R: 0.5 * (1 + np.tanh(a_ * R - b_))
    Fpp = lambda R: 0.5 * a_ / np.cosh(a_ * R - b_)**2
    g = lambda R: R * Fp(R) - 2 * F(R)
    grid = np.linspace(1e-6, 4 * b_ / a_ + 5 / a_, 400001)
    v = g(grid)
    idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    roots = [brentq(g, grid[i], grid[i + 1], xtol=1e-14) for i in idx]
    R1 = max(roots)
    return {'a': a_, 'b': b_, 'F(0)': float(F(0.0)), "F'(0) = 1/(1+e^{2b})": float(Fp(0.0)), "F''(0)": float(Fpp(0.0)),
            "G_Minkowski/G_dS = F'(R_1)/F'(0)": float(Fp(R1) / Fp(0.0)),
            "Minkowski scalaron m^2 = F'/(3F'') at 0": float(Fp(0.0) / (3 * Fpp(0.0))),
            'dS/vacuum roots R (excl. 0)': roots, 'R_1': R1, "F'(R_1)": float(Fp(R1)),
            "stable dS: F'(R_1) > R_1 F''(R_1)": bool(Fp(R1) > R1 * Fpp(R1))}


def threshold_check():
    """B4: q-realizations of an Appleby-Battye f (compressed units a = 2, b = 3) under eps -> eps + c with P_q re-selection.
    Realization (mu = 0):  eps(q) = [K(q) R(q) - F(R(q))]/2, so f = K R - 2 eps = F and the q-equation eps' = K' R/2 holds
    identically.  (A) linear K(q) = F'(0) + K1 (q - q0), R(q) = (F')^{-1}(K(q));  (B) affine gauge R(q) = lam (q - q0),
    K(q) = F'(R(q)).  After eps -> eps + c, P_q selects (q0', mu'):  eps'(q0') = mu', eps(q0') + c - mu' q0' = 0, and
        f_new(R) = K(q_R) R - 2 (eps(q_R) + c) + 2 mu' q_R,   eps'(q_R) - K'(q_R) R/2 = mu'.
    Compared with the first-order predictions  F(R + dR) - dR F'(0), dR = 2c/(q0 K1)  (A)  and  F(R) + 2cR/(lam q0)  (B)."""
    a_, b_ = 2.0, 3.0
    F = lambda R: R / 2 + np.log(np.cosh(a_ * R) - np.tanh(b_) * np.sinh(a_ * R)) / (2 * a_)
    Fp = lambda R: 0.5 * (1 + np.tanh(a_ * R - b_))
    Rinv = lambda K: (np.arctanh(2 * K - 1) + b_) / a_            # (F')^{-1}
    q0 = 1.0
    out = {}
    for gauge in ['linearK', 'affine']:
        if gauge == 'linearK':
            K1, lam = 1.0, None
            Kq = lambda q: Fp(0.0) + K1 * (q - q0)
            dKq = lambda q: K1
            Rq = lambda q: Rinv(Kq(q))
            qlo, qhi = q0 - Fp(0.0) + 1e-9, q0 + (1 - Fp(0.0)) - 1e-9
        else:
            lam = 1.0
            Rq = lambda q: lam * (q - q0)
            Kq = lambda q: Fp(Rq(q))
            dKq = lambda q: lam * 0.5 * a_ / np.cosh(a_ * Rq(q) - b_)**2
            qlo, qhi = q0 - 3.0, q0 + 6.0
        eps = lambda q: (Kq(q) * Rq(q) - F(Rq(q))) / 2
        deps = lambda q: dKq(q) * Rq(q) / 2                        # exact for this realization (mu = 0)
        h = 1e-6
        qq = np.linspace(qlo + 1e-3, qhi - 1e-3, 7)
        res = [(eps(x + h) - eps(x - h)) / (2 * h) - deps(x) for x in qq]
        out[gauge + ': q-equation residual (finite-difference eps\' - K\'R/2, max)'] = float(np.max(np.abs(res)))
        for c in ([1e-3, 1e-2] if gauge == 'linearK' else [1e-7, 1e-6]):
            gsel = lambda q: eps(q) + c - deps(q) * q
            # nearest root to q0 on a fine grid, then brentq
            grid = np.linspace(max(qlo, q0 - 0.2), min(qhi, q0 + 0.2), 40001)
            v = np.array([gsel(x) for x in grid])
            idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
            j = idx[np.argmin(np.abs(grid[idx] - q0))]
            q0n = brentq(gsel, grid[j], grid[j + 1], xtol=1e-15)
            mun = deps(q0n)
            Rs = np.linspace(0.0, 4.0 if gauge == 'linearK' else 3.5, 81)
            fnew, fpred = [], []
            for R in Rs:
                gq = lambda q: deps(q) - dKq(q) * R / 2 - mun
                if gauge == 'linearK':
                    qR = brentq(gq, qlo + 1e-7, qhi - 1e-7, xtol=1e-15)
                    dR = 2 * c / (q0 * K1)
                    fpred.append(F(R + dR) - dR * Fp(0.0))
                else:
                    gg = np.linspace(q0 + R / lam - 2.0, q0 + R / lam + 2.0, 4001)
                    vv = gq(gg)
                    ii = np.where(np.sign(vv[:-1]) != np.sign(vv[1:]))[0]
                    jj = ii[np.argmin(np.abs(gg[ii] - (q0 + R / lam)))]
                    qR = brentq(gq, gg[jj], gg[jj + 1], xtol=1e-15)
                    fpred.append(F(R) + 2 * c * R / (lam * q0))
                fnew.append(Kq(qR) * R - 2 * (eps(qR) + c) + 2 * mun * qR)
            fnew, fpred = np.array(fnew), np.array(fpred)
            dev = np.max(np.abs(fnew - F(Rs)))
            out['%s c=%g: (q0\' - q0), mu\' vs first order c/(q0 eps\'\'), c/q0' % (gauge, c)] = [float(q0n - q0), float(mun), float(c / q0)]
            out['%s c=%g: max|f_new - F|' % (gauge, c)] = float(dev)
            out['%s c=%g: max|f_new - prediction| / max|f_new - F|' % (gauge, c)] = float(np.max(np.abs(fnew - fpred)) / dev)
            # where is the de Sitter vacuum of f_new? (R f' = 2 f)
            fn = lambda R: np.interp(R, Rs, fnew)
            Rg = Rs[1:-1]
            fprime = np.gradient(fnew, Rs)
            gvac = Rs * fprime - 2 * fnew
            k = np.where(np.sign(gvac[1:-1][:-1]) != np.sign(gvac[1:-1][1:]))[0]
            out['%s c=%g: dS roots of f_new (grid)' % (gauge, c)] = [float(Rg[i]) for i in k]
    gv = lambda R: R * Fp(R) - 2 * F(R)
    out['unperturbed AB dS roots'] = [float(brentq(gv, x0, x1)) for x0, x1 in [(0.5, 2.0), (2.0, 3.5)] if np.sign(gv(x0)) != np.sign(gv(x1))]
    return out


def main():
    t0 = time.time()
    out = {'conventions': 'HS units m^2 = kappa^2 rho_m0/3 = 1; f_HS is the modification (action (R + f_HS)/2kappa^2)',
           'calibration': 'c1/c2 = 6 (1 - Om)/Om with Om = %.3f (Planck 2018) [calibrated: the effective Lambda]' % OM_T}
    # ---------------- B2 vacua
    vac = {}
    for fR0 in [1e-2, 1e-4, 1e-6]:
        vac['HS n=1 |f_R0|=%g' % fR0] = vacua_HS(HS(-fR0, 1))
    vac['HS n=2 |f_R0|=1e-4'] = vacua_HS(HS(-1e-4, 2))
    vac['Starobinsky n=2 lambda=1'] = starobinsky_vacua(1.0, 2)
    vac['Starobinsky n=2 lambda=2'] = starobinsky_vacua(2.0, 2)
    vac['Appleby-Battye a=2,b=3 (compressed)'] = appleby_battye(2.0, 3.0)
    vac['Appleby-Battye a=10,b=10 (compressed)'] = appleby_battye(10.0, 10.0)
    out['B2_vacua'] = vac
    for k, v in vac.items():
        print(k, {kk: vv for kk, vv in v.items() if not isinstance(vv, list)})
    # ---------------- B3 background
    Ngrid = np.linspace(np.log(1 / 3.0), 0.0, 400)
    bg = {}
    for fR0 in [1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        model = HS(-fR0, 1)
        row = {}
        for N_i in ([np.log(1 / 11.), np.log(1 / 21.)] if fR0 >= 1e-3 else ([np.log(1 / 11.)] if fR0 >= 1e-4 else [])):
            try:
                sol = solve_hs(model, N_i=N_i)
                yH = sol.sol(Ngrid)[0]
                w, w0, wa, wloc0, waloc = w_and_cpl(Ngrid, yH)
                R_today = sol.sol(0.0)[1] + 3
                row['N_i=%.2f' % N_i] = {'status': int(sol.status), 'w0_CPL': w0, 'wa_CPL': wa, 'w(z=0)': wloc0,
                                         'wa_local': waloc, 'w(z=1)': float(np.interp(np.log(0.5), Ngrid, w)),
                                         'f_R today (actual)': float(model.fR(R_today)),
                                         'Omega_m true = 1/(1 + y_H(0))': float(1 / (1 + sol.sol(0.0)[0])),
                                         'nfev': int(sol.nfev)}
            except Exception as e:
                row['N_i=%.2f' % N_i] = {'error': repr(e)}
        yHa = analytic_first_order(model, Ngrid)
        w, w0, wa, wloc0, waloc = w_and_cpl(Ngrid, yHa)
        row['first_order_analytic'] = {'w0_CPL': w0, 'wa_CPL': wa, 'w(z=0)': wloc0, 'w(z=1)': float(np.interp(np.log(0.5), Ngrid, w)),
                                       'Omega_m': float(1 / (1 + yHa[-1]))}
        bg['|f_R0|=%g' % fR0] = row
        print('B3', fR0, {k: ({kk: (round(vv, 7) if isinstance(vv, float) else vv) for kk, vv in v.items()}) for k, v in row.items()})
    out['B3_background'] = bg
    # w_eff(z) table (phantom crossing; HS 2007 state w_eff < -1 at high z for these models)
    wz = {}
    for fR0 in [1e-2, 1e-4]:
        model = HS(-fR0, 1)
        Nw = np.linspace(np.log(1 / 11.), 0.0, 800)
        yH = solve_hs(model).sol(Nw)[0]
        ww = -1 - np.gradient(np.log(yH), Nw) / 3
        wz['|f_R0|=%g (full)' % fR0] = {('z=%g' % z): float(np.interp(-np.log1p(z), Nw, ww)) for z in [0, 0.5, 1, 1.5, 2, 3, 5]}
        zz = np.exp(-Nw) - 1
        sel = (zz < 4)          # start-up transients near z_i = 10 excluded
        s = np.where(np.sign(ww[sel][:-1] + 1) != np.sign(ww[sel][1:] + 1))[0]
        wz['|f_R0|=%g (full)' % fR0]['phantom-divide crossing z'] = [float(zz[sel][i]) for i in s]
    out['B3_w_of_z'] = wz
    print('wz', wz)
    # scaling of (1 + w0, wa) with f_R0 from the analytic first-order formula, and values at the published bounds
    bounds = {'Schmidt-Vikhlinin-Hu 2009 clusters (95%)': 1.3e-4, 'Cataneo+ 2015 clusters (95.4%)': 10**-4.79,
              'DESI DR1 full shape, Gonzalez+ 2026 (95%)': 10**-4.59, 'BOSS DR12 full shape, Aviles 2024 (95%)': 1.53e-5,
              'Hu-Sawicki 2007 Milky Way static thin shell (v=300 km/s)': 2e-6,
              'Jain-Vikram-Sakstein 2013 distance indicators (95%)': 5e-7,
              'Desmond-Ferreira 2020 galaxy morphology (1 sigma)': 1.4e-8, 'Desmond-Ferreira 2020 (2 sigma)': 2.8e-8}
    at_b = {}
    for name, fb in bounds.items():
        model = HS(-fb, 1)
        yHa = analytic_first_order(model, Ngrid)
        w, w0, wa, wloc0, waloc = w_and_cpl(Ngrid, yHa)
        at_b[name] = {'|f_R0|': fb, '1 + w0 (CPL)': 1 + w0, 'wa (CPL)': wa, '1 + w(z=0)': 1 + wloc0,
                      'Compton wavelength today [Mpc/h] ~ 1/m_s, m_s^2 = (1+f_R)/(3 f_RR) at R0':
                          float(2997.92458 / np.sqrt(OM_T) / np.sqrt((1 + model.fR(3 + 2 * model.c1_over_c2)) / (3 * model.fRR(3 + 2 * model.c1_over_c2))))}
    out['B3_at_published_fR0_bounds (first-order formula; validated against full integration above)'] = at_b
    for k, v in at_b.items():
        print('bound', k, v)
    # DESI DR2 comparison (ledger row 8 numbers): distance of the model point from each DESI posterior, diagonal Gaussian
    desi = {'DESI+CMB': (-0.42, 0.21, -1.75, 0.58), 'DESI+CMB+PantheonPlus': (-0.838, 0.055, -0.62, 0.205),
            'DESI+CMB+Union3': (-0.667, 0.088, -1.09, 0.29), 'DESI+CMB+DESY5': (-0.752, 0.057, -0.86, 0.215)}
    comp = {}
    for fname, (fb) in [('|f_R0| = 1e-6 (first-order formula)', 1e-6), ('|f_R0| = 1e-4 (full)', 1e-4),
                        ('|f_R0| = 1e-2 (full)', 1e-2), ('|f_R0| = 1e-1 (full; excluded by every bound)', 1e-1)]:
        model = HS(-fb, 1)
        yH = analytic_first_order(model, Ngrid) if fb < 1e-5 else solve_hs(model).sol(Ngrid)[0]
        w, w0, wa, _, _ = w_and_cpl(Ngrid, yH)
        comp[fname] = {'w0': w0, 'wa': wa}
        for dn, (mw0, sw0, mwa, swa) in desi.items():
            comp[fname]['chi (diag) vs ' + dn] = float(np.hypot((w0 - mw0) / sw0, (wa - mwa) / swa))
    for dn, (mw0, sw0, mwa, swa) in desi.items():
        comp.setdefault('LCDM point', {})['chi (diag) vs ' + dn] = float(np.hypot((-1 - mw0) / sw0, (0 - mwa) / swa))
    comp['note'] = ('diagonal-Gaussian distance in (w0, wa) ignoring the DESI w0-wa correlation; indicative only. DESI DR2 quotes '
                    'the LCDM preference as 2.8-4.2 sigma from Delta chi^2_MAP (ledger row 8) [Q].')
    out['B3_vs_DESI_DR2'] = comp
    print('DESI', comp)
    # future attractor
    fut = {}
    for fR0 in [1e-2, 1e-4]:
        model = HS(-fR0, 1)
        sol = solve_ivp(hs_rhs(model), (np.log(1 / 11.), np.log(200.0)), hs_initial(model, np.log(1 / 11.)), method='Radau',
                        rtol=1e-7, atol=1e-12, jac=hs_jac(model), dense_output=True)
        Rf = sol.sol(np.log(200.0))[1] + 3 * 200.0**-3
        vac1 = vacua_HS(model)
        fut['|f_R0|=%g' % fR0] = {'R(a=200)/m^2': float(Rf), 'R_1 dS root': vac1.get('R_1 (largest dS root)'),
                                  'relative difference': float(Rf / vac1['R_1 (largest dS root)'] - 1),
                                  "1 + f_R(a=200)": float(1 + model.fR(Rf)), 'status': int(sol.status)}
    out['B3_future_attractor'] = fut
    print('future', fut)
    # growth
    gr = {}
    model_l = HS(-1e-6, 1)
    N_i = np.log(1 / 31.)
    sl = growth(model_l, None, 1.0, N_i, lcdm=True)
    for fR0 in [1e-4, 1e-5, 1e-6, 1e-7, 1.4e-8]:
        model = HS(-fR0, 1)
        sol = None
        row = {}
        for k_hMpc in [0.02, 0.1, 0.2]:
            k_over_m = k_hMpc * 2997.92458 / np.sqrt(OM_T)
            sg = growth(model, sol, k_over_m, N_i)
            Ne = np.log(1 / 1.5)
            d, dp = sg.sol(Ne)
            dl, dpl = sl.sol(Ne)
            row['k=%g h/Mpc' % k_hMpc] = {'D/D_LCDM at z=0.5': float(d / dl), 'f/f_LCDM at z=0.5': float((dp / d) / (dpl / dl)),
                                          'f D/(f D)_LCDM at z=0.5 (f sigma8 proxy)': float(dp / dpl),
                                          'D/D_LCDM at z=0': float(sg.sol(0.0)[0] / sl.sol(0.0)[0])}
        gr['|f_R0|=%g' % fR0] = row
        print('growth', fR0, row)
    out['B3_linear_growth_quasistatic'] = gr
    # ---------------- B4
    out['B4_threshold_reselection'] = threshold_check()
    print('B4', out['B4_threshold_reselection'])
    out['runtime_s'] = time.time() - t0
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'hu_sawicki.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
