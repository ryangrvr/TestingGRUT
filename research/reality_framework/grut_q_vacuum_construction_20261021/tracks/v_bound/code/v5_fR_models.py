"""v_bound: independent checks of bound part B (published disappearing-Lambda f(R) models as q-theories).

Written from scratch (nothing imported from ../../bound/code).  Formulation differs from bound's (HS reduced y_H, y_R):
the Hu-Sawicki background is integrated in scalar-tensor form, phi = f_R, with
   H^2 = [rho_m + (phi R - f)/2] / (3 (1 + phi + psi)),   psi = dphi/dN,        R = R(phi) (explicit inverse, n = 1)
   psi' = [rho_m + phi R - R - 2 f]/(3 H^2) - psi (R/(6 H^2) + 1)
(the f(R) trace equation and 00 equation, V1a), units m^2 = kappa^2 rho_m0/3 = 1, kappa^2 = 1, rho_m = 3 a^-3.
Checks:
 B3  vacua: HS n=1 (f'(0) = 1 - c1, R_* where f' = 0, dS roots and stability), HS n=2, Starobinsky n=2 (lambda = 1, 2),
     Appleby-Battye F'(0).
 B5  w_eff(z), CPL (w0, wa) with three fit weightings, Omega_m(true), phantom crossing, future attractor, small-field
     scaling (first-order formula, own implementation), Compton wavelengths at the published f_R0 bounds, DESI diagonal
     distances (both symmetric-average and correct-side asymmetric errors), linear growth incl. a run on the true f(R)
     background.
 B4  exact threshold re-selection in the linear-K realization: f_new(R) = f(R + delta) - delta f'(0) + 2 (mu' q0 - c),
     delta = 2 mu'/K1 (derived here; verified numerically on an HS-type f with c1 < 1).
Run: python3 v5_fR_models.py   (writes ../results/v5_fR_models.json)
"""
import json
import os
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OMT = 0.315
OUT = {}


class HS1:
    """total F(R) = R + f, f = -c1 R/(c2 R + 1) (n = 1, m^2 = 1)."""

    def __init__(self, fR0abs, Om=OMT):
        self.cc = 6 * (1 - Om) / Om                 # c1/c2 (calibrated, effective Lambda)
        x0 = 12 / Om - 9
        self.c1c22 = fR0abs * x0**2                 # c1/c2^2 from f_R0 = -(c1/c2^2) x0^-2 (HS approximation)
        self.c2 = self.cc / self.c1c22
        self.c1 = self.cc * self.c2

    def f(self, R):
        return -self.c1 * R / (self.c2 * R + 1)

    def fR(self, R):
        return -self.c1 / (self.c2 * R + 1)**2

    def fRR(self, R):
        return 2 * self.c1 * self.c2 / (self.c2 * R + 1)**3

    def R_of_phi(self, phi):
        return (np.sqrt(self.c1 / (-phi)) - 1) / self.c2


def hs_rhs(mod):
    def rhs(N, y):
        phi, psi = y
        R = mod.R_of_phi(phi)
        f = mod.f(R)
        rm = 3 * np.exp(-3 * N)
        H2 = (rm + (phi * R - f) / 2) / (3 * (1 + phi + psi))
        dpsi = (rm + phi * R - R - 2 * f) / (3 * H2) - psi * (R / (6 * H2) + 1)
        return [psi, dpsi]
    return rhs


def hs_H2(mod, N, phi, psi):
    R = mod.R_of_phi(phi)
    return (3 * np.exp(-3 * N) + (phi * R - mod.f(R)) / 2) / (3 * (1 + phi + psi))


def hs_init(mod, Ni):
    """quasi-static: trace without derivatives, rho_m + phi R - R - 2 f = 0 with phi = f_R(R); psi by differencing."""
    def Rqs(N):
        rm = 3 * np.exp(-3 * N)
        g = lambda R: rm + mod.fR(R) * R - R - 2 * mod.f(R)
        return brentq(g, rm * 0.5, rm * 2 + 10 * mod.cc, xtol=1e-14, rtol=1e-15)
    h = 1e-4
    phi = mod.fR(Rqs(Ni))
    psi = (mod.fR(Rqs(Ni + h)) - mod.fR(Rqs(Ni - h))) / (2 * h)
    return [phi, psi]


def solve_bg(mod, zi=10.0, Nf=0.0, rtol=1e-9):
    Ni = -np.log1p(zi)
    return solve_ivp(hs_rhs(mod), (Ni, Nf), hs_init(mod, Ni), method='Radau', rtol=rtol, atol=1e-14, dense_output=True)


def w_of(Ng, yH):
    return -1 - np.gradient(np.log(yH), Ng) / 3


def cpl(Ng, w, mode):
    a = np.exp(Ng)
    if mode == 'uniform in N (bound)':
        sel = a >= 1 / 3
        aa, ww = a[sel], w[sel]
    else:
        if mode == 'uniform in a':
            aa = np.linspace(1 / 3, 1, 400)
        else:                                     # uniform in z
            aa = 1 / (1 + np.linspace(0, 2, 400))
        ww = np.interp(np.log(aa), Ng, w)
    A = np.vstack([np.ones_like(aa), 1 - aa]).T
    w0, wa = np.linalg.lstsq(A, ww, rcond=None)[0]
    return float(w0), float(wa)


def first_order_yH(mod, Ng):
    """own implementation of rho_DE,eff = Lam + (R F' - F)/2 - 3 F' H^2 - 3 H F'' Rdot on GR + Lam (F = f + 2 Lam)."""
    Lam = mod.cc / 2
    rm = 3 * np.exp(-3 * Ng)
    R = rm + 4 * Lam
    H2 = (rm + Lam) / 3
    F = mod.f(R) + 2 * Lam
    Fp, Fpp = mod.fR(R), mod.fRR(R)
    HRdot = H2 * (-3 * rm)                       # H dR/dt = H^2 dR/dN, dR/dN = -3 rho_m
    rho = Lam + (R * Fp - F) / 2 - 3 * Fp * H2 - 3 * Fpp * HRdot
    return rho / 3


# ------------------------------------------------------------------------------------------------ B3 vacua
vac = {}
for fr in [1e-2, 1e-4, 1e-6]:
    mod = HS1(fr)
    Ft = lambda R: R + mod.f(R)
    Fp = lambda R: 1 + mod.fR(R)
    g = lambda R: R * Fp(R) - 2 * Ft(R)
    Rg = np.logspace(-4, 3, 50001)
    v = g(Rg)
    roots = [brentq(g, Rg[i], Rg[i + 1], xtol=1e-14) for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]]
    R1 = max(roots)
    Rstar = brentq(Fp, 1e-12, R1)
    vac['HS n=1 |f_R0|=%g' % fr] = {'c1': mod.c1, "f'(0) = 1 - c1": 1 - mod.c1, 'dS roots': roots, 'R_1': R1,
                                     'R_1/(2 c1/c2)': R1 / (2 * mod.cc), "R_* (f' = 0)": Rstar, 'R_*/R_1': Rstar / R1,
                                     "f'(R_1)/(R_1 f''(R_1))": Fp(R1) / (R1 * mod.fRR(R1))}
# HS n = 2: f = -c1 R^2/(c2 R^2 + 1)
fr = 1e-4
x0 = 12 / OMT - 9
cc = 6 * (1 - OMT) / OMT
c1c22 = fr * x0**3 / 2
c2 = cc / c1c22
c1 = cc * c2
f2 = lambda R: -c1 * R**2 / (c2 * R**2 + 1)
f2pp0 = -2 * c1
vac['HS n=2 |f_R0|=1e-4'] = {"f''(0) = -2 c1": f2pp0, 'c1': c1}
# Starobinsky n = 2: f = R + lam [(1 + R^2)^-2 - 1] (R0 = 1)
for lam in [1.0, 2.0]:
    Ft = lambda R: R + lam * ((1 + R * R)**-2 - 1)
    Fp = lambda R: 1 - 4 * lam * R * (1 + R * R)**-3
    Fpp = lambda R: -4 * lam * ((1 + R * R)**-3 - 6 * R * R * (1 + R * R)**-4)
    g = lambda R: R * Fp(R) - 2 * Ft(R)
    Rg = np.linspace(1e-3, 10, 100001)
    v = g(Rg)
    roots = [brentq(g, Rg[i], Rg[i + 1], xtol=1e-14) for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]]
    x1 = max(roots)
    n = 2
    vac['Starobinsky n=2 lambda=%g' % lam] = {
        'roots': roots, 'x1': x1, "f''(0)": Fpp(0.0), 'm^2(0) = f\'/(3 f\'\') [R0]': 1 / (3 * Fpp(0.0)),
        "stability f'(R1) > R1 f''(R1)": bool(Fp(x1) > x1 * Fpp(x1)),
        'Starobinsky ineq (1+x^2)^(n+2) > 1+(n+2)x^2+(n+1)(2n+1)x^4': bool((1 + x1**2)**(n + 2) > 1 + (n + 2) * x1**2 + (n + 1) * (2 * n + 1) * x1**4),
        'lambda from his eq (dS) at x1': x1 * (1 + x1**2)**(n + 1) / (2 * ((1 + x1**2)**(n + 1) - 1 - (n + 1) * x1**2))}
# Appleby-Battye, compressed units: vacua R F' = 2 F and their stability (F' > R F'' means m^2 > 0 at a dS point)
for (aA, bA) in [(2.0, 3.0), (10.0, 10.0)]:
    FA = lambda R: R / 2 + np.log(np.cosh(aA * R) - np.tanh(bA) * np.sinh(aA * R)) / (2 * aA)
    FAp = lambda R: 0.5 * (1 + np.tanh(aA * R - bA))
    FApp = lambda R: 0.5 * aA / np.cosh(aA * R - bA)**2
    g = lambda R: R * FAp(R) - 2 * FA(R)
    Rg = np.linspace(1e-6, 3 * bA / aA + 3, 300001)
    v = g(Rg)
    roots = [brentq(g, Rg[i], Rg[i + 1], xtol=1e-14) for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]]
    vac['Appleby-Battye a=%g b=%g' % (aA, bA)] = {
        "F'(0)": FAp(0.0), '1/(1+e^{2b})': 1 / (1 + np.exp(2 * bA)), 'nonzero vacua': roots,
        'stable (F\' > R F\'\') at each': [bool(FAp(r) > r * FApp(r)) for r in roots],
        "G_Mink/G_dS = F'(R_1)/F'(0)": FAp(max(roots)) / FAp(0.0)}
OUT['B3_vacua'] = vac

# ------------------------------------------------------------------------------------------------ B5 background
Ng = np.linspace(np.log(1 / 3), 0, 400)
bg = {}
for fr in [1e-1, 3e-2, 1e-2, 1e-3, 1e-4]:
    mod = HS1(fr)
    row = {}
    for zi in ([10.0, 20.0] if fr >= 1e-3 else [10.0]):
        sol = solve_bg(mod, zi)
        phi, psi = sol.sol(Ng)
        H2 = hs_H2(mod, Ng, phi, psi)
        yH = H2 - np.exp(-3 * Ng)
        w = w_of(Ng, yH)
        fits = {mode: cpl(Ng, w, mode) for mode in ['uniform in N (bound)', 'uniform in a', 'uniform in z']}
        row['z_i=%g' % zi] = {'status': int(sol.status), 'CPL (w0, wa)': fits, 'w(z=0)': float(w[-1]),
                              'w(z=1)': float(np.interp(np.log(0.5), Ng, w)), 'f_R today': float(phi[-1]),
                              'Omega_m true': float(1 / H2[-1])}
    yHa = first_order_yH(mod, Ng)
    row['first-order formula'] = {'CPL uniform in N': cpl(Ng, w_of(Ng, yHa), 'uniform in N (bound)')}
    bg['|f_R0|=%g' % fr] = row
    print('bg', fr, row['z_i=10']['CPL (w0, wa)'], row['z_i=10']['Omega_m true'], flush=True)
OUT['B5_background'] = bg
# w(z) and the phantom crossing
wz = {}
for fr in [1e-2, 1e-4]:
    mod = HS1(fr)
    Nw = np.linspace(-np.log(11), 0, 1600)
    sol = solve_bg(mod, 10.0)
    phi, psi = sol.sol(Nw)
    yH = hs_H2(mod, Nw, phi, psi) - np.exp(-3 * Nw)
    w = w_of(Nw, yH)
    z = np.exp(-Nw) - 1
    sel = z < 4
    idx = np.where(np.sign(w[sel][:-1] + 1) != np.sign(w[sel][1:] + 1))[0]
    wz['|f_R0|=%g' % fr] = {'w(z)': {('z=%g' % zz): float(np.interp(-np.log1p(zz), Nw, w)) for zz in [0, 0.5, 1, 1.5, 2, 3]},
                            'phantom crossing z': [float(z[sel][i]) for i in idx]}
OUT['B5_w_of_z'] = wz
# future attractor
fut = {}
for fr in [1e-2, 1e-4]:
    mod = HS1(fr)
    sol = solve_bg(mod, 10.0, Nf=np.log(200.0), rtol=1e-8)
    Rf = mod.R_of_phi(sol.y[0, -1])
    R1 = vac['HS n=1 |f_R0|=%g' % fr]['R_1']
    fut['|f_R0|=%g' % fr] = {'R(a=200)/R_1 - 1': float(Rf / R1 - 1), '1 + f_R(a=200)': float(1 + sol.y[0, -1])}
OUT['B5_future'] = fut
# small-field scaling and published bounds (first-order formula, own implementation)
sc = {}
for fr in [1e-5, 1e-6]:
    w0, wa = cpl(Ng, w_of(Ng, first_order_yH(HS1(fr), Ng)), 'uniform in N (bound)')
    sc['|f_R0|=%g' % fr] = {'(1+w0)/|f_R0|': (1 + w0) / fr, 'wa/|f_R0|': wa / fr}
OUT['B5_small_field_scaling'] = sc
bounds = {'SVH 2009': 1.3e-4, 'Cataneo+ 2015': 10**-4.79, 'DESI DR1 FS': 10**-4.59, 'BOSS Aviles': 1.53e-5,
          'HS Milky Way': 2e-6, 'JVS 2013': 5e-7, 'Desmond-Ferreira 1sig': 1.4e-8}
bb = {}
for k, fr in bounds.items():
    mod = HS1(fr)
    w0, wa = cpl(Ng, w_of(Ng, first_order_yH(mod, Ng)), 'uniform in N (bound)')
    R0 = 3 + 2 * mod.cc
    lam_simple = np.sqrt(3 * mod.fRR(R0) / (1 + mod.fR(R0)))
    lam_full = np.sqrt(3 * mod.fRR(R0) / (1 + mod.fR(R0) - R0 * mod.fRR(R0)))
    to_Mpch = 2997.92458 / np.sqrt(OMT)          # 1/m in Mpc/h, m = sqrt(Om) H0
    bb[k] = {'|f_R0|': fr, '1+w0': 1 + w0, 'wa': wa, 'Compton 1/m_s [Mpc/h]': lam_simple * to_Mpch,
             'with R f_RR term [Mpc/h]': lam_full * to_Mpch, 'in Mpc (h = 0.674)': lam_simple * to_Mpch / 0.674}
OUT['B5_at_bounds'] = bb
# DESI DR2 diagonal distances (ledger row 8)
desi = {'DESI+CMB': (-0.42, 0.21, 0.21, -1.75, 0.58, 0.58), 'DESI+CMB+PantheonPlus': (-0.838, 0.055, 0.055, -0.62, 0.22, 0.19),
        'DESI+CMB+Union3': (-0.667, 0.088, 0.088, -1.09, 0.31, 0.27), 'DESI+CMB+DESY5': (-0.752, 0.057, 0.057, -0.86, 0.23, 0.20)}
dd = {}
for k, (w0c, sp_, sm_, wac, ap, am) in desi.items():
    w0m, wam = -1.0, 0.0
    sw0 = sp_ if w0m > w0c else sm_
    swa = ap if wam > wac else am
    dd[k] = {'symmetric average (bound)': float(np.hypot((w0m - w0c) / sw0, (wam - wac) / ((ap + am) / 2))),
             'correct-side asymmetric': float(np.hypot((w0m - w0c) / sw0, (wam - wac) / swa))}
OUT['B5_DESI_LCDM_point'] = dd


# growth
def growth_run(mod, k_hMpc, lcdm=False, true_bg=None, zi=30.0):
    Lam = mod.cc / 2
    kk = k_hMpc * 2997.92458 / np.sqrt(OMT)

    def H2_of(N):
        if true_bg is None:
            return np.exp(-3 * N) + Lam / 3
        phi, psi = true_bg.sol(N)
        return hs_H2(mod, N, phi, psi)

    def R_of(N):
        if true_bg is None:
            return 3 * np.exp(-3 * N) + 4 * Lam
        return mod.R_of_phi(true_bg.sol(N)[0])

    def rhs(N, y):
        d, dp = y
        h = 1e-5
        dlnH = (np.log(H2_of(N + h)) - np.log(H2_of(N - h))) / (4 * h)
        if lcdm:
            ge = 1.0
        else:
            R = R_of(N)
            u = (kk * np.exp(-N))**2 * mod.fRR(R) / (1 + mod.fR(R))
            ge = (1 + 4 * u) / ((1 + 3 * u) * (1 + mod.fR(R)))
        return [dp, -(2 + dlnH) * dp + 1.5 * np.exp(-3 * N) / H2_of(N) * ge * d]
    Ni = -np.log1p(zi)
    return solve_ivp(rhs, (Ni, 0.0), [np.exp(Ni), np.exp(Ni)], method='LSODA', rtol=1e-10, atol=1e-14, dense_output=True)


gr = {}
ref = growth_run(HS1(1e-6), 0.1, lcdm=True)
Ne = -np.log(1.5)
for fr in [1e-4, 1e-5, 1e-6, 1e-7, 1.4e-8]:
    mod = HS1(fr)
    row = {}
    for k in [0.02, 0.1, 0.2]:
        s = growth_run(mod, k)
        row['k=%g' % k] = float(s.sol(Ne)[1] / ref.sol(Ne)[1] - 1)
    gr['|f_R0|=%g' % fr] = row
# true background (|f_R0| = 1e-4), and its own LCDM reference with the same true Omega_m: check of the neglected shift
mod = HS1(1e-4)
tb = solve_bg(mod, 30.0)
s_true = growth_run(mod, 0.1, true_bg=tb)
s_appr = growth_run(mod, 0.1)
gr['|f_R0|=1e-4, k=0.1: fD on true f(R) background / fD on GR+Lam background - 1'] = float(s_true.sol(Ne)[1] / s_appr.sol(Ne)[1] - 1)
OUT['B5_growth_fD_ratio_minus_1_at_z0.5'] = gr
print('growth', gr, flush=True)

# ------------------------------------------------------------------------------------------------ B4 exact threshold
c1t, c2t = 0.5, 1.0                       # HS-type f with c1 < 1: f' > 0, f'' > 0 everywhere (compressed units)
ft = lambda R: R - c1t * R / (c2t * R + 1)
fpt = lambda R: 1 - c1t / (c2t * R + 1)**2
Rinv = lambda Kv: (np.sqrt(c1t / (1 - Kv)) - 1) / c2t
q0, K1 = 1.0, 1.0
Kq = lambda q: fpt(0.0) + K1 * (q - q0)
Rq = lambda q: Rinv(Kq(q))
eps = lambda q: (Kq(q) * Rq(q) - ft(Rq(q))) / 2           # mu = 0
deps = lambda q: K1 * Rq(q) / 2                            # exact: eps' = K' R(q)/2
res = {}
for c in [1e-4, 1e-3, 1e-2]:
    q0n = brentq(lambda q: eps(q) + c - deps(q) * q, q0 - 0.3, q0 + 0.3 * (1 - fpt(0.0)), xtol=1e-15)
    mun = deps(q0n)
    delta = 2 * mun / K1
    Rs = np.linspace(0, 5, 51)
    fnew = []
    for R in Rs:
        qR = brentq(lambda q: deps(q) - K1 * R / 2 - mun, q0 - fpt(0.0) + 1e-12, q0 + (1 - fpt(0.0)) - 1e-12, xtol=1e-15)
        fnew.append(Kq(qR) * R - 2 * (eps(qR) + c) + 2 * mun * qR)
    fnew = np.array(fnew)
    exact = ft(Rs + delta) - delta * fpt(0.0) + 2 * (mun * q0 - c)
    first = ft(Rs + 2 * c / (q0 * K1)) - 2 * c / (q0 * K1) * fpt(0.0)
    res['c=%g' % c] = {'mu\' q0/c': mun * q0 / c, 'delta/(2c/(q0 K1))': delta / (2 * c / (q0 * K1)),
                       'max|f_new - exact form|': float(np.max(np.abs(fnew - exact))),
                       'max|f_new - first-order form| / max|f_new - f|': float(np.max(np.abs(fnew - first)) / np.max(np.abs(fnew - ft(Rs))))}
OUT['B4_threshold_linearK_exact'] = res
print('B4', res)

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v5_fR_models.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=float)
