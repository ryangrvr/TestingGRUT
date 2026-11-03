#!/usr/bin/env python3
"""e2e_sm: N-decay leptogenesis on the e2e_kinetic post-waterfall history (model choices throughout).

Unflavoured Boltzmann equations in cosmic time (H_i units), comoving densities X = n a^3:
  dX_th,j/dt = -G_j <K1/K2>_j (X_th,j - X_eq,j)                                   (thermal N_j, inverse decays)
  dX_L/dt    = sum_j eps_j [G_j <K1/K2>_j (X_th,j - X_eq,j) + G_j (M_j/E_j) X_NT,j]  (decays: thermal + non-thermal)
               - (sum_j G_W,j + G_DL2) X_L                                        (washout)
  G_W,j = G_j z_j^2 K1(z_j)/4  (inverse decays, z = M/T; Buchmuller-Di Bari-Plumacher hep-ph/0401240 form [Q])
  G_DL2 = 0.186 (M/1e10 GeV)(mbar^2/eV^2) H(z=1)/z^3 = 2.62e-29 T^3 mbar^2 GeV  (Delta L = 2, same ref. [Q])
X_NT,j (axially produced / chi-remnant N_j) is read from the e2e_kinetic history: f_Nj rho a^3 / E_j with
E_j = max(3/4 A_j m_parent, M_j, m_parent/2) [E]. T is the e2e_kinetic thermalised-bath temperature. After the
history ends: pure radiation (T ~ 1/a, entropy conserved). Final Y_{B-L} = X_L/(s a^3); Y_B = (28/79) Y_{B-L};
eta_B = 7.04 Y_B. DL=1 scatterings, thermal masses, spectator processes and flavour are not included [E].

Usage: python3 -B sm_lepto.py HIST.json SEESAW.json SCEN H_i_GeV OUT.json
"""
import json, sys, math
import numpy as np
from scipy.special import kn, kve

MP_GEV = 2.435e18
GSTAR = 106.75
SUM_M2 = 2.513e-3 + 7.49e-5      # eV^2 (m1 = 0, NO)
ETA_OBS = 6.12e-10               # Planck 2018 omega_b h^2 = 0.02237 -> eta = 273.8e-10 * omega_b h^2 [R]
PARENT_M = {'sigma': 599.99, 'chi': 249.995}


def K1K2(z):
    z = np.clip(z, 1e-8, 600.0)
    return kve(1, z) / kve(2, z)


def z2K1(z):
    z = np.clip(z, 1e-8, 600.0)
    return z * z * kve(1, z) * np.exp(-z)


def neq(M, T):
    """equilibrium number density of a 2-dof Majorana fermion (Boltzmann), H_i^3 units."""
    z = np.clip(M / np.maximum(T, 1e-30), 1e-8, 600.0)
    T = M / z
    return (2.0 / (2 * math.pi**2)) * M * M * T * kve(2, z) * np.exp(-z)


def load_hist(path):
    d = json.load(open(path))
    h = d['history']
    t = np.array(h['t']); a = np.array(h['a']); H = np.array(h['H'])
    T = np.array([x if x is not None else 0.0 for x in h['T_thermal']])
    fN = np.array([h['f_N%d' % (j + 1)] for j in range(4)])
    A = np.array([[x if x is not None else 0.0 for x in h['A_N%d' % (j + 1)]] for j in range(4)])
    return d, t, a, H, T, fN, A


def run(hist_path, seesaw_path, scen, H_GeV, delta_ref=1e-2, out_path=None, t_after=200.0, washout_scale=1.0,
        nt_scale=1.0):
    d, t, a, Hh, T, fN, A = load_hist(hist_path)
    ss = json.load(open(seesaw_path))
    sc = ss['scenarios'][scen]
    row = [r for r in sc['rows'] if abs(r['delta'] - delta_ref) < 1e-12][0]
    M = np.array(ss['M_over_Hi'])
    # rates and CP asymmetries rescaled to H_GeV (y^2 = mt M/v^2 -> G/H_i ~ H_i; eps ~ M ~ H_i)
    sH = H_GeV / ss['H_i_GeV']
    G = np.array(row['Gamma_rest_over_Hi']) * sH
    eps_unit = np.array(row['eps']) * sH / delta_ref          # eps_j per unit delta (linear regime)
    mpar = np.array([PARENT_M[p] for p in ss['parent']])
    MP = MP_GEV / H_GeV
    rho = 3 * Hh**2 * MP**2
    E = np.maximum(np.maximum(0.75 * A * mpar[:, None], M[:, None]), mpar[:, None] / 2)
    XNT = nt_scale * fN * rho[None, :] * a[None, :]**3 / E
    rNT = G[:, None] * (M[:, None] / E) * XNT                  # NT decays per unit time (comoving)
    lt = np.log(t)
    t_end, a_end, T_end = t[-1], a[-1], T[-1]
    t0 = t[0]
    tmax = t_end * t_after
    # exponential (exact-per-step) integrator on a fine log-time grid: coefficients frozen at step midpoints
    NG = 60000
    ue = np.linspace(math.log(t0), math.log(tmax), NG + 1)
    um = 0.5 * (ue[1:] + ue[:-1]); tm = np.exp(um); dt = np.diff(np.exp(ue))
    inside = tm <= t_end
    am = np.where(inside, np.exp(np.interp(um, lt, np.log(a))), a_end * np.sqrt(tm / t_end))
    Tm = np.where(inside, np.exp(np.interp(um, lt, np.log(np.maximum(T, 1e-12)))), T_end * a_end / am)
    srcm = np.array([np.where(inside, np.interp(um, lt, rNT[j]), 0.0) for j in range(4)])
    z = M[:, None] / Tm[None, :]
    k = K1K2(z)
    GW = 0.25 * G[:, None] * z2K1(z) * washout_scale
    GDL2 = 2.62e-29 * (Tm * H_GeV)**3 * SUM_M2 / H_GeV * washout_scale
    W = GW.sum(0) + GDL2
    Xeq = neq(M[:, None], Tm[None, :]) * am[None, :]**3
    r = G[:, None] * k
    Xth = np.zeros(4); L = 0.0; Lth = 0.0; Wc = 0.0
    Wcum = np.zeros(NG); Ls = np.zeros(NG); Lths = np.zeros(NG)
    for i in range(NG):
        x = r[:, i] * dt[i]
        ex = np.exp(-x)
        dev0 = Xth - Xeq[:, i]
        avgdev = dev0 * np.where(x > 1e-12, (1 - ex) / np.maximum(x, 1e-300), 1.0)
        Sth = np.dot(eps_unit, r[:, i] * avgdev)
        Snt = np.dot(eps_unit, srcm[:, i])
        Xth = Xeq[:, i] + dev0 * ex
        w = W[i] * dt[i]
        ew = math.exp(-w)
        f = (1 - ew) / W[i] if w > 1e-12 else dt[i]
        L = L * ew + (Sth + Snt) * f
        Lth = Lth * ew + Sth * f
        Wc += w
        Wcum[i] = Wc; Ls[i] = L; Lths[i] = Lth
    aa, TT = am[-1], Tm[-1]
    s_com = (2 * math.pi**2 / 45) * GSTAR * TT**3 * aa**3
    YBL_unit = L / s_com
    YBL_th_unit = Lth / s_com
    etaB_unit = 7.04 * (28 / 79) * YBL_unit
    Wtot = Wc
    surv = {}
    for tt in [0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4]:
        if tt < t0 or tt > tmax:
            continue
        j = int(np.searchsorted(tm, tt))
        surv['%g' % tt] = float(math.exp(-(Wtot - Wcum[min(j, NG - 1)])))
    samp = []
    for tt in [1, 10, 100, 1000, 3000, 1e4, 3e4, 1e5, 1e6]:
        if tt < t0 or tt > tmax:
            continue
        j = min(int(np.searchsorted(tm, tt)), NG - 1)
        sc_ = (2 * math.pi**2 / 45) * GSTAR * Tm[j]**3 * am[j]**3
        samp.append(dict(t=float(tt), T_over_Hi=float(Tm[j]), z4=float(M[3] / Tm[j]), W_over_H=float(W[j] * tm[j]),
                         L_over_s_bath_unit=float(Ls[j] / sc_), L_thermal_only_unit=float(Lths[j] / sc_)))
    delta_needed = ETA_OBS / etaB_unit if etaB_unit != 0 else None
    res = dict(label='C on [Q]/[E] rates, model choices; unflavoured, decays + inverse decays + DL=2 only',
               history=hist_path.split('tracks/')[-1], seesaw=scen, H_i_GeV=H_GeV, M_GeV=(M * H_GeV).tolist(),
               Gamma_over_Hi=G.tolist(), eps_per_unit_delta=eps_unit.tolist(),
               etaB_per_unit_delta=etaB_unit, etaB_thermal_only_per_unit_delta=7.04 * (28 / 79) * YBL_th_unit,
               YBL_per_unit_delta=YBL_unit, eta_obs=ETA_OBS, delta_needed=delta_needed,
               eps_at_delta_needed=(eps_unit * delta_needed).tolist() if delta_needed else None,
               etaB_at_delta_0p1=0.1 * etaB_unit, etaB_at_delta_1_linear=etaB_unit,
               washout_survival_from_t=surv, W_total=float(Wtot), samples=samp,
               washout_scale=washout_scale, nt_scale=nt_scale, t_final=tmax, T_final_over_Hi=TT)
    if out_path:
        json.dump(res, open(out_path, 'w'), indent=1)
    return res


if __name__ == '__main__':
    hp, sp, scen, H = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    out = sys.argv[5] if len(sys.argv) > 5 else None
    r = run(hp, sp, scen, H, out_path=out)
    print(json.dumps({k: r[k] for k in ['etaB_per_unit_delta', 'etaB_thermal_only_per_unit_delta', 'delta_needed',
                                        'eps_at_delta_needed', 'washout_survival_from_t', 'W_total']}, indent=1))
    for s in r['samples']:
        print(s)
