"""Verifier v_reach: independent recomputation of the reach track's cosmological numbers.
Own constants, own sub-10-MeV g*(T) (e+e- annihilation computed from Fermi-Dirac integrals with entropy-conserving
T_nu/T), own free-streaming integral. [R] = recalled, order-of-magnitude reference values only.
Usage: python3 v_cosmo.py <qball_min_summary.json or ''> <out.json>
"""
import json, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

MP = 2.435323e18; GEVCM = 1.973270e-14; GEVG = 1.782662e-24; MSUN = 1.115416e57; MPC = 3.085678e24
T0 = 2.7255 * 8.617333e-14
S0 = 2 * np.pi**2 / 45 * 3.909 * (2.7255 * 8.617333e-5 / 1.973270e-5) ** 3   # cm^-3
RHOC = 1.05368e-5                                                          # GeV cm^-3 h^-2
OMC, OMB = 0.120, 0.02233                                                  # [RETRO], [Q]
YB = OMB * RHOC / (0.93827 * S0)  # baryon yield with proton-mass convention (reach used 8.674e-11)
h0 = 0.674; H0 = 2.13312e-42 * h0; OMM = 0.315
v, gX, lphi, lX = 1000.0, 0.1, 0.01, 0.001
mX = gX * v; mh = np.sqrt(2 * lphi) * v; a_, b_ = gX**2 / lphi, lX / lphi
mstar = b_ ** 0.25 / np.sqrt(a_)
sigma = 2 * np.sqrt(2 * lphi) / 3 * v**3
out = {}

# ---------------- g*(T): computed below 20 MeV, step table above [R]
me = 0.000511


def _fd(m_over_T, kind):
    if kind == 'rho':
        f = lambda y: y**2 * np.sqrt(y**2 + m_over_T**2) / (np.exp(np.sqrt(y**2 + m_over_T**2)) + 1)
        ref = 7 * np.pi**4 / 120
    else:  # entropy: (rho + P)/T
        f = lambda y: y**2 * (np.sqrt(y**2 + m_over_T**2) + y**2 / (3 * np.sqrt(y**2 + m_over_T**2))) / (np.exp(np.sqrt(y**2 + m_over_T**2)) + 1)
        ref = 7 * np.pi**4 / 90
    return quad(f, 0, 60)[0] / ref


def g_low(T):
    ge_r = 3.5 * _fd(me / T, 'rho'); ge_s = 3.5 * _fd(me / T, 's')
    Tnu = ((2 + ge_s) / 5.5) ** (1 / 3) if T < 2e-3 else 1.0   # instantaneous nu decoupling at 2 MeV
    return 2 + ge_r + 5.25 * Tnu**4, 2 + ge_s + 5.25 * Tnu**3


_hi = np.array([[0.02, 10.75], [0.1, 10.75], [0.15, 14.25], [0.2, 17.25], [0.3, 61.75], [1.5, 61.75], [3.0, 72.25],
                [6.0, 75.75], [15, 86.25], [100, 86.25], [180, 96.25], [250, 106.75], [1e20, 106.75]])
_Tl = np.logspace(-6, np.log10(0.02), 120)
_gl = np.array([g_low(T) for T in _Tl])


def gstar(T):
    if T <= 0.02:
        return float(np.interp(np.log(T), np.log(_Tl), _gl[:, 0]))
    return float(np.interp(np.log(T), np.log(_hi[:, 0]), _hi[:, 1]))


def gstars(T):
    if T <= 0.02:
        return float(np.interp(np.log(T), np.log(_Tl), _gl[:, 1]))
    return gstar(T)


def H(T):
    return np.sqrt(np.pi**2 * gstar(T) / 90) * T**2 / MP


def s_of(T):
    return 2 * np.pi**2 / 45 * gstars(T) * T**3


def a_of(T):
    return T0 / T * (3.909 / gstars(T)) ** (1 / 3)


def T_of_a(aa):
    T = T0 / aa
    for _ in range(6):
        T = T0 / aa * (3.909 / gstars(T)) ** (1 / 3)
    return T


def H_of_a(aa):
    T = T_of_a(aa)
    rho_r = np.pi**2 / 30 * gstar(T) * T**4
    rhoc0 = 3 * H0**2 * MP**2
    return np.sqrt((rho_r + rhoc0 * (OMM / aa**3 + (1 - OMM))) / (3 * MP**2))


out['g_check'] = {str(T): [gstar(T), gstars(T)] for T in (3.1e-4, 1e-3, 3.06e-2)}

# ---------------- 1. calibration
coef = YB * S0 / RHOC
EQr = OMC / coef
fam = None
if len(sys.argv) > 2 and sys.argv[1]:
    fam = json.load(open(sys.argv[1]))
e_cases = {'free X (e=1)': 1.0, 'homogeneous condensate at Codex density (e=0.9442) [Q charge_phase]': 0.9442,
           'archived soliton e=0.91109': 0.91108652, 'Q-matter floor e=m*': mstar}
if fam:
    for lab, val in fam['e_at'].items():
        e_cases['Q-ball ' + lab] = val
out['1_calibration'] = dict(YB=YB, s0=S0, coef_per_GeV=coef, EQ_r_GeV=EQr,
                            r_allDM_v1TeV={k: EQr / (e * mX) for k, e in e_cases.items()},
                            omega_asym_per_r_at_mstar=coef * mstar * mX,
                            maser=dict(mX_r=[EQr / 1.0, EQr / mstar], bound=63.93, H0_floor=100 * np.sqrt(OMB + OMC)),
                            T_n_TeV={k: 198.2 * (EQr / (e * mX)) ** (-1 / 3) for k, e in [('e=1', 1.0), ('e=m*', mstar)]})

# ---------------- 2. fragmentation floor (gamma_max = 0.45 x sqrt(lam) v = N H)


def x_of(T, r, vv):
    return r * YB * s_of(T) / (np.sqrt(lphi) * vv**3)


def T_floor(r, vv, N=30.0):
    f = lambda lT: np.log(0.45 * x_of(np.exp(lT), r, vv) * np.sqrt(lphi) * vv / (N * H(np.exp(lT))))
    return float(np.exp(brentq(f, np.log(1e-10), np.log(1e9))))


fl = {}
for vv in (1000.0, 100.0, 89.7, 50.4):
    for lab, e in (('e=1', 1.0), ('e=m*', mstar)):
        r = EQr / (e * gX * vv)
        Tf = T_floor(r, vv)
        xf = x_of(Tf, r, vv)
        fl[f'v={vv} {lab}'] = dict(r=r, T_floor_30H=Tf, x=xf, Q_frag=290.0 * xf ** -0.5 / lphi)
out['2_fragmentation_floor'] = fl

# ---------------- 3. free streaming


def lam_fs(Ti, p_over_m):
    ai = a_of(Ti)

    def integrand(lna):
        aa = np.exp(lna)
        pm = p_over_m * ai / aa
        return pm / np.sqrt(1 + pm**2) / (aa * H_of_a(aa))
    return quad(integrand, np.log(ai), 0.0, limit=500)[0] * GEVCM / MPC


def pm(beta):
    return beta / np.sqrt(1 - beta**2)


fs = {}
for k, row in fl.items():
    fs[k] = {f'beta={bb}': lam_fs(row['T_floor_30H'], pm(bb)) for bb in (0.01, 0.5, 0.9)}
fs['T_Q=4.4e5 GeV'] = {f'beta={bb}': lam_fs(4.4e5, pm(bb)) for bb in (0.5, 0.9)}
deb = {}
for f, T in ((1.0, 0.0306), (0.1, 3.3e-3), (0.01, 3.1e-4)):
    deb[f'f={f}, T_ann={T}'] = lam_fs(T, 1.0 / f)
out['3_free_streaming_Mpc'] = dict(cold_floor=fs, biased_debris_p_over_m_1_over_f=deb,
                                   lya_reference='[R] ~0.03-0.1 Mpc (thermal WDM 3-5 keV)')

# ---------------- 4. GW: HKS/Saikawa amplitude, relic cap, spectral tail into PTA band
eps, A = 0.7, 0.8
Y_DM = OMC * RHOC / S0


def today(T, Om):
    Omr = 2.473e-5 * (3.909 / 2) ** 0 * 1.0  # photons h^2
    # Omega_rad,0 h^2 including neutrinos = 2.473e-5 * (1 + 0.6813) ; g-redshift factor
    return 2.473e-5 * 1.6813 * (gstar(T) / 3.363) * (3.909 / gstars(T)) ** (4 / 3) * Om


def fpk(T):
    return a_of(T) * H(T) / 6.582120e-25


def cap(T, f):
    x = Y_DM / f * s_of(T) / (np.pi**2 / 30 * gstar(T) * T**4)
    return today(T, 3 * eps / (32 * np.pi) * x**2)


def bench(T):
    x = 2 * A * sigma * H(T) / (np.pi**2 / 30 * gstar(T) * T**4)
    return today(T, 3 * eps / (32 * np.pi) * x**2)


def spec(fr, fp, Omp):
    return Omp * ((fr / fp) ** 3 if fr < fp else fp / fr)


gw = dict(cap_coefficient_at_1MeV=cap(1e-3, 1.0), cap_at_T={})
for T in (3.1e-4, 1e-3, 1e-2, 0.1):
    gw['cap_at_T'][str(T)] = dict(cap_f1=cap(T, 1.0), cap_f001=cap(T, 0.01), fpeak=fpk(T))
gw['benchmark_v1TeV'] = {str(T): dict(Om=bench(T), fpeak=fpk(T), cap_f=cap(T, f)) for f, T in ((1.0, 0.0306), (0.1, 3.3e-3), (0.01, 3.1e-4))}
# max in-band amplitude over T_ann >= 0.3 MeV (relic cap, any v) at PTA frequencies, including the f^-1 tail
inband = {}
Ts = np.logspace(np.log10(3e-4), 1, 400)
for f in (1.0, 0.1, 0.01):
    for fr in (2e-9, 1e-8, 3e-8):
        vals = [spec(fr, fpk(T), cap(T, f)) for T in Ts]
        i = int(np.argmax(vals))
        inband[f'f={f}, freq={fr}'] = dict(max_Omega_h2=vals[i], at_T_ann=Ts[i], peak_only_cap=max(cap(T, f) for T in Ts if 0.5 * fr < fpk(T) < 2 * fr) if any(0.5 * fr < fpk(T) < 2 * fr for T in Ts) else None)
gw['in_PTA_band'] = inband
gw['benchmark_f001_tail_at_2e-9Hz'] = spec(2e-9, fpk(3.1e-4), bench(3.1e-4))
gw['reference'] = '[R] NANOGrav-15 ~1e-9 (Omega h^2 near 3e-9 Hz); SKA-era PTA forecasts quoted in the literature range roughly 1e-16 to 1e-13 (uncertain)'
out['4_GW'] = gw

# ---------------- 5. lensing: causal ceiling at the fragmentation floor
caus = {}
for k, row in fl.items():
    if not k.startswith('v=1000'):
        continue
    T = row['T_floor_30H']
    Mh = Y_DM * s_of(T) * 4 * np.pi / 3 / H(T) ** 3
    caus[k] = dict(T=T, M_DM_in_horizon_Msun=Mh / MSUN, M_frag_Msun=row['Q_frag'] * mstar * mX / MSUN)
out['5_lensing_ceilings'] = caus

# ---------------- 6. SIDM, PBH, r_t, portal, M_cut
ell = 1 / (np.sqrt(lphi) * v) * GEVCM
sidm = {}
for lab, Q, e, R in (('Q_b', 1.77e4, 1.0, 5.21), ('archived', 39944, 0.91109, 5.1034)):
    M = e * mX * Q * GEVG
    sidm[lab] = np.pi * (2 * R * ell) ** 2 / M
MXX = abs(4 * gX**4 / lphi - 4 * lX)
sXX = MXX**2 / (32 * np.pi * 4 * mX**2) * GEVCM**2 / (mX * GEVG)
Mhh = 24 * lphi
shh = Mhh**2 / (32 * np.pi * 4 * mh**2) * GEVCM**2 / (mh * GEVG)
out['6_misc'] = dict(sigma_over_m_contact=sidm, sigma_over_m_XX=sXX, sigma_over_m_hh=shh,
                     PBH_Mmin_Msun=4 * np.pi * MP**4 / sigma / MSUN,
                     r_t={str(Hi): 2 * Hi**2 / (np.pi**2 * MP**2 * 2.1e-9) for Hi in (210.0, 560.0)},
                     kappa_th=np.sqrt(H(1.56 * v) / (2e-3 * 1.56 * v)),
                     Mcut_Msun_lambda6em6Mpc=4 * np.pi / 3 * OMM * 2.775e11 * h0**2 * (3e-6) ** 3,
                     Neff_full_conversion_T1TeV=(4 / 3) * EQr * YB / 1e3 * (106.75 / 2) * (3.909 / 106.75) ** (4 / 3) / (7 / 8 * (4 / 11) ** (4 / 3)),
                     rho_ratio_1MeV=Y_DM * s_of(1e-3) / (np.pi**2 / 30 * gstar(1e-3) * 1e-12))
json.dump(out, open(sys.argv[-1], 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
