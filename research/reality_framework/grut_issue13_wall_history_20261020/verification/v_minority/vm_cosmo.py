"""v_minority: independent thermal-history bookkeeping (instant reheating, N_0, re-entry map, relic bound).

Own SM g_rho(T), g_s(T): ideal gases with masses (Bose/Fermi integrals), quark-gluon plasma above T_QCD and a pion
gas below (tanh blend at 160 MeV, width 15 MeV), instantaneous neutrino decoupling at 2 MeV with separate gamma+e
entropy conservation afterwards.  Only the low-T limits (keV) and the high-T limit (1e10 GeV, 106.75) matter for the
numbers checked here; the MeV-GeV band is crude and is not used for any verdict.
"""
import json, math, os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

OUT = os.path.dirname(os.path.abspath(__file__))

MP = 2.435323e18            # reduced Planck mass, GeV
v = 1000.0
lam = 0.01
sigma = (2 * math.sqrt(2 * lam) / 3) * v ** 3     # GeV^3, BRIEF
T0 = 2.7255 * 8.617333262e-14                      # GeV
h_H0 = 0.6736                                       # Planck 2018 (calibration only)
H0 = h_H0 * 100 / 2.99792458e5 / 3.0856775814913673e19 * 6.582119569e-25 * 1e0  # 100 km/s/Mpc in GeV
# 100 km/s/Mpc = 100/(3.0857e19 km) s^-1 ; hbar = 6.582e-25 GeV s
H0 = h_H0 * (100.0 / 3.0856775814913673e19) * 6.582119569e-25
omega_m = 0.1430                                    # Omega_m h^2 (Planck 2018 calibration)
rho_c_over_h2 = 3 * (H0 / h_H0) ** 2 * MP ** 2

# ---------------- g_* model ----------------
def _fx(x, fermion):
    """returns (rho, s) per dof normalized to massless boson (rho: pi^2/30 T^4; s: 2pi^2/45 T^3)."""
    if x > 60:
        return 0.0, 0.0
    sgn = 1.0 if fermion else -1.0
    # rho = (1/2pi^2) ∫ E p^2 /(e^E ± 1) ; p = T q
    fr = lambda q: q * q * math.sqrt(q * q + x * x) / (math.exp(math.sqrt(q * q + x * x)) + sgn)
    fp = lambda q: q ** 4 / math.sqrt(q * q + x * x) / (math.exp(math.sqrt(q * q + x * x)) + sgn) / 3.0
    r = quad(fr, 0, 80, limit=200)[0] / (2 * math.pi ** 2)
    p = quad(fp, 0, 80, limit=200)[0] / (2 * math.pi ** 2)
    s = r + p
    return r / (math.pi ** 2 / 30), s / (2 * math.pi ** 2 / 45)

# (mass GeV, dof, fermion)
SM_EW = [  # always present (leptons, EW bosons, Higgs, photon)
    (0.0, 2, False),          # photon
    (0.000511, 4, True), (0.10566, 4, True), (1.77686, 4, True),
    (80.38, 6, False), (91.19, 3, False), (125.1, 1, False),
]
QUARKS = [(0.0022, 12, True), (0.0047, 12, True), (0.095, 12, True), (1.275, 12, True), (4.18, 12, True),
          (172.7, 12, True), (0.0, 16, False)]   # 3 colours x 2 spins x 2 (q, qbar); gluons 16
HADRONS = [(0.13957, 2, False), (0.13498, 1, False), (0.4937, 4, False), (0.5479, 1, False), (0.7755, 9, False)]
NU = (0.0, 6, True)
T_QCD, W_QCD, T_DEC = 0.160, 0.015, 0.002


def _sum(T, lst):
    r = s = 0.0
    for m, g, fer in lst:
        a, b = _fx(m / T, fer)
        r += g * a; s += g * b
    return r, s


def g_eq(T):
    """all SM species in equilibrium at T (above neutrino decoupling)."""
    r, s = _sum(T, SM_EW)
    wq = 0.5 * (1 + math.tanh((T - T_QCD) / W_QCD))
    rq, sq = _sum(T, QUARKS) if wq > 1e-12 else (0.0, 0.0)
    rh, sh = _sum(T, HADRONS) if wq < 1 - 1e-12 else (0.0, 0.0)
    r += wq * rq + (1 - wq) * rh; s += wq * sq + (1 - wq) * sh
    rn, sn = 6 * 7 / 8, 6 * 7 / 8
    return r + rn, s + sn


_GE_DEC = None
def g_gamma_e(T):
    r = 2.0; s = 2.0
    a, b = _fx(0.000511 / T, True)
    return r + 4 * a, s + 4 * b


def g_star(T):
    """(g_rho, g_s) with instantaneous nu decoupling at T_DEC."""
    global _GE_DEC
    if T >= T_DEC:
        return g_eq(T)
    if _GE_DEC is None:
        # at T_DEC: gamma+e entropy dof (with muons negligible there)
        _GE_DEC = g_gamma_e(T_DEC)[1]
    rge, sge = g_gamma_e(T)
    tnu = (sge / _GE_DEC) ** (1 / 3)
    return rge + 6 * 7 / 8 * tnu ** 4, sge + 6 * 7 / 8 * tnu ** 3

# tabulate
_lT = np.linspace(math.log(1e-10), math.log(1e12), 1101)
_G = np.array([g_star(math.exp(x)) for x in _lT])
def gr(T): return float(np.interp(math.log(T), _lT, _G[:, 0]))
def gs(T): return float(np.interp(math.log(T), _lT, _G[:, 1]))
gs0 = gs(1e-10)


def s_of_T(T): return 2 * math.pi ** 2 / 45 * gs(T) * T ** 3

def rho_rad(T): return math.pi ** 2 / 30 * gr(T) * T ** 4

def a_over_a0(T): return (gs0 / gs(T)) ** (1 / 3) * T0 / T

def H_of_T(T, matter=True):
    rho = rho_rad(T)
    if matter:
        rho += omega_m * rho_c_over_h2 / a_over_a0(T) ** 3
    return math.sqrt(rho / 3) / MP

def nprime(T):
    """e-folds after today's Hubble scale exits for a scale re-entering at T (a R = 1/H with c_R = 1)."""
    return math.log(a_over_a0(T) * H_of_T(T) / H0)

def T_RH(Hinf):
    """instant reheating: 3 MP^2 H^2 = rho_rad(T)."""
    f = lambda lT: math.log(rho_rad(math.exp(lT))) - math.log(3 * MP ** 2 * Hinf ** 2)
    return math.exp(brentq(f, math.log(1e6), math.log(1e14)))

def N0_instant(Hinf):
    T = T_RH(Hinf)
    return math.log(a_over_a0(T) * Hinf / H0)

def T_dom_iso():
    f = lambda lT: math.log(sigma / (MP ** 2 * H_of_T(math.exp(lT))))
    return math.exp(brentq(f, math.log(1e-7), math.log(1e-3)))


if __name__ == '__main__':
    R = {}
    R['inputs'] = dict(MP=MP, sigma=sigma, T0=T0, h=h_H0, H0_GeV=H0, omega_m=omega_m)
    R['g_low'] = dict(g_rho=gr(1e-6), g_s=gs(1e-6), gs0=gs0)
    R['g_high'] = dict(g_rho=gr(1e10), g_s=gs(1e10))
    s0 = s_of_T(T0)
    R['relic_bound_rho_over_s_GeV'] = 0.12 * rho_c_over_h2 / s0
    R['Y_bound_m100'] = R['relic_bound_rho_over_s_GeV'] / 100
    R['Y_bound_m141.42'] = R['relic_bound_rho_over_s_GeV'] / 141.42
    tab = []
    for Hi in [130, 150, 198, 241, 250, 300, 302, 500, 1000, 1600, 2000]:
        T = T_RH(Hi)
        N0 = N0_instant(Hi)
        rho_inf = 3 * MP ** 2 * Hi ** 2
        # kination (w=1) from end of inflation until radiation takes over at T_x = 4 MeV
        rho_x = rho_rad(4e-3)
        N0_kin = N0 + math.log(rho_inf / rho_x) / 12
        # matter-dominated reheating to T_RH = 10 MeV
        rho_rh = rho_rad(1e-2)
        N0_md = N0 - math.log(rho_inf / rho_rh) / 12
        tab.append(dict(H=Hi, T_RH=T, N0_instant=N0, N0_kination_to_4MeV=N0_kin, N0_MD_reheat_10MeV=N0_md,
                        DeltaN_50=50 - N0, DeltaN_60=60 - N0))
    R['history'] = tab
    # direct numerical cross-check of the kination formula (integrate ln(aH) through a w=1 era)
    Hi = 1600.0
    rho_inf = 3 * MP ** 2 * Hi ** 2; rho_x = rho_rad(4e-3)
    # ln(a_end H_end) - ln(a_x H_x) for w: (1+3w)/(6(1+w)) ln(rho_inf/rho_x)
    lx = math.log(a_over_a0(4e-3) * math.sqrt(rho_x / 3) / MP / H0)
    R['kination_direct_1600'] = dict(N0=lx + (1 / 3) * math.log(rho_inf / rho_x),
                                     N0_radiation=lx + 0.25 * math.log(rho_inf / rho_x))
    # re-entry map
    Td = T_dom_iso()
    R['T_dom_iso_keV'] = Td * 1e6
    R['sigmaH_over_s_at_Tdom'] = sigma * H_of_T(Td) / s_of_T(Td)
    R['nprime'] = {lab: nprime(T) for lab, T in [('T_rec_0.256eV', 0.256e-9), ('T_eq_0.80eV', 0.80e-9),
                   ('T_dom_net_6.37keV', 6.37e-6), ('T_dom_iso', Td), ('70keV', 70e-6), ('1MeV', 1e-3)]}
    # T_re examples at 250 GeV, m = 10, 30, 40 (e-folds before end), c_R = 1
    Hi = 250.0
    lend = math.log(a_over_a0(T_RH(Hi)) * Hi / H0)
    ex = {}
    for m in [10, 30, 40]:
        target = lend - m
        f = lambda lT: nprime(math.exp(lT)) - target
        ex[f'm={m}'] = math.exp(brentq(f, math.log(1e-11), math.log(1e12)))
    R['T_re_examples_250GeV'] = ex
    R['g_factor_high_over_low'] = (gr(1e10) ** 0.5 / gs(1e10) ** (1 / 3)) / (gr(1e-6) ** 0.5 / gs(1e-6) ** (1 / 3))
    # sigma/(MP^2 H) at 1 MeV and 42 keV
    R['wall_fraction'] = {k: sigma / (MP ** 2 * H_of_T(T)) for k, T in [('1MeV', 1e-3), ('42keV', 42e-6)]}
    with open(os.path.join(OUT, 'results_cosmo.json'), 'w') as f:
        json.dump(R, f, indent=1)
    print(json.dumps(R, indent=1))
