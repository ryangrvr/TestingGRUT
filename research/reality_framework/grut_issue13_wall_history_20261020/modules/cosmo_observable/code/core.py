"""Core inputs, unit conversions, SM degrees of freedom and the background history with a wall component.

Module cosmo_observable of the issue #13 wall-history packet. numpy/scipy only, deterministic, no random numbers.
Natural units (hbar = c = k_B = 1), energies in GeV, times in GeV^-1 unless a name says otherwise.

Every input carries a provenance tag:
  ACTION   - from the unchanged Candidate 1 action / BRIEF.md
  CALIB    - calibration input taken from published data (stated, retrospective)
  LIT      - number quoted from literature (arXiv id given)
  FREE     - free input of the history, varied explicitly
"""
import math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

# ----------------------------------------------------------------------------------------------------------------
# Physical constants (CODATA 2018 exact/recommended values) and unit conversions
# ----------------------------------------------------------------------------------------------------------------
HBAR_GEV_S = 6.582119569e-25          # GeV s
HBARC_GEV_M = 1.973269804e-16         # GeV m
GEV_J = 1.602176634e-10               # J per GeV (exact)
C_M_S = 299792458.0                   # m/s (exact)
KB_GEV_K = 8.617333262e-14            # GeV per K (exact)
MPC_M = 3.0856775814913673e22         # m
JULIAN_YEAR_S = 365.25 * 86400.0
GYR_S = 1e9 * JULIAN_YEAR_S
KMSMPC_S = 1e3 / MPC_M                # 1 km/s/Mpc in s^-1

M_P = 2.435323e18                     # reduced Planck mass, GeV (ACTION / BRIEF)
G_N = 1.0 / (8.0 * math.pi * M_P**2)  # Newton constant, GeV^-2

# 1 GeV^3 (energy per area in natural units) -> SI
GEV3_J_M2 = GEV_J / HBARC_GEV_M**2            # J/m^2 per GeV^3
GEV3_KG_M2 = GEV3_J_M2 / C_M_S**2             # kg/m^2 per GeV^3 (mass per area)

def kmsmpc_to_gev(h_kmsmpc):
    return h_kmsmpc * KMSMPC_S * HBAR_GEV_S

H100_GEV = kmsmpc_to_gev(100.0)
RHO_C100 = 3.0 * M_P**2 * H100_GEV**2         # critical density for h = 1, GeV^4

# ----------------------------------------------------------------------------------------------------------------
# Action inputs (BRIEF.md; GRUT_CANDIDATE_1.md sections 2-3)
# ----------------------------------------------------------------------------------------------------------------
LAMBDA_PHI = 0.01      # ACTION benchmark
V_BENCH = 1000.0       # GeV, ACTION benchmark
G_X = 0.1              # ACTION benchmark
M_H = math.sqrt(2 * LAMBDA_PHI) * V_BENCH          # 141.42 GeV, radial mode
M_X = G_X * V_BENCH                                # 100 GeV
M_D = math.sqrt((G_X * V_BENCH)**2 + (0.01 * V_BENCH)**2)   # 100.50 GeV (m_D0 = g_X v, h_D = 0.01)

def sigma_tree(v, lam=LAMBDA_PHI):
    """Planar tree tension (2 sqrt(2 lam)/3) v^3 in GeV^3 (BRIEF; Codex DERIVATION.md at 135d535)."""
    return 2.0 * math.sqrt(2.0 * lam) / 3.0 * v**3

SIGMA_BENCH = sigma_tree(V_BENCH)                   # 9.428090416e7 GeV^3

def v_from_sigma(sigma, lam=LAMBDA_PHI):
    return (sigma / (2.0 * math.sqrt(2.0 * lam) / 3.0)) ** (1.0 / 3.0)

# ----------------------------------------------------------------------------------------------------------------
# Calibration inputs (retrospective, stated)
# ----------------------------------------------------------------------------------------------------------------
T_CMB_K = 2.7255                         # CALIB: Fixsen 2009 arXiv:0911.1955 (2.72548 +- 0.00057 K); BRIEF value
T0 = T_CMB_K * KB_GEV_K                  # GeV
OMEGA_B_H2 = 0.02237                     # CALIB: Planck 2018 arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing
OMEGA_C_H2 = 0.1200                      # CALIB: same column
OMEGA_M_H2 = OMEGA_B_H2 + OMEGA_C_H2     # massless neutrinos assumed (sum m_nu = 0.06 eV would add 6.4e-4)
H_PLANCK = 0.6736                        # CALIB (only used to define a Planck-like rho_Lambda variant)
Z_STAR = 1089.92                         # CALIB: Planck 2018 Table 2 z_* (photon decoupling); used only to fix T_rec
AGE_PLANCK_GYR = 13.797                  # retrospective comparison only
H0_MEAS_RANGE = (67.0, 74.0)             # retrospective known data: Planck 2018 67.36+-0.54 (LCDM-inferred),
                                         # Riess et al. 2022 arXiv:2112.04510 73.04+-1.04 (distance ladder)
OMEGA_DM_H2_MAX = 0.12                   # retrospective known: Planck 2018 Omega_c h^2 = 0.1200 +- 0.0012

# Wall scaling parameter (LIT): Hiramatsu-Kawasaki-Saikawa arXiv:1309.5001 eq. (4) A = rho_wall t / sigma_wall,
# eq. (29) A = 0.8 +- 0.1 (radiation-era Z2 lattice simulations; sigma_wall = 2 sqrt(2 lam) eta^3/3, same as ours).
A_CENTRAL, A_LOW, A_HIGH = 0.8, 0.7, 0.9
# Lazanu-Martins-Shellard arXiv:1505.03673 eq. (27): rho = (1 - lambda)(A eta/V) sigma/t, matter era A eta/V = 1.93
A_MATTER_LMS = (1.0 - 2.0 / 3.0) * 1.93     # = 0.643 in the HKS convention rho t / sigma

# ----------------------------------------------------------------------------------------------------------------
# SM effective degrees of freedom: Saikawa & Shirai arXiv:1803.01038 Appendix C fitting functions (C.1)-(C.9)
# ----------------------------------------------------------------------------------------------------------------
_A = [1.0, 1.11724e0, 3.12672e-1, -4.68049e-2, -2.65004e-2, -1.19760e-3, 1.82812e-4, 1.36436e-4, 8.55051e-5,
      1.22840e-5, 3.82259e-7, -6.87035e-9]
_B = [1.43382e-2, 1.37559e-2, 2.92108e-3, -5.38533e-4, -1.62496e-4, -2.87906e-5, -3.84278e-6, 2.78776e-6,
      7.40342e-7, 1.17210e-7, 3.72499e-9, -6.74107e-11]
_C = [1.0, 6.07869e-1, -1.54485e-1, -2.24034e-1, -2.82147e-2, 2.90620e-2, 6.86778e-3, -1.00005e-3, -1.69104e-4,
      1.06301e-5, 1.69528e-6, -9.33311e-8]
_D = [7.07388e1, 9.18011e1, 3.31892e1, -1.39779e0, -1.52558e0, -1.97857e-2, -1.60146e-1, 8.22615e-5, 2.02651e-2,
      -1.82134e-5, 7.83943e-5, 7.13518e-5]
_ME, _MMU, _MPI0, _MPIC = 511e-6, 0.1056, 0.135, 0.140
_M1, _M2, _M3, _M4 = 0.5, 0.77, 1.2, 2.0

def _poly(c, t):
    return sum(ci * t**i for i, ci in enumerate(c))

def _frho(x): return np.exp(-1.04855 * x) * (1 + 1.03757 * x + 0.508630 * x**2 + 0.0893988 * x**3)
def _brho(x): return np.exp(-1.03149 * x) * (1 + 1.03317 * x + 0.398264 * x**2 + 0.0648056 * x**3)
def _fs(x):   return np.exp(-1.04190 * x) * (1 + 1.03400 * x + 0.456426 * x**2 + 0.0595248 * x**3)
def _bs(x):   return np.exp(-1.03365 * x) * (1 + 1.03397 * x + 0.342548 * x**2 + 0.0506182 * x**3)
def _Sfit(x): return 1 + 7.0 / 4.0 * np.exp(-1.0419 * x) * (1 + 1.034 * x + 0.456426 * x**2 + 0.0595249 * x**3)

def gstar_rho(T):
    """g_*rho(T), T in GeV; Saikawa-Shirai (C.1) for T >= 0.12 GeV, (C.3) below."""
    if T >= 0.12:
        t = math.log(T)
        return _poly(_A, t) / _poly(_B, t)
    x = lambda m: m / T
    return (2.030 + 1.353 * _Sfit(x(_ME))**(4.0 / 3.0) + 3.495 * _frho(x(_ME)) + 3.446 * _frho(x(_MMU))
            + 1.05 * _brho(x(_MPI0)) + 2.08 * _brho(x(_MPIC)) + 4.165 * _brho(x(_M1)) + 30.55 * _brho(x(_M2))
            + 89.4 * _brho(x(_M3)) + 8209 * _brho(x(_M4)))

def gstar_s(T):
    """g_*s(T), T in GeV; Saikawa-Shirai (C.2) for T >= 0.12 GeV, (C.4) below."""
    if T >= 0.12:
        t = math.log(T)
        return gstar_rho(T) / (1.0 + _poly(_C, t) / _poly(_D, t))
    x = lambda m: m / T
    return (2.008 + 1.923 * _Sfit(x(_ME)) + 3.442 * _fs(x(_ME)) + 3.468 * _fs(x(_MMU))
            + 1.034 * _bs(x(_MPI0)) + 2.068 * _bs(x(_MPIC)) + 4.16 * _bs(x(_M1)) + 30.55 * _bs(x(_M2))
            + 90 * _bs(x(_M3)) + 6209 * _bs(x(_M4)))

G_RHO0 = gstar_rho(T0)
G_S0 = gstar_s(T0)

def dlngs_dlnT(T, h=1e-4):
    return (math.log(gstar_s(T * math.exp(h))) - math.log(gstar_s(T * math.exp(-h)))) / (2 * h)

def rho_rad(T):
    return math.pi**2 / 30.0 * gstar_rho(T) * T**4

def s_ent(T):
    return 2 * math.pi**2 / 45.0 * gstar_s(T) * T**3

def a_ratio(T):
    """a(T)/a0 from SM entropy conservation (walls are a separate, non-interacting component)."""
    return (G_S0 / gstar_s(T))**(1.0 / 3.0) * T0 / T

def rho_matter(T, omega_m_h2=OMEGA_M_H2):
    return omega_m_h2 * RHO_C100 / a_ratio(T)**3

RHO_R0 = rho_rad(T0)
RHO_M0 = OMEGA_M_H2 * RHO_C100
RHO_L_PLANCK = H_PLANCK**2 * RHO_C100 - RHO_R0 - RHO_M0   # closes a flat no-wall universe at H0 = 67.36

# ----------------------------------------------------------------------------------------------------------------
# Neutrino temperature (instantaneous decoupling, entropy transfer from e+-), for Delta N_eff equivalents
# ----------------------------------------------------------------------------------------------------------------
def g_s_electrons(x):
    """Entropy degrees of freedom of e+ e- (g = 4) at x = m_e/T, in units where s = (2 pi^2/45) g T^3."""
    if x < 1e-8:
        return 7.0 / 2.0
    def integrand(u):  # u = p/T
        e = math.sqrt(u * u + x * x)
        if e > 700:
            return 0.0
        return (u * u * e + u**4 / (3 * e)) / (math.exp(e) + 1.0)
    val, _ = quad(integrand, 0, 60 + 2 * x, limit=200)
    s_over_T3 = 4.0 / (2 * math.pi**2) * val
    return s_over_T3 / (2 * math.pi**2 / 45.0)

def Tnu_over_T(T):
    x = _ME / T
    return ((2.0 + g_s_electrons(x)) / (2.0 + 3.5))**(1.0 / 3.0)

def rho_one_nu(T):
    """Energy density of one neutrino species (nu + nubar), (7/8)(pi^2/15) T_nu^4."""
    tn = Tnu_over_T(T) * T
    return 7.0 / 8.0 * math.pi**2 / 15.0 * tn**4

# ----------------------------------------------------------------------------------------------------------------
# Background history with a wall component
# ----------------------------------------------------------------------------------------------------------------
def hubble(T, rho_w, rho_L=0.0):
    return math.sqrt((rho_rad(T) + rho_matter(T) + rho_w + rho_L) / (3 * M_P**2))

def history(sigma, A=A_CENTRAL, n_post=1.0, f_switch=1.0, rho_L=0.0, T_start=100.0, T_end=None,
            switch_on='rad', freeze_at_acceleration=False, backreaction=True, rtol=1e-10, convention='t'):
    """Integrate t(T) from T_start down to T_end (default T0).

    Phase 1 (scaling, backreaction included): rho_w = A sigma / t with t solved self-consistently from
             dt = dln a / H(T, rho_w). Phase 1 ends when rho_w / rho_X = f_switch, where rho_X = rho_rad
             (switch_on='rad') or rho_rad + rho_m + rho_L (switch_on='all'); if freeze_at_acceleration, phase 1
             also ends where the total deceleration parameter first becomes negative.
    Phase 2 (post-domination): rho_w = rho_w,switch (a_switch/a)^n_post (n_post = 1 frozen, w = -2/3;
             n_post = 2 v_rms^2 = 1/3 network, w = -1/3; 3, 4 conversion stress tests).
    convention='t': rho_w = A sigma / t (HKS eq. 4, cosmic time t); convention='H': rho_w = A_H sigma H with
             A_H = 2A (identical to 't' in an exact radiation era, H = 1/(2t)), H solved self-consistently from
             3 M_P^2 H^2 = rho_others + A_H sigma H. The two differ only through backreaction near domination.
    Returns a dict with switch temperature, H0, age, wall fractions at requested temperatures.
    """
    if T_end is None:
        T_end = T0
    x0, x1 = math.log(T_start), math.log(T_end)
    rho_tot0 = rho_rad(T_start) + rho_matter(T_start) + rho_L
    t_init = 1.0 / (2.0 * math.sqrt(rho_tot0 / (3 * M_P**2)))   # radiation era start (walls negligible there)

    def rw1(T, t):
        if convention == 't':
            return A * sigma / t
        ro = rho_rad(T) + rho_matter(T) + rho_L
        AH = 2.0 * A
        if backreaction:
            H = (AH * sigma + math.sqrt((AH * sigma)**2 + 12 * M_P**2 * ro)) / (6 * M_P**2)
        else:
            H = math.sqrt(ro / (3 * M_P**2))
        return AH * sigma * H

    def rhs1(x, y):
        T = math.exp(x)
        t = math.exp(y[0])
        rw = rw1(T, t)
        H = hubble(T, rw if backreaction else 0.0, rho_L)
        dlna = -1.0 - dlngs_dlnT(T) / 3.0
        return [dlna / (H * t)]

    def ev_switch(x, y):
        T = math.exp(x)
        rw = rw1(T, math.exp(y[0]))
        ref = rho_rad(T) if switch_on == 'rad' else rho_rad(T) + rho_matter(T) + rho_L
        return math.log(rw / ref) - math.log(f_switch)
    ev_switch.terminal = True
    ev_switch.direction = 1

    def ev_accel(x, y):
        T = math.exp(x)
        t = math.exp(y[0])
        rw = A * sigma / t
        H = hubble(T, rw, rho_L)
        w_w = 1.0 / (3.0 * H * t) - 1.0
        return 2 * rho_rad(T) + rho_matter(T) + (1 + 3 * w_w) * rw - 2 * rho_L   # rho + 3p
    ev_accel.terminal = True
    ev_accel.direction = -1

    events = [ev_switch] + ([ev_accel] if freeze_at_acceleration else [])
    sol1 = solve_ivp(rhs1, (x0, x1), [math.log(t_init)], method='DOP853', rtol=rtol, atol=rtol * 1e-2,
                     events=events, dense_output=True)
    out = {'sigma_GeV3': sigma, 'A': A, 'n_post': n_post, 'f_switch': f_switch, 'rho_L_GeV4': rho_L,
           'switch_on': switch_on, 'T_start_GeV': T_start, 'backreaction': backreaction, 'convention': convention,
           'rw1': rw1}
    switched = any(len(e) > 0 for e in sol1.t_events)
    if not switched:
        # scaling all the way to T_end
        t_end = math.exp(sol1.y[0, -1])
        rw_end = rw1(T_end, t_end)
        out.update({'switched': False, 'T_switch_GeV': None, 't_end_GeV-1': t_end, 'rho_w_end': rw_end,
                    'sol1': sol1, 'sol2': None})
        return out
    xs = sol1.t[-1]
    Tsw = math.exp(xs)
    tsw = math.exp(sol1.y[0, -1])
    rw_sw = rw1(Tsw, tsw)
    a_sw = a_ratio(Tsw)

    def rho_w2(T):
        return rw_sw * (a_sw / a_ratio(T))**n_post

    def rhs2(x, y):
        T = math.exp(x)
        t = math.exp(y[0])
        H = hubble(T, rho_w2(T), rho_L)
        dlna = -1.0 - dlngs_dlnT(T) / 3.0
        return [dlna / (H * t)]

    sol2 = solve_ivp(rhs2, (xs, x1), [math.log(tsw)], method='DOP853', rtol=rtol, atol=rtol * 1e-2, dense_output=True)
    out.update({'switched': True, 'T_switch_GeV': Tsw, 't_switch_GeV-1': tsw, 'rho_w_switch': rw_sw,
                'a_switch_over_a0': a_sw, 'sol1': sol1, 'sol2': sol2, 'rho_w2': rho_w2,
                't_end_GeV-1': math.exp(sol2.y[0, -1]), 'rho_w_end': rho_w2(T_end)})
    return out

def wall_density_at(h, T):
    """Wall energy density at temperature T for a history dict returned by history()."""
    x = math.log(T)
    if h['switched'] and T <= h['T_switch_GeV']:
        return h['rho_w2'](T)
    t = math.exp(h['sol1'].sol(x)[0])
    return h['rw1'](T, t)

def time_at(h, T):
    x = math.log(T)
    if h['switched'] and T <= h['T_switch_GeV']:
        return math.exp(h['sol2'].sol(x)[0])
    return math.exp(h['sol1'].sol(x)[0])
