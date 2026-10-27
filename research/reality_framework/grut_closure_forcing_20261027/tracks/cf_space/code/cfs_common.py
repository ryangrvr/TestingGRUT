"""Shared constants and cosmology helpers for track cf_space (closure in space and capacity).

Units: natural units, GeV. Every input carries its source label:
  [Q] quoted from a source read in this pass, [R] recalled, [S] supplied by the GRUT branch packets.
"""
import math

# ---------------------------------------------------------------- physical constants
M_P = 2.435e18            # reduced Planck mass, GeV [R, PDG]
M_PL = 1.22089e19         # non-reduced Planck mass, GeV [R, PDG]
GEV_PER_INV_S = 6.582119569e-25   # hbar in GeV s [R, PDG]
MPC_KM = 3.0856775814913673e19    # km per Mpc [R]
KELVIN_GEV = 8.617333262e-14      # GeV per K [R]
SEC_PER_YR = 3.15576e7

T0_K = 2.7255             # CMB temperature today [R, Fixsen 2009]
T0 = T0_K * KELVIN_GEV    # GeV
G_S0 = 43.0 / 11.0        # entropy dof today (photons + 3 nu, instantaneous decoupling) [R]
G_SM = 106.75             # SM relativistic dof at high T [R]
G_GRUT = 106.75 + 10.0    # + phi(1) + X,Y,Sigma (6 real) + massive hidden vector (3) [D, counting]

# ---------------------------------------------------------------- cosmological inputs (sources in NOTES.md)
H0_PLANCK = 67.36         # km/s/Mpc, Planck 2018 base LCDM TT,TE,EE+lowE+lensing [R; abstract quotes 67.4 +- 0.5]
OMEGA_L_PLANCK = 0.6847   # [R, Planck 2018 VI Table 2]
H0_DESI_CMB = 68.17       # km/s/Mpc, DESI DR2 + CMB, LCDM, arXiv:2503.14738 eq. (21) [Q]
OMEGA_M_DESI_CMB = 0.3027 # arXiv:2503.14738 eq. (21) [Q]
A_S = 2.099e-9            # Planck 2018: ln(1e10 A_s) = 3.044 [R]
K_PIVOT_MPC = 0.05        # Mpc^-1, Planck pivot [R]
R_BK18 = 0.036            # BICEP/Keck 2018, 95% CL, arXiv:2110.00483 [Q via search]
R_TRISTRAM = 0.032        # Planck PR4 + BK18 (Tristram et al. 2022, arXiv:2112.07961) [Q via search]
T_BBN_MIN = 5e-3          # GeV, lowest reheating temperature compatible with BBN (several MeV) [R]

# ---------------------------------------------------------------- GRUT branch inputs [S]
HEAVY = dict(H_i=1.3335e11, M_Y=30 * 1.3335e11, h=40 * 1.3335e11, kappa=3.4587e-3, e=1e-7,
             T_RH_lo=2e14, T_RH_hi=5e14,   # REVIEW_REHEATING.md: reheating within ~1 Hubble time
             c6=1.2504e6, c6_ceiling=1.486e13,   # hw_gates NOTES gate 5 (1% ceiling)
             eps=0.146, lambda_X=1e-3, c_H=10.0)
GAUGEPOINT_H = 4.7e13     # REVIEW_GAUGEPOINT.md / REVIEW_XCHARGE.md H_inf
SCANPOINT_H = (3e11, 6e11)  # REVIEW_SCANPOINT.md regular low-scale branch
CAND1 = dict(g=0.1, v=0.924089637375362, lambda_X=1e-3, lambda_phi=0.03,
             mX=0.09240896373753621, V1loop_over_v4=-1.2265e-6)


def H0_GeV(h0_kms=H0_PLANCK):
    """Hubble constant in GeV."""
    return h0_kms / MPC_KM * GEV_PER_INV_S


def a0H0_per_Mpc(h0_kms=H0_PLANCK):
    """Comoving Hubble wavenumber today, Mpc^-1 (a0 = 1, c = 1)."""
    return h0_kms / 299792.458


def rho_crit(h0_kms=H0_PLANCK):
    return 3.0 * M_P**2 * H0_GeV(h0_kms)**2


def rho_lambda(h0_kms=H0_PLANCK, omega_l=OMEGA_L_PLANCK):
    return omega_l * rho_crit(h0_kms)


def H_lambda(h0_kms=H0_PLANCK, omega_l=OMEGA_L_PLANCK):
    """Asymptotic de Sitter Hubble rate of the observed Lambda."""
    return H0_GeV(h0_kms) * math.sqrt(omega_l)


def S_dS(H):
    """Gibbons-Hawking entropy pi/(G H^2) = 8 pi^2 M_P^2 / H^2 (k_B = 1)."""
    return 8.0 * math.pi**2 * M_P**2 / H**2


def V_of_H(H):
    """Inflationary energy density, Friedmann with slow roll: V ~ 3 M_P^2 H^2."""
    return 3.0 * M_P**2 * H**2


def H_of_V(V):
    return math.sqrt(V / 3.0) / M_P


def r_of_H(H, As=A_S):
    """Tensor-to-scalar ratio given H at horizon exit and the observed A_s: r = 2 H^2/(pi^2 M_P^2 A_s)."""
    return 2.0 * H**2 / (math.pi**2 * M_P**2 * As)


def H_of_r(r, As=A_S):
    return math.pi * M_P * math.sqrt(As * r / 2.0)


def T_of_rho(rho, g=G_SM):
    return (30.0 * rho / (math.pi**2 * g)) ** 0.25


def rho_of_T(T, g=G_SM):
    return math.pi**2 * g * T**4 / 30.0


def s_of_T(T, gs=G_SM):
    return 2.0 * math.pi**2 / 45.0 * gs * T**3


def ln_aH_end_over_a0H0(H_end, rho_end=None, T_rh=None, w=0.0, g_rh=G_SM, h0_kms=H0_PLANCK):
    """ln(a_end H_end / (a_0 H_0)) for: inflation ends at H_end with density rho_end (default 3 M_P^2 H_end^2),
    an intermediate phase with constant equation of state w until reheating at T_rh (default: instant),
    then adiabatic radiation/matter expansion to today with entropy conservation.
    This is the minimal number of inflationary e-folds that puts today's Hubble scale inside the
    horizon at the start of inflation (taking H constant during inflation, conservative)."""
    if rho_end is None:
        rho_end = 3.0 * M_P**2 * H_end**2
    if T_rh is None:
        T_rh = T_of_rho(rho_end, g_rh)
    rho_rh = rho_of_T(T_rh, g_rh)
    if rho_rh > rho_end * (1 + 1e-12):
        raise ValueError("T_rh above the end-of-inflation density")
    ln_a_end_over_a_rh = -math.log(rho_end / rho_rh) / (3.0 * (1.0 + w))
    ln_a_rh_over_a0 = math.log(T0 / T_rh) + math.log(G_S0 / g_rh) / 3.0
    return ln_a_end_over_a_rh + ln_a_rh_over_a0 + math.log(H_end / H0_GeV(h0_kms))


def N_of_k(k_Mpc, H_end, **kw):
    """e-folds before the end of inflation at which comoving k exits (H constant approximation)."""
    h0 = kw.get("h0_kms", H0_PLANCK)
    return ln_aH_end_over_a0H0(H_end, **kw) - math.log(k_Mpc / a0H0_per_Mpc(h0))
