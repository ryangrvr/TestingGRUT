"""Shared inputs for track gp_leftovers (own code; values quoted from the packets as labelled)."""
import numpy as np

# Light target (7737487 / b5174e5 stated point) [Q]
G = 0.1            # g_X
LX = 1e-3          # lambda_X
LP = 0.03          # lambda_phi
V = 0.924089637375362  # v (GeV), calibrated at Q=200000 (portal-debris FIXED_YIELD)
HD = 0.01          # h_D  (portal-debris widths.py) [Q]
MX = G * V
MS = np.sqrt(2 * LP) * V
MD = np.sqrt(G**2 + HD**2) * V   # m_D0 = 0.1 v [Q]
Q_REF = 200000
OMEGA_REF = 0.8611951308217225 * MX      # refined profile [Q]
E_PER_FREE_MASS = 0.9233161353726057     # E(Q)/(Q m_X), refined [Q]
KAPPA = 9.006401966916295e-10            # joint-screen Higgs portal [Q]
VH, MH = 246.0, 125.0

# Abundance calibration [Q]
YQ = 4.8294009076988175e-9
BDM = 4.4e-10      # GeV per visible entropy
B_PER_CHARGE = BDM / YQ   # GeV

# Gauge model (7737487) [Q]
E_GAUGE = 1e-7
H_SIGMA = 1e11
M_AP = 6 * E_GAUGE * H_SIGMA   # vector mass 6 e h
M_Y = 1e12

MRED = 2.435e18
MPL = 1.22089e19
HBAR = 6.582119569e-25
GEV2_CM3S = 1.167e-17   # 1 GeV^-2 * c in cm^3/s  (0.3894 mb * 2.9979e10 cm/s)


def gstar(T):
    """Coarse SM g* (= g*s) step table, adequate at the 10% level (same table as tracks/affleck_dine)."""
    tab = [(1e-6, 3.36), (5e-4, 3.36), (1e-3, 10.75), (0.1, 10.75), (0.15, 17.25), (0.3, 61.75),
           (1.0, 72.25), (5.0, 86.25), (100.0, 96.25), (200.0, 106.75), (1e19, 106.75)]
    Ts, gs = zip(*tab)
    return float(np.interp(np.log(T), np.log(Ts), gs))


def gstar_s(T):
    g = gstar(T)
    return 3.91 if g < 3.5 else g


def entropy(T):
    return 2 * np.pi**2 / 45 * gstar_s(T) * T**3


def hubble(T):
    return np.sqrt(np.pi**2 * gstar(T) / 90) * T**2 / MRED
