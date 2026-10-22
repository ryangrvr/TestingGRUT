"""v_dissipation, script 7: independent recomputation of every dimensionful number in ../../dissipation/NOTES.md
sections 0, 5.1, 5.2, 5.4, 5.5, 6.2 (closed forms; the formulas are re-derived in v1_symbolic.py, checks C and D).

Inputs (own constants, not read from the track's code):
  M = 2.435e18 GeV, rho_Lambda = 2.52e-47 GeV^4                    BRIEF.md
  gamma = 1/(32 pi^2)                                                KV 1601.00601 eq (5.5) [Q]
  H0 = 67.4 km/s/Mpc, Omega_m = 0.315, t0 = 13.787 Gyr, T_CMB = 2.7255 K, z_eq = 3400    Planck 2018 [R]
  hbar = 6.582119569e-25 GeV s, k_B = 8.617333262e-14 GeV/K, Mpc = 3.0856776e19 km, Gyr = 3.15576e16 s  [R]
  g*: 106.75 (T > m_t), 10.75 at 1 MeV, g*s0 = 3.909                 [R]
  Delta_EW = m_H^2 v^2/8 with m_H = 125.25, v = 246.22 GeV; Delta_QCD = 4.4e-4 GeV^4    ../../equations/NOTES.md sec 7
Track values are copied by hand from ../../dissipation/NOTES.md for the comparison (column 'track').
Output: ../results/v7_numbers.json.  Run: python3 -I v7_numbers.py
"""
import json, os
import numpy as np
from scipy.integrate import quad

M = 2.435e18
rhoL = 2.52e-47
gam = 1 / (32 * np.pi**2)
hbar = 6.582119569e-25
kB = 8.617333262e-14
Mpc_km = 3.0856776e19
Gyr_s = 3.15576e16
H0 = 67.4 / Mpc_km * hbar
Om = 0.315
t0 = 13.787 * Gyr_s / hbar
Tcmb = 2.7255 * kB
gs0 = 3.909
mH, vev = 125.25, 246.22
D_EW = mH**2 * vev**2 / 8
D_QCD = 4.4e-4
z_eq = 3400.0


def rho_rad(T, g):
    return np.pi**2 / 30 * g * T**4


def H_rad(T, g):
    return np.sqrt(rho_rad(T, g) / (3 * M**2))


rows = []


def rec(name, mine, track, note=''):
    rel = (mine - track) / track if track not in (0, None) else None
    rows.append({'quantity': name, 'verifier': mine, 'track': track, 'rel_diff': rel, 'note': note})


# ---------------------------------------------------------------- 5.1
N_L = 9 * M**4 / (4 * gam * rhoL)
rec('N_Lambda = 9M^4/(4 gamma rho_L) [e-folds]', N_L, 9.9e122)
N_std = np.log(M * (106.75 / gs0)**(1 / 3) / Tcmb)
rec('N_std (T=M -> T_CMB, entropy conserved)', N_std, 72.5)
floor = 9 * M**4 / (4 * gam * N_std)
rec('floor 9M^4/(4 gamma N_std) [GeV^4]', floor, 3.4e74)
rec('floor / M^4', floor / M**4, 9.8)
rec('floor / rho_L', floor / rhoL, 1.4e121)
# drain after BBN when rho_V <= rho_r(1 MeV): relative change (4 gamma/9)(rho_V/M^4) N, N = 72.5 (track) or ln(1MeV/T_CMB)
rr1 = rho_rad(1e-3, 10.75)
rec('relative drain after BBN, N = 72.5 (track bound)', 4 * gam / 9 * rr1 / M**4 * N_std, 1.0e-86)
N_bbn = np.log(1e-3 / Tcmb * (10.75 / gs0)**(1 / 3))
rec('relative drain after BBN, N = ln(a0/a_BBN)', 4 * gam / 9 * rr1 / M**4 * N_bbn, 1.0e-86, 'tighter span; bound holds')

# ---------------------------------------------------------------- 5.2 attractor at t0
X = 6 * gam * np.sqrt(3) * M * t0
rhoV_t0 = 9 * M**4 * X**(-2 / 3)
H_t0 = np.sqrt(3) * M * X**(-1 / 3)
N_t0 = (3 * np.sqrt(3) * M * t0 / (4 * np.sqrt(gam)))**(2 / 3)
rec('attractor rho_V(t0) [GeV^4]', rhoV_t0, 2.24e35)
rec('attractor rho_V(t0)/rho_L', rhoV_t0 / rhoL, 8.9e81)
rec('attractor H(t0) [GeV]', H_t0, 0.112)
rec('attractor H(t0)/H0', H_t0 / H0, 7.8e40)
rec('attractor H(t0) t0', H_t0 * t0, 7.4e40)
rec('attractor e-folds N(t0)', N_t0, 1.1e41)
rec('attractor Omega_M = 1/(4N)', 1 / (4 * N_t0), 2.2e-42)
t_H0 = M**2 / (2 * gam * H0**3)
rec('attractor age when H = H0, in 1/H0', t_H0 * H0, 4.5e122)
rec('attractor age when H = H0 [Gyr]', t_H0 * hbar / Gyr_s, 6.6e123)
# observed H0 t0 (retrospective): from data and from flat LCDM integral
OL = 1 - Om
Or = 2.473e-5 / 0.674**2 * (1 + 0.2271 * 3.046)
age_int = quad(lambda aa: 1 / (aa * np.sqrt(Om / aa**3 + Or / aa**4 + OL)), 1e-12, 1)[0]
rec('observed H0 t0 (data)', H0 * t0, 0.950)
rec('flat LCDM H0 t0 integral', age_int, 0.951)

# ---------------------------------------------------------------- 5.4
rec('Delta_EW = mH^2 v^2/8 [GeV^4]', D_EW, 1.19e8)
N_end = 9 * M**4 / (4 * gam * (D_EW + D_QCD + rhoL))
rec('N_end (rho_V,end = D_EW + D_QCD + rho_L)', N_end, 2.1e68)
rec('fractional tuning rho_L/(D_EW + D_QCD)', rhoL / (D_EW + D_QCD), 2.1e-55)
H_end = np.sqrt((D_EW + D_QCD) / (3 * M**2))
rec('H at N_end [GeV]', H_end, 2.6e-15)
rec('radiation fraction 1/(4 N_end)', 1 / (4 * N_end), 1.2e-69)
rec('radiation density at N_end [GeV^4]', (D_EW + D_QCD) / (4 * N_end), 1.4e-61)
rec('rho_r(160 GeV), g*=106.75 [GeV^4]', rho_rad(160, 106.75), 2.3e10)
V = (1e16)**4
dI = -(4 * gam / 9) * V**2 / M**4 * 60
rec('inflaton drain V^1/4=1e16, 60 e-folds [GeV^4]', dI, -2.4e53)
rec('inflaton drain / rho_L', dI / rhoL, -9.5e99)
for T_RH, tv in ((1e16, (1.5e93, 1.5e95)), (1e9, (1.5e37, 1.5e39)), (1e3, (1.5e-11, 1.5e-9)), (160.0, (6.6e-18, 6.6e-16))):
    for kap, trv in zip((1e-3, 1e-2), tv):
        D = gam * kap**2 * H_rad(T_RH, 106.75)**4 / 32
        rec('trace drain T_RH=%g, kappa=%g, /rho_L' % (T_RH, kap), D / rhoL, trv)
D_qcd = gam * 0.3**2 * H_rad(0.2, 61.75)**4 / 32
rec('QCD epoch trace drain (kappa=0.3, T=0.2 GeV, g*=61.75) / rho_L', D_qcd / rhoL, 1.2e-36)
D_qcd2 = gam * 0.3**2 * H_rad(0.2, 17.25)**4 / 32
rec('QCD epoch with g*=17.25 / rho_L', D_qcd2 / rhoL, 1.2e-36, 'g* choice: order of magnitude only')
H_eq_tot = H0 * np.sqrt(2 * Om * (1 + z_eq)**3)
rec('matter era from equality, H_1 = H_eq(total) / rho_L', gam / 24 * H_eq_tot**4 / rhoL, 1.4e-104)
# sharper: R = rho_m/M^2 with rho_m = half the total at equality; integrate exactly over the matter+radiation era
rm_eq = Om * 3 * M**2 * H0**2 * (1 + z_eq)**3
def integrand(Nn):  # N from equality, a/a_eq = e^N ; drain per e-fold (gamma/36) (rho_m/M^2)^2
    return gam / 36 * (rm_eq * np.exp(-3 * Nn) / M**2)**2
D_eq_exact = quad(integrand, 0, 50)[0]
rec('matter era from equality, exact rho_m trace / rho_L', D_eq_exact / rhoL, 1.4e-104, 'track bound is 4x high (uses total H_eq)')
# KV/AHHKM scale
T126 = 1.26e3
rec('H(1.26 TeV)^4 / rho_L (g*=106.75)', H_rad(T126, 106.75)**4 / rhoL, 0.98)

# ---------------------------------------------------------------- 5.5 late-time observables
rho_c = 3 * M**2 * H0**2
rho_m0 = Om * rho_c
rho_g0 = np.pi**2 / 15 * Tcmb**4
R_today = (rho_m0 + 4 * rhoL) / M**2
inj_today = gam / 36 * R_today**2
rec('injection per e-fold today / rho_gamma', inj_today / rho_g0, 1.6e-119)
rec('w+1 today = drain/(3 rho_V)', inj_today / (3 * rhoL), 4e-124)
z = 2e6
R_z = rho_m0 * (1 + z)**3 / M**2
rec('injection per e-fold at z=2e6 / rho_gamma', gam / 36 * R_z**2 / (rho_g0 * (1 + z)**4), 6.7e-109)
# graviton fraction per e-fold = (gamma/36) R^2 / rho_r, R = kappa rho_r / M^2
for T, g, kap, tv in ((1e-3, 10.75, 1e-2, 8.8e-94), (1e16, 106.75, 1e-2, 8.8e-17)):
    rr = rho_rad(T, g)
    rec('graviton fraction per e-fold T=%g, kappa=%g' % (T, kap), gam / 36 * (kap * rr / M**2)**2 / rr, tv)
# physical trace fraction of e+- at T = 1 MeV (high-T expansion, fermions: rho - 3p = g m^2 T^2/48 per species, g=4)
me = 0.511e-3
kap_e = (4 * me**2 * (1e-3)**2 / 48) / rho_rad(1e-3, 10.75)
rec('e+- trace fraction at T=1 MeV (estimate)', kap_e, 1e-2, '[E] the track used 1e-2; physical ~6e-3 (and m/T=0.5 is not << 1)')

out = {'inputs': {'M': M, 'rho_Lambda': rhoL, 'gamma': gam, 'H0_GeV': H0, 't0_GeVinv': t0, 'T_CMB_GeV': Tcmb,
                  'Delta_EW': D_EW, 'Delta_QCD': D_QCD, 'z_eq': z_eq},
       'rows': rows}
bad = [r_ for r_ in rows if r_['rel_diff'] is not None and abs(r_['rel_diff']) > 0.06]
out['rows_differing_by_more_than_6pct'] = [r_['quantity'] for r_ in bad]
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'v7_numbers.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for r_ in rows:
    print('%-62s %12.4g %12.4g %s' % (r_['quantity'], r_['verifier'], r_['track'],
                                     '' if r_['rel_diff'] is None else '%+.3f' % r_['rel_diff']))
print('differ >6%:', out['rows_differing_by_more_than_6pct'])
