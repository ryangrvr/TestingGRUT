#!/usr/bin/env python3
"""gp_critic: cross-track checks at the 7737487 gauge point (own code; no Codex code run).

Each block recomputes a number that is load-bearing for the consolidation or that two sources state differently.
Labels follow the brief: [D] derived, [C] computed, [Q] quoted, [E] estimate.
Usage: python3 critic_checks.py <path to gp_matching/results/release.json> <out.json>
"""
import json
import sys

import numpy as np

MRED = 2.435e18          # GeV, reduced Planck mass
MPL = 1.2209e19          # GeV
GEV_INV_M = 1.97327e-16  # metres per GeV^-1
T0 = 2.348e-13           # GeV, CMB temperature today
YEAR = 3.15576e7

# Stated point [Q]
M_Y, h, mu, kap, e = 1e12, 1e11, 5e11, 1e-3, 1e-7
lam_Y = lam_S = 1.0
lam_X, g_X, lam_phi, v = 1e-3, 0.1, 0.03, 0.924089637
m_X = g_X * v            # m_X = g_X v = 0.0924 GeV (Candidate 1 light point)
H_inf = 4.7e13
Y_Q = 4.8294009e-9       # required net yield per visible entropy [Q]
B_DM = 4.4e-10           # GeV per visible entropy [Q]

out = {}

# ---------------------------------------------------------------- 1. matched coefficients and survival
Sig = h / np.sqrt(2)
Dm = M_Y**4 - 4 * (mu * Sig)**2
alpha6 = mu * kap**2 * Sig / Dm                 # angular X^6 coefficient [D] (7737487 CONSTRUCTION)
beta6 = -M_Y**2 * kap**2 / Dm                    # compulsory radial companion [D]
r = 2 * abs(alpha6) / abs(beta6)                 # = sqrt2 mu h / M_Y^2 < 1 whenever Dm > 0
c6 = alpha6 * MRED**2
# unit-coefficient d=6 rates: Codex b5174e5 PAIR_RESULTS summaries [Q]; our 2e9f296 review gives 1.04e-46 total [Q]
G_coh, G_pair = 3.91958726382564e-50, 1.0399701293581432e-46
t0 = 13.8e9 * YEAR
rate = c6**2 * (G_coh + G_pair)
out['survival'] = dict(
    alpha6_GeV_m2=alpha6, beta6_GeV_m2=beta6, companion_ratio_r=r, abs_beta6_over_2abs_alpha6=1 / r,
    c6_reduced_planck=c6, codex_c6=210683.17611461342,
    rate_per_s=rate, codex_rate_per_s=4.6178968874546585e-36, rate_times_age=rate * t0,
    codex_rate_times_age=2.01107045085444e-18, c6_ceiling_1pct=np.sqrt(0.01 / (t0 * (G_coh + G_pair))),
    note='[D] c6 from mediator parameters; rates per unit coefficient [Q]. Agreement with Codex checks arithmetic only.')

# Robustness of the late X^6 coefficient to the operators enumerated by gp_phase (c = 1, Lambda = M_red) [D]
Lam = MRED
X_core = 1.71  # GeV, |X0| at the reference core centre [Q 7737487]
out['survival']['late_alpha6_shifts'] = {
    'SigmaStar_Y_X3_over_Lambda (dim 5, Y integrated out)': Dm / (mu * kap * M_Y**2 * Lam),
    'SigmaStar_X6_over_Lambda3 (dim 7)': (Sig / Lam**3) / alpha6,
    'SigmaStar_Y2_absX2_over_Lambda (dim 5) at core': X_core**2 / (mu * Lam),
    'YStar_X3_absX2_over_Lambda2 (dim 6) at core, relative kappa shift': X_core**2 / Lam**2,
}

# ---------------------------------------------------------------- 2. remainder law used by sibling tracks
def f_quoted(eps):
    return 0.9 * (eps**(-4 / 3) - 1)

def f_nr(eps):
    return 1.01427 * (1 / eps - 1)   # gp_leftovers NR count [Q]

out['remainder_law'] = {f'{x:g}': dict(quoted=f_quoted(x), nr_count=f_nr(x)) for x in (0.027, 0.03, 0.11, 0.17, 0.216)}

# ---------------------------------------------------------------- 3. wall annihilation temperature
sigma_w = 2 * np.sqrt(2) / 3 * np.sqrt(lam_phi) * v**3
Delta = 1e-22

def T_ann(A, C, gstar):
    t = C * A * sigma_w / Delta
    H = 1 / (2 * t)
    return np.sqrt(H * MRED / np.sqrt(np.pi**2 * gstar / 90))

out['walls'] = dict(sigma_GeV3=sigma_w,
                    codex_min_gstar61p75_A2_C10_MeV=1e3 * T_ann(2, 10, 61.75),
                    min_gstar10p75_A2_C10_MeV=1e3 * T_ann(2, 10, 10.75),
                    central_A0p8_C5_MeV=1e3 * T_ann(0.8, 5, 10.75),
                    central_A0p9_C5_MeV=1e3 * T_ann(0.9, 5, 10.75),
                    central_A0p7_C2_MeV=1e3 * T_ann(0.7, 2, 10.75))

# ---------------------------------------------------------------- 4. longitudinal vector (GMR scaling) and vector mass in inflation
m_A = 6 * e * h
XY_S = (1.119736e15, 8.889516e13, 1.254787e13)   # inflationary |X|,|Y|,|Sigma| [C gp_matching, Q Codex]
mA_inf = np.sqrt(2 * e**2 * (XY_S[0]**2 + 9 * XY_S[1]**2 + 36 * XY_S[2]**2))
rho_over_s_A = B_DM * np.sqrt(m_A * 1e9 / 6e-6) * (H_inf / 1e14)**2   # Omega_A/Omega_cdm = sqrt(m/6e-6 eV)(H/1e14)^2
e_for_mA_eq_H = H_inf / np.sqrt(2 * (XY_S[0]**2 + 9 * XY_S[1]**2 + 36 * XY_S[2]**2))
out['vector'] = dict(m_A_vacuum=m_A, m_A_inflation=mA_inf, m_A_inf_over_H=mA_inf / H_inf,
                     GMR_rho_over_s_GeV=rho_over_s_A, pairs_per_net_charge_if_all_to_XXbar=rho_over_s_A / m_A / Y_Q,
                     e_needed_for_m_A_inf_ge_H=e_for_mA_eq_H,
                     note='[E] GMR scaling assumes a constant Stueckelberg-like mass; pairs per net charge = (rho/s)/m_A / Y_Q '
                          '(one A decay -> one X Xbar pair).')

# ---------------------------------------------------------------- 5. hidden-sector energy, thermal restoration, Y condensate
rel = json.load(open(sys.argv[1]))
L = rel['long_run_pi_over_4']['at']
a_pts = np.array([float(k) for k in L])
EX = np.array([L[k]['E_X'] for k in L]); EY = np.array([L[k]['E_Y'] for k in L]); NX = np.array([L[k]['N_X'] for k in L])
sel = a_pts >= 600
slope = np.polyfit(a_pts[sel], EY[sel], 1)[0]                 # comoving a^4 rho_Y grows linearly (matter vs radiation)
EX_late = EX[a_pts == 3000.0][0]; NX_late = abs(NX[a_pts == 3000.0][0])
rho_rad_com = 3 * (MRED / H_inf)**2                            # a^4 rho_rad in H_inf^4 units (a=1 at release)
Gamma_Y = 3 * kap**2 * M_Y / (256 * np.pi**3)                  # [E, gp_matching]
a_dec = np.sqrt(H_inf / Gamma_Y)
EY_dec = slope * a_dec
gstar, gh = 106.75, 7.0                                        # hidden: X(2), phi(1), D(1), A'(3) relativistic
frac_dec = EY_dec / rho_rad_com
Th_over_T = (frac_dec * gstar / gh)**0.25
s_vis_coef = 2 * np.pi**2 / 45 * gstar
s_h_over_s = (2 * np.pi**2 / 45 * gh * Th_over_T**3) / s_vis_coef
Y_net_uniform = NX_late * rel['Y_per_unit_N']
sigma_per_charge = s_h_over_s / Y_net_uniform
# non-thermal: Y quanta -> 3 X each; comoving Y number = E_Y/(M_Y a) in H_inf units
NY_dec = EY_dec / ((M_Y / H_inf) * a_dec)
pairs_nonthermal = 1.5 * NY_dec / NX_late                      # 3 quanta per Y = 1.5 pairs, per net charge
f_nonthermal = 2 * pairs_nonthermal * m_X * Y_Q / B_DM
# X condensate thermalised at release
T_R = (30 * rho_rad_com * H_inf**4 / (np.pi**2 * gstar))**0.25
rho_X0 = EX[0] * H_inf**4
Th_release = (30 * rho_X0 / (np.pi**2 * 10.0))**0.25          # g_h ~ 10 incl. Y, Sigma relativistic
out['hidden_energy'] = dict(
    EX_over_rho_rad=EX_late / rho_rad_com, EY_slope_per_a=slope, Gamma_Y_GeV=Gamma_Y, a_Y_decay=a_dec,
    EY_over_EX_at_decay=EY_dec / EX_late, EY_over_rho_vis_at_decay=frac_dec, Th_over_T=Th_over_T,
    entropy_per_net_charge_uniform_sign=sigma_per_charge, Y_net_uniform_sign=Y_net_uniform,
    Y_net_over_Y_Q=Y_net_uniform / Y_Q,
    pairs_per_net_charge_nonthermal=pairs_nonthermal, f_extra_nonthermal=f_nonthermal,
    codex_pair_budget_1pct=0.0049296,
    T_R=T_R, Th_if_condensate_thermalised_at_release=Th_release,
    portal_note='[E] SM-hidden portal kappa ~ 9e-10: Gamma ~ kappa^2 T << H at all T > h, so visible reheating cannot '
                'restore the hidden U(1); restoration needs the hidden sector\'s own energy or the release resonance.',
    note='[E] linear extrapolation of [C] gp_matching E_Y(a) from a = 600-3000 to Y decay; uniform-sign net charge '
         '(best case); gh = 7.')

# ---------------------------------------------------------------- 6. separate-universe domain size
gs_R, gs_0 = 106.75, 3.909
a0_over_ar = (gs_R / gs_0)**(1 / 3) * T_R / T0
dom = {}
for a_c in (100.0, 150.0, 300.0):
    com = a_c / H_inf                                         # comoving Hubble radius 1/(aH) = a/H_i in radiation
    T_on = 86.0
    phys_on = com * (gs_R / 95.0)**(1 / 3) * T_R / T_on
    dom[f'a={a_c:g}'] = dict(today_m=com * a0_over_ar * GEV_INV_M, at_onset_86GeV_GeVinv=phys_on)
s86 = 2 * np.pi**2 / 45 * 95.0 * 86.0**3
spacing = (2e5 / (Y_Q * s86))**(1 / 3)
out['domains'] = dict(by_release_a=dom, object_spacing_at_86GeV_GeVinv=spacing, object_r90_GeVinv=87.9,
                      note='[A/E] chaotic sensitivity (1e-9 rad) << inflationary phase spread, so super-horizon patches '
                           'at a ~ 100-300 get uncorrelated charge signs (separate-universe argument).')

# ---------------------------------------------------------------- 7. compulsory companion and the one-field window
# Hubble-held minimum of -cH H^2 R^2 + lam R^4 - b R^6 exists iff b X0^2/lam <= 1/6, X0^2 = cH H^2/(2 lam)  [D]
# with b = |beta6| - 2|alpha6| cos(phi0); eta = |alpha6| X0^2/lam (tracks/affleck_dine normalisation)
X_edge = np.sqrt(lam_X / (6 * abs(beta6)))                     # phi0 = pi/2
eta_edge = abs(alpha6) * X_edge**2 / lam_X                     # = r/12
eps_cap_lin = 4.5 * r / (12 * np.sqrt(1 - r**2))               # max over phi0 of 4.5 eta sin(phi0) at existence edge
D7 = 2 * np.pi**2 / 45 * gstar * (90 / (np.pi**2 * gstar))**0.75
A6 = {1: 12.97, 10: 8.66, 100: 17.5}                           # [Q] tracks/affleck_dine radiation release
win = {}
for cH in (1, 10, 100):
    H_edge = X_edge / np.sqrt(cH / (2 * lam_X))
    Y_edge = A6[cH] * c6 * cH**2.5 * (H_edge / MRED)**3.5 / (4 * lam_X**3 * D7)
    H_gmr = 2e10                                                # [E gp_leftovers] ceiling for 1% after thermalisation
    Y_gmr = Y_edge * (min(H_gmr, H_edge) / H_edge)**3.5
    win[f'c_H={cH}'] = dict(H_edge_GeV=H_edge, kappaX2_over_MY2_at_edge=kap * X_edge**2 / M_Y**2,
                            Y_at_edge_sin1=Y_edge, Y_at_edge_over_Y_Q=Y_edge / Y_Q,
                            Y_at_GMR_ceiling_over_Y_Q=Y_gmr / Y_Q)
c_cap_scaled = 1.61e3 * eta_edge / 0.3                          # affleck_dine c_cap at eta_safe=0.3 rescaled [D/E]
r_needed_case_H = 0.19 / 0.375 / np.sqrt(1 + (0.19 / 0.375)**2)
out['companion'] = dict(
    r=r, X_edge_GeV=X_edge, eta_edge=eta_edge, eta_edge_equals_r_over_12=r / 12,
    eps_cap_linear_SUPERSEDED=eps_cap_lin,
    one_field_ccap_rescaled_cH1_H4p7e13=c_cap_scaled, isocurvature_floor_c_heavy=137.0,
    r_needed_for_eps_0p19_linear=r_needed_case_H, windows=win,
    note='[D] existence edge of the Hubble-held minimum with the phase-independent radial companion. The linear-law '
         'eps cap and the A_6-based yields here ignore the amplitude shift at the edge (y0 up to 2) and are SUPERSEDED '
         'by the numerical release in companion_ad.py (eps ~ 0.10, larger q).')

json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
