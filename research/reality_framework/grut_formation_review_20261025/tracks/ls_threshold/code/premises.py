#!/usr/bin/env python3
"""ls_threshold (own code): order-of-magnitude premise numbers for the case-L benchmark.

  1. Thermal corrections at the onset (one-loop high-T masses of the light potential
     U = lp/4 (phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lX |X|^4): T_c and m_X(T) at the onset T_h.
  2. The on-shell light-scalar route: Gamma(s -> X Xbar), K(T) = n_s/(n_X n_Xbar), and its rate relative to the
     threshold 4->2 rate at equal densities.
  3. Extra light hidden degrees of freedom (s, D) in the onset bath: entropy-per-charge factor.
  4. Timing of the delayed vector decays relative to the light-scale onset (inputs quoted from b448f25 [Q]).
Usage: python3 -B premises.py <out.json>
"""
import json, sys
import numpy as np
import lsrelic as LR

G, LP, LX, V, M, MS = LR.GX, LR.LP, LR.LX, LR.V, LR.M, LR.MS
out = {}

# 1. thermal masses (T in units of m_X)
Pi_phi_full = LP / 4 + G**2 / 6           # phi and X in the bath
Pi_phi_Xonly = G**2 / 6
Tc_full = np.sqrt(LP / Pi_phi_full) * V / M
Tc_Xonly = np.sqrt(LP / Pi_phi_Xonly) * V / M
rows = []
for Th in (11.99, 12.40, 12.68, 13.05):   # T_h0/m at the onset for the four conservative/central ledgers (this track)
    r = {}
    for lab, Tc, PiX in (('phi_and_X_in_bath', Tc_full, LX / 3 + G**2 / 12), ('X_only_bath', Tc_Xonly, LX / 3)):
        v2 = max(1 - (Th / Tc)**2, 0)
        r[lab] = dict(Tc_over_m=Tc, v_T_over_v=np.sqrt(v2), mX_T_over_m=np.sqrt(v2 + PiX * Th**2))
    rows.append(dict(Th_over_m=Th, **r))
out['thermal_masses'] = dict(rows=rows, note='one-loop high-T estimate [E]; corrections fall as T_h^2 and are <1% '
                             'for T_h < 1.5 m, before the 4->2 freeze-out (T_h ~ 0.4 m)')

# 2. s route
lam3 = 2 * G**2 * V
Gs = lam3**2 / (16 * np.pi * MS) * np.sqrt(1 - 4 * M**2 / MS**2)
tab = []
for t in (10, 3, 1, 0.3, 0.1, 0.05, 0.03, 0.02, 0.015, 0.01):
    K = LR.K_s(t) / M**3
    ratio = 4 * LR.SV_SS * K * K / LR.AA      # (route rate)/(A_a 4->2 rate) at n+ ~ n-, small K n
    tab.append(dict(Th_over_m=t, K_GeV_minus3=K, route_over_threshold42=ratio))
out['s_route'] = dict(m_s=MS, Gamma_s_to_XXbar_GeV=Gs, Gamma_over_H_at_onset=Gs / LR.hubble(26.4),
                      sigma_v_ss_GeV2=LR.SV_SS, sigma_v_ss_cm3s=LR.SV_SS * 1.167e-17, table=tab,
                      note='threshold 4->2 coefficients evaluate the s propagator at p^2 = 4 m_X^2 < m_s^2; in a plasma '
                           'pairs with p^2 = m_s^2 exist, giving the resonant (on-shell) route')

# 3. extra dof
out['extra_dof'] = dict(g_X_Xbar=2, g_with_s=3, g_with_s_and_D=4,
                        sigma_factor_with_s=(3 / 2)**0.25, sigma_factor_with_s_and_D=(4 / 2)**0.25,
                        equivalent_R_factor_with_s_and_D=(4 / 2)**(1 / 3),
                        note='at fixed onset energy, s ~ g^(1/4) rho^(3/4); sigma ~ (1+R)^(3/4) so a sigma factor F equals '
                             'a (1+R) factor F^(4/3); this lowers R_L, it cannot rescue')

# 4. vector decay timing [Q inputs: e = 1e-7, h = 1e11, m_A = 6 e h, H_i = 4.348677982066179e11 GeV]
e, h, Hi = 1e-7, 1e11, 4.348677982066179e11
mA = 6 * e * h
GA = e**2 * mA / (48 * np.pi)
a_mean = np.sqrt(np.pi * Hi / (2 * GA))
H_dec = Hi / a_mean**2
from scipy.optimize import brentq
T_dec = brentq(lambda T: LR.hubble(T) - H_dec, 1e-3, 1e8)
out['vector_timing'] = dict(m_A=mA, Gamma_A=GA, mean_a=a_mean, H_at_mean=H_dec, T_vis_at_mean_decay=T_dec,
                            onset_T_vis_raw=[26.42, 'cand10'], note='vectors decay at T_vis ~ TeV, before the light-scale '
                            'onset (T_vis 21-46 GeV), so injecting R at the onset is consistent')
# 5. what an annihilation sink would need: Gamma = n_X sigma v > H when the cannibal p_eq reaches 0.005
for T in (1.3e-8, 1e-6, 1e-4):
    nq = LR.YQ * 5.69 * LR.entropy(T)
    out.setdefault('sink_needed', []).append(dict(T_vis=T, n_Q=nq, sigma_v_needed_GeV2=LR.hubble(T) / nq,
                                                  sigma_v_needed_cm3s=LR.hubble(T) / nq * 1.167e-17))
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
