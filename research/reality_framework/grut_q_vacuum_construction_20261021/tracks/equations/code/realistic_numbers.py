"""Order-of-magnitude numbers for the q-theory formulas at realistic scales.
Inputs (exposed data / recalled values, see NOTES section 7 for labels):
  M = 2.435e18 GeV (reduced Planck mass), rho_L,obs = 2.52e-47 GeV^4 (brief, Planck 2018),
  m_H = 125.25 GeV, v = 246.22 GeV (recalled PDG), QCD bag constant B^(1/4) = 0.145-0.235 GeV (recalled range),
  H0 = 67.4 km/s/Mpc -> 1.44e-42 GeV, Omega_m = 0.315 (recalled Planck 2018).
Natural q-theory scales (estimate [E], NOT derived): q0^2 eps''(q0) ~ M^4, q0 K'(q0) ~ M^2, K(q0) = M^2.
"""
import json, os
import numpy as np

M = 2.435e18
K = M**2
rhoL = 2.52e-47
mH, v = 125.25, 246.22
H0 = 67.4 * 1e3 / 3.0857e22 * 6.582e-25  # s^-1 -> GeV
Om = 0.315
rho_crit = 3 * K * H0**2
rho_m0 = Om * rho_crit
out = {'M_GeV': M, 'K_GeV2': K, 'rhoL_obs': rhoL, 'H0_GeV': H0, 'rho_crit0_GeV4': rho_crit, 'rho_m0_GeV4': rho_m0}

# natural mass of the q-mode: m^2 = (2/3) K eps''/K'^2 with q0^2 eps'' = x M^4, q0 K' = y M^2
for x, y in [(1, 1), (1, 0.1), (1e-4, 1)]:
    m2 = (2 / 3) * K * (x * M**4) / (y * M**2)**2
    out['m_q_GeV (q0^2 eps\'\'=%g M^4, q0 K\'=%g M^2)' % (x, y)] = float(np.sqrt(m2))
m_nat = np.sqrt((2 / 3) * K)

# transition depths
D_EW = mH**2 * v**2 / 8                      # tree-level SM: V(0) - V(v) = lambda v^4/4 = m_H^2 v^2/8
D_QCD = [0.145**4, 0.235**4]
out['Delta_EW_GeV4 (tree level m_H^2 v^2/8)'] = D_EW
out['Delta_QCD_GeV4 (bag range)'] = D_QCD
out['first_order_remnant_EW_over_rhoL (pre-transition mu)'] = -D_EW / rhoL
out['first_order_remnant_QCD_over_rhoL (pre-transition mu)'] = [-d / rhoL for d in D_QCD]
sec_EW = (4 / 3) * D_EW**2 / (K * m_nat**2)
out['second_order_(4/3)Delta_EW^2/(K m_q^2)_GeV4 (m_q=sqrt(2/3)M)'] = sec_EW
out['second_order_EW_over_rhoL'] = sec_EW / rhoL
# m_q for which (4/3) Delta_EW^2/(K m^2) would equal rho_L (illustration only; no mechanism produces it as a constant)
out['m_q_that_would_make_(4/3)D_EW^2/(K m^2)=rhoL_GeV'] = float(np.sqrt((4 / 3) * D_EW**2 / (K * rhoL)))
# adiabaticity of the EW crossover: x = pi m tau/2 with tau ~ 1/H(T=100 GeV)
gstar = 106.75
T = 100.0
rho_rad = np.pi**2 / 30 * gstar * T**4
H_EW = np.sqrt(rho_rad / (3 * K))
out['H(T=100GeV)_GeV'] = H_EW
out['m_q/H_EW (natural m_q)'] = m_nat / H_EW
# matter-induced vacuum energy rho_V ~ rho_m^2/(12 K m^2) today
out['rho_V_matter_induced_today_natural_m_GeV4'] = rho_m0**2 / (12 * K * m_nat**2)
m_needed = np.sqrt(rho_m0**2 / (12 * K * rhoL))
out['m_q_for_rho_m^2/(12 K m^2) = rhoL_GeV'] = m_needed
out['that_m_q_over_H0'] = m_needed / H0
out['that_Yukawa_range_over_Hubble_radius'] = (1 / m_needed) / (1 / H0)
# fifth-force bound (recalled): Yukawa strength alpha = 1/3 requires range lambda <~ 5e-5 m (order of magnitude)
hbarc = 1.97327e-16  # GeV m
out['m_q_lower_bound_from_lambda<50um_eV'] = hbarc / 5e-5 * 1e9
# Cassini |gamma-1| < 2.3e-5 for light q needs exp(-m r) < ~3.5e-5 at r ~ R_sun
Rsun = 6.96e8
out['m_q_lower_bound_from_Cassini_eV (m R_sun > 10.3)'] = 10.3 * hbarc / Rsun * 1e9
# relative shift of K (Newton constant) by matter today: dK/K = rho_m/(3 K m^2)
out['dG/G_today_matter_induced_natural_m'] = rho_m0 / (3 * K * m_nat**2)
# required precision of mu (offset for rho_L): mu - mu0[eps_inf] = -rho_L/q0, relative to mu0 ~ M^4/q0
out['relative_precision_of_mu_for_rhoL (q0 mu0 ~ M^4)'] = rhoL / M**4
for k_, v_ in out.items():
    print(k_, ':', v_)
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'realistic_numbers.json'), 'w') as f:
    json.dump(out, f, indent=1, default=float)
