"""Track m1_transport, item 4 [A/E]: what sets the gas fraction and its velocity, from Codex's own ledgers.

Independent arithmetic on numbers quoted from Codex's packets (labelled [Q]); no Codex code used.
  ed73791 (dynamics): negative charge emitted / Q for 5/10/20/40 % pulses (receiver lower bound)
  bf8ab3b (scattering): direct Xbar/Q0 from the spectral split, 5/10/20/40 %
  c481530 (kinetics): 10 % pulse source 94.859 X+, 1695.010 Xbar; at tau=1000: 11629.389 X+, 0.559 Xbar
  6dc9556 (core matching): through t=20000: N+ = 13189.13, Q_core = 186810.87, gas KE 31.04 m_X at a/a0=10
  94f660d: v_rms = 0.06838 c at 2 MeV, 35.5 pc; 4fe963f: clock 20 MeV, g*=10.75, tau(10)=566135
"""
import json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
Q = 200000.0
pulses = [0.05, 0.10, 0.20, 0.40]
emitted_neg_receiver = [0.00200, 0.00799, 0.0304, 0.0949]     # [Q ed73791]
direct_xbar = [0.002088, 0.008448, 0.03447, 0.1459]           # [Q bf8ab3b]
out = {}
# scaling of the antiparticle emission with pulse amplitude
sl = [math.log(direct_xbar[i+1]/direct_xbar[i])/math.log(pulses[i+1]/pulses[i]) for i in range(3)]
out['xbar_emission_scaling'] = {'pulses': pulses, 'direct_xbar_over_Q0': direct_xbar, 'receiver_negative_over_Q': emitted_neg_receiver,
                                'local_log_slopes_d_ln_Xbar_d_ln_A': sl, 'note': 'Xbar emission ~ A^2.0-2.1 over 5-40 %; nothing in M1 fixes A'}
# catalytic multiplication on cores for the 10 % pulse
src_plus, src_minus = 94.859, 1695.010                        # [Q c481530]
tau1000_plus, tau1000_minus = 11629.389, 0.559                # [Q c481530]
final_plus = 13189.13146074248                                # [Q 6dc9556] (through t = 20000 source window)
Qcore_after_pulse = Q + src_minus - src_plus
out['chain_multiplication_10pct'] = {
    'source_X_plus': src_plus, 'source_Xbar': src_minus, 'Q_core_after_pulse': Qcore_after_pulse,
    'free_X_at_tau_1000': tau1000_plus, 'free_X_final_checkpoint': final_plus,
    'free_X_if_only_one_charge_exchange_per_Xbar': src_plus + src_minus,
    'multiplication_factor_final_over_one_exchange': final_plus / (src_plus + src_minus),
    'core_charge_released_beyond_one_exchange': final_plus - (src_plus + src_minus),
    'free_charge_fraction_final': final_plus / Q, 'free_charge_fraction_one_exchange_only': (src_plus + src_minus) / Q,
    'eta_equivalent_final': 1 - final_plus / Q,
    'note': 'The 6.6 % is 0.85 % emitted antiparticles times a ~7.4x catalytic chain (X-core neutral production, h->X Xbar, Xbar charge exchange) run at the supplied 20 MeV density for tau <= 1000 plus the 10000-20000 source window; the chain is density-clock dependent.'}
# mass-identity check (d4747cb frozen mass test 0.0667805(1-eta) + 2.13356 p <= 0.01)
eta = 1 - final_plus / Q
out['frozen_mass_test_with_this_gas'] = {'eta': eta, 'p': 0.0, 'lhs': 0.0667805 * (1 - eta), 'budget': 0.01,
                                        'codex_stationary_cold_excess_pct_at_checkpoint': 0.8222566687221544, 'codex_total_pct': 0.8944364105298863}
# velocity: gas KE per particle at a/a0=10 and what sets it
KE_x10 = 31.038545373806038                                   # m_X [Q 6dc9556 / 94f660d]
KE_tau1000_x1 = 2526.289                                      # m_X at 20 MeV, 11629 particles [Q c481530]
out['velocity_origin'] = {
    'KE_per_particle_at_2MeV_mX': KE_x10 / final_plus, 'v_rms_from_KE_NR': math.sqrt(2 * KE_x10 / final_plus),
    'codex_v_rms': 0.06837528835889453,
    'KE_per_particle_at_20MeV_tau1000_mX': KE_tau1000_x1 / tau1000_plus, 'v_rms_at_20MeV_NR_estimate': math.sqrt(2 * KE_tau1000_x1 / tau1000_plus),
    'pulse_energy_excess_pct': 2.305, 'pulse_energy_excess_mX': 0.02305 * 187479.9981052459,
    'note': 'v_rms(2 MeV) = v_rms(20 MeV)/10 by redshift; v_rms(20 MeV) ~ 0.66 c is set by the pulse energy shared among the chain-produced particles. Changing the formation clock T_f changes v at 2 MeV as T_f/20 MeV and lambda_fs in proportion.'}
# dependence of the forecast scales on the clock
I_eq_pc_per_unit_v = 580.46                                   # own background (transport_check.json)
mean_v = math.sqrt(8 / (3 * math.pi)) * 0.0684
out['clock_scaling'] = {}
for Tf in (2., 20., 200., 2000.):
    v2 = 0.0684 * Tf / 20.
    lam = mean_v * (Tf / 20.) * I_eq_pc_per_unit_v
    # BGK half-suppression k ~ 5e4/Mpc scales as 20 MeV/T_f; half-mode mass M = (4 pi/3) rho_m (pi/k)^3 ~ (T_f/20 MeV)^3
    out['clock_scaling'][f'T_form_{Tf:g}_MeV'] = {'v_rms_at_2MeV': v2, 'lambda_fs_mean_pc_approx': lam, 'k_half_approx_Mpc^-1_2.7e4_scaled': 2.7e4 * 20. / Tf,
                                                   'M_half_mode_Msun_approx': 4 * math.pi / 3 * 39687737661.0 * (math.pi / (5e4 * 20. / Tf))**3}
# gas fraction scaling with pulse amplitude assuming the same chain factor
for A in pulses:
    xb = direct_xbar[pulses.index(A)]
    out.setdefault('gas_fraction_vs_pulse_same_chain', {})[f'{A:g}'] = {'Xbar_emitted_over_Q': xb, 'f_free_if_same_chain_factor': xb * final_plus / src_minus,
                                                                         'P_suppression_plateau_(1-f)^2_approx': (1 - xb * final_plus / src_minus)**2}
# mass inside the free-streaming length and inside 2pi/k for the BGK nodes
rho_m = 39687737661.0  # Msun/Mpc^3
lam_fs_Mpc = mean_v * I_eq_pc_per_unit_v * 1e-6
out['mass_scales'] = {'lambda_fs_mean_pc': lam_fs_Mpc * 1e6, 'M_in_lambda_fs_sphere_Msun': 4 * math.pi / 3 * rho_m * lam_fs_Mpc**3,
                      'M_in_lambda_fs_cube_Msun': rho_m * lam_fs_Mpc**3,
                      'BGK_half_suppression_k_Mpc^-1_interp': None}
# interpolate the BGK curve to the half-suppression point ((1 + (1-f)^2)/2) in log k
f = 0.058767411685782034
bgk = {1e4: 0.997879, 3e4: 0.979835, 1e5: 0.893210, 3e5: 0.849806, 1e6: 0.764169}
target = (1 + (1 - f)**2) / 2
ks = sorted(bgk); import numpy as np
lk = np.log(ks); r = np.array([bgk[k] for k in ks])
k_half = float(np.exp(np.interp(-target, -r, lk)))
out['mass_scales']['BGK_half_suppression_k_Mpc^-1_interp'] = k_half
out['mass_scales']['BGK_half_suppression_target_ratio'] = target
out['mass_scales']['BGK_half_suppression_M_half_mode_Msun'] = 4 * math.pi / 3 * rho_m * (math.pi / k_half)**3
out['mass_scales']['k_half_thermal_relic_same_vrms_Mpc^-1'] = 27627.76
out['mass_scales']['ratio_BGK_khalf_over_thermal_khalf'] = k_half / 27627.76
# antiparticles in the gas
out['antiparticle_content'] = {'Xbar_per_core_at_tau_1000': 0.559, 'fraction_of_gas_number': 0.559 / tau1000_plus,
                               'Xbar_emitted_in_10000_20000_window_per_core': 1.503, 'h_emitted_in_window_per_core': 3.238,
                               'checkpoint_solver_Nminus_residual': -1.377686108445133e-13,
                               'Gamma_Xbar_core_over_H_at_20MeV_codex': 1.3e6,
                               'note': 'the source is antiparticle emission; the checkpoint gas is >99.99 % positive X because Xbar convert on cores within tau ~ 100 (Gamma/H ~ 1e6); late emission (1.5 Xbar + 3.2 h per core per 1e4/m_X) continues at a low level and is not in the transport.'}
(ROOT / 'results' / 'source_scalings.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
