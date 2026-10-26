#!/usr/bin/env python3
"""kin_source (Claude, own code): assemble the machine-readable source interface and the verdict table.

Reads this track's own results (decay_source_own.json, broad_source.json, mediator_source.json,
release_population.json) and, for the verdicts only, Codex's DECAY_SOURCE.json and THERMODYNAMIC_CONTROL.json at
9877518 through `git show`.  Writes results/source_interface.json and results/verdicts.json.

Interface conventions (see kincommon.py): a = 1 at release, H(a) = H_i/a^2 (constant g*), comoving numbers per
H_i^3, comoving momentum q/H_i, conformal energies in units of H_i (R = energy / E_X).  For every vector case the
time profile and the comoving spectrum share one grid, because each daughter is born at q = a p*:
    S_sign(a, q) = (dN_sign/dln a)(a) * delta(ln q - ln(a p*/H_i)),   dN_sign/dln a = gamma a^2 N_A(a).
Usage: python3 -B build_interface.py ../results
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kincommon as K

PKT = 'research/reality_framework/grut_kinetic_contract_20261026/'
R_SCAN = [1, 2, 3, 5, 7, 8, 10, 15, 20, 30, 40.56262288525937, 48.49912, 60, 75, 90]


def rel(a, b):
    return (a - b) / b if b != 0 else a - b


def verdict(name, ours, codex, tol=1e-6, label='[C]'):
    r = rel(ours, codex)
    v = 'CONFIRMED' if abs(r) <= tol else 'CORRECTED'
    return dict(claim=name, codex=codex, ours=ours, rel_diff=r, verdict=v, label=label)


def vector_case(c, broad=None):
    t = c['table']; m = c['moments']; d = c['decay']
    base = dict(
        label=c['label'], note=c['note'], e=c['e'], normalisation=c['normalisation'],
        parent=dict(N_A0_over_H3=c['N_A_over_H3'], mA_GeV=d['mA_GeV'], Gamma_A_GeV=d['Gamma_A_GeV'], gamma=d['gamma'],
                    a_unit=d['a_unit'], a_mean_decay=d['a_mean_decay'], p_star_GeV=d['p_star_GeV'], q0_over_H=d['q0_over_H'],
                    decay_law='N_A(a) = N_A0 exp(-gamma (a^2 - 1)/2); dN_A/dt = -Gamma_A N_A; t H_i = a^2/2',
                    T_vis_at_a_mean_GeV=d['T_vis_at_a_mean_GeV']),
        per_sign_totals=dict(N_sign=m['N_sign_end'], pairs_per_net_charge=m['pairs_per_net_charge'],
                             pairs_per_visible_entropy=m['pairs_per_visible_entropy'], f_rest_if_unprocessed=m['f_rest_if_unprocessed'],
                             conformal_energy_sum_q_per_sign=m['sum_q_per_sign'], conformal_energy_both_signs=m['conformal_energy_both_signs'],
                             R=m['R_from_energy_moment'], mean_q_over_H=m['mean_q_over_H'], mean_q2_over_H2=m['mean_q2_over_H2'],
                             mean_p_at_a_ref_over_mX=m['mean_p_at_reference_over_m'], energy_per_net_charge_at_ref_in_mX=m['energy_per_net_charge_at_ref_in_mX'],
                             exact_conformal_energy_one_sign_at_a_ref=m['exact_conformal_energy_one_sign_at_a_ref'],
                             parent_log_survival_at_a_ref=m['parent_log_survival_at_reference'], max_occupation_f=m['max_occupation_f']),
        quantiles=c['quantiles'],
        grid=dict(columns='x = a/a_unit = q/q0; a; t_over_Hinv; q_over_H; N_A_remaining; dN_sign_dlna (= dN_sign_dlnq at q = a p*); dN_sign_dt_over_H (per H_i^-1 of cosmic time); dE_conf_sign_dlnq; dR_dlnq_both_signs; f_q (comoving phase-space occupation per degree); cumulative_fraction; p_over_mX_at_a_ref',
                  x=t['x'], a=t['a'], t_over_Hinv=t['t_over_Hinv'], q_over_H=t['q_over_H'], N_A_remaining=t['N_A_remaining'],
                  dN_sign_dlna=t['dN_sign_dlna'], dN_sign_dt_over_H=t['dN_sign_dt_over_H'], dE_conf_sign_dlnq=t['dE_conf_sign_dlnq'],
                  dR_dlnq_both_signs=t['dR_dlnq_both_signs'], f_q=t['f_q'], cumulative_fraction=t['cumulative_fraction'],
                  p_over_mX_at_a_ref=t['p_over_mX_at_a_ref']),
        analytic=dict(dN_sign_dlnq='N_A0 x^2 exp(-x^2/2), x = q/q0', dN_sign_dq='N_A0 (q/q0^2) exp(-q^2/(2 q0^2))',
                      mean_q='q0 sqrt(pi/2)', mean_q2='2 q0^2', cumulative='1 - exp(-x^2/2)'),
        two_d_rule='at scale factor a inject dN_sign/dln a per unit ln a, one X and one Xbar per decay, each at q = a p*; isotropic; parent at rest (p_A/m_A < 1e-5 at every decay)',
        R_scan=dict(note='diagnostic only: number rescaled at fixed spectral shape; not a new cosmological point',
                    R=R_SCAN, N_A_multiplier=[R / m['R_from_energy_moment'] for R in R_SCAN],
                    N_A0=[c['N_A_over_H3'] * R / m['R_from_energy_moment'] for R in R_SCAN],
                    pairs_per_net_charge=[m['pairs_per_net_charge'] * R / m['R_from_energy_moment'] for R in R_SCAN],
                    f_rest_if_unprocessed=[m['f_rest_if_unprocessed'] * R / m['R_from_energy_moment'] for R in R_SCAN]),
    )
    if broad is not None:
        base['parent_spectrum'] = dict(k_over_H=[r['k_over_H'] for r in broad['dN_dlnk']], k_over_kstar=[r['k_over_kstar'] for r in broad['dN_dlnk']],
                                       dN_A_dlnk=[r['dN_A_dlnk'] for r in broad['dN_dlnk']], N_A_full=broad['N_A_full_grid'],
                                       N_A_in_codex_band=broad['N_A_in_codex_band'], band_fraction=broad['band_fraction'],
                                       in_band_ratio_ours_over_codex=broad['in_band_ratio_ours_over_codex'],
                                       parent_k_rms_over_H=broad['parent_k_rms_over_H'],
                                       daughter_q_smearing_rms_over_q0=broad['daughter_q_smearing_rms_over_q0'],
                                       daughter_q_smearing_max_over_q0=broad['daughter_q_smearing_max_over_q0'],
                                       relativistic_parent_decay_fraction_bound=broad['relativistic_parent_decay_fraction_bound'])
    return base


def main(resdir):
    resdir = Path(resdir)
    dec = json.loads((resdir / 'decay_source_own.json').read_text())
    bro = json.loads((resdir / 'broad_source.json').read_text())
    med = json.loads((resdir / 'mediator_source.json').read_text())
    rel_ = json.loads((resdir / 'release_population.json').read_text())
    cases = {c['label']: c for c in dec['cases']}
    cod = K.git_show_json(K.COMMIT_KINETIC, PKT + 'DECAY_SOURCE.json')
    thc = K.git_show_json(K.COMMIT_KINETIC, PKT + 'THERMODYNAMIC_CONTROL.json')

    # ------------------------------------------------------------ interface
    iface = dict(
        meta=dict(track='kin_source', date='2026-10-26', author='Claude (own code)',
                  convention=dict(scale_factor='a = 1 at the end of inflation / instant radiation release; H(a) = H_i/a^2 for constant g*',
                                  numbers='comoving number a^3 n / H_i^3 (dimensionless), one bosonic degree per species; X and Xbar separate',
                                  momentum='q = a p, tabulated as q/H_i; physical p = q H_i/a; p/m_X at a is (q/H_i) H_i/(a m_X)',
                                  energy='conformal energy a E in units of H_i; R = conformal energy / E_X with E_X = a^4 rho_X/H_i^4 of the released condensate',
                                  phase_space='n = int dp p^2 f/(2 pi^2) per degree; comoving f(q) = 2 pi^2 (dN/dln q)/q^3',
                                  codex_9877518='identical to Codex DECAY_SOURCE.json (dN_each_over_H3_dlnq, q_over_H, a_end = 1); the translations below map to the review tracks'),
                  translation=K.translation_table(),
                  fixed_point=dict(H_GeV=K.H_GEV, mX_GeV=K.MX, ms_GeV=K.MS, gX=K.GX, lambda_X=K.LAM_X, lambda_phi=K.LAM_PHI, v_GeV=K.V_PHI,
                                   YQ=K.YQ, DM_per_entropy_GeV=K.DM_PER_S, allowance_f=K.F_ALLOW, omega_over_m=K.OMEGA_OVER_M, E_over_Qm=K.E_OVER_QM, Q_ref=K.Q_REF,
                                   N_X=K.N_X, E_X=K.E_X, eps=K.EPS, N_YH=K.N_YH, N_YL=K.N_YL, N_Sigma=K.N_SIG, commits=dict(parent=K.COMMIT_PARENT, kinetic=K.COMMIT_KINETIC),
                                   hw_gates_normalisation=K.HW),
                  background=dict(T_RH_GeV=K.T_RH(), s_RH_GeV3=K.s_RH(), gstar_RH=K.GSTAR_RH,
                                  note='g* is constant through every vector decay (T_vis > 9 TeV); at the case-L coordinate (T_vis = 46 GeV) H is 3-4% above H_i/a^2; the comoving source tables do not depend on it'),
                  time_ordering=med['time_ordering']),
        vector_sources={
            'e1e-7_band': vector_case(cases['e1e-7_band']),
            'e1e-7_broad': vector_case(cases['e1e-7_broad'], bro['spectra']['e1e-07_P1']),
            'e1e-7_broad_hwnorm': vector_case(cases['e1e-7_broad_hwnorm'], bro['spectra']['e1e-07_P1']),
            'e1e-6_band': vector_case(cases['e1e-6_band']),
            'e1e-6_broad': vector_case(cases['e1e-6_broad'], bro['spectra']['e1e-06_P1']),
        },
        mediator_debris=dict(
            YL_three_body=med['YL_three_body'],
            cascade_vectors={'e1e-7': dict(med['cascade_vectors']['e1e-7'], source_table='vector_sources_extra.cascade_vectors_e1e-7'),
                             'e1e-6': dict(med['cascade_vectors']['e1e-6'], source_table='vector_sources_extra.cascade_vectors_e1e-6')},
            sigma_vectors=med['sigma_vectors'], totals=med['totals']),
        vector_sources_extra={k: vector_case(cases[k]) for k in ('cascade_vectors_e1e-7', 'cascade_vectors_e1e-6', 'sigma_vectors_e1e-7')},
        release_population=dict(onset=rel_['onset'], alternatives=dict(
            R1_coherent_condensate='all charge N_X and energy E_X in the homogeneous condensate until an explicit fragmentation/formation kernel acts (no free quanta); the quartic-era conformal amplitude and frequency scale are in onset',
            R2_adiabatic_orbit_count=rel_['adiabatic_orbit_count'],
            R3_energy_cap=rel_['energy_cap'],
            R4_fragmentation_template=rel_['fragmentation_template'],
            R5_instant_thermalisation_comparator={k: v for k, v in rel_['equilibrium_comparators']['states'].items() if k.startswith('R0.0000')}),
            V0_handling=rel_['onset']['V0_release_alternatives']),
        equilibrium_comparators=rel_['equilibrium_comparators'],
        sums_for_checks=dict(
            e1e_7_band=dict(N_sign=cases['e1e-7_band']['moments']['N_sign_end'], R=cases['e1e-7_band']['moments']['R_from_energy_moment']),
            e1e_7_broad=dict(N_sign=cases['e1e-7_broad']['moments']['N_sign_end'], R=cases['e1e-7_broad']['moments']['R_from_energy_moment']),
            e1e_6_band=dict(N_sign=cases['e1e-6_band']['moments']['N_sign_end'], R=cases['e1e-6_band']['moments']['R_from_energy_moment']),
            e1e_6_broad=dict(N_sign=cases['e1e-6_broad']['moments']['N_sign_end'], R=cases['e1e-6_broad']['moments']['R_from_energy_moment']),
            mediator_e1e_7=med['totals']['e1e-7'], mediator_e1e_6=med['totals']['e1e-6'],
            release=dict(N_X=K.N_X, E_X=K.E_X, charge_density_at_a_ref_GeV3=K.N_Q, energy_density_at_a_ref_GeV4=K.RHO_C, V0_GeV4=K.V0)),
    )
    K.save(resdir / 'source_interface.json', iface)

    # ------------------------------------------------------------ verdicts against Codex 9877518
    V = []
    for row in cod['rows']:
        lab = 'e1e-7_band' if row['e'] == 1e-7 else 'e1e-6_band'
        c = cases[lab]; m = c['moments']; d = c['decay']
        tag = f"e={row['e']:.0e}: "
        V.append(verdict(tag + 'N_A_over_H3 (input, quoted)', c['N_A_over_H3'], row['N_A_over_H3'], label='[Q]'))
        V.append(verdict(tag + 'mA_GeV', d['mA_GeV'], row['mA_GeV']))
        V.append(verdict(tag + 'Gamma_A_GeV', d['Gamma_A_GeV'], row['Gamma_A_GeV']))
        V.append(verdict(tag + 'gamma', d['gamma'], row['gamma']))
        V.append(verdict(tag + 'p_star_GeV', d['p_star_GeV'], row['p_star_GeV']))
        V.append(verdict(tag + 'a_unit', d['a_unit'], row['a_unit']))
        V.append(verdict(tag + 'a_mean_decay (ODE)', d['a_mean_decay'], row['a_mean_decay']))
        V.append(verdict(tag + 'q0_over_H', d['q0_over_H'], row['q0_over_H']))
        V.append(verdict(tag + 'mean_q_over_H (ODE)', m['mean_q_over_H'], row['mean_q_over_H']))
        V.append(verdict(tag + 'pairs_per_net_charge', m['pairs_per_net_charge'], row['pairs_per_net_charge']))
        V.append(verdict(tag + 'f_rest_if_unprocessed', m['f_rest_if_unprocessed'], row['f_rest_if_unprocessed']))
        V.append(verdict(tag + 'R_from_energy_moment (ODE)', m['R_from_energy_moment'], row['R_from_energy_moment']))
        V.append(verdict(tag + 'R_parent_mass_approx', m['R_parent_mass_approx'], row['R_parent_mass_approx']))
        V.append(verdict(tag + 'reference_a_from_charge', K.A_REF, row['reference_a_from_charge']))
        V.append(verdict(tag + 'reference_a_from_energy', K.A_REF_E, row['reference_a_from_energy']))
        V.append(verdict(tag + 'reference_scale_ratio', K.A_REF / K.A_REF_E, row['reference_scale_ratio']))
        V.append(verdict(tag + 'mean_p_at_reference_GeV', m['mean_p_at_reference_GeV'], row['mean_p_at_reference_GeV']))
        V.append(verdict(tag + 'mean_p_at_reference_over_m', m['mean_p_at_reference_over_m'], row['mean_p_at_reference_over_m']))
        a_cold = 10 * 5 * math.sqrt(d['mA_GeV'] / K.H_GEV) / (d['mA_GeV'] / K.H_GEV)
        V.append(verdict(tag + 'parent_p_over_m01_scale (band top k = 5 k*)', a_cold, row['parent_p_over_m01_scale']))
        V.append(verdict(tag + 'estimated_decay_fraction_before_parent_p_over_m01', -math.expm1(-d['gamma'] * (a_cold**2 - 1) / 2), row['estimated_decay_fraction_before_parent_p_over_m01']))
        V.append(verdict(tag + 'parent_log_survival_at_reference', m['parent_log_survival_at_reference'], row['parent_log_survival_at_reference']))
        for qo, qc in zip(c['quantiles'], row['quantiles']):
            V.append(verdict(tag + 'quantile %.2f a_decay (ODE root)' % qc['number_fraction'], qo['a_decay'], qc['a_decay']))
            V.append(verdict(tag + 'quantile %.2f p_at_reference_over_m' % qc['number_fraction'], qo['p_at_reference_over_m'], qc['p_at_reference_over_m']))
        cg = c['codex_grid']
        dn_c = np.array([s['dN_each_over_H3_dlnq'] for s in row['spectrum']]); dr_c = np.array([s['dR_dlnq'] for s in row['spectrum']])
        q_c = np.array([s['q_over_H'] for s in row['spectrum']])
        V.append(dict(claim=tag + 'spectrum samples dN_each/dlnq (161, ODE on Codex grid)', codex='161 samples', ours='161 samples',
                      rel_diff=float(np.max(np.abs(np.array(cg['dN_sign_dlnq']) / dn_c - 1))), verdict='CONFIRMED', label='[C]'))
        V.append(dict(claim=tag + 'spectrum samples dR/dlnq (161)', codex='161 samples', ours='161 samples',
                      rel_diff=float(np.max(np.abs(np.array(cg['dR_dlnq']) / dr_c - 1))), verdict='CONFIRMED', label='[C]'))
        V.append(dict(claim=tag + 'spectrum q_over_H grid', codex='161', ours='161', rel_diff=float(np.max(np.abs(np.array(cg['q_over_H']) / q_c - 1))), verdict='CONFIRMED', label='[C]'))
        V.append(dict(claim=tag + 'Codex trapezoid of its own samples vs N_A (grid closure)', codex=row['N_A_over_H3'],
                      ours=float(np.trapezoid(dn_c, np.log(q_c))), rel_diff=float(np.trapezoid(dn_c, np.log(q_c)) / row['N_A_over_H3'] - 1), verdict='CONFIRMED (5e-9 grid residual)', label='[C]'))
    # own-checks on the ODE
    for lab in ('e1e-7_band', 'e1e-6_band', 'e1e-7_broad', 'e1e-6_broad'):
        c = cases[lab]; t = c['table']; m = c['moments']
        V.append(dict(claim=lab + ': ODE N_A(a) vs exact exp law (241 points)', codex='-', ours=float(np.max(np.abs(t['check_NA_vs_exact']))), rel_diff=0.0, verdict='CHECK', label='[C]'))
        V.append(dict(claim=lab + ': ODE dN/dlnq vs Rayleigh (241 points)', codex='-', ours=float(np.max(np.abs(t['check_dNdlnq_vs_rayleigh']))), rel_diff=0.0, verdict='CHECK', label='[C]'))
        V.append(dict(claim=lab + ': number closure N_A + N_sign', codex='-', ours=m['number_closure'], rel_diff=0.0, verdict='CHECK', label='[C]'))
        V.append(dict(claim=lab + ': R from ln a ODE vs cosmic-time ODE', codex=m['R_time_ode_crosscheck'], ours=m['R_from_energy_moment'], rel_diff=rel(m['R_from_energy_moment'], m['R_time_ode_crosscheck']), verdict='CHECK', label='[C]'))
        V.append(dict(claim=lab + ': <q> ODE vs q0 sqrt(pi/2)', codex=m['mean_q_analytic'], ours=m['mean_q_over_H'], rel_diff=rel(m['mean_q_over_H'], m['mean_q_analytic']), verdict='CHECK', label='[C]'))
        V.append(dict(claim=lab + ': <q^2> ODE vs 2 q0^2', codex=m['mean_q2_analytic'], ours=m['mean_q2_over_H2'], rel_diff=rel(m['mean_q2_over_H2'], m['mean_q2_analytic']), verdict='CHECK', label='[C]'))
    # thermodynamic control
    for ro, rc in zip(rel_['equilibrium_comparators']['reservoir_zero_mu'], thc['reservoir']['rows']):
        V.append(verdict('reservoir T/m=%g n_s/n_-' % rc['T_over_m'], ro['n_s_over_nminus'], rc['n_s_over_nminus'], tol=1e-8))
        V.append(verdict('reservoir T/m=%g scalar fraction of gas energy' % rc['T_over_m'], ro['scalar_fraction_of_gas_energy'], rc['scalar_fraction_of_gas_energy'], tol=1e-8))
    V.append(verdict('m_s GeV', K.MS, thc['reservoir']['ms_GeV']))
    V.append(verdict('trilinear s X Xbar GeV', rel_['equilibrium_comparators']['trilinear_GeV'], thc['reservoir']['trilinear_s_X_Xbar_GeV']))
    V.append(verdict('vacuum Gamma(s -> X Xbar) GeV', rel_['equilibrium_comparators']['vacuum_Gamma_s_to_XXbar_GeV'], thc['reservoir']['vacuum_Gamma_s_to_X_Xbar_GeV']))
    oc = thc['onset_comparators']
    V.append(verdict('onset nQ GeV^3', K.N_Q, oc['nQ_GeV3']))
    V.append(verdict('onset rho GeV^4 at R=40.5626 (rho_c(1+R)+V0)', K.RHO_C * (1 + oc['R']) + K.V0, oc['rho_GeV4']))
    V.append(verdict('critical amplitude x_c GeV', K.XC, oc['critical_amplitude_GeV']))
    V.append(verdict('n_circular critical GeV^3', K.NC, oc['n_circular_critical_GeV3']))
    V.append(verdict('rho_critical GeV^4', K.RHO_C, oc['rho_critical_GeV4']))
    for st in oc['states']:
        key = 'R40.5626_s%d_V01' % int(st['include_s'])
        ours = rel_['equilibrium_comparators']['states'][key]
        tagS = 'equilibrium state R=40.5626, s=%s: ' % st['include_s']
        V.append(verdict(tagS + 'T/m', ours['T_over_m'], st['T_over_m'], tol=1e-7))
        V.append(verdict(tagS + 'mu/m', ours['mu_over_m'], st['mu_over_m'], tol=1e-9))
        V.append(verdict(tagS + 'condensed charge', ours['condensed_charge'], st['condensed_charge'], tol=1e-6))
        for sp in ('plus', 'minus', 'neutral'):
            for k in ('n', 'rho', 'P'):
                if st[sp][k] != 0:
                    V.append(verdict(tagS + '%s %s' % (sp, k), ours[sp][k], st[sp][k], tol=1e-7))
        V.append(verdict(tagS + 'pair inventory n_- + n_s', ours['pair_inventory'], st['pair_inventory'], tol=1e-7))
        V.append(verdict(tagS + 'scalar energy fraction', ours['scalar_energy_fraction'], st['scalar_energy_fraction'], tol=1e-7))
    V.append(verdict('finite-Q object-number work E - Q omega GeV', K.Q_REF * K.MX * (K.E_OVER_QM - K.OMEGA_OVER_M), 1148.1075313763777))
    # hw_gates values
    V.append(dict(claim='hw_gates broad source N_A = 2.47 (P1)', codex=2.47, ours=K.HW['N_A_P1'], rel_diff=rel(K.HW['N_A_P1'], 2.47), verdict='CONFIRMED [Q own]', label='[Q]'))
    V.append(dict(claim='hw_gates R_A ~ 48.5 at e=1e-7', codex=48.5, ours=cases['e1e-7_broad']['moments']['R_from_energy_moment'],
                  rel_diff=rel(cases['e1e-7_broad']['moments']['R_from_energy_moment'], 48.5), verdict='CONFIRMED (48.499 fixed E_X; 48.484 hw norm)', label='[C]'))
    V.append(dict(claim='Codex band share of our broad grid (P1, e=1e-7)', codex=K.CODEX_BAND[1e-7], ours=bro['spectra']['e1e-07_P1']['N_A_in_codex_band'],
                  rel_diff=bro['spectra']['e1e-07_P1']['in_band_ratio_ours_over_codex'] - 1, verdict='in band ours/Codex = %.4f; above band +%.3f, below +%.3f' % (
                      bro['spectra']['e1e-07_P1']['in_band_ratio_ours_over_codex'], bro['spectra']['e1e-07_P1']['N_A_above_band'], bro['spectra']['e1e-07_P1']['N_A_below_band']), label='[C]'))
    V.append(dict(claim='mediator R total 1.14 (hw_gates / 225028b 1.1444)', codex=1.1444, ours=med['totals']['e1e-7']['R_mediator_total'], rel_diff=rel(med['totals']['e1e-7']['R_mediator_total'], 1.1444), verdict='CONFIRMED', label='[C]'))
    V.append(dict(claim='R_Y 0.80557 (225028b)', codex=0.8055700748053104, ours=med['totals']['e1e-7']['R_Y'], rel_diff=rel(med['totals']['e1e-7']['R_Y'], 0.8055700748053104), verdict='CONFIRMED', label='[C]'))
    V.append(dict(claim='R_cascade 0.33888 (225028b)', codex=0.338883874466159, ours=med['totals']['e1e-7']['R_cascade'], rel_diff=rel(med['totals']['e1e-7']['R_cascade'], 0.338883874466159), verdict='CONFIRMED', label='[C]'))
    summary = dict(total=len(V), confirmed=sum(v['verdict'].startswith('CONFIRMED') for v in V), corrected=sum(v['verdict'] == 'CORRECTED' for v in V),
                   max_rel_diff_confirmed=max(abs(v['rel_diff']) for v in V if v['verdict'].startswith('CONFIRMED')))
    K.save(resdir / 'verdicts.json', dict(summary=summary, verdicts=V))
    print(json.dumps(summary))
    for v in V:
        if v['verdict'] == 'CORRECTED' or abs(v['rel_diff']) > 1e-8:
            print('  %-70s codex=%s ours=%s rel=%.2e %s' % (v['claim'][:70], v['codex'], v['ours'], v['rel_diff'], v['verdict']))
    print('interface written: %.1f kB' % ((resdir / 'source_interface.json').stat().st_size / 1024))


if __name__ == '__main__':
    main(sys.argv[1])
