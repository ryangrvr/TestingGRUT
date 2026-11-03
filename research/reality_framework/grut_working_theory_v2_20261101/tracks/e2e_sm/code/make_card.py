#!/usr/bin/env python3
"""Assemble results/sm_card.json (every SM-side parameter of M2 with its status) from seesaw/lepto/checks."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
ss = json.load(open(os.path.join(RES, 'seesaw_H1.3e11.json')))
lp = json.load(open(os.path.join(RES, 'lepto.json')))['cases']
ck = json.load(open(os.path.join(RES, 'checks.json')))
S1 = ss['scenarios']['S1_mt050']
r0 = [r for r in S1['rows'] if r['delta'] == 0.0][0]
r1 = [r for r in S1['rows'] if r['delta'] == 0.01][0]
base = lp['S1_HFa_1.3e11']
P = lambda v, st, note='', **kw: dict(value=v, status=st, note=note, **kw)

card = dict(
    label='e2e_sm parameter card for M2 (SM side). Status: supplied = model choice; fitted = to a named datum; '
          'derived/computed. Nothing here is fixed by closure; nothing is a GRUT prediction.',
    H_i_reference_GeV=P(1.3e11, 'supplied (e2e_kinetic history scale); M2 points at 8.4-8.9e10 rescale',
                        'rates ~ H_i; see lepto case S1_HFa_rescaled_8.4e10'),
    light_neutrinos=dict(
        ordering=P('normal', 'supplied', 'IO needs sum m >= 0.0995 eV, above DESI DR2+CMB 0.0642 (retrospective)'),
        m1_eV=P(0.0, 'supplied', 'window: m1 <~ 4.2e-3 eV keeps sum m < 0.0642 eV (sum = 0.0641 at m1 = 4.2 meV)'),
        dm21_eV2=P(7.49e-5, 'fitted', 'NuFIT 6.0 IC24 (arXiv:2410.05380)'),
        dm31_eV2=P(2.513e-3, 'fitted', 'NuFIT 6.0 IC24 NO (arXiv:2410.05380)'),
        masses_eV=P(ss['light_masses_eV'], 'derived'),
        sum_mnu_eV=P(ss['sum_mnu_eV'], 'derived', 'Planck 2018+BAO < 0.12 (1807.06209); DESI DR2+CMB < 0.0642 '
                     '(2503.14744, LCDM, 95%); retrospective: below both'),
        PMNS=P(dict(th12_deg=33.68, th13_deg=8.56, th23_deg=43.3, delta_deg=212, majorana='0, 0 (irrelevant, '
                    'unflavoured)'), 'fitted (NuFIT 6.0 IC24) / supplied (Majorana)',
               'drop out of h^dag h in the unflavoured regime [D]')),
    singlets=dict(
        content=P('four axial Majorana N1..N4 (mf_hybrid), masses -(M_j/2) Nbar e^{i g5 phi_j/f_j} N, + massless N_s',
                  'supplied'),
        M_over_Hi=P(ss['M_over_Hi'], 'supplied (mf_hybrid)'), M_GeV=P(ss['M_GeV'], 'derived at H_i ref'),
        f_over_Hi=P([1200, 2000, 2000, 2000], 'supplied (mf_hybrid)', 'UV: see uv block'),
        mtilde_eV=P(r0['mtilde_eV'], 'supplied (scenario S1)', 'all >= 0.05 eV: mf_hybrid rates and the '
                    'mf_higgs-skeptic N-gas medium mass are kept or raised (N2 x1.82)'),
        CP_parities_at_delta0=P(ss['ETA'], 'supplied', 'CP-conserving R requires sum_j eta_j mt_j = sum m_i [D]'),
        CI_R0_re=P(r0['R_re'], 'computed (least-norm solution of the m-tilde conditions)'),
        CI_R0_im=P(r0['R_im'], 'computed'),
        CI_A_odd=P(ss['A_ODD'], 'supplied (CP-odd direction; R = R0 exp(delta A))'),
        tuning_sum_mt_over_sum_m=P(S1['tuning_sum_mt_over_sum_m'], 'derived', 'cancellation in m_nu at the x4 level'),
        y2=P(r0['y2_hdagh_jj'], 'derived'), Gamma_rest_over_Hi=P(r0['Gamma_rest_over_Hi'], 'derived'),
        K_washout=P(r0['K_washout'], 'derived (mt/m*, m* = 1.08e-3 eV)'),
        one_loop_mnu=P(r0['one_loop'], 'estimate [E]', 'universal 22% absorbed in the fit; non-universal '
                       '<= 7e-4 eV (1.4% of m3) -> R re-tuned at that level'),
        alternative_S2=P(ss['scenarios']['S2_mt030']['mtilde_target_eV'], 'supplied alternative',
                         'tuning x2, rates x0.6 (needs an e2e_kinetic rerun)')),
    CP_phase=dict(
        delta_explicit=P(base['delta_needed'], 'fitted to eta_B (Planck 2018, 6.12e-10 [R])',
                         'window over washout x0.5..x2 and H_i 1.3e11..8.4e10: %.1e .. %.1e' % (
                             min(lp[k]['delta_needed'] for k in lp if k.startswith('S1') and 'no_non' not in k),
                             max(lp[k]['delta_needed'] for k in lp if k.startswith('S1') and 'no_non' not in k))),
        eps_per_unit_delta=P(base['eps_per_unit_delta'], 'computed (Covi-Roulet-Vissani, hep-ph/9605319)'),
        eps_at_fit=P(base['eps_at_delta_needed'], 'derived'),
        vacuum_phase_chi=P(ck['uv']['vacuum_phase_chi'], 'derived from the hybrid (v/f mod 2pi = 1.31 rad)',
                           'enters the relative phase of N1 vs N2-4: an O(1) contribution that the Yukawa phases '
                           'must cancel to ~1e-6 for the fit [A]'),
        relation_to_DM_sign=P('independent', 'derived [D/A]', 'the baryon sign is set by delta (sigma/chi/N sector); the '
                              'DM sign by alpha + sign(s_*) (holder); they are uncorrelated supplied inputs')),
    baryogenesis=dict(
        etaB_per_unit_delta=P(base['etaB_per_unit_delta'], 'computed [C on Q/E rates]'),
        etaB_thermal_only_per_unit_delta=P(base['etaB_thermal_only_per_unit_delta'], 'computed'),
        etaB_at_delta_0p1=P(base['etaB_at_delta_0p1'], 'computed', 'overproduced x%.0e at an O(0.1) phase' % (
            base['etaB_at_delta_0p1'] / base['eta_obs'])),
        dominant_source=P('non-thermal N4 from the late chi-remnant decays (chi -> N4 N4) at T ~ T_RH ~ M4/5',
                          'computed'),
        washout_survival=P(base['washout_survival_from_t'], 'computed', 'asymmetry present at t (H_i^-1) survives by this'),
        spontaneous_rotating_phases=P('erased: survival from t <= 30 H_i^-1 <= 1e-45; no tree-level theta W Wtilde '
                                      'coupling in the mass-phase basis', 'computed (survival) / argued (basis)'),
        sphaleron_conversion=P(28 / 79, 'quoted [R] (Harvey-Turner)'),
        flavour=P('unflavoured (T ~ 2-3e12 GeV, at the edge of tau equilibration)', 'estimate [E]')),
    Ns=dict(
        choice=P('massless dark radiation: no soft mass, no seesaw Yukawa; N_s -> i g5 N_s exact', 'supplied'),
        DeltaNeff=P(ck['Ns_dark_radiation'], 'computed on the e2e_kinetic history',
                    'Planck 2018 2.99 +- 0.17 (1807.06209); ACT DR6 2.86 +- 0.13 (2503.14454)'),
        soft_variant=P(ck['Ns_soft_variant'], 'estimate [E]', 'not adopted: adds a second CP phase and an s-odd piece')),
    uv=dict(
        completion=P('Phi_a = (rho_a/sqrt2) e^{i theta_a}, sigma = f1 theta_1, chi = f2 theta_2; M_j = y_j <rho_a>/sqrt2; '
                     'the hybrid potential is a monodromy (multi-branched) potential on the angles', 'supplied [A]'),
        radial_test=P(ck['uv'], 'computed', 'weakly coupled (lambda <= 4pi) completion is not frozen at the slice '
                      'amplitude: theta-dot/m_rho = 1.5-2.2; strong coupling (m_rho = 4 pi f) gives 0.6-0.9')),
    other=dict(
        lifetimes_s=P(ck['timing'], 'computed'),
        BBN=P('no new state lives past ~1e-29 s; N_s is the only new light state (dark radiation)', 'derived'),
        higgs_medium=P('N-gas y^2 per species >= mf_higgs-skeptic values (S1)', 'derived')))
json.dump(card, open(os.path.join(RES, 'sm_card.json'), 'w'), indent=1)
print('written', os.path.getsize(os.path.join(RES, 'sm_card.json')))
