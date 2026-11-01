#!/usr/bin/env python3
"""m1_ledger_skeptic: assemble the per-claim verdict table (skeptic_summary.json) from skeptic_algebra.json and
skeptic_born.json.  Verdict semantics: refuted=true only with a concrete recomputed reason."""
import json
from pathlib import Path
R = Path(__file__).resolve().parents[1] / 'results'
A = json.loads((R / 'skeptic_algebra.json').read_text())['results']
B = json.loads((R / 'skeptic_born.json').read_text())['results']
f = A['frozen']; o = A['old_target']; lc = A['lambda_chain']; ck = A['clock_ambiguity']; ex = A['exponent_from_codex_rows']
g = A['gamma_t0']; cc = A['cold_completion']; pr = A['pulse_rows']['rows']; bm = B['M1']['native_quad']; bo = B['old']['native_quad']
checks = [
 dict(claim_id='1', refuted=False, recomputed='c_eta %.10f, c_p %.9f, eta_min %.10f, p_max %.13f (mpmath, from v g M Q D)' % (f['c_eta'], f['c_p'], f['eta_min'], f['p_max']),
      reason='exact identity; agrees with Codex CONSEQUENCES.json formation_failure_test to 2e-15 and with m1_ledger to all quoted digits'),
 dict(claim_id='2', refuted=False, recomputed='Y_Q = %.15e (rel to Codex %.1e); Y_Q eps_Q = D exactly' % (f['YQ_new'], f['YQ_new_rel_codex']), reason='D/eps_Q with eps_Q = 17324.832346421026/200000 GeV'),
 dict(claim_id='3', refuted=False, recomputed='D/Y_Q,old = %.10f GeV = %.6f m_X = d71cc01 late mass per charge (identity %.0e); = object 0.923316 m_X + 2 x %.6f m_X released pairs; old fraction %.7f, multiplier %.7f, kin f coefficient %.6f' % (o['D_over_YQ_old_GeV'], o['over_mX'], o['identity_YQ_old_times_late_mass_over_D'], o['released_pairs_per_charge'], o['old_object_fraction'], o['old_multiplier'], o['kin_f_coefficient']),
      reason='confirmed and sharpened: the v-calibrated 0.98593 m_X is object + 0.031306 scalar-release pairs per charge (d71cc01 source_scalar_number 6261.2; ABUNDANCE_LEDGER required_pair_number 0.031306)'),
 dict(claim_id='4', refuted=False, recomputed='lambda = %.13f from Codex b3 baseline row (Codex %.13f); Codex rows obey Y = |N_X| H^3/(a^3 s), g* = 106.75, to 4e-16 (no engine); clock band 1.00484-1.00834' % (lc['lambda_M1_b3'], lc['codex']),
      reason='the chain is arithmetic on the baseline yield 1.040807 Y_Q,old; the 8-digit lambda is uncertain at -0.22/+0.13% absolute from the clock choice alone and -0.40% if D is taken at the Planck 2018 value'),
 dict(claim_id='5', refuted=False, recomputed='1.0172235041614 (b3 old object), 1.0162835956038 (b5 old object; m1_ledger wrote 1.016277, a 7e-6 slip), 1.006078579112 (b5 M1)' , reason='Codex rows b3_g0.5_matched 1.0172235042, b5_g1_matched 1.0162835956; superseded old-object targets as m1_ledger says'),
 dict(claim_id='6', refuted=False, recomputed='Codex scale_0.5/scale_2.0 control rows (factor 4 in H_i, not used by m1_ledger): exponent_Y %.10f, exponent_T %.10f; N_X and a unchanged to %.1e / %.1e' % (ex['Y_scale2_over_scale05']['exponent'], ex['T_scale2_over_scale05']['exponent'], abs(ex['NX_spread_over_factor4']['s2']), abs(ex['a_spread_over_factor4']['s2'])),
      reason='confirmed by a different route (Codex\'s own factor-4 controls + the row-level yield identity) rather than our engine; the "ODE rerun" is a unit conversion of an unchanged dimensionless state'),
 dict(claim_id='7', refuted=False, recomputed='Y recomputed from Codex row (N_X, a, T, H): 5.026475516056737e-9 and 5.03345022205823e-9 (4e-16)', reason='not re-run here (engine); the m1_ledger/Codex 1e-11 engine agreement is a same-formula (g* = 106.75, tau >= 0.8 end window) two-engine comparison, which is what it claims'),
 dict(claim_id='8', refuted=False, recomputed='rp_reheating extra.json yr_b3_g0.5: end 80, late N_X -8156.7346 (window 64-80), late90 -8153.30; 8156.734642/8120.740082 = %.7f' % ck['m1_ledger_window_ratio_tau80_over_tau800'],
      reason='the tau = 80 origin is verified; note m1_ledger pairs this window ratio with Codex\'s instantaneous clocks (-8147.1/-8116.1 = %.5f), a different estimator of the same drift' % ck['codex_instantaneous_tau80_over_tau800']),
 dict(claim_id='9', refuted=False, recomputed='pulse rows with charge-conserving eta = 1 + (Nbar-Nplus)/Q0: terms %.5f/%.5f/%.5f/%.4f (m1_ledger %.5f/%.5f/%.5f/%.4f, 3%% high); kin rows and crossing R 5.824 -> 5.768 reproduced; mask rows are not eta by Codex\'s own statement' % (pr['5%']['term_charge_conserving'], pr['10%']['term_charge_conserving'], pr['20%']['term_charge_conserving'], pr['40%']['term_charge_conserving'], pr['5%']['term_m1_ledger_eta1'], pr['10%']['term_m1_ledger_eta1'], pr['20%']['term_m1_ledger_eta1'], pr['40%']['term_m1_ledger_eta1']),
      reason='verdict (no formation-type population passes unconditionally) stands; the pulse rows should use the coefficient c_p - c_eta = 2.0668 on Nbar/Q0 (core charge rises), and the frozen formula also ignores the pulsed core\'s 2.3% excitation energy (c481530)'),
 dict(claim_id='10', refuted=False, recomputed='excess from (E(Q_f) + N_free m_X)/M_ref - 1: %.6e/%.6e/%.6e/%.6e (rel to Codex <= 2e-13); frozen/actual 0.529-0.551' % tuple(cc[k]['excess_from_masses'] for k in ('5%', '10%', '20%', '40%')), reason='Codex COLD_COMPLETION_RESULTS reproduced; Q_f + N_free = Q0 to 1e-11'),
 dict(claim_id='11', refuted=False, recomputed='A from M1 parameters with alpha6 = kappa^2 mu <Sigma>/(M_Y^4 - 4 mu^2 <Sigma>^2), <Sigma> = h/sqrt2: %.15e (rel %.1e); A_old,full/lambda^2 = Codex A_new to %.1e; r^2 = 0.2025' % (A['A_coupling']['A_from_M1_parameters'], A['A_coupling']['rel'], A['A_coupling']['rel_to_codex_new']),
      reason='own evaluation from the heavy parameters, not only the lambda^-2 scaling'),
 dict(claim_id='12', refuted=False, recomputed='I_M1 = %.7f (native grid dr 0.05, Simpson transform, adaptive quad with (qmax-q)^{1/2} weight; rel to Codex %.1e), spline/trapezoid variant rel %.1e; I_old = %.6f (%.1e); Parseval closure %.0e; Born c1 rates %.6e / %.6e' % (bm['I'], bm['rel_I_codex'], B['M1']['spline_trapezoid']['rel_I_codex'], bo['I'], bo['rel_I_m1'], abs(bm['parseval_rel']), bm['born_fractional_core_rate_c1'], bo['born_fractional_core_rate_c1']),
      reason='prefactor 900/(32 pi^3) re-derived (two-body factor (pi/2) beta, identical-particle 1/2, chi^2 contraction 2, C(6,2) = 15); charge loss 4 per event'),
 dict(claim_id='13', refuted=False, recomputed='Gamma t0: new %.7e (rel %.1e), old %.7e = e160356 7.08367857848652e-17 (m1_ledger\'s 7.0836800e-17 used the 6-digit A, 2e-7 off); ratio %.6f = profile %.7f x coupling %.7f' % (B['Gamma_t0']['M1'], B['Gamma_t0']['rel_M1'], B['Gamma_t0']['old'], g['ratio_new_over_old'], g['profile_factor'], g['coupling_factor']),
      reason='coherent channel quoted from Codex (6.4e-5 of the total); c6 = A M_red^2; age 4.3549488e17 s'),
]
new_issues = [
 'N1 [E,R]: the target D = 4.4e-10 GeV is a rounded value; Omega_c h^2 = 0.1200 +- 0.0012 (Planck 2018), rho_c/h^2 = 1.05368e-5 GeV cm^-3 and s0 = 2891.2 cm^-3 give D = %.4e GeV (Codex 0.61%% high, +-1%% band); at the Planck value lambda shifts by %.2f%% and the band is +-%.2f%%: the +0.7009%% calibration step is inside the target\'s own rounding and uncertainty. Not raised by Codex or m1_ledger.' % (A['target_D']['D_planck2018'], 100 * A['target_D']['lambda_shift_if_D_planck'], 100 * A['target_D']['lambda_band_from_Och2']),
 'N2 [C]: lambda = 1.00700905 is quoted to 1e-8 but Codex\'s own charge_clock entries (tau 20/50/80/200/800 instantaneous N_X) give lambda = %.5f-%.5f (window average %.5f): -0.22/+0.13%% absolute, i.e. a third of the calibration step; m1_ledger gave the 0.44%% Y drift but not its lambda band.' % (ck['lambda_min'], ck['lambda_max'], ck['lambda_window_640_800']),
 'N3 [D]: m1_ledger\'s pulse rows read the emitted Xbar as p at eta = 1; charge conservation makes eta = 1 + (Nbar - Nplus)/Q0 > 1 (ed73791), so the frozen term is 2.0668 Nbar/Q0 + 0.0668 Nplus/Q0, 3% below the table (10%: 0.01749 vs 0.01802; x1.75 not x1.8).',
 'N4 [D]: m1_ledger\'s (5,1) old-object repair "1.016277" is a transcription slip for 1.0162836 (Codex row b5_g1_matched uniform_scale 1.0162835956); superseded number, +1.628% unchanged.',
 'N5 [C]: Codex\'s scale_0.5 and scale_2.0 control rows (not used by m1_ledger) already fix the exponent over a factor 4: 1.50000000026 (Y), 0.4999999996 (T_vis), with N_X and a unchanged to 5e-10; m1_ledger\'s own-engine measurement adds nothing beyond this and the "ODE rerun" cannot fail.',
 'N6 [A]: the frozen formula applied to pulsed states uses eps_Q for an excited core; c481530 gives 2.305% excess energy in the full 10% pulse and 0.2237% unprocessed remainder, so the eta-row readings of pulses are lower bounds on the mass excess, not only the kinetics row.',
 'N7 [A]: T_vis at tau = 800 (7.7e12 GeV) exceeds M_Y = 4.03e12 GeV; g* = 106.75 is the SM count and holds only because the conservative candidate has zero SM portals; the "separately tested portal extension" would need g* and the yield recomputed (m1_ledger lists g* as an assumption; the reason is this).',
]
out = dict(track='m1_ledger_skeptic', date='2026-11-01', scope='adversarial recomputation of tracks/m1_ledger; algebra and quadrature only, no engine runs, no Codex code run; packets read via git show',
           independence=dict(own_code=['skeptic/code/skeptic_algebra.py', 'skeptic/code/skeptic_born.py'], reproduction=[], codex_code_run=False, engine_runs=0),
           checks=checks, new_issues=new_issues, refuted_count=sum(c['refuted'] for c in checks))
(R / 'skeptic_summary.json').write_text(json.dumps(out, indent=1))
for c in checks: print(c['claim_id'], 'REFUTED' if c['refuted'] else 'stands', '|', c['recomputed'][:150])
for n in new_issues: print(n[:160])
