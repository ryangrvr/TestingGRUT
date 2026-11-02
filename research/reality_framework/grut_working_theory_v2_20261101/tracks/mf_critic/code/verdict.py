#!/usr/bin/env python3
"""mf_critic (Claude): assemble the benchmark verdict table -> results/verdict.json.

Sources: this track's results (settle_checks.json, m1r_c0.json, hyb_c0.json) [C/D/E here]; the reviewed tracks and skeptics
[Q]: mf_gemini (+ skeptic), mf_hybrid (+ skeptic), mf_ideal, mf_chain (+ skeptic), mf_higgs (+ skeptic).  Every entry is
conditional on that column's model choices (architecture, reheating closure, holder amplitude c0, rescaling rule,
renormalisation conditions, ledger, p map, survival criteria).  Data pulls are retrospective.
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / 'results'
sc = json.load(open(RES / 'settle_checks.json'))
m1 = json.load(open(RES / 'm1r_c0.json'))
hy = json.load(open(RES / 'hyb_c0.json'))


def pull(ns):
    return dict(planck=(ns - 0.9649) / 0.0042, act=(ns - 0.9743) / 0.0034)


def chain(s):
    a, b = s['e1e-6'], s['e3e-6']
    return dict(SQ_1e6=a['S_per_Q'], SQ_3e6=b['S_per_Q'], SQ_inf=s['S_per_Q_e_inf'], branch_1e6=a['branch'], branch_3e6=b['branch'],
                e_crit_bulk=s['e_crit_bulk'], e_crit_Q2e5=s['e_crit_Q2e5'], M_1e6=a['M'], H_req_GeV=a['H_i_req_GeV'], K=s['K'],
                eps=s['eps'], N_YL=s['N_YL'], N_YH=s['N_YH'], R_vec=a['R_vec'], R_Y=a['R_Y'], R_casc=a['R_casc'])


T = sc['K7']
cols = {
 'proposal_as_written': dict(
    label='[Q mf_gemini + skeptic; D here]',
    architecture='V = lam/4 (chi^2 - v^2)^2 + g^2 sigma^2 chi^2/2 + V_SR; m_chi = 3 H0, lam v^4/4 = V0 -> v = 1.633 M_P, lam = 4.81e-15 (H0 = 1.3e11)',
    n_s='convex V_SR: 1.004-1.04 (blue); hilltop |eta| = 0.0136-0.0166; log family with the waterfall attached >= 0.972 (in 0.967-0.973 only for Delta >~ 0.2-0.3)',
    r=2.75e-7, alpha_s='hilltop ~ -4e-9; log ~ -8e-4', H_i='1.25-1.55e11 GeV asserted; not self-consistent (no history closes)',
    waterfall='mild: >= 17 e-folds (16.9-41), O(1) delta-N tail, PBH overclosure f_PBH >= 1e10; real chi walls at 0.99 of the topological-inflation threshold',
    reheating='y_N chi NN (blocked x4e7; CW x2e13; loop cap y <= 3-5.5e-4 -> Gamma <= 2e-8 H) + Higgs portal (tachyonic Higgs or Landau pole 9e12 GeV); best preheating K >= 2.2 or never',
    K_needed='<= 1.38-1.41 for Q-matter at 1e-6 on M1r\'s own release', chain='not run: no consistent history',
    higgs='portal as written: m_H^2 = -9e7 H^2 during inflation -> fails', act=pull(0.972),
    verdict='fails its own requirements'),
 'HFa_exported': dict(
    label='[Q mf_hybrid, mf_chain, skeptics]',
    architecture='heavy small-field hybrid, tilted hilltop clock + chi loop, biased real waterfall; m_chi = 250, m_sig = 600 H_i; lam 2.32e-7, v 0.0196 M_P, g 1.63e-3, sigma_c 0.0058 M_P, kappa = 10; spectator m_s = g_s chi; axial N (4 seesaw singlets)',
    n_s=0.9649, r=2.75e-7, alpha_s=1.07e-4, H_i='1.30e11 design; self-consistent 1.28e11 (track map) / ~1.11e11 (gap-filled map)',
    waterfall='fast (eps_H = 1 at 0.19 H_i^-1); no PBH; walls excluded only with unblocked friction',
    reheating='axial N -> LH, UNBLOCKED (Pauli blocking ignored: refuted)', K=2.048, T=T['HFa exported (unblocked, c0-independent background)'],
    spectator=dict(s_star_GeV=7.3e15, beta_iso='6.9e-3 (track map) / 0.017-0.022 (gap-filled)', f_after_turn_on=1.35e-5, f_at_decay=0.013,
                   decay='y_s N_s, Gamma_s = 1.3e-6 H_i; needs soft mass + seesaw or N_s is dark radiation'),
    chain=dict(SQ_1e6=14.69, SQ_3e6=12.25, e_crit=None, branch='free X', M='1.58 (track map) / 1.77 (gap-filled)'),
    higgs=dict(Mt_crit_grav=172.88, with_media='>= 174.6'), act=pull(0.9649)),
 'HFb_exported': dict(
    label='[Q]', architecture='as HFa, clock retuned (c1 1.75e4, m2 0.0098); spectator m_s^2 = 9 - (3/4) R (xi_s = -3/4)',
    n_s=0.9743, r=2.75e-7, alpha_s=1.43e-4, H_i='1.28e11 self-consistent (track map)', K=2.047,
    T=T['HFb exported'], chain=dict(SQ_1e6=14.58, SQ_3e6=12.15, e_crit=None, branch='free X', M=1.57),
    higgs=dict(Mt_crit_grav=172.88, with_media='>= 174.6'), act=pull(0.9743)),
 'HFa_Pauli_P0_selfconsistent_c0_1': dict(
    label='[Q mf_chain; skeptic]', reheating='axial N, Pauli-blocked fluid (Gamma_th = 0, tails free)', K=7.84, H_i='2.40e11 GeV',
    r=9.37e-7, n_s='0.9649 retuned; CMB capacity 48.6 < N_* at kappa = 10 (needs kappa = 15-50); HFb passes (53.7)',
    T=T['HFa Pauli P0 sc 2.40e11 (mf_chain)'],
    chain=dict(SQ_1e6=61.5, SQ_3e6=30.2, e_crit=None, branch='free X', M=6.82),
    spectator=dict(beta_iso='<= 4.5e-3', f_at_decay='0.1-0.2%'), higgs=dict(Mt_crit_grav=172.74, with_media='>= 174.6')),
 'HFa_Pauli_P0_strong_holder': dict(
    label='[C here: c0 = 1.40 reproduced to all digits on a rebuilt history; Q skeptic for c0 = 1.36 and robustness]',
    holder='c0 = 1.36 / 1.40 at the 10i phase (|g_h| s_* = 16.7 / 17.2 H_i)', K=6.50, H_i='8.4-8.9e10 GeV (c0 1.40: H_req/H_i = 0.998; c0 1.36: 1.05)',
    r=sc['K5']['rows']['strong holder P0']['r'], axial_EFT_ratio=sc['K5']['rows']['strong holder P0']['axial_EFT_ratio'],
    T=T['HFa Pauli P0 8.40e10 (rebuilt here)'],
    chain={'1.40': chain(hy['crP0H084']['runs']['1.400'])},
    chain_c0_1p36='S/Q 4.16 / 3.28, e_crit 2.18e-6 (bulk) / 7.7e-6 (Q = 2e5), free X M = 1.000 [Q skeptic]',
    relic_1p40='Q-matter, M = 0.807 m_X; <Q> ~ 1e20 and object ~0.6 M_Pl transferred from qm_micro (+ skeptic) [Q, not recomputed]',
    spectator=dict(D=9.2, s_star_min_GeV=1.36e16, leak='c_L <~ 1/16 pi^2 or vanish at M_Y'),
    higgs=dict(Mt_crit_grav_E=sc['K6']['Mt_crit']['8.4e+10'], with_media='>= 174.6'), cmb='capacity 99.8 >= N_* 51.9 at kappa = 10'),
 'HFa_Pauli_Gth10_strong_holder': dict(
    label='[C here]', K={'8.4e10': hy['crP10H084']['meta']['reheating']['K'], '5.0e10': hy['crP10H050']['meta']['reheating']['K']},
    T={'8.4e10': T['HFa Pauli Gamma_th=10 8.40e10 (here)'], '5.0e10': T['HFa Pauli Gamma_th=10 5.0e10 (here)']},
    runs={tag: {k: chain(v) for k, v in blk['runs'].items()} for tag, blk in hy.items() if tag.startswith('crP10')},
    axial_EFT_ratio_at_5p2e10=sc['K5']['rows']['strong holder Gamma_th = 10']['axial_EFT_ratio'],
    higgs=dict(Mt_crit_grav_E=sc['K6']['Mt_crit']['5.2e+10'], with_media='>= 174.6'),
    reading='self-consistent M1r-like point at c0 ~ 1.25, H_i ~ 5.2e10 (S/Q 3.64, e_crit 1.3e-6); Q-matter at 1e-6 needs c0 ~ 1.27-1.30 at H_i ~ 4.3-5e10; axial EFT ratio 1.05-1.5 there'),
 'M1r_supplied_reference': dict(
    label='[C here: M1r (3, 0.5) static slice, mf_ideal bridge table]', K=1.98581, T=T['M1r (3, 0.5) supplied'],
    runs={k: chain(v) for k, v in m1.items()}, c0_threshold_at_1e6=sc['K8']['c0_threshold_at_1e6'],
    dln_ecrit_dln_c0=sc['K8']['dln_ecrit_dln_c0_near_1']),
}
json.dump(dict(note=__doc__, columns=cols), open(RES / 'verdict.json', 'w'), indent=1, default=float)
print('written', RES / 'verdict.json')
for k, v in cols['M1r_supplied_reference']['runs'].items():
    print('M1r c0', k, '%.3f %.3f' % (v['SQ_1e6'], v['SQ_3e6']), v['e_crit_bulk'], v['branch_1e6'], '%.3g' % v['H_req_GeV'])
