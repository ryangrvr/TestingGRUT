#!/usr/bin/env python3
"""inf_critic (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): assemble the per-benchmark verdict
table of the inflaton-sector pass from the reviewed tracks (READ-ONLY) and this track's checks.

Sources: tracks/inf_chain/results/summary.json [C there]; tracks/inf_chain/skeptic/results/sk_higgs_TB.json [C there];
inf_higgs NOTES/SKEPTIC numbers [Q, typed in below with their source]; ../results/settle_checks.json and
../results/quiet_scan.json (this track).  Every entry is conditional on the model choices of its track.
Output ../results/verdict.json
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]
RES = HERE.parent / 'results'
S = json.load(open(TRK / 'inf_chain' / 'results' / 'summary.json'))['rows']
K = json.load(open(RES / 'settle_checks.json'))
Q = json.load(open(RES / 'quiet_scan.json'))

POT = {
    'HT6t': 'tilted p=6 hilltop V0(1-beta x-x^6)^2/f_top^2, v=0.34694 M_P, beta=1.5529e-5 (no phi_I -> -phi_I symmetry)',
    'HT6t_cal1': 'same family, v=0.122614 M_P, beta=1.23528e-6 (re-tuned at 2.805e10 GeV)',
    'HT6t_Send': 'HT6t, holder calibrated at eps_H=1 (bracket)',
    'T1': 'T-model V0 tanh^2(phi/sqrt(6 alpha)), alpha=6.0968e-5',
    'T1_cal1': 'T-model alpha=3.7855e-7 (re-tuned at 1.043e10 GeV)',
    'TB1': 'T-model alpha=0.3',
    'TB2': 'T-model alpha=1',
}
CH = {
    'HT6t': 'Gamma_I=0.5 H_i (engine convention); kinematically blocked for y phi_I NbarN (needs y~0.13-0.54, y|phi|/(m/2)~1.8e4): no realisable channel',
    'T1': 'Gamma_I=0.5 H_i; blocked (y|phi|/(m/2)~3.1e3); no realisable channel; fragments at ~0.35-0.43 H_i^-1',
    'TB1': 'Majorana -(y/2)phi_I NbarN, y=4.1e-4 (S1, inside y_nat 4-9e-4) or axial S2 with f~2-9e15 GeV; Gamma_I=7.5e-9 H_i fixed by Omega_DM; unblocked ~3 e-folds after S1, decay after ~12',
    'TB2': 'as TB1, y=2.6e-4 or S2 f~4e15-1.4e16 GeV; Gamma_I=2.1e-9 H_i',
}
CH['HT6t_cal1'] = CH['HT6t']; CH['HT6t_Send'] = CH['HT6t']; CH['T1_cal1'] = CH['T1']
HIGGS = {
    'HT6t': 'at 1.29e11 [Q inf_higgs + skeptic]: inflation safe (xi_min 0.070 tree at 172.57; fails >=173.75 tree/174.15 loop); '
            'post-inflation gravitational-only M_t,crit 172.92 (MD)/172.875 (Gamma=0.5) = +1.2/+1.05 sigma; at the DM-consistent 2.7e10 '
            'extrapolated ~173.4 [E]; S1 plasma reach capped by inflation, needs |A|<~1.9 TeV; lam_IH>=3.4e-11 safe',
    'T1': 'as HT6t (T1/E1 bracket 172.87-172.92 for the fragmented fluid); lam_IH>=3.1e-10 (T1), 1.2e-9 (E1) safe',
    'TB1': 'inflation, no portal: M_t,crit ~173.10-173.13 GeV (+1.8 sigma) [C, K7; inf_chain skeptic xi_min 0.0865/0.0938/0.1006]; '
           'post-inflation gravitational-only extrapolated ~172.1 GeV (fails at central M_t) [E, K8]; with lam_IH>=lam_min=3.1e-10 '
           '(>=~2.5 lam_min if S1 loop renormalised at M_N): log10<N> <= -3.4e4 to 174.6 GeV [C, K7], post-inflation transfer of inf_higgs portal result [E]',
    'TB2': 'inflation, no portal: xi_min 0.0871/0.0944/0.1011 -> M_t,crit ~173.1 [Q inf_chain skeptic]; post-inflation ~172.1 [E]; lam_IH>=3.0e-10 safe [E]',
}
HIGGS['HT6t_cal1'] = 'not computed (H_i 2.8e10); extrapolated grav-only M_t,crit ~173.3 [E]'
HIGGS['HT6t_Send'] = HIGGS['HT6t']
HIGGS['T1_cal1'] = 'not computed (H_i 1.0e10); extrapolated grav-only M_t,crit ~173.6 [E]'
LEPTO = {
    'TB1': 'non-thermal from phi_I -> NN: eps_1 = 1.35e-8 needed (Y_L = 1.5 eps T_RH/m_I); M_N = 2-7e12 GeV gives 0.9-4.5e-5 of the '
           'Davidson-Ibarra bound with negligible inverse-decay washout (z_RH >= 10); M_N = 1e12 is washed out (exponent 9.6) [E]',
    'TB2': 'eps_1 = 2.1e-8 needed; M_N >= 1e12 works (washout exponent 1.3 at 1e12) [E]',
}
OPEN = {
    'HT6t': 'not DM-consistent at its design H_i (ratio 0.21; secant fixed point ~2e9 GeV, r~7e-11, not run); no reheating channel; no R symmetry (CP bookkeeping fails); ACT n_s -2.8 sigma',
    'T1': 'not DM-consistent (0.076); near fixed point ~1.2e10 GeV (T1_cal1) is fragmenting and needs Gamma=0.5 (blocked); ACT -3.7 sigma',
    'TB1': 'p(S/Q~600) recipe (sets M, Gamma_cal ~ M^-2, T_RH ~ M^-1); |A_I Sigma| <~ 1e-6 H_i (else leak R~6e2-9e5 re-calibrates Gamma); lam_IH added (an R-even SM portal: model choice); '
           'Higgs post-inflation at H~7e12 and plasma timing not computed; thin x2.4 margin of "no e_crit at any e"; ACT n_s -3.7 sigma',
    'TB2': 'as TB1 (margin x5.5); r = 4.2e-3',
}
OPEN['HT6t_cal1'] = OPEN['HT6t']; OPEN['HT6t_Send'] = 'charge sign flips vs S1 calibration (CP sign depends on holder calibration point); H_i,req 1.3e9 (ratio 0.010)'
OPEN['T1_cal1'] = OPEN['T1']

rows = {}
for k, r in S.items():
    d = dict(potential=POT[k], reheating=CH[k], n_s=r['n_s'], r=r['r'], H_star_GeV=r['H_star_GeV'], H_i_GeV=r['H_i_potential_GeV'],
             m_I_over_Hi=r['m_I_over_Hi'], Gamma_I_over_Hi=r['Gamma_I_over_Hi'], N_star=r['N_star'],
             T_RH_GeV=K['K4_reheating_temperatures'][k]['T_RH_GeV'],
             T_RH_note='realisable' if k.startswith('TB') else 'nominal only (Gamma=0.5 kinematically blocked)',
             K=(r.get('reheating', {}).get('K') or r['release']['K']), K_norm='a=1 at S1 (M1r: 1.986)',
             S_per_Q_e1e6=r['S_per_Q_e1e6'], S_per_Q_e_inf=r['S_per_Q_e_inf'], branch_e1e6=r['branch_e1e6'],
             e_crit=r['e_crit_bulk'], relic='free X, mostly symmetric', p=r['p'], M_per_charge_mX=r['M'],
             asymmetric_share=1 / r['M'], H_i_required_GeV=r['H_i_req_GeV'], consistency_ratio=r['consistency_ratio'],
             higgs=HIGGS[k], leptogenesis=LEPTO.get(k, 'not defined (no realisable reheating channel)'), open=OPEN[k],
             ns_pull_ACT_DR6=K['K9_ns_pulls'][k]['pull_ACT_DR6'])
    rows[k] = d
rows['E1'] = dict(potential='E-model V0(1-exp(-phi/sqrt(3 alpha/2)))^2, alpha=6.0969e-5 (same inflation as T1)', m_I_over_Hi=256,
                  n_s=0.9618, r=2.71e-7, H_i_GeV=1.29e11, reheating=CH['T1'],
                  S_per_Q_e1e6=285, note='[Q inf_potential skeptic: own MODE-4 release N_YH 9.2 (x883 M1r), dq_heat ledger; not run through the chain]',
                  branch_e1e6='free X', e_crit=None)
out = dict(method=__doc__, rows=rows, iteration=K['K2_iteration_secant'], quiet_threshold_T=Q['quiet_threshold_T'],
           fast_closure=Q['fast_closure'], quiet_closure=Q['quiet_closure'])
json.dump(out, open(RES / 'verdict.json', 'w'), indent=1)
for k, d in rows.items():
    print(k, {x: d.get(x) for x in ('n_s', 'r', 'H_i_GeV', 'T_RH_GeV', 'K', 'S_per_Q_e1e6', 'M_per_charge_mX', 'consistency_ratio')})
