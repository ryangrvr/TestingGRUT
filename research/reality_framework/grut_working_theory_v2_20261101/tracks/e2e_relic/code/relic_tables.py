#!/usr/bin/env python3
"""e2e_relic: writes ../results/relic_table.json (relic vs S/Q and Y_Q), ../results/relic_observables.json (benchmark
observables and regression) and ../results/relic_runs.json (the Y_Q runs of relic_runs.py, copied from scratch)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import relic as R
RES = os.path.join(HERE, '..', 'results')
SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_relic'

rows = []
for s in [1.84, 2.0, 2.12, 2.5, 2.88, 2.99, 3.0, 3.2, 3.28, 3.4, 3.4716, 3.5, 4.0, 4.16, 4.955, 6.0, 8.0, 10.0, 12.3, 14.7,
          17.1, 21.8, 30.2, 61.5, 155, 987]:
    for f in (None, 0.1, 1.0, 10.0):
        Y = None if f is None else f * (R.YQ_REF_QM if s <= R.SQ_BULK else R.YQ_M1)
        r = R.relic(s, Y_Q=Y)
        rows.append({k: r[k] for k in r if k != 'inputs'} | dict(S_per_Q=s, YQ_factor_input=f))
json.dump(dict(note='relic(S/Q, Y_Q): YQ_factor_input None = abundance-calibrated; otherwise x Y_Q,ref (6.0653e-9 for Q-matter, '
                    '5.0794e-9 for free X)', rows=rows), open(os.path.join(RES, 'relic_table.json'), 'w'), indent=1, default=float)

bench = {}
for lab, s, H in [('M1r_Qmatter_e3e-6_SQ2.99', 2.99, 1.49e11), ('HFa_c0_1.40_e1e-6_SQ2.88', 2.88, 8.40e10),
                  ('HFa_c0_1.36_e3e-6_SQ3.28', 3.28, 8.85e10), ('HFa_c0_1.40_e3e-6_SQ2.12', 2.12, 8.40e10),
                  ('M1r_freeX_e1e-6_SQ4.955', 4.955, 1.29e11), ('HFa_c0_1.36_e1e-6_SQ4.16', 4.16, 8.85e10),
                  ('HF_exported_c0_1_SQ14.7', 14.7, 1.11e11), ('HFa_Pauli_c0_1_SQ61.5', 61.5, 2.40e11)]:
    r = R.relic(s, H_i=H, e=1e-6)
    bench[lab] = dict(relic=r, observables=R.observables(r))
json.dump(dict(regression=R.regression(), benchmarks=bench), open(os.path.join(RES, 'relic_observables.json'), 'w'),
          indent=1, default=float)
runs = {}
for f in ('1', '0.1', '10'):
    try:
        runs[f] = json.load(open(os.path.join(SCR, 'yq%s' % f, 'relic_%s.json' % f)))
    except Exception as e:
        runs[f] = repr(e)
json.dump(runs, open(os.path.join(RES, 'relic_runs.json'), 'w'), indent=1, default=float)
for k, v in bench.items():
    r, o = v['relic'], v['observables']
    if r['branch'] == 'Q-matter':
        print('%-28s Q=%.3e [%.2e,%.2e] fv=%.4f/%.4f M=%.4f-%.4f m_obj=%.2f ug (%.2f M_Pl) R=%.3f nm flux=%.3f /m2/yr si=%.1e meff=%.0f-%.0f keV' % (
            k, r['Q_mean'], *r['Q_mean_range'], r['f_v']['perfect'], r['f_v']['beliaev'], *r['M_range'], o['object_mass_ug'],
            o['object_mass_over_M_Planck'], o['radius_nm'], o['flux_per_m2_per_yr'], o['self_interaction']['sigma_over_m_cm2_g'],
            min(o['vapour']['m_eff_Lya_keV']), max(o['vapour']['m_eff_Lya_keV'])))
    else:
        print('%-28s free X p=%.4e (aref %s) M=%.4f Xbar/X=%.3e n=%.2f/cm3 si=%.2e cm2/g' % (
            k, r['p'], '%.4e' % r['p_aref'] if r['p_aref'] else '-', r['M_central'], r['Xbar_over_X'], o['n_local_per_cm3'],
            o['self_interaction']['sigma_over_m_mix_cm2_g']))
