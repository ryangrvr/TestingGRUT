#!/usr/bin/env python3
"""inf_potential (Claude, own code): assemble the benchmark table (results/summary.json) from bench_*.json,
exit_factor.json, yy_harmonics.json, modes_*.json, engine_test.json and scan_H.json, and the rough [E] re-ledger of the
hidden heat on the engine smoke-test release (see NOTES section 6; NOT a release analysis)."""
import json, math
from pathlib import Path
RES = Path(__file__).resolve().parent.parent / 'results'
C_HART = 6.977568          # E_c = (N/c)^{4/3} [Q wt_critic chain_consistency / dq_critic]
EX_M1R, NX_M1R, K_M1R = 188706.29, 8045.485140, 1.98581
DEP_VEC_M1R = 1.7805 * EX_M1R      # per vector at e = 1e-6 (K = 1.986) [Q dq_heat]
NA_M1R = 2.026                      # vectors per H_i^3 at e = 1e-6 [Q dq_heat]
GL = 1.006e-7; ML = 22.249


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def main():
    S = dict(note='every number is a model choice; labels in NOTES.md', benchmarks={})
    ef = json.load(open(RES / 'exit_factor.json'))['results']
    yy = json.load(open(RES / 'yy_harmonics.json'))['runs']
    eng = json.load(open(RES / 'engine_test.json'))
    for name in ('T1', 'E1', 'HT6t', 'HT6', 'HT8', 'HT4', 'PL2'):
        d = json.load(open(RES / ('bench_%s.json' % name)))
        rows = []
        for r in d['by_Hi']:
            h = {round(t['N_before_end'], 3): t for t in r['history']}
            row = dict(H_i_GeV=r['H_i_GeV'], params=r['params'], V0_quarter_GeV=r['V0_quarter_GeV'], Nstar=r['Nstar'], As=r['As'],
                       ns=r['ns'], ns_pull=r['ns_pull_sigma'], r=r['r'], alpha_s=r['alpha_s'], n_t=r['n_t'],
                       Hstar_GeV=r['Hstar_GeV'], Hend_over_Hi=r['Hend_over_Hi'], m_I_GeV=r['m_I_GeV'],
                       m_I_over_Hi=r['m_I_over_Hi'], m_I_over_Hend=r['m_I_over_Hend'], phi_min_MP=r['phi_min_MP'],
                       phiI_S1_MP=r['phiI_S1_MP'], phiI_end_MP=r['phiI_end_MP'], phidot_S1=r['phidot_S1_over_HiMP'],
                       g_I_S1=r['g_I_S1'], g_I_end=r['g_I_end'], alpha_I_S1=r['alpha_I_S1'], K_end_gamma05=r['K_end_gamma05'],
                       s_S1_at=dict((str(k), h[k]['s_S1']) for k in (2.0, 3.8, 4.8, 10.0) if k in h),
                       s_end_at=dict((str(k), h[k]['s_end']) for k in (1.0, 3.8, 10.0) if k in h),
                       Nstar_alt={'%s' % (a['history'],): (a.get('Nstar'), a.get('ns')) for a in r['Nstar_alternatives']})
            if 'osc_gamma0' in r:
                o = r['osc_gamma0']
                row['osc'] = dict(first_turn_phiI_MP=o['turning'][0]['phiI_turn'], second_turn_phiI_MP=o['turning'][1]['phiI_turn'],
                                  V_first_turn_over_V0=o['turning'][0]['V_over_V0'], omega_first_over_mI=o['omega_first_over_mI'],
                                  omega_min_over_Hi=o['omega_min_over_Hi'], omega_max_over_Hi=o['omega_max_over_Hi'], t_max_Hi=o['t_max_Hi'])
            rows.append(row)
        ent = dict(by_Hi=rows)
        for conv in ('S1', 'Send', 'Cexit'):
            k = '%s_%s' % (name, conv)
            if k in ef:
                ent.setdefault('exit_factor_Hi1.29e11', {})[conv] = dict(f_eff=ef[k]['f_eff_e1e6'], SQ_e1e6=ef[k]['SQ_e1e6'],
                                                                          branch=ef[k]['branch_e1e6'], e_crit_bulk=ef[k]['e_crit_bulk'])
        for run in yy:
            if run['benchmark'] == name:
                ent.setdefault('yy_harmonics_Hi1.29e11', {})['gamma%.1f' % run['gamma_over_Hi']] = {
                    lab: dict(g_I=v['g_I'], N_YL=v['N_YL'], N_YH=v['N_YH'], R_YY_MR1units=v['R_YY_with_tail']) for lab, v in run['res'].items()}
        mf = RES / ('modes_%s.json' % name)
        if mf.exists():
            md = json.load(open(mf))
            if 'floquet' in md:
                ent['floquet_first'] = md['floquet'][0]
            if 'resonance' in md:
                ent['fragmentation'] = {('gamma%.1f' % x['gamma_over_Hi']): x['first_ratio_0p1'] for x in md['resonance']}
            if 'ycontrol' in md:
                ent['ycontrol'] = md['ycontrol'][0]['N_YL_LL_linear']
        if name in eng.get('smoke', {}):
            sm = eng['smoke'][name]
            ent['engine_smoke'] = sm
            # rough re-ledger [E]: absolute deposits on the smoke-test release (end normalisation, a = 1 at eps_H = 1)
            K = rows[0]['K_end_gamma05']
            fe = ef.get('%s_S1' % name, {}).get('f_eff_e1e6', 1.0)
            dep_YL = ML * math.sqrt(math.pi * K / (2 * GL))
            vec = NA_M1R * (K / K_M1R) ** 1.5 * fe * DEP_VEC_M1R * math.sqrt(K / K_M1R)
            hom = sm['NYL'] / math.e ** 3 * dep_YL
            yyr = [run for run in yy if run['benchmark'] == name and run['gamma_over_Hi'] == 0.5]
            yyd = yyr[0]['res']['S1']['R_YY_with_tail'] * EX_M1R * (dep_YL / (ML * 5569.0)) if yyr else 0.0
            NX = abs(sm['NX_endnorm']); EXv = sm['EX_endnorm']
            Ec = (NX / C_HART) ** (4 / 3)
            ent['rough_reledger_E'] = dict(K=K, E_X=EXv, N_X=NX, E_c=Ec, D_vectors_e1e6=vec, D_homYL=hom, D_YY=yyd,
                                           SQ_vectors_only=SQ((EXv + vec) / Ec), SQ_plus_homY=SQ((EXv + vec + hom) / Ec),
                                           SQ_all=SQ((EXv + vec + hom + yyd) / Ec), SQ_all_e_inf=SQ((EXv + hom + yyd) / Ec))
        S['benchmarks'][name] = ent
    S['scan_H'] = json.load(open(RES / 'scan_H.json'))['rows']
    S['engine_regression_mode0_identical'] = eng['mode0_identical']
    S['engine_mode4_vs_python'] = {k: (v.get('max_rel_dH'), v.get('max_rel_dphihat')) for k, v in eng['mode4'].items()}
    json.dump(S, open(RES / 'summary.json', 'w'), indent=1, default=float)
    for name, ent in S['benchmarks'].items():
        r = ent['by_Hi'][0]
        print('%-5s %s V0^1/4=%.3e N*=%.2f ns=%.4f(%+.1fs) r=%.2e as=%.1e H*/Hi=%.4f Hend/Hi=%.3f mI=%.2e GeV mI/Hi=%.1f phiI(S1)=%.4f phiI(end)=%.4f gI=%.2e/%.2e K=%s' % (
            name, {k: '%.4g' % v for k, v in r['params'].items()}, r['V0_quarter_GeV'], r['Nstar'], r['ns'], r['ns_pull'], r['r'], r['alpha_s'],
            r['Hstar_GeV'] / r['H_i_GeV'], r['Hend_over_Hi'], r['m_I_GeV'], r['m_I_over_Hi'], r['phiI_S1_MP'], r['phiI_end_MP'],
            r['g_I_S1'], r['g_I_end'], r['K_end_gamma05']))
        if 'osc' in r:
            print('      osc: turn1 %.4f turn2 %.4f V/V0 %.3f w1/mI %.3f w range %.1f-%.1f H_i' % (r['osc']['first_turn_phiI_MP'], r['osc']['second_turn_phiI_MP'],
                  r['osc']['V_first_turn_over_V0'], r['osc']['omega_first_over_mI'], r['osc']['omega_min_over_Hi'], r['osc']['omega_max_over_Hi']))
        for k in ('exit_factor_Hi1.29e11', 'yy_harmonics_Hi1.29e11', 'floquet_first', 'fragmentation', 'rough_reledger_E'):
            if k in ent:
                print('      %s: %s' % (k, json.dumps(ent[k], default=lambda x: float(x))[:600]))


if __name__ == '__main__':
    main()
