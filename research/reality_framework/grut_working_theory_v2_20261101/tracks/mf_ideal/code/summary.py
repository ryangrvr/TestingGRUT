#!/usr/bin/env python3
"""mf_ideal (Claude, own code): assemble the map S/Q(K, m_s, t_on, e) and the headline numbers from results/chain/*.json,
results/isocurvature.json, results/regression.json, results/hartree_check.json and results/bench/*.json.
Also the [D] arithmetic checks of the brief's estimates G1, G2 and G5 that this track can verify.
Usage: python3 -B summary.py   -> ../results/sq_map.json, ../results/summary.json
"""
import glob, json, math, os
import numpy as np
import mfi_common as mc
import chain_mfi as cm

ES = (1e-7, 1e-6, 3e-6, 1e-5, 1e-4)


def row(o):
    sq = {('%.0e' % r['e']): r['S_per_Q'] for r in o['rows']}
    L = cm.ledger(o['tag'])
    sqe = {}
    for e in ES:
        c = L['comps'](e); sqe['%g' % e] = cm.SQ(L['E_Ec'] * (1 + c['R_total']))
    sqe['inf'] = o['S_per_Q_e_inf']
    sp = o['spec']
    return dict(tag=o['tag'], K=o['K'], m_s=sp['m_s'], t_on=sp['t_on'], Gamma_s=sp['Gamma_s'], c0=sp['c0'], s_star_GeV=sp['s_star_GeV'],
                xi_s=sp.get('xi_s', 0.0), bg=sp['bg'].split('/')[-1], N_X=o['N_X'], eps=o['eps'], E_over_Ec=o['E_over_Ec'], N_YL=o['N_YL'],
                N_YH=o['N_YH'], N_sigma=o['N_sigma'], N_A_1e6=o['N_A_1e6'], R_vec_1e6=o['e1e-6']['R_vec'], R_Y=o['R_Y'],
                R_casc_1e6=o['e1e-6']['R_casc'], W_1e6=o['e1e-6']['W'], S_per_Q=sqe, branch_1e6=o['e1e-6']['branch'],
                branch_3e6=o['e3e-6']['branch'], e_crit_bulk=o['e_crit_bulk'], e_crit_Q2e5=o['e_crit_Q2e5'],
                relic_1e6=dict(M=o['e1e-6']['M'], p=o['e1e-6'].get('p'), H_i_req_GeV=o['e1e-6']['H_i_req_GeV'], r=o['e1e-6']['r']),
                relic_3e6=dict(M=o['e3e-6']['M'], p=o['e3e-6'].get('p'), H_i_req_GeV=o['e3e-6']['H_i_req_GeV'], r=o['e3e-6']['r']),
                Y_over_YQ_at_HM1=o['Y_over_YQ_at_HM1'], checks=o['checks'])


def group(tag):
    if tag.startswith('M1r'):
        return 'regression'
    if tag.startswith('smoke'):
        return 'bench_smoke'
    if tag.startswith('f_osc'):
        return 'oscillator_background'
    if tag.startswith('f_send'):
        return 'smooth_end'
    if tag.endswith('_xi'):
        return 'xi_s_tracking'
    if '_mirror' in tag:
        return 'mirror_sign'
    if '_s1e' in tag:
        return 's_star'
    if '_g0.05' in tag:
        return 'small_decay'
    if '_g0.5' in tag:
        return 'fast_decay'
    if '_c' in tag:
        return 'holder_amplitude'
    from grid import grid
    if tag in [j['tag'] for j in grid('P')]:
        return 'primary'
    return 'fine_K1'


def g_checks():
    """[D] arithmetic of G1, G2, G5 (the parts this track can check)."""
    HMP = mc.H_M1 / mc.MP
    out = {}
    out['G1_amplitude_over_MP_at_m3'] = math.sqrt(6) / 3                         # sqrt(6) H M_P / m
    v = math.sqrt(8 / 3); lam = 9 * HMP ** 2 / (2 * v ** 2)                      # m = sqrt(2 lam) v = 3H, lam v^4/4 = 3 M_P^2 H^2
    out['G2'] = dict(v_chi_over_MP=v, lambda_chi=lam, y2_loop_over_lam=(2 ** 4 / (16 * math.pi ** 2)) / lam,
                     y_v_over_half_m=2 * v / (1.5 * HMP), tachyonic_m2_over_H2=-lam * v ** 2 / HMP ** 2, H_i_GeV_used=mc.H_M1)
    amp = math.sqrt(cm.BETA_ISO * cm.AS / (1 - cm.BETA_ISO))
    out['G5'] = dict(s_star_min_per_unit_D_GeV_at_1p3e11=1.3e11 / (2 * math.pi * amp), energy_fraction_coeff_m3=9 / 6,
                     xi_s_for_m3=-9 / 12, note='s_* > |D| H_i/(2 pi sqrt(beta A_s/(1-beta))), beta = 0.038, uncorrelated; '
                                               'rho_s/rho = (m_s^2/6)(s_*/M_P)^2 at the slice')
    return out


def main():
    rows = []
    for f in sorted(glob.glob(str(mc.RES / 'chain' / 'chain_*.json'))):
        o = json.load(open(f))
        r = row(o); r['group'] = group(o['tag']); rows.append(r)
    by = {}
    for r in rows:
        by.setdefault(r['group'], []).append(r)
    iso = json.load(open(mc.RES / 'isocurvature.json'))
    reg = json.load(open(mc.RES / 'regression.json'))
    hc = json.load(open(mc.RES / 'hartree_check.json'))
    prim = by['primary']
    q6 = [r for r in rows if r['S_per_Q']['1e-06'] < cm.SQB]
    best = sorted([r for r in rows if r['group'] in ('primary', 'fine_K1', 'small_decay', 'xi_s_tracking', 'oscillator_background', 'smooth_end')],
                  key=lambda r: r['S_per_Q']['1e-06'])[:6]
    m1 = [r for r in rows if r['tag'] == 'M1r_bridge'][0]
    summ = dict(
        regression=dict(mode0_byte_identical_vs_rpdq=reg['R1']['byte_identical_vs_rpdq'], mode0_byte_identical_vs_rpchain=reg['R1']['byte_identical_vs_rpchain'],
                        NX_400_1000=reg['R1']['mine']['NX_400_1000'], table_bridge_rel_diff_dt1e3=reg['R2']['dt0.001']['rel_diff_vs_mode0'],
                        chain_M1r=dict(S_per_Q_1e6=m1['S_per_Q']['1e-06'], e_crit_bulk=m1['e_crit_bulk'], e_crit_Q2e5=m1['e_crit_Q2e5'],
                                       R_vec_1e6=m1['R_vec_1e6'], R_Y=m1['R_Y'])),
        primary_ranges=dict(S_per_Q_1e6=[min(r['S_per_Q']['1e-06'] for r in prim), max(r['S_per_Q']['1e-06'] for r in prim)],
                            eps=[min(r['eps'] for r in prim), max(r['eps'] for r in prim)],
                            E_over_Ec=[min(r['E_over_Ec'] for r in prim), max(r['E_over_Ec'] for r in prim)],
                            N_YL_ton0p1=[min(r['N_YL'] for r in prim if r['t_on'] == 0.1), max(r['N_YL'] for r in prim if r['t_on'] == 0.1)],
                            N_YL_ton1=[min(r['N_YL'] for r in prim if r['t_on'] == 1.0), max(r['N_YL'] for r in prim if r['t_on'] == 1.0)],
                            E_X=[min(json.load(open(mc.RES / 'release' / ('release_%s.json' % r['tag'])))['result']['E_X'] for r in prim),
                                 max(json.load(open(mc.RES / 'release' / ('release_%s.json' % r['tag'])))['result']['E_X'] for r in prim)],
                            N_A_1e6_by_K={('%g' % K): float(np.mean([r['N_A_1e6'] for r in prim if abs(r['K'] - K) < 1e-6])) for K in (1.0, 1.3, 1.6, 2.0, 3.0)},
                            n_free_X_1e6=sum(r['branch_1e6'] == 'free X' for r in prim), n=len(prim),
                            with_e_crit=[(r['tag'], r['e_crit_bulk']) for r in prim if r['e_crit_bulk']]),
        Q_matter_at_1e6=[(r['tag'], r['S_per_Q']['1e-06'], r['e_crit_bulk'], r['e_crit_Q2e5'], r['group']) for r in q6],
        best_at_c0_1=[(r['tag'], r['S_per_Q']['1e-06'], r['e_crit_bulk'], r['group']) for r in best],
        isocurvature={k: dict(D_1e6=v['e1e-6']['D'], D_charge_only=v['e1e-6']['D_charge_only'], s_star_min_over_Hi=v['e1e-6']['s_star_min_GeV_at_Hi']['per_H_i'],
                              D_3e6=v['e3e-6']['D'], branches_3e6=[v['e3e-6']['pm'][x]['branch'] for x in v['e3e-6']['pm']]) for k, v in iso.items()},
        hartree_check=dict(max_rel_dev=max(abs(x['rel']) for v in hc.values() for x in v if 'rel' in x),
                           W_max_checked=max(x['W'] for v in hc.values() for x in v if 'rel' in x)),
        G_checks=g_checks())
    bfiles = sorted(glob.glob(str(mc.RES / 'bench' / 'bench_run_*.json')))
    summ['bench_smoke'] = {os.path.basename(f)[10:-5]: {k: v for k, v in json.load(open(f)).items() if k in ('exit_factor', 's_check', 'table')}
                           for f in bfiles}
    json.dump(dict(note='mf_ideal map: every history is a model choice; S/Q from the closed form; branch thresholds 3.4716 (bulk), 2.9284 (Q = 2e5)',
                   groups=by), open(mc.RES / 'sq_map.json', 'w'), indent=1, default=float)
    json.dump(summ, open(mc.RES / 'summary.json', 'w'), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float)[:6000])


if __name__ == '__main__':
    main()
