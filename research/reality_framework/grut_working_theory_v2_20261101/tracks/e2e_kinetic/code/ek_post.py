#!/usr/bin/env python3
"""e2e_kinetic post-processing (resumed run): scan table, reference numbers, the condensate transfer rate column.
  * adds Gamma_cond = -d ln(rho_cond a^3)/dt (condensate -> N + SM energy transfer rate, redshift removed) to
    ../results/fluid_HFa.json and fluid_HFb.json (computed from their own rows; centred differences in ln t);
  * ../results/scan.json: settings and summaries of every run in $SP (fx_* = after the shell factor-2 fix; others are
    pre-fix and kept for the record), plus the LCUT = 80 run in $SP/lcut80;
  * ../results/reference.json: decisive numbers, the f_rad / Gamma table and the m~ and H_i power laws.
Usage: python3 -B ek_post.py
"""
import json, math
import numpy as np
import ek_common as C


def add_gamma_cond(fn):
    d = json.load(open(fn))
    h = d['history']
    t = np.array(h['t'], float); a = np.array(h['a'], float); H = np.array(h['H'], float)
    fc = np.array(h['f_sig'], float) + np.array(h['f_chi'], float)
    rho = 3 * (C.MP_GEV / d['H_i_GeV']) ** 2 * H ** 2
    Ec = np.maximum(fc * rho * a ** 3, 1e-300)
    G = np.full(len(t), np.nan)
    ok = fc > 1e-6
    lt, lE = np.log(t), np.log(Ec)
    G[1:-1] = -(lE[2:] - lE[:-2]) / (t[2:] - t[:-2])
    G = np.where(ok, G, np.nan)
    h['Gamma_cond'] = [None if not np.isfinite(x) else float('%.6g' % x) for x in G]
    d['history_columns']['Gamma_cond'] = ('-d ln(rho_sig + rho_chi) a^3 / dt (H_i): condensate -> N + SM transfer rate, '
                                          'cycle-mean, redshift removed; None once the condensates are gone')
    C.dump(fn, d)
    return d


def table(d):
    h = d['history']; t = np.array(h['t'], float)
    rows = []
    for x in (0.19, 0.25, 0.3, 0.45, 1, 2.2, 3, 6, 10, 30, 100, 300, 1000, 3000, 5000, 1e4, 3e4, 1e5):
        i = min(int(np.searchsorted(t, x)), len(t) - 1)
        g = lambda k: h[k][i]
        rows.append(dict(t=t[i], f_rad=g('f_rad'), f_SM=g('f_SM'), f_N=g('f_N'), f_SM_thermal=g('f_SM_thermal'), f_sig=g('f_sig'),
                         f_chi=g('f_chi'), Gamma_th=g('Gamma_th'), Gamma_cond=g('Gamma_cond'), K_running=g('K_running'),
                         T_SM_GeV=g('T_SM') * d['H_i_GeV']))
    return rows


def plaw(x, y):
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


if __name__ == '__main__':
    fl = {nm: add_gamma_cond(C.RES / ('fluid_%s.json' % nm)) for nm in ('HFa', 'HFb')}
    runs = {}
    for fn in sorted(C.SP.glob('run_*.json')) + sorted((C.SP / 'lcut80').glob('run_*.json')):
        d = json.load(open(fn))
        inf = d['info']
        tag = fn.stem[4:] + ('_lcut80' if fn.parent.name == 'lcut80' else '')
        runs[tag] = dict(settings={k: inf.get(k) for k in ('name', 'H_GeV', 'mnu', 't0', 't_coh', 'lam_group', 'KLscale', 'cBBscale',
                                                          'G12', 'dshell', 'nqlog', 'dk_m', 'cdt', 'winv', 'lcut')},
                         summary={k: v for k, v in d['summary'].items() if k not in ('Gamma_th_at',)},
                         prefix=('after the shell factor-2 fix' if tag.startswith('fx_') else
                                 'pre-fix (shell rate x 1/2), superseded; kept for the record'),
                         restart_from=inf.get('restart_from'))
    C.dump(C.RES / 'scan.json', dict(label='C on the stated rates (model choices)', note=__doc__, runs=runs))
    S = lambda tag: runs[tag]['summary']
    mn = [(0.015, 'fx_natw_mnu0.015'), (0.05, 'fx_natw_nq6400'), (0.15, 'fx_natw_mnu0.15')]
    Hs = [(5.0e10, 'fx_natw_H5.0e10'), (1.3e11, 'fx_natw_nq6400'), (2.4e11, 'fx_natw_H2.4e11')]
    ref = dict(label='C on the stated rates (model choices; NOTES.md)', reference_tag='fx_natw_nq6400 (HFa), fx_natw_HFb (HFb)',
               HFa={k: S('fx_natw_nq6400')[k] for k in ('K', 'T_RH_GeV', 't_SM_half', 't_rad_half', 't_SM_99', 't_thermal',
                                                        't_bath_thermal', 'T_max_GeV', 't_Tmax', 'T_K_GeV', 't_SM_1em3')},
               HFb={k: S('fx_natw_HFb')[k] for k in ('K', 'T_RH_GeV', 't_SM_half', 't_rad_half', 't_SM_99', 't_thermal',
                                                     't_bath_thermal', 'T_max_GeV', 't_Tmax', 'T_K_GeV', 't_SM_1em3')},
               grid_check={t: {k: S(t)[k] for k in ('K', 'T_RH_GeV', 't_SM_half')} for t in ('fx_natw_nq6400', 'fx_natw_nq12800')},
               mnu_scaling=dict(points={str(m): {k: S(t)[k] for k in ('K', 'T_RH_GeV', 't_SM_half')} for m, t in mn},
                                exponent_K=plaw([m for m, _ in mn], [S(t)['K'] for _, t in mn]),
                                exponent_T_RH=plaw([m for m, _ in mn], [S(t)['T_RH_GeV'] for _, t in mn]),
                                exponent_t_SM_half=plaw([m for m, _ in mn], [S(t)['t_SM_half'] for _, t in mn]),
                                note='all four singlets at a common m~ (h_j^2 = m~ M_j / v^2), H_i = 1.3e11 GeV'),
               H_dependence={str(Hv): {k: S(t)[k] for k in ('K', 'T_RH_GeV', 't_SM_half')} for Hv, t in Hs},
               sensitivity={t: {k: S(t)[k] for k in ('K', 'T_RH_GeV', 't_SM_half')} for t in
                            ('fx_natw_x2', 'fx_g0.1', 'fx_g0.3', 'fx_natw_tcoh10', 'fx_dshell0') if t in runs},
               tables={nm: table(d) for nm, d in fl.items()})
    C.dump(C.RES / 'reference.json', ref)
    print(json.dumps({k: ref[k] for k in ('HFa', 'HFb', 'mnu_scaling', 'H_dependence', 'sensitivity')}, indent=1, default=float))
    for r in ref['tables']['HFa']:
        print(' '.join('%s %.4g' % (k, v) if isinstance(v, float) else '%s %s' % (k, v) for k, v in r.items()))
