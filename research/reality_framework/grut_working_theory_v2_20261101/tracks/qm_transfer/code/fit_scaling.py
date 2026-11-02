"""qm_transfer: parametric summary for evaluation at other tracks' (f_v, zeta, T_f) (own, Claude 2 Nov 2026).  [C fit]
(1) sigma_p (1D comoving velocity dispersion today, units of c) against zeta T_f (thermo and sudden, s = 1);
(2) effective thermal mass per probe, m_eff = A (sigma_p/1e-8)^b (f_v/0.1)^c, fitted over models with 1.5 < m_eff < 80 keV;
(3) the resulting zeta_max(f_v, T_f) for the current bounds and the 10 / 20 keV reach.
usage: fit_scaling.py  (reads ../results/constraints.json, transfer_grid.json; writes ../results/scaling_fit.json)"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
G = {m['name']: m for m in json.load(open(os.path.join(RES, 'transfer_grid.json')))['models']}
Cn = {m['name']: m for m in json.load(open(os.path.join(RES, 'constraints.json')))['models']}
out = {}
# (1) sigma_p vs zeta * T_f
for fam, mode in (('thermo', ('P', 'X')), ('sudden', ('S',))):
    x, y = [], []
    for n, m in G.items():
        if m['fam'] in mode and m['zeta'] and m['spec'].get('s', 1) == 1:
            Tf = m['spec'].get('Tf') or m['T_f_eq']
            x.append(math.log(m['zeta'] * Tf)); y.append(math.log(m['sigma_p']))
    b, a = np.polyfit(x, y, 1)
    out['sigma_p_fit_' + fam] = dict(form='sigma_p = A (zeta T_f)^b', A=math.exp(a), b=b,
                                     rms_dex=float(np.std(np.array(y) - (a + b * np.array(x))) / math.log(10)), n=len(x))
# (2) m_eff fits
for probe in ('lya', 'sat', 'lens'):
    X, Y = [], []
    for n, m in Cn.items():
        me = m['m_eff_keV'][probe]
        if me is None or not (1.5 < me < 80) or G[n]['fam'] not in ('P', 'S', 'E', 'X'):
            continue
        X.append([1.0, math.log(m['sigma_p'] / 1e-8), math.log(m['f_v_mass'] / 0.1)]); Y.append(math.log(me))
    X, Y = np.array(X), np.array(Y)
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    res = Y - X @ coef
    out['m_eff_fit_' + probe] = dict(form='m_eff[keV] = A (sigma_p/1e-8)^b (f_v/0.1)^c', A=math.exp(coef[0]), b=coef[1], c=coef[2],
                                     rms_dex=float(np.std(res) / math.log(10)), max_dex=float(np.max(abs(res)) / math.log(10)), n=len(Y))
# (3) zeta_max(f_v, T_f) (thermo, s = 1)
sp = out['sigma_p_fit_thermo']
tab = []
for fv in (0.01, 0.03, 0.06, 0.1, 0.2, 0.35, 0.5):
    for Tf in (0.004, 0.006, 0.009, 0.012):
        row = dict(f_v_mass=fv, T_f=Tf)
        for probe, mnow in (('lya', 5.7), ('sat', 6.5), ('lens', 6.1)):
            f = out['m_eff_fit_' + probe]
            for tag, mth in (('now', mnow), ('10keV', 10.0), ('20keV', 20.0)):
                sig = 1e-8 * (mth / (f['A'] * (fv / 0.1)**f['c']))**(1 / f['b'])
                row['zeta_%s_%s' % (probe, tag)] = (sig / sp['A'])**(1 / sp['b']) / Tf
        row['zeta_max_now'] = min(row['zeta_lya_now'], row['zeta_sat_now'], row['zeta_lens_now'])
        row['zeta_detect_20keV'] = min(row['zeta_lya_20keV'], row['zeta_sat_20keV'], row['zeta_lens_20keV'])
        row['zeta_detect_10keV'] = min(row['zeta_lya_10keV'], row['zeta_sat_10keV'], row['zeta_lens_10keV'])
        tab.append(row)
out['zeta_limits_from_fit'] = tab
out['note'] = ('fits valid for 1.5 < m_eff < 80 keV, s = 1 thermo kinetics; the physical M1r line has T_f ~ 0.004-0.009; '
               'models with zeta < zeta_detect_20keV are CDM-like for every probe considered')
json.dump(out, open(os.path.join(RES, 'scaling_fit.json'), 'w'), indent=1)
for k, v in out.items():
    if k != 'zeta_limits_from_fit':
        print(k, v)
for r in tab:
    if r['T_f'] in (0.004, 0.009):
        print('f_v=%.2f T_f=%.3f zeta_max_now=%.2e (lya %.2e sat %.2e lens %.2e)  detect10=%.2e detect20=%.2e' % (
            r['f_v_mass'], r['T_f'], r['zeta_max_now'], r['zeta_lya_now'], r['zeta_sat_now'], r['zeta_lens_now'],
            r['zeta_detect_10keV'], r['zeta_detect_20keV']))
