"""qm_transfer summary (own, Claude 2 Nov 2026): collects the scan, applies the constraint criteria, maps the allowed
region in (f_v, zeta) and the observable band, and writes results/*.json.  [C]/[E]
usage: summary.py SCAN_DIR VAL_TWOFLUID.json CONTROLS.json
Outputs (track results/): transfer_grid.json (all P(k), z = 0 and z = 3), constraints.json (metrics, effective
thermal masses per probe, verdicts), allowed_region.json (zeta limits per family), validation.json (WDM, two-fluid,
controls), qm_micro_eval.json.
"""
import glob, json, math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import constraints as C

RES = os.path.join(HERE, '..', 'results')


def m_eff(value, key, increasing_with_m):
    """thermal mass (Viel fit) with the same probe metric; None if CDM-like beyond 200 keV or below 0.5 keV"""
    f = lambda lm: C.metrics(C.viel_T2(math.exp(lm), C.KG))[key] - value
    lo, hi = math.log(0.5), math.log(200.0)
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None if abs(fhi) < abs(flo) else 0.5
    return math.exp(brentq(f, lo, hi, xtol=1e-4))


def main():
    scan, valfn, ctlfn = sys.argv[1:4]
    refs = C.references(scan)
    W = {}
    models = []
    for fn in sorted(glob.glob(os.path.join(scan, 'model_*.json'))):
        d = json.load(open(fn))
        if 'rows' not in d:
            continue
        ks = [r['k'] for r in d['rows'] if 'P' in r]
        Ps = [r['P']['1.0'] for r in d['rows'] if 'P' in r]
        P3 = [r['P']['0.25'] for r in d['rows'] if 'P' in r]
        if d['spec']['fam'] == 'W':
            W[d['spec']['m_keV']] = dict(k=ks, P=Ps, sigma_p=d['summary']['sigma_p'])
            continue
        r = C.ratio_on_grid(ks, Ps)
        mt = C.metrics(r)
        S = d['summary']
        models.append(dict(name=d['spec']['name'], fam=d['spec']['fam'], spec=d['spec'], inversion=d.get('inversion'),
                           f_v_mass=S['f_v_mass'], f_v_charge=S.get('f_v_charge'), sigma_p=S['sigma_p'],
                           a_f=S.get('a_f_kin') or S.get('a_f_eq'), a_coll_end=S.get('a_coll_end'), T_f_eq=S.get('T_f_eq'),
                           T_f_kin=d.get('T_f_kin'), zeta=S.get('zeta'), zeta_eff_f=S.get('zeta_eff_f'), Q=S.get('Q'),
                           SQ=S.get('SQ'), g=S.get('g'), k=ks, P_z0=Ps, P_z3=P3, metrics=mt,
                           k_half=C.half_mode(ks, Ps), judge_viel=C.judge(mt, refs, 'viel'), judge_own=C.judge(mt, refs, 'own'),
                           mixed=C.mixed_dm_crosscheck(S['f_v_mass'], S['sigma_p'])))
    # effective thermal masses per probe (Viel-fit scale, the scale of the published bounds)
    for m in models:
        mt = m['metrics']
        m['m_eff_keV'] = dict(lya=m_eff(mt['dA'], 'dA', False), sat=m_eff(mt['N_gt_1e8'], 'N_gt_1e8', True),
                              lens=m_eff(mt['dn_1e7p5'], 'dn_1e7p5', True))
    json.dump(dict(k_Mpc=models[0]['k'] if models else [], note='P = P_model/P_LCDM, same primordial amplitude; z = 0 and z = 3',
                   models=[{x: m[x] for x in ('name', 'fam', 'spec', 'f_v_mass', 'f_v_charge', 'sigma_p', 'a_f', 'T_f_eq',
                                              'T_f_kin', 'zeta', 'zeta_eff_f', 'Q', 'SQ', 'g', 'P_z0', 'P_z3')} for m in models]),
              open(os.path.join(RES, 'transfer_grid.json'), 'w'), indent=1, default=float)
    json.dump(dict(scope=C.__doc__.strip().splitlines()[0], references={k: {x: y for x, y in v.items() if x != 'P'} for k, v in refs.items()},
                   models=[{x: m[x] for x in ('name', 'fam', 'f_v_mass', 'sigma_p', 'zeta', 'metrics', 'k_half', 'm_eff_keV',
                                              'judge_viel', 'judge_own', 'mixed')} for m in models]),
              open(os.path.join(RES, 'constraints.json'), 'w'), indent=1, default=float)
    # allowed region: for each family line (fixed SQ or fixed (f_c, T_f)), interpolate m_eff in log zeta
    lines = {}
    for m in models:
        sp = m['spec']
        if m['fam'] in ('P', 'S', 'E'):
            key = '%s_SQ%g_s%g' % (m['fam'], sp['SQ'], sp['s'])
        elif m['fam'] == 'X':
            key = 'X_fc%g_Tf%g' % (sp['fc'], sp['Tf'])
        else:
            continue
        lines.setdefault(key, []).append(m)
    region = {}
    thr = dict(lya_now=5.7, lya_2017=5.3, lya_cons=3.5, sat_now=6.5, lens_now=6.1, lens_2020=5.2, fut10=10.0, fut20=20.0)
    for key, ms in sorted(lines.items()):
        ms.sort(key=lambda x: x['zeta'])
        lz = np.log10([x['zeta'] for x in ms])

        def crossing(probe, mth):
            # largest zeta with m_eff >= mth (allowed); m_eff decreases with zeta
            vals = np.array([(x['m_eff_keV'][probe] if x['m_eff_keV'][probe] is not None else 1e3) for x in ms])
            ok = vals >= mth
            if ok.all():
                return dict(zeta_max='> %.0e' % ms[-1]['zeta'])
            if not ok.any():
                return dict(zeta_max='< %.0e' % ms[0]['zeta'])
            j = np.where(~ok)[0][0]
            if j == 0:
                return dict(zeta_max='< %.0e' % ms[0]['zeta'])
            x = np.interp(math.log(mth), [math.log(vals[j]), math.log(vals[j - 1])], [lz[j], lz[j - 1]])
            return dict(zeta_max=float(10**x))
        reg = dict(f_v_mass=[x['f_v_mass'] for x in ms], zeta=[x['zeta'] for x in ms], sigma_p=[x['sigma_p'] for x in ms],
                   m_eff_lya=[x['m_eff_keV']['lya'] for x in ms], m_eff_sat=[x['m_eff_keV']['sat'] for x in ms],
                   m_eff_lens=[x['m_eff_keV']['lens'] for x in ms])
        for probe in ('lya', 'sat', 'lens'):
            for nm, mth in (('now', {'lya': 5.7, 'sat': 6.5, 'lens': 6.1}[probe]), ('10keV', 10.0), ('20keV', 20.0)):
                reg['zeta_%s_%s' % (probe, nm)] = crossing(probe, mth)['zeta_max']
        region[key] = reg
    json.dump(dict(note='zeta_<probe>_now: largest zeta still allowed by the current bound (Lya 5.7 keV area criterion, '
                        'MW satellites 6.5 keV N(>1e8), lensing 6.1 keV dn/dlnM(10^7.5)); zeta_<probe>_10keV/20keV: '
                        'smallest zeta a probe reaching 10/20 keV would see (models with zeta above it are detectable). '
                        'm_eff: Viel-fit thermal mass with the same probe metric.', thresholds_keV=thr, lines=region),
              open(os.path.join(RES, 'allowed_region.json'), 'w'), indent=1, default=float)
    # validation file
    val = dict()
    if os.path.exists(valfn):
        val['fx_history_vs_two_fluid'] = json.load(open(valfn))
    if os.path.exists(ctlfn):
        val['controls'] = json.load(open(ctlfn))
    wd = {}
    for mk, w in sorted(W.items()):
        k = np.array(w['k'])
        al = 0.0437 * mk**(-1.188) / C.H
        va = (1 + (al * k)**(2 * 1.049))**(-10 / 1.049)
        wd['%g' % mk] = dict(k=w['k'], P_own=w['P'], P_viel=list(C.viel_T2(mk, k)), P_VA22_CLASSfit=list(va),
                             k_half_own=C.half_mode(w['k'], w['P']), k_half_viel=C.half_mode(list(k), list(C.viel_T2(mk, k))),
                             k_half_VA22=C.half_mode(list(k), list(va)), sigma_p=w['sigma_p'],
                             dA_own=refs.get('own_%g' % mk, {}).get('dA'), dA_viel=refs.get('viel_%g' % mk, {}).get('dA'))
    val['thermal_wdm_own_vs_fits'] = wd
    val['murgia_dA_reproduction'] = dict(viel_3p5=refs['viel_3.5']['dA'], viel_5p3=refs['viel_5.3']['dA'],
                                         murgia_quoted=dict(m3p5=0.38, m5p3=0.21))
    json.dump(val, open(os.path.join(RES, 'validation.json'), 'w'), indent=1, default=float)
    # qm_micro evaluation
    qm = [m for m in models if m['fam'] == 'M']
    json.dump(dict(note='tracks/qm_micro points (NOTES section 3), energy-dependent sticking from stick_thermal.json',
                   models=[{x: m[x] for x in ('name', 'spec', 'f_v_mass', 'f_v_charge', 'sigma_p', 'a_f', 'T_f_eq', 'T_f_kin',
                                              'zeta_eff_f', 'P_z0', 'metrics', 'm_eff_keV', 'judge_viel')} for m in qm]),
              open(os.path.join(RES, 'qm_micro_eval.json'), 'w'), indent=1, default=float)
    # print
    for key, reg in region.items():
        print('%-22s' % key, ' '.join('%s=%s' % (p, ('%.2e' % reg[p]) if isinstance(reg[p], float) else reg[p])
                                      for p in sorted(reg) if p.startswith('zeta_')))
    for m in models:
        me = m['m_eff_keV']
        print('%-28s fv=%.3f sigp=%.2e  m_eff lya/sat/lens = %s/%s/%s keV  now=%d' % (
            m['name'], m['f_v_mass'], m['sigma_p'], *(('%.1f' % me[x]) if me[x] else '>200' for x in ('lya', 'sat', 'lens')),
            m['judge_viel']['allowed_now']))
    for mk, w in wd.items():
        print('W %s keV: k_half own %.1f viel %.1f VA22 %.1f' % (mk, w['k_half_own'] or -1, w['k_half_viel'] or -1, w['k_half_VA22'] or -1))


if __name__ == '__main__':
    main()
