"""qm_transfer SKEPTIC: constraint metrics of the K_rsa = 800 / 1600 reruns against the production (K_rsa = 200)
grid, using the track's own constraints.py (read-only) so that only the transfer ratios change.  [C]
usage: sk_krsa_metrics.py OUT.json gridK_*.json"""
import json, os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
import constraints as C
RES = os.path.join(HERE, '..', '..', 'results')
G = {m['name']: m for m in json.load(open(os.path.join(RES, 'transfer_grid.json')))['models']}
KS = [1.5, 3.0, 5.0, 8.0, 12.0, 17.0, 24.0, 34.0, 48.0, 68.0, 96.0, 136.0]
refs = C.references(None)
masses = (3.5, 5.2, 5.3, 5.7, 6.1, 6.5, 10.0, 20.0)
def m_eff(val, key, inc):
    # Viel-fit mass with the same metric (log interpolation over a fine mass grid)
    ms = np.exp(np.linspace(np.log(1.0), np.log(150.0), 120))
    vs = np.array([C.metrics(C.viel_T2(m, C.KG))[key] for m in ms])
    if inc:   # metric increases with m (abundance ratios)
        if val >= vs[-1]: return None
        return float(np.exp(np.interp(val, vs, np.log(ms))))
    if val <= vs[-1]: return None
    return float(np.exp(np.interp(-val, -vs, np.log(ms))))
out = dict(scope=__doc__.strip().splitlines()[0], thresholds={k: {x: refs[k][x] for x in ('dA', 'N_gt_1e8', 'dn_1e7p5')}
                                                               for k in ('viel_5.7', 'viel_6.5', 'viel_6.1', 'viel_10', 'viel_20')}, runs={})
for fn in sys.argv[2:]:
    d = json.load(open(fn))
    for key, row in d.items():
        if key.startswith('cdm'):
            continue
        name, K = key.rsplit('_K', 1)
        Ps = [row[str(k)] for k in KS]
        P0 = G[name]['P_z0']
        r = {}
        for lab, P in (('prod_K200', P0), ('K' + K, Ps)):
            mt = C.metrics(C.ratio_on_grid(KS, P))
            r[lab] = dict(P=dict(zip(map(str, KS), P)), metrics=mt,
                          m_eff=dict(lya=m_eff(mt['dA'], 'dA', False), sat=m_eff(mt['N_gt_1e8'], 'N_gt_1e8', True),
                                     lens=m_eff(mt['dn_1e7p5'], 'dn_1e7p5', True)),
                          allowed_now=bool(mt['dA'] <= refs['viel_5.7']['dA'] and mt['N_gt_1e8'] >= refs['viel_6.5']['N_gt_1e8']
                                           and mt['dn_1e7p5'] >= refs['viel_6.1']['dn_1e7p5']))
        r['zeta'] = G[name]['zeta']; r['f_v_mass'] = G[name]['f_v_mass']
        out['runs'][key] = r
        a, b = r['prod_K200'], r['K' + K]
        print('%-26s K=%s  dA %.4f -> %.4f   N8 %.4f -> %.4f   dn %.4f -> %.4f   m_eff(lya) %s -> %s  sat %s -> %s  P(48) %.3f -> %.3f  P(96) %.3f -> %.3f' % (
            name, K, a['metrics']['dA'], b['metrics']['dA'], a['metrics']['N_gt_1e8'], b['metrics']['N_gt_1e8'],
            a['metrics']['dn_1e7p5'], b['metrics']['dn_1e7p5'], a['m_eff']['lya'] and round(a['m_eff']['lya'], 2),
            b['m_eff']['lya'] and round(b['m_eff']['lya'], 2), a['m_eff']['sat'] and round(a['m_eff']['sat'], 2),
            b['m_eff']['sat'] and round(b['m_eff']['sat'], 2), a['P']['48.0'], b['P']['48.0'], a['P']['96.0'], b['P']['96.0']))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
