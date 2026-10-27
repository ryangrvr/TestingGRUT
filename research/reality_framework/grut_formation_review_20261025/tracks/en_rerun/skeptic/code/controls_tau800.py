"""en_rerun_skeptic -- own projections/floors at tau 800 for all seven controls (en_rerun's reproduced final
snapshots, Codex engine), compared with Codex's PUBLISHED FIELD_RESULTS reference projections. Own code.
usage: python3 -I controls_tau800.py <FIELD_RESULTS.json> <early_workdir> <out_json> tag ..."""
import sys, json, os
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import proj_own as P
from own_engine import read_snapshot

F = json.load(open(sys.argv[1]))['cases']; work = sys.argv[2]; outp = sys.argv[3]
out = {}
for tag in sys.argv[4:]:
    snap = json.load(open(os.path.join(work, tag, 'run_result.json')))['snapshots'][-1]
    N, L, s, q, p = read_snapshot(snap)
    r = P.diagnostics(q, p, L, P.fd4_symbol_r(N, L), s, 800., proj=True)
    pub = {x['name']: x['raw_nminus_over_net_current'] for x in F[tag]['spectra']['reference_projections']}
    out[tag] = dict(N=N, L=L, massless=r['refs']['massless']['total'], hartree=r['refs']['hartree_codex']['total'],
                    published_massless=pub['massless'], published_hartree=pub['isotropic_Hartree_reference'],
                    massless_k_lt_1=r['refs']['massless']['k_lt_1'], massless_k_lt_1p8=r['refs']['massless']['k_lt_1p8'],
                    hartree_k_lt_1p8=r['refs']['hartree_codex']['k_lt_1p8'],
                    hartree_variants={k: v['total'] for k, v in r['refs'].items() if k != 'massless'},
                    floors=r['floors'], mean_r2=r['mean_r2'], X2=r['X2'])
    o = out[tag]
    print(tag, 'massless %.4f (pub %.4f) hartree %.4f (pub %.4f) k<1 %.3f floors w %.3f wg %.3f mode_wg %.3f' % (
        o['massless'], o['published_massless'], o['hartree'], o['published_hartree'], o['massless_k_lt_1'],
        o['floors']['fd4_shell_w'], o['floors']['fd4_shell_wg'], o['floors']['mode_wg']),
        {k: round(v, 3) for k, v in o['hartree_variants'].items()})
allH = [v for o in out.values() for v in o['hartree_variants'].values()]
res = dict(cases=out, hartree_all_variants_range=[min(allH), max(allH)],
           shell_floor_w_range=[min(o['floors']['fd4_shell_w'] for o in out.values()), max(o['floors']['fd4_shell_w'] for o in out.values())],
           shell_floor_wg_range=[min(o['floors']['fd4_shell_wg'] for o in out.values()), max(o['floors']['fd4_shell_wg'] for o in out.values())])
json.dump(res, open(outp, 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'cases'})
