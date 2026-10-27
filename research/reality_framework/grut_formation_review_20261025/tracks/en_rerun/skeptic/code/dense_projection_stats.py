"""en_rerun_skeptic -- own analysis of the dense (dtau = .5) projection history from own_engine.py, tau 600-800, N48.
Checks: (i) own-engine trace vs reproduced Codex trace; (ii) own tau-800 fields vs Codex-engine snapshot;
(iii) range/mean of the reference projections and floors with dense sampling vs en_rerun's 20-tau snapshots;
(iv) whole-orbit means (pericentre to pericentre of |chibar|^2).
usage: python3 -I dense_projection_stats.py <codex_state.tsv> <codex_tau800_snapshot> <own_800.bin> <out_json> dense1.json [dense2.json ...]"""
import sys, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from own_engine import read_snapshot

tsv, snap800, own800, outp = sys.argv[1:5]
rows = []
for f in sys.argv[5:]:
    r = json.load(open(f))['rows']
    rows += r if not rows else r[1:]
t = np.genfromtxt(tsv, names=True, delimiter='\t')
keys = ['rho', 'nQ', 'KX', 'Kphi', 'GX', 'Gphi', 'UX', 'Uphi', 'Umix', 'X2', 'phi2', 'charge_RMS', 'q0', 'q1', 'p0', 'p1']
mx = {k: 0. for k in keys}
for r in rows:
    i = np.argmin(abs(t['tau'] - r['tau']))
    assert abs(t['tau'][i] - r['tau']) < 1e-8
    for k in keys:
        mx[k] = max(mx[k], abs(r[k] - t[k][i])/max(abs(t[k][i]), 1e-300))
N, L, s1, q1, p1 = read_snapshot(snap800)
_, _, s2, q2, p2 = read_snapshot(own800)
field_diff = float(max(np.max(abs(q1 - q2))/np.max(abs(q1)), np.max(abs(p1 - p2))/np.max(abs(p1))))

P = [r for r in rows if 'refs' in r]
tau = np.array([r['tau'] for r in P])
ser = dict(massless=np.array([r['refs']['massless']['total'] for r in P]),
           hartree=np.array([r['refs']['hartree_codex']['total'] for r in P]),
           hartree_no_zero=np.array([r['refs']['hartree_no_zero_mode']['total'] for r in P]),
           radial=np.array([r['refs']['radial']['total'] for r in P]),
           massless_k_ge_1p8=np.array([r['refs']['massless']['k_ge_1p8'] for r in P]),
           hartree_k_ge_1p8=np.array([r['refs']['hartree_codex']['k_ge_1p8'] for r in P]),
           floor_fd4_shell_w=np.array([r['floors']['fd4_shell_w'] for r in P]),
           floor_fd4_shell_wg=np.array([r['floors']['fd4_shell_wg'] for r in P]),
           floor_mode_w=np.array([r['floors']['mode_w'] for r in P]),
           floor_mode_wg=np.array([r['floors']['mode_wg'] for r in P]),
           hartree_m2=np.array([r['refs']['hartree_codex']['m2'] for r in P]))
mr2 = np.array([r['mean_r2'] for r in P])
# pericentres of |chibar|^2 on the dense grid
peri = [i for i in range(1, len(mr2) - 1) if mr2[i] < mr2[i - 1] and mr2[i] <= mr2[i + 1] and mr2[i] == mr2[max(0, i - 20):i + 21].min()]
apo = [i for i in range(1, len(mr2) - 1) if mr2[i] > mr2[i - 1] and mr2[i] >= mr2[i + 1] and mr2[i] == mr2[max(0, i - 20):i + 21].max()]
snap_taus = np.arange(600., 800.001, 20.)
sidx = [int(np.argmin(abs(tau - x))) for x in snap_taus]
res = dict(own_engine_vs_codex_trace_maxrel_600_800=mx, own_vs_codex_tau800_field_maxrel=field_diff,
           n_projection_samples=len(P), projection_dtau=float(np.median(np.diff(tau))),
           pericentres=[float(tau[i]) for i in peri], apocentres=[float(tau[i]) for i in apo],
           mean_r2_range_600_800=[float(mr2.min()), float(mr2.max())])
summ = {}
for k, v in ser.items():
    sv = v[sidx]
    d = dict(dense_min=float(v.min()), dense_max=float(v.max()), dense_mean=float(v.mean()),
             tau_at_min=float(tau[v.argmin()]), tau_at_max=float(tau[v.argmax()]),
             snap20_min=float(sv.min()), snap20_max=float(sv.max()), snap20_mean=float(sv.mean()),
             at_800=float(v[-1]))
    if len(peri) >= 2:
        whole = (tau >= tau[peri[0]]) & (tau < tau[peri[-1]])
        d['whole_orbits_mean'] = float(v[whole].mean()); d['whole_orbits_span'] = [float(tau[peri[0]]), float(tau[peri[-1]])]
        d['per_orbit_means'] = [float(v[(tau >= tau[a]) & (tau < tau[b])].mean()) for a, b in zip(peri[:-1], peri[1:])]
    # correlation with orbit phase
    d['corr_with_mean_r2'] = float(np.corrcoef(v, mr2)[0, 1])
    summ[k] = d
res['series'] = summ
ratio = ser['massless']/ser['hartree']
res['massless_over_hartree'] = dict(dense_min=float(ratio.min()), dense_max=float(ratio.max()), dense_mean=float(ratio.mean()),
                                    snap20=[float(ratio[i]) for i in sidx])
res['hartree_over_shell_floor_w'] = dict(min=float((ser['hartree']/ser['floor_fd4_shell_w']).min()), max=float((ser['hartree']/ser['floor_fd4_shell_w']).max()))
res['hartree_over_shell_floor_wg'] = dict(min=float((ser['hartree']/ser['floor_fd4_shell_wg']).min()), max=float((ser['hartree']/ser['floor_fd4_shell_wg']).max()))
res['dense_rows'] = [dict(tau=float(tau[i]), mean_r2=float(mr2[i]), **{k: float(v[i]) for k, v in ser.items()}) for i in range(len(P))]
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('series', 'dense_rows')}, indent=0))
for k, d in summ.items():
    print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in d.items() if kk not in ('per_orbit_means', 'whole_orbits_span')},
          'per-orbit', [round(x, 3) for x in d.get('per_orbit_means', [])])
