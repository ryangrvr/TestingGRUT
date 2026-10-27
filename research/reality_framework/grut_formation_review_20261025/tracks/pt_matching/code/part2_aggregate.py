"""Aggregate the transition scan into results/transition_scan.json (summary per orbit, eps, seed)."""
import sys, os, json, glob, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from pm_common import *

OUT = os.path.join(SCRATCH, 'scan')


def load(pattern='*.json'):
    rows = []
    for fn in sorted(glob.glob(os.path.join(OUT, pattern))):
        r = json.load(open(fn)); r['file'] = os.path.basename(fn); rows.append(r)
    return rows


def stats(a):
    a = np.asarray(a, float)
    return dict(n=len(a), mean=a.mean(), std=a.std(ddof=1) if len(a) > 1 else 0.0, sem=a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0,
                min=a.min(), max=a.max(), median=np.median(a), q10=np.quantile(a, 0.1), q90=np.quantile(a, 0.9))


def summarize(rows, verbose=True):
    out = []
    for r in rows:
        r['variant'] = r['file'].split('_seed')[1].split('_', 1)[1].replace('.json', '')
    keys = sorted(set((r['orbit'], r['growth'], r['variant'], r['eps']) for r in rows), key=lambda k: (k[0], k[1], k[2], -k[3]))
    for orb, gr, var, eps in keys:
        sel = [r for r in rows if r['orbit'] == orb and r['growth'] == gr and r['eps'] == eps and r['variant'] == var]
        if not sel: continue
        J0 = sel[0]['J0_per_q']
        ent = dict(orbit=orb, growth=['exponential', 'RD_linear'][gr], variant=var, eps=eps, J0_per_q=J0, A0_per_q=sel[0]['A0_per_q'],
                   nminus=stats([r['nminus'] for r in sel]), Nphi=stats([r['Nphi'] for r in sel]),
                   nminus_plus_Nphi=stats([r['nminus'] + r['Nphi'] for r in sel]),
                   ratio_nminus_to_J0=stats([r['nminus'] / J0 for r in sel]),
                   M_escape_over_Mperi=stats([r['M_escape_over_Mperi'] for r in sel if r['M_escape_over_Mperi'] is not None]),
                   max_qrel=max(r['qrel'] for r in sel), max_work_rel=max(r['work_rel'] for r in sel),
                   by_seed={})
        for sd in sorted(set(r['seed'] for r in sel)):
            ss = [r for r in sel if r['seed'] == sd]
            ent['by_seed'][f'{sd:.0e}'] = dict(nminus=stats([r['nminus'] for r in ss]), Nphi=stats([r['Nphi'] for r in ss]),
                                               M_escape_over_Mperi=stats([r['M_escape_over_Mperi'] for r in ss if r['M_escape_over_Mperi'] is not None]),
                                               by_phase=[dict(phase=r['phase'], nminus=r['nminus'], Nphi=r['Nphi']) for r in sorted(ss, key=lambda x: x['phase'])])
        out.append(ent)
        if verbose:
            n = ent['nminus']; N = ent['Nphi']; e = ent['M_escape_over_Mperi']
            print(f"{orb:6s} {ent['growth'][:3]} {var:12s} eps {eps:7.0e} n={n['n']:3d} J0 {J0:.4f} | n-/q mean {n['mean']:.4f} sem {n['sem']:.4f} [{n['min']:.4f},{n['max']:.4f}] "
                  f"| Nphi/q mean {N['mean']:.4f} [{N['min']:.4f},{N['max']:.4f}] | sum {ent['nminus_plus_Nphi']['mean']:.4f} | Mesc/Mperi {e.get('mean', float('nan')):.3f} [{e.get('min', float('nan')):.3f},{e.get('max', float('nan')):.3f}]")
            for sd, v in ent['by_seed'].items():
                print(f"        seed {sd}: n- {v['nminus']['mean']:.4f} +- {v['nminus']['sem']:.4f}  Nphi {v['Nphi']['mean']:.4f}  Mesc {v['M_escape_over_Mperi'].get('mean', float('nan')):.3f}")
    return out


if __name__ == '__main__':
    rows = load()
    summ = summarize(rows)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        save_json('transition_scan.json', dict(summary=summ))
