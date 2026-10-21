"""Compare a rerun's scalar CSV and topology JSONL with an archived run (or summarize a run alone).
Usage:
  python3 -I compare.py full   RERUN_PREFIX ARCHIVE_PREFIX OUT.json   # same-configuration rerun
  python3 -I compare.py single RERUN_PREFIX OUT.json                  # reduced-cost run summary
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reader import load_csv, load_topo  # noqa: E402  (own reader)


def summary(prefix):
    d = load_csv(prefix + '.csv')
    tp = load_topo(prefix + '.topology.jsonl') if Path(prefix + '.topology.jsonl').exists() else []
    both = (d['fplus_05'] >= .1) & (d['fminus_05'] >= .1)
    ix = np.flatnonzero(both)
    q = d['Q_com']
    H, W = d['H_conf'], d['work_conf']
    res = H - H[0] - W
    norm = np.maximum.reduce([np.full(len(H), abs(H[0])), np.abs(H), np.abs(W)])
    out = {'a_form': float(d['a'][ix[0]]) if ix.size else None,
           'a_end': float(d['a'][-1]),
           'end_fplus_05': float(d['fplus_05'][-1]), 'end_fminus_05': float(d['fminus_05'][-1]),
           'end_genuine_plus': float(d['genuine_plus'][-1]), 'end_genuine_minus': float(d['genuine_minus'][-1]),
           'max_fminus_05': float(d['fminus_05'].max()), 'a_at_max_fminus_05': float(d['a'][np.argmax(d['fminus_05'])]),
           'max_fminus_04': float(d['fminus_04'].max()),
           'charge_rel_drift_max': float(np.max(np.abs(q / q[0] - 1))),
           'budget_pointwise': float(np.max(np.abs(res) / norm)),
           'budget_over_H0': float(np.max(np.abs(res)) / abs(H[0])),
           'w_end': float(d['p'][-1] / d['rho'][-1]),
           'v2_core_proxy_end': float(d['v2_core_proxy'][-1]),
           'rho_end': float(d['rho'][-1]),
           'area_core_end': float(d['area_core_com_over_V'][-1])}
    if tp:
        last = tp[-1]
        out['end_winding_plus_any'] = last['plus']['any_winding']
        out['end_winding_minus_any'] = last['minus']['any_winding']
        out['end_winding_plus_largest'] = last['plus']['largest_winding']
        out['end_winding_minus_largest'] = last['minus']['largest_winding']
        out['end_components'] = [last['plus']['components'], last['minus']['components']]
        if out['a_form'] is not None:
            post = [r for r in tp if r['a'] >= out['a_form'] - 1e-9]
            out['post_formation_records'] = len(post)
            out['both_wind_all_every_post_record'] = bool(post) and all(
                r['plus']['any_winding'] == 7 and r['minus']['any_winding'] == 7 for r in post)
            out['post_windings_any'] = [(r['plus']['any_winding'], r['minus']['any_winding']) for r in post]
    out['fraction_history'] = [(round(float(a), 4), float(fp), float(fm)) for a, fp, fm in
                               zip(d['a'], d['fplus_05'], d['fminus_05'])][::2]
    return out


def compare(re, ar):
    d1, d2 = load_csv(re + '.csv'), load_csv(ar + '.csv')
    n = min(len(d1['t']), len(d2['t']))
    out = {'rows_rerun': len(d1['t']), 'rows_archive': len(d2['t'])}
    same_t = np.array_equal(d1['t'][:n], d2['t'][:n]) and np.array_equal(d1['a'][:n], d2['a'][:n])
    out['identical_time_grid'] = bool(same_t)
    exact_cols, diffs = [], {}
    for k in d1:
        x, y = d1[k][:n], d2[k][:n]
        if np.array_equal(x, y):
            exact_cols.append(k)
        else:
            den = np.maximum(np.abs(y), 1e-300)
            diffs[k] = {'max_abs': float(np.max(np.abs(x - y))), 'max_rel': float(np.max(np.abs(x - y) / den))}
    out['bitwise_identical_columns'] = exact_cols
    out['differing_columns'] = diffs
    t1, t2 = load_topo(re + '.topology.jsonl'), load_topo(ar + '.topology.jsonl')
    out['topology_records'] = [len(t1), len(t2)]
    out['topology_records_identical'] = t1 == t2
    mism = []
    for r1, r2 in zip(t1, t2):
        if r1 != r2:
            mism.append({'a': r1['a'], 'rerun': r1, 'archive': r2})
    out['topology_mismatches'] = mism[:5]
    out['rerun_summary'] = summary(re)
    out['archive_summary'] = summary(ar)
    for k in ('fraction_history',):
        out['rerun_summary'].pop(k)
        out['archive_summary'].pop(k)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'full':
        r = compare(sys.argv[2], sys.argv[3])
        dest = sys.argv[4]
    else:
        r = summary(sys.argv[2])
        dest = sys.argv[3]
    Path(dest).write_text(json.dumps(r, indent=1) + '\n')
    print(json.dumps({k: v for k, v in r.items() if k not in ('fraction_history',)}, indent=1)[:4000])
