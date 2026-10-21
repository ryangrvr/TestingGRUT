"""Independent reader for Codex's archived 3D diagnostics (a69bd86).

Does NOT import or execute any Codex analysis script. Reads raw/<run>.csv and
raw/<run>.topology.jsonl with numpy/json only and recomputes, per run:
  a_form   first scalar output with fplus_05 >= 0.10 and fminus_05 >= 0.10
  end f+/- fplus_05 / fminus_05 at the last scalar row
  charge   max_t |Q_com(t)/Q_com(0) - 1|
  budget   work residual  max_t |H(t)-H(0)-W(t)| / max(max_t H, H(0))   (variant A, as stated in
           Codex validate(); we also report variants B: / H(0), C: / max|H|, D: / abs_work)
  topology last-record plus/minus any_winding and largest_winding, and the first post-formation
           record where both phases wind all three axes, and whether this holds at every
           post-formation topology record.
Usage: python3 -I reader.py RAW_DIR RUN_TABLE.md OUT.json
"""
import json
import sys
from pathlib import Path

import numpy as np


def load_csv(path):
    with open(path) as f:
        header = f.readline().strip().split(',')
    data = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    return {h: data[:, i] for i, h in enumerate(header)}


def load_topo(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def analyse(raw, name, meta=None):
    d = load_csv(raw / f'{name}.csv')
    out = {'name': name, 'rows': int(len(d['t'])), 'a_end': float(d['a'][-1]), 't_end': float(d['t'][-1])}
    both = (d['fplus_05'] >= 0.10) & (d['fminus_05'] >= 0.10)
    idx = np.flatnonzero(both)
    out['a_form'] = float(d['a'][idx[0]]) if idx.size else None
    out['a_form_prev_output'] = float(d['a'][idx[0] - 1]) if idx.size and idx[0] > 0 else None
    out['end_fplus_05'] = float(d['fplus_05'][-1])
    out['end_fminus_05'] = float(d['fminus_05'][-1])
    out['end_genuine_plus'] = float(d['genuine_plus'][-1])
    out['end_genuine_minus'] = float(d['genuine_minus'][-1])
    out['max_fminus_05_any_time'] = float(d['fminus_05'].max())
    out['min_fplus_05_after_a1_5'] = float(d['fplus_05'][d['a'] >= 1.5].min())
    q = d['Q_com']
    out['Q_initial'] = float(q[0])
    out['charge_rel_drift_max'] = float(np.max(np.abs(q / q[0] - 1)))
    out['Q_negative_end_over_Q'] = float(d['Q_negative'][-1] / q[0])
    H, W = d['H_conf'], d['work_conf']
    res = H - H[0] - W
    out['budget_A_over_max_H_H0'] = float(np.max(np.abs(res)) / max(H.max(), H[0]))
    out['budget_B_over_H0'] = float(np.max(np.abs(res)) / H[0])
    out['budget_C_over_max_absH'] = float(np.max(np.abs(res)) / np.max(np.abs(H)))
    out['budget_argmax_a'] = float(d['a'][np.argmax(np.abs(res))])
    # variant E: pointwise denominator max(|H0|,|H(t)|,|W(t)|) (the definition Codex's analyze.py uses;
    # identified only AFTER our variant A disagreed in the 3rd significant digit)
    norm = np.maximum.reduce([np.full(len(H), abs(H[0])), np.abs(H), np.abs(W)])
    out['budget_E_pointwise'] = float(np.max(np.abs(res) / norm))
    out['H_max_over_H0'] = float(H.max() / H[0])
    if meta is not None and 'abs_work' in meta:
        out['budget_D_over_abs_work'] = float(np.max(np.abs(res)) / meta['abs_work'])
    out['H_end_over_H0'] = float(H[-1] / H[0])
    out['w_end_rho_p'] = float(d['p'][-1] / d['rho'][-1])
    out['v2_core_proxy_end'] = float(d['v2_core_proxy'][-1])
    out['v2_smooth_proxy_end'] = float(d['v2_smooth_proxy'][-1])
    topo_path = raw / f'{name}.topology.jsonl'
    if topo_path.exists() and topo_path.stat().st_size:
        tp = load_topo(topo_path)
        last = tp[-1]
        out['topo_records'] = len(tp)
        out['topo_last_a'] = float(last['a'])
        for s in ('plus', 'minus'):
            out[f'end_{s}_any_winding'] = int(last[s]['any_winding'])
            out[f'end_{s}_largest_winding'] = int(last[s]['largest_winding'])
            out[f'end_{s}_components'] = int(last[s]['components'])
            out[f'end_{s}_largest_fraction'] = float(last[s]['largest_fraction'])
            out[f'end_{s}_occupied_fraction'] = float(last[s]['occupied_fraction'])
        out['end_interface_links'] = int(last['interface_links'])
        out['end_basin_area_over_V'] = float(last['basin_link_area_comoving_over_V'])
        if out['a_form'] is not None:
            post = [r for r in tp if r['a'] >= out['a_form'] - 1e-12]
            out['post_formation_topo_records'] = len(post)
            out['post_formation_topo_a'] = [round(float(r['a']), 6) for r in post]
            out['both_largest_wind_all_at_every_post_record'] = bool(post) and all(
                r['plus']['largest_winding'] == 7 and r['minus']['largest_winding'] == 7 for r in post)
            out['both_any_wind_all_at_every_post_record'] = bool(post) and all(
                r['plus']['any_winding'] == 7 and r['minus']['any_winding'] == 7 for r in post)
            out['post_record_windings'] = [(r['plus']['any_winding'], r['minus']['any_winding'],
                                            r['plus']['largest_winding'], r['minus']['largest_winding'])
                                           for r in post]
    return out


def parse_table(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith('|')]
    head = [h.strip() for h in lines[0].strip('|').split('|')]
    rows = []
    for l in lines[2:]:
        cells = [c.strip() for c in l.strip('|').split('|')]
        rows.append(dict(zip(head, cells)))
    return rows


def rel(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def main():
    raw, table, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    rows = parse_table(table)
    results = []
    for r in rows:
        name = r['name']
        meta = json.loads((raw / f'{name}.meta.json').read_text())
        mine = analyse(raw, name, meta)
        mine['N_meta'] = meta['N']
        mine['L_meta'] = meta['N'] * meta['dx']
        cmp = {}
        # a_form: table prints 2 decimals or 'none'
        ta = None if r['a_form'] == 'none' else float(r['a_form'])
        cmp['a_form'] = {'table': r['a_form'], 'mine': mine['a_form'],
                         'agree': (ta is None and mine['a_form'] is None) or
                                  (ta is not None and mine['a_form'] is not None and abs(round(mine['a_form'], 2) - ta) < 1e-9)}
        for k_tab, k_mine in (('end_fplus', 'end_fplus_05'), ('end_fminus', 'end_fminus_05')):
            tv = float(r[k_tab])
            mv = mine[k_mine]
            # table prints 6 significant digits
            cmp[k_tab] = {'table': tv, 'mine': mv, 'agree_6sig': float(f'{mv:.6g}') == tv}
        for k_tab, k_mine in (('end_plus_winding', 'end_plus_any_winding'), ('end_minus_winding', 'end_minus_any_winding')):
            cmp[k_tab] = {'table': int(r[k_tab]), 'mine_any': mine.get(k_mine),
                          'mine_largest': mine.get(k_mine.replace('any', 'largest')),
                          'agree_any': int(r[k_tab]) == mine.get(k_mine)}
        tb = float(r['budget'])
        cmp['budget'] = {'table': tb, 'mine_A': mine['budget_A_over_max_H_H0'],
                         'agree_A_3sig': float(f"{mine['budget_A_over_max_H_H0']:.3g}") == tb,
                         'mine_E': mine['budget_E_pointwise'],
                         'agree_E_3sig': float(f"{mine['budget_E_pointwise']:.3g}") == tb,
                         'mine_B': mine['budget_B_over_H0'], 'mine_D': mine.get('budget_D_over_abs_work')}
        cmp['N'] = {'table': int(r['N']), 'mine': mine['N_meta'], 'agree': int(r['N']) == mine['N_meta']}
        cmp['L'] = {'table': float(r['L']), 'mine': mine['L_meta'], 'agree': abs(float(r['L']) - mine['L_meta']) < 1e-9}
        results.append({'run': name, 'mine': mine, 'compare': cmp})
    allagree = {k: all(x['compare'][k].get('agree', x['compare'][k].get('agree_6sig', x['compare'][k].get('agree_any'))) for x in results)
                for k in ('a_form', 'end_fplus', 'end_fminus', 'end_plus_winding', 'end_minus_winding', 'N', 'L')}
    allagree['budget_variantA_max_over_maxH_3sig'] = all(x['compare']['budget']['agree_A_3sig'] for x in results)
    allagree['budget_variantE_pointwise_3sig'] = all(x['compare']['budget']['agree_E_3sig'] for x in results)
    summary = {
        'max_charge_drift_all_rows': max(x['mine']['charge_rel_drift_max'] for x in results),
        'max_budget_A_all_rows': max(x['mine']['budget_A_over_max_H_H0'] for x in results),
        'max_budget_E_all_rows': max(x['mine']['budget_E_pointwise'] for x in results),
        'max_budget_B_over_H0_all_rows': max(x['mine']['budget_B_over_H0'] for x in results),
        'budget_E_agrees_3sig_all_rows': all(x['compare']['budget']['agree_E_3sig'] for x in results),
        'all_columns_agree': allagree,
    }
    Path(out).write_text(json.dumps({'summary': summary, 'rows': results}, indent=1) + '\n')
    print(json.dumps(summary, indent=1))
    for x in results:
        m, c = x['mine'], x['compare']
        print(f"{x['run']:22s} a_form={m['a_form']!s:6s}[{c['a_form']['table']}] f+={m['end_fplus_05']:.6g}[{c['end_fplus']['table']}] "
              f"f-={m['end_fminus_05']:.6g}[{c['end_fminus']['table']}] wind={m.get('end_plus_any_winding')}/{m.get('end_minus_any_winding')}"
              f"(lg {m.get('end_plus_largest_winding')}/{m.get('end_minus_largest_winding')})[{c['end_plus_winding']['table']}/{c['end_minus_winding']['table']}] "
              f"budA={m['budget_A_over_max_H_H0']:.4g} budE={m['budget_E_pointwise']:.4g}[{c['budget']['table']}] Hmax/H0={m['H_max_over_H0']:.3g} budD={m.get('budget_D_over_abs_work',0):.3g} budB={m['budget_B_over_H0']:.4g} dQ={m['charge_rel_drift_max']:.3g}")


if __name__ == '__main__':
    main()
