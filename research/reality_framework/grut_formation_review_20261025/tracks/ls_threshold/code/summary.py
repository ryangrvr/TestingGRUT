#!/usr/bin/env python3
"""ls_threshold: collect the verdict numbers from results/*.json into results/summary.json (own code)."""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / 'results'
L = lambda n: json.loads((RES / n).read_text())
CODEX = {  # quoted from b448f25 CLAUDE_HANDOFF.md / THRESHOLD_RESULTS.json / BENCHMARK_RESULTS.json [Q]
    'R_L': {'10_1': [91.174, 91.246], '10_raw': [102.956, 103.037], '17_1': [83.115, 83.181], '17_raw': [93.393, 93.467]},
    'f': {'10_R_lower_raw': 10.8493, '17_R_lower_raw': 13.9423, '10_R_central_raw': 12.0800, '17_R_central_raw': 15.2658,
          '10_R_Y_raw': 0.570, '17_R_Y_raw': 5.750, '10_bound': 9.90, '17_bound': 12.80}}

out = dict(codex=CODEX)
pb = L('proxy_benchmark.json')
out['proxy_thresholds_own'] = {f"{t['candidate']}_{'1' if t['ymult'] == 1 else 'raw'}": t['bracket'] for t in pb['thresholds']}
out['proxy_f_own'] = {f"{p['candidate']}_{p['ledger']}_{'1' if p['ymult'] == 1 else 'raw'}_{'s' if p['s_channel'] else 'bench'}_N{p['N']}": p['f']
                      for p in pb['points']}
ours = L('ours_threshold.json')
out['our_earlier_code'] = dict(thresholds={f"{t['candidate']}_{'1' if t['ymult'] == 1 else 'raw'}": t['R_L_bracket'] for t in ours['thresholds']},
                               points={f"{p['candidate']}_{p['ledger']}_raw": p['f_extra'] for p in ours.get('points', [])},
                               note='endpoint T_vis = 1e-6 GeV, 260 steps, explicit update (code unchanged)')
for name in ('sensitivity.json', 'ke_consistent.json'):
    try:
        out[name.replace('.json', '')] = [{k: v for k, v in r.items() if k not in ('history', 'method', 'kw')} for r in L(name)['rows']]
    except FileNotFoundError:
        pass
(RES / 'summary.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps({k: v for k, v in out.items() if k in ('proxy_thresholds_own', 'our_earlier_code')}, indent=1))
