"""Build results/summary.json from the per-run result files and the raw run outputs in scratch.
Usage: python3 -I aggregate.py SCRATCH_RUNS_DIR RESULTS_DIR
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import summary  # noqa: E402

runs, res = Path(sys.argv[1]), Path(sys.argv[2])
out = {'archived_recompute': None, 'full_reruns': {}, 'reduced_and_new_runs': {}, 'own_topology': {}, 'wall_seconds': {}}
arch = json.loads((res / 'archived_recompute.json').read_text())
out['archived_recompute'] = arch['summary']
COUNT = ['fplus_04', 'fminus_04', 'fplus_05', 'fminus_05', 'fplus_06', 'fminus_06', 'genuine_plus', 'genuine_minus', 'min_S', 'max_S']
for f in sorted(runs.glob('*.cmd.json')):
    name = f.name[:-len('.cmd.json')]
    out['wall_seconds'][name] = round(json.loads(f.read_text())['wall_seconds'], 1)
    rf = res / f'rerun_{name}.json'
    if name.startswith('full_') and rf.exists():
        r = json.loads(rf.read_text())
        s = r['rerun_summary']
        out['full_reruns'][name[5:]] = {
            'count_columns_bitwise_identical': all(c in r['bitwise_identical_columns'] for c in COUNT),
            'topology_records_identical': r['topology_records_identical'],
            'max_rel_diff_reduction_columns': max([v['max_rel'] for k, v in r['differing_columns'].items()
                                                   if k not in ('mean_S', 'linear_chi_projection')] + [0.0]),
            'a_form': s['a_form'], 'end_fplus_05': s['end_fplus_05'], 'end_fminus_05': s['end_fminus_05'],
            'end_winding': [s.get('end_winding_plus_any'), s.get('end_winding_minus_any')],
            'charge_drift_rerun': s['charge_rel_drift_max'],
            'charge_drift_archive': r['archive_summary']['charge_rel_drift_max'],
            'budget_rerun': s['budget_pointwise'], 'budget_archive': r['archive_summary']['budget_pointwise']}
    elif not name.startswith('full_'):
        try:
            s = summary(str(runs / name))
        except Exception as e:  # incomplete run
            out['reduced_and_new_runs'][name] = {'error': str(e)}
            continue
        s.pop('fraction_history', None)
        s.pop('post_windings_any', None)
        out['reduced_and_new_runs'][name] = s
for f in sorted(res.glob('topo_*.json')):
    r = json.loads(f.read_text())
    out['own_topology'][f.name[5:-5]] = {k: {'plus': [v['plus']['any_winding'], v['plus']['largest_winding'], round(v['plus']['occupied_fraction'], 5), round(v['plus']['largest_fraction'], 5), v['plus']['components']],
                                             'minus': [v['minus']['any_winding'], v['minus']['largest_winding'], round(v['minus']['occupied_fraction'], 5), round(v['minus']['largest_fraction'], 5), v['minus']['components']]}
                                         for k, v in r.items() if isinstance(v, dict)}
    out['own_topology'][f.name[5:-5]]['a'] = r['a']
    if 'solver_record_match' in r:
        out['own_topology'][f.name[5:-5]]['solver_record_match'] = r['solver_record_match']
out['legend_own_topology'] = '[any_winding, largest_winding, occupied_fraction, largest_fraction, components] per threshold S>|cut| & F2<cut2'
(res / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1)[:6000])
