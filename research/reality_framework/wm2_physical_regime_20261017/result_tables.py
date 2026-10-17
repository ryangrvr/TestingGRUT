"""Recorded finite costs; probability display endpoints round down."""
from pathlib import Path
from math import floor
import json

BASE=Path(__file__).resolve().parent


def run():
    lines=['# Frozen finite recovery outcomes','',
        'Model target time is 100 in receiver-hopping units. Worker wall/CPU',
        'seconds are computational resources. A stopped run earns no target-time',
        'probability. The ledger column describes accumulated costs only.','',
        '| Remaining counts / original weight | Status | Completed model time | Wall s | CPU s | Peak MiB | Ledger |',
        '|---|---|---:|---:|---:|---:|---:|']
    rows=[]
    for i,weight in enumerate(['3/5','2/5']):
        p=BASE/f'results/S9_N5_root{i}.json'
        if not p.exists():
            lines.append(f'| Root {i} / {weight} | Pending | — | — | — | — | — |');continue
        r=json.loads(p.read_text());rows.append(r)
        lines.append(f"| {r['nonroot_type_counts']} / {weight} | {r['status']} | {r['completed_physical_steps_time']:.6f} | {r['whole_worker_wall_seconds']:.3f} | {r['whole_worker_cpu_seconds']:.3f} | {r['peak_resident_MiB']:.3f} | {r['analytic_state_error_account']:.6f} |")
    lines+=['','| Counts | Orthogonal loss | Linear in-range defect | Strang | Field product | Stop |',
        '|---|---:|---:|---:|---:|---|']
    for r in rows:
        e=r['error_components']
        lines.append(f"| {r['nonroot_type_counts']} | {e['schmidt_truncation']:.8f} | {e['in_range_defect']:.3g} | {e['strang']:.8f} | {e['field_factorization']:.8f} | {r['stop_reason'] or 'Complete'} |")
    lines+=['','The worker limit is 1200 seconds, checked at complete physical steps;',
        'overshoot is retained in each raw JSON. Rank ceiling 96, RSS ceiling',
        '2048 MiB, and all error gates remain frozen. Retained tensor entries',
        'are not total memory. No equal-error/budget speed claim is made.','']
    for r in rows:
        interval=r['physical_probability_intervals']
        if interval:
            a=floor(interval['macro'][0]*1e6)/1e6;b=floor(interval['full'][0]*1e6)/1e6
            lines.append(f"Complete counts {r['nonroot_type_counts']}: macro lower {a:.6f}, full lower {b:.6f} (rounded down).")
    if (BASE/'summary.json').exists():
        s=json.loads((BASE/'summary.json').read_text());A=s['A']
        lines+=['',f"Full original-weight recovery: **{A['status']}**."]
        if A['event_intervals'] is None:lines.append('No combined T=100 family interval is claimed.')
        lines.append('No fallback or six-event target belongs to S9.')
    (BASE/'RESULT_TABLES.md').write_text('\n'.join(lines).rstrip()+'\n')


if __name__=='__main__':run()
