"""Generate factual cost/status tables from archived worker output only."""
from pathlib import Path
from math import floor
import json

BASE=Path(__file__).resolve().parent


def read(path):return json.loads(path.read_text())


def main():
    text=['# Recorded primary outcomes','',
        'Times are completed whole physical steps. A mid-sweep stop does not',
        'provide a stored state at that time or a target-time probability.',
        'Retained tensors/coefficients are not total resident memory.','',
        '| Protocol | Remaining counts | Status | Completed / required | Wall s | CPU s | Peak MiB | Accumulated ledger |',
        '|---|---|---|---:|---:|---:|---:|---:|']
    for prefix,label in [('s8','Original S8'),('g1','S8-G1'),('s7','Fresh S7')]:
        for root in range(2):
            path=BASE/f'results/{prefix}_N5_root{root}.json'
            if not path.exists():continue
            r=read(path);time=r['completed_physical_steps_time'];count=r.get('nonroot_type_counts')
            state=r['status']+(' (mid-sweep)' if not r.get('whole_step_state',True) else '')
            text.append(f"| {label} | {count} | {state} | {time:.6f} / 100 | {r['whole_worker_wall_seconds']:.3f} | {r['whole_worker_cpu_seconds']:.3f} | {r['peak_resident_MiB']:.3f} | {r['analytic_state_error_account']:.6f} |")
    text+=['','All limits, allocation choices and error components remain in the raw JSON.',
        'Original S8 has a 16-per-charge cap and effective total <=80 at N=5;',
        'G1 has a global cap 96. S8/G1 have 900 seconds; S7 has 600 seconds.',
        'These are different numerical protocols, not an equal-rank/equal-error',
        'speed comparison. A stopped run does not demonstrate a cheaper target.',
        'A mid-sweep ledger does not certify its partially updated state against',
        'physical H evolution at the last completed-step time.','']
    if (BASE/'summary.json').exists():
        s=read(BASE/'summary.json')
        text+=['| Complete N=5 family | Macro lower | Full-assembly lower | Recovery gate |',
            '|---|---:|---:|---|']
        for name,key in [('Original S8','original_S8'),('S8-G1','G1'),('Fresh S7','fresh_S7')]:
            row=s.get(key);interval=row.get('event_intervals') if row else None
            if interval:
                a=floor(interval['macro'][0]*1e6)/1e6;b=floor(interval['full'][0]*1e6)/1e6
                text.append(f'| {name} | {a:.6f} | {b:.6f} | {row.get("gate_pass","reference")} |')
            else:text.append(f'| {name} | No target-time bound | No target-time bound | Not passed |')
        text+=['','Displayed probability lower endpoints are rounded downward; exact',
            'endpoints and any upper bounds are in summary.json.','']
    (BASE/'RESULT_TABLES.md').write_text('\n'.join(text).rstrip()+'\n')


if __name__=='__main__':main()
