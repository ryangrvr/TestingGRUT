"""Frozen component launches and outward, root-weighted result assembly."""
from fractions import Fraction as Q
from pathlib import Path
import json
import os
import subprocess
import sys
import bounds as E

BASE=Path(__file__).resolve().parent
S7=BASE.parent/'wm2_positive_residual_20261017'


def read(name):return json.loads((BASE/name).read_text())


def combined(rows,weights):
    if not all(r['status']=='COMPLETE' and r['physical_probability_intervals'] for r in rows):return None
    return {event:tuple(E.OLD.rational_bounds(sum(w*Q.from_float(r['physical_probability_intervals'][event][j])
        for w,r in zip(weights,rows)))[j] for j in range(2)) for event in ['macro','full']}


def root_rows(prefix):return [read('results/'+prefix+'_N5_root'+str(i)+'.json') for i in range(2)]


def recovery(prefix):
    rows=root_rows(prefix);bounds=combined(rows,[Q(3,5),Q(2,5)])
    gate=bool(bounds and all(r['finite_instance_gate'] for r in rows) and bounds['macro'][0]>=.5)
    return dict(gate_pass=gate,event_intervals=bounds,
        completed_components=sum(r['status']=='COMPLETE' for r in rows),
        endpoint_times=[r['completed_physical_steps_time'] for r in rows],
        statuses=[r['status'] for r in rows],resources=resources(rows))


def resources(rows):
    return dict(total_sequential_worker_wall_seconds=sum(r['whole_worker_wall_seconds'] for r in rows),
        total_sequential_worker_CPU_seconds=sum(r['whole_worker_cpu_seconds'] for r in rows),
        maximum_worker_peak_MiB=max(r['peak_resident_MiB'] for r in rows))


def launch(kind,counts,name):
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
    if kind=='s7':
        command=[sys.executable,str(S7/'worker.py'),'--method','s7','--counts',counts,
            '--output',str(BASE/'replay'/name)]
        cwd=S7
    else:
        command=[sys.executable,str(BASE/('worker_global.py' if kind=='g1' else 'worker.py')),
            '--counts',counts,'--output',str(BASE/'replay'/name)]
        cwd=BASE
    subprocess.run(command,cwd=cwd,env=env,check=True)


def gate_six(directory='results'):
    rows=[read(directory+'/g1_N5_root'+str(i)+'.json') for i in range(2)]
    interval=combined(rows,[Q(3,5),Q(2,5)])
    assert interval and interval['macro'][0]>=.5
    assert all(r['finite_instance_gate'] and r['analytic_state_error_account']<=.60 for r in rows)
    return True


def summary():
    result=dict(parent='a1dd220ee28cf8fd985a7c144922e0bb6c299902',
        original_charter='0e42a5d3023c7ec45029f38ef918c47f62b6c603',
        G1_charter='b15fcf0618ac53e2a5a69f216ac346e16fa3f6ff',
        original_S8=recovery('s8'),G1=recovery('g1'))
    if (BASE/'results/s7_N5_root0.json').exists() and (BASE/'results/s7_N5_root1.json').exists():
        rows=root_rows('s7');interval=combined(rows,[Q(3,5),Q(2,5)])
        result['fresh_S7']=dict(event_intervals=interval,resources=resources(rows),
            completed_components=sum(r['status']=='COMPLETE' for r in rows),equal_error_speed_claim=False,
            worker_wall_budget_seconds=600)
    six=BASE/'results/g1_N6_root0.json'
    if six.exists():
        gate_six();r=json.loads(six.read_text())
        result['six_event']=dict(status=r['status'],target_time=144,
            completed_physical_steps_time=r['completed_physical_steps_time'],
            stop_reason=r['stop_reason'],event_intervals=r['physical_probability_intervals'],
            resources=resources([r]),original_equal_root_weights_and_proved_mirror=True)
    else:result['six_event']=dict(status='NOT_RUN',reason='Requires completed G1 recovery')
    result.update(full_independent_disorder_scaling=False,dimension_selected=False,
        original_CP12_matter_universality=False,distinctive_empirical_prediction=False,
        external_independent_validation=False,
        arithmetic='roundoff-aware under declared IEEE/reduction model; not formal software certification')
    (BASE/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':print(json.dumps(summary(),indent=2))
