#!/usr/bin/env python3
"""Verify the saved packet, or explicitly repeat frozen measured workers.

No primary trajectories are launched by default. Full replay can be expensive.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

BASE=Path(__file__).resolve().parent


def read(name):return json.loads((BASE/name).read_text())


def verify():
    manifest=read('MANIFEST.json')
    for name,expected in manifest['sha256'].items():
        assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==expected,name
    for name,expected in manifest['dependencies'].items():
        assert hashlib.sha256((BASE.parent/name).read_bytes()).hexdigest()==expected,name
    for name,expected in read('EXECUTED_SOURCES.json')['sha256'].items():
        assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==expected,name
    controls=read('controls.json');assert controls['status']=='PASS'
    p=controls['physical_576_reference']
    assert p['physical_receiver_modes']==576 and p['physical_dimension']==333505
    assert p['full_wavefunction_error']<=p['complete_error_account']
    for row in controls['physical_full_wavefunctions']:
        assert row['full_wavefunction_error']<=row['complete_error_account']
    for row in controls['native_against_S5']:
        assert row['independent_interval_coefficient_error_upper']<=row['fast_coefficient_error_upper']
    summary=read('summary.json')
    components=[read('results/s7_N5_root0.json'),read('results/s7_N5_root1.json')]
    passed=all(r['status']=='COMPLETE' and r['analytic_state_error_account']<=.60 for r in components)
    if passed:
        lower=sum(w*r['physical_probability_intervals']['macro'][0] for w,r in zip([.6,.4],components))
        assert abs(lower-summary['recovery']['macro_probability_interval'][0])<2e-14
        passed=lower>=.5
    assert passed==summary['recovery']['gate_pass']
    six=BASE/'results/s7_N6_root0.json'
    if six.exists():assert passed,'Six-event experiment requires completed recovery'
    print(json.dumps(dict(status='PASS',verified_files=len(manifest['sha256']),
        recovery_gate=passed,saved_results_only=True),indent=2))


def run(method,counts,name):
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
    subprocess.run([sys.executable,'worker.py','--method',method,'--counts',counts,
        '--output','replay/'+name],cwd=BASE,env=env,check=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--controls',action='store_true')
    p.add_argument('--recovery',action='store_true');p.add_argument('--baselines',action='store_true')
    p.add_argument('--six',action='store_true');args=p.parse_args()
    if args.controls:
        import controls
        (BASE/'replay').mkdir(exist_ok=True)
        (BASE/'replay/controls.json').write_text(json.dumps(controls.all_checks(),indent=2)+'\n')
    if args.recovery:
        run('s7','2,2','s7_N5_root0.json');run('s7','3,1','s7_N5_root1.json')
    if args.baselines:
        run('s5','2,2','s5_N5_root0.json');run('s5','3,1','s5_N5_root1.json')
    if args.six:
        directory='replay' if args.recovery else 'results'
        rows=[read(directory+'/s7_N5_root0.json'),read(directory+'/s7_N5_root1.json')]
        assert all(x['status']=='COMPLETE' for x in rows)
        assert sum(w*x['physical_probability_intervals']['macro'][0] for w,x in zip([.6,.4],rows))>=.5
        run('s7','2,3','s7_N6_root0.json')
    if not any(vars(args).values()):verify()


if __name__=='__main__':main()
