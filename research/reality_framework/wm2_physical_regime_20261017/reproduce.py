"""Saved checks by default; explicit replays preserve all archived outcomes."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

BASE=Path(__file__).resolve().parent


def read(name):return json.loads((BASE/name).read_text())


def verify():
    manifest=read('MANIFEST.json')
    for name,value in manifest['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==value,name
    for name,value in manifest['dependencies'].items():assert hashlib.sha256((BASE.parent/name).read_bytes()).hexdigest()==value,name
    for name,value in read('EXECUTED_SOURCES.json')['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==value,name
    for name in ['SWEEP_CONTROLS.json','FORMATION_CONTROLS.json','TRANSPORT_CONTROLS.json']:
        assert read(name)['status']=='PASS'
    controls=read('SWEEP_CONTROLS.json')['independent_physical_and_arithmetic']
    for row in controls['full_wavefunctions']+[controls['physical_576_reference']]:
        assert row['full_wavefunction_error']<=row['complete_error_account']
    assert controls['physical_576_reference']['physical_dimension']==333505
    import formation_controls
    assert formation_controls.constants()['claimed_uniform_lower']=='1/2'
    summary=read('summary.json')
    assert summary['B']['N0']==2**28 and not summary['external_validation']
    assert not list((BASE/'results').glob('*N6*'))
    for i in range(2):
        p=BASE/f'results/S9_N5_root{i}.json'
        if not p.exists():continue
        row=json.loads(p.read_text());assert row['nonroot_type_counts']==[[2,2],[3,1]][i]
        assert row['physical_receiver_modes']==1600 and row['observation_time']==100
        assert row['wall_limit_seconds']==1200
        if row['status']=='COMPLETE':
            assert row['completed_steps']==row['steps'] and row['analytic_state_error_account']<=.60
        else:
            assert row['physical_probability_intervals'] is None
            assert row['computed_normalized_event_intervals'] is None
    print(json.dumps(dict(verification_status='PASS',saved_files=len(manifest['sha256']),
        default_launches_trajectories=False,numerical_gate=summary['A']['status'],
        analytic_status=summary['B']['status']),indent=2))


def main():
    p=argparse.ArgumentParser()
    for name in ['sweep-controls','analytic-controls','primary']:p.add_argument('--'+name,action='store_true')
    args=p.parse_args();replay=BASE/'replay';replay.mkdir(exist_ok=True)
    if any(vars(args).values()):verify()
    if args.sweep_controls:
        import sweep_controls
        (replay/'SWEEP_CONTROLS.json').write_text(json.dumps(sweep_controls.run(),indent=2)+'\n')
    if args.analytic_controls:
        import formation_controls,transport_controls
        for name,function in [('FORMATION_CONTROLS.json',formation_controls.run),('TRANSPORT_CONTROLS.json',transport_controls.run)]:
            (replay/name).write_text(json.dumps(function(),indent=2)+'\n')
    if args.primary:
        for i,counts in enumerate(['2,2','3,1']):
            output=replay/f'S9_N5_root{i}.json';assert not output.exists(),'Refuse to overwrite '+str(output)
            subprocess.run([sys.executable,str(BASE/'worker.py'),'--counts',counts,'--output',str(output),
                '--seconds','1200'],cwd=BASE,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)
    if not any(vars(args).values()):verify()


if __name__=='__main__':main()
