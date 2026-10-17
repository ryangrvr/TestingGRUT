"""Saved checks by default; explicit, gated fresh-worker replays only."""
import argparse
import hashlib
import json
from pathlib import Path
import campaign

BASE=Path(__file__).resolve().parent


def read(name):return json.loads((BASE/name).read_text())


def verify():
    manifest=read('MANIFEST.json')
    for name,value in manifest['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==value,name
    for name,value in manifest['dependencies'].items():assert hashlib.sha256((BASE.parent/name).read_bytes()).hexdigest()==value,name
    for record in ['EXECUTED_SOURCES.json','G1_EXECUTED_SOURCES.json']:
        for name,value in read(record)['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==value,name
    assert read('controls.json')['status']=='PASS' and read('G1_CONTROLS.json')['status']=='PASS'
    import global_controls
    global_controls.source_check()
    for record in ['controls.json','G1_CONTROLS.json']:
        r=read(record);r=r.get('physical_and_arithmetic',r)
        for row in r['full_wavefunctions']:
            assert row['full_wavefunction_error']<=row['complete_error_account']
        p=r['physical_576_reference'];assert p['physical_dimension']==333505
        assert p['full_wavefunction_error']<=p['complete_error_account']
    six=BASE/'results/g1_N6_root0.json'
    if six.exists():campaign.gate_six()
    auxiliary=read('SWEEP_CONTROLS.json')
    assert auxiliary['status']=='PASS' and not auxiliary['used_in_S8_or_G1']
    assert auxiliary['coherent_return_counterexample']['actual_error']>auxiliary['coherent_return_counterexample']['invalid_all_time_square_sum']
    summary=read('summary.json') if (BASE/'summary.json').exists() else None
    if summary:
        for key,prefix in [('original_S8','s8'),('G1','g1')]:
            assert summary[key]['gate_pass']==campaign.recovery(prefix)['gate_pass']
            assert summary[key]['event_intervals']==campaign.recovery(prefix)['event_intervals']
        if not summary['G1']['gate_pass']:assert not six.exists()
        assert summary['six_event']['status']=='NOT_RUN'
        for prefix in ['s8','g1','s7']:
            rows=campaign.root_rows(prefix)
            assert [r['nonroot_type_counts'] for r in rows]==[[2,2],[3,1]]
            for row in rows:
                assert row['observation_time']==100
                assert row.get('physical_receiver_modes',row.get('receiver_modes'))==1600
                if row['status']=='COMPLETE':
                    assert row['completed_physical_steps_time']==100
                    assert row['finite_instance_gate'] and row['analytic_state_error_account']<=.60
                    for interval in row['physical_probability_intervals'].values():
                        assert 0<=interval[0]<=interval[1]<=1
                else:
                    assert row['physical_probability_intervals'] is None
                    assert not row['finite_instance_gate']
                    assert row['completed_physical_steps_time']<100
        reference=campaign.combined(campaign.root_rows('s7'),[campaign.Q(3,5),campaign.Q(2,5)])
        assert summary['fresh_S7']['event_intervals']=={k:list(v) for k,v in reference.items()}
    print(json.dumps(dict(verification_status='PASS',files=len(manifest['sha256']),saved_records_only=True,
        original_recovery_gate=(summary['original_S8']['gate_pass'] if summary else 'pending'),
        G1_recovery_gate=(summary['G1']['gate_pass'] if summary else 'pending')) ,indent=2))


def main():
    parser=argparse.ArgumentParser()
    for flag in ['controls','original','global','baselines','six']:parser.add_argument('--'+flag,action='store_true')
    args=parser.parse_args();(BASE/'replay').mkdir(exist_ok=True)
    if args.controls:
        import controls,global_controls
        (BASE/'replay/controls.json').write_text(json.dumps(controls.run(),indent=2)+'\n')
        (BASE/'replay/G1_CONTROLS.json').write_text(json.dumps(global_controls.run(),indent=2)+'\n')
    if args.original:
        campaign.launch('s8','2,2','s8_N5_root0.json');campaign.launch('s8','3,1','s8_N5_root1.json')
    if getattr(args,'global'):
        campaign.launch('g1','2,2','g1_N5_root0.json');campaign.launch('g1','3,1','g1_N5_root1.json')
    if args.baselines:
        campaign.launch('s7','2,2','s7_N5_root0.json');campaign.launch('s7','3,1','s7_N5_root1.json')
    if args.six:
        campaign.gate_six('replay' if getattr(args,'global') else 'results')
        campaign.launch('g1','2,3','g1_N6_root0.json')
    if not any(vars(args).values()):verify()


if __name__=='__main__':main()
