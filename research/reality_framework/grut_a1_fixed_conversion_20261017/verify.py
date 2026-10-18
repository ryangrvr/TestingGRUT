"""Integrity and bounded replay, not a proof checker."""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-checks',action='store_true')
    args=ap.parse_args()
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    for f,d in manifest['sha256'].items():
        assert hashlib.sha256((HERE/f).read_bytes()).hexdigest()==d,f
    executed=json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for f,d in executed['sha256'].items():
        assert hashlib.sha256((HERE/f).read_bytes()).hexdigest()==d,f
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for f,d in pins['sources'].items():
        assert hashlib.sha256((HERE.parents[2]/f).read_bytes()).hexdigest()==d,f
    observations=json.loads((HERE/'CHECKS.json').read_text())
    assert observations['resources']['wall_seconds']<60
    assert observations['resources']['peak_rss_MiB']<1024
    assert observations['resources']['written_bytes']<100*1024**2
    assert observations['coupled_cells']['many_cell_formation_lower']>.946
    assert observations['coupled_cells']['many_cell_return_upper']<.215<.25
    f=json.loads((HERE/'FINDINGS.json').read_text())
    assert f['decision']=='refuted by explicit permitted counterexample'
    assert not f['independent_validation'] and f['experimental_observations']==0
    assert f['physical_variants_run']==1 and f['rescue_variant_launched']==False
    replay=None
    if args.run_checks:
        from checks import run
        replay=run()
        assert replay['coupled_cells']['dimension']==observations['coupled_cells']['dimension']
        assert replay['moving_geometry']['exact_resource_dimension']==observations['moving_geometry']['exact_resource_dimension']
        assert replay['resources']['wall_seconds']+observations['resources']['wall_seconds']<60
    print(json.dumps({'integrity':'pass','owned_files':len(manifest['sha256']),
                      'source_pins':len(pins['sources']),'replayed':args.run_checks,
                      'replay_seconds':replay['resources']['wall_seconds'] if replay else None,
                      'independent_proof_validation':False}))


if __name__=='__main__': main()
