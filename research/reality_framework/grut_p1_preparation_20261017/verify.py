"""Packet integrity and optional bounded replay, not mathematical validation."""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-checks',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    for name,digest in manifest['sha256'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    executed=json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for name,digest in executed['revision']['sha256'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name,digest in pins['sources'].items():
        assert hashlib.sha256((HERE.parents[2]/name).read_bytes()).hexdigest()==digest,name
    observations=json.loads((HERE/'CHECKS.json').read_text())
    assert observations['source_pins_unchanged']==len(pins['sources'])
    assert observations['resources']['wall_seconds']<60
    assert observations['resources']['peak_rss_MiB']<1024
    assert observations['tensor']['constraint_commutator_max']<1e-12
    assert observations['fixed_path']['full_state_error']<1e-10
    g=observations['complete_small_geometry']
    assert g['receiver_modes']==2 and g['all_physical_number_sectors_retained']
    assert g['preparation_state_error']<=g['preparation_bound']
    assert g['operation_state_error']<=g['operation_bound']
    f=json.loads((HERE/'FINDINGS.json').read_text())
    assert not f['independent_validation'] and not f['autonomous_preparation']
    assert f['experimental_observations']==0
    if args.run_checks:
        from checks import run
        result=run()
        # run() contains the full assertions and deliberately does not write
        # the archived observations. Runtime/RSS and harmless roundoff may vary.
        assert result['fixed_path']['dimension']==observations['fixed_path']['dimension']
        assert result['complete_small_geometry']['dimension']==g['dimension']
    print(json.dumps({'integrity':'pass','owned_files':len(manifest['sha256']),
                      'baseline_pins':len(pins['sources']),'replayed':args.run_checks,
                      'independent_math_validation':False}))


if __name__=='__main__':
    main()
