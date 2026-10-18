"""Integrity and declared execution checks, not proof validation."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE/'MANIFEST.json').read_text())['sha256']
    actual = {p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
    assert actual == set(manifest)
    for name, expected in manifest.items():
        assert sha(HERE/name) == expected, name
    for name in ('SOURCE_PINS.json','EXECUTED_SOURCES.json'):
        for path, expected in json.loads((HERE/name).read_text())['sha256'].items():
            assert sha(ROOT/path) == expected, path
    primary = json.loads((HERE/'PRIMARY_WORKER.json').read_text())
    replay = json.loads((HERE/'REPLAY.json').read_text())
    metrics = json.loads((HERE/'EXECUTION.json').read_text())
    data = json.loads((HERE/'CONTROLS.json').read_text())
    assert primary['exit_code'] == replay['exit_code'] == 0
    assert set(replay['data_files']) == {'CONTROLS.json','BULK_CONTROL.json',
                                       'SPECTRUM_8.json','SPECTRUM_12.json','SPECTRUM_16.json'}
    for name, check in replay['data_files'].items():
        assert check['byte_identical'] and sha(HERE/name) == check['sha256']
    assert replay['combined_whole_worker_seconds'] < 60
    assert max(metrics['peak_rss_MiB'],replay['metrics']['peak_rss_MiB']) < 1024
    assert data['no_source_time_evolution'] and data['no_operating_probabilities']
    assert data['source_and_proof_pins_unchanged']
    assert data['full_source_scaling_verdict'] == 'MATHEMATICALLY UNRESOLVED'
    assert [x['length'] for x in data['finite_source_controls']] == [8,12,16]
    total = sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())
    assert total < 25*1024**2
    print(json.dumps({'status':'PASS: integrity and protocol only',
                      'manifest_covered_files':len(manifest),
                      'total_files_including_manifest':len(manifest)+1,
                      'bytes':total,'scaling_source_verdict':'MATHEMATICALLY UNRESOLVED'}))


if __name__ == '__main__':
    main()
