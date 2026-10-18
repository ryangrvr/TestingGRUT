"""Integrity/protocol checks; not a validation of analytic proofs."""
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
    for name, value in manifest.items():
        assert sha(HERE/name) == value, name
    for name in ('SOURCE_PINS.json', 'EXECUTED_SOURCES.json'):
        for path, value in json.loads((HERE/name).read_text())['sha256'].items():
            assert sha(ROOT/path) == value, path
    data = json.loads((HERE/'CONTROLS.json').read_text())
    replay = json.loads((HERE/'REPLAY.json').read_text())
    primary = json.loads((HERE/'PRIMARY_WORKER.json').read_text())
    metrics = json.loads((HERE/'EXECUTION.json').read_text())
    assert primary['exit_code'] == replay['exit_code'] == 0
    assert replay['data_byte_identical'] and replay['data_sha256'] == sha(HERE/'CONTROLS.json')
    assert replay['combined_whole_worker_seconds'] < 60
    assert max(metrics['peak_rss_MiB'], replay['metrics']['peak_rss_MiB']) < 1024
    assert [x['length'] for x in data['matrix_embedding_controls']] == [4, 6, 8]
    assert data['source_and_proof_pins_unchanged']
    assert data['scaling_operating_gate'].startswith('UNRESOLVED')
    total = sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())
    assert total < 25*1024**2
    print(json.dumps({'status':'PASS: integrity/protocol checks only',
                      'manifest_covered_files':len(manifest),
                      'total_files_including_manifest':len(manifest)+1,
                      'bytes':total,'data_sha256':sha(HERE/'CONTROLS.json'),
                      'scaling_operating_gate':'MATHEMATICALLY UNRESOLVED'}))


if __name__ == '__main__':
    main()
