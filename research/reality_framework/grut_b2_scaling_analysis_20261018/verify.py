"""Integrity/provenance only; does not validate the mathematical proofs."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    actual = sorted(p.name for p in HERE.iterdir()
                    if p.is_file() and p.name != 'MANIFEST.json')
    assert actual == sorted(manifest['sha256']), 'File inventory differs'
    for name, value in manifest['sha256'].items():
        assert sha(HERE / name) == value, name
    for name in ('SOURCE_PINS.json', 'EXECUTED_SOURCES.json'):
        for path, value in json.loads((HERE / name).read_text())['sha256'].items():
            assert sha(ROOT / path) == value, (name, path)
    c = json.loads((HERE / 'CONTROLS.json').read_text())
    replay = json.loads((HERE / 'REPLAY.json').read_text())
    worker = json.loads((HERE / 'PRIMARY_WORKER.json').read_text())
    metrics = json.loads((HERE / 'EXECUTION.json').read_text())
    assert worker['exit_code'] == replay['exit_code'] == 0
    assert replay['CONTROLS_byte_identical'] is True
    assert replay['CONTROLS_sha256'] == sha(HERE / 'CONTROLS.json')
    assert replay['combined_whole_worker_seconds'] < 60
    assert max(metrics['peak_rss_MiB'], replay['metrics']['peak_rss_MiB']) < 1024
    assert all(c['pins_unchanged'].values())
    assert [x['length'] for x in c['small_matrix_controls']] == [4, 6, 8]
    assert all(x['ordered_color_entry_defect'] == 0
               for x in c['small_matrix_controls'])
    total = sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())
    assert total < 25 * 1024 * 1024
    print(json.dumps({'status': 'PASS: integrity and declared execution checks only',
                      'manifest_covered_files': len(actual),
                      'files_including_manifest': len(actual) + 1,
                      'total_bytes': total,
                      'primary_data_sha256': sha(HERE / 'CONTROLS.json'),
                      'no_scaling_operating_verdict': True}))


if __name__ == '__main__':
    main()
