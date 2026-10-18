"""Integrity/accounting verification only; no trajectory execution."""
import hashlib, json, pathlib

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        assert digest(HERE/name)==expected, name
    frozen=json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for name,expected in frozen.items(): assert digest(HERE/name)==expected, name
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    for name,expected in pins.items(): assert digest(ROOT/name)==expected, name
    checks=json.loads((HERE/'CHECKS.json').read_text())
    replay=json.loads((HERE/'REPLAY.json').read_text())
    assert replay['exact_same_nonperformance_json']
    assert replay['total_worker_wall_seconds']<60
    assert max(replay['primary']['peak_rss_MiB'],replay['replay']['peak_rss_MiB'])<1024
    assert all(checks['baseline_pins_unchanged'].values())
    assert checks['scalar']['matter_uniform_lower']>.5
    assert checks['scalar']['arrival_lower']>.6
    assert sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())<100*1024**2
    print(json.dumps({'packet_files':len(manifest['files']),
        'frozen_sources':len(frozen),'historical_source_pins':len(pins),
        'integrity':'PASS','scope':'Hash/accounting verification, not independent mathematical validation.'},indent=2))

if __name__=='__main__': main()
