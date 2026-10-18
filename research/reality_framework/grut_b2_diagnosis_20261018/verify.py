"""Hash/accounting checks only; do not rerun physical trajectories."""
import hashlib,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    manifest=json.loads((HERE/'MANIFEST.json').read_text())['files']
    for name,sha in manifest.items():assert digest(HERE/name)==sha,name
    frozen=json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for name,sha in frozen.items():assert digest(HERE/name)==sha,name
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    for name,sha in pins.items():assert digest(ROOT/name)==sha,name
    run=json.loads((HERE/'EXECUTION.json').read_text())
    for name,sha in run['data_sha256'].items():assert digest(HERE/name)==sha,name
    replay=json.loads((HERE/'REPLAY.json').read_text())
    assert replay['all_ten_data_files_identical']
    assert replay['whole_worker_wall_total_seconds']<60
    assert max(replay['primary']['peak_rss_MiB'],replay['replay']['peak_rss_MiB'])<1024
    s=json.loads((HERE/'SECTORS.json').read_text());row=s['integer_time_points'][-1]
    old=json.loads((HERE.parent/'grut_b2_collision_channel_20261018/CHECKS.json').read_text())
    orig=old['two_composites']['time_points'][-1]
    assert abs(row['p_output_R']-orig['outgoing_probability'])<1e-12
    assert abs(row['p_q2']-orig['two_pairs_probability'])<1e-12
    assert abs(row['p_q0']+row['p_q1']+row['p_q2']-1)<1e-10
    assert abs(row['p_fuel11']-row['p_q0'])<1e-12
    assert abs(row['fuel_left']-row['p_fuel11']-row['p_fuel10'])<1e-12
    assert row['p_output_R']<.1
    b=json.loads((HERE/'CHANNEL_ACCOUNT.json').read_text())
    assert b['gamma_modulus_interval'][1]<1
    files_bytes=sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())
    assert files_bytes<100*1024**2
    print(json.dumps({'integrity':'PASS','packet_files':len(manifest),'frozen_sources':len(frozen),
        'historical_pins':len(pins),'bytes':files_bytes,'B2_failed_endpoint':'unchanged',
        'scope':'Hash/accounting verification, not external validation.'},indent=2))

if __name__=='__main__':main()
