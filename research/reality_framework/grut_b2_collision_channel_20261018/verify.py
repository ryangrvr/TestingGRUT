"""Verify frozen sources, historical pins, packet and resource accounting."""
import hashlib,json,pathlib
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    manifest=json.loads((HERE/'MANIFEST.json').read_text())['files']
    for n,h in manifest.items():assert digest(HERE/n)==h,n
    frozen=json.loads((HERE/'EXECUTED_SOURCES.json').read_text())
    for n,h in frozen.items():assert digest(HERE/n)==h,n
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    for n,h in pins.items():assert digest(ROOT/n)==h,n
    d=json.loads((HERE/'CHECKS.json').read_text())
    r=json.loads((HERE/'REPLAY.json').read_text())
    assert r['same_nonperformance_json']
    assert r['whole_worker_wall_total_seconds']<60
    assert max(r['primary']['peak_rss_MiB'],r['replay']['peak_rss_MiB'])<1024
    g=d['two_composites']['gates']
    assert g['outgoing']=='FAIL' and g['encounter']=='PASS' and g['contact_sensitivity']=='FAIL'
    assert d['two_composites']['outgoing_probability_enclosure'][1]<.1
    assert d['two_composites']['colored_matrix_entry_defect']==0
    assert d['interior_fuel_exchange']['shortest_primitive_word']==8
    assert all(d['baseline_pins_unchanged'].values())
    size=sum(p.stat().st_size for p in HERE.iterdir() if p.is_file())
    assert size<100*1024**2
    print(json.dumps({'integrity':'PASS','packet_files':len(manifest),
        'frozen_sources':len(frozen),'historical_pins':len(pins),'bytes':size,
        'physical_gate':'FAIL as frozen','scope':'Integrity and accounting, not independent validation.'},indent=2))

if __name__=='__main__':main()
