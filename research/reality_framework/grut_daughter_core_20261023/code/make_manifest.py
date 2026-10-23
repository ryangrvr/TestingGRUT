#!/usr/bin/env python3
"""Pin delivered data, including exact trace shards and consumed inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PARENT='../grut_portal_debris_window_20261023/'
INPUTS=['FIXED_YIELD_RESULTS.json','SOURCE_INPUTS.json','code/widths.py','code/profiles.py',
        'data/profile_Q200000_lambda0.03_refined.csv','data/hipso_provenance.json']
INPUTS += ['data/hipso_'+p+'_'+m+'.csv' for p in ['pi','K'] for m in ['DGL','BTPZ']]
INPUTS += ['data/foresee_'+x+'.txt' for x in ['ctau','br_pi_pi','br_K_K','br_e_e','br_mu_mu']]


def fingerprint(path):return dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def payload():
    paths=[p for p in ROOT.iterdir() if p.is_file() and p.name!='MANIFEST.json' and not p.name.startswith('.')]
    paths += [p for p in (ROOT/'code').iterdir() if p.is_file() and p.suffix in ['.py','.cpp']]
    runs=json.loads((ROOT/'CORE_RESULTS.json').read_text())['runs']
    for run in runs:
        traces=[]
        for part in run['trace_parts']:
            with np.load(ROOT/part) as saved:
                trace=saved['trace']
                assert trace.ndim==2 and trace.shape[1]==10 and len(trace)>0,part
                assert np.isfinite(trace).all(),part
                traces.append(trace)
        trace=np.concatenate(traces)
        assert len(trace)==run['trace_samples'],run['name']
        assert hashlib.sha256(trace.tobytes()).hexdigest()==run['trace_float64_sha256'],run['name']
        paths += [ROOT/run[k] for k in ['history_file','final_file']]
        paths += [ROOT/p for p in run['trace_parts']]
        paths += [ROOT/'data'/f"{run['name']}_controlled_particle_spectrum.csv"]
    assert len(paths)==len(set(paths))
    return dict(parent_commit='d71cc012f0e27ab947ca942675b441fe8979f9f6',
                packet_files={str(p.relative_to(ROOT)):fingerprint(p) for p in sorted(paths)},
                consumed_parent_files={PARENT+p:fingerprint(ROOT/PARENT/p) for p in INPUTS})


def main(check):
    target=ROOT/'MANIFEST.json';actual=payload()
    if check:assert json.loads(target.read_text())==actual,'Manifest mismatch'
    else:target.write_text(json.dumps(actual,indent=2)+'\n')
    print('Manifest matches' if check else 'Manifest written',len(actual['packet_files']),
          'packet files plus',len(INPUTS),'consumed parent files; manifest itself additional.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
