"""Electric circulation can coexist with holes; measure it separately."""
from pathlib import Path
import argparse,json,time
import numpy as np
from loops import P,UP,model,fixed_triangle
from references import energy_basis,inclusion
from persistence import case


def run(snapshot_dir):
    start=time.monotonic();full,energies,V,weights,beta=energy_basis()
    source=model(4,4,256,source=True);J,checks=inclusion(full,source)
    D,cycle=fixed_triangle(4,source[3]);trajectory=json.loads((UP/'N4_Q4_B256_T20.json').read_text())
    cases=[]
    for wait in [5.,10.]:
        v=np.load(Path(snapshot_dir)/f'N4_Q4_B256_T{wait:g}.npz')['wave']
        err=min(trajectory['rows'],key=lambda a:abs(a['T']-wait))['state_error_upper']
        cases.append(case(f'neutral_sea_wait_{wait:g}',(J@v)[:,None],source,D,err,trajectory,wait,
            initial_channels=['occupied_flux','vacuum_flux','electric_circulation'],additional_detectors=['electric_circulation']))
    for label,w in weights.items():
        X=J@(V*np.sqrt(w)[None,:])
        cases.append(case(label,X,source,D,0.,trajectory,0.,energies=energies,
            initial_channels=['occupied_flux','vacuum_flux','electric_circulation'],additional_detectors=['electric_circulation']))
    out=dict(N=4,Q=4,budget=256,named_triangle=[1,2,3],cases=cases,seconds=time.monotonic()-start,
        scope='A directed cycle of opposite-sign nonzero half-leg fields carries an extractable integer circulation, even on a hole edge. This is stricter than protected support and broader than entirely occupied-edge flux. Conditional probabilities have the same diagnostic filtering and resource caveats as PERSISTENCE.json; not stationary currents or a causal cone.')
    (P/'CIRCULATION_SURVIVAL.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(saved='CIRCULATION_SURVIVAL.json',seconds=out['seconds'])),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
