"""One-sided positive-observable dual lower bound for closed circulation."""
from pathlib import Path
import argparse,json,time
import numpy as np
from loops import P,UP,model,observables
from evolution import Evolution


def run(snapshot_dir):
    started=time.monotonic();sub,H,B,states,sizes,psi=model(6,12,10)
    names,O,mm=observables(6,12,states)
    detector=O[:,names.index('directed_flux_cycle_any')]
    v=np.load(Path(snapshot_dir)/'N6_Q12_B10_T1.npz')['wave']
    value=float(abs(v)**2@detector);w=detector*v;d=0.
    E=Evolution(H,B,6.,2)
    for j in range(10):
        w,_,de=E(w,.1,back=True);d+=de
        print(json.dumps(dict(stage='dense_dual',step=j+1,seconds=time.monotonic()-started)),flush=True)
    e=json.loads((UP/'N6_Q12_B10_T1.json').read_text())['rows'][-1]['state_error_upper']
    lower=max(0.,value-e*d)**2/value if value else 0.
    upper=min(1.,value+2*e*d+e*e)
    basic=[max(0.,np.sqrt(value)-e)**2,min(1.,(np.sqrt(value)+e)**2)]
    out=dict(N=6,Q=12,budget=10,T=1.,detector='directed_flux_cycle_any',
        projected_probability=value,forward_error_upper=e,backward_error_upper=d,
        positive_dual_interval=[lower,upper],
        positive_amplitude_interval=basic,
        combined_interval=[max(lower,basic[0]),min(upper,basic[1])],
        positive_cross_term_lower=max(0.,value-2*e*d),
        scope='For a positive contraction O, |<Ov,psi-v>|<=e*d and Cauchy-Schwarz in the O form yield p_exact >= max(0,p_projected-e*d)^2/p_projected. Upper bound is p_projected+2ed+e^2. Same full-law leakage matrix as before; ordinary floating point, not formal roundoff intervals.',
        seconds=time.monotonic()-started)
    (P/'DENSE_CIRCULATION_BOUND.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
