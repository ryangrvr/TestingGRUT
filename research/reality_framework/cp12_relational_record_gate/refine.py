"""Detector-specific bound for the existing dense N6,Q12,T1 state.

The TV witness is chosen after inspecting the projected distribution, then
held fixed as an ordinary bounded observable. Its expectation bounds the
full-law TV from below, with no statistical fitting or physical intervention.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from analyze import P,UP,model,resources,setup
from evolution import Evolution


def run(snapshot_dir):
    started=time.monotonic();n=6;q=12;budget=10;end=1.
    sub,H,B,states,sizes,psi=model(n,q,budget)
    r,k,nv,e2=resources(n,q,sub,states);S=setup(n,q,states,sizes,r,nv,e2)
    v=np.load(Path(snapshot_dir)/'N6_Q12_B10_T1.npz')['wave']
    p=abs(v)**2
    hist=np.bincount(S['classes'],weights=p,minlength=len(S['catalog']['keys']))
    kernel,inverse=S['controls']['gauss_parent_resources']
    mass=np.bincount(inverse,weights=p,minlength=len(kernel))
    null=mass@kernel;sign=np.sign(hist-null)
    F=.5*(sign[S['classes']]-(kernel@sign)[inverse])
    tv=float(abs(hist-null).sum()/2)
    assert abs(p@F-tv)<1e-12 and np.max(abs(F))<=1
    names=['core_has_cycle','core_connected','gauss_parent_TV_witness']
    O=np.column_stack([S['actual'][:,S['names'].index('has_cycle')],
                       S['actual'][:,S['names'].index('connected')],F])
    W=O*v[:,None];d=np.zeros(len(names));E=Evolution(H,B,q/2,2)
    for j in range(10):
        W,de,_=E(W,.1,back=True);d+=de
        print(json.dumps(dict(stage='backward',step=j+1,seconds=time.monotonic()-started)),flush=True)
    old=json.loads((UP/'N6_Q12_B10_T1.json').read_text())
    e=old['rows'][-1]['state_error_upper']
    values=p@O;bound=2*e*d+e*e
    rows=[dict(detector=name,value=float(x),backward_error_upper=float(b),
               expectation_error_upper=float(err),
               interval=[max(0.,float(x-err)),min(1.,float(x+err))])
          for name,x,b,err in zip(names,values,d,bound)]
    rows[-1]['interpretation']='Interval for the fixed witness expectation; the lower endpoint is also a lower bound on the true core-shape TV from the declared Gauss/parent/resource null. It is not a two-sided refined TV interval.'
    out=dict(N=n,Q=q,budget=budget,T=end,forward_error_upper=e,
        numerical_snapshot_norm_drift=float(abs(p.sum()-1)),backward_leakage_subintervals=2,
        rows=rows,projected_TV=tv,TV_witness_class_signs=sign.tolist(),
        core_shape_canonical_masks=S['catalog']['keys'],
        bound='2 e_forward d_backward + e_forward^2 for Hermitian contractions commuting with the retained projector; outgoing full-law B retained. Analytic bound evaluated in ordinary floating point, not interval arithmetic.',
        null_extension='For invariant physical configurations outside the represented retained basis, the Gauss/resource null is the identity. Their contribution cancels in the actual-minus-null class distribution and the TV witness is zero there.',
        seconds=time.monotonic()-started)
    (P/'REFINED_DENSE_GEOMETRY.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
