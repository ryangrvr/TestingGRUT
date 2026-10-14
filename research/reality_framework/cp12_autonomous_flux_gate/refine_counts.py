"""Detector-specific source-free count bounds in the dense sector.

For O a contraction commuting with the retained-space projector, the real
expectation error is <= 2*e_forward*d_backward + e_forward^2. Backward
evolution is an error estimator, not an intervention on the physical state.
"""
from pathlib import Path
import json,sys,time
import numpy as np

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_unmasked_formation'))
sys.path.insert(0,str(P.parent/'cp12_distance_history'))
from run import model
from evolution import Evolution
from dynamics import resources


def run(n,q,budget,end,steps=8,only_flux=False):
    start=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget)
    r,k,nv,e2=resources(n,q,sub,states)
    E=Evolution(H,B,q/2,steps);v=psi.copy();e=0.
    for _ in range(round(end/.1)):
        v,_,de=E(v,.1);e+=de
    O=np.vstack([k/q] if only_flux else [r/q,k/q,nv/q])
    names=['flux_protected_fraction'] if only_flux else ['removable_fraction','flux_protected_fraction','hole_protected_fraction']
    W=O.T*v[:,None];d=np.zeros(len(names))
    for _ in range(round(end/.1)):
        W,dd,_=E(W,.1,back=True);d+=dd
    values=O@abs(v)**2
    return dict(N=n,Q=q,budget=budget,T=end,dimension=len(states),forward_error=e,leakage_subintervals=steps,
                rows=[dict(detector=name,value=float(val),backward_error=float(dual),
                           expectation_error_upper=float(2*e*dual+e*e),
                           interval=[max(0.,float(val-2*e*dual-e*e)),min(1.,float(val+2*e*dual+e*e))])
                      for name,val,dual in zip(names,values,d)],
                seconds=time.monotonic()-start)


def main():
    models=[run(6,12,6,.4),run(6,12,10,1.,steps=2,only_flux=True)]
    out=dict(status='CP12_DENSE_AUTONOMOUS_COUNT_REFINEMENT_PASS',models=models,
             scope='Unforced state, no external graph record. Analytically justified forward/backward detector bounds evaluated in ordinary floating point; no formal roundoff certificate and no long-time phase claim.')
    (P/'REFINED_COUNTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
