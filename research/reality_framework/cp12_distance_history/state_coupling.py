"""Actual state's direct-edge coupling vs full configuration dephasing.

This is an operator-resource diagnostic, not a finite-time Taylor prediction
or an energy-matched preparation comparison.
"""
from pathlib import Path
from math import log2
import argparse
import json,numpy as np
from basis import model,distance
from evolution import Evolution

P=Path(__file__).resolve().parent

def annotate(model):
    n=model['N'];q=model['Q'];M=n*(n-1)//2;ew=model['wait_state_error_upper'];pe=2*ew+ew**2
    for r in model['records']:
        p=r['class_probability'];lo=p-pe;hi=p+pe
        g=r['global_operator_norm_square']**.5
        num=r['coherent_cross_norm_square']*p;error=2*num**.5*g*ew+(g*ew)**2
        dep=r['dephased_configuration_cross_norm_square']*p;dep_error=q/(100*M*M)*pe
        if lo>0:
            r['coherent_cross_norm_square_interval']=[max(0.,(num-error)/hi),(num+error)/lo]
            r['dephased_cross_norm_square_interval']=[max(0.,(dep-dep_error)/hi),(dep+dep_error)/lo]
    return model

def run(n,budget):
    q=n-1;sub,H,B,states,sizes,psi=model(n,q,budget)
    v,_,ew=Evolution(H,B,q/2,16)(psi,.1)
    D=distance(n,states,sub.es);present=(D==1).astype(float);records=[]
    removable=[];rawdiag=[]
    for g,bits,field in states:
        r=sum(bool(g>>e&1) and bool(bits>>(n+e)&1) and field[2*e]==field[2*e+1]==0 for e in range(sub.M))
        removable.append(r)
        rawdiag.append(q/2+2*(bits&((1<<n)-1)).bit_count()+sum(x*x for x in field)/8+r*(sub.M-q)/(10*sub.M))
    removable=np.array(removable,float);rawdiag=np.array(rawdiag)
    g2=q*(sub.M-q)/(100*sub.M**2)
    for r in sorted(set(D.tolist())):
        if r==1:continue
        a=(D==r)*v;pr=float(np.linalg.norm(a)**2);weights=abs(a)**2/pr
        amp=present*(H@a);coh=float(np.linalg.norm(amp)**2/pr)
        dep=float(weights@removable)/(100*sub.M**2)
        assert coh<=g2+1e-10
        assert dep<=q/(100*sub.M**2)+1e-12
        assert coh<=dep*(sub.M-q)+1e-10
        entropy=float(-sum(w*log2(w/size) for w,size in zip(weights,sizes) if w>0))
        energy=float(np.vdot(a,H@a).real/pr);dephenergy=float(weights@rawdiag)
        records.append(dict(initial_distance=r,disconnected=(r==n),class_probability=pr,
                            coherent_cross_norm_square=coh,dephased_configuration_cross_norm_square=dep,
                            amplification=coh/dep if dep else None,
                            removable_edge_fraction=float(weights@removable/q),
                            global_operator_norm_square=g2,
                            conditional_pure_energy=energy,configuration_dephased_energy=dephenergy,
                            configuration_dephased_entropy_bits=entropy))
    return annotate(dict(N=n,Q=q,budget=budget,waiting_time=.1,wait_state_error_upper=ew,records=records,
                scope='Direct edge-addition matrix-element norm in the supplied state, not a probability at finite time. Complete physical-configuration dephasing changes energy and entropy; no free coherence benefit is inferred.')
                )

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--annotate-only',action='store_true');a=ap.parse_args()
    if a.annotate_only:
        out=json.loads((P/'STATE_COUPLING.json').read_text());out['models']=[annotate(m) for m in out['models']]
    else:out=dict(status='CP12_STATE_COUPLING_AUDIT_PASS',models=[run(4,10),run(5,10),run(6,8)])
    (P/'STATE_COUPLING.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
