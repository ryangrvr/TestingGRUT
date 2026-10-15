"""Original-neutral formation and complete core-record resources."""
from pathlib import Path
from math import factorial,sqrt
import argparse,json,time
import numpy as np
from shared import *


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--waits',default='2,5');ap.add_argument('--dt',type=float,default=1/32)
    ap.add_argument('--batch',type=int,default=4)
    ap.add_argument('--backend',choices=['matrix_free','sparse'],default='matrix_free')
    args=ap.parse_args();start=time.monotonic();n,q,b=args.n,args.q,args.budget
    m=quotient_model(n,q,b);H,B=m['H'],m['B']
    waits=sorted(float(t) for t in args.waits.split(','))
    propagator=MatrixFreePropagator if args.backend=='matrix_free' else Propagator
    result=propagator(H,{'B':B},q/2,dt=args.dt,batch=args.batch).run(m['psi'],waits)
    save_results(m,n,q,b,result,dict(dt=args.dt,batch=args.batch,backend=args.backend),start)


def save_results(m,n,q,b,result,meta,start):
    H,B=m['H'],m['B'];C=np.array([flow(n,s)[0] for s in m['states']],np.int64)
    yes={c:is_spanning_cycle(n,int(c)) for c in set(C)};mask=np.array([yes[c] for c in C])
    nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in m['states']]);count=factorial(n-1)//2
    rows=[]
    for wait,(w,er) in result.items():
        error=float(er['B'][0]);X=mask*w;p=float(np.linalg.norm(X)**2);pI=probability_interval(p,error)
        filling=float(np.sum(nv*abs(X)**2)/n);fillI=probability_interval(filling,error)
        mean=filling/p if p else None
        meanI=[fillI[0]/pI[1],min(1.,fillI[1]/pI[0])] if pI[0]>0 else [0.,1.]
        variance=mean*(1-mean) if mean is not None else None
        varianceI=[min(a*(1-a) for a in meanI),.25 if meanI[0]<=.5<=meanI[1] else max(a*(1-a) for a in meanI)]
        energy=-2*float(np.vdot(X,H@((1-mask)*w)).real)
        workerr=2*sqrt(2*q*min(n,q))*(2*error+error*error)+ALLOWANCE
        key=f'N{n}_Q{q}_B{b}_W{wait:g}'
        np.savez(cache_root()/f'{key}.npz',w=w,mask=mask,core=C,wait=wait,state_error=error)
        row=dict(wait=wait,cache_key=key,formation_probability=p,formation_interval=pI,
            state_error=error,source_mean=mean,source_mean_interval=meanI,source_variance=variance,
            source_variance_interval=varianceI,selected_parent_orbits=int(mask.sum()),
            physical_states_per_fixed_core=int(sum(m['sizes'][mask])//count),
            complete_record_outcomes=count+1,record_energy_change=energy,
            record_energy_interval=[energy-workerr,energy+workerr],
            record_entropy_bits=binary_entropy(p)+p*np.log2(count),
            record_entropy_interval_bits=complete_entropy_interval(pI,count))
        rows.append(row);print(json.dumps(dict(stage='core_formation',N=n,Q=q,**row)),flush=True)
    data=dict(status='PASS',N=n,Q=q,budget=b,dimension=H.shape[0],outside_dimension=B.shape[0],
        neutral_preparation='Original equal coherent neutral sea over all labelled fixed-Q hosts.',
        core_event='The entire directed-cycle edge mask is exactly one spanning simple cycle. Outside active connections, including one-leg charged ones, are not removed or filtered.',
        complete_record='One projector for each spanning core mask, plus a single failure projector. The core pointer is fixed before transport; the host is not separately measured.',
        record_work_identity='Different exact spanning-cycle records have no one-step H matrix element. Rewiring changes zero-field edges and preserves the cycle; a physical hop can preserve the core, destroy its cycle, or add an extra cyclic chord, but cannot change one spanning simple cycle directly into another. Complete-record work therefore equals binary success/failure work.',
        source='Gauge-invariant exp(-i theta n0). The rooted cycle symmetry makes the conditional source mean <nv>/N. Variance intervals include original waiting error.',
        rows=rows,**meta,seconds=time.monotonic()-start)
    (P/f'NEUTRAL_N{n}_Q{q}_B{b}.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
