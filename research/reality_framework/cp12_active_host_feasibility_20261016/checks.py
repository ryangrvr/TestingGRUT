"""Independent labelled controls for records, persistence and graph constraints."""
from fractions import Fraction as F
from math import lcm,sqrt
import json,time
import numpy as np
from scipy.sparse import csr_matrix
from common import *
from graphs import metrics


def main():
    start=time.monotonic();n,q,b=4,5,32;full=physical_model(n,q,b)
    H,B,states=(full[k] for k in ['H','B','states']);M=6;host=states[0][0]
    chosen=[i for i,s in enumerate(states) if s[0]==host];scale=lcm(8,10*M)
    destinations={};cross_checks=0;maximum=F(0)
    for col in chosen:
        s=states[col];dest=[]
        removable=sum(bool(s[0]>>e&1 and s[1]>>(n+e)&1 and s[2][2*e]==s[2][2*e+1]==0) for e in range(M))
        for t,c in transitions(n,s,scale):
            if t[0]==host:continue
            gauss(n,t);assert t not in destinations;destinations[t]=col;dest.append(F(c,scale));cross_checks+=1
        norm=sum(c*c for c in dest);expected=F(removable*(M-q),100*M*M)
        assert norm==expected;maximum=max(maximum,norm)
    assert maximum==F(q*(M-q),100*M*M)
    rng=np.random.default_rng(61);v=rng.normal(size=len(states))+1j*rng.normal(size=len(states));v/=np.linalg.norm(v)
    nv=np.array([(s[1]&15).bit_count() for s in states]);parity=(-1.)**nv
    diag=np.array([q/2+2*a+sum(x*x for x in s[2])/8+
        sum(bool(s[0]>>e&1 and s[1]>>(n+e)&1 and s[2][2*e]==s[2][2*e+1]==0) for e in range(M))*(M-q)/(10*M)
        for s,a in zip(states,nv)])
    recorded=np.column_stack([v*np.array([s[0]==g for s in states]) for g in sorted({s[0] for s in states})])
    direct=float(np.vdot(recorded,H@recorded).real-np.vdot(v,H@v).real)
    parity_work=float(np.vdot(v,diag*v).real-.5*(np.vdot(v,H@v).real+np.vdot(parity*v,H@(parity*v)).real))
    assert abs(direct-parity_work)<1e-11
    neutral=(host,sum(1<<(n+e) for e in range(M)),(0,)*(2*M));j=states.index(neutral)
    v=np.zeros(len(states),complex);v[j]=1
    ids=np.array(chosen);Hg=H[ids,:][:,ids];Bg=B[:,ids];vv=v[ids]
    full_trajectory=Propagator(H,{'electric':B},q/2,dt=1/64,batch=8).run(v,[.125,.25,.5])
    fixed_trajectory=Propagator(Hg,{'electric':Bg},q/2,dt=1/64,batch=8).run(vv,[.125,.25,.5])
    rows=[];c=sqrt(q*(M-q))/(10*M)
    for t,(x,er) in full_trajectory.items():
        y,ey=fixed_trajectory[t];embedded=np.zeros(len(states),complex);embedded[ids]=y
        difference=float(np.linalg.norm(x-embedded));error=float(er['electric'][0]+ey['electric'][0])
        assert difference<=c*t+error+1e-10
        retention=float(np.linalg.norm(x[ids])**2)
        assert probability_interval(retention,float(er['electric'][0]))[1]>=1-(c*t)**2
        rows.append(dict(lag=t,retention=retention,universal_retention_lower=1-(c*t)**2,
            full_vs_frozen_state_difference=difference,state_difference_bound=c*t,combined_electric_error=error))
    structural=[]
    for nn,qq in [(5,10),(6,12)]:
        all_metrics=[metrics(nn,sum(1<<e for e in es)) for es in combinations(range(nn*(nn-1)//2),qq)]
        assert all(a['connected'] and a['diameter']<=2 for a in all_metrics)
        structural.append(dict(N=nn,Q=qq,hosts=len(all_metrics),diameters=sorted({a['diameter'] for a in all_metrics}),
            mean_balls=all_metrics[0]['balls'],mean_pair_distance=all_metrics[0]['mean_finite_pair_distance']))
    nn=64;edges=list(combinations(range(nn),2))
    def mask(pairs):return sum(1<<edges.index(tuple(sorted(p))) for p in set(tuple(sorted(p)) for p in pairs))
    ring=mask([(i,(i+d)%nn) for i in range(nn) for d in [1,2]])
    square=mask([(8*x+y,8*((x+1)%8)+y) for x in range(8) for y in range(8)]+
        [(8*x+y,8*x+(y+1)%8) for x in range(8) for y in range(8)])
    aa,bb=metrics(nn,ring,3),metrics(nn,square,3)
    assert aa['degrees']==bb['degrees']==[4]*nn and ring.bit_count()==square.bit_count()==2*nn
    assert aa['balls']==[1.,5.,9.,13.] and bb['balls']==[1.,5.,13.,25.]
    data=dict(status='PASS',labelled_dimension=len(states),raw_rewiring_checks=cross_checks,
        exact_cross_host_norm_square=str(maximum),host_record_work_identity_difference=abs(direct-parity_work),
        persistence_controls=rows,primary_fixed_edge_constraints=structural,
        equal_degree_counterexamples=dict(N=nn,Q=2*nn,ring_next_neighbor_balls=aa['balls'],square_balls=bb['balls']),
        interpretation='Supplied labelled control states and frozen-host comparisons verify bounds only; they do not count as formed spatial phases or full-law locality successes.',
        seconds=time.monotonic()-start)
    (P/'CHECKS.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)


if __name__=='__main__':main()
