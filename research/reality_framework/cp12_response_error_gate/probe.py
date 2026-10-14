"""Same CP12 unmasked preparation, observable-specific error gate."""
from pathlib import Path
from itertools import combinations
import argparse,hashlib,json,sys,time
import numpy as np
from control import response
from run import model,graphstats

P=Path(__file__).resolve().parent

def detectors(n,states,es):
    distances={}
    for g,_,_ in states:
        if g in distances:continue
        adj=[set() for _ in range(n)]
        for e,(a,b) in enumerate(es):
            if g>>e&1:adj[a].add(b);adj[b].add(a)
        dist=[n]*n;dist[0]=0;todo=[0]
        for a in todo:
            for b in adj[a]:
                if dist[b]==n:dist[b]=dist[a]+1;todo.append(b)
        distances[g]=dist
    return [np.array([sum((bits>>j)&1 for j in range(1,n) if distances[g][j]==r)/(n-1) for g,bits,E in states]) for r in [1,2]]


def run(n,budget,steps):
    start=time.monotonic();q=2*n;sub,H,B,states,sizes,psi=model(n,q,budget,True)
    source=np.array([bits&1 for g,bits,E in states],float)
    O=detectors(n,states,sub.es)
    total=np.array([sum((bits>>j)&1 for j in range(1,n))/(n-1) for g,bits,E in states])
    assert np.max(abs(O[0]+O[1]-total))<1e-14  # Q=2N, N=5 or 6.
    O.append(total)
    rows,Z,v=response(H,B,psi,source,O,.1,.1,steps,q/2)
    norm_error=max(abs(float(np.linalg.norm(v)**2)-1),abs(float(np.linalg.norm(Z[:,0])**2)-1),abs(float(np.linalg.norm(Z[:,1])**2)-float(np.linalg.norm(source*v)**2)))
    energy_error=abs(float(np.vdot(Z[:,0],H@Z[:,0]).real)-q/2)
    assert norm_error<1e-10 and energy_error<1e-9
    for j,row in enumerate(rows):
        row['channel']=['final_distance_1','final_distance_2','mean_other_vertex_density'][j]
        row['final_distance']=j+1 if j<2 else None
        row['response_interval']=[row['projected_response']-row['observable_absolute_error_upper'],row['projected_response']+row['observable_absolute_error_upper']]
        row['sign_resolved']=abs(row['projected_response'])>row['observable_absolute_error_upper']
        print(json.dumps(row),flush=True)
    gs={g:graphstats(n,g,sub.es) for g in {s[0] for s in states}}
    G=np.array([gs[g] for g,_,_ in states],float).T
    names=['connected','triangles','squares','short_cycle_packing','degree_square_sum']+[f'ball_{r}' for r in range(n)]+[f'shell_{r}' for r in range(n)]
    assert abs(rows[0]['projected_response']+rows[1]['projected_response']-rows[2]['projected_response'])<1e-12
    out=dict(N=n,Q=q,electric_square_budget=budget,wait=.1,lag=.1,integration_steps=steps,
             retained_quotient_dimension=H.shape[0],outgoing_quotient_orbits=B.shape[0],
             responses=rows,geometric_observables_at_wait=dict(zip(names,map(float,G@abs(v)**2))),
             geometric_observables_at_wait_plus_lag=dict(zip(names,map(float,G@abs(Z[:,0])**2))),
             initial_energy=q/2,initial_variance=2*q,initial_min_entropy_bits=0,
             numerical_norm_error=norm_error,numerical_baseline_energy_error=energy_error,
             floating_point_bounds=True,scope='Final graph-distance conditional density response. A local source phase kick; includes graph redistribution. No fitted dimension, photon identification, or autonomous state selection.',seconds=time.monotonic()-start)
    (P/f'N{n}_B{budget}_S{steps}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--budget',type=int,required=True);ap.add_argument('--steps',type=int,default=64);a=ap.parse_args();run(a.n,a.budget,a.steps)
