"""Unfiltered original-neutral circulation connectivity and intrinsic balls."""
import argparse,json,time
import numpy as np
from shared import *


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--wait',type=float,default=5);a=ap.parse_args();start=time.monotonic()
    z=np.load(cache_root()/f'N{a.n}_Q{a.q}_B{a.budget}_W{a.wait:g}.npz')
    cores=z['core'];weights=abs(z['w'])**2;e=float(z['state_error']);n=a.n
    masks,inv=np.unique(cores,return_inverse=True);probabilities=np.bincount(inv,weights=weights)
    any_p=0.;span_p=0.;cover_p=0.;cycle_p=0.;balls=np.zeros(n);by_edges={}
    for C,p in zip(masks,probabilities):
        C=int(C);ds=[distances(n,C,i) for i in range(n)]
        connected=all(v is not None for v in ds[0])
        covers=all(sum(v==1 for v in d)>0 for d in ds)
        if C:any_p+=p
        if covers:cover_p+=p
        if is_spanning_cycle(n,C):cycle_p+=p
        if connected:
            span_p+=p
            balls+=p*np.mean([[sum(v is not None and v<=r for v in d) for r in range(n)] for d in ds],axis=0)
        count=C.bit_count();by_edges[count]=by_edges.get(count,0.)+p
    pI=probability_interval(float(span_p),e)
    ball_rows=[]
    for r,moment in enumerate(balls):
        lo=max(0.,np.sqrt(moment)-np.sqrt(n)*e)**2
        hi=min(float(n),(np.sqrt(moment)+np.sqrt(n)*e)**2)
        ball_rows.append(dict(radius=r,conditional_mean=float(moment/span_p) if span_p else None,
            conditional_mean_interval=[lo/pI[1],min(float(n),hi/pI[0])] if pI[0]>0 else [0.,float(n)]))
    data=dict(status='PASS',N=n,Q=a.q,budget=a.budget,wait=a.wait,state_error=e,
        any_circulation_probability=float(any_p),any_circulation_probability_interval=probability_interval(float(any_p),e),
        all_vertices_covered_probability=float(cover_p),all_vertices_covered_interval=probability_interval(float(cover_p),e),
        connected_spanning_circulation_probability=float(span_p),connected_spanning_circulation_interval=pI,
        exact_spanning_cycle_probability=float(cycle_p),
        circulation_edge_count_distribution={str(k):dict(probability=float(v),interval=probability_interval(float(v),e)) for k,v in sorted(by_edges.items())},
        all_root_mean_balls_given_spanning_connected_circulation=ball_rows,
        interpretation='All original-neutral amplitudes are retained. Connectivity and balls use the entire directed-cycle edge mask, without selecting rings. Conditional ball averages explicitly condition on spanning connected circulation. These few small finite graphs do not define a continuum dimension; coordination is not converted to dimension.',
        seconds=time.monotonic()-start)
    (P/f'GEOMETRY_N{n}_Q{a.q}_B{a.budget}_W{a.wait:g}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data),flush=True)


if __name__=='__main__':main()
