"""Electric resources needed to represent spanning, cycle-rich circulation."""
from functools import lru_cache
from itertools import combinations
import argparse,json,time
import numpy as np
from audit import P,REL,model,positive_interval
from loops import incidence_masks,directed_support
from frozen_model import build


@lru_cache(None)
def support(n,g):
    vertices=set();parent=list(range(n));edges=0
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    for e,(a,b) in enumerate(combinations(range(n),2)):
        if g>>e&1:
            vertices.update([a,b]);edges+=1;aa,bb=root(a),root(b);parent[aa]=bb
    components=len({root(v) for v in vertices})
    rank=edges-len(vertices)+components
    return edges,len(vertices),components,rank


def check(n,states):
    margins=[];spanning=0;maxrank=0;maxspan=0;events=[]
    for s in states:
        through,sign=incidence_masks(n,s);C=directed_support(n,through,sign)
        edges,vertices,components,rank=support(n,C)
        nv=(s[1]&((1<<n)-1)).bit_count()
        hc=(C&~(s[1]>>n)).bit_count()
        e2=sum(x*x for x in s[2]);margin=e2-(2*edges+nv+2*hc)
        assert margin>=0
        spans=vertices==n and components==1
        if spans:
            assert rank<=e2/2-n+1-nv/2-hc
            spanning+=1;maxspan=max(maxspan,rank)
        margins.append(margin);maxrank=max(maxrank,rank)
        events.append((bool(C),spans,spans and rank>=2))
    return dict(configurations_checked=len(states),minimum_inequality_margin=min(margins),
        spanning_circulation_configurations=spanning,maximum_cycle_rank=maxrank,
        maximum_spanning_cycle_rank=maxspan),np.array(events,float)


def run(snapshot_dir):
    start=time.monotonic();HF,physical,_,_=build(4,4,3)
    independent,_=check(4,physical);rows=[]
    for n,q,budget,t in [(4,4,256,10),(5,5,128,10),(6,6,128,10),(6,12,10,1)]:
        sub,H,B,states,sizes,psi=model(n,q,budget)
        controls,events=check(n,states)
        old=json.loads((REL/f'RELATIONAL_N{n}_Q{q}_B{budget}.json').read_text())
        e=min(old['rows'],key=lambda r:abs(r['T']-t))['state_error_upper']
        v=np.load(str(snapshot_dir)+f'/N{n}_Q{q}_B{budget}_T{t:g}.npz')['wave']
        values=(abs(v)**2)@events
        names=['circulation_any','circulation_spanning','circulation_spanning_cycle_rank_at_least_2']
        rows.append(dict(N=n,Q=q,budget=budget,T=t,controls=controls,
            probabilities={name:dict(projected=float(value),full_law_interval=positive_interval(float(value),e))
                           for name,value in zip(names,values)},
            input_state_error_upper=e,
            necessary_neutral_spanning_circulation_budget_lower_bound=2*n,
            necessary_neutral_spanning_rank_2_budget_lower_bound=2*n+2,
            necessary_neutral_full_degree_4_budget_lower_bound=4*n))
    out=dict(status='CP12_CIRCULATION_CAPACITY_PASS',independent_physical_check=independent,cases=rows,
        seconds=time.monotonic()-start,
        exact_inequality='E_squared >= 2*e_C + n_v + 2*h_C, where C is the union of all directed circulation cycles and h_C counts its hole edges. For a connected spanning C, beta_C <= E_squared/2 - N + 1 - n_v/2 - h_C.',
        scope='Necessary resource and numerical-representation conditions, not a proof that degree four is sufficient or necessary for three dimensions. A degree-four 3D network is an example with extensive cycle rank. Q=N spanning circulation still has only one independent cycle. Numerical projected zeros outside the budget are not full-law absences.')
    (P/'CYCLE_CAPACITY.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],cases=[dict(N=r['N'],Q=r['Q'],checked=r['controls']['configurations_checked']) for r in rows],seconds=out['seconds'])),flush=True)
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
