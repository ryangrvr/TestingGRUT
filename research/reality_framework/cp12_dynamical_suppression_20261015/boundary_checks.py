"""Independent Python canonicalizer checks of sampled large adaptive columns."""
from collections import defaultdict
from math import lcm,sqrt
from pathlib import Path
import argparse,json,time
import numpy as np
from boundary import load_model,transitions,gauss,hopping_edge,outside_weight,flow
from invariant_evolution import Invariant

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--model',default='N6_SHELL2')
    ap.add_argument('--samples',type=int,default=40);args=ap.parse_args();start=time.monotonic()
    m=load_model(6,12,12,args.model);sub=m['sub'];sc=lcm(8,150)
    H,V,B,BV=(m[k] for k in ['H','V','B','BV']);states=m['states'];outer=m['outer_states']
    samples=set(int(j) for j in np.linspace(0,len(states)-1,args.samples))
    for j in [0,715984,715985,940984,940985,1421196,1421197,len(states)-1]:
        if j<len(states):samples.add(j)
    maximum=0.;raw_hops=0;rewiring_exits=0
    for col in sorted(samples):
        s=states[col];gauss(6,s);info=flow(6,s);expectedH=defaultdict(float);expectedV=defaultdict(float)
        sub.cache.clear()
        for t,c in transitions(6,s,sc):
            gauss(6,t);edge=hopping_edge(6,s,t);raw_hops+=edge is not None
            key,sign,size=Invariant.canonical(sub,t)
            if not sign:continue
            a=c/sc*sign*sqrt(m['sizes'][col]/size)
            expectedH[key]+=a
            if edge is not None:expectedV[key]+=outside_weight(6,s,t,edge,info)*a
        actualH={};actualV={}
        for matrix,rows,dest in [(H,states,actualH),(V,states,actualV),(B,outer,actualH),(BV,outer,actualV)]:
            z=matrix[:,col].tocoo()
            for r,a in zip(z.row,z.data):dest[rows[r]]=dest.get(rows[r],0.)+a
        for actual,expected in [(actualH,expectedH),(actualV,expectedV)]:
            maximum=max(maximum,max((abs(actual.get(s,0.)-expected.get(s,0.)) for s in set(actual)|set(expected)),default=0.))
        z=B[:,col].tocoo()
        rewiring_exits+=sum(t[1]==s[1] and t[0]!=s[0] for t in (outer[r] for r in z.row))
    assert maximum<1e-10
    result=dict(status='PASS',model=args.model,dimension=len(states),outgoing_dimension=len(outer),
        sampled_columns=len(samples),raw_Gauss_hops=raw_hops,maximum_full_column_entry_error=maximum,
        sampled_outgoing_rewiring_entries=rewiring_exits,
        comparison='Regenerate unregulated raw destinations, independently canonicalize with the original pure Python implementation, compare complete retained plus outgoing H and outside hopping columns. Samples include all extension boundaries.',
        seconds=time.monotonic()-start)
    (P/'BOUNDARY_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
