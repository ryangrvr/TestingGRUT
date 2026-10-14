"""Independent small-graph oracles and physical orbit checks."""
from fractions import Fraction as F
from itertools import combinations,permutations
from math import comb
import json,time
import numpy as np
from analyze import P,catalog,features,distances,core_mask,model,resources,setup
from exact_destination import transitions,gauss


def oracle(n,g):
    es=list(combinations(range(n),2));A=np.zeros((n,n),dtype=np.int64)
    for e,(a,b) in enumerate(es):
        if g>>e&1:A[a,b]=A[b,a]=1
    powers=[np.eye(n,dtype=np.int64)]
    for _ in range(n):powers.append(powers[-1]@A)
    D=np.full((n,n),n,dtype=int)
    for r in range(n-1,-1,-1):D[powers[r]>0]=r
    reach=D<n;C={tuple(row.tolist()) for row in reach}
    rank=g.bit_count()-n+len(C)
    oncycle=set();short=[]
    for length in range(3,n+1):
        for verts in combinations(range(n),length):
            for tail in permutations(verts[1:]):
                order=(verts[0],)+tail
                if order[1]>order[-1]:continue
                if all(A[order[j],order[(j+1)%length]] for j in range(length)):
                    oncycle.update(order)
                    if length<=4:short.append(set(order))
    pack=0
    for subset in combinations(short,2):
        if not subset[0]&subset[1]:pack=2
    pack=max(pack,int(bool(short)))
    deg=A.sum(axis=0)
    triangle=int(np.trace(powers[3]))//6
    square=(int(np.trace(powers[4]))-2*g.bit_count()-4*sum(comb(int(d),2) for d in deg))//8
    vals=[float(len(C)==1),float(rank>0),max(reach.sum(axis=1))/n,
        sum(deg>0)/n,len(oncycle)/n,triangle/comb(n,3),square/(3*comb(n,4)),
        pack/(n//3),rank/(len(es)-n+1),sum(deg**2)/(n*(n-1)**2),
        (reach.sum()-n)/(n*(n-1))]
    vals += [np.mean((D<=r).sum(axis=1))/n for r in range(n)]
    return np.array(vals),D


def run():
    start=time.monotonic();graphs_checked=0;maxerr=0.;catalogs=[]
    for n,expected in [(4,11),(5,34),(6,156)]:
        C=catalog(n);M=C['M'];assert len(C['keys'])==expected
        masks=list(range(1<<M)) if n<=5 else C['keys']
        for g in masks:
            f,D=oracle(n,g)
            maxerr=max(maxerr,float(np.max(abs(f-features(n,g)))))
            assert np.array_equal(D,distances(n,g))
            graphs_checked+=1
        for k in range(M+1):
            tri=float(F(comb(k,3),comb(M,3))) if k>=3 else 0.
            sq=float(F(comb(k,4),comb(M,4))) if k>=4 else 0.
            assert abs(C['uniform_scores'][k,5]-tri)<1e-12
            assert abs(C['uniform_scores'][k,6]-sq)<1e-12
            assert abs(C['uniform_scores'][k,11+1]-(1+2*k/n)/n)<1e-12
        catalogs.append(dict(N=n,labelled_graphs=1<<M,unlabelled_graphs=len(C['keys'])))
    assert maxerr<1e-12
    # Independent raw-permutation unfolding at N4, including fermionic signs.
    n=q=4;sub,H,B,states,sizes,psi=model(n,q,256)
    r,k,nv,e2=resources(n,q,sub,states);S=setup(n,q,states,sizes,r,nv,e2)
    from scipy.sparse.linalg import expm_multiply
    v=expm_multiply(-10j*H,psi,traceA=-10j*float(H.diagonal().sum()));p=abs(v)**2
    edge=np.zeros(sub.M);prot=edge.copy();matter=np.zeros(n)
    scores=np.zeros(len(S['names']));rewire_targets=0;expanded=0
    for j,s in enumerate(states):
        gauss(n,s);seen={}
        for perm in sub.group:
            target,sign=sub.transform(s,perm)
            if target in seen:assert seen[target]==sign
            seen[target]=sign
        assert len(seen)==sizes[j]
        for target in seen:
            g,bits,fields=target;core=core_mask(n,target);w=p[j]/sizes[j]
            edge+=w*np.array([g>>e&1 for e in range(sub.M)])
            prot+=w*np.array([core>>e&1 for e in range(sub.M)])
            matter+=w*np.array([bits>>a&1 for a in range(n)])
            scores+=w*features(n,core);expanded+=1
        for target,c in transitions(n,s,120,hopping=False):
            assert core_mask(n,s)==core_mask(n,target)
            rewire_targets+=1
    edgeerr=float(np.max(abs(edge-q/sub.M)))
    proterr=float(np.max(abs(prot-(q-p@r)/sub.M)))
    mattererr=float(np.max(abs(matter-(p@nv)/n)))
    foldederr=float(np.max(abs(scores-p@S['actual'])))
    assert max(edgeerr,proterr,mattererr,foldederr)<1e-11
    # Exact rational energy gap on the full allowed resource range.
    triples=0
    for n in range(3,10):
        M=comb(n,2)
        for q in range(1,M):
            a=F(M-q,10*M)
            for nv0 in range(min(n,q)+1):
                for r0 in range(q-nv0+1):
                    k0=q-r0-nv0;e2min=2*k0+nv0
                    gap=2*nv0+F(e2min,8)+a*r0
                    assert gap>=a*q
                    triples+=1
            assert a*q>0
    out=dict(status='CP12_RELATIONAL_INDEPENDENT_CHECKS_PASS',
        graph_oracles_checked=graphs_checked,maximum_graph_feature_error=maxerr,
        exact_graph_catalogs=catalogs,expanded_N4_physical_states=expanded,
        raw_N4_nonhopping_transitions=rewire_targets,
        N4_labelled_edge_marginal_error=edgeerr,N4_protected_marginal_error=proterr,
        N4_matter_marginal_error=mattererr,N4_folded_observable_error=foldederr,
        rational_energy_gap_resource_cases=triples,seconds=time.monotonic()-start)
    (P/'INDEPENDENT_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':run()
