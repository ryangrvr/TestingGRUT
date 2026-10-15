"""Unregulated full-H actions on capacity-valid charged circulation controls.

Every graph is explicitly supplied. These controls test energetic availability
and failure of energy-only dimension selection, never autonomous formation.
"""
from collections import defaultdict
from itertools import combinations
from math import sqrt,lcm
from pathlib import Path
import json,sys,time
import numpy as np
from bounds import P,dressed_moments

FROZEN=P.parent/'cp12_ground_connectivity'
sys.path.insert(0,str(FROZEN))
from exact_destination import transitions,hop,gauss
from frozen_model import hop as independent_hop,gauss as independent_gauss


def euler(n,es):
    adjacency=[[] for _ in range(n)]
    for e,(a,b) in enumerate(es):adjacency[a].append((b,e));adjacency[b].append((a,e))
    assert all(len(x)==4 for x in adjacency)
    used=set();stack=[0];edge_stack=[];tour=[]
    while stack:
        v=stack[-1]
        while adjacency[v] and adjacency[v][-1][1] in used:adjacency[v].pop()
        if adjacency[v]:
            w,e=adjacency[v].pop();used.add(e);stack.append(w);edge_stack.append(e)
        else:
            stack.pop()
            if edge_stack:tour.append((edge_stack.pop(),stack[-1],v))
    assert len(used)==len(es)
    positive={e:a for e,a,b in tour}
    assert all(sum(positive[e]==v for e,pair in enumerate(es) if v in pair)==2 for v in range(n))
    return positive


def graph(n,kind):
    if kind=='complete':return list(combinations(range(n),2))
    if kind=='octahedron':return [(a,b) for a,b in combinations(range(n),2) if a//2!=b//2]
    if kind=='bipartite':return [(a,b) for a in range(n//2) for b in range(n//2,n)]
    if kind=='ring_power':return sorted({tuple(sorted((a,(a+s)%n))) for a in range(n) for s in [1,2]})
    raise ValueError(kind)


def state(n,es,eta,phase=0.,reverse=False):
    all_edges=list(combinations(range(n),2));index={x:i for i,x in enumerate(all_edges)}
    positive=euler(n,es)
    if reverse:positive={e:(b if positive[e]==a else a) for e,(a,b) in enumerate(es)}
    g=sum(1<<index[e] for e in es);M=len(all_edges)
    bits=sum(1<<(n+e) for e in range(M));E=[0]*(2*M)
    assigned=[[] for _ in range(n)]
    for local,(a,b) in enumerate(es):
        e=index[a,b];p=positive[local];leg=0 if p==a else 1
        E[2*e+leg]=1;E[2*e+1-leg]=-1;assigned[p].append((e,leg,1 if leg==0 else -1))
    wave={(g,bits,tuple(E)):1.};c=sqrt(1-eta);s=sqrt(eta/2);phase_factor=np.exp(1j*phase)
    for v in range(n):
        nxt=defaultdict(complex)
        for (g,b,E),amp in wave.items():
            nxt[g,b,E]+=c*amp
            for e,leg,orientation in assigned[v]:
                nb,phase=hop(b,v,n+e);F=list(E);F[2*e+leg]+=1
                nxt[g,nb,tuple(F)]+=-s*phase_factor*orientation*phase*amp
        wave=dict(nxt)
    for key in wave:gauss(n,key);assert not np.any(independent_gauss(n,key))
    assert abs(sum(abs(a)**2 for a in wave.values())-1)<1e-11
    return wave


def action(n,wave,independent=False):
    all_edges=list(combinations(range(n),2));M=len(all_edges);out=defaultdict(float)
    scale=lcm(8,10*M)
    for key,amp in wave.items():
        g,b,E=key
        if not independent:
            for target,value in transitions(n,key,scale):out[target]+=value/scale*amp
        else:
            q=g.bit_count();nv=(b&((1<<n)-1)).bit_count()
            out[key]+=(q/2+2*nv+sum(x*x for x in E)/8)*amp
            for e,(i,j) in enumerate(all_edges):
                if not(g>>e&1):continue
                assert E[2*e]*E[2*e+1]<0
                for leg,v,orientation in [(0,i,1),(1,j,-1)]:
                    if b>>(n+e)&1 and not(b>>v&1):dst,src,shift=v,n+e,1
                    elif b>>v&1 and not(b>>(n+e)&1):dst,src,shift=n+e,v,-1
                    else:continue
                    nb,phase=independent_hop(b,dst,src);F=list(E);F[2*e+leg]+=shift
                    out[g,nb,tuple(F)]+=orientation*phase*amp
    return {s:a for s,a in out.items() if a!=0}


def run():
    start=time.monotonic();a=19/8
    eta0=(60-32*sqrt(2))/291;etag=(1-a/sqrt(a*a+8))/2
    cases=[(5,'complete',eta0),(5,'complete',etag),(5,'complete',.3),
           (6,'octahedron',eta0),(8,'ring_power',eta0),(8,'bipartite',eta0)]
    rows=[];worst=0.
    for n,kind,eta in cases:
        es=graph(n,kind);assert len(es)==2*n
        wave=state(n,es,eta);hv=action(n,wave);other=action(n,wave,True)
        difference=sqrt(sum(abs(hv.get(s,0)-other.get(s,0))**2 for s in hv.keys()|other.keys()))
        assert difference<1e-10
        mean=float(sum(np.conj(amp)*hv.get(s,0) for s,amp in wave.items()).real)
        variance=sum(abs(hv.get(s,0)-mean*wave.get(s,0))**2 for s in hv.keys()|wave.keys())
        loss=0.;protected_loss=0.
        for (g,b,E),amp in hv.items():
            gauss(n,(g,b,E))
            present=[e for e in range(len(E)//2) if g>>e&1]
            if not all(E[2*e]*E[2*e+1]<0 for e in present):loss+=abs(amp)**2
            if not all(E[2*e] or E[2*e+1] for e in present):protected_loss+=abs(amp)**2
        prediction=dressed_moments(eta)
        errors=[abs(mean/n-prediction['mean_energy_per_vertex']),abs(variance/n-prediction['variance_per_vertex']),
                abs(loss/n-prediction['full_circulation_loss_T2_per_vertex']),protected_loss]
        worst=max(worst,*errors);assert max(errors)<1e-10
        rows.append(dict(N=n,Q=len(es),supplied_graph=kind,matter_density=eta,
            support_configurations=len(wave),one_full_H_action_configurations=len(hv),
            full_mean_energy=mean,full_energy_variance=variance,
            full_circulation_loss_T2_coefficient=loss,full_protection_loss_T2_coefficient=protected_loss,
            electric_square_maximum_in_input=max(sum(x*x for x in s[2]) for s in wave),
            electric_square_maximum_after_full_H_action=max(sum(x*x for x in s[2]) for s in hv),
            independent_H_action_norm_difference=difference,maximum_formula_difference=max(errors),
            prediction_per_vertex=prediction))
        print(json.dumps(dict(stage='charged_comparator',N=n,graph=kind,eta=eta,energy=mean,variance=variance,error=max(errors))),flush=True)
    n=5;eta=.25;h=.05;es=graph(n,'complete');bright=sqrt(3/2)
    wa=state(n,es,eta,float(np.arccos((35/32-h)/bright)))
    wb=state(n,es,eta,float(np.arccos((35/32+h)/bright)),reverse=True)
    ha=action(n,wa);hb=action(n,wb)
    assert wa.keys().isdisjoint(wb) and ha.keys().isdisjoint(hb)
    pair={s:amp/sqrt(2) for w in [wa,wb] for s,amp in w.items()}
    hp=action(n,pair);mean=float(sum(np.conj(a)*hp.get(s,0) for s,a in pair.items()).real)
    variance=sum(abs(hp.get(s,0)-mean*pair.get(s,0))**2 for s in hp.keys()|pair.keys())
    predicted=n*(31/16+(n-1)*h*h)
    assert abs(mean-n)<1e-10 and abs(variance-predicted)<1e-10
    pair_control=dict(N=n,eta=eta,energy_offset_per_vertex=h,full_mean_energy=mean,
        full_variance=variance,predicted_full_variance=predicted,
        opposite_orientation_supports_orthogonal=True,one_action_ranges_orthogonal=True,
        scope='Small check of the general two-orientation moment formula. Exact mean-and-variance matching to the neutral ray uses h^2=33/[16(N-1)] and is feasible for N>=128, not claimed for this N5 check.')
    result=dict(status='CP12_CHARGED_CIRCULATION_CONTROLS_PASS',cases=rows,small_pair_control=pair_control,
        maximum_formula_difference=worst,seconds=time.monotonic()-start,
        scope='Supplied four-regular Euler-oriented hosts only. All input configurations are fully circulating, r=0, with matter and finite-support fields. Both Hamiltonian actions include every unbounded-field outgoing transition, without an electric cutoff. No physical full-time evolution, autonomous geometry selection, or propagation cone is inferred.')
    (P/'CHARGED_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':run()
