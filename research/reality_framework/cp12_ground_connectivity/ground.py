"""Full frozen CP12 Q=N ground connectivity, with explicit symmetry scope.

The invariant run is an exact compression, not a proof that the global ground
state has this symmetry. Full runs independently test that question at small N.
Hard electric cutoffs are reported and varied; they are not changes to the law.
"""
from pathlib import Path
import argparse
import json
import time
import math
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import eigsh, LinearOperator
from invariant_evolution import Invariant
from frozen_model import build

P=Path(__file__).resolve().parent

def topology(n,g,es):
    adj=[set() for _ in range(n)]
    for e,(a,b) in enumerate(es):
        if g>>e&1: adj[a].add(b);adj[b].add(a)
    components=[];remaining=set(range(n))
    while remaining:
        todo=[min(remaining)];part=set(todo);remaining-=part
        while todo:
            i=todo.pop()
            for j in adj[i]&remaining:
                remaining.remove(j);part.add(j);todo.append(j)
        components.append(part)
    ds=[len(a) for a in adj]
    bridges=0
    if len(components)==1:
        # In Q=N the unique cycle is what survives leaf pruning.
        leaf=[i for i,d in enumerate(ds) if d==1];pruned=set();dd=ds.copy()
        while leaf:
            i=leaf.pop();pruned.add(i)
            for j in adj[i]-pruned:
                dd[j]-=1
                if dd[j]==1:leaf.append(j)
        bridges=len(pruned)
    return [len(components)==1,len(components),max(map(len,components))/n,
            sum(d==0 for d in ds),sum(abs(d-2) for d in ds),all(d==2 for d in ds),
            bridges,max(ds),sum(d*(d-1)/2 for d in ds)]

def lowest(H,k=2):
    if H.shape[0]<=k+1:
        vals,vec=np.linalg.eigh(H.toarray());return vals[:k],vec[:,:k]
    rng=np.random.default_rng(1234567)
    vals,vec=eigsh(H,k=k,which='SA',tol=2e-11,maxiter=30000,
                   v0=rng.standard_normal(H.shape[0]))
    ix=np.argsort(vals);return vals[ix],vec[:,ix]

def run(n,L,full=False,blocks=True):
    start=time.monotonic();Q=n
    if full:
        H,states,_,_=build(n,Q,L)
        import itertools
        es=list(itertools.combinations(range(n),2));fd=len(states)
        sizes=np.ones(fd,int)
        scope='complete regulated physical Hilbert space'
    else:
        sub=Invariant(n,Q,L);H=sub.matrix();states=sub.states
        es=sub.es;fd=sub.full_dimension;sizes=sub.sizes
        scope='neutral-background S_N-invariant compression; global minimum not assumed'
    print(json.dumps(dict(stage='matrix',N=n,Q=Q,L=L,full=full,dimension=H.shape[0],
                          full_dimension=fd,nnz=H.nnz)),flush=True)
    vals,vec=lowest(H);g=vec[:,0];E=float(vals[0]);prob=g*g
    flags=np.array([topology(n,s[0],es) for s in states],float).T
    names=['connected','components','largest_component_fraction','isolates','degree_defect',
           'degree_two','bridges_on_connected_graph','maximum_degree','wedges']
    obs=dict(zip(names,(flags@prob).tolist()))
    boundary=np.array([max(map(abs,s[2]))==L for s in states])
    nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in states])
    fields=np.array([sum(x*x for x in s[2]) for s in states])
    graph_weights={}
    for i,s in enumerate(states):graph_weights[s[0]]=graph_weights.get(s[0],0.)+float(prob[i])
    shape_weights={}
    for i,s in enumerate(states):
        # Distribution of unlabeled topology features, not graph identities.
        f=flags[:,i]
        key=f"components={int(f[1])};isolates={int(f[3])};D2={int(f[4])};maxdegree={int(f[7])}"
        shape_weights[key]=shape_weights.get(key,0.)+float(prob[i])
    records=dict(status='CP12_QN_GROUND_CONNECTIVITY_PASS',N=n,Q=Q,L=L,full=full,
        dimension=H.shape[0],full_physical_dimension=fd,scope=scope,
        energies=vals.tolist(),gap=float(vals[1]-vals[0]),energy_per_vertex=E/n,
        eigen_residuals=[float(np.linalg.norm(H@vec[:,i]-vals[i]*vec[:,i])) for i in range(len(vals))],
        observables=obs,cutoff_boundary_probability=float(boundary@prob),
        charge_density=float(nv@prob/n),electric_square_per_vertex=float(fields@prob/n),
        topology_distribution=shape_weights,
        graph_class_probability=[dict(graph_mask=int(mask),probability=p)
                                 for mask,p in sorted(graph_weights.items())])
    # Compressing a full matrix must keep the boundary diagonal; rebuild on
    # connected masks would instead delete outgoing rewiring degree, a new law.
    connected=flags[0].astype(bool)
    if blocks and connected.any() and (~connected).any():
        ic=np.flatnonzero(connected);idc=np.flatnonzero(~connected)
        HC=H[ic,:][:,ic];HD=H[idc,:][:,idc];B=H[ic,:][:,idc]
        ec,vc=lowest(HC,k=1);ed,vd=lowest(HD,k=1)
        cross=LinearOperator((B.shape[0],)*2,matvec=lambda v:-(B@(B.T@v)),dtype=float)
        bs,bv=lowest(cross,k=1);b=float(math.sqrt(max(0,-bs[0])))
        pc=obs['connected'];d=float(ec[0]-E)
        rhs=B@g[idc];lhs=(HC-E*eye(HC.shape[0]))@g[ic]
        records['sector_comparison']=dict(connected_energy=float(ec[0]),disconnected_energy=float(ed[0]),
            connected_energy_minus_ground=d,cross_norm=b,
            resolvent_probability_upper=b*b/(d*d+b*b) if d>0 else 1.,
            projected_eigenequation_residual=float(np.linalg.norm(lhs+rhs)),
            conditional_connected_energy=float(np.vdot(g[ic],HC@g[ic]).real/pc),
            conditional_disconnected_energy=float(np.vdot(g[idc],HD@g[idc]).real/(1-pc)),
            graph_mixing_energy=float(2*np.vdot(g[ic],B@g[idc]).real),
            energy_decomposition_error=float(abs(E-np.vdot(g[ic],HC@g[ic]).real-
                np.vdot(g[idc],HD@g[idc]).real-2*np.vdot(g[ic],B@g[idc]).real)))
        assert pc <= records['sector_comparison']['resolvent_probability_upper']+1e-9
        assert records['sector_comparison']['projected_eigenequation_residual']<1e-8
    assert max(records['eigen_residuals'])<1e-8
    records['seconds']=time.monotonic()-start
    output=P/f"ground_N{n}_L{L}_{'full' if full else 'invariant'}.json"
    output.write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps({k:records[k] for k in ('N','L','dimension','energies','observables','cutoff_boundary_probability','seconds')}),flush=True)
    if 'sector_comparison' in records:print(json.dumps(records['sector_comparison']),flush=True)
    return records

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--N',type=int,required=True)
    a.add_argument('--L',type=int,required=True);a.add_argument('--full',action='store_true')
    a.add_argument('--skip-blocks',action='store_true');a=a.parse_args()
    run(a.N,a.L,a.full,not a.skip_blocks)
