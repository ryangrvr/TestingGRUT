"""Topology-neutral CP12 formation test with explicit electric-budget leakage.

The budget is a numerical projection, never a change to the frozen law.
Complete closure inside the budget; all outgoing full-law transitions retained
in a separate residual matrix. Signed permutation quotient is exact.
"""
from pathlib import Path
from collections import deque
from itertools import combinations, permutations
from math import comb, sqrt, lcm
import argparse, hashlib, json, sys, time, os, pickle
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
S=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((S/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
sys.path.insert(0,str(S))
from invariant_evolution import Invariant
from exact_destination import transitions, gauss


def model(n,q,budget,source=False):
    cache=Path(os.environ['CP12_FORMATION_CACHE']) if os.environ.get('CP12_FORMATION_CACHE') else None
    cache_name=f'n{n}_q{q}_b{budget}_s{int(source)}_v1.pkl'
    if cache and (cache/cache_name).exists():
        with (cache/cache_name).open('rb') as f:return pickle.load(f)
    sub=Invariant.__new__(Invariant)
    sub.n=n;sub.Q=q;sub.es=list(combinations(range(n),2));sub.M=len(sub.es)
    group=( [(0,)+p for p in permutations(range(1,n))] if source else permutations(range(n)) )
    sub.group=[sub.permutation(p) for p in group]
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)}
    sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for a,b in enumerate(p['vertex']): inv[b]=a
        sub.inverse.append(gid[tuple(inv)])
    sub.gcache={};sub.cache={}
    masks=[sum(1<<e for e in c) for c in combinations(range(sub.M),q)]
    reps=sorted({sub.graph_info(g)[0] for g in masks})
    sea=sum(1<<(n+e) for e in range(sub.M));zero=(0,)*(2*sub.M)
    states=[];index={};sizes=[];queue=deque()
    def add(s):
        key,sign,size=sub.canonical(s)
        if sign and key not in index:
            index[key]=len(states);states.append(key);sizes.append(size);queue.append(key)
        return key,sign,size
    for g in reps: add((g,sea,zero))
    sc=lcm(8,10*sub.M);columns=[];outer={};br=[];bc=[];bv=[]
    while queue:
        s=queue.popleft();col=index[s];gauss(n,s);entries=[]
        for target,c in transitions(n,s,sc):
            inside=sum(x*x for x in target[2])<=budget
            key,sign,size=add(target) if inside else sub.canonical(target)
            if not sign: continue
            value=c/sc*sign*sqrt(sizes[col]/size)
            if inside: entries.append((index[key],value))
            else:
                if key not in outer:outer[key]=len(outer)
                br.append(outer[key]);bc.append(col);bv.append(value)
        columns.append(entries)
        if col and col%10000==0:print(json.dumps(dict(stage='closure',N=n,budget=budget,processed=col,queued=len(queue))),flush=True)
    rr=[];cc=[];vv=[]
    for j,col in enumerate(columns):
        for i,v in col:rr.append(i);cc.append(j);vv.append(v)
    d=len(states)
    H=coo_matrix((vv,(rr,cc)),shape=(d,d)).tocsr()
    B=coo_matrix((bv,(br,bc)),shape=(len(outer),d)).tocsr()
    assert np.linalg.norm((H-H.T).data)<1e-9
    psi=np.array([sqrt(sizes[i]/comb(sub.M,q)) if bits==sea and not any(E) else 0 for i,(g,bits,E) in enumerate(states)],complex)
    assert abs(np.linalg.norm(psi)-1)<1e-12
    assert abs(np.vdot(psi,H@psi).real-q/2)<1e-10
    assert abs(np.linalg.norm(H@psi)**2-(q*q/4+2*q))<1e-9
    print(json.dumps(dict(stage='matrix',N=n,Q=q,budget=budget,dimension=d,outside=len(outer),graphs=len(reps),nnz=H.nnz)),flush=True)
    result=(sub,H,B,states,np.array(sizes),psi)
    if cache:
        cache.mkdir(parents=True,exist_ok=True)
        sub.cache.clear()
        with (cache/cache_name).open('wb') as f:pickle.dump(result,f)
    return result


def error_evolve(H,B,V,end,steps,shift):
    """Analytic Duhamel upper bound, evaluated in floating point.

    On each midpoint interval, Taylor-expand the leakage vector to first
    order and bound its remainder by ||B|| ||(H-shift)^2 V|| u^2/2.
    Every input column is treated independently and need not be normalized.
    """
    from scipy.sparse import eye
    A=H-shift*eye(H.shape[0],format='csr');dt=end/steps
    matrix=V.ndim==2
    if not matrix:V=V[:,None]
    bnorm=sqrt(float(np.max(np.asarray(abs(B).sum(axis=1))))*float(np.max(np.asarray(abs(B).sum(axis=0))))) if B.shape[0] else 0.
    h2=np.linalg.norm(A@(A@V),axis=0);bound=np.zeros(V.shape[1]);v=V.copy()
    trace=float(A.diagonal().sum())
    for _ in range(steps):
        mid=expm_multiply(-.5j*dt*A,v,traceA=-.5j*dt*trace)
        bound+=dt*np.linalg.norm(B@mid,axis=0)+dt*dt/4*np.linalg.norm(B@(A@mid),axis=0)+dt**3/24*bnorm*h2
        v=expm_multiply(-.5j*dt*A,mid,traceA=-.5j*dt*trace)
    return (v if matrix else v[:,0]),(bound if matrix else float(bound[0]))


def graphstats(n,g,es):
    adj=[set() for _ in range(n)]
    for e,(a,b) in enumerate(es):
        if g>>e&1:adj[a].add(b);adj[b].add(a)
    cycles=[]
    for r in [3,4]:
        for vertices in combinations(range(n),r):
            a=min(vertices)
            for rest in permutations([v for v in vertices if v!=a]):
                order=(a,)+rest
                if order[1]>order[-1]:continue
                if all(order[(i+1)%r] in adj[order[i]] for i in range(r)):
                    cycles.append((r,sum(1<<v for v in order)))
    def packing(available,start=0):
        return max([0]+[1+packing(available^m,i+1) for i,(_,m) in enumerate(cycles[start:],start) if m&available==m])
    vols=np.zeros(n);shells=np.zeros(n);connected=True
    for a in range(n):
        seen={a};front={a};vols[0]+=1
        for r in range(1,n):
            nxt=set().union(*(adj[v] for v in front))-seen if front else set()
            seen|=nxt;front=nxt;vols[r]+=len(seen);shells[r]+=len(nxt)
        connected=connected and len(seen)==n
    # No fitted dimension from these tiny diameters.
    ds=[len(x) for x in adj]
    return [float(connected),sum(r==3 for r,_ in cycles),sum(r==4 for r,_ in cycles),packing((1<<n)-1),sum(d*d for d in ds)]+list(vols/n)+list(shells/n)


def run(n,q,budget,end,steps):
    start=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget)
    stats={g:graphstats(n,g,sub.es) for g,_,_ in states}
    names=['connected','triangles','squares','short_cycle_packing','degree_square_sum']+[f'ball_{r}' for r in range(n)]+[f'shell_{r}' for r in range(n)]
    O=np.array([stats[g] for g,_,_ in states],float).T
    E2=np.array([sum(x*x for x in E) for g,b,E in states])
    nv=np.array([(b&((1<<n)-1)).bit_count() for g,b,E in states])
    names+=['electric_square','vertex_fermions'];O=np.vstack([O,E2,nv])
    # Residual-integral upper bound via midpoint Taylor inequality.
    # Shift by initial mean energy; this leaves norms/dynamics unchanged.
    from scipy.sparse import eye
    A=H-q/2*eye(H.shape[0],format='csr')
    bnorm=sqrt(float(np.max(np.asarray(abs(B).sum(axis=1))))*float(np.max(np.asarray(abs(B).sum(axis=0))))) if B.shape[0] else 0.
    h2norm=float(np.linalg.norm(A@(A@psi)))
    dt=end/steps;residual_bound=0.;rows=[];v=psi.copy()
    for j in range(steps+1):
        if j in {0,steps//4,steps//2,steps}:
            p=abs(v)**2
            rows.append(dict(T=j*dt,observables=dict(zip(names,map(float,O@p))),norm_error=float(abs(p.sum()-1)),energy_error=float(abs(np.vdot(v,H@v).real-q/2)),state_error_upper=residual_bound,graph_total_variation_upper=min(1.,residual_bound)))
            print(json.dumps(dict(stage='time',N=n,budget=budget,**rows[-1])),flush=True)
        if j==steps:break
        mid=expm_multiply(-.5j*dt*A,v,traceA=-.5j*dt*float(A.diagonal().sum()))
        f=float(np.linalg.norm(B@mid));f1=float(np.linalg.norm(B@(A@mid)))
        residual_bound+=dt*f+dt*dt/4*f1+dt**3/24*bnorm*h2norm
        v=expm_multiply(-.5j*dt*A,mid,traceA=-.5j*dt*float(A.diagonal().sum()))
    out=dict(N=n,Q=q,budget=budget,dimension=len(states),retained_full_orbit_dimension=int(sum(sizes)),outside_orbits=B.shape[0],initial_energy=q/2,initial_variance=2*q,initial_min_entropy_bits=0,neutral_reference_deficit_bits=float(np.log2(float(comb(sub.M,q)))),graph_representatives=len(stats),rows=rows,boundary_operator_norm_upper=bnorm,second_shifted_moment_vector_norm=h2norm,integration_steps=steps,seconds=time.monotonic()-start,scope='Complete reachable electric-budget projection in the signed invariant sector, with full-law outgoing residual; floating-point evaluation of analytic error bound. No inferred spatial dimension.')
    # Sequential projective loop-class persistence in the same preparation.
    C=(O[3]>=n//3).astype(float)
    w=min(.25,end/4);lag=min(.25,end/4)
    v,ew=error_evolve(H,B,psi,w,32,q/2)
    projected=C*v;pc=float(np.linalg.norm(projected)**2)
    z,ez=error_evolve(H,B,projected,lag,32,q/2)
    joint=float(np.linalg.norm(C*z)**2)
    out['loop_persistence']=dict(wait=w,lag=lag,packing_threshold=n//3,initial_event_probability=pc,sequential_joint_probability=joint,conditional_survival=joint/pc if pc else None,wait_state_error_upper=ew,postmeasurement_subnormalized_error_upper=ez,sequential_joint_absolute_error_upper=min(1.,2*(ew+ez)))
    (P/f'N{n}_Q{q}_B{budget}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def probe(n,q,budget,wait,lag,steps=256):
    start=time.monotonic();sub,H,B,states,sizes,psi=model(n,q,budget,source=True)
    v,ew=error_evolve(H,B,psi,wait,steps,q/2)
    source=np.array([b&1 for g,b,E in states],float)
    detectors={}
    for g,_,_ in states:
        if g in detectors:continue
        adj=[set() for _ in range(n)]
        for e,(a,b) in enumerate(sub.es):
            if g>>e&1:adj[a].add(b);adj[b].add(a)
        dist=[n]*n;dist[0]=0;todo=[0]
        for a in todo:
            for b in adj[a]:
                if dist[b]==n:dist[b]=dist[a]+1;todo.append(b)
        detectors[g]=dist
    # A_r is a gauge-invariant detector density averaged over N-1 potential
    # targets, sorted by FINAL graph distance. Graph changes can contribute.
    O=np.array([[sum(((bits>>j)&1) for j in range(1,n) if detectors[g][j]==r)/(n-1) for g,bits,E in states] for r in range(1,n+1)])
    # Exact linear response: -i< [A_r(lag), n_0] >. No finite kick bias.
    gs={g:graphstats(n,g,sub.es) for g,_,_ in states}
    G=np.array([gs[g] for g,_,_ in states],float).T
    gn=['connected','triangles','squares','short_cycle_packing','degree_square_sum']+[f'ball_{r}' for r in range(n)]+[f'shell_{r}' for r in range(n)]
    C=(G[3]>=n//3).astype(float)
    V=np.column_stack([v,source*v,C*v])
    Z,ez=error_evolve(H,B,V,lag,steps,q/2)
    chi=2*np.imag(np.sum(Z[:,0].conj()[None,:]*O*Z[:,1][None,:],axis=1))
    # Both output columns inherit wait-state error (n_0 is a contraction).
    e0=ew+float(ez[0]);e1=ew+float(ez[1])
    # Exact response vector: projected expectation error <=2(e0+e1+e0*e1).
    bound=2*(e0+e1+e0*e1)
    pc=float(np.linalg.norm(C*v)**2);joint=float(np.linalg.norm(C*Z[:,2])**2)
    out=dict(N=n,Q=q,budget=budget,wait=wait,lag=lag,integration_steps=steps,dimension=len(states),retained_full_orbit_dimension=int(sum(sizes)),source_density=float(source@abs(v)**2),response_by_final_distance={str(r):float(x) for r,x in enumerate(chi,1)},wait_state_error_upper=ew,lag_column_error_upper=ez.tolist(),each_response_absolute_error_upper=bound,geometric_diagnostics_at_wait=dict(zip(gn,map(float,G@abs(v)**2))),geometric_diagnostics_at_wait_plus_lag=dict(zip(gn,map(float,G@abs(Z[:,0])**2))),loop_persistence=dict(packing_threshold=n//3,event_probability=pc,joint_probability=joint,conditional_survival=joint/pc if pc else None,joint_absolute_error_upper=min(1.,2*(ew+float(ez[2])))),source_gauge_invariant=True,detectors_gauge_invariant=True,seconds=time.monotonic()-start,interpretation='Linear response to exp(-i eta n_0); detector shell follows final graph. Includes graph redistribution and cannot alone certify matter transport or a causal cone.')
    (P/f'probe_N{n}_Q{q}_B{budget}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int);ap.add_argument('--budget',type=int,required=True);ap.add_argument('--end',type=float,default=2.);ap.add_argument('--steps',type=int,default=64)
    ap.add_argument('--probe',action='store_true');ap.add_argument('--wait',type=float,default=.1);ap.add_argument('--lag',type=float,default=.1)
    a=ap.parse_args()
    if a.probe:probe(a.n,a.q or 2*a.n,a.budget,a.wait,a.lag,a.steps)
    else:run(a.n,a.q or 2*a.n,a.budget,a.end,a.steps)
