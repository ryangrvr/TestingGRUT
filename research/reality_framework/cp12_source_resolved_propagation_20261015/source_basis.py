"""Signed rooted-core lift and full-law numerical Galerkin bases."""
from collections import defaultdict,deque
from math import sqrt,lcm,factorial
import gc,json,pickle,time
import numpy as np
from scipy.sparse import bmat,csr_matrix,vstack,hstack
from shared import *


def save_cache(path,model):
    temporary=path.with_suffix('.partial')
    with temporary.open('wb') as f:pickle.dump(model,f)
    temporary.replace(path)


def lift(n,q,budget,wait,tail_relative=1e-11,chord=None):
    start=time.monotonic();m=quotient_model(n,q,budget);sub=m['sub'];C=cycle_mask(n)
    key=f'N{n}_Q{q}_B{budget}_W{wait:g}';z=np.load(cache_root()/f'{key}.npz')
    w=z['w'];mask=z['mask'].copy();cores=z['core'];ew=float(z['state_error'])
    all_cycle_probability=float(np.linalg.norm(mask*w)**2);host=None
    if chord is not None:
        assert q==n+1 and 2<=chord<=n//2
        edge=list(combinations(range(n),2)).index((0,chord));host=C|(1<<edge)
        for col in np.flatnonzero(mask):
            s=m['states'][col];extra=s[0]^int(cores[col]);assert extra.bit_count()==1
            i,j=m['sub'].es[extra.bit_length()-1]
            if distances(n,int(cores[col]),i)[j]!=chord:mask[col]=False
    root=rooted_group(n,q,C,host);maps={};values={};old_ids={}
    for col in np.flatnonzero(mask):
        s=m['states'][col];cs=int(cores[col])
        mapkey=(cs,s[0] if host is not None else None)
        if mapkey not in maps:
            maps[mapkey]=[p for p in sub.group
                if sum(1<<p['edge'][e] for e in range(sub.M) if cs>>e&1)==C and
                (host is None or sum(1<<p['edge'][e] for e in range(sub.M) if s[0]>>e&1)==host)]
            assert len(maps[mapkey])==(2*n if host is None else 4 if 2*chord==n else 2)
        for p in maps[mapkey]:
            t,phase=sub.transform(s,p);r,sign,size=root.canonical(t);assert sign
            value=w[col]*phase*sign*sqrt(size/m['sizes'][col])
            if r in values:assert abs(values[r]-value)<1e-13
            else:values[r]=value;old_ids[r]=int(col)
    states=list(values);X=np.array(list(values.values()),complex)
    count=factorial(n-1)//2
    if chord is not None:count*=n//2 if 2*chord==n else n
    p_total=float(np.linalg.norm(mask*w)**2);a=float(np.linalg.norm(X))
    assert abs(a*a-p_total/count)<1e-13
    order=np.argsort(abs(X)**2,kind='stable');cumulative=np.cumsum(abs(X[order])**2)
    drop=int(np.searchsorted(cumulative,(tail_relative*a)**2,side='right'))
    omitted=float(np.linalg.norm(X[order[:drop]]));keep=sorted(order[drop:])
    states=[states[i] for i in keep];X=X[keep]
    initial_error=ew/sqrt(count)+omitted
    a_ret=float(np.linalg.norm(X));rho=initial_error/a_ret
    shape_error=(sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho))) if rho<1 else 2.
    # A global phase is free in a conditional pure ray. A ball of radius e
    # about X_hat subtends sin(angle)<=e/||X_hat|| after phase alignment.
    X/=a_ret
    root.cache.clear()
    info=dict(neutral_cache=key,neutral_state_error=ew,core_mask=C,cycle_records=count,
        rooted_group_order=len(root.group),initial_untrimmed_orbits=len(values),
        record_type='core_only' if host is None else 'core_and_initial_host',
        initial_host_mask=host,recorded_chord_core_distance=chord,
        all_spanning_cycle_probability=all_cycle_probability,
        equivalent_record_count=count,
        fixed_source_endpoint_class_probability=None if host is None else 2*p_total/n,
        initial_retained_orbits=len(states),discarded_initial_norm=omitted,
        original_total_formation_probability=p_total,
        original_total_formation_interval=probability_interval(p_total,ew),
        probability_per_core=p_total/count,original_projected_error_per_core=ew/sqrt(count),
        total_unnormalized_initial_error=initial_error,
        normalized_phase_aligned_initial_error=shape_error,initial_approximate_norm=a_ret,
        lift_probability_error=abs(a*a-p_total/count),seconds=time.monotonic()-start)
    return root,states,X,info


def build_seed(n,q,budget,wait,tag,chord=None):
    path=cache_root()/f'{tag}_L0.pkl'
    if path.exists():
        with path.open('rb') as f:return pickle.load(f)
    sub,states,X,info=lift(n,q,budget,wait,chord=chord);index={s:i for i,s in enumerate(states)}
    sizes=np.array([sub.canonical(s)[2] for s in states]);outside={};hh=Triplets();bb=Triplets()
    sc=lcm(8,10*sub.M)
    for col,s in enumerate(states):
        gauss(n,s)
        for t,c in transitions(n,s,sc):
            r,sign,size=sub.canonical(t)
            if not sign:continue
            a=c/sc*sign*sqrt(sizes[col]/size)
            if r in index:hh.add(index[r],col,a)
            else:bb.add(outside.setdefault(r,len(outside)),col,a)
    H=hh.csr((len(states),len(states)));B=bb.csr((len(outside),len(states)))
    err=float(np.linalg.norm((H-H.T).data));assert err<1e-9
    sub.cache.clear()
    model=dict(sub=sub,states=states,outer_states=list(outside),sizes=sizes,H=H,B=B,
        X=X,initial=info,level=0,build_checks=dict(H_hermiticity=err))
    save_cache(path,model)
    print(json.dumps(dict(stage='source_seed',tag=tag,dimension=len(states),outside=len(outside),**info)),flush=True)
    return model


def extend(model,n,q,selected,tag):
    level=model['level']+1;path=cache_root()/f'{tag}_L{level}.pkl'
    if path.exists():
        with path.open('rb') as f:return pickle.load(f)
    start=time.monotonic();sub=model['sub'];old=model['states'];outer=model['outer_states']
    selected=np.sort(np.asarray(selected,np.int64));taken=set(map(int,selected))
    extra=[outer[i] for i in selected];nb=len(old);ne=len(extra);states=old+extra
    index={s:i for i,s in enumerate(states)}
    remaining=[i for i in range(len(outer)) if i not in taken]
    outside=[outer[i] for i in remaining];oi={s:i for i,s in enumerate(outside)}
    sizes=np.r_[model['sizes'],np.array([sub.canonical(s)[2] for s in extra])]
    hh=Triplets();bb=Triplets();cross=Triplets();sc=lcm(8,10*sub.M);last=start
    for col,s in enumerate(extra):
        gauss(n,s)
        for t,c in transitions(n,s,sc):
            r,sign,size=sub.canonical(t)
            if not sign:continue
            a=c/sc*sign*sqrt(sizes[nb+col]/size);row=index.get(r)
            if row is not None:
                if row<nb:cross.add(row,col,a)
                else:hh.add(row-nb,col,a)
            else:bb.add(oi.setdefault(r,len(oi)),col,a)
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='source_extend',tag=tag,level=level,processed=col,
                total=ne,outside=len(oi))),flush=True);last=time.monotonic()
    # Append new outside states in their exact assigned order.
    outside=[None]*len(oi)
    for s,i in oi.items():outside[i]=s
    C=model['B'][selected,:].T.tocsr();ce=float(np.linalg.norm((cross.csr((nb,ne))-C).data));assert ce<1e-9
    H=bmat([[model['H'],C],[C.T,hh.csr((ne,ne))]],format='csr')
    B0=vstack([model['B'][remaining,:],csr_matrix((len(outside)-len(remaining),nb))],format='csr')
    B=hstack([B0,bb.csr((len(outside),ne))],format='csr')
    he=float(np.linalg.norm((H-H.T).data));assert he<1e-8;sub.cache.clear()
    result=dict(sub=sub,states=states,outer_states=outside,sizes=sizes,H=H,B=B,
        X=np.r_[model['X'],np.zeros(ne,complex)],initial=model['initial'],level=level,
        build_checks=dict(H_hermiticity=he,independent_new_to_old_error=ce,
            old_dimension=nb,added=ne,dimension=H.shape[0],outside=B.shape[0],seconds=time.monotonic()-start))
    save_cache(path,result)
    print(json.dumps(dict(stage='source_extended',tag=tag,level=level,**result['build_checks'])),flush=True)
    return result
