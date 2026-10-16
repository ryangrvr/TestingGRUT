"""Full-law phase response after a host-only record from the neutral ray."""
from math import sqrt,lcm
import argparse,gc,json,pickle,time
import numpy as np
from common import *
from graphs import metrics
from response_run import importance


def seed(n,q,b,wait,tag):
    target=cache_root()/f'{tag}_L0.pkl'
    if target.exists():
        with target.open('rb') as f:return pickle.load(f)
    m=quotient_model(n,q,b);sub=m['sub'];cat=json.loads((P/f'CATALOG_N{n}_Q{q}.json').read_text())
    choice=cat['source_choice'];canonical=choice['canonical_mask'];host=choice['recorded_host'];L=choice['labelled_count']
    z=np.load(neutral_path(n,q,b,wait));w=z['w'];e=float(z['state_error'])
    selected=np.array([s[0]==canonical for s in m['states']]);class_probability=float(np.linalg.norm(w*selected)**2)
    maps=[p for p in sub.group if sum(1<<p['edge'][j] for j in range(sub.M) if canonical>>j&1)==host]
    root=rooted_group(n,q,host);values={}
    for col in np.flatnonzero(selected):
        s=m['states'][col]
        for p in maps:
            t,phase=sub.transform(s,p);r,sign,size=root.canonical(t);assert sign
            value=w[col]*phase*sign*sqrt(size/m['sizes'][col])
            if r in values:assert abs(value-values[r])<1e-13
            else:values[r]=value
    states=list(values);X=np.array(list(values.values()),complex);norm=float(np.linalg.norm(X))
    assert abs(norm*norm-class_probability/L)<1e-12
    order=np.argsort(abs(X)**2,kind='stable');cumulative=np.cumsum(abs(X[order])**2)
    drop=int(np.searchsorted(cumulative,(norm*1e-10)**2,side='right'));omitted=float(np.linalg.norm(X[order[:drop]]))
    keep=sorted(order[drop:]);states=[states[i] for i in keep];X=X[keep];a=float(np.linalg.norm(X))
    error=e/sqrt(L)+omitted;rho=error/a;assert rho<1
    delta=sqrt(2)*rho/sqrt(1+sqrt(1-rho*rho));X/=a
    initial=dict(host=host,canonical_host=canonical,labelled_hosts_in_class=L,
        class_probability=class_probability,class_probability_interval=probability_interval(class_probability,e),
        one_labelled_host_probability=class_probability/L,
        fixed_source_peripheral_class_probability=class_probability*choice['source_orbit_size']/n,
        source_orbit_size=choice['source_orbit_size'],neutral_state_error=e,
        retained_orbits=len(states),untrimmed_orbits=len(values),discarded_initial_norm=omitted,
        normalized_initial_error=delta,rooted_group_order=len(root.group),lift_probability_error=abs(norm*norm-class_probability/L))
    del m,w,z;gc.collect();index={s:i for i,s in enumerate(states)};outside={};hh=Triplets();bb=Triplets()
    sizes=np.array([root.canonical(s)[2] for s in states]);scale=lcm(8,10*root.M)
    for col,s in enumerate(states):
        gauss(n,s)
        for t,c in transitions(n,s,scale):
            r,sign,size=root.canonical(t)
            if not sign:continue
            amplitude=c/scale*sign*sqrt(sizes[col]/size)
            if r in index:hh.add(index[r],col,amplitude)
            else:bb.add(outside.setdefault(r,len(outside)),col,amplitude)
    H=hh.csr((len(states),len(states)));B=bb.csr((len(outside),len(states)))
    he=float(np.linalg.norm((H-H.T).data));assert he<1e-8;root.cache.clear()
    model=dict(sub=root,states=states,outer_states=list(outside),sizes=sizes,H=H,B=B,X=X,
        initial=initial,level=0,build_checks=dict(H_hermiticity=he))
    save_cache(target,model);print(json.dumps(dict(stage='host_source_seed',**initial,dimension=len(states),outside=B.shape[0])),flush=True)
    return model


def measure(model,n,q,lags):
    H,B,X=model['H'],model['B'],model['X'];host=model['initial']['host'];ds=metrics(n,host)['distances'][0]
    bits=np.array([s[1] for s in model['states']]);n0=(bits&1).astype(float);Y=n0*X
    delta=model['initial']['normalized_initial_error'];mu=float(np.linalg.norm(Y)**2);muI=probability_interval(mu,delta)
    var=mu*(1-mu);varI=[min(v*(1-v) for v in muI),.25 if muI[0]<=.5<=muI[1] else max(v*(1-v) for v in muI)]
    assert varI[0]>0
    S=np.array([s[0]==host for s in model['states']],float)
    overlap=np.array([(s[0]&host).bit_count()/q for s in model['states']])
    es=list(combinations(range(n),2));neighborhood=sum(1<<j for j,(u,v) in enumerate(es) if host>>j&1 and 0 in (u,v))
    all_incident=sum(1<<j for j,(u,v) in enumerate(es) if 0 in (u,v))
    unchanged=np.array([s[0]&all_incident==neighborhood for s in model['states']],float)
    trajectory=MatrixFreePropagator(H,{'all_outgoing':B},q/2,dt=1/64,batch=4).run(np.column_stack([X,Y]),lags)
    rows=[];M=n*(n-1)//2;c=sqrt(q*(M-q))/(10*M)
    for t,(z,er) in trajectory.items():
        x,y=z[:,0],z[:,1];ex=delta+float(er['all_outgoing'][0]);ey=delta+float(er['all_outgoing'][1])
        pairs={}
        for j in range(1,n):
            orbit=sorted({p['vertex'][j] for p in model['sub'].group})
            O=np.mean([(bits>>k)&1 for k in orbit],axis=0)
            pairs[str(j)]=dict(distance=ds[j],symmetry_orbit=orbit,response=response(x,y,O,ex,ey))
        shells={}
        for d in range(1,max(ds)+1):
            vertices=[j for j in range(n) if ds[j]==d];O=np.mean([(bits>>j)&1 for j in vertices],axis=0)
            shells[str(d)]=dict(vertices=vertices,response=response(x,y,O,ex,ey))
        retention=float(np.vdot(x,S*x).real);edge=float(np.vdot(x,overlap*x).real)
        neigh=float(np.vdot(x,unchanged*x).real)
        row=dict(lag=t,receivers=pairs,shells=shells,host_retention=retention,
            host_retention_interval=probability_interval(retention,ex),
            retained_edge_fraction=edge,retained_edge_fraction_interval=probability_interval(edge,ex),
            source_neighborhood_retention=neigh,source_neighborhood_retention_interval=probability_interval(neigh,ex),
            universal_host_retention_lower_bound=max(0.,1-(c*t)**2),
            full_state_error=ex,source_vector_error=ey,
            frozen_host_full_state_difference_upper_bound=c*t,
            full_vs_frozen_host_response_difference_upper_bound=4*c*t+2*(c*t)**2)
        rows.append(row);print(json.dumps(dict(stage='host_source_response',N=n,lag=t,retention=row['host_retention_interval'],
            shells={d:v['response']['interval'] for d,v in shells.items()})),flush=True)
    result=dict(status='PASS',N=n,Q=q,dimension=H.shape[0],outside=B.shape[0],level=model['level'],
        initial=model['initial'],source_mean=mu,source_variance=var,source_variance_interval=varI,
        initial_host_distances=ds,initial_host_balls=metrics(n,host)['balls'],rows=rows,
        build_checks=model['build_checks'],
        evolution='Complete frozen H: all hopping, escape diagonals and graph rewiring retained. No circulating-core condition.',
        interpretation='All results condition on one host record chosen by the state-independent protocol; complete record work/entropy and original class probabilities are separate. Initial-host distances are fixed. Sampled retention is not first-passage survival.',
        error_scope='Every unretained physical row in all_outgoing; normalized initial error and analytic residual upper sums; ordinary floating point with reporting allowances.')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int,required=True)
    ap.add_argument('--budget',type=int,required=True);ap.add_argument('--wait',type=float,default=5)
    ap.add_argument('--level',type=int,default=6);ap.add_argument('--growth',type=int,default=100000)
    ap.add_argument('--lags',default='.125,.25,.5');a=ap.parse_args();start=time.monotonic()
    tag=f'HOST_N{a.n}_Q{a.q}_B{a.budget}_W{a.wait:g}_G{a.growth}'
    lags=sorted(float(t) for t in a.lags.split(','));model=seed(a.n,a.q,a.budget,a.wait,tag)
    for level in range(1,a.level+1):
        cached=cache_root()/f'{tag}_L{level}.pkl'
        if cached.exists():
            del model;gc.collect()
            with cached.open('rb') as f:model=pickle.load(f)
            continue
        count=min(a.growth,len(model['outer_states']),1200000-len(model['states']))
        if count<=0:break
        selected=np.arange(count) if count==len(model['outer_states']) else np.argsort(-importance(model,a.n,max(lags)),kind='stable')[:count]
        model=extend(model,a.n,a.q,selected,tag)
    result=measure(model,a.n,a.q,lags);result.update(wait=a.wait,budget=a.budget,seconds=time.monotonic()-start)
    (P/f'HOST_RESPONSE_N{a.n}_Q{a.q}_B{a.budget}_L{model["level"]}.json').write_text(json.dumps(result,indent=2)+'\n')
