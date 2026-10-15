"""Source-resolved propagation at distances fixed by the initial core record."""
from pathlib import Path
from math import sqrt
import argparse,gc,json,time,pickle
import numpy as np
from scipy.sparse.linalg import expm_multiply
from shared import *
from source_basis import build_seed,extend
from linear_propagator import shifted_operators


def ratio(I,J):
    if J[0]<=0:return None
    values=[a/b for a in I for b in J];return [min(values),max(values)]


def product(I,J):
    values=[a*b for a in I for b in J];return [min(values),max(values)]


def absolute(I):
    return [0. if I[0]<=0<=I[1] else min(abs(a) for a in I),max(abs(a) for a in I)]


def response(x,y,O,ex,ey):
    value=float(2*np.vdot(x,O*y).imag)
    radius=2*(ex*np.linalg.norm(O*y)+ey*np.linalg.norm(O*x)+ex*ey)+ALLOWANCE
    interval=[value-radius,value+radius]
    plus=x+(-1j-1)*y;minus=x+(1j-1)*y
    finite=float((np.vdot(plus,O*plus).real-np.vdot(minus,O*minus).real)/2)
    assert abs(finite-value)<1e-11
    return dict(value=value,interval=interval,absolute_interval=absolute(interval),
        zero_excluded=bool(abs(value)>radius),error_bound=float(radius),
        exact_finite_phase_identity_error=abs(finite-value))


def stats(model,n):
    C=model['initial']['core_mask'];d=distances(n,C);bits=np.array([s[1] for s in model['states']])
    n0=(bits&1).astype(float)
    O={r:np.mean(np.array([(bits>>j)&1 for j in range(n) if d[j]==r]),axis=0)
        for r in range(1,max(d)+1)}
    S=np.array([ring_present(n,s,C) for s in model['states']],float)
    return n0,O,S,d


def importance(model,n,T):
    X=model['X'];n0=np.array([s[1]&1 for s in model['states']],float)
    V=np.column_stack([X,n0*X]);A,G,tr=shifted_operators(model['H'],model['sub'].Q/2)
    grid=expm_multiply(G,V,start=0,stop=T,num=3,traceA=tr)
    score=np.zeros(model['B'].shape[0])
    for z in grid:score+=np.sum(abs(model['B']@z)**2,axis=1)
    return score


def measure(model,n,q,lags,dt):
    H,B,X=model['H'],model['B'],model['X'];n0,Os,S,d=stats(model,n);Y=n0*X
    e0=model['initial']['normalized_phase_aligned_initial_error']
    mu=float(np.linalg.norm(Y)**2);muI=probability_interval(mu,e0)
    var=mu*(1-mu);varI=[min(a*(1-a) for a in muI),
        .25 if muI[0]<=.5<=muI[1] else max(a*(1-a) for a in muI)]
    assert varI[0]>0,'No resolved active density-phase source.'
    results=MatrixFreePropagator(H,{'all_outgoing':B},q/2,dt=dt,batch=4).run(np.column_stack([X,Y]),lags)
    rows=[]
    for lag,(z,errors) in results.items():
        x,y=z[:,0],z[:,1];ex=e0+float(errors['all_outgoing'][0]);ey=e0+float(errors['all_outgoing'][1])
        sp=float(np.vdot(x,S*x).real);spI=probability_interval(sp,ex)
        sr=response(x,y,S,ex,ey);shells={}
        for r,O in Os.items():
            raw=response(x,y,O,ex,ey);joint=response(x,y,O*S,ex,ey)
            density=float(np.vdot(x,O*S*x).real);densityI=probability_interval(density,ex)
            if spI[0]>0:
                first=ratio(joint['interval'],spI)
                second=ratio(product(densityI,sr['interval']),[spI[0]**2,spI[1]**2])
                conditional=[first[0]-second[1],first[1]-second[0]]
                cv=joint['value']/sp-density*sr['value']/sp**2
            else:conditional=None;cv=None
            shells[str(r)]=dict(vertices=[j for j in range(n) if d[j]==r],
                density_response=raw,density_and_same_core_response=joint,
                density_given_same_core_response=cv,
                density_given_same_core_response_interval=conditional,
                final_density_and_core=density,final_density_and_core_interval=densityI)
        near=shells['1']['density_response']['absolute_interval']
        far_ratios={r:ratio(v['density_response']['absolute_interval'],near) for r,v in shells.items() if r!='1'}
        pairs={}
        bits=np.array([s[1] for s in model['states']])
        for j in range(1,n):
            orbit=sorted({p['vertex'][j] for p in model['sub'].group})
            O=np.mean(np.array([(bits>>k)&1 for k in orbit]),axis=0)
            raw=response(x,y,O,ex,ey);joint=response(x,y,O*S,ex,ey)
            pairs[str(j)]=dict(core_distance=d[j],symmetry_orbit=orbit,
                density_response=raw,density_and_same_core_response=joint)
        row=dict(lag=lag,shells=shells,fixed_receiver_pairs=pairs,same_core_probability=sp,
            same_core_probability_interval=spI,same_core_response=sr,
            far_to_neighbor_absolute_response_intervals=far_ratios,
            state_error_X=ex,state_error_source=ey,lag_errors=errors['all_outgoing'].tolist())
        rows.append(row)
        print(json.dumps(dict(stage='distance_response',N=n,Q=q,level=model['level'],lag=lag,
            source_variance_interval=varI,same_core=spI,
            responses={r:v['density_response']['interval'] for r,v in shells.items()},
            far_neighbor=far_ratios)),flush=True)
    # Expectations of initial host balls are diagnostics; no host measurement.
    host_data={g:distances(n,g) for g in {s[0] for s in model['states']}}
    host_balls=[];host_ball_intervals=[]
    for radius in range(max(d)+1):
        counts={g:sum(v is not None and v<=radius for v in ds) for g,ds in host_data.items()}
        obs=np.array([counts[s[0]] for s in model['states']],float)
        value=float(np.vdot(X,obs*X).real);host_balls.append(value)
        host_ball_intervals.append([max(0.,sqrt(value)-sqrt(n)*e0)**2,
            min(float(n),(sqrt(value)+sqrt(n)*e0)**2)])
    result=dict(status='PASS',N=n,Q=q,level=model['level'],dimension=H.shape[0],outside_dimension=B.shape[0],
        H_nnz=H.nnz,initial=model['initial'],source_mean=mu,source_mean_interval=muI,
        source_variance=var,source_variance_interval=varI,
        fixed_core_mask=model['initial']['core_mask'],fixed_core_distances=d,
        fixed_core_ball_volumes=[sum(v<=r for v in d) for r in range(max(d)+1)],
        initial_mean_active_host_ball_volumes=host_balls,initial_active_host_ball_intervals=host_ball_intervals,
        build_checks=model['build_checks'],
        source='Gauge-invariant exp(-i theta n0), applied after recording the core; variance demonstrably nonzero.',
        receivers='Distance-shell means and separately fixed receivers. For each fixed receiver, average only over its exact initial-core/host/source stabilizer orbit, which has identical response by symmetry. No final-distance binning.',
        evolution='Complete frozen Hamiltonian including all active external hopping and all graph rewiring.',
        conditioning='All quantities are conditional on the initial core record, with separately reported original formation probability, complete record work/entropy and normalized state error. No final-core selection is hidden: both joint and final-conditional responses are reported, and the final-conditional derivative subtracts its survival-denominator response.',
        error_scope='All outgoing physical rows, including rewiring; analytic residual upper sums; phase-aligned conditional initial-state bound; bilinear response intervals. Ordinary floating point with5e-10 reporting allowance, not outward rounding.',
        rows=rows,dt=dt)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int,required=True)
    ap.add_argument('--budget',type=int,required=True);ap.add_argument('--wait',type=float,required=True)
    ap.add_argument('--level',type=int,default=3);ap.add_argument('--growth',type=int,default=200000)
    ap.add_argument('--chord',type=int,help='Also record the initial host, with the extra edge from0 to this core distance; Q=N+1 only.')
    ap.add_argument('--lags',default='.125,.25,.5');ap.add_argument('--dt',type=float,default=1/64)
    args=ap.parse_args();start=time.monotonic();n,q=args.n,args.q
    tag=f'N{n}_Q{q}_B{args.budget}_W{args.wait:g}_G{args.growth}'
    if args.chord is not None:tag+=f'_F{args.chord}'
    lags=sorted(float(x) for x in args.lags.split(','));model=build_seed(n,q,args.budget,args.wait,tag,chord=args.chord)
    for level in range(1,args.level+1):
        cached=cache_root()/f'{tag}_L{level}.pkl'
        if cached.exists():
            del model;gc.collect()
            with cached.open('rb') as f:model=pickle.load(f)
            continue
        count=min(args.growth,len(model['outer_states']))
        if count==len(model['outer_states']):selected=np.arange(count)
        else:selected=np.argsort(-importance(model,n,max(lags)),kind='stable')[:count]
        model=extend(model,n,q,selected,tag)
    result=measure(model,n,q,lags,args.dt);result['wait']=args.wait
    result['neutral_budget']=args.budget;result['seconds']=time.monotonic()-start
    out=P/f'RESPONSE_{tag}_L{args.level}.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='source_done',file=out.name,seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
