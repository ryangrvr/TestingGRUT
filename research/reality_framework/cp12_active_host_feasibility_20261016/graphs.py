"""Exact finite graphs, rooted source choice and uniform fixed-Q controls."""
from collections import Counter
from itertools import combinations,permutations
from math import comb,factorial,exp,log1p,sqrt,log
import argparse,json
import numpy as np
from common import P,FastInvariant


def metrics(n,g,max_radius=None):
    adj=[[] for _ in range(n)]
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if g>>e&1:adj[i].append(j);adj[j].append(i)
    ds=[];components=[];unused=set(range(n))
    for root in range(n):
        d=[None]*n;d[root]=0;todo=[root]
        for i in todo:
            for j in adj[i]:
                if d[j] is None:d[j]=d[i]+1;todo.append(j)
        ds.append(d)
        if root in unused:
            c={j for j,v in enumerate(d) if v is not None};components.append(len(c));unused-=c
    connected=len(components)==1;finite=[v for d in ds for v in d if v is not None]
    pairs=[v for i,d in enumerate(ds) for j,v in enumerate(d) if i<j and v is not None]
    return dict(connected=connected,components=sorted(components,reverse=True),
        diameter=max(finite) if connected else None,largest_component=max(components),
        maximum_component_diameter=max(finite),
        mean_finite_pair_distance=float(np.mean(pairs)) if pairs else 0.,
        disconnected_pair_fraction=1-len(pairs)/(n*(n-1)/2),
        balls=[sum(sum(v is not None and v<=r for v in d) for d in ds)/n for r in range(n if max_radius is None else max_radius+1)],
        degrees=list(map(len,adj)),distances=ds)


def aggregate(n,weighted):
    connected=0.;diameters={};balls=np.zeros(n);cb=np.zeros(n);pair=0.;dp=0.;degree=np.zeros(n)
    for g,p in weighted:
        a=metrics(n,g);balls+=p*np.array(a['balls']);dp+=p*a['disconnected_pair_fraction']
        for d in a['degrees']:degree[d]+=p/n
        if a['connected']:
            connected+=p;diameters[a['diameter']]=diameters.get(a['diameter'],0.)+p
            cb+=p*np.array(a['balls']);pair+=p*a['mean_finite_pair_distance']
    return dict(connectivity_probability=connected,disconnected_probability=1-connected,
        diameter_distribution={str(d):v for d,v in sorted(diameters.items())},mean_balls=balls.tolist(),
        mean_balls_given_connected=(cb/connected).tolist() if connected else None,
        mean_pair_distance_given_connected=pair/connected if connected else None,
        disconnected_pair_probability=dp,degree_distribution=degree.tolist())


def full_group(n,q):
    sub=FastInvariant.__new__(FastInvariant);sub.n=n;sub.Q=q;sub.M=n*(n-1)//2
    sub.es=list(combinations(range(n),2));sub.group=[sub.permutation(p) for p in permutations(range(n))]
    lookup={tuple(p['vertex']):i for i,p in enumerate(sub.group)};sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for i,j in enumerate(p['vertex']):inv[j]=i
        sub.inverse.append(lookup[tuple(inv)])
    sub.gcache={};sub.cache={};return sub


def catalog(n,q):
    sub=full_group(n,q);counts=Counter(sub.graph_info(sum(1<<e for e in es))[0]
        for es in combinations(range(sub.M),q));rows=[]
    for g,count in sorted(counts.items()):
        a=metrics(n,g);a.pop('distances');aut=sub.graph_info(g)[1]
        assert count*len(aut)==factorial(n)
        rows.append(dict(mask=g,labelled_count=count,automorphism_order=len(aut),**a))
    connected=[r for r in rows if r['connected']];D=max(r['diameter'] for r in connected)
    chosen=min((r for r in connected if r['diameter']==D),key=lambda r:(-r['labelled_count'],r['mask']))
    g=chosen['mask'];a=metrics(n,g);aut=sub.graph_info(g)[1]
    peripheral={i for i,d in enumerate(a['distances']) if max(d)==D};orbits=[]
    while peripheral:
        i=min(peripheral);o=sorted({sub.group[k]['vertex'][i] for k in aut});peripheral-=set(o);orbits.append(o)
    orbit=min(orbits,key=lambda o:(-len(o),o[0]));source=orbit[0];perm=[None]*n;perm[source]=0
    for label,i in enumerate(j for j in range(n) if j!=source):perm[i]=label+1
    p=sub.permutation(tuple(perm));target=sum(1<<p['edge'][e] for e in range(sub.M) if g>>e&1)
    total=comb(sub.M,q)
    data=dict(status='PASS',N=n,Q=q,labelled_hosts=total,classes=rows,
        baseline=aggregate(n,[(r['mask'],r['labelled_count']/total) for r in rows]),
        source_choice=dict(canonical_mask=g,recorded_host=target,diameter=D,
            labelled_count=chosen['labelled_count'],source_orbit=orbit,source_orbit_size=len(orbit),
            distances=metrics(n,target)['distances'][0],
            prior_class_probability=chosen['labelled_count']/total,
            prior_fixed_source_class_probability=chosen['labelled_count']/total*len(orbit)/n))
    (P/f'CATALOG_N{n}_Q{q}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(N=n,Q=q,classes=len(rows),baseline=data['baseline'],source=data['source_choice'])),flush=True)
    return data


def isolate_moments(n,q):
    M=n*(n-1)//2
    def avoid(k):return exp(sum(log1p(-k/(M-j)) for j in range(q))) if q<=M-k else 0.
    p1=avoid(n-1);p2=avoid(2*n-3);mean=n*p1;var=mean+n*(n-1)*p2-mean*mean
    return dict(N=n,Q=q,expected_isolated_vertices=mean,isolated_count_variance=var,
        connected_probability_upper_bound=min(1.,max(0.,var)/(mean*mean)) if mean else 1.)


def baselines():
    rows=[];rng=np.random.default_rng(20261016);samples=512
    # Each row's connected fraction and normalized ball means are bounded
    # variables. Union bound over rows and listed quantities at 95% confidence.
    sizes=[8,12,16,24,32,64,128];quantities=len(sizes)*2*8
    eps=sqrt(log(2*quantities/.05)/(2*samples))
    for family in ['N+1','2N']:
        for n in sizes:
            q=n+1 if family=='N+1' else 2*n;M=n*(n-1)//2
            values=[];connected=0;diam=Counter()
            for _ in range(samples):
                g=sum(1<<int(e) for e in rng.choice(M,q,replace=False));a=metrics(n,g,4)
                values.append(a['balls'][:5]);connected+=a['connected']
                if a['connected']:diam[a['diameter']]+=1
            balls=np.mean(values,axis=0);p=connected/samples
            rows.append(dict(family=family,N=n,Q=q,samples=samples,connected=p,
                connectivity_interval=[max(0.,p-eps),min(1.,p+eps)],
                mean_balls=balls.tolist(),ball_intervals=[[max(1.,v-n*eps),min(n,v+n*eps)] for v in balls],
                connected_diameter_counts={str(k):v for k,v in sorted(diam.items())},
                isolate_moments=isolate_moments(n,q)))
            print(json.dumps(dict(stage='uniform_baseline',family=family,N=n,connected=p,balls=balls.tolist())),flush=True)
    data=dict(status='PASS',kind='Uniform fixed-Q combinatorial baseline, not CP-12 dynamics',seed=20261016,
        simultaneous_confidence=.95,hoeffding_radius=eps,rows=rows,
        exact_isolate_moments=[isolate_moments(n,q) for n in [128,512,2048,8192] for q in [n+1,2*n]])
    (P/'LARGE_N_BASELINE.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int);ap.add_argument('--q',type=int)
    ap.add_argument('--sample',action='store_true');a=ap.parse_args()
    if a.sample:baselines()
    else:catalog(a.n,a.q)
