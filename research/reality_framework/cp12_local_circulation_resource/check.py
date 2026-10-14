"""Exact finite graph controls and explicit asymptotic CP-12 bounds.
No graph preference or coupling is added to the frozen law.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial, log, pi, sqrt
from pathlib import Path
import json,os

P=Path(__file__).resolve().parent

def cycles(n, edges, R):
    adj=[set() for _ in range(n)]
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    out=[]
    for root in range(n):
        def visit(path):
            last=path[-1]
            if len(path)>=3 and root in adj[last] and path[1]<last:
                out.append(tuple(path))
            if len(path)==R:return
            for v in sorted(adj[last]):
                if v>root and v not in path:visit(path+[v])
        visit([root])
    return out

def packings(cs,k):
    masks=[sum(1<<v for v in c) for c in cs]
    def count(start,used,left):
        if left==0:return 1
        return sum(count(i+1,used|masks[i],left-1) for i in range(start,len(cs)) if not(used&masks[i]))
    return count(0,0,k)

def expectation_exact(n,q,R,k):
    M=comb(n,2); ans=F(0)
    for lens in product(range(3,R+1),repeat=k):
        ell=sum(lens)
        if ell>n or ell>q:continue
        ordered=F(factorial(n),factorial(n-ell))
        for r in lens:ordered/=2*r
        ordered/=factorial(k)
        inclusion=F(comb(M-ell,q-ell),comb(M,q))
        ans+=ordered*inclusion
    return ans

def union_bound(n,q,R,k):
    M=comb(n,2)
    A=sum((F(n*q,M)**r/F(2*r) for r in range(3,R+1)),F(0))
    return A**k/factorial(k)

def torus(lengths):
    points=list(product(*(range(a) for a in lengths))); ids={x:i for i,x in enumerate(points)}
    es=set();faces=[]
    for x in points:
        for d in range(len(lengths)):
            y=list(x);y[d]=(y[d]+1)%lengths[d]
            es.add(tuple(sorted((ids[x],ids[tuple(y)]))))
        for a,b in combinations(range(len(lengths)),2):
            xa=list(x);xa[a]=(xa[a]+1)%lengths[a]
            xb=list(x);xb[b]=(xb[b]+1)%lengths[b]
            xab=xa.copy();xab[b]=(xab[b]+1)%lengths[b]
            faces.append((ids[x],ids[tuple(xa)],ids[tuple(xab)],ids[tuple(xb)]))
    return len(points),sorted(es),faces

def rational_rank(rows):
    pivots={}
    for row in rows:
        r={i:F(x) for i,x in enumerate(row) if x}
        while r:
            j=min(r)
            if j not in pivots:
                v=r[j];pivots[j]={i:x/v for i,x in r.items()};break
            a=r[j]
            for i,x in pivots[j].items():
                z=r.get(i,F(0))-a*x
                if z:r[i]=z
                else:r.pop(i,None)
    return len(pivots)

def loop_vectors(es,cs):
    index={e:i for i,e in enumerate(es)};out=[]
    for c in cs:
        row=[0]*len(es)
        for a,b in zip(c,c[1:]+c[:1]):row[index[tuple(sorted((a,b)))]]+=1 if a<b else -1
        out.append(row)
    return out

def asymptotic_bound(logN,nu,beta=1.,R=4,alpha=.25):
    # Conservative upper constants at N>=3: NQ/M <= 3nu.
    A=sum((3*nu)**r/(2*r) for r in range(3,R+1))
    # log C(N+Q,Q) <= N[(1+nu)log(1+nu)-nu log nu].
    theta_upper=1+sqrt(pi/(beta*.25/2))
    f=(1+nu)*log(1+nu)-nu*log(nu)+2*nu*log(theta_upper)
    ell_per_N=alpha*(logN+log(alpha)-1-log(A))-f
    omega_per_N=4*nu+2*.1*nu
    return dict(logN=logN,nu=nu,R=R,k_over_N=alpha,
                conservative_log_inverse_reference_bound_per_N=ell_per_N,
                gibbs_probability_upper=(min(1.,beta*omega_per_N/ell_per_N+1e-12) if ell_per_N>0 else 1.),
                asymptotic_only=True)

def main():
    exact=[]
    for n,q,R,k in [(5,5,4,1),(6,6,4,2),(6,8,4,2),(7,14,4,2)]:
        es=list(combinations(range(n),2));hist=Counter();total_packings=0
        for picked in combinations(es,q):
            a=packings(cycles(n,picked,R),k);hist[bool(a)]+=1;total_packings+=a
        denom=comb(len(es),q);actual=F(total_packings,denom);event=F(hist[True],denom)
        expected=expectation_exact(n,q,R,k);upper=union_bound(n,q,R,k)
        assert actual==expected and event<=expected<=upper
        row=dict(N=n,Q=q,R=R,k=k,graphs=denom,event_count=hist[True],
                 packing_sum=total_packings,packing_expectation=str(actual),
                 exact_event_probability=str(event),union_upper=str(upper))
        exact.append(row);print(json.dumps(row),flush=True)
    tori=[]
    for shape in [(6,6),(5,5,5)]:
        n,es,faces=torus(shape);rank=rational_rank(loop_vectors(es,faces));d=len(shape)
        assert rank==len(es)-n+1-d
        row=dict(shape=list(shape),N=n,Q=len(es),cycle_rank=len(es)-n+1,
                 elementary_faces=len(faces),rational_local_face_rank=rank,
                 remaining_winding_rank=d)
        tori.append(row);print(json.dumps(row),flush=True)
    # Disjoint plaquettes covering all vertices: same physical density in each dimension.
    tile=[]
    for shape in [(8,8),(4,4,4)]:
        n,es,faces=torus(shape);coords=list(product(*(range(a) for a in shape)))
        ids={x:i for i,x in enumerate(coords)};chosen=[]
        for x in coords:
            if x[0]%2 or x[1]%2:continue
            a=list(x);a[0]+=1;b=list(x);b[1]+=1;ab=a.copy();ab[1]+=1
            chosen.append((ids[x],ids[tuple(a)],ids[tuple(ab)],ids[tuple(b)]))
        assert len(chosen)==n//4 and len({v for c in chosen for v in c})==n
        assert all(set(tuple(sorted(e)) for e in zip(c,c[1:]+c[:1]))<=set(es) for c in chosen)
        tile.append(dict(shape=list(shape),N=n,Q=len(es),vertex_disjoint_four_cycles=len(chosen)))
    out=dict(status='CP12_LOCAL_CIRCULATION_GRAPH_CONTROLS_PASS',enumerations=exact,
             tori=tori,disjoint_plaquette_controls=tile,
             illustrative_asymptotic_bounds=[asymptotic_bound(x,nu) for nu in [2.,3.] for x in [50.,100.,300.]],
             scope='Exact graph counts and rational incidence ranks; prepared tori are benchmarks, not selected states.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'graph_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':main()
