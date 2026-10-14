"""Exact phase-mask topology response in one main sector Q=3N.

Coefficients through T^3 are full frozen-Hamiltonian coefficients for invariant
graph-diagonal observables, using the proven projection identities. No large-N
many-body time evolution is asserted.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb,factorial,acos,sqrt
from pathlib import Path
import hashlib,json,os,sys
import numpy as np

P=Path(__file__).resolve().parent
OLD=P.parent/'cp12_energy_matched_preparation'
assert hashlib.sha256((OLD/'preparation.py').read_bytes()).hexdigest()=='def88f44d6ab15a23f585eabf448352887de355b23b66adcfab31b6fc267feef'
sys.path.insert(0,str(OLD))
from preparation import motif_colors,falling

A6={(0,1),(1,2),(2,3),(3,4),(4,5),(0,5),(0,3)}
B6={(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3)}


def matching_order(n):
    assert n%2==0
    ring=list(range(n-1));fixed=n-1;out=[]
    for _ in range(n-1):
        out.extend([(fixed,ring[0])]+[(ring[j],ring[-j]) for j in range(1,n//2)])
        ring=[ring[-1]]+ring[:-1]
    out=[tuple(sorted(e)) for e in out]
    assert len(set(out))==comb(n,2)
    return out


def masks(n):
    assert n%12==0
    blocks=n//6;M=comb(n,2)
    target=int(F(M,2)-7*blocks+18)//36
    assert 0<=target<=comb(blocks,2)
    coarse=set(matching_order(blocks)[:target])
    pair=[]
    for base in [A6,B6]:
        S=np.zeros((n,n),np.int64)
        for b in range(blocks):
            for i,j in base:S[6*b+i,6*b+j]=S[6*b+j,6*b+i]=1
        for b,c in coarse:
            S[6*b:6*b+6,6*c:6*c+6]=1
            S[6*c:6*c+6,6*b:6*b+6]=1
        pair.append(S)
    A,B=pair;a=int(A.sum())//2
    assert a==int(B.sum())//2
    assert np.array_equal(A.sum(axis=1),B.sum(axis=1))
    assort=[]
    for X in pair:
        d=X.sum(axis=1)
        assort.append(sum(int((d[i]-1)*(d[j]-1)) for i,j in combinations(range(n),2) if X[i,j]))
    assert assort[0]==assort[1]
    assert -1<=1-F(M*(M-1),2*a*(M-a))<=1
    return A,B,dict(N=n,M=M,blocks=blocks,marked_edges=a,coarse_edges=target,
                   common_marked_edge_degree_product=assort[0])


def probability_polynomials(M,Q,a):
    b=M-a;g=F(1,10*M);c=1-F(M*(M-1),2*a*b);ss=1-c*c
    assert ss>0
    def add(x,y):return x[0]+y[0],x[1]+y[1]
    def mul(x,y):return x[0]*y[0]-ss*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
    def scale(t,x):return t*x[0],t*x[1]
    @lru_cache(None)
    def u(order,k):
        if order==0:return F(1),F(0)
        p=(Q-k)*(a-k);m=k*(b-Q+k)
        out=scale(F(p+m),u(order-1,k))
        out=add(out,scale(F(-p),mul((c,F(1)),u(order-1,k+1))))
        out=add(out,scale(F(-m),mul((c,F(-1)),u(order-1,k-1))))
        return scale(g,out)
    def q(order,k):
        out=F(0)
        for j in range(order+1):
            x=u(j,k);y=u(order-j,k);prod=mul((x[0],-x[1]),y)
            phase=(2*j-order)%4
            term=prod[0] if phase==0 else -prod[0] if phase==2 else -prod[1] if phase==1 else prod[1]
            out+=term/F(factorial(j)*factorial(order-j))
        return out
    return q,c,ss


def hypergeometric_average(fun,degree,M,Q,a):
    vals=[fun(k) for k in range(degree+2)];out=F(0)
    for j in range(degree+1):
        out+=vals[0]*F(falling(Q,j)*falling(a,j),factorial(j)*falling(M,j))
        vals=[v-u for u,v in zip(vals,vals[1:])]
    assert vals==[0],('degree certificate',degree,vals)
    return out


def coefficients(M,Q,a,counts):
    q,c,ss=probability_polynomials(M,Q,a);b=M-a;r=len(counts)-1
    def obs(k):
        return sum((F(num*falling(k,j)*falling(Q-k,r-j),falling(a,j)*falling(b,r-j)) for j,num in enumerate(counts) if num),F(0))
    ans=[]
    for j in range(4):
        ans.append(hypergeometric_average(lambda k:q(j,k)*obs(k),r+j,M,Q,a))
    assert ans[1]==0
    return ans,c,ss


def energy_variance(M,Q,a):
    b=M-a;g=F(1,10*M);c=1-F(M*(M-1),2*a*b);ss=1-c*c
    def square(k):
        p=(Q-k)*(a-k);m=k*(b-Q+k)
        return g*g*((p+m)**2*(1-c)**2+ss*(p-m)**2)
    Rmean=g*Q*(M-Q)
    return 2*Q+hypergeometric_average(square,4,M,Q,a)-Rmean**2


def cubic_leader(M,c):
    # Symbolically verified by symbolic_check.py; q3 is divided by sin(theta).
    return F(1,10*M)**3*F(4,3)*(1-c)*(8*(1-3*c)-2*M*(1-13*c)-3*M*M*(1+3*c)+M**3*(1+c))


def run(n):
    A,B,meta=masks(n);a=meta['marked_edges'];M=meta['M'];Q=3*n
    counts=[motif_colors(X) for X in [A,B]]
    assert counts[0]['triangles'][-1]-counts[1]['triangles'][-1]==-2*meta['blocks']
    assert counts[0]['squares'][-1]-counts[1]['squares'][-1]==2*meta['blocks']
    assert counts[0]['wedges'][-1]==counts[1]['wedges'][-1]
    obs={}
    for name in ['wedges','triangles','squares']:
        co=[]
        for row in counts:
            cert,c,ss=coefficients(M,Q,a,row[name]);co.append(cert)
        assert co[0][:3]==co[1][:3]
        delta=co[1][3]-co[0][3]
        obs[name]=dict(initial=str(co[0][0]),common_T2=str(co[0][2]),
                       A_T3_over_sin=str(co[0][3]),B_T3_over_sin=str(co[1][3]),
                       B_minus_A_T3_over_sin=str(delta),B_minus_A_T3=float(delta)*sqrt(float(ss)),
                       phase_mask_counts=[row[name][-1] for row in counts])
    assert obs['triangles']['B_minus_A_T3']!=0 and obs['squares']['B_minus_A_T3']!=0
    factor=F(6*falling(Q,3)*falling(M-Q,3),falling(M,6))
    delta=cubic_leader(M,c)*(2*meta['blocks'])*factor
    assert F(obs['triangles']['B_minus_A_T3_over_sin'])==delta
    assert F(obs['squares']['B_minus_A_T3_over_sin'])==-F(3*(Q-3),M-6)*delta
    assert obs['wedges']['B_minus_A_T3_over_sin']=='0'
    out=dict(**meta,Q=Q,phase_cosine=str(c),phase_sine_squared=str(ss),
             mean_energy=str(F(Q,2)+F(Q*(M-Q),10*M)),
             full_energy_variance=str(energy_variance(M,Q,a)),
             min_entropy=0,coherence_bits=__import__('math').log2(comb(M,Q)),
             marked_degrees=list(map(int,A.sum(axis=1))),observables=obs)
    print(json.dumps(dict(N=n,Q=Q,marked_edges=a,topology_response=obs),indent=2),flush=True)
    return out


def main():
    out=dict(status='CP12_DEGREE_MATCHED_TOPOLOGY_DERIVATIVES_PASS',models=[run(n) for n in [24,48,96,192]],
             scope='Q=3N throughout. Exact full-H graph-observable coefficients through T^3, by analytical projection identities. No full finite-time simulation, uniform Taylor remainder, loop-packing formation or propagation claim.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'topology_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])

if __name__=='__main__':main()
