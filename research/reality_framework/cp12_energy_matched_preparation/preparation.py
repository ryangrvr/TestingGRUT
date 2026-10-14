"""Equal-mean-energy neutral preparations and exact initial graph response.
The large-N derivatives are full-CP12 coefficients through T^2 only.
Later-time two-color evolution is NOT substituted for full CP12 dynamics.
"""
from itertools import combinations
from math import comb,acos,log,sqrt,factorial
from fractions import Fraction as F
from pathlib import Path
import json,os
import numpy as np

P=Path(__file__).resolve().parent

def pattern(n,kind,seed=20261014):
    es=list(combinations(range(n),2));M=len(es);a=M//2
    if kind=='random':
        rng=np.random.default_rng(seed);idx=sorted(rng.choice(M,a,replace=False).tolist())
    elif kind=='concentrated':idx=list(range(a))
    elif kind=='spread':
        # Round-robin one-factorization, independent of any spatial lattice.
        assert n%2==0
        ring=list(range(n-1));fixed=n-1;picked=[]
        for _ in range(n-1):
            pairs=[(fixed,ring[0])]+[(ring[j],ring[-j]) for j in range(1,n//2)]
            picked.extend(tuple(sorted(e)) for e in pairs)
            ring=[ring[-1]]+ring[:-1]
        assert len(set(picked))==M
        ids={e:i for i,e in enumerate(es)};idx=sorted(ids[e] for e in picked[:a])
    else:raise ValueError(kind)
    S=np.zeros((n,n),dtype=np.int64)
    for i in idx:u,v=es[i];S[u,v]=S[v,u]=1
    return S

def theta_match(M,a):
    cosine=1-M*(M-1)/(2*a*(M-a))
    assert -1<=cosine<=1
    return acos(cosine)

def quotient(n,Q):
    M=comb(n,2);a=M//2;b=M-a;k=np.arange(max(0,Q-b),min(Q,a)+1)
    # Mode-centered probability recurrence avoids cancellation of large log-gammas.
    w=np.zeros(len(k));mode=min(max(((Q+1)*(a+1))//(M+2),k[0]),k[-1]);j=int(mode-k[0]);w[j]=1.
    for l in range(j,len(k)-1):
        x=k[l];w[l+1]=w[l]*(a-x)*(Q-x)/((x+1)*(b-Q+x+1))
    for l in range(j,0,-1):
        x=k[l]-1;w[l-1]=w[l]*(x+1)*(b-Q+x+1)/((a-x)*(Q-x))
    w/=w.sum()
    theta=theta_match(M,a)
    phases=np.r_[1.+0j,np.cumprod(np.full(len(k)-1,np.exp(1j*theta),complex))]
    psi=np.sqrt(w)*phases
    target=.1*Q*(M-Q)/M
    diag=.1/M*(k*(b-Q+k)+(Q-k)*(a-k))
    off=-.1/M*np.sqrt((k[:-1]+1)*(a-k[:-1])*(Q-k[:-1])*(b-Q+k[:-1]+1))
    def act(v,center=0.):
        out=(diag-center)*v;out=out.astype(complex)
        out[:-1]+=off*v[1:];out[1:]+=off*v[:-1];return out
    assert abs(np.vdot(psi,act(psi)).real-target)<2e-10*max(1,target)
    h1=act(psi,target);h2=act(h1,target)
    p1=2*np.imag(psi.conj()*h1)
    p2=abs(h1)**2-np.real(psi.conj()*h2)
    assert abs(p1.sum())<1e-9 and abs(p2.sum())<1e-9
    if a==b:
        direct=4*.1**2*(M-1)**2/M**4*(4*(M-1)*(k-Q/2)**2-Q*(M-Q))
        assert np.max(abs(p2-w*direct))<2e-10
    return dict(n=n,Q=Q,M=M,a=a,b=b,k=k,w=w,psi=psi,theta=theta,
                p1=p1,p2=p2,rewiring_energy=target,pure_rewiring_variance=float(np.vdot(h1,h1).real),
                mixed_rewiring_variance=.01*Q*(M-Q)/M**2)

def motif_colors(S):
    n=len(S);B=np.ones_like(S)-np.eye(n,dtype=np.int64)-S
    S2=S@S;B2=B@B;SB=S@B+B@S
    tr=lambda X,Y:int(np.sum(X*Y.T,dtype=np.int64))
    tri=[tr(B2,B)//6,tr(B2,S)//2,tr(S2,B)//2,tr(S2,S)//6]
    # Trace of (B+zS)^4, then remove weighted retracing walks.
    coeff=[0]*5;cs=[B2,SB,S2]
    for i in range(3):
        for j in range(3):coeff[i+j]+=tr(cs[i],cs[j])
    ds=S.sum(axis=1);db=B.sum(axis=1);a=int(ds.sum())//2;b=int(db.sum())//2
    wedges=[int(np.sum(db*(db-1)//2)),int(np.sum(ds*db)),int(np.sum(ds*(ds-1)//2))]
    coeff[0]-=2*b+4*wedges[0];coeff[2]-=4*wedges[1];coeff[4]-=2*a+4*wedges[2]
    assert all(x>=0 and x%8==0 for x in coeff)
    squares=[x//8 for x in coeff]
    assert sum(tri)==comb(n,3) and sum(squares)==3*comb(n,4)
    assert sum(wedges)==n*comb(n-1,2)
    return dict(wedges=wedges,triangles=tri,squares=squares)

def ratio(k,a,r):
    if not r:return np.ones_like(k,dtype=float)
    if a<r:return np.zeros_like(k,dtype=float)
    out=np.ones_like(k,dtype=float)
    for j in range(r):out*=(k-j)/(a-j)
    out[k<r]=0
    return out

def conditional(q,counts):
    r=len(counts)-1;k=q['k'];out=np.zeros_like(k,dtype=float)
    for j,c in enumerate(counts):out+=c*ratio(k,q['a'],j)*ratio(q['Q']-k,q['b'],r-j)
    return out


def falling(x,r):
    ans=1
    for j in range(r):ans*=x-j
    return ans

def exact_expectation(n,Q,counts):
    """Independent rational factorial-moment route; no quotient eigenvectors.
    q2(k) and the conditional motif are polynomials, so finite differences
    evaluate their hypergeometric average exactly.
    """
    M=comb(n,2);a=M//2;b=M-a;r=len(counts)-1;g=F(1,10*M)
    c=1-F(M*(M-1),2*a*b);ss=1-c*c
    def AB(k):
        plus=(Q-k)*(a-k);minus=k*(b-Q+k)
        return g*(plus+minus)*(1-c),-g*(plus-minus)
    def q2(k):
        plus=(Q-k)*(a-k);minus=k*(b-Q+k)
        A,B=AB(k);Ap,Bp=AB(k+1);Am,Bm=AB(k-1)
        return A*A+ss*B*B-g*((plus+minus)*A-plus*(c*Ap-ss*Bp)-minus*(c*Am+ss*Bm))
    def obs(k):
        return sum((F(num*falling(k,j)*falling(Q-k,r-j),falling(a,j)*falling(b,r-j)) for j,num in enumerate(counts) if num),F(0))
    def avg(fun,degree):
        vals=[fun(k) for k in range(degree+2)];out=F(0)
        for j in range(degree+1):
            out+=vals[0]*F(falling(Q,j)*falling(a,j),factorial(j)*falling(M,j))
            vals=[v-u for u,v in zip(vals,vals[1:])]
        assert vals==[0],('polynomial degree certificate',vals)
        return out
    initial=avg(obs,r);first_over_sqrt=avg(lambda k:2*AB(k)[1]*obs(k),r+2)
    second=avg(lambda k:q2(k)*obs(k),r+4)
    assert first_over_sqrt==0
    return initial,second

def exact_color_check():
    for n in [4,6,8]:
        for kind in ['random','concentrated','spread']:
            S=pattern(n,kind);computed=motif_colors(S);count=[0]*5
            for a,b,c,d in combinations(range(n),4):
                for cyc in [(a,b,c,d),(a,b,d,c),(a,c,b,d)]:
                    j=sum(S[x,y] for x,y in zip(cyc,cyc[1:]+cyc[:1]));count[j]+=1
            assert count==computed['squares']

def main():
    exact_color_check();rows=[]
    # One declared main connection density. Large-N entries are only derivatives at T=0.
    for n in [8,16,32,64,128,256]:
        Q=3*n;q=quotient(n,Q)
        for kind in ['random','concentrated','spread']:
            S=pattern(n,kind);cs=motif_colors(S);obs={}
            for name,counts in cs.items():
                O=conditional(q,counts);r=len(counts)-1
                expected=sum(counts)*np.prod([(Q-j)/(q['M']-j) for j in range(r)])
                initial=float(q['w']@O);first=float(q['p1']@O);second=float(q['p2']@O)
                assert abs(initial-expected)<1e-9*max(1,expected)
                assert abs(first)<1e-8*max(1,expected)
                exact0,exact2=exact_expectation(n,Q,counts)
                assert abs(initial-float(exact0))<1e-9*max(1,expected)
                assert abs(second-float(exact2))<2e-9
                obs[name]=dict(initial=float(exact0),T_coefficient=0.,T2_coefficient=float(exact2),
                               exact_initial=str(exact0),exact_T2_coefficient=str(exact2),
                               quotient_T2_error=abs(second-float(exact2)))
            row=dict(N=n,Q=Q,pattern=kind,seed=(20261014 if kind=='random' else None),
                     phase_mask_wedges=cs['wedges'][-1],phase_mask_degree_variance=float(np.var(S.sum(axis=1))),
                     phase_angle=q['theta'],mean_full_energy=Q/2+q['rewiring_energy'],
                     pure_full_energy_variance=2*Q+q['pure_rewiring_variance'],
                     mixture_full_energy_variance=2*Q+q['mixed_rewiring_variance'],
                     observables=obs)
            rows.append(row);print(json.dumps(row),flush=True)
    out=dict(status='CP12_EQUAL_ENERGY_PREPARATION_INITIAL_RESPONSE_PASS',rows=rows,
             scope='Full-law exact graph-observable derivatives through T^2 for supplied flat-amplitude phase states. No finite-time approximation or formation claim.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'preparation_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':main()
