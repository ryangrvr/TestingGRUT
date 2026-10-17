"""Exact small fixed-path A1 matrices and unchanged-P1 neutral controls.

Fixed-resource reduction is an exact invariant, not numerical truncation.
All field-compatible charge positions and reverse transitions are retained.
"""
from itertools import combinations
import math
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

JP, JM, KAPPA = 1/64, 1/96, 1/4


def field(plus, minus, length):
    charge=0; result=[]
    for x in range(length-1):
        charge+=int(x in plus)-int(x in minus)
        result.append(charge)
    assert charge+int(length-1 in plus)-int(length-1 in minus)==0
    return tuple(result)


def moved(slots, old, new):
    return tuple(sorted(new if x==old else x for x in slots))


class PathModel:
    def __init__(self, length, matching=(), resource=1, omega=None, lam=None):
        self.length,self.matching,self.resource=length,tuple(matching),resource
        k=len(matching)
        self.omega=np.array(omega if omega is not None else [13/4]*k)
        self.lam=np.array(lam if lam is not None else [1.]*k)
        self.states=[]; self.fields=[]
        for f in range(1<<k):
            n=resource-f.bit_count()
            if n<0 or n>length: continue
            for plus in combinations(range(length),n):
                for minus in combinations(range(length),n):
                    e=field(plus,minus,length)
                    if max(map(abs,e),default=0)<=1:
                        self.states.append((plus,minus,f)); self.fields.append(e)
        self.index={s:i for i,s in enumerate(self.states)}
        self.fields=np.array(self.fields,int)
        self.dim=len(self.states)
        diagonal=[]; hops=[]; converts=[]
        def put(target,col,val,entries): entries.append((self.index[target],col,val))
        for col,(plus,minus,f) in enumerate(self.states):
            density=sum(int(x in plus and x+1 in minus)+int(x in minus and x+1 in plus)
                        for x in range(length-1))
            diagonal.append(2*len(plus)+sum(e*e for e in self.fields[col])
                            +KAPPA*density+sum(om for a,om in enumerate(self.omega) if f>>a&1))
            for species,slots,jhop in [(0,plus,JP),(1,minus,JM)]:
                for old in slots:
                    for new in (old-1,old+1):
                        if not 0<=new<length or new in slots: continue
                        pp,mm=(moved(plus,old,new),minus) if species==0 else (plus,moved(minus,old,new))
                        e2=field(pp,mm,length)
                        if max(map(abs,e2),default=0)>1: continue
                        changes=[x for x,(a,b) in enumerate(zip(self.fields[col],e2)) if a!=b]
                        assert changes==[min(old,new)]
                        expected=(-1 if new>old else 1)*(1 if species==0 else -1)
                        assert e2[changes[0]]-self.fields[col][changes[0]]==expected
                        put((pp,mm,f),col,-jhop,hops)
            for a,((i,j),coupling) in enumerate(zip(matching,self.lam)):
                assert j==i+1
                if not f>>a&1 or i in plus or j in minus: continue
                pp,mm=tuple(sorted(plus+(i,))),tuple(sorted(minus+(j,)))
                e2=field(pp,mm,length)
                if max(map(abs,e2),default=0)>1: continue
                assert e2[i]-self.fields[col][i]==1
                target=(pp,mm,f^(1<<a))
                put(target,col,coupling,converts)
                converts.append((col,self.index[target],coupling))
        def matrix(entries):
            if not entries: return coo_matrix((self.dim,self.dim),dtype=complex).tocsr()
            row,col,val=zip(*entries)
            return coo_matrix((np.array(val,complex),(row,col)),shape=(self.dim,self.dim)).tocsr()
        self.T=matrix(hops); self.C=matrix(converts)
        self.D=diags(diagonal,format='csr',dtype=complex)
        self.H=self.D+self.T+self.C
        assert np.max(abs((self.H-self.H.getH()).data),initial=0)<1e-13
        self.particles=np.array([2*len(p) for p,m,f in self.states],float)
        self.fcount=np.array([f.bit_count() for p,m,f in self.states],float)
        assert np.all(self.particles==2*(resource-self.fcount))

    def initial(self):
        v=np.zeros(self.dim,complex)
        v[self.index[((),(),(1<<len(self.matching))-1)]]=1
        return v

    def evolve(self,v,t): return expm_multiply((-1j*t)*self.H,v)


def exact_cell_matrix():
    # Independently written order F,D+,D-,M_left,M_right.
    h=np.diag([13/4,13/4,13/4,2,2]).astype(complex)
    for a,b,value in [(0,1,1),(1,3,-JM),(1,4,-JP),(2,3,-JP),(2,4,-JM)]:
        h[a,b]=h[b,a]=value
    return h


def log_tail_upper(k,x):
    # Chernoff bound tau_k(x)<=exp(e*x-k). Used only as a conservative
    # floating formula evaluation. Keep a positive floor below float range.
    return max(-700.,math.e*x-k)


def geometry_bound(n,t,a=None,L=None):
    a=n//4 if a is None else a
    L=3*n/16-1 if L is None else L
    c=4/(1-math.exp(-1)); v=12*(math.e-1)
    # For small-order diagnostics sum the positive tail directly.
    if a<50:
        term=(1.1*t)**(a+1)/math.factorial(a+1); tail=term
        for r in range(a+2,a+200):
            term*=1.1*t/r; tail+=term
            if term<max(1e-300,1e-17*tail): break
    else:
        tail=math.exp(log_tail_upper(a+1,1.1*t))
    bterm=14*math.sqrt(c)/v*math.exp(max(-700.,-L/2))*math.expm1(v*t/2)
    return min(2.,tail+bterm)


def neutral_constants():
    tau=2*(JP+JM); gap=5/4
    rem0=(tau*tau/gap)*(tau/gap)**2/(1-(tau/gap)**2)
    d0=rem0+tau**4/(gap-tau)**3
    g1=.75; c1=g1-tau; b1=tau*tau/(g1-2*tau)
    rem1=tau*tau/g1*(tau/g1)**2/(1-(tau/g1)**2)
    d1=rem1+tau*tau*b1/(c1*(c1-b1))
    tau2=2*tau; b2=tau2*tau2/(.5-tau2)
    return {'tau':tau,'delta0':gap,'ground_remainder_at_E2':rem0,'ground_dispersion_error':d0,
            'ground_hopping':1/3840,'excited_dispersion_error':d1,
            'excited_even_hopping':-1/23040,'excited_odd_hopping':1/4608,
            'low_band_overlap_probability_upper':1/529,
            'two_pair_self_energy_upper':b2,'collision_diagonal_difference_lower':.5-2*b2}
