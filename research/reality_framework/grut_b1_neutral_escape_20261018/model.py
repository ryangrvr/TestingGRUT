"""B1 exact physical matrices. No relative-distance or field-energy cutoff.

Only the exactly conserved resource sector is reduced. Color factorization
is used for the primary one-resource controls, not assumed at density.
"""
from itertools import combinations
import numpy as np
from scipy.sparse import coo_matrix, diags

JP, JM, KAPPA = 1/64, 1/96, 1/4
GAP = 5/4
SELF = (JP*JP+JM*JM)/GAP
JPAIR = 2*JP*JM/GAP
OMEGA = 2-2*SELF


def field(plus, minus, length):
    charge = 0; result = []
    for x in range(length-1):
        charge += int(x in plus)-int(x in minus)
        result.append(charge)
    assert charge+int(length-1 in plus)-int(length-1 in minus) == 0
    return tuple(result)


class PathModel:
    def __init__(self, length, sites=(0,), resource=1, omega=OMEGA, lam=JPAIR):
        self.length, self.sites, self.resource = length, tuple(sites), resource
        self.states = []; self.fields = []
        for fuels in range(1 << len(sites)):
            pairs = resource-fuels.bit_count()
            if not 0 <= pairs <= length: continue
            for plus in combinations(range(length), pairs):
                for minus in combinations(range(length), pairs):
                    e = field(plus, minus, length)
                    if max(map(abs,e), default=0) <= 1:
                        self.states.append((plus, minus, fuels)); self.fields.append(e)
        self.index = {s:i for i,s in enumerate(self.states)}
        self.dim = len(self.states)
        diagonal = []; hops = []; conversion = []
        def put(target, col, value, dest): dest.append((self.index[target],col,value))
        for col,(plus,minus,fuels) in enumerate(self.states):
            density = sum(int(x in plus and x+1 in minus)+int(x in minus and x+1 in plus)
                          for x in range(length-1))
            correction = sum((2-int(x>0)-int(x<length-1)) for x in set(plus)&set(minus))
            diagonal.append(2*len(plus)+sum(x*x for x in self.fields[col])
                            +KAPPA*density-SELF*correction+omega*fuels.bit_count())
            for species, slots, coupling in [(0,plus,JP),(1,minus,JM)]:
                for old in slots:
                    for new in (old-1,old+1):
                        if not 0 <= new < length or new in slots: continue
                        moved = tuple(sorted(new if x==old else x for x in slots))
                        pp,mm = (moved,minus) if species==0 else (plus,moved)
                        e2 = field(pp,mm,length)
                        if max(map(abs,e2),default=0)>1: continue
                        change = [x for x,(a,b) in enumerate(zip(self.fields[col],e2)) if a!=b]
                        assert change == [min(old,new)]
                        put((pp,mm,fuels),col,-coupling,hops)
            for a,x in enumerate(sites):
                if not fuels>>a&1 or x in plus or x in minus: continue
                pp = tuple(sorted(plus+(x,))); mm = tuple(sorted(minus+(x,)))
                assert field(pp,mm,length) == self.fields[col]
                target = (pp,mm,fuels^(1<<a))
                put(target,col,-lam,conversion)
                conversion.append((col,self.index[target],-lam))
        def sparse(entries):
            if not entries: return coo_matrix((self.dim,self.dim),dtype=float).tocsr()
            rows,cols,values = zip(*entries)
            return coo_matrix((values,(rows,cols)),shape=(self.dim,self.dim)).tocsr()
        self.D = diags(diagonal,format='csr'); self.T = sparse(hops)
        self.C = sparse(conversion); self.H = self.D+self.T+self.C
        self.pairs = np.array([len(p) for p,m,f in self.states])
        self.fuels = np.array([f.bit_count() for p,m,f in self.states])
        assert np.all(self.pairs+self.fuels==resource)
        assert np.max(abs((self.H-self.H.T).data),initial=0)<1e-15

    def initial(self):
        v = np.zeros(self.dim,complex)
        v[self.index[((),(),(1<<len(self.sites))-1)]] = 1
        return v

    def intertwiner(self):
        assert self.resource==1 and self.sites==(0,)
        ids = [self.index[((),(),1)]]+[self.index[((x,),(x,),0)] for x in range(self.length)]
        other = [i for i in range(self.dim) if i not in set(ids)]
        h = self.H.toarray()-2*np.eye(self.dim)
        ap = h[np.ix_(ids,ids)]; b = h[np.ix_(other,ids)]
        qh = h[np.ix_(other,other)]; diagonal = np.diag(qh)
        inv = 1/diagonal
        tq = qh-np.diag(diagonal)
        k = ap-b.T@(inv[:,None]*b)
        x1 = -inv[:,None]*b
        x2 = -inv[:,None]*(tq@x1)
        x3 = -inv[:,None]*(tq@x2)
        xk = -(inv*inv)[:,None]*(b@k)
        x = x1+x2+x3+xk
        z = np.zeros((self.dim,len(ids)))
        z[ids,:] = np.eye(len(ids)); z[other,:] = x
        rp = b.T@(x3+xk)
        rq = tq@(x3+xk)-x2@k-(x3+xk)@k
        residual = h@z-z@k
        formula = np.zeros_like(residual); formula[ids,:] = rp; formula[other,:] = rq
        return {'ids':ids,'other':other,'K':k,'Z':z,'residual':residual,'formula':formula,
                'parity_zero':b.T@x2}
