"""Complete normalized Fock contractions between unequal receiver frames."""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, sqrt
import numpy as np


def two_mode(t, n):
    out = np.zeros((n+1, n+1), complex)
    for old in range(n+1):
        for new in range(n+1):
            for x in range(max(0, new-(n-old)), min(old, new)+1):
                y = new-x
                out[new, old] += (comb(old,x)*comb(n-old,y)
                    *t[0,0]**x*t[1,0]**(old-x)*t[0,1]**y*t[1,1]**(n-old-y)
                    *sqrt(factorial(new)*factorial(n-new)/(factorial(old)*factorial(n-old))))
    return out


class Space:
    def __init__(self, rank, photons):
        self.rank, self.k = rank, photons
        self.occ = list(combinations_with_replacement(range(rank),photons))
        self.index = {o:i for i,o in enumerate(self.occ)}
        self.size = len(self.occ)
        self.annihilation = []
        if photons:
            lower = space(rank,photons-1)
            for mode in range(rank):
                source, target, factor = [], [], []
                for j,o in enumerate(self.occ):
                    if mode in o:
                        old = list(o); old.remove(mode)
                        source.append(j); target.append(lower.index[tuple(old)])
                        factor.append(sqrt(o.count(mode)))
                self.annihilation.append((np.array(source,int),np.array(target,int),np.array(factor)))
        self.groups = {}

    def ann(self,c,x):
        out=np.zeros(space(self.rank,self.k-1).size,complex)
        for mode,(src,dst,factor) in enumerate(self.annihilation):
            out[dst] += c[mode].conjugate()*factor*x[src]
        return out

    def create(self,c,x):
        upper=space(self.rank,self.k+1)
        out=np.zeros(upper.size,complex)
        for mode,(dst,src,factor) in enumerate(upper.annihilation):
            out[dst] += c[mode]*factor*x[src]
        return out

    def density(self,x):
        if not self.k:
            return np.zeros((self.rank,self.rank),complex)
        reduced=np.zeros((space(self.rank,self.k-1).size,self.rank),complex)
        for mode,(src,dst,factor) in enumerate(self.annihilation):
            reduced[dst,mode]=factor*x[src]
        return reduced.T@reduced.conj()

    def rotation_groups(self,i,j):
        if (i,j) not in self.groups:
            groups={}
            for index,o in enumerate(self.occ):
                other=tuple(a for a in o if a not in (i,j))
                n=self.k-len(other)
                if n:
                    groups.setdefault((n,other),np.empty(n+1,int))[o.count(i)]=index
            self.groups[i,j]={n:np.stack([v for (nn,_),v in groups.items() if n==nn])
                              for n in range(1,self.k+1) if any(nn==n for nn,_ in groups)}
        return self.groups[i,j]

    def beam(self,x,i,j,t,matrices=None):
        out=x.copy()
        for n,index in self.rotation_groups(i,j).items():
            matrix=two_mode(t,n) if matrices is None else matrices[n]
            out[index]=out[index]@matrix.T
        return out


@lru_cache(maxsize=128)
def space(rank,photons):
    return Space(rank,photons)


def symmetric_power(g,k):
    """Dense exact polynomial; rectangular g need not be unitary."""
    rd,rs=g.shape
    if k==0:
        return np.ones((1,1),complex)
    source,dest=space(rs,k),space(rd,k)
    previous=symmetric_power(g,k-1)
    oldsrc=space(rs,k-1)
    last=np.array([o[-1] for o in source.occ],int)
    prev=np.array([oldsrc.index[o[:-1]] for o in source.occ],int)
    divisor=np.sqrt([o.count(o[-1]) for o in source.occ])
    out=np.zeros((dest.size,source.size),complex)
    for mode,(rows,oldrows,factor) in enumerate(dest.annihilation):
        out[rows,:] += (factor[:,None]*previous[np.ix_(oldrows,prev)]
                        *(g[mode,last]/divisor)[None,:])
    return out


def unitary_plan(u,k,prepare_rotations=True):
    """Gamma(u) through complex Givens rotations, with no dense Fock lift."""
    n=u.shape[0]
    a=u.copy(); rotations=[]
    for col in range(n):
        for row in range(n-1,col,-1):
            x,y=a[col,col],a[row,col]
            r=sqrt(abs(x)**2+abs(y)**2)
            if r==0:
                continue
            t=np.array([[x.conjugate(),y.conjugate()],[-y,x]])/r
            a[[col,row],:]=t@a[[col,row],:]
            inverse=t.conj().T
            rotations.append((col,row,inverse,
                              {m:two_mode(inverse,m) for m in range(1,k+1)}
                              if prepare_rotations else None))
    phases=np.array([np.prod(np.diag(a)[list(o)]) for o in space(n,k).occ])
    return phases,list(reversed(rotations))


def apply_unitary(plan,x,rank,k):
    phases,rotations=plan
    out=phases*x
    f=space(rank,k)
    for i,j,t,matrices in rotations:
        out=f.beam(out,i,j,t,matrices)
    return out


class Map:
    def __init__(self,g,k,dense_limit=65536):
        self.g,self.k=g,k
        self.rd,self.rs=g.shape
        self.shape=(space(self.rd,k).size,space(self.rs,k).size)
        self.dense_entries=0
        self.matrix=None
        if self.shape[0]*self.shape[1]<=dense_limit:
            self.matrix=symmetric_power(g,k)
            self.dense_entries=self.matrix.size
        else:
            u,s,vh=np.linalg.svd(g,full_matrices=True)
            self.left=unitary_plan(u,k)
            self.right=unitary_plan(vh,k)
            self.diagonal=np.zeros(min(self.rs,self.rd))
            self.diagonal[:]=s
            self.source=[]; self.dest=[]; self.scale=[]
            target=space(self.rd,k)
            for j,o in enumerate(space(self.rs,k).occ):
                if not o or max(o)<len(s):
                    self.source.append(j); self.dest.append(target.index[o])
                    self.scale.append(np.prod(s[list(o)]))
            self.source=np.array(self.source,int)
            self.dest=np.array(self.dest,int)
            self.scale=np.array(self.scale)

    def __call__(self,x):
        if self.matrix is not None:
            return self.matrix@x
        v=apply_unitary(self.right,x,self.rs,self.k)
        y=np.zeros(self.shape[0],complex)
        y[self.dest]=self.scale*v[self.source]
        return apply_unitary(self.left,y,self.rd,self.k)

    def adjoint(self,x):
        if self.matrix is not None:
            return self.matrix.conj().T@x
        # The plans are algebraic lifts; the adjoint lift is the lift of g^dagger.
        return Map(self.g.conj().T,self.k,dense_limit=0)(x)
