#!/usr/bin/env python3
"""Independent full Gauss-law ring basis and two-star CP12 operator check.

Original construction; no import of the 32-state local certificate.
Requires NumPy and SciPy. Uses exact frozen global Fock orientation signs.
"""
import itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

def fields(n, vertex, holes, z):
    left=[z]*n
    for i in range(1,n):
        left[i]=left[i-1]+((vertex>>i)&1)-((holes>>(i-1))&1)
    assert left[0]+((holes>>(n-1))&1)-left[n-1]==(vertex&1)
    right=[((holes>>i)&1)-left[i] for i in range(n)]
    return left,right

def fock_sign(n, v, h, i, e):
    bits=v+((((1<<n)-1)^h)<<n)
    src=n+e;dst=i
    tmp=bits^(1<<src)
    return (-1)**((bits&((1<<src)-1)).bit_count()+(tmp&((1<<dst)-1)).bit_count())

def calculate(n,L,a=19/200):
    assert n>=4 and L>=1
    states=[]
    for k in range(n+1):
        for vi in itertools.combinations(range(n),k):
            v=sum(1<<x for x in vi)
            for hi in itertools.combinations(range(n),k):
                h=sum(1<<x for x in hi)
                for z in range(-L,L+1):
                    left,right=fields(n,v,h,z)
                    if all(abs(x)<=L for x in left+right):states.append((v,h,z))
    index={x:i for i,x in enumerate(states)}
    rr=[];cc=[];vv=[];spectators=[]
    for j,(v,h,z) in enumerate(states):
        left,right=fields(n,v,h,z)
        neutral=[int(not (h>>e&1) and left[e]==right[e]==0) for e in range(3)]
        electrical=(left[0]**2+right[0]**2+2*(left[1]**2+right[1]**2)+left[2]**2+right[2]**2)/16
        diagonal=1+2*((v>>1&1)+(v>>2&1))+electrical+a*(neutral[0]/2+neutral[1]+neutral[2]/2)
        rr.append(j);cc.append(j);vv.append(diagonal)
        spectators.append((v & ~6,h & ~7,z,left[0],right[2]))
        for e,i,orientation in ((0,1,-1),(1,1,1),(1,2,-1),(2,2,1)):
            if (h>>e&1) or (v>>i&1):continue
            target=(v^(1<<i),h^(1<<e),z)
            t=index.get(target)
            if t is None:continue
            sign=orientation*fock_sign(n,v,h,i,e)
            rr.extend((j,t));cc.extend((t,j));vv.extend((sign,sign))
    H=coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr()
    assert (H-H.T).nnz==0
    for i,j in zip(*H.nonzero()):
        if i!=j:assert spectators[i]==spectators[j]
    lowest=float(eigsh(H,k=1,which="SA",tol=1e-11,return_eigenvectors=False)[0])
    assert lowest>-27/125
    return dict(N=n,L=L,dimension=len(states),lowest_pair_energy=lowest,
                margin_above_minus_27_over_125=lowest+27/125,
                independent_spectator_blocks=len(set(spectators)))

if __name__=="__main__":
    for n,L in ((4,2),(4,3),(4,4),(5,2),(5,3)):
        print(calculate(n,L),flush=True)
    print("CP12_INDEPENDENT_PHYSICAL_PAIR_EMBEDDING_PASS")
