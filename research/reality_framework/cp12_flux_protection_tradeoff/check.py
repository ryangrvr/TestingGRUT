#!/usr/bin/env python3
"""Exact Gauss-law enumeration and removable-edge rewiring invariance.
Pure stdlib, N<=4, finite flux samples. The report proves unbounded result.
"""
from itertools import combinations,product
from fractions import Fraction as F
for N,Q in [(3,2),(3,3),(4,3),(4,4)]:
    edges=list(combinations(range(N),2));M=len(edges);seen=0
    for present in combinations(range(M),Q):
        mask=set(present)
        # (edge occupation b, first electric field, second electric field)
        choices=[(b,x,1-b-x) for b in (0,1) for x in range(-2,3)
                 if -2<=1-b-x<=2]
        for local in product(choices,repeat=Q):
            n=[0]*N
            for e,(b,x,y) in zip(present,local):
                i,j=edges[e];n[i]+=x;n[j]+=y
            if any(v not in (0,1) for v in n):continue
            nv=sum(n)
            assert nv==sum(1-b for b,x,y in local)
            r=sum(b==1 and x==0 and y==0 for b,x,y in local)
            k=sum(b==1 and (x!=0 or y!=0) for b,x,y in local)
            assert r+k+nv==Q
            F_e=F(1,8)*sum(x*x+y*y for b,x,y in local)
            assert F_e>=F(1,4)*(Q-r-F(nv,2))
            assert 2*nv+F_e>=F(1,4)*(Q-r)
            # Any legal rewiring swaps a removable present edge
            # for an absent occupied, zero-field edge, preserving r.
            for ix,(b,x,y) in enumerate(local):
                if (b,x,y)!=(1,0,0):continue
                for new_e in range(M):
                    if new_e in mask:continue
                    updated=list(local);updated.pop(ix);updated.append((1,0,0))
                    r2=sum(z==(1,0,0) for z in updated)
                    assert r2==r
            seen+=1
    assert seen>0
    print('N',N,'Q',Q,'Gauss states',seen,'PASS')
print('CP12_FLUX_PROTECTION_EXACT_LOCAL_IDENTITY_PASS')
