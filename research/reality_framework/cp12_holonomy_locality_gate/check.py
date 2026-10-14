#!/usr/bin/env python3
"""Independent CP12 graph-cycle and Gauss-incidence rank controls; stdlib."""
from collections import Counter
from itertools import combinations
from math import factorial

def parts(n,ed):
    p=list(range(n))
    def root(x):
        while p[x]!=x:
            p[x]=p[p[x]];x=p[x]
        return x
    for a,b in ed:
        a,b=root(a),root(b)
        p[a]=b
    return len({root(i) for i in range(n)})

def unique_cycle_length(n,ed):
    adj=[[] for _ in range(n)];d=[0]*n
    for u,v in ed:
        adj[u].append(v);adj[v].append(u);d[u]+=1;d[v]+=1
    deleted=[x<=1 for x in d]
    queue=[v for v in range(n) if deleted[v]];idx=0
    while idx<len(queue):
        v=queue[idx];idx+=1
        for u in adj[v]:
            if not deleted[u]:
                d[u]-=1
                if d[u]<=1: deleted[u]=True;queue.append(u)
    return n-sum(deleted)

def count_formula(n,ell):
    numerator=factorial(n)*n**(n-ell)
    denominator=2*factorial(n-ell)*n
    assert numerator%denominator==0
    return numerator//denominator

def mat_rank(A):
    A=[list(map(float,row)) for row in A]
    if not A:return 0
    m,n=len(A),len(A[0]);r=0
    for j in range(n):
        pivot=next((i for i in range(r,m) if abs(A[i][j])>1e-9),None)
        if pivot is None: continue
        A[r],A[pivot]=A[pivot],A[r]
        v=A[r][j];A[r]=[x/v for x in A[r]]
        for i in range(r+1,m):
            a=A[i][j]
            if abs(a)>1e-9:A[i]=[x-a*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==m:break
    return r

def circulant(n,k):
    return sorted({tuple(sorted((v,(v+s)%n))) for v in range(n)
                   for s in range(1,k+1)})

def main():
    for n in (5,6,7):
        ed=list(combinations(range(n),2));counts=Counter();connected=0
        for ids in combinations(range(len(ed)),n):
            chosen=[ed[i] for i in ids]
            if parts(n,chosen)!=1:continue
            connected+=1
            counts[unique_cycle_length(n,chosen)]+=1
        expected={k:count_formula(n,k) for k in range(3,n+1)}
        assert dict(counts)==expected,(n,counts,expected)
        print('count',n,connected,dict(sorted(counts.items())))
    for n,ed in [(8,circulant(8,k)) for k in (1,2,3)]+[
        (6,[(0,1),(0,2),(1,2),(3,4),(3,5),(4,5)])]:
        q=len(ed);c=parts(n,ed)
        B=[[0]*q for _ in range(n)]
        for k,(a,b) in enumerate(ed):B[a][k]=1;B[b][k]=-1
        rank=mat_rank(B)
        assert rank==n-c and q-rank==q-n+c
        print('incidence',n,q,c,'free_integer_circulations',q-rank)
    print('CP12_HOLONOMY_LOCALITY_GRAPH_CHECK_PASS')
if __name__=='__main__': main()
