"""Inherited frozen-law tools and relational-core definitions."""
from pathlib import Path
from itertools import combinations,permutations
from math import factorial
import os,sys,tempfile
import numpy as np

P=Path(__file__).resolve().parent;RF=P.parent;ROOT=RF.parent.parent
sys.path.insert(0,str(RF/'cp12_dynamical_suppression_20261015'))
from boundary import quotient_model,Triplets,transitions,gauss
from operators import flow,hopping_edge,physical_model
from fast_symmetry import FastInvariant
from invariant_evolution import Invariant
from linear_propagator import MatrixFreePropagator
from integrals import Propagator
from stream import probability_interval
from resources import binary_entropy

ALLOWANCE=5e-10


def cache_root():
    path=Path(os.environ.get('GRUT_SOURCE_CACHE',str(Path(tempfile.gettempdir())/'grut_source_cache')))
    path.mkdir(parents=True,exist_ok=True);return path


def cycle_mask(n,order=None):
    order=list(range(n)) if order is None else list(order);es=list(combinations(range(n),2))
    return sum(1<<es.index(tuple(sorted((i,j)))) for i,j in zip(order,order[1:]+order[:1]))


def cycle_masks(n):
    return sorted({cycle_mask(n,(0,)+p) for p in permutations(range(1,n)) if p[0]<p[-1]})


def distances(n,C,source=0):
    adjacency=[set() for _ in range(n)]
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if C>>e&1:adjacency[i].add(j);adjacency[j].add(i)
    d=[None]*n;d[source]=0;todo=[source]
    for i in todo:
        for j in adjacency[i]:
            if d[j] is None:d[j]=d[i]+1;todo.append(j)
    return d


def is_spanning_cycle(n,C):
    if C.bit_count()!=n:return False
    degree=[0]*n
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if C>>e&1:degree[i]+=1;degree[j]+=1
    return degree==[2]*n and None not in distances(n,C)


def ring_present(n,s,C):
    if s[0]&C!=C:return False
    outgoing=[0]*n;incoming=[0]*n
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if not(C>>e&1):continue
        a,b=s[2][2*e:2*e+2]
        if a*b>=0:return False
        u,v=(i,j) if a>0 else (j,i);outgoing[u]+=1;incoming[v]+=1
    return outgoing==incoming==[1]*n


def rooted_group(n,q,C,host=None):
    sub=FastInvariant.__new__(FastInvariant)
    sub.n=n;sub.Q=q;sub.es=list(combinations(range(n),2));sub.M=len(sub.es)
    sub.group=[]
    for p in permutations(range(n)):
        if p[0]!=0:continue
        rep=sub.permutation(p)
        target=sum(1<<rep['edge'][e] for e in range(sub.M) if C>>e&1)
        host_ok=host is None or sum(1<<rep['edge'][e] for e in range(sub.M) if host>>e&1)==host
        if target==C and host_ok:sub.group.append(rep)
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)};sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for a,b in enumerate(p['vertex']):inv[b]=a
        sub.inverse.append(gid[tuple(inv)])
    sub.gcache={};sub.cache={}
    return sub


def complete_entropy_interval(pI,count):
    def entropy(p):return binary_entropy(p)+p*np.log2(count)
    lo,hi=pI;points=[entropy(lo),entropy(hi)]
    maximum=count/(count+1)
    if lo<=maximum<=hi:points.append(entropy(maximum))
    return [min(entropy(lo),entropy(hi)),max(points)]
