#!/usr/bin/env python3
"""CP12 exact Q=N phase-mask topology-blindness at T^2, N=6.

Uses a separate Q(i sqrt(15)) Johnson-amplitude calculation and enumerates
all 5005 labelled six-edge graphs. Transfer of these neutral amplitudes to
the full-H T^2 coefficient is justified in parent's DERIVATION.md §2.
"""
from itertools import combinations, permutations
from collections import Counter
from fractions import Fraction as F
from math import comb
N=Q=6
EDGES=list(combinations(range(N),2)); M=len(EDGES); D=comb(M,Q)
A={tuple(sorted(e)) for e in [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,3)]}
B={tuple(sorted(e)) for e in [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3)]}
a=7;b=M-a;g=F(1,10*M)
z=(F(-7,8),F(1,8))
def plus(x,y):return x[0]+y[0],x[1]+y[1]
def scale(x,c):return x[0]*c,x[1]*c
def mul(x,y):return x[0]*y[0]-15*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def conj(x):return x[0],-x[1]
def p(k):return (Q-k)*(a-k)
def m(k):return k*(b-Q+k)
def apply(u):
 out={}
 for k in range(Q+1):
  v=scale(u[k],F(p(k)+m(k)))
  if p(k):v=plus(v,scale(mul(z,u[k+1]),F(-p(k))))
  if m(k):v=plus(v,scale(mul(conj(z),u[k-1]),F(-m(k))))
  out[k]=scale(v,g)
 return out
u0={k:(F(1),F(0)) for k in range(Q+1)}
u1=apply(u0);u2=apply(u1)
q2={k:mul(u1[k],conj(u1[k]))[0]-u2[k][0] for k in u0}
assert all(q2[k]==F(203*k*k-1133*k+1386,24000) for k in u0)
assert sum(comb(a,k)*comb(b,Q-k)*q2[k] for k in u0)==0

def degrees(E):
 d=[0]*N
 for x,y in E:d[x]+=1;d[y]+=1
 return d
def cyc(E,r):
 count=0
 for vs in combinations(range(N),r):
  for seq in permutations(vs[1:]):
   if seq[0]>seq[-1]:continue
   path=(vs[0],)+seq
   if all(tuple(sorted((path[j],path[(j+1)%r]))) in E for j in range(r)):count+=1
 return count
def connected(E):
 adj=[set() for _ in range(N)]
 for x,y in E:adj[x].add(y);adj[y].add(x)
 seen={0};stack=[0]
 while stack:
  for w in adj[stack.pop()]-seen:seen.add(w);stack.append(w)
 return len(seen)==N
def obs(E):
 d=degrees(E)
 return dict(triangles=cyc(E,3),four_cycles=cyc(E,4),
             degree_wedges=sum(comb(x,2) for x in d),
             connected=int(connected(E)),degree_two=int(all(x==2 for x in d)))
def evaluate(S):
 moments=Counter();init=Counter()
 for edges in combinations(EDGES,Q):
  E=set(edges);k=len(E&S)
  for name,v in obs(E).items():
   moments[name]+=v*q2[k];init[name]+=v
 return ({name:v/D for name,v in moments.items()},
         {name:F(v,D) for name,v in init.items()})
def main():
 assert len(A)==len(B)==7
 assert sorted(degrees(A))==sorted(degrees(B))==[2,2,2,2,3,3]
 assert sum(comb(d,2) for d in degrees(A))==10
 assert (cyc(A,3),cyc(A,4))==(0,2) and (cyc(B,3),cyc(B,4))==(2,0)
 coefs,init=evaluate(A);coefs2,init2=evaluate(B)
 assert coefs==coefs2 and init==init2
 assert coefs=={'triangles':F(-29,35750),'four_cycles':F(-87,357500),
                'degree_wedges':F(-29,13000),'connected':F(29,66000),
                'degree_two':F(203,4290000)}
 assert init=={'triangles':F(80,91),'four_cycles':F(45,91),
               'degree_wedges':F(60,7),'connected':F(732,1001),
               'degree_two':F(2,143)}
 print('CP12_PHASE_MOTIF_TOPOLOGY_BLINDNESS_EXACT_PASS')
 for key in coefs:print(key,'initial',init[key],'common_t2',coefs[key])
if __name__=='__main__':main()
