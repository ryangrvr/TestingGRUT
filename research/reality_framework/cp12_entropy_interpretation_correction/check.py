#!/usr/bin/env python3
"""CP12 exact sparse graph isolation and regular entropy controls; stdlib."""
from itertools import combinations
from math import comb,lgamma,exp,log
def lc(n,k):
 if k<0 or k>n:return float('-inf')
 return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
def probabilities(n,q):
 m=comb(n,2)
 return (exp(lc(comb(n-1,2),q)-lc(m,q)),
         exp(lc(comb(n-2,2),q)-lc(m,q)))
def exact(n,q):
 edges=list(combinations(range(n),2))
 total=first=second=0
 for g in combinations(edges,q):
  deg=[0]*n
  for a,b in g:deg[a]+=1;deg[b]+=1
  i=deg.count(0);first+=i;second+=i*(i-1);total+=1
 p1,p2=probabilities(n,q)
 assert total==comb(comb(n,2),q)
 assert first==round(total*n*p1)
 assert second==round(total*n*(n-1)*p2)
 return (n,q,total,first/total)
def reg(n,k):
 q=k*n//2
 logreg=(lgamma(n*k+1)-lgamma(q+1)-q*log(2)
         -n*lgamma(k+1)-(k*k-1)/4)
 logall=lc(comb(n,2),q)
 return ((logreg-logall)/n,logreg/(n*log(n)))
def main():
 print('finite exhaustive',*[exact(n,q) for n,q in [(6,6),(7,7),(7,14)]],sep='\n')
 for n in (100,1000,10000,100000):
  for nu in (1,2,3):
   p,p2=probabilities(n,nu*n)
   var=n*p*(1-p)+n*(n-1)*(p2-p*p)
   assert var>=-1e-3
   assert abs(p-exp(-2*nu))<1/n**.4
   if n==10000:print('nu',nu,'isolates',n*p)
 for k in (2,4,6):
  c,s=reg(10000,k)
  assert c<0 and s>0
  print('regular degree',k,'log fraction per N',c,'entropy exponent est',s)
 print('CP12_ENTROPY_INTERPRETATION_CORRECTION_PASS')
if __name__=='__main__':main()
