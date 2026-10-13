#!/usr/bin/env python3
"""CP12 Q=N coherent matching trial and exact regular-graph cut certificates.
Outward 90-digit arithmetic; NO many-body diagonalization or electric cutoff.
"""
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from math import comb,isqrt
def dec(x,rounding):
 with localcontext() as c:
  c.prec=90;c.rounding=rounding
  return Decimal(x.numerator)/Decimal(x.denominator)
def certificate(n):
 M=comb(n,2);q=n
 center=round(n/7);width=isqrt(n)
 ids=list(range(max(0,center-width),min(n,center+width)+1))
 w=[width+1-abs(k-center) for k in ids]
 norm=sum(i*i for i in w)
 diag=F(n,2)+F(17,8)*F(sum(k*v*v for k,v in zip(ids,w)),norm)
 # H_k recurrence from exact matching counts m_k=binom(n,k)*H_k.
 h=[1,n-1]
 for j in range(1,ids[-1]+1):h.append((n-1)*h[-1]-j*h[-2])
 with localcontext() as c:
  c.prec=90;c.rounding=ROUND_FLOOR
  term=Decimal(0)
  for j,k in enumerate(ids[:-1]):
   c2=F((k+1)*(n-k)*(q-k)*h[k+1],(M-k)*h[k])
   root=dec(c2,ROUND_FLOOR).sqrt().next_minus()
   term+=Decimal(w[j]*w[j+1])*root
  off=2*term/Decimal(norm)
 with localcontext() as c:
  c.prec=90;c.rounding=ROUND_CEILING
  E=dec(diag,ROUND_CEILING)-off
  gap=(-Decimal(27)*Decimal(n)/Decimal(250))-E
  b2=F(3*q*(M-q),100*M*M)
  p=min(Decimal(1),dec(b2,ROUND_CEILING)/gap**2) if gap>0 else Decimal(1)
 assert gap>0
 return n, E/Decimal(n),p
if __name__=="__main__":
 for n in (42,64,128,256,512,1024,2048):
  N,e,p=certificate(n)
  print("N=",N,"trial_E_per_N_upper=",e,"ground_degree2_P_upper=",p)
 print("CP12_QN_OUTWARD_RITZ_AND_GROUND_BOUND_PASS")
