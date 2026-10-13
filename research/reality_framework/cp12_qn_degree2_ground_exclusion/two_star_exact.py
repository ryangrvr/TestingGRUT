#!/usr/bin/env python3
"""Exact rational CP12 two-star inequality at a=19/200; full integer rotor tail.
Vertex Gauss at two central vertices, edge-center Gauss at three edges;
outer vertex Gauss deliberately relaxed. No finite-electric-cutoff premise.
"""
from fractions import Fraction as F
A=F(19,200)
LOWER=-F(27,125)

def fields(bits,z):
 n1=bits&1;n2=(bits>>1)&1
 h0=1-((bits>>2)&1);h1=1-((bits>>3)&1);h2=1-((bits>>4)&1)
 e0o=z; e0i=h0-z
 e1l=n1-e0i; e1r=h1-e1l
 e2i=n2-e1r; e2o=h2-e2i
 return (e0o,e0i,e1l,e1r,e2i,e2o)

def hop(bits,src,dst):
 if not (bits>>src)&1 or (bits>>dst)&1:return None
 mid=bits^(1<<src)
 s=(-1)**((bits&((1<<src)-1)).bit_count()+(mid&((1<<dst)-1)).bit_count())
 return mid^(1<<dst),s

def block(z):
 H=[[F(0) for _ in range(32)] for _ in range(32)]
 for bits in range(32):
  e=fields(bits,z)
  holes=[1-((bits>>(i+2))&1) for i in range(3)]
  neutral=[int(not holes[i] and e[2*i]==e[2*i+1]==0) for i in range(3)]
  mass=2*((bits&1)+((bits>>1)&1))
  electric=F(1,16)*(e[0]**2+e[1]**2+2*e[2]**2+2*e[3]**2+e[4]**2+e[5]**2)
  H[bits][bits]=F(1)+mass+electric+A*(F(neutral[0],2)+neutral[1]+F(neutral[2],2))-LOWER
  for src,dst,orientation in ((2,0,-1),(3,0,1),(3,1,-1),(4,1,1)):
   r=hop(bits,src,dst)
   if r is None:continue
   target,s=r
   H[bits][target]+=s*orientation;H[target][bits]+=s*orientation
 assert all(H[i][j]==H[j][i] for i in range(32) for j in range(32))
 return H

def exact_positive(H):
 # Positive pivot Schur elimination with Fraction, equivalent to rational LDL.
 for k in range(len(H)):
  pivot=H[k][k]
  assert pivot>0,("nonpositive pivot",k,pivot)
  for i in range(k+1,len(H)):
   aik=H[i][k]
   if aik:
    for j in range(i,len(H)):
     H[j][i]-=aik*H[j][k]/pivot
     if j!=i:H[i][j]=H[j][i]
  for i in range(k+1,len(H)):H[i][k]=H[k][i]=F(0)

if __name__=="__main__":
 for z in range(-5,6):
  for bits in range(32):
   assert all(abs(v-s*z)<=3 for v,s in zip(fields(bits,z),(1,-1,1,-1,1,-1)))
  exact_positive(block(z))
 # At |z|>=6 each of six E legs differs from +/-z by at most 3;
 # electric energy >=(1/2)(|z|-3)^2 >= 9/2; cost 1; hopping norm <=4.
 assert F(9,2)+1-4>LOWER
 print("EXACT_CP12_TWO_STAR_POSITIVE_LDL_PASS")
