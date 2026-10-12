"""Controlled polynomial resolvent remainder in the full-geometry Schur problem.
Small complete finite-host benchmark; not a large-N phase certificate.
Dense complement spectra are used for this validation. On a larger system,
independent spectral enclosures and matrix-polynomial applications are needed.
"""
from pathlib import Path
import json,math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator,eigsh
from mixing_fibers import run
P=Path(__file__).resolve().parent

def cheb_inverse(vals,z,lo,hi,K):
 a=lo-z;b=hi-z;assert a>0
 d=(hi+lo)/2;h=(hi-lo)/2;root=math.sqrt(a*b)
 r=(math.sqrt(b)-math.sqrt(a))/(math.sqrt(b)+math.sqrt(a))
 x=(vals-d)/h;T0=np.ones_like(x);T1=x.copy();p=T0/root
 if K>=1:p+=2*(-r)*T1/root
 for j in range(2,K+1):
  T2=2*x*T1-T0;p+=2*(-r)**j*T2/root;T0,T1=T1,T2
 tail=2/root*r**(K+1)/(1-r)
 assert max(abs(p-1/(vals-z)))<=tail+2e-13
 return p,tail,r

def main(n=4,q=3,L=3):
 H,U,Q,A,info=run(n,q,L,False,True)
 C=(Q.T@H@Q).toarray();B=(Q.T@H@U).toarray();A=A.toarray();c,W=eigh(C)
 # Residual-padding records floating-point numerical controls. Not an interval
 # arithmetic proof of an infinite-cutoff physical spectrum.
 orth=np.linalg.norm(W.T@W-np.eye(len(c)),ord='fro');res=np.linalg.norm(C-W@np.diag(c)@W.T,ord='fro')
 pad=res+max(abs(c))*orth+1e-10;lo=float(c[0]-pad);hi=float(c[-1]+pad);B=W.T@B;b2=float(np.linalg.norm(B,ord='fro')**2)
 upper=info['uniform_fiber_variational_energy'];lower=info['best_graph_block_energy']-info['offdiag_norm_exact']-.01
 assert upper<lo
 def smallest(z,K=None):
  if K is None:g=1/(c-z)
  else:g,_,_=cheb_inverse(c,z,lo,hi,K)
  F=A-B.T@(g[:,None]*B)-z*np.eye(A.shape[0])
  return float(eigh(F,subset_by_index=[0,0],eigvals_only=True)[0])
 assert smallest(lower)>0 and smallest(upper)<0
 exact_root=brentq(lambda z:smallest(z),lower,upper,xtol=1e-13)
 assert abs(exact_root-info['full_energy'])<1e-10
 rows=[]
 for K in [16,32,64,128,256]:
  # For small K the polynomial may not bracket its approximate root; retain
  # the complete resolvent's root as a benchmark evaluation point instead.
  fl=smallest(lower,K);fu=smallest(upper,K)
  if fl*fu<0:point=brentq(lambda z:smallest(z,K),lower,upper,xtol=1e-13)
  else:point=exact_root
  _,tail,r=cheb_inverse(c,point,lo,hi,K);defect=abs(smallest(point,K));eta=b2*tail+defect+1e-10
  low=point-eta;high=point+eta;valid=high<lo
  if valid:
   # F is strictly decreasing at speed at least1. The controlled uniform
   # approximation error therefore transfers into this energy enclosure.
   assert low<=info['full_energy']<=high
  rows.append(dict(degree=K,point=point,ratio=r,resolvent_uniform_tail_bound=tail,schur_norm_error_bound=b2*tail,energy_halfwidth=eta,energy_lower=low,energy_upper=high,below_complement_spectrum=bool(valid),actual_energy_error=abs(point-info['full_energy'])))
 assert rows[-1]['below_complement_spectrum'] and rows[-1]['energy_halfwidth']<1e-8
 out=dict(status='SCHUR_POLYNOMIAL_REMAINDER_CONTROL_PASS',N=n,Q=q,L=L,full_dimension=H.shape[0],kernel_dimension=U.shape[1],complement_dimension=Q.shape[1],full_energy=info['full_energy'],exact_schur_root=exact_root,complement_spectral_bracket=[lo,hi],spectral_padding=pad,complement_eigenbasis_residual=res,coupling_frobenius_square=b2,rows=rows,scope='Conditional analytic resolvent bound with floating-point finite-cutoff validation; no certified large-N spectral enclosure or graph-cluster tail.')
 (P/f'schur_N{n}_Q{q}_L{L}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
