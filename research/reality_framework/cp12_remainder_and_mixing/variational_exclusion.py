"""Gauss-law star lower bound and a coherent incidence-matching trial state.
No coupling is varied. Decimal arithmetic gives outward-rounded Rayleigh
upper bounds and regular-sector lower bounds at Q=2N, frozen t=m=1.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import itertools,math,json,sys
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.linalg import eigh_tridiagonal
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,hop

def matching_counts(n):
 z=n-1;h=[1,z]
 for k in range(2,n+1):h.append(z*h[-1]-(k-1)*h[-2])
 return [math.comb(n,k)*h[k] for k in range(n+1)]

def squared_couplings(n,q):
 M=n*(n-1)//2;m=matching_counts(n);kmax=min(n,q)
 return [F(q-k,M-k)*(k+1)**2*F(m[k+1],m[k]) for k in range(kmax)]

def frac_dec(x,rounding):
 with localcontext() as c:
  c.prec=80;c.rounding=rounding
  return Decimal(x.numerator)/Decimal(x.denominator)

def sqrt_lower(x):
 with localcontext() as c:
  c.prec=80
  return frac_dec(x,ROUND_FLOOR).sqrt().next_minus()

def sqrt_upper(x):
 with localcontext() as c:
  c.prec=80
  return frac_dec(x,ROUND_CEILING).sqrt().next_plus()

def regular_lower(n,d=4):
 q=n*d//2;M=n*(n-1)//2;a=F(M-q,10*M);delta=F(2)+F(1,16)-a/2
 baseline=F(d,4)+a*d/2
 with localcontext() as c:
  c.prec=80;c.rounding=ROUND_FLOOR
  # Subtract an upper square root from a lower rational term.
  root=sqrt_upper(delta*delta+4*d)
  with localcontext() as up:
   up.prec=80;up.rounding=ROUND_CEILING;half=root/2
  return frac_dec(baseline+delta/2,ROUND_FLOOR)-half

def finite_certificate(n):
 d=4;q=2*n;M=n*(n-1)//2;c2=squared_couplings(n,q);center=n//4;width=math.isqrt(n)
 ids=list(range(max(0,center-width),min(n,center+width)+1));w=[width+1-abs(k-center) for k in ids];norm=sum(x*x for x in w)
 diag=F(q,2)+F(17,8)*F(sum(k*x*x for k,x in zip(ids,w)),norm)
 with localcontext() as c:
  c.prec=80;c.rounding=ROUND_FLOOR
  term=sum(Decimal(w[i]*w[i+1])*sqrt_lower(c2[ids[i]]) for i in range(len(ids)-1))
  subtraction=Decimal(2)*term/Decimal(norm)
  with localcontext() as up:
   up.prec=80;up.rounding=ROUND_CEILING
   upper=frac_dec(diag,ROUND_CEILING)-subtraction
  lower_density=regular_lower(n);gap=Decimal(n)*lower_density-upper
  b2=F((d+1)*q*(M-q),100*M*M)
  if gap>0:
   gap_square=gap*gap
   with localcontext() as up:
    up.prec=80;up.rounding=ROUND_CEILING
    probability=min(Decimal(1),frac_dec(b2,ROUND_CEILING)/gap_square)
  else:probability=Decimal(1)
  # Ritz optimum is a useful diagnostic, not the primary certificate.
  diag_float=np.arange(n+1)*17/8+q/2;off=np.sqrt([float(x) for x in c2]);ritz=float(eigh_tridiagonal(diag_float,off,select='i',select_range=(0,0),eigvals_only=True)[0])
  assert ritz<=float(upper)+1e-9
  with localcontext() as up:
   up.prec=80;up.rounding=ROUND_CEILING
   upper_density=upper/Decimal(n)
  return dict(N=n,Q=q,a=str(F(M-q,10*M)),trial_support=[ids[0],ids[-1]],trial_integer_weight_norm=norm,rayleigh_upper_total=str(upper),rayleigh_upper_density=str(upper_density),regular_lower_density=str(lower_density),regular_minus_trial_gap_lower=str(gap),regular_ground_probability_upper=str(probability),optimized_trial_ritz_density=ritz/n,scope='Analytic inequalities with outward Decimal arithmetic; no finite electric cutoff, cluster truncation or omitted graph mixing.')

def validate_small(n,q,L=2):
 H,states,lookup,degrees=build(n,q,L);M=n*(n-1)//2;m=matching_counts(n);kmax=min(n,q);cols=[];rows=[];data=[];level_counts=[0]*(kmax+1)
 sea=sum(1<<(n+e) for e in range(M));edges=list(itertools.combinations(range(n),2))
 for i,(mask,bits,E) in enumerate(states):
  if any(x not in (0,1) for x in E):continue
  k=(bits&((1<<n)-1)).bit_count()
  if k>kmax:continue
  # Nonnegative 0/1 fields plus Gauss's law force local incidence dimers.
  phase=1;created=sea
  for e,(a,b) in enumerate(edges):
   if E[2*e] or E[2*e+1]:
    leg=0 if E[2*e] else 1;vertex=a if leg==0 else b
    new,sg=hop(created,vertex,n+e);created=new;phase*=sg*(1 if leg==0 else -1)
  assert created==bits
  expected=m[k]*math.comb(M-k,q-k)
  rows.append(i);cols.append(k);data.append(phase/math.sqrt(expected));level_counts[k]+=1
 U=coo_matrix((data,(rows,cols)),shape=(len(states),kmax+1)).tocsr();assert all(level_counts[k]==m[k]*math.comb(M-k,q-k) for k in range(kmax+1))
 assert np.linalg.norm((U.T@U).toarray()-np.eye(kmax+1))<1e-10
 Hsmall=(U.T@H@U).toarray();predicted=np.diag(np.arange(kmax+1)*17/8+q/2);c=np.sqrt([float(x) for x in squared_couplings(n,q)])
 predicted+=np.diag(c,1)+np.diag(c,-1);error=float(np.linalg.norm(Hsmall-predicted));assert error<1e-10
 # The graph Laplacian annihilates every trial column.
 coo=H.tocoo();masks=np.array([s[0] for s in states]);sel=masks[coo.row]!=masks[coo.col]
 R=coo_matrix((coo.data[sel],(coo.row[sel],coo.col[sel])),shape=H.shape).tocsr()+diags(.1/M*degrees)
 err=float(np.linalg.norm((R@U).data));assert err<1e-10
 trial=float(np.linalg.eigvalsh(Hsmall)[0]);from scipy.sparse.linalg import eigsh
 ground=float(eigsh(H,k=1,which='SA',return_eigenvectors=False,tol=1e-12)[0]);assert ground<=trial+1e-9
 return dict(N=n,Q=q,L=L,full_dimension=len(states),matching_levels=level_counts,projected_matrix_error=error,rewiring_annihilation_error=err,trial_energy=trial,full_energy=ground)

if __name__=='__main__':
 small=[validate_small(3,2),validate_small(4,3),validate_small(4,4)]
 certificates=[finite_certificate(n) for n in [16,32,64,128,256,512,1024,2048]]
 out=dict(status='COHERENT_MATCHING_VARIATIONAL_CONTROL_PASS',small_hamiltonian_controls=small,certificates=certificates)
 (P/'variational_exclusion_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
