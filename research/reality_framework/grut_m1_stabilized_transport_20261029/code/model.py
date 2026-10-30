"""Dimension-five cubic curvature holder, fixed-curvature scalar sector.

Units: complex fields/H_i, time H_i t, curvature R/H_i², potential/H_i⁴.
Kinetic term is sum |dz|². Gravity is prescribed in the spectator evolution.
"""
from pathlib import Path
import json,numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
P=json.loads((RF/'grut_matched_candidate_20261027/CANDIDATE_M1.json').read_text())
H=P['heavy']['H_GeV'];M=P['heavy']['M_Y_GeV']/H;h=P['heavy']['h_GeV']/H;mu=P['heavy']['mu_GeV']/H;k=P['heavy']['kappa']
lx=P['heavy']['lambda_X'];ly=P['heavy']['lambda_Y'];ls=P['heavy']['lambda_Sigma'];xi=P['reheating']['xiX'];e=P['heavy']['e'];eta=H/2.435e18
LOCK=(10j-mu)/12
T=np.zeros((6,6))
for i,q in enumerate([1,3,6]):T[2*i,2*i+1]=-q;T[2*i+1,2*i]=q

def real(z):return np.sqrt(2)*np.column_stack((z.real,z.imag)).ravel()
def complexvec(r):return (r[::2]+1j*r[1::2])/np.sqrt(2)
def evaluate(z,R=12.,lock=LOCK):
 x,y,s=z;u=mu+lock*R
 V=lx*abs(x)**4+ly*abs(y)**4+M*M*abs(y)**2+ls*(abs(s)**2-h*h/2)**2+2*np.real(k*np.conj(y)*x**3+u*np.conj(s)*y*y)+xi*R*abs(x)**2
 F=np.array([2*lx*abs(x)**2*x+3*k*y*np.conj(x)**2+xi*R*x,
  (M*M+2*ly*abs(y)**2)*y+k*x**3+2*np.conj(u)*s*np.conj(y),
  2*ls*(abs(s)**2-h*h/2)*s+u*y*y])
 VR=xi*abs(x)**2+2*np.real(lock*np.conj(s)*y*y)
 return V,F,VR

def hessian(z,R=12.,lock=LOCK):
 x,y,s=z;u=mu+lock*R;A=np.zeros((3,3),complex);B=A.copy()
 A[0,0]=4*lx*abs(x)**2+xi*R;B[0,0]=2*lx*x*x+6*k*y*np.conj(x)
 A[0,1]=3*k*np.conj(x)**2;A[1,0]=np.conj(A[0,1])
 A[1,1]=4*ly*abs(y)**2+M*M;B[1,1]=2*ly*y*y+2*np.conj(u)*s
 A[1,2]=2*np.conj(u)*np.conj(y);A[2,1]=np.conj(A[1,2])
 A[2,2]=4*ls*abs(s)**2-ls*h*h;B[2,2]=2*ls*s*s
 D=np.empty((6,6));D[::2,::2]=(A+B).real;D[::2,1::2]=-(A-B).imag;D[1::2,::2]=(A+B).imag;D[1::2,1::2]=(A-B).real
 return D

def initial():
 u=abs(mu+12*LOCK)
 def eq(a):
  x,y,s=a
  return [2*lx*x**3-3*k*y*x*x+12*xi*x,(M*M+2*ly*y*y)*y-k*x**3-2*u*s*y,2*ls*(s*s-h*h/2)*s-u*y*y]
 r=root(eq,[88,6,28.4],tol=1e-12)
 if np.linalg.norm(r.fun)>1e-7:raise RuntimeError(r.message)
 x,y,s=r.x;return np.array([x,-y,s*np.exp(1j*(np.pi+np.angle(mu+12*LOCK)))])

def lower_bound(R,lock=LOCK,reserve=.2):
 A=reserve*lx;B=ly-27*k**4/(16*((1-reserve)*lx)**3);u=abs(mu+lock*R);ax=max(-xi*R,0)
 assert B>0
 return ls*h**4/4-ax**2/(4*A)-(ls*h*h+u*u/B)**2/(4*ls),A,B

def main():
 z=initial();V,F,VR=evaluate(z);lo,A,B=lower_bound(12);vals=np.linalg.eigvalsh(hessian(z));out={'parent_H_GeV':H,'lock_complex_over_H':[LOCK.real,LOCK.imag],
 'mu_eff_initial_over_H':[float((mu+12*LOCK).real),float((mu+12*LOCK).imag)],'xi_muR':float(abs(LOCK)*M),'alpha_R':float(np.angle(LOCK)),
 'initial_fields_over_H':[[float(x.real),float(x.imag)] for x in z],'initial_amplitudes_over_H':abs(z).tolist(),'initial_phase_sigma_minus_6x':float((np.angle(z[2])-6*np.angle(z[0])+np.pi)%(2*np.pi)-np.pi),
 'initial_force_norm':float(np.linalg.norm(F)),'initial_potential_over_H4':float(V),'canonical_mass2_over_H2':vals.tolist(),'quartic_reserved_X':A,'quartic_reduced_Y':B,
 'fixed_R12_global_lower_bound_over_H4':lo,'initial_effective_Planck_ratio':float(1-2*eta**2*VR),
 'scope':'Fixed-curvature scalar boundedness with dimension-five curvature cubic. Desired initial phase encoded in chosen complex Wilson coefficient; no UV completion, unique state law, or global gravitational positivity.'}
 (ROOT/'holder_model.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
