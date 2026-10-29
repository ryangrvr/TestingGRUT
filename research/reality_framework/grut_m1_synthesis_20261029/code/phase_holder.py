"""Gauge-invariant curvature phase holder: dimensional and local-vacuum test.

All complex fields are z/H_i, kinetic L = |partial z|^2, and all potentials
are V/H_i^4. Positive Wilson coefficient yields a holder minimum at
arg(Sigma)-6 arg(X)=-pi/2 if it acts alone. This tests a local truncated EFT;
the operator by itself makes the large-field potential unbounded.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
p=json.loads((RF/'grut_matched_candidate_20261027'/'CANDIDATE_M1.json').read_text());heavy=p['heavy'];H=heavy['H_GeV'];cut=p['EFT_cutoff_GeV_declared'];ell=cut/H
lx,ly,ls=.001,1.,1.;M=heavy['M_Y_GeV']/H;h=heavy['h_GeV']/H;mu=heavy['mu_GeV']/H;k=heavy['kappa'];xiX=-5/6
z0=np.array([76.39978266076902,-2.111509854283515-0.9393290404896636j,-28.293142125868524j])

def pot_force(z,xi=0.,R=12.):
 x,y,s=z;c=R*xi/ell**5
 base=lx*abs(x)**4+ly*abs(y)**4+M*M*abs(y)**2+ls*(abs(s)**2-h*h/2)**2+2*np.real(k*np.conj(y)*x**3+mu*np.conj(s)*y*y)+xiX*R*abs(x)**2
 lock=2*np.real(1j*c*np.conj(s)*x**6)
 f=np.array([2*lx*abs(x)**2*x+3*k*y*np.conj(x)**2+xiX*R*x-6j*c*s*np.conj(x)**5,
             (M*M+2*ly*abs(y)**2)*y+k*x**3+2*mu*s*np.conj(y),
             2*ls*(abs(s)**2-h*h/2)*s+mu*y*y+1j*c*x**6])
 return base+lock,f,base,lock

def pack(z):return np.array([z[0].real,z[1].real,z[1].imag,z[2].real,z[2].imag])
def unpack(a):return np.array([a[0],a[1]+1j*a[2],a[3]+1j*a[4]])
def grad5(a,xi):
 _,f,_,_=pot_force(unpack(a),xi);return 2*np.array([f[0].real,f[1].real,f[1].imag,f[2].real,f[2].imag])
def wrap(x):return float((x+np.pi)%(2*np.pi)-np.pi)
def hessian6(z,xi):
 a=np.column_stack((z.real,z.imag)).ravel();D=np.zeros((6,6))
 def g(a):
  w=a[::2]+1j*a[1::2];_,f,_,_=pot_force(w,xi);return 2*np.column_stack((f.real,f.imag)).ravel()
 for j in range(6):
  hstep=2e-5*max(1,abs(a[j]));v=np.zeros(6);v[j]=hstep;D[:,j]=(g(a+v)-g(a-v))/(2*hstep)
 return (D+D.T)/4 # /2 for canonical real mass since kinetic is |dot z|²

def info(z,xi,res):
 V,f,b,l=pot_force(z,xi);eig=np.linalg.eigvalsh(hessian6(z,xi));return {'xi':xi,'fields_over_H':[[float(v.real),float(v.imag)] for v in z],
  'amplitudes_over_H':[float(abs(v)) for v in z],'delta_sigma_minus_6x':wrap(np.angle(z[2])-6*np.angle(z[0])),
  'delta_y_minus_3x':wrap(np.angle(z[1])-3*np.angle(z[0])),
  'force_norm_H3':float(np.linalg.norm(f)),'root_residual':float(res),'potential_over_H4':float(V),'base_over_H4':float(b),'holder_over_H4':float(l),
  'canonical_scalar_mass2_over_H2':eig.tolist(),'X_over_cutoff':float(abs(z[0])/ell),'Sigma_over_cutoff':float(abs(z[2])/ell)}

def main():
 f2=2/(1/abs(z0[2])**2+36/abs(z0[0])**2);m2unit=24*abs(z0[2])*abs(z0[0])**6/ell**5/f2;threshold=9/m2unit
 results=[];seed=np.array([80.2,-3.44,0.,-28.31,0.])
 # Continue from the aligned holder-free minimum in small fixed steps. The
 # truncated operator is unbounded at large field, so this tracks a local
 # stationary branch only, not a global de Sitter ground state.
 wanted=[0.,1e8,2e8,3e8,4e8,threshold,5e8,5.6e8]
 previous=0.
 for xi in wanted:
  points=np.linspace(previous,xi,max(2,int(np.ceil((xi-previous)/1e7))+1))
  for value in points:
   s=root(lambda x:grad5(x,value),seed,tol=1e-11)
   if np.linalg.norm(s.fun)>1e-6:raise RuntimeError(f'continuation failed xi={value}: {np.linalg.norm(s.fun)}')
   seed=s.x
  z=unpack(seed);item=info(z,xi,np.linalg.norm(s.fun));item['root_success']=bool(s.success);results.append(item);previous=xi
 initial=[info(z0,xi,float('nan')) for xi in (0.,threshold)]
 for a in initial:a['root_residual']=None
 # exact U(1) dimensions in four spacetime dimensions
 out={'operator_dimensions':{'Sigma_star_Y2':3,'R_Sigma_star_Y2':5,'correct_denominator_power_for_R_Sigma_star_Y2':1,
        'Sigma_star_X6':7,'R_Sigma_star_X6':9,'correct_denominator_power_for_R_Sigma_star_X6':5},
      'unit_charge_check':{'Sigma_star_Y2':-6+2*3,'Sigma_star_X6':-6+6},
      'prepared_state':initial,'holder_only_relative_phase_f2_over_H2':f2,'holder_only_mass2_per_xi':m2unit,'holder_only_xi_for_9H2':threshold,
      'inflation_H_GeV':H,'declared_cutoff_GeV':cut,'cutoff_over_H':ell,'mediator_over_H':M,'initial_X_over_mediator':float(abs(z0[0])/M),
      'initial_Sigma_over_mediator':float(abs(z0[2])/M),'effective_cutoff_at_xi_for_9_GeV':float(cut/threshold**.2),
      'states':results,'scope':'Local fixed-de-Sitter test with a new higher-curvature operator and original spectator potential. The trial operator is not a UV completion; no unique initial state, reheating trajectory, perturbation transfer or isocurvature bound follows from these roots.'}
 (ROOT/'phase_holder.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
