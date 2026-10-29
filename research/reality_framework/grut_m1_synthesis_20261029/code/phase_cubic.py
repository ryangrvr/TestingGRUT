"""Second holder: local de Sitter minimum for R Sigma* Y² / M_Y.

At fixed curvature this replaces mu by mu+i alpha. For kappa,mu>0,
phase minimization gives theta_Y=pi and theta_Sigma=-pi+arg(mu+i alpha)
when X is real and positive. This is a local stationary-point test.
"""
from pathlib import Path
import json,math
import numpy as np
from scipy.optimize import root
import phase_holder as ph
ROOT=Path(__file__).resolve().parents[1]

def equations(v,mu_eff):
 x,y,s=v
 return np.array([2*ph.lx*x**3-3*ph.k*y*x*x-10*x,
  (ph.M**2+2*ph.ly*y*y)*y-ph.k*x**3-2*mu_eff*s*y,
  2*ph.ls*(s*s-ph.h**2/2)*s-mu_eff*y*y])
def full_force(z,mu_complex):
 x,y,s=z
 return np.array([2*ph.lx*abs(x)**2*x+3*ph.k*y*np.conj(x)**2-10*x,
  (ph.M**2+2*ph.ly*abs(y)**2)*y+ph.k*x**3+2*np.conj(mu_complex)*s*np.conj(y),
  2*ph.ls*(abs(s)**2-ph.h**2/2)*s+mu_complex*y*y])
def mass_eigenvalues(z,mu_complex):
 a=np.column_stack((z.real,z.imag)).ravel();D=np.zeros((6,6))
 def grad(a):
  f=full_force(a[::2]+1j*a[1::2],mu_complex);return 2*np.column_stack((f.real,f.imag)).ravel()
 for i in range(6):
  d=np.zeros(6);d[i]=2e-5*max(1,abs(a[i]));D[:,i]=(grad(a+d)-grad(a-d))/(2*d[i])
 return np.linalg.eigvalsh((D+D.T)/4)
def main():
 seed=np.array([80.2020754,3.44245486,28.31074753]);rows=[];last=0
 for xi in (0.,1.,5.,10.,20.,30.,50.,100.,200.,500.,2000.):
  for xii in np.linspace(last,xi,max(2,int(np.ceil((xi-last)/2))+1)):
   alpha=12*xii/ph.M;mu_eff=np.hypot(ph.mu,alpha);r=root(lambda v:equations(v,mu_eff),seed,tol=1e-11)
   if np.linalg.norm(r.fun)>1e-6:raise RuntimeError(f'root failure {xii} {r.fun}')
   seed=r.x
  alpha=12*xi/ph.M;mu_c=ph.mu+1j*alpha;x,y,s=seed;angle=ph.wrap(np.pi+np.angle(mu_c));z=np.array([x,-y,s*np.exp(1j*angle)]);ev=mass_eigenvalues(z,mu_c)
  rows.append({'xi':xi,'alpha_over_H':alpha,'abs_mu_eff_over_H':abs(mu_c),'delta_sigma_minus_6x':angle,
      'deviation_from_minus_pi_over_2':ph.wrap(angle+np.pi/2),'amplitudes_over_H':[float(x),float(y),float(s)],
      'force_norm_over_H3':float(np.linalg.norm(full_force(z,mu_c))),
      'smallest_eigenvalues_m2_over_H2':ev[:3].tolist(),'fixed_sigma_Y_quadratic_min_over_H2':float(ph.M**2-2*abs(mu_c)*s),
      'X_over_M_Y':float(x/ph.M),'Y_over_M_Y':float(y/ph.M)})
  last=xi
 s0=abs(ph.z0[2]);crit=math.sqrt(max((ph.M**2/(2*s0))**2-ph.mu**2,0))*ph.M/12
 out={'rows':rows,'xi_at_zero_quadratic_Y_mass_if_sigma_amplitude_frozen':crit,
      'scope':'All rows are fixed-de-Sitter stationary points in a spectator scalar potential with the supplied curvature operator, not a UV completion or a new cosmological prediction. The quadratic Y criterion freezes Sigma and is not a full stability test. X/MY>1 invalidates an automatic expansion in the proposed mediator scale.'}
 (ROOT/'phase_cubic.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
