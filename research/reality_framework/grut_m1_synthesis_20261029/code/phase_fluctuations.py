"""Illustrative linear de-Sitter phase power for the local holder minimum.

Canonical constant-H, constant-Hessian, independent Bunch-Davies scalar modes.
No gauge-field, metric, curved field-space, reheating or isocurvature transfer.
All of those omissions preclude a physical isocurvature verdict.
"""
from pathlib import Path
import json, numpy as np
from scipy.integrate import solve_ivp
import phase_holder as ph
ROOT=Path(__file__).resolve().parents[1]

def modes(m2,x0,points):
 n=len(m2);N0=-np.log(x0);b=m2-2.;w=np.sqrt(1+b/x0**2);dw=-b/(w*x0**3)
 y=1/np.sqrt(w)+0j;yx=(1j*w-.5*dw/w)*y
 S=x0*y;dS=-x0*(y+x0*yx);initial=np.r_[S,dS]
 def rhs(N,Y):return np.r_[Y[n:],-3*Y[n:]-(np.exp(-2*N)+m2)*Y[:n]]
 sol=solve_ivp(rhs,(N0,max(points)),initial,method='DOP853',rtol=2e-9,atol=2e-11,t_eval=points)
 if not sol.success:raise RuntimeError(sol.message)
 return sol.y[:n,:],sol.y[n:,:]

def main():
 j=json.loads((ROOT/'phase_holder.json').read_text());st=next(x for x in j['states'] if x['xi']==4e8);z=np.array([complex(*x) for x in st['fields_over_H']]);vals,U=np.linalg.eigh(ph.hessian6(z,4e8))
 grad=np.zeros(6);grad[0]=6*z[0].imag/(np.sqrt(2)*abs(z[0])**2);grad[1]=-6*z[0].real/(np.sqrt(2)*abs(z[0])**2)
 grad[4]=-z[2].imag/(np.sqrt(2)*abs(z[2])**2);grad[5]=z[2].real/(np.sqrt(2)*abs(z[2])**2)
 proj=grad@U;keep=vals>1e-4;m2=vals[keep];proj=proj[keep];assert min(m2)>2.25
 points=np.array([0.,1.,2.,4.,8.]);S,D=modes(m2,1000.,points);S2,D2=modes(m2,2000.,points)
 rows=[]
 for i,N in enumerate(points):rows.append({'N_after_horizon_exit':float(N),'phase_sqrt_power':float(np.sqrt(np.sum(proj**2*abs(S[:,i])**2))/(2*np.pi)),
     'start_depth_1000_vs_2000_fraction':float(np.sqrt(np.sum(proj**2*abs(S[:,i])**2)/np.sum(proj**2*abs(S2[:,i])**2))-1)})
 for N in (20.,50.):
  d=N-8.;mu=np.sqrt(m2-2.25);late=np.exp(-1.5*d)*(S[:,-1]*np.cos(mu*d)+(D[:,-1]+1.5*S[:,-1])/mu*np.sin(mu*d))
  rows.append({'N_after_horizon_exit':N,'phase_sqrt_power':float(np.sqrt(np.sum(proj**2*abs(late)**2))/(2*np.pi)),'start_depth_1000_vs_2000_fraction':None})
 test,_=modes(np.array([0.]),1000.,np.array([0.,4.]));expected=np.sqrt(1+np.exp(-2*np.array([0.,4.])))
 out={'xi':4e8,'mass2_over_H2':m2.tolist(),'relative_phase_projection_squared':(proj**2).tolist(),
      'goldstone_projection_squared':float((grad@U[:,~keep])@(grad@U[:,~keep])),
      'massless_analytic_relative_error':float(np.max(abs(abs(test[0])-expected)/expected)),
      'assumed_massless_horizon_exit_sqrt_power':float(np.linalg.norm(grad)*np.sqrt(2)/(2*np.pi)),
      'benchmark_phase_sqrt_power_limit_old_transfer':4.140502195878647e-5,'rows':rows,
      'scope':'Fixed de Sitter and local constant scalar Hessian, independent Bunch-Davies test fields. N is the assumed duration after a mode exits. Gauge/metric perturbations, sources and mode mixing during the real transition, state choice, and the new holder-to-charge transfer are not computed. The old Planck-inspired limit is not automatically applicable.'}
 (ROOT/'phase_fluctuations.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
