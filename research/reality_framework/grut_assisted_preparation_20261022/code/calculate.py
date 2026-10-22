"""Same-action normalization and a bounded fixed-inflaton stochastic probe.

No fitted density, complete reheating model or general empirical exclusion.
"""
from pathlib import Path
import json,math,platform,sys,time
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.special import hankel1
P=Path(__file__).resolve().parents[1]
I=dict(base_commit='4f4b4b04f67294213be15ec3cb51ff0a6c90243d',M_GeV=2.435e18,v_GeV=246.,lam=.129,As=2.1e-9,ns_snapshot=.9649,ns_snapshot_sigma=.0042,r_snapshot_upper=.036,rho_reference_GeV4=2.52e-47,g=1.5,c_values=[.5,1.,2.5,3.],N_values=[50.,55.,60.],probe_c=10.,probe_b=3.,probe_calibration_N=55.,probe_times=[10.,100.],grid_configs=[[.2,80.],[.1,80.],[.05,80.],[.1,100.]],source='1905.00045v2 literal PDF eq 5.2, HTML eq 20',exposure='exploratory; all observational snapshots already known; no held-out comparison')
def dump(name,x):(P/name).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def rho_relation(ns,r,As=I['As']):
 V0=I['lam']*I['v_GeV']**4/4
 return 64*V0*V0/(3*I['g']*math.pi**2*As*I['M_GeV']**4*(8*(1-ns)-r)**2)
def calibration(c,N):
 k=I['g']*c*c;xe=math.log1p(math.sqrt(2/k))
 x=brentq(lambda x:k/2*(math.exp(x)-x-math.exp(xe)+xe)-N,xe,20)
 z=math.exp(x);epsilon=2/(k*(z-1)**2);eta=2*(2-z)/(k*(z-1)**2)
 ns=1-6*epsilon+2*eta;r=16*epsilon;f=(-math.expm1(-x))**2
 alpha2=24*math.pi**2*I['As']*epsilon/f
 V0=I['lam']*I['v_GeV']**4/4;M=I['M_GeV']
 rho=c*c*V0*V0/(alpha2*M**4);pred=rho_relation(ns,r)
 # Canonical finite derivatives are independent of the closed slow-roll form.
 q=math.sqrt(k)*x;h=1e-3
 def fun(q):return math.expm1(-q/math.sqrt(k))**2
 f0=fun(q);fp=(fun(q-2*h)-8*fun(q-h)+8*fun(q+h)-fun(q+2*h))/(12*h)
 fpp=(-fun(q+2*h)+16*fun(q+h)-30*f0+16*fun(q-h)-fun(q-2*h))/(12*h*h)
 efd=.5*(fp/f0)**2;etafd=fpp/f0
 return dict(c=c,N=N,x_pivot=x,epsilon=epsilon,eta=eta,ns=ns,r=r,alpha=math.sqrt(alpha2),rho_direct_GeV4=rho,rho_eliminated_GeV4=pred,relative_elimination_error=abs(pred/rho-1),relative_epsilon_fd_error=abs(efd/epsilon-1),eta_fd_error=abs(etafd-eta))
def bernoulli(x):
 out=np.ones_like(x);mask=abs(x)>1e-7
 out[mask]=x[mask]/np.expm1(x[mask]);out[~mask]=1-x[~mask]/2+x[~mask]**2/12
 return out
def fp(dx,L,drift,times,name):
 x=np.arange(-round(L/dx),round(L/dx)+1)*dx;phi=drift*np.maximum(x,0)**2
 d=np.diff(phi);f=.5/dx**2*bernoulli(d);b=.5/dx**2*bernoulli(-d)
 diagonal=np.zeros(len(x));diagonal[:-1]-=f;diagonal[1:]-=b
 Q=diags([f,diagonal,b],[-1,0,1],format='csc')
 p0=np.zeros(len(x));p0[len(x)//2]=1.;vectors=[];rows=[]
 for t in times:
  prob=expm_multiply(t*Q,p0,traceA=t*np.sum(diagonal));vectors.append(prob)
  pos=np.maximum(x,0)
  rows.append(dict(N=t,norm=float(sum(prob)),minimum=float(min(prob)),edge_mass=float(sum(prob[:5])+sum(prob[-5:])),positive_probability=float(sum(prob[x>0])+.5*prob[len(x)//2]),positive_first_moment=float(np.dot(pos,prob)),positive_second_moment=float(np.dot(pos*pos,prob))))
 np.savetxt(P/(name+'.csv'),np.column_stack([x]+vectors),delimiter=',',header='canonical_over_Hdiv2pi,'+','.join('prob_N%g'%t for t in times),comments='')
 return dict(dx=dx,L=L,drift=drift,points=len(x),rows=rows)
def mode(c,xinitial):
 ratio=math.sqrt(6/I['g'])/c;nu=math.sqrt(2.25-ratio*ratio)
 targets=np.array([.1,.01]);times=np.log(xinitial/targets)
 v0=np.exp(1j*xinitial)/math.sqrt(2)
 def rhs(t,y):
  x=xinitial*np.exp(-t);v,vp=y
  return [vp,-vp-(x*x-2+ratio*ratio)*v]
 sol=solve_ivp(rhs,(0,times[-1]),[v0,-1j*xinitial*v0],t_eval=times,method='DOP853',rtol=1e-10,atol=1e-12)
 if not sol.success:raise RuntimeError(sol.message)
 rows=[]
 for j,x in enumerate(targets):
  v,vp=sol.y[:,j];direct=x*x*abs(v)**2/(2*math.pi**2);exact=x**3*abs(hankel1(nu,x))**2/(8*math.pi)
  w=(v*np.conj(vp)-np.conj(v)*vp)/(1j*x)
  rows.append(dict(x=float(x),power_over_H2_ode=float(direct),power_over_H2_hankel=float(exact),relative_error=float(abs(direct/exact-1)),wronskian_relative=float(abs(w-1))))
 return dict(c=c,mass_over_H=ratio,nu=nu,xinitial=xinitial,nfev=sol.nfev,rows=rows)
def main():
 start=time.perf_counter();dump('INPUTS.json',I);checks=[]
 def check(name,ok,value):
  checks.append(dict(name=name,passed=bool(ok),value=value))
  if not ok:raise AssertionError((name,value))
 cal=[calibration(c,N) for c in I['c_values'] for N in I['N_values']]
 check('literal canonical normalization',abs(.5*I['g']-.75)<1e-15,I['g'])
 check('same-action elimination identity',max(q['relative_elimination_error'] for q in cal)<1e-12,max(q['relative_elimination_error'] for q in cal))
 check('potential finite-derivative slow-roll',max(q['relative_epsilon_fd_error'] for q in cal)<1e-6 and max(q['eta_fd_error'] for q in cal)<1e-7,dict(epsilon=max(q['relative_epsilon_fd_error'] for q in cal),eta=max(q['eta_fd_error'] for q in cal)))
 np.savetxt(P/'calibration.csv',[[q[k] for k in ['c','N','ns','r','alpha','rho_direct_GeV4','rho_eliminated_GeV4']] for q in cal],delimiter=',',header='c,N,ns,r,alpha,rho_direct_GeV4,rho_eliminated_GeV4',comments='')
 ns=I['ns_snapshot'];rr=I['r_snapshot_upper']
 corners=[dict(ns=n,r=r,rho_GeV4=rho_relation(n,r)) for n in [ns-2*I['ns_snapshot_sigma'],ns+2*I['ns_snapshot_sigma']] for r in [0.,rr]]
 central=[dict(r=r,rho_GeV4=rho_relation(ns,r)) for r in [0.,rr]]
 V0=I['lam']*I['v_GeV']**4/4;M=I['M_GeV']
 root=math.sqrt(64*V0*V0/(3*I['g']*math.pi**2*I['As']*M**4*I['rho_reference_GeV4']))
 snapshot=dict(central=central,rectangular_sensitivity_corners=corners,scope='dated parameter box, not a combined confidence region or likelihood; no held-out data',ns_required_for_reference=[dict(r=r,ns=1-(r+root)/8) for r in [0.,rr]],max_ratio_to_reference=max(q['rho_GeV4'] for q in corners)/I['rho_reference_GeV4'])
 cp=calibration(I['probe_c'],I['probe_calibration_N']);c=I['probe_c'];k=I['g']*c*c;alpha=cp['alpha'];S=alpha*M/(2*math.pi*math.sqrt(3));s0=c*V0/(alpha*alpha*M**4)
 xearly=math.log(I['probe_b']**2/(k*s0*s0))+20;pole=I['probe_b']*math.exp(-xearly/2)
 drift=2/(k);H=alpha*M/math.sqrt(3)
 probe=dict(c=c,b=I['probe_b'],alpha=alpha,H_GeV=H,kappa=k,sigma0_over_M=s0,fixed_inflaton_over_M=xearly,pole_amplitude=pole,pole_metric_ratio_at_minimum=pole*pole/(k*s0*s0),minimum_canonical_scale_over_noise=math.sqrt(k)*s0*M/S,pole_canonical_scale_over_noise=pole*M/S,mass_over_H=math.sqrt(6/k),drift=drift,calibration_r=cp['r'],calibration_ns=cp['ns'],scope='fixed early plateau/light-field stochastic limit; c=10 control is not a claimed CMB-compatible point')
 runs=[]
 for dx,L in I['grid_configs']:
  name='fp_dx%g_L%g'%(dx,L);q=fp(dx,L,drift,I['probe_times'],name);runs.append(q);print('DONE',name,flush=True)
  for row in q['rows']:
   check(name+' N%g probability'%row['N'],abs(row['norm']-1)<1e-8 and row['minimum']>=-1e-12,row)
   check(name+' N%g boundaries'%row['N'],row['edge_mass']<1e-8,row['edge_mass'])
 for j,t in enumerate(I['probe_times']):
  base=runs[1]['rows'][j];fine=runs[2]['rows'][j];large=runs[3]['rows'][j]
  for label,other,tol in [('refinement',fine,.005),('volume',large,1e-5)]:
   errors={key:abs(base[key]-other[key])/max(abs(other[key]),1e-300) for key in ['positive_probability','positive_second_moment']}
   check(label+' N%g'%t,max(errors.values())<tol,errors)
 brown=fp(.1,80,0.,I['probe_times'],'brownian_control')
 for row in brown['rows']:check('Brownian N%g'%row['N'],abs(row['positive_probability']-.5)<1e-10 and abs(row['positive_second_moment']/(row['N']/2)-1)<1e-8,row)
 moments=[]
 for row in runs[2]['rows']:
  sigma2=(S/M)**2/k*row['positive_second_moment']
  moments.append(dict(N=row['N'],positive_probability=row['positive_probability'],sigma_rms_over_M=math.sqrt(sigma2),sigma_rms_over_selected_sigma0=math.sqrt(sigma2)/s0,coarse_grained_sigma_energy_over_inflation_energy=sigma2,scope='only this fixed-inflaton asymptotic stochastic probe; not post-reheating relic energy'))
 modes=[mode(c,initial) for c in [2.5,10.] for initial in [1000.,2000.]]
 check('mode ODE vs Bunch-Davies Hankel',max(row['relative_error'] for q in modes for row in q['rows'])<2e-5,max(row['relative_error'] for q in modes for row in q['rows']))
 check('mode Wronskian',max(row['wronskian_relative'] for q in modes for row in q['rows'])<1e-4,max(row['wronskian_relative'] for q in modes for row in q['rows']))
 finite=[]
 for eps in [.01,.1]:
  for t in [1.,10.]:
   nu=math.sqrt(2.25-probe['mass_over_H']**2)
   val=quad(lambda ell:(eps*math.exp(-ell))**3*abs(hankel1(nu,eps*math.exp(-ell)))**2/(8*math.pi),0,t,epsabs=0,epsrel=1e-10)[0]
   finite.append(dict(epsilon=eps,efold_band=t,linear_canonical_variance_over_H2=val,linear_rms_over_selected_sigma0=math.sqrt(val/k)*H/(s0*M),scope='local linear diagnostic; not a Gaussian continuation across the pole'))
 dump('RESULTS.json',dict(calibration=cal,snapshot=snapshot,probe=probe,fp_runs=runs,brownian_control=brown,probe_moments=moments,mode_controls=modes,linear_finite_band=finite,checks=checks,status='PASS' if all(q['passed'] for q in checks) else 'FAIL'))
 dump('ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=__import__('scipy').__version__,platform=platform.platform(),seconds=time.perf_counter()-start,command='python3 code/calculate.py'))
 print('PASS',len(checks),'internal checks',flush=True)
if __name__=='__main__':main()
