"""Isotropic NR elastic XX evolution in comoving kinetic energy.

No physical-momentum redshift remapping: epsilon=(a/a_i)^2 p^2/(2m²).
Constant isotropic threshold cross section, classical dilute gas, no sources,
core scattering, annihilation, Bose factors or perturbations. Full angular
redistribution conserves number and comoving NR kinetic energy.
"""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp,cumulative_trapezoid,quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq,root
from scipy.special import gammainc,gammaln,gamma
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent;sys.path.insert(0,str(RF/'grut_m1_synthesis_20261029/code'))
import free_stream as fs
AMP=4*.001-(.2**2)/3;SIGMA=AMP*AMP/(128*np.pi)

def project(grid,values,weights):
 j=np.clip(np.searchsorted(grid,values)-1,0,len(grid)-2);f=(values-grid[j])/(grid[j+1]-grid[j]);out=np.zeros(len(grid));np.add.at(out,j,weights*(1-f));np.add.at(out,j+1,weights*f);return out

def maxwell(grid,temp):
 a=grid[:-1];b=grid[1:];mass=gammainc(1.5,b/temp)-gammainc(1.5,a/temp);eng=1.5*temp*(gammainc(2.5,b/temp)-gammainc(2.5,a/temp));w=np.zeros(len(grid));w[:-1]+=(b*mass-eng)/(b-a);w[1:]+=(eng-a*mass)/(b-a);return w/w.sum()

def operator(grid,angles):
 mu,wm=leggauss(angles);ii,jj=np.triu_indices(len(grid));D=np.zeros((len(grid),len(ii)));K=np.zeros(len(ii));omit=K.copy();tic=time.monotonic()
 for c,(i,j) in enumerate(zip(ii,jj)):
  total=grid[i]+grid[j];rootprod=np.sqrt(grid[i]*grid[j]);vrel=np.sqrt(np.maximum(2*total-4*rootprod*mu,0));weight=wm/2*vrel;K[c]=sum(weight)
  if total==0:continue
  spread=.5*np.sqrt(np.maximum(total*total-4*grid[i]*grid[j]*mu*mu,0));lo=total/2-spread;hi=total/2+spread
  lower=np.maximum(lo,total-grid[-1]);upper=np.minimum(hi,grid[-1]);frac=np.maximum(upper-lower,0)/(hi-lo);rate=weight@frac;omit[c]=K[c]-rate
  left=np.maximum(lower[:,None],grid[None,:-1]);right=np.minimum(upper[:,None],grid[None,1:]);width=np.maximum(right-left,0);middle=(left+right)/2
  v=weight[:,None]*width/(hi-lo)[:,None]/np.diff(grid)[None,:]
  D[:-1,c]+=2*np.sum(v*(grid[None,1:]-middle),axis=0);D[1:,c]+=2*np.sum(v*(middle-grid[None,:-1]),axis=0);D[i,c]-=rate;D[j,c]-=rate
 return ii,jj,D,K,omit,{'build_seconds':time.monotonic()-tic,'number_column_error':float(max(abs(D.sum(axis=0)))),'energy_column_error':float(max(abs(grid@D)))}

def clock(N):
 _,_,ai,_,_=fs.background(N=50);aa,E,_,_,bg=fs.background(N=1201,a_end=10*ai);x=aa/ai
 HgeV=fs.H0/3.0856775814913673e19*6.582119569e-25*E
 point=json.loads((RF/'grut_matched_candidate_20261027/CANDIDATE_M1.json').read_text());m=point['light_renormalized_targets']['mX_GeV'];nbar=json.loads((RF/'grut_m1_core_matching_20261029/SCENARIO_AUDIT.json').read_text())['ncore_over_mX3']
 factor=SIGMA*N*nbar*m/1000/HgeV/x**4
 exposure=cumulative_trapezoid(factor,np.log(x),initial=0)
 temps=np.geomspace(.05,2.,250);gs=np.array([2+fs.e_dof(t)[1] for t in temps]);xt=2/temps*(gs[-1]/gs)**(1/3);T=PchipInterpolator(np.log(xt[::-1]),np.log(temps[::-1]))(np.log(x))
 return x,np.exp(T),factor,exposure,{'a_initial':ai,'H_initial_GeV':float(HgeV[0]),'sigma_threshold_mX_minus2':SIGMA,'amplitude_threshold':AMP,'density_inherited_a0_10':True}

def fit(grid,w,xend):
 p=np.sqrt(2*grid);mean2=w@(p*p);mean4=w@(p**4);ratio=mean4/mean2**2
 beta=brentq(lambda b:np.exp(gammaln(7/b)+gammaln(3/b)-2*gammaln(5/b))-ratio,.15,30)
 p0=np.sqrt(mean2*np.exp(gammaln(3/beta)-gammaln(5/beta)))
 speed2=(p/xend)**2/(1+(p/xend)**2);v2=w@speed2;v4=w@(speed2**2)
 def moments(vals):
  scale,b=np.exp(vals);shape=3/b
  def fun(y,power):
   if y==0:return 0.
   pp=scale*y**(1/b)/xend;vv=pp*pp/(1+pp*pp)
   return np.exp((shape-1)*np.log(y)-y-gammaln(shape))*vv**power
  return np.array([quad(lambda y:fun(y,power),0,150,epsabs=1e-16,epsrel=1e-9)[0] for power in [1,2]])
 rr=root(lambda vals:np.log(moments(vals)/[v2,v4]),np.log([p0,beta]),tol=1e-9)
 if not rr.success and np.linalg.norm(rr.fun)>1e-7:raise RuntimeError('fit failure '+rr.message)
 p0,beta=np.exp(rr.x);cdf=gammainc(3/beta,(p/p0)**beta);emp=np.cumsum(w)-.5*w
 q=np.geomspace(1e-5,25,2400);q0=p0/.05;f=beta/(q0**3*gamma(3/beta))*np.exp(-(q/q0)**beta);f+=1e-30*np.exp(-q);f/=np.trapezoid(q*q*f,q)
 return {'beta':float(beta),'p0_initial_units_mX':float(p0),'q0':float(q0),'velocity2_target':float(v2),'velocity4_target':float(v4),'fit_velocity_relative_error':(moments(rr.x)/[v2,v4]-1).tolist(),
 'max_node_midpoint_CDF_difference':float(max(abs(cdf-emp))),'kurtosis_p':float(ratio),'scope':'Three-parameter generalized Maxwell shape; number and final exact v²/v⁴ matched, full shape and perturbations are not fixed by these moments.'},np.c_[q,f]

def run(nodes=128,angles=16,source=96,emax=.05):
 tic=time.monotonic();tag=f'elastic_s{source}_n{nodes}_a{angles}';z=np.load(RF/'grut_m1_synthesis_20261029/free_stream.npz');u=z[f'p_at_2MeV_{source}']**2/2;w=z[f'weights_{source}'];N=float(w.sum());keep=u<=emax
 lost={'number_fraction':float(w[~keep].sum()/N),'energy_fraction':float(w[~keep]@u[~keep]/(w@u)),'energy2_fraction':float(w[~keep]@(u[~keep]**2)/(w@(u*u)))}
 grid=emax*(np.arange(nodes)/(nodes-1))**2;initial=project(grid,u[keep],w[keep]);initial/=initial.sum();ii,jj,D,K,omit,checks=operator(grid,angles);sym=np.where(ii==jj,.5,1.)
 def rhs(s,y):
  pairs=sym*y[ii]*y[jj];return np.r_[D@pairs,2*K@pairs,2*omit@pairs]
 x,temp,factor,exposure,bg=clock(N);t_eval=exposure[::5];xx=x[::5];tt=temp[::5];ff=factor[::5]
 sol=solve_ivp(rhs,(0,exposure[-1]),np.r_[initial,0,0],method='DOP853',rtol=1e-8,atol=1e-13,t_eval=t_eval)
 if not sol.success:raise RuntimeError(sol.message)
 pop=sol.y[:nodes];means=grid@pop;number=pop.sum(axis=0);events=np.array([2*K@(sym*v[ii]*v[jj]) for v in pop.T]);rates=events*ff;visc=2/3*rates
 dec=float(np.exp(np.interp(0,np.log(visc[::-1]),np.log(xx[::-1]))));tdec=float(np.interp(np.log(dec),np.log(xx),tt));final=pop[:,-1];fitted,psd=fit(grid,final,10.)
 thermal=maxwell(grid,2*means[0]/3);rr=sym*thermal[ii]*thermal[jj];therm_res=float(np.sum(abs(D@rr))/(4*(K@rr)))
 checks.update({'number_relative_drift':float(max(abs(number-number[0]))),'comoving_NR_energy_relative_drift':float(max(abs(means/means[0]-1))),
  'min_probability':float(pop.min()),'maxwell_stationarity_L1_per_four_events':therm_res,'initial_energy_projection_relative_error':float(means[0]/(w[keep]@u[keep]/w[keep].sum())-1)})
 out={'tag':tag,'nodes':nodes,'angles':angles,'source_bins':source,'energy_max':emax,'source_tail_removed':lost,'background':bg,'checks':checks,'number_per_reference_core':N,
  'initial_collision_rate_over_H':float(rates[0]),'final_collision_rate_over_H':float(rates[-1]),'initial_viscosity_proxy_over_H':float(visc[0]),'final_viscosity_proxy_over_H':float(visc[-1]),
  'viscosity_proxy_crossing_x':dec,'viscosity_proxy_crossing_temperature_MeV':tdec,'temperature_at_x10_MeV':float(tt[-1]),'collisions_per_particle_integrated':float(sol.y[-2,-1]),'omitted_collisions_per_particle_integrated':float(sol.y[-1,-1]),
  'initial_probability_zero_node':float(initial[0]),'final_probability_zero_node':float(final[0]),'mean_comoving_kinetic':float(means[-1]),'initial_kurtosis_p':float((grid*grid)@initial/means[0]**2),
  'final_kurtosis_p':float((grid*grid)@final/means[-1]**2),'fit':fitted,'nfev':sol.nfev,'seconds':time.monotonic()-tic,
  'scope':'Classical isotropic NR elastic gas with threshold s-wave amplitude and inherited homogeneous density. Viscosity rate is a proxy, not a perturbation eigenvalue. No cosmic formation, core reactions, metric perturbations or acoustic-to-free-streaming transfer is solved by this homogeneous calculation.'}
 np.savez_compressed(ROOT/(tag+'.npz'),energy=grid,x=xx,T_MeV=tt,probability=pop,collision_rate_over_H=rates,viscosity_proxy_over_H=visc,collisions=sol.y[-2],initial_maxwell=thermal)
 np.savetxt(ROOT/(tag+'_psd.dat'),psd,fmt='%.14e');(ROOT/(tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2),flush=True)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--nodes',type=int,default=128);p.add_argument('--angles',type=int,default=16);p.add_argument('--source',type=int,default=96);a=p.parse_args();run(a.nodes,a.angles,a.source)
