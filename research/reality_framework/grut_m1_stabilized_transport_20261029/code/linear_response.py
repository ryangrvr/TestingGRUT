"""Conditional parsec-scale response with a conserving BGK gas hierarchy.

This is an explicitly simplified test cosmology, NOT the full M1 power
spectrum: prescribed adiabatic formation matching at 2 MeV; all other matter
(including baryons) cold; radiation a single barotropic perfect fluid;
Maxwell or fitted generalized-Maxwell gas with the measured moments; relaxation-time collision
operator matched to the homogeneous viscosity proxy, not the exact XX
linearized collision integral. Radiation oscillations are averaged after a
controlled k/(aH) threshold. Phi=Psi; gas gravitational anisotropic stress
and O(v²) background-energy corrections are neglected.
"""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline,PchipInterpolator
from scipy.special import roots_genlaguerre,gamma,sici
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent;sys.path.insert(0,str(RF/'grut_m1_synthesis_20261029/code'))
import free_stream as fs

def collision_term(P,v,q,w,deriv,tbasis):
 """Unit-rate relaxation preserving the discrete number, energy and flow."""
 delta=w@P[:,0];ug=w@(v*P[:,1]);thetaT=(w@(q*q*P[:,0]))/(w@(q*q))-delta
 C=-P.copy();C[:,0]+=delta+tbasis*thetaT;C[:,1]-=deriv*ug/(3*v)
 return C

class Background:
 def __init__(self):
  a,E,ai,ae,meta=fs.background(4500,a_end=1.)
  self.ai=ai;self.loga=np.log(a);self.E=CubicSpline(self.loga,np.log(E));rho=E*E-fs.OM/a**3-fs.OL;self.lograd=CubicSpline(self.loga,np.log(rho));self.maxN=np.log(1/ai)
  gas=json.loads((ROOT/'elastic_s96_n512_a16.json').read_text());self.sigma2=2*gas['mean_comoving_kinetic']/3;self.rate0=gas['initial_viscosity_proxy_over_H'];self.rateN=gas['final_viscosity_proxy_over_H']
  z=np.load(ROOT/'elastic_s96_n512_a16.npz');self.rates=PchipInterpolator(np.log(z['x']),np.log(z['viscosity_proxy_over_H']));self.H10=float(self.E(np.log(ai*10)))
  self.fgas=json.loads((RF/'grut_m1_synthesis_20261029/free_stream.json').read_text())['populations']['96']['mass_fraction_of_cold_plus_held']*.12/(.12+.0224)
 def at(self,N):
  lna=np.log(self.ai)+N;a=np.exp(lna);E=np.exp(self.E(lna));Or=np.exp(self.lograd(lna))/E**2;Om=fs.OM/a**3/E**2
  if a<1e-7:
   w=-1-self.lograd(lna,1)/3;dw=-self.lograd(lna,2)/3;cs=w-dw/(3*(1+w))
  else:w=cs=1/3;dw=0.
  hubble=a*fs.H0/fs.C*E
  rate=np.exp(self.rates(N)) if N<=np.log(10) else self.rateN*np.exp(-4*(N-np.log(10))+self.H10-self.E(lna))
  return hubble,Or,Om,w,dw,cs,rate

def rd_initial(K):
 y=K/np.sqrt(3)
 if abs(y)<.02:
  phi=1-y*y/10+y**4/280-y**6/15120;dphi=-y*y/5+y**4/70-y**6/2520;integ=y*y/12-y**4/480+y**6/30240;uc=K/2*(1-y*y/20+y**4/840)
 else:
  phi=3*(np.sin(y)-y*np.cos(y))/y**3;dphi=3*np.sin(y)/y-3*phi;integ=np.log(y)+np.euler_gamma-1-sici(y)[1]+np.sin(y)/y;uc=3*np.sqrt(3)/y*(1-np.sin(y)/y)
 ur=K/2*(dphi+phi);dr=-2*K*K*phi/3-4*ur/K;dc=3*phi-4.5-9*integ
 return np.array([dr,ur,dc,uc])

def evolve(k,bg,nq=24,L=80,rtol=2e-7,collision=1.,switch=150.,thermal_seed=2/3,fgas=None,shape="maxwell"):
 tic=time.monotonic()
 if shape=='fit':
  fit=json.loads((ROOT/'elastic_s96_n512_a16.json').read_text())['fit'];beta=fit['beta'];q0=fit['p0_initial_units_mX'];nodes,w=roots_genlaguerre(nq,3/beta-1);w=w/gamma(3/beta);q=q0*nodes**(1/beta);deriv=-beta*nodes
 else:
  nodes,w=roots_genlaguerre(nq,.5);w=w/gamma(1.5);q=np.sqrt(2*bg.sigma2*nodes);deriv=-2*nodes
 energymean=w@(q*q);tbasis=-.5*(3+deriv);tbasis/=w@(q*q*tbasis)/energymean;fg=bg.fgas if fgas is None else fgas
 K0=k/bg.at(0)[0];cold=rd_initial(K0);psi=np.zeros((nq,L+1));psi[:,0]=cold[2]*(1+tbasis*thermal_seed);psi[:,1]=-deriv*cold[3]/(3*q)
 init=np.r_[cold,cold,psi.ravel()];ell=np.arange(2,L+1);lastI=8+nq*(L+1)
 def rhs(N,y,rad_on):
  hub,Or,Om,wr,dwr,csr,rate=bg.at(N);K=k/hub;v=q*np.exp(-N);P=y[8:].reshape(nq,L+1);delta=w@P[:,0];ug=w@(v*P[:,1]);thetaT=(w@(q*q*P[:,0]))/energymean-delta
  def fluid(d,gas):
   dr,ur,dc,uc=d;dm=(1-fg)*dc+fg*delta if gas else dc;um=(1-fg)*uc+fg*ug if gas else uc
   if not rad_on:dr=ur=0.
   rho=Or*dr+Om*dm;mom=(1+wr)*Or*ur+Om*um
   phi=-1.5/K**2*(rho+3*mom/K);phid=-phi+1.5*mom/K
   dd=np.array([-(1+wr)*(K*ur-3*phid)-3*(csr-wr)*dr,
      -(1-3*wr+dwr/(1+wr))*ur+K*(csr*dr/(1+wr)+phi),-K*uc+3*phid,-uc+K*phi])
   if not rad_on:dd[:2]=0.
   return dd,phi,phid
  dr,phi,phid=fluid(y[:4],True);dc,_,_=fluid(y[4:8],False);G=np.zeros_like(P);stream=K*v
  G[:,0]=-stream*P[:,1]-phid*deriv
  G[:,1]=stream/3*(P[:,0]-2*P[:,2])-K*phi/(3*v)*deriv
  G[:,2:]=stream[:,None]/(2*ell+1)[None,:]*(ell[None,:]*P[:,1:L]-(ell+1)[None,:]*np.c_[P[:,3:],np.zeros(nq)])
  if collision:G+=collision*rate*collision_term(P,v,q,w,deriv,tbasis)
  return np.r_[dr,dc,G.ravel()]
 # The radiation tail is an explicit approximation; switch refinement tests
 # its effect. All perturbations of the gas remain evolved after the switch.
 from scipy.optimize import brentq
 if K0>=switch:split=0.
 else:split=brentq(lambda N:k/bg.at(N)[0]-switch,0,min(bg.maxN,np.log(.001/bg.ai)))
 samples=np.linspace(0,bg.maxN,241);parts=[];last=init;total_eval=0
 if split>0:
  first=solve_ivp(lambda N,y:rhs(N,y,True),(0,split),last,method='DOP853',rtol=rtol,atol=rtol*1e-3,dense_output=True)
  if not first.success:raise RuntimeError(first.message)
  mask=samples<split;parts.append((samples[mask],first.sol(samples[mask])));last=first.y[:,-1];total_eval+=first.nfev
 second=solve_ivp(lambda N,y:rhs(N,y,False),(split,bg.maxN),last,method='DOP853',rtol=rtol,atol=rtol*1e-3,dense_output=True)
 if not second.success:raise RuntimeError(second.message)
 mask=samples>=split;parts.append((samples[mask],second.sol(samples[mask])));total_eval+=second.nfev;Y=np.concatenate([p[1] for p in parts],axis=1);Ps=Y[8:].reshape(nq,L+1,len(samples));dg=np.einsum('q,qt->t',w,Ps[:,0]);dm=(1-fg)*Y[2]+fg*dg;ratio=(dm/Y[6])**2
 return {'k_Mpc_inverse':k,'shape':shape,'nq':nq,'Lmax':L,'collision_multiplier':collision,'radiation_average_K':switch,'initial_thermal_seed':thermal_seed,'rtol':rtol,'fgas_of_matter':fg,
    'T2_z0':float(ratio[-1]),'T2_equality':float(np.interp(np.log(fs.OR/fs.OM/bg.ai),samples,ratio)),
    'initial_k_over_aH':K0,'nfev':total_eval,'seconds':time.monotonic()-tic},np.c_[samples,dm,Y[6],dg,ratio]

def main():
 p=argparse.ArgumentParser();p.add_argument('--tag',default='linear');p.add_argument('--ks',default='1000,3000,10000,30000,100000,300000,1000000');p.add_argument('--nq',type=int,default=24);p.add_argument('--L',type=int,default=100);p.add_argument('--rtol',type=float,default=2e-7);p.add_argument('--collision',type=float,default=1);p.add_argument('--switch',type=float,default=150);p.add_argument('--seed',type=float,default=2/3);p.add_argument('--shape',choices=['maxwell','fit'],default='maxwell');a=p.parse_args()
 bg=Background();rows=[];arrays={}
 for k in map(float,a.ks.split(',')):
  o,curve=evolve(k,bg,a.nq,a.L,a.rtol,a.collision,a.switch,a.seed,shape=a.shape);rows.append(o);arrays[f'k_{k:g}']=curve;print(o,flush=True)
 out={'rows':rows,'scope':__doc__,'background_parameters':{'h':fs.h,'Omega_m':fs.OM,'Omega_r':fs.OR,'a_initial':bg.ai,'sigma2_comoving_initial_mX':bg.sigma2,'gas_fraction_of_total_matter':bg.fgas},'columns':['ln(a/a_initial)','delta_mixture','delta_cold_comparator','delta_gas','T2']}
 np.savez_compressed(ROOT/(a.tag+'.npz'),**arrays);(ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
