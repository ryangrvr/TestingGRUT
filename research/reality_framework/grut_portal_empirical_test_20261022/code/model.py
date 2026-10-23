"""Tree portal recoil model and fixed public-source conventions."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.special import erf,spherical_jn

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_unified_working_theory_20261022'
G,LAM,HD=.1,.01,.01
MH,VH,MN,FN=125.,246.,.939,.308
CKMS,CCMS,HBARC=299792.458,29979245800.,.1973269804
GEV2CM2,AMUG,YEAR=.3893793721e-27,1.66053906660e-24,365.25*86400
V0,VE,VESC,RHO=220.,232.,544.,.3
MPLANCK=2.435e18
ISOTOPES=np.array([
 [124,123.9058920,.000952],[126,125.9042983,.000890],
 [128,127.9035310,.019102],[129,128.9047808611,.264006],
 [130,129.903509349,.040710],[131,130.90508406,.212324],
 [132,131.9041550856,.269086],[134,133.90539466,.104357],
 [136,135.907214484,.088573]])
NZ=erf(VESC/V0)-2*(VESC/V0)*np.exp(-(VESC/V0)**2)/np.sqrt(np.pi)


def eta(vmin):
 x=np.asarray(vmin)/V0;y=VE/V0;z=VESC/V0
 a=erf(x+y)-erf(x-y)-4*y*np.exp(-z*z)/np.sqrt(np.pi)
 b=erf(z)-erf(x-y)-2*(z+y-x)*np.exp(-z*z)/np.sqrt(np.pi)
 return np.where(x<z-y,a,np.where(x<z+y,b,0))/(2*NZ*VE)


def helm(q,a,legacy=False):
 skin=.9
 r=np.sqrt((1.2*a**(1/3))**2-5*skin**2) if legacy else np.sqrt((1.23*a**(1/3)-.6)**2+7*np.pi**2*.52**2/3-5*skin**2)
 x=q*r/HBARC
 return 3*spherical_jn(1,x)/np.where(x==0,1,x)*np.exp(-.5*(q*skin/HBARC)**2)+np.where(x==0,1,0)


class Object:
 def __init__(self,charge):
  self.charge=charge
  old=json.loads((PARENT/'RESULTS.json').read_text())
  self.row=next(r for r in old['profiles'] if r['Q']==charge)
  data=np.loadtxt(PARENT/'data'/f'profile_Q{charge}.csv',delimiter=',',skiprows=1)
  self.rho=data[:,0];self.s=data[:,3]
  self.weight=self.rho**2*(1-self.s**2)
  self.norm=simpson(self.weight,x=self.rho)
  qh=np.linspace(0,5,2501)
  vals=[]
  for qs in np.array_split(qh,25):
   vals.extend(simpson(self.weight[None,:]*np.sinc(qs[:,None]*self.rho[None,:]/np.pi),x=self.rho,axis=1)/self.norm)
  self.form=CubicSpline(qh,vals,extrapolate=False)

 def mass(self,v):return G*v*self.row['mass_over_gv']

 def sigma_n(self,v,kappa):
  mass=self.mass(v);mu=mass*MN/(mass+MN)
  c=kappa*FN*MN*mass/(2*LAM*v*v*MH*MH)
  return mu*mu*c*c/np.pi*GEV2CM2

 def rate(self,e,v,kappa=1e-9,point=False,natural=True,legacy=False,standard_mass=None,standard_sigma=None):
  e=np.asarray(e);mass=self.mass(v) if standard_mass is None else standard_mass
  sn=self.sigma_n(v,kappa) if standard_sigma is None else standard_sigma
  mun=mass*MN/(mass+MN)
  isotopes=ISOTOPES if natural else np.array([[131,131.,1.]])
  mean_u=np.sum(isotopes[:,1]*isotopes[:,2]);ntarget=1e6/(mean_u*AMUG)
  ans=np.zeros_like(e,dtype=float)
  for a,u,f in isotopes:
   ma=u*.93149410242-(54*.00051099895 if natural else 0)
   mua=mass*ma/(mass+ma);q=np.sqrt(2*ma*e*1e-6)
   sa=sn*a*a*(mua/mun)**2
   fobj=np.ones_like(e) if point or standard_mass is not None else self.form(q/(G*v))
   vmin=CKMS*np.sqrt(ma*e*1e-6/(2*mua*mua))
   ans+=f*ntarget*RHO/mass*CCMS*ma*sa/(2*mua*mua)*1e-6*CKMS*eta(vmin)*YEAR*helm(q,a,legacy)**2*fobj**2*(MH*MH/(MH*MH+q*q))**2
  return ans


def bin_integrals(fun,edges,order=8):
 x,w=np.polynomial.legendre.leggauss(order)
 center=(edges[:-1]+edges[1:])/2;width=np.diff(edges)/2
 ee=center[:,None]+width[:,None]*x
 return width*np.sum(fun(ee)*w[None,:],axis=1)
