#!/usr/bin/env python3
"""Continuum radial Jost/Fredholm determinant remainder, same fixed action.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python spectral.py
MSbar first/second Born graphs are outside this worker's scope.
"""
import json,sys,time,platform,hashlib,csv
from pathlib import Path
import numpy as np
import scipy
from scipy.special import ive,kve,gammaln
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss

OUT=Path(__file__).resolve().parent
REF=OUT.parent/'grut_candidate_1_20261018/code/results_profile.npz'
p=np.load(REF);SPL=scipy.interpolate.CubicSpline(p['rho'],p['y'],axis=1);OMEGA=float(p['Omega'])
U=np.array([[1,1,0],[1j,-1j,0],[0,0,np.sqrt(2)]],complex)/np.sqrt(2)
PROTOCOL={
 'reference':str(REF),'reference_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),'Omega':OMEGA,'a':1.,'b':.1,'D_mass2_vacuum':1.01,
 'operator':'P(zeta)=-d²+ell(ell+1)/r²+Klocal+zeta²+zeta*G, transformed radial measure dr',
 'free_kappa':'sqrt(1+(zeta±iOmega)²),sqrt(2+zeta²),sqrt(1.01+zeta²), positive real branches',
 'method':'continuum regular matrix ODE with QR stabilization; decaying-free Wronskian Jost determinant; exact continuum free Green-function Volterra first/second Born coefficients; scalar D Riccati/Born control',
 'truncated_potential_radius':20.,'origin_radius':.005,'QR_segment_length':1.,'ODE_method':'DOP853','rtol':2e-7,'atol':2e-10,
 'angular_cutoff':40,'zeta_intervals':[0.,.25,1.,3.,10.],'quadrature_orders':[6,8],
 'collective_factors':'For ell0 and1 subtract log[zeta²/(zeta²+1)] from numerical integrand and restore its exact integral0..Z=Z log[Z²/(Z²+1)]-2atan(Z). Full determinant and full vacuum UV degree retained.',
 'controls':{'vacuum':[[0,.3],[2,.7],[12,2.]],'weak_potential':{'ell':2,'zeta':.7,'scales':[.2,.1,.05]},'Dscalar':'independent scalar Riccati, first/second coefficients and weak cubic scaling recorded','radial_tolerance_probes':[[0,.005],[0,.1],[0,1.],[1,.005],[1,.1],[2,.7],[12,2.],[40,5.]],'tight_rtol':3e-8,'domain_probes':{'ell':[0,1,2,6],'zeta':[.02,.5,2.],'R':[15.,20.,25.]},'IR_zero_probes':{'ell':[0,1],'zeta':[.002,.005,.01,.02]}},
 'UV_and_angular_tails':'No fitted tails. Root separately supplies norm/HS bounds; record finite-cutoff remainder and explicit tail inputs. Potential endpoint controls are empirical unless separately bounded.',
 'budget_seconds':540,'budget_behavior':'stop starting physical points after budget, preserve completed points and report missing pieces',
 'units':'Energy remainder=(g*v)/(2pi) integral_0^infinity dzeta sum_ell(2ell+1) logdet_Born3plus. Numerical energy coefficient multiplies g*v.',
 'scope':'Tree-background coupled continuum determinant remainder only; renormalized low-order graphs, total fixed-Q correction, pole threshold comparison and any lifetime are not calculated here.',
 'frozen_before_physical_execution':True,
}
ROWS=[];START=None

def serial(o):
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 if isinstance(o,complex):return {'real':o.real,'imag':o.imag}
 raise TypeError(type(o).__name__)
def dump(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,default=serial)+'\n')
def log(obj):
 print(json.dumps(obj,default=serial),flush=True)
 with (OUT/'spectral_execution.jsonl').open('a') as f:f.write(json.dumps(obj,default=serial)+'\n')

def potential(r,scale=1.):
 F,_,S,_=SPL(r)
 real=np.array([[S*S+.3*F*F-1,0,2*F*S],[0,S*S+.1*F*F-1,0],[2*F*S,0,3*(S*S-1)+F*F]],complex)
 return scale*(U.conj().T@real@U),scale*.01*(S*S-1)


def free_data(ell,k,r):
 nu=ell+.5;z=k*r
 iv=ive(nu,z);kv=kve(nu,z)
 diag=r*iv*kv*np.exp(z.real-z)
 E=-ell/r-k*kve(nu-1,z)/kv
 return diag,E

def run_point(ell,zeta,R=20.,tol=2e-7,scale=1.,r0=.005,label='physical'):
 begin=time.perf_counter();k=np.sqrt(np.array([1+(zeta+1j*OMEGA)**2,1+(zeta-1j*OMEGA)**2,2+zeta*zeta],complex));kd=np.sqrt(1.01+zeta*zeta)
 C=np.diag(k*k);V0,vd0=potential(0.,scale);Veps,vdeps=potential(.001,scale);V2=(Veps-V0)/1e-6;vd2=(vdeps-vd0)/1e-6
 A=(C+V0)/(2*(2*ell+3));B=((C+V0)@A+V2)/(4*(2*ell+5))
 F=np.eye(3)+A*r0*r0+B*r0**4;Fp=2*A*r0+4*B*r0**3
 # Accurate near-origin first Born and Volterra state; free regular-function ratios avoid overflow.
 x,w=leggauss(24);rr=(x+1)*r0/2;ww=w*r0/2
 g0,e0=free_data(ell,k,r0);H=np.zeros((3,3),complex);born1=0j
 iv0=ive(ell+.5,k*r0)
 for t,weight in zip(rr,ww):
  V,_=potential(t,scale);gg,_=free_data(ell,k,t)
  born1+=weight*np.dot(np.diag(V),gg)
  ratio=np.sqrt(t/r0)*ive(ell+.5,k*t)/iv0*np.exp((k*t).real-(k*r0).real)
  H+=weight*V.T*(g0[:,None]*g0[None,:])*(ratio[:,None]*ratio[None,:])
 born2=-np.trace(V0@V0)*r0**4/(4*(2*ell+1)**2*(2*ell+3))
 den=(2*ell+3)**2*(2*ell+5);aD=vd0*r0/(2*ell+3)+(vd2/(2*ell+5)-2*kd*kd*vd0/den)*r0**3;bD=-vd0*vd0*r0**3/den;dD=aD+bD
 y=np.r_[F.ravel(),Fp.ravel(),H.ravel(),born1,born2,dD,aD,bD,0j]
 accum=0j;nfev=0;success=True;messages=[]
 def rhs(r,y):
  F=y[:9].reshape(3,3);Fp=y[9:18].reshape(3,3);HH=y[18:27].reshape(3,3);V,vd=potential(r,scale)
  gg,ee=free_data(ell,k,r)
  d,a,b=y[29:32]
  z=kd*r;R0=(ell+1)/r+kd*ive(ell+1.5,z)/ive(ell+.5,z)
  out=np.empty(33,complex);out[:9]=Fp.ravel();out[9:18]=((C+V)@F-2*(ell+1)/r*Fp).ravel()
  out[18:27]=((ee[:,None]+ee[None,:])*HH+(gg[:,None]*gg[None,:])*V.T).ravel()
  out[27]=np.dot(np.diag(V),gg);out[28]=-np.sum(V*HH)
  out[29]=vd-2*R0*d-d*d;out[30]=vd-2*R0*a;out[31]=-2*R0*b-a*a;out[32]=d-a-b
  return out
 edges=np.r_[r0,np.arange(1.,R,1.),R]
 for low,high in zip(edges[:-1],edges[1:]):
  sol=solve_ivp(rhs,(low,high),y,method='DOP853',rtol=tol,atol=tol*1e-3)
  nfev+=sol.nfev
  if not sol.success:success=False;messages.append(sol.message);break
  y=sol.y[:,-1];stack=np.vstack([y[:9].reshape(3,3),y[9:18].reshape(3,3)])
  Q,tri=np.linalg.qr(stack,mode='reduced');accum+=np.log(np.diag(tri).astype(complex)).sum();y[:9]=Q[:3].ravel();y[9:18]=Q[3:].ravel()
 if not success:
  rec={'ell':ell,'zeta':zeta,'R':R,'tol':tol,'scale':scale,'label':label,'success':False,'messages':messages,'seconds':time.perf_counter()-begin};ROWS.append(rec);return rec
 F=y[:9].reshape(3,3);Fp=y[9:18].reshape(3,3);gg,ee=free_data(ell,k,R)
 W=Fp+((ell+1)/R-ee)[:,None]*F
 sign,lnabs=np.linalg.slogdet(W)
 logb=.5*np.log(R)+np.log(kve(ell+.5,k*R))-k*R
 lognorm=gammaln(ell+1.5)+(ell+.5)*np.log(2/k)
 logJ=accum+lnabs+np.log(complex(sign))+np.sum(logb+(ell+1)*np.log(R)-lognorm)
 phase=np.angle(np.exp(1j*logJ.imag));B1=y[27];B2=y[28];rem=logJ.real-B1.real-B2.real
 d,a,b=y[29:32];z=kd*R;R0=(ell+1)/R+kd*ive(ell+1.5,z)/ive(ell+.5,z);E=-ell/R-kd*kve(ell-.5,z)/kve(ell+.5,z);T=R0-E
 drem=(y[32]+np.log1p(d/T)-a/T-b/T+.5*(a/T)**2).real
 rec={'ell':ell,'zeta':zeta,'R':R,'tol':tol,'scale':scale,'label':label,'success':True,'logJost_real':logJ.real,'logJost_phase_mod_2pi':phase,'Born1_real':B1.real,'Born1_imag':B1.imag,'Born2_real':B2.real,'Born2_imag':B2.imag,'coupled_remainder':rem,'D_remainder':drem,'total_remainder':rem+drem,'nfev':nfev,'seconds':time.perf_counter()-begin}
 ROWS.append(rec);return rec


def quadrature(n):
 x,w=leggauss(n);out=[]
 for low,high in zip(PROTOCOL['zeta_intervals'][:-1],PROTOCOL['zeta_intervals'][1:]):
  out.extend(zip(low+(x+1)*(high-low)/2,w*(high-low)/2))
 return out

def checkpoint(result):
 result['points']=ROWS;dump('spectral_results.json',result)
 fields=sorted(set().union(*(r.keys() for r in ROWS)))
 with (OUT/'spectral_points.csv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(ROWS)

def main():
 global START
 START=time.perf_counter();(OUT/'spectral_execution.jsonl').write_text('');dump('spectral_protocol.json',PROTOCOL)
 result={'protocol':PROTOCOL,'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'controls':[],'partial_integrals':[],'skipped':[]}
 for ell,zeta in PROTOCOL['controls']['vacuum']:
  rec=run_point(ell,zeta,scale=0.,label='vacuum');result['controls'].append(rec);log(rec)
 if max(abs(r['total_remainder']) for r in result['controls'])>2e-5:
  result['error']='vacuum control failed';checkpoint(result);return
 for scale in [.2,.1,.05]:
  rec=run_point(2,.7,scale=scale,label='weak_potential');result['controls'].append(rec);log(rec)
 for ell,zeta in PROTOCOL['controls']['IR_zero_probes']['zeta'] if False else []:pass
 for n in PROTOCOL['quadrature_orders']:
  integrals=[];complete=True
  for ell in range(41):
   subtotal=0.;dtotal=0.;maxphase=0.
   for zeta,weight in quadrature(n):
    if time.perf_counter()-START>540:complete=False;result['skipped'].append({'quadrature_order':n,'ell':ell,'zeta':zeta,'reason':'frozen time budget'});break
    rec=run_point(ell,float(zeta),label=f'quadrature_n{n}')
    if not rec['success']:complete=False;continue
    zero=np.log(zeta*zeta/(zeta*zeta+1)) if ell in [0,1] else 0.
    subtotal+=weight*(rec['coupled_remainder']-zero);dtotal+=weight*rec['D_remainder'];maxphase=max(maxphase,abs(rec['logJost_phase_mod_2pi']))
   if ell in [0,1]:
    Z=10.;subtotal+=Z*np.log(Z*Z/(Z*Z+1))-2*np.arctan(Z)
   entry={'ell':ell,'quadrature_order':n,'coupled_integral':subtotal,'D_integral':dtotal,'weighted_energy_coefficient':(2*ell+1)*(subtotal+dtotal)/(2*np.pi),'max_Jost_phase':maxphase,'complete':complete}
   integrals.append(entry);log(entry);result['partial_integrals'].append(entry);checkpoint(result)
   if not complete:break
  result[f'energy_n{n}']=sum(r['weighted_energy_coefficient'] for r in integrals)
  result[f'complete_n{n}']=complete and len(integrals)==41
  if not complete:break
 # Declared radial/domain/IR controls, bounded by the same global budget.
 for ell,zeta in PROTOCOL['controls']['radial_tolerance_probes']:
  for tol in [2e-7,3e-8]:
   if time.perf_counter()-START>540:result['skipped'].append({'control':'radial','ell':ell,'zeta':zeta});continue
   rec=run_point(ell,zeta,tol=tol,label='radial_tolerance');result['controls'].append(rec);log(rec)
 for ell in [0,1,2,6]:
  for zeta in [.02,.5,2.]:
   for R in [15.,20.,25.]:
    if time.perf_counter()-START>540:result['skipped'].append({'control':'domain','ell':ell,'zeta':zeta,'R':R});continue
    rec=run_point(ell,zeta,R=R,label='domain');result['controls'].append(rec);log(rec)
 for ell in [0,1]:
  for zeta in [.002,.005,.01,.02]:
   if time.perf_counter()-START>540:continue
   rec=run_point(ell,zeta,tol=3e-8,label='IR_zero');rec['logJost_minus_log_zeta_squared']=rec.get('logJost_real',np.nan)-2*np.log(zeta);result['controls'].append(rec);log(rec)
 result['elapsed_seconds']=time.perf_counter()-START;result['completed']=True;checkpoint(result);log({'completed':True,'elapsed_seconds':result['elapsed_seconds'],'energy_n6':result.get('energy_n6'),'energy_n8':result.get('energy_n8'),'skipped_count':len(result['skipped'])})

if __name__=='__main__':main()
