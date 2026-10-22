"""Exploratory action-level kinetic-pole construction; no empirical fit.

Restored from the recorded first execution after a transient workspace reset.
Same positive action, parameters, numerical gates and controls. All integrations
are repeated; old raw files are not claimed to have survived the reset.
"""
from pathlib import Path
import gzip, hashlib, io, json, math, platform, sys, tarfile, time
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

P=Path(__file__).resolve().parents[1]
A=math.sqrt(2/3)
INPUTS={
 'base_commit':'10e66b7c19a42aaf60fad9111d0d46644fa27774',
 'model_status':'proposed QV3-KP0/KP1 hybrid extensions; not adopted empirical targets',
 'credit':'Dimopoulos-Markkanen 1807.04359v3 eq 2.3 potential plus separate eq 3.6/3.8 kinetic illustrations; not their full eq 3.6',
 'full':dict(M=1.,m=10.,v=.2,lam=1.,c=1.,b=3.,rad_to_Higgs_V0=100.,h_over_v_initial=1e-6,N_end=3.),
 'numerics':dict(rtol=1e-9,atol=1e-12,max_step=.02,refined_rtol=1e-10,refined_atol=1e-13,refined_max_step=.01,samples=1501,method='Radau'),
 'physical':dict(M_GeV=2.435e18,v_GeV=246.,lam=.129,alpha=9.97e-6,mh_GeV=125.20,SM_Higgs_width_GeV=.0041,ATLAS_Binv95=.107,rho_Lambda_reference_GeV4=2.52e-47,H0_reference_kms_Mpc=67.4,Omega_m_reference=.315,Omega_r_reference=9.2e-5,hbar_GeV_s=6.582119569e-25,Mpc_km=3.085677581491367e19),
 'physical_provenance':'Rounded historical Higgs/inflation inputs; alpha calibrated within an inflation model, not a universal measured constant; cosmology already known',
 'late':dict(a_initial=1e-8,b_values=[1.,3.,10.],c_values=[1.,2.5],rtol=1e-10,atol=1e-12,max_step=.05,refined_max_step=.025,samples=1001),
 'c_status':'c=1 supplied benchmark; c=2.5 is already-dark-energy-calibrated literature choice, not independent prediction',
 'b_status':'supplied sensitivity values, no fit',
 'scope':'homogeneous classical action, compressed transition not physical SM thermal crossover; no absolute quantum vacuum selection, noise transport, held-out test or positive identification'
}

def dump(name,x):
 (P/name).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')

def phi_minimum(p):
 M,m,v,lam,c=[p[k] for k in ['M','m','v','lam','c']]
 V0=lam*v**4/4;eps=V0/(.75*m*m*M*M)
 def f(k):
  y=eps*k
  return math.exp(-y)*(-math.expm1(-y))/eps-c*math.exp(-2*c*y)
 k=brentq(f,.25*c,4*c,xtol=1e-14,rtol=1e-14)
 return M/A*(eps*k)

def potential(u,h,p,phi0):
 M,m,v,lam,c=[p[k] for k in ['M','m','v','lam','c']]
 phi=phi0*np.exp(u);x=A*phi/M;em=np.exp(-x);f=-np.expm1(-x);C=np.exp(-2*c*x)
 Vh=lam*(h*h-v*v)**2/4;Vp=.75*m*m*M*M*f*f
 U=Vp+C*Vh
 Uphi=1.5*m*m*M*A*em*f-(2*c*A/M)*C*Vh
 Uh=C*lam*h*(h*h-v*v)
 return U,phi*Uphi,Uh,phi,Vp,Vh

def metric(u,h,p,phi0,kind):
 phi=phi0*np.exp(u);B=p['b']*p['M']/p['v'];a=phi*phi+B*B*h*h
 if kind=='original':return a,0.,1.,a
 z=u+np.log(phi0/p['M']);c=B*B*h*z;d=1+B*B*z*z
 return a,c,d,phi*phi*d+B*B*h*h

def kinetic(u,h,up,hp,p,phi0,kind):
 phi=phi0*np.exp(u);B=p['b']*p['M']/p['v']
 extra=B*h*up if kind=='original' else B*(h*up+(u+np.log(phi0/p['M']))*hp)
 return (phi*up)**2+hp*hp+extra*extra

def full_rhs(kind,p,phi0,rr0):
 M=p['M'];B=p['b']*M/p['v'];L=np.log(phi0/M)
 def fun(N,y):
  u,h,up,hp,logH,W=y;H2=np.exp(2*logH)
  U,Uu,Uh,phi,_,_=potential(u,h,p,phi0)
  aa,cc,dd,det=metric(u,h,p,phi0,kind)
  kin=kinetic(u,h,up,hp,p,phi0,kind);rr=rr0*np.exp(-4*N)
  eps=-kin/(2*M*M)-2*rr/(3*M*M*H2)
  cu=phi*phi*up*up+2*B*B*h*up*hp
  ch=-B*B*h*up*up if kind=='original' else 2*B*B*(u+L)*up*hp
  fu=cu+Uu/H2;fh=ch+Uh/H2
  return [up,hp,-(3+eps)*up-(dd*fu-cc*fh)/det,-(3+eps)*hp-(-cc*fu+aa*fh)/det,eps,3*H2*kin]
 return fun

def full_run(kind,refined):
 p=INPUTS['full'];o=INPUTS['numerics'];phi0=phi_minimum(p)
 h0=p['v']*p['h_over_v_initial'];U0=potential(0,h0,p,phi0)[0]
 rr0=p['rad_to_Higgs_V0']*p['lam']*p['v']**4/4
 H0=np.sqrt((rr0+U0)/(3*p['M']**2));N=np.linspace(0,p['N_end'],o['samples'])
 start=time.perf_counter()
 sol=solve_ivp(full_rhs(kind,p,phi0,rr0),(0,p['N_end']),[0,h0,0,0,np.log(H0),0],t_eval=N,method=o['method'],rtol=o['refined_rtol'] if refined else o['rtol'],atol=o['refined_atol'] if refined else o['atol'],max_step=o['refined_max_step'] if refined else o['max_step'])
 name=kind+('_refined' if refined else '')
 if not sol.success:
  dump('FAILED_'+name+'.json',dict(message=sol.message,last_N=float(sol.t[-1])))
  raise RuntimeError((name,sol.message))
 u,h,up,hp,logH,W=sol.y;H2=np.exp(2*logH);rr=rr0*np.exp(-4*N)
 U,Uu,Uh,phi,Vp,Vh=potential(u,h,p,phi0)
 kin=kinetic(u,h,up,hp,p,phi0,kind);rs=.5*H2*kin+U;total=rs+rr
 fc=np.abs(3*p['M']**2*H2-total)/total;wc=np.abs(rs+W-U0)/U0
 speed=np.sqrt(kin);length=np.trapezoid(speed,N)
 z=p['b']*p['M']*h/p['v']*(u+np.log(phi0/p['M']))
 acc=1-kin/(2*p['M']**2)-2*rr/(3*p['M']**2*H2)
 data=np.column_stack([N,u,h/p['v'],phi/phi0,logH,rs,rr,U,Vp,Vh,W,fc,wc,z,speed,acc,up,hp])
 np.savetxt(P/(name+'.csv'),data,delimiter=',',header='N,log_phi_ratio,h_over_v,phi_ratio,log_H,rho_scalar,rho_r,U,V_phi,V_h,expansion_work,Friedmann_relative,work_relative,z,field_speed_per_N,ddot_a_over_aH2,u_per_N,h_per_N',comments='')
 r=dict(kind=kind,refined=refined,success=True,nfev=sol.nfev,seconds=time.perf_counter()-start,phi0=phi0,endpoint=dict(h_over_v=float(h[-1]/p['v']),log_phi_ratio=float(u[-1]),phi_ratio=float(phi[-1]/phi0),radiation_fraction=float(rr[-1]/total[-1]),scalar_fraction=float(rs[-1]/total[-1]),ddot_a_over_aH2=float(acc[-1])),max_Friedmann_relative=float(np.max(fc)),max_scalar_work_relative=float(np.max(wc)),length_over_M=float(length/p['M']),endpoint_z_displacement_over_M=float(abs(z[-1]-z[0])/p['M']),Hubble_length_upper_bound_over_M=float(math.sqrt(6)*p['N_end']))
 dump(name+'.json',r);print('DONE',name,flush=True);return r

def physical_scale(c=1.):
 p=INPUTS['physical'];M=p['M_GeV'];alpha=p['alpha'];v=p['v_GeV'];lam=p['lam']
 m=np.sqrt(4/3)*alpha*M;phi0=phi_minimum(dict(M=M,m=m,v=v,lam=lam,c=c))
 rho=.75*m*m*M*M*(-np.expm1(-A*phi0/M))**2
 approx=c*c*(lam*v**4/4)**2/(alpha*alpha*M**4)
 return dict(c=c,phi0_GeV=phi0,phi0_over_M=phi0/M,log_phi0_over_M=np.log(phi0/M),scalaron_m_GeV=m,rho_initial_GeV4=rho,rho_approx_GeV4=approx,energy_scale_meV=1e12*rho**.25,ratio_to_known_reference=rho/p['rho_Lambda_reference_GeV4'])

def late_run(b,c,refined):
 p=INPUTS['physical'];o=INPUTS['late'];M=p['M_GeV']
 Href=p['H0_reference_kms_Mpc']/p['Mpc_km']*p['hbar_GeV_s'];C=3*M*M*Href*Href
 amp=physical_scale(c)['rho_initial_GeV4']/C;om=p['Omega_m_reference'];orr=p['Omega_r_reference'];N0=np.log(o['a_initial'])
 def rhs(N,y):
  f,u,W=y;V=amp*np.exp(2*f/b);dm=om*np.exp(-3*N);rad=orr*np.exp(-4*N)
  E2=(dm+rad+V)/(1-u*u/6);eps=-u*u/2-(3*dm+4*rad)/(2*E2)
  return [u,-(3+eps)*u-(6/b)*V/E2,E2*u*u]
 N=np.linspace(N0,0,o['samples']);factor=.1 if refined else 1.
 sol=solve_ivp(rhs,(N0,0),[0,0,0],method='DOP853',t_eval=N,rtol=o['rtol']*factor,atol=o['atol']*factor,max_step=o['refined_max_step'] if refined else o['max_step'])
 if not sol.success:raise RuntimeError(sol.message)
 f,u,W=sol.y;V=amp*np.exp(2*f/b);dm=om*np.exp(-3*N);rad=orr*np.exp(-4*N)
 E2=(dm+rad+V)/(1-u*u/6);K=E2*u*u/6;w=(K-V)/(K+V);op=(K+V)/E2
 work=np.abs(K+V+W-amp)/amp
 r=dict(b=b,c=c,refined=refined,initial_displacement='prepared physical minimum, initial canonical velocity zero',rho_phi_today_GeV4=float(C*(K[-1]+V[-1])),w_today=float(w[-1]),Omega_phi_today=float(op[-1]),H_today_kms_Mpc=float(np.sqrt(E2[-1])*p['H0_reference_kms_Mpc']),log_phi_ratio_today=float(f[-1]/b),fraction_potential_retained_today=float(V[-1]/amp),max_frozen_change_before_a1e_3=float(np.max(abs(f[N<np.log(1e-3)]))),max_scalar_work_relative=float(np.max(work)),numerical_nfev=sol.nfev)
 name='late_b%g_c%g%s'%(b,c,'_refined' if refined else '')
 np.savetxt(P/(name+'.csv'),np.column_stack([N,np.exp(N),f,u,w,op,C*(K+V),np.sqrt(E2)*p['H0_reference_kms_Mpc'],W,work]),delimiter=',',header='N,a,canonical_displacement_over_M,velocity_per_N,w_phi,Omega_phi,rho_phi_GeV4,H_kms_Mpc,expansion_work_over_reference_critical,work_relative',comments='')
 dump(name+'.json',r);print('DONE',name,flush=True);return r

def analytical():
 p=INPUTS['physical'];M=p['M_GeV'];v=p['v_GeV'];mh=p['mh_GeV'];ss=[physical_scale(c) for c in [1.,2.5]]
 path=[]
 for b in [1.,3.,10.]:
  d=b*abs(ss[0]['log_phi0_over_M'])
  path.append(dict(b=b,minimum_length_over_M=d,minimum_efolds_all_positive_energy=d/math.sqrt(6),minimum_efolds_scalar_fraction_at_most_0p01=d/math.sqrt(.06)))
 width=mh**3/(32*math.pi*v*v);BR=width/(width+p['SM_Higgs_width_GeV'])
 ew=[dict(E_GeV=E,rho_GeV4=E**8/M**4,energy_scale_meV=1e12*E*E/M,ratio_to_reference=(E**8/M**4)/p['rho_Lambda_reference_GeV4']) for E in [160.,246.,1000.]]
 H=p['H0_reference_kms_Mpc']/p['Mpc_km']*p['hbar_GeV_s']
 window=[dict(c=s['c'],b_max_for_completed_path_within_one_efold=math.sqrt(6)/abs(s['log_phi0_over_M']),b_min_for_initial_effective_mass_at_most_reference_H=2*math.sqrt(s['rho_initial_GeV4'])/(M*H)) for s in ss]
 return dict(scales=ss,pure_E8_over_M4=ew,derivative_completion_path_bounds=path,completed_rapid_transition_vs_late_mass=window,original_portal=dict(width_GeV=width,reference_total_SM_width_GeV=p['SM_Higgs_width_GeV'],branching_fraction=BR,ATLAS_bound=p['ATLAS_Binv95'],assumptions='SM-like Higgs production, invisible light scalar and SM reference visible width'))

def controls():
 records=[]
 def ck(name,ok,value):
  records.append(dict(name=name,passed=bool(ok),evidence=value))
  if not ok:raise AssertionError((name,value))
 p=INPUTS['full'];phi0=phi_minimum(p);U,uu,uh,*_=potential(0.,0.,p,phi0)
 ck('prepared phi stationary before roll',abs(uu)/U<1e-12,float(abs(uu)/U))
 for kind in ['original','completed']:
  q=np.array([.03,.11]);velocity=np.array([.07,-.02]);step=1e-5
  def gm(q):
   a,c,d,det=metric(q[0],q[1],p,phi0,kind);return np.array([[a,c],[c,d]])
  der=np.array([(gm(q+np.eye(2)[j]*step)-gm(q-np.eye(2)[j]*step))/(2*step) for j in range(2)])
  cov=np.array([sum((der[j,i,k]+der[k,i,j]-der[i,j,k])/2*velocity[j]*velocity[k] for j in range(2) for k in range(2)) for i in range(2)])
  ph=phi0*np.exp(q[0]);B=p['b']/p['v'];z=q[0]+np.log(phi0)
  expected=np.array([ph*ph*velocity[0]**2+2*B*B*q[1]*velocity[0]*velocity[1],-B*B*q[1]*velocity[0]**2 if kind=='original' else 2*B*B*z*velocity[0]*velocity[1]])
  err=float(np.max(abs(cov-expected))/max(np.max(abs(expected)),1e-300))
  ck(kind+' metric variation connection',err<1e-8,err)
  det=metric(q[0],q[1],p,phi0,kind)[-1]
  ck(kind+' positive metric',det>0 and gm(q)[0,0]>0,float(det))
 mh=INPUTS['physical']['mh_GeV'];v=INPUTS['physical']['v_GeV']
 direct=quad(lambda cosine:(mh*mh/v)**2/(64*math.pi*mh),-1,1)[0];closed=mh**3/(32*math.pi*v*v)
 ck('portal phase space normalization',abs(direct/closed-1)<1e-14,float(direct))
 x=physical_scale();ck('physical remnant small-displacement limit',abs(x['rho_initial_GeV4']/x['rho_approx_GeV4']-1)<1e-13,x)
 return records

def archive_raw():
 """Retain every raw CSV, with deterministic container metadata."""
 manifest={};buf=io.BytesIO()
 with tarfile.open(fileobj=buf,mode='w') as archive:
  for path in sorted(P.glob('*.csv')):
   data=path.read_bytes();info=tarfile.TarInfo(path.name)
   info.size=len(data);info.mode=0o644;info.mtime=0
   archive.addfile(info,io.BytesIO(data))
   manifest[path.name]=dict(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
 (P/'RAW_DATA.tar.gz').write_bytes(gzip.compress(buf.getvalue(),mtime=0))
 dump('RAW_MANIFEST.json',manifest)
 for path in P.glob('*.csv'):path.unlink()

def main():
 dump('INPUTS.json',INPUTS);records=controls()
 full=[full_run(k,r) for k in ['original','completed'] for r in [False,True]]
 late=[];ref=[]
 for b in INPUTS['late']['b_values']:
  for c in INPUTS['late']['c_values']:
   late.append(late_run(b,c,False));ref.append(late_run(b,c,True))
 gates=[]
 for r in full:
  gates.append(dict(name=r['kind']+str(r['refined'])+' Friedmann',passed=r['max_Friedmann_relative']<1e-5,value=r['max_Friedmann_relative']))
  gates.append(dict(name=r['kind']+str(r['refined'])+' scalar work',passed=r['max_scalar_work_relative']<1e-6,value=r['max_scalar_work_relative']))
 for x,y in [(full[0],full[1]),(full[2],full[3])]:
  dh=abs(x['endpoint']['h_over_v']-y['endpoint']['h_over_v']);dp=abs(x['endpoint']['log_phi_ratio']-y['endpoint']['log_phi_ratio'])
  gates.append(dict(name=x['kind']+' refinement',passed=max(dh,dp)<1e-4,value=dict(h=dh,log_phi=dp)))
 for r in late+ref:gates.append(dict(name='late b%g c%g refined=%s scalar work'%(r['b'],r['c'],r['refined']),passed=r['max_scalar_work_relative']<1e-6,value=r['max_scalar_work_relative']))
 for x,y in zip(late,ref):
  d={k:abs(x[k]-y[k]) for k in ['w_today','Omega_phi_today','H_today_kms_Mpc']}
  gates.append(dict(name='late b%g c%g refinement'%(x['b'],x['c']),passed=all(d[k]<1e-6*abs(x[k])+1e-9 for k in d),value=d))
 r=dict(analytical=analytical(),full_action_runs=full,conditional_late_runs=late,conditional_late_refinements=ref,analytic_controls=records,numerical_gates=gates,status='PASS' if all(x['passed'] for x in records+gates) else 'FAIL')
 dump('RESULTS.json',r)
 dump('ENVIRONMENT.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=__import__('scipy').__version__,command='python3 code/build.py'))
 archive_raw()
 print('STATUS',r['status'],'CONTROLS',len(records),'GATES',len(gates))
 if r['status']!='PASS':sys.exit(1)

if __name__=='__main__':main()
