#!/usr/bin/env python3
"""Fixed-coupling finite-charge branch and explicit fragmentation candidates.

Run: python branch_scan.py (from any working directory).
Physical parameters are never fit or changed. All attempted numerical starts are
logged, including convergence failures and vacuum convergence. No global minimum
or unsampled-partition assertion is inferred from a finite scan.
"""
from __future__ import annotations
import csv,json,sys,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_bvp,simpson,cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar,root_scalar

OUT=Path(__file__).resolve().parent
REFERENCE=OUT.parent/'grut_candidate_1_20261018/code/results_profile.npz'
PROTOCOL={
 'physical':{'a':1.,'b':.1,'g':.1,'lambda':.01,'lambda_X':.001},
 'reference_profile':str(REFERENCE),
 'main_omega':[round(x,2) for x in np.arange(.80,1.,.01)]+[.995,.9975,.999],
 'main_grid':{'L':'max(40,14/sqrt(1-Omega^2))','initial_spacing_max':.1,'tol':1e-6,'max_nodes':30000},
 'refinement_grid':{'L':'max(50,18/sqrt(1-Omega^2))','initial_spacing_max':.06,'tol':1e-8,'max_nodes':40000},
 'alternate_starts':[{'Omega':o,'A':a,'width':w,'kind':'gaussian'} for o in [.90,.97] for a,w in [(.3,8.),(2.3,3.5)]],
 'turning_point_refinement':{'method':'bounded minimize_scalar of Q(Omega) in adjacent-sample bracket','maxiter':16,'xatol':2e-7},
 'binding_point_refinement':{'method':'brentq on E/Q-1 in sampled sign-change bracket','maxiter':16,'xtol':2e-8},
 'fragments':{'equal_counts':[2,3,4,5],'unequal_two_fractions':[.20,.25,1/3,.40,.50],'lump_plus_free_fractions':[.25,.5,.75]},
 'target_resolution':'Fixed-Q BVP with Omega as unknown and integral charge auxiliary equation; try each sampled monotone branch containing target Q.',
 'acceptance':'BVP status0, positive Omega<1, nodeless F up to1e-7, 0<S<=1 up to1e-6, F0>1e-4; vacuum and failed starts retained.',
 'scope':'Found nodeless spherical tree-level flat-space states and positive-charge fragmentation candidates plus free X carriers; no global minimization, no proof of all partitions.',
 'frozen_before_first_scan':True,
}

def serial(o):
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 raise TypeError(type(o).__name__)

def dump(path,data):path.write_text(json.dumps(data,indent=2,default=serial)+'\n')
def log(rec):
 print(json.dumps(rec,default=serial),flush=True)
 with (OUT/'branch_execution.jsonl').open('a') as f:f.write(json.dumps(rec,default=serial)+'\n')

def bc_fixed(ya,yb):return np.array([ya[1],ya[3],yb[0],yb[2]-1])
def rhs(r,y,O):
 F,Fp,S,Sp=y[:4]
 return np.array([Fp,(S*S+.1*F*F-O*O)*F,Sp,((S*S-1)+F*F)*S])

def grid(O,refine=False):
 k=np.sqrt(max(1e-9,1-O*O));L=max(50 if refine else 40,(18 if refine else 14)/k)
 return np.linspace(0,L,int(np.ceil(L/(.06 if refine else .1)))+1)

class Saved:
 def __init__(self,r,y,O):self.x=r;self.sol=CubicSpline(r,y,axis=1);self.O=O

def observables(sol,O):
 r=np.linspace(0,sol.x[-1],max(12001,int(sol.x[-1]/.004)+1));F,Fp,S,Sp=sol.sol(r)[:4]
 U=.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4
 grad=.5*(Fp*Fp+Sp*Sp);tim=.5*O*O*F*F
 integ=lambda f:4*np.pi*simpson(r*r*f,x=r)
 Q=O*integ(F*F);E=integ(grad+tim+U);K=integ(grad);V=integ(U);T=integ(tim)
 cum=4*np.pi*cumulative_trapezoid(r*r*(grad+tim+U),r,initial=0)
 return {'Omega':float(O),'Qbar':Q,'Ebar':E,'E_over_Q':E/Q if Q else np.inf,'F0':F[0],'S0':S[0],'R90bar':np.interp(.9*cum[-1],cum,r),'gradient_energy':K,'virial_relative':(K+3*(V-T))/E if E else np.inf,'E_minus_OmegaQ':E-O*Q,'virial_identity_E_minus_OmegaQ_minus_2K_over3':E-O*Q-2*K/3,'L':r[-1]}

ATTEMPTS=[];STATES=[];SOLUTIONS=[]
def solve_fixed(O,seed,label,refine=False,guess=None,domain_override=None,tol_override=None):
 r=grid(O,refine) if domain_override is None else np.linspace(0,domain_override,int(np.ceil(domain_override/.05))+1)
 if guess is None:
  y=seed.sol(np.minimum(r,seed.x[-1]))[:4]
  outside=r>seed.x[-1];y[:,outside]=np.array([0.,0.,1.,0.])[:,None]
 else:
  F=guess['A']*np.exp(-(r/guess['width'])**2)
  S=1-.99*np.exp(-(r/guess['width'])**2)
  y=np.array([F,-2*r/guess['width']**2*F,S,2*r/guess['width']**2*(1-S)])
 start=time.perf_counter()
 sol=solve_bvp(lambda r,y:rhs(r,y,O),bc_fixed,r,y,S=np.diag([0.,-2.,0.,-2.]),tol=(tol_override if tol_override is not None else (1e-8 if refine else 1e-6)),max_nodes=40000 if refine else 30000)
 rr=np.linspace(0,r[-1],12001);yy=sol.sol(rr)
 accepted=sol.success and yy[0,0]>1e-4 and yy[0].min()>-1e-7 and yy[2].min()>0 and yy[2].max()<1.000001
 rec={'label':label,'Omega':O,'refined':refine,'guess':guess,'status':sol.status,'message':sol.message,'mesh_nodes':len(sol.x),'max_residual':np.max(sol.rms_residuals),'minF':yy[0].min(),'minS':yy[2].min(),'F0':yy[0,0],'accepted':bool(accepted),'seconds':time.perf_counter()-start}
 ATTEMPTS.append(rec);log(rec)
 if accepted:
  obs=observables(sol,O);obs['label']=label;obs['refined']=refine
  STATES.append(obs);SOLUTIONS.append((O,sol,obs));return sol,obs
 return None,None

def nearest(O):return min(SOLUTIONS,key=lambda x:abs(x[0]-O))[1]
def save_progress(extra=None):
 data={'protocol':PROTOCOL,'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'attempts':ATTEMPTS,'states':STATES}
 if extra:data.update(extra)
 dump(OUT/'branch_results.json',data)


def target_solve(target,seed,O,label,refine=False):
 # Auxiliary J=r-integral of charge. Six BC for five ODE plus unknown Omega.
 r=grid(O,refine);y4=seed.sol(np.minimum(r,seed.x[-1]))[:4];outside=r>seed.x[-1];y4[:,outside]=np.array([0.,0.,1.,0.])[:,None]
 j=4*np.pi*O*cumulative_trapezoid(r*r*y4[0]**2,r,initial=0);j*=target/max(j[-1],1e-12)
 y=np.vstack([y4,j])
 def fun(r,y,p):return np.vstack([rhs(r,y,p[0]),4*np.pi*p[0]*r*r*y[0]**2])
 def bc(ya,yb,p):return np.array([ya[1],ya[3],yb[0],yb[2]-1,ya[4],yb[4]-target])
 start=time.perf_counter()
 sol=solve_bvp(fun,bc,r,y,p=np.array([O]),S=np.diag([0.,-2.,0.,-2.,0.]),tol=1e-8 if refine else 1e-6,max_nodes=40000)
 oo=float(sol.p[0]);rr=np.linspace(0,r[-1],12001);yy=sol.sol(rr)
 accepted=sol.success and 0<oo<1 and yy[0,0]>1e-4 and yy[0].min()>-1e-7 and yy[2].min()>0 and yy[2].max()<1.000001
 rec={'label':label,'target_Qbar':target,'Omega':oo,'refined':refine,'status':sol.status,'message':sol.message,'mesh_nodes':len(sol.x),'max_residual':np.max(sol.rms_residuals),'minF':yy[0].min(),'minS':yy[2].min(),'accepted':bool(accepted),'seconds':time.perf_counter()-start}
 ATTEMPTS.append(rec);log(rec)
 if accepted:
  obs=observables(sol,oo);obs.update({'label':label,'refined':refine,'target_Qbar':target,'charge_residual':obs['Qbar']-target})
  STATES.append(obs);return sol,obs
 return None,None


def main():
 (OUT/'branch_execution.jsonl').write_text('')
 dump(OUT/'branch_protocol.json',PROTOCOL)
 p=np.load(REFERENCE);seed=Saved(p['rho'],p['y'],float(p['Omega']))
 ref=observables(seed,float(p['Omega']));Qstar=ref['Qbar'];Estar=ref['Ebar']
 maincurve=[]
 for O in PROTOCOL['main_omega']:
  sol,obs=solve_fixed(O,seed,'continuation')
  if sol is not None:seed=sol;maincurve.append((O,sol,obs))
  save_progress({'reference':ref})
 for guess in PROTOCOL['alternate_starts']:
  solve_fixed(guess['Omega'],nearest(guess['Omega']),'alternate_seed',guess=guess)
  save_progress({'reference':ref})
 # Locate turning points by minimum Q in sampled curve, only if bracketed.
 curve=sorted(maincurve);q=np.array([s[2]['Qbar'] for s in curve]);imin=int(np.argmin(q));thresholds={}
 if 0<imin<len(curve)-1:
  bracket=(curve[imin-1][0],curve[imin+1][0])
  def objective(O):
   sol,obs=solve_fixed(float(O),nearest(O),'turning_refinement')
   return obs['Qbar'] if obs else 1e50
  opt=minimize_scalar(objective,bounds=bracket,method='bounded',options={'maxiter':16,'xatol':2e-7})
  ss,oo=solve_fixed(float(opt.x),nearest(opt.x),'turning_confirmation',refine=True)
  thresholds['turning']={'bracket':bracket,'optimizer_success':bool(opt.success),'optimizer_message':str(opt.message),'observation':oo}
 else:thresholds['turning']={'status':'minimum not interior to sampled branch'}
 # Every sampled sign change gets its own root and refinement; no hidden omission.
 bindings=[]
 for low,high in zip(curve[:-1],curve[1:]):
  if (low[2]['E_over_Q']-1)*(high[2]['E_over_Q']-1)<0:
   def binding(O):
    ss,oo=solve_fixed(float(O),nearest(O),'binding_refinement')
    if oo is None:raise RuntimeError('Binding refinement BVP failed')
    return oo['E_over_Q']-1
   root=root_scalar(binding,bracket=(low[0],high[0]),method='brentq',xtol=2e-8,maxiter=16)
   ss,oo=solve_fixed(root.root,nearest(root.root),'binding_confirmation',refine=True)
   bindings.append({'bracket':[low[0],high[0]],'root_converged':root.converged,'observation':oo})
 thresholds['binding']=bindings
 save_progress({'reference':ref,'thresholds':thresholds})
 # Split sampled connected curve at its minimum Q. Exact-charge solves try both.
 segments=[('lower_Omega',curve[:imin+1]),('higher_Omega',curve[imin:])]
 fractions=sorted(set([1/n for n in PROTOCOL['fragments']['equal_counts']]+PROTOCOL['fragments']['unequal_two_fractions']+[1-x for x in PROTOCOL['fragments']['unequal_two_fractions']]+PROTOCOL['fragments']['lump_plus_free_fractions']))
 targets={};turnq=thresholds.get('turning',{}).get('observation',{}).get('Qbar',float(np.min(q)))
 for frac in fractions:
  target=frac*Qstar;found=[]
  for name,seg in segments:
   if len(seg)<2:continue
   qlo=min(a[2]['Qbar'] for a in seg);qhi=max(a[2]['Qbar'] for a in seg)
   if not qlo<=target<=qhi:continue
   O,ss,oo=min(seg,key=lambda a:abs(a[2]['Qbar']-target))
   candidate,obs=target_solve(target,ss,O,'target_'+name)
   if obs is not None:
    fine,ff=target_solve(target,candidate,obs['Omega'],'target_'+name+'_refined',refine=True)
    found.append({'segment':name,'coarse':obs,'fine':ff,'energy_refinement_abs':abs(ff['Ebar']-obs['Ebar']) if ff else None,'charge_refinement_abs':abs(ff['Qbar']-obs['Qbar']) if ff else None})
  targets[str(frac)]={'fraction':frac,'target_Qbar':target,'found':found,'below_found_turning_charge':target<turnq}
  save_progress({'reference':ref,'thresholds':thresholds,'targets':targets})
 def choose(frac):
  rec=min(targets.values(),key=lambda a:abs(a['fraction']-frac))
  if abs(rec['fraction']-frac)>1e-12:return None
  vals=[a['fine'] or a['coarse'] for a in rec['found']]
  return min(vals,key=lambda a:a['Ebar']) if vals else None
 comparisons=[]
 comparisons.append({'kind':'all_free_carriers','Ebar':Qstar,'delta_Ebar':Qstar-Estar,'resolved':True})
 for n in PROTOCOL['fragments']['equal_counts']:
  s=choose(1/n);comparisons.append({'kind':'equal_fragments','count':n,'resolved':s is not None,'fragment':s,'Ebar':n*s['Ebar'] if s else None,'delta_Ebar':n*s['Ebar']-Estar if s else None})
 for frac in PROTOCOL['fragments']['unequal_two_fractions']:
  s,t=choose(frac),choose(1-frac);ok=s is not None and t is not None
  comparisons.append({'kind':'unequal_two_fragments','fraction':frac,'resolved':ok,'fragments':[s,t],'Ebar':s['Ebar']+t['Ebar'] if ok else None,'delta_Ebar':s['Ebar']+t['Ebar']-Estar if ok else None})
 for frac in PROTOCOL['fragments']['lump_plus_free_fractions']:
  s=choose(frac);e=s['Ebar']+(1-frac)*Qstar if s else None
  comparisons.append({'kind':'one_lump_plus_free','fraction_in_lump':frac,'resolved':s is not None,'fragment':s,'Ebar':e,'delta_Ebar':e-Estar if e else None})
 # Numerical first-law check, and the conditional continuum comparison identity.
 sortedq=sorted([a[2] for a in curve[:imin+1]],key=lambda s:s['Qbar'])
 low_envelope_checks={'sampled_lower_Omega_E_over_Q_decreases_with_Q':all(b['E_over_Q']<a['E_over_Q'] for a,b in zip(sortedq,sortedq[1:])),'minimum_sampled_E_over_Q_for_Q_le_Qstar':min(s['E_over_Q'] for _,_,s in curve if s['Qbar']<=Qstar*(1+1e-6)),'reference_E_over_Q':Estar/Qstar,'identity':'On each smooth stationary branch dE/dQ=Omega. Virial gives E-Omega Q=(2/3)K>0, so d(E/Q)/dQ=-(2K)/(3Q^2)<0. If every allowed lower-charge component in the found continuous envelope has E(q)/q >= E(Qstar)/Qstar and free carriers cost1, every positive-charge partition drawn from that envelope costs at least E(Qstar). Missing states/branches are not excluded.'}
 save_progress({'reference':ref,'thresholds':thresholds,'targets':targets,'comparisons':comparisons,'conditional_envelope':low_envelope_checks})
 # Export all stationary states and sampled profiles for independent checking.
 with (OUT/'branch_states.csv').open('w',newline='') as f:
  fields=list(STATES[0]);writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(s.keys() for s in STATES))))
  writer.writeheader();writer.writerows(STATES)
 arrays={}
 for i,(O,s,obs) in enumerate(curve):
  r=np.linspace(0,s.x[-1],4001);arrays[f'rho_{i:03}']=r;arrays[f'y_{i:03}']=s.sol(r);arrays[f'Omega_{i:03}']=O
 np.savez(OUT/'branch_profiles.npz',**arrays)
 log({'complete':True,'accepted_states':len(STATES),'attempts':len(ATTEMPTS),'reference_Qbar':Qstar,'reference_Ebar':Estar,'thresholds':thresholds,'comparisons':[{'kind':c['kind'],'count':c.get('count'),'fraction':c.get('fraction',c.get('fraction_in_lump')),'resolved':c['resolved'],'delta_Ebar':c['delta_Ebar']} for c in comparisons]})



def refine_followup():
 """Predeclared, bounded second pass; first pass remains archived unchanged."""
 from scipy.sparse import diags,bmat
 from scipy.sparse.linalg import eigsh,LinearOperator
 global ATTEMPTS,STATES,SOLUTIONS
 data=json.loads((OUT/'branch_results.json').read_text())
 if not (OUT/'branch_first_pass_results.json').exists():dump(OUT/'branch_first_pass_results.json',data)
 refinement_protocol={
  'reason':'Initial 16-call turning optimizer hit its declared cap; confirm turning by linear BVP sensitivity and extend explicit multi-component comparisons.',
  'turning':{'method':'root of dQ/dOmega from differentiated fixed-frequency ODE','bracket':[.958,.962],'domain':80,'tol':1e-9,'xtol':2e-9,'maxiter':12,'confirmation_domain':100,'confirmation_tol':3e-10},
  'binding':{'method':'brentq E/Q-1 with refined BVP','bracket':[.90649,.90659],'domain':60,'tol':1e-9,'xtol':2e-10,'maxiter':12,'confirmation_domain':80,'confirmation_tol':3e-10},
  'additional_exact_charge_fractions':[.375,.45,.55,.625,1.],
  'additional_multicomponent_candidates':[[.375,.375,'free .25'],[.4,.4,'free .2'],[.45,.45,'free .1']],
  'amplitude_l0_checks':{'Omega':[.95,.97],'N':[800,1200],'L':70,'lowest_count':2,'fixed_Q_rank_one':True},
  'first_pass_preserved':True,'physical_parameters_unchanged':True,'frozen_before_second_pass':True,
 }
 dump(OUT/'branch_refinement_protocol.json',refinement_protocol)
 ATTEMPTS=data['attempts'];STATES=data['states']
 pp=np.load(OUT/'branch_profiles.npz');i=0
 while f'Omega_{i:03}' in pp:
  O=float(pp[f'Omega_{i:03}']);ss=Saved(pp[f'rho_{i:03}'],pp[f'y_{i:03}'],O)
  SOLUTIONS.append((O,ss,observables(ss,O)));i+=1
 sensitivity_records=[]
 def slope(O):
  ss,obs=solve_fixed(float(O),nearest(O),'turning_sensitivity_profile',refine=True,domain_override=80,tol_override=1e-9)
  if ss is None:raise RuntimeError('Turning profile failed')
  r=np.linspace(0,80,1601)
  def fun(r,z):
   F,_,S,_=ss.sol(r);u,up,w,wp=z
   return np.array([up,(S*S+.3*F*F-O*O)*u+2*F*S*w-2*O*F,wp,2*F*S*u+(3*S*S-1+F*F)*w])
  def bc(ya,yb):return np.array([ya[1],ya[3],yb[0],yb[2]])
  zz=solve_bvp(fun,bc,r,np.zeros((4,len(r))),S=np.diag([0.,-2.,0.,-2.]),tol=1e-9,max_nodes=20000)
  rr=np.linspace(0,80,30001);F=ss.sol(rr)[0];u=zz.sol(rr)[0]
  derivative=4*np.pi*simpson(rr*rr*(F*F+2*O*F*u),x=rr)
  rec={'Omega':O,'dQbar_dOmega':derivative,'sensitivity_status':zz.status,'sensitivity_max_residual':np.max(zz.rms_residuals),'Qbar':obs['Qbar']}
  sensitivity_records.append(rec);log({'sensitivity':rec})
  return derivative
 rt=root_scalar(slope,bracket=(.958,.962),method='brentq',xtol=2e-9,maxiter=12)
 st,ot=solve_fixed(rt.root,nearest(rt.root),'turning_final_domain80',refine=True,domain_override=80,tol_override=1e-9)
 st2,ot2=solve_fixed(rt.root,st,'turning_final_domain100',refine=True,domain_override=100,tol_override=3e-10)
 def binding(O):
  ss,oo=solve_fixed(float(O),nearest(O),'binding_high_accuracy',refine=True,domain_override=60,tol_override=1e-9)
  if oo is None:raise RuntimeError('Binding profile failed')
  return oo['E_over_Q']-1
 rb=root_scalar(binding,bracket=(.90649,.90659),method='brentq',xtol=2e-10,maxiter=12)
 sb,ob=solve_fixed(rb.root,nearest(rb.root),'binding_final_domain60',refine=True,domain_override=60,tol_override=1e-9)
 sb2,ob2=solve_fixed(rb.root,sb,'binding_final_domain80',refine=True,domain_override=80,tol_override=3e-10)
 thresholds={
  'turning':{'root_converged':rt.converged,'domain80':ot,'domain100':ot2,'Qbar_domain_difference':abs(ot2['Qbar']-ot['Qbar']),'Ebar_domain_difference':abs(ot2['Ebar']-ot['Ebar']),'sensitivity_evaluations':sensitivity_records},
  'binding':{'root_converged':rb.converged,'domain60':ob,'domain80':ob2,'Qbar_domain_difference':abs(ob2['Qbar']-ob['Qbar']),'Ebar_domain_difference':abs(ob2['Ebar']-ob['Ebar'])},
 }
 targets=data['targets'];Qstar=data['reference']['Qbar'];Estar=data['reference']['Ebar']
 for frac in refinement_protocol['additional_exact_charge_fractions']:
  found=[]
  for name,pred in [('lower_Omega',lambda O:O<rt.root),('higher_Omega',lambda O:O>rt.root)]:
   choices=[x for x in SOLUTIONS if pred(x[0])]
   if not choices:continue
   target=frac*Qstar;O,ss,oo=min(choices,key=lambda x:abs(x[2]['Qbar']-target))
   candidate,coarse=target_solve(target,ss,O,'additional_'+name)
   if coarse is not None:
    fine,ff=target_solve(target,candidate,coarse['Omega'],'additional_'+name+'_refined',refine=True)
    found.append({'segment':name,'coarse':coarse,'fine':ff,'energy_refinement_abs':abs(ff['Ebar']-coarse['Ebar']) if ff else None,'charge_refinement_abs':abs(ff['Qbar']-coarse['Qbar']) if ff else None})
  targets[str(frac)]={'fraction':frac,'target_Qbar':frac*Qstar,'found':found,'below_found_turning_charge':frac*Qstar<ot2['Qbar']}
 def choose(frac):
  v=min(targets.values(),key=lambda x:abs(x['fraction']-frac));vals=[f['fine'] or f['coarse'] for f in v['found']]
  return min(vals,key=lambda x:x['Ebar']) if vals else None
 extra=[]
 for frac in [.375,.4,.45]:
  s=choose(frac);E=2*s['Ebar']+(1-2*frac)*Qstar if s else None
  extra.append({'kind':'two_lumps_plus_free','lump_fraction_each':frac,'free_fraction':1-2*frac,'resolved':s is not None,'fragment':s,'Ebar':E,'delta_Ebar':E-Estar if E else None})
 # Explicit branch combinations for resolved equal and unequal splits.
 branch_combinations=[]
 for frac in [.4,.5,.375,.45]:
  v=min(targets.values(),key=lambda x:abs(x['fraction']-frac));free=0 if frac in [.4,.5] else 1-2*frac
  other=1-frac if free==0 else frac
  w=min(targets.values(),key=lambda x:abs(x['fraction']-other))
  if abs(w['fraction']-other)>1e-10:continue
  for first in v['found']:
   for second in w['found']:
    a=first['fine'] or first['coarse'];b=second['fine'] or second['coarse'];E=a['Ebar']+b['Ebar']+free*Qstar
    branch_combinations.append({'first_fraction':frac,'second_fraction':other,'free_fraction':free,'segments':[first['segment'],second['segment']],'Ebar':E,'delta_Ebar':E-Estar})
 # Fixed-Q radial Hessian on either side of the turning point.
 stability=[]
 for O in [.95,.97]:
  ss,oo=solve_fixed(O,nearest(O),'radial_check_profile',refine=True,domain_override=70,tol_override=1e-9)
  I=oo['Qbar']/O
  for N in [800,1200]:
   L=70.;h=L/(N+1);r=h*np.arange(1,N+1);F,_,S,_=ss.sol(r)
   A=diags([-np.ones(N-1)/h**2,2/h**2+S*S+.3*F*F-O*O,-np.ones(N-1)/h**2],[-1,0,1]);B=diags([-np.ones(N-1)/h**2,2/h**2+3*S*S-1+F*F,-np.ones(N-1)/h**2],[-1,0,1]);C=diags(2*F*S)
   H=bmat([[A,C],[C,B]],format='csr');q=np.r_[np.sqrt(4*np.pi)*r*F,np.zeros(N)];factor=4*O*O/I
   op=LinearOperator(H.shape,matvec=lambda x:H@x+factor*h*q*np.dot(q,x),dtype=float)
   ev=np.sort(eigsh(op,k=2,which='SA',tol=2e-9,return_eigenvectors=False))
   rec={'Omega':O,'N':N,'L':L,'h':h,'lowest_eigenvalues':ev};stability.append(rec);log({'fixedQ_radial':rec})
 data.update({'initial_thresholds':data['thresholds'],'refined_thresholds':thresholds,'targets':targets,'additional_comparisons':extra,'all_resolved_branch_combinations':branch_combinations,'branch_radial_stability':stability,'attempts':ATTEMPTS,'states':STATES,'refinement_protocol':refinement_protocol})
 add_postprocessing(data)
 dump(OUT/'branch_results.json',data)
 dump(OUT/'branch_thresholds_refined.json',thresholds)
 with (OUT/'branch_states.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(s.keys() for s in STATES))));writer.writeheader();writer.writerows(STATES)
 log({'second_pass_complete':True,'turning':ot2,'binding':ob2,'additional_comparisons':extra})

def add_postprocessing(data):
 """Algebraic comparisons only; no additional solve or parameter selection."""
 Q=data['reference']['Qbar'];E=data['reference']['Ebar'];targets=data['targets']
 allfree=[];allpairs=[]
 def found(v):return [(x['segment'],x['fine'] or x['coarse']) for x in v['found']]
 for v in targets.values():
  f=v['fraction']
  if not 0<f<1:continue
  for segment,state in found(v):
   cost=state['Ebar']+(1-f)*Q
   allfree.append({'lump_fraction':f,'segment':segment,'Ebar':cost,'delta_Ebar':cost-E,'lump_E_over_Q':state['E_over_Q']})
  if f>.5:continue
  other=min(targets.values(),key=lambda v:abs(v['fraction']-(1-f)))
  if abs(other['fraction']-(1-f))>1e-10:continue
  for aseg,a in found(v):
   for bseg,b in found(other):
    cost=a['Ebar']+b['Ebar']
    allpairs.append({'fractions':[f,1-f],'segments':[aseg,bseg],'Ebar':cost,'delta_Ebar':cost-E})
 data['all_resolved_lump_plus_free']=allfree;data['all_resolved_two_lump_pairs']=allpairs
 minimum=data['refined_thresholds']['turning']['domain100']['Qbar']
 lower_found=[s for s in data['states'] if 1e-7<s['Qbar']<=Q*(1+1e-9)]
 data['conditional_envelope'].update({
  'minimum_resolved_stationary_Qbar':minimum,
  'three_times_found_minimum_Qbar':3*minimum,
  'three_found_lumps_exceed_reference_charge':3*minimum>Q,
  'smallest_resolved_positive_partition_delta_Ebar':min([x['delta_Ebar'] for x in allfree+allpairs if x['delta_Ebar']>1e-6]),
  'minimum_resolved_lower_charge_E_over_Q':min(s['E_over_Q'] for s in lower_found),
  'conditional_statement':'For positive-charge components from the continuous branches found here with 0<q<=Qstar, if no unobserved lower-energy state intervenes, e(q)=E(q)/q >= e(Qstar). Free same-sign X carriers have e=1. Therefore any such partition sums to energy >= E(Qstar). The finite scan does not prove completeness of these branches or control arbitrary field configurations, nonspherical states, opposite-charge creation, or quantum channels.',
  'no_uniform_gap_for_arbitrarily_small_free_charge':'A lump of charge Qstar-dq plus free charge dq approaches the original energy as dq->0; first order cost is (1-Omega_star)*dq.',
 })
 curve=[s for s in data['states'] if s['label']=='continuation']
 identity=[abs(s['virial_identity_E_minus_OmegaQ_minus_2K_over3']) for s in data['states'] if s['refined']]
 data['identity_checks']={'maximum_refined_abs_E_minus_OmegaQ_minus_2K_over3':max(identity),'analytic_first_law':'dE/dQ=Omega at fixed couplings on a smooth stationary branch','analytic_e_over_q_derivative':'d(E/Q)/dQ=(Omega*Q-E)/Q^2=-2*K/(3*Q^2)<0','K_definition':'K=4pi integral rho^2 (Fprime^2+Sprime^2)/2 drho'}

if __name__=='__main__':
 if '--refine' in sys.argv:refine_followup()
 else:main()
