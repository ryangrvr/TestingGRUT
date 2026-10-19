#!/usr/bin/env python3
"""Frozen flat-space MSbar scalar Coleman-Weinberg diagnostic, v=1.

No finite counterterm/matching value is fitted. This zero-derivative constant-
background potential is not the quantum energy or lifetime of a rotating Q-ball.
Run: python effective_potential.py
"""
from __future__ import annotations
import itertools,json,platform,sys,time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import root_scalar,minimize_scalar,minimize

OUT=Path(__file__).resolve().parent
PARAM={'lambda_phi':.01,'lambda_X':.001,'c':.01,'d':.0001,'lambda_D':.001,'alpha':0.,'mphi2':-.01,'mX2':0.,'mD2':.01,'C':.0025,'mu':.1,'v':1.,'v_over_reduced_Mpl':1e-5}
PROTOCOL={
 'parameters':PARAM,'scheme':'flat-space MSbar, scalar modes only, finite constant -3/2, fixed mu=.1v',
 'potential':'C+mphi2*phi²/2+lambda_phi*phi⁴/4+mX2*(x²+y²)/2+lambda_X*(x²+y²)²/4+c*phi²*(x²+y²)/2+mD2*D²/2+lambda_D*D⁴/4+d*phi²*D²/2+alpha*D²*(x²+y²)/2',
 'field_order':['phi','x','y','D'],'background':'x=F,y=D=0, X=(x+i*y)/sqrt(2)',
 'CW':'sum eigenvalue²*(log(eigenvalue/mu²)-3/2)/(64*pi²), eigenvalue=m_a²',
 'negative_eigenvalues':'retain complex principal logarithm log(negative)=log(abs)+i*pi; report imaginary part and unstable-background flag; no clipping',
 'zero_eigenvalue':'continuous t² log(t) -> 0 for potential only; analytic curvature routines require strictly positive tree masses',
 'vacuum':{'tree':[1.,0.,0.,0.],'root_bracket_phi':[.98,1.02],'xtol':1e-14,'report':'strict order-hbar shifts and separately minimum of V0+V1'},
 'bulk':{'tree_F':10**.25,'local_rectangle_phi':[-.15,.15],'local_rectangle_F':[1.4,2.2],'grid':[31,101],'one_dimensional_xatol':1e-13,'two_dimensional_starts':[[0.,10**.25],[.1,10**.25-.15],[-.1,10**.25+.15]],'optimizer_maxiter':100,'report':'strict one-loop shift plus truncated-potential numerical minimum; no derivative terms or finite-charge energy'},
 'derivatives':{'method':'analytic spectral Frechet derivatives of trace f(H), including degeneracies','finite_difference_steps':[1e-3,3e-4,1e-4,3e-5,1e-5],'points':[[1.,0.,0.,0.],[0.,10**.25,0.,0.],[.95,.08,.04,.03]]},
 'global_background_scan':{'phi':[0,1.2,61],'F':[0,2.4,121]},
 'profile_background_scan':{'source':'../grut_candidate_1_20261018/code/results_profile.npz','rho':[0,20,1001],'interpretation':'Only flags constant-background Hessian signs on classical-profile values; no local-CW integration as quantum soliton energy'},
 'gravity':'Raw fixed-C vacuum shift reported. Subtraction for relative energies does not calibrate gravity or solve the cosmological-constant problem.',
 'scope':'Scalar constant-background one-loop potential only. No pole masses, running-scale claim, nonlocal determinant, rotating-phase derivative correction, full soliton one-loop energy, or quantum lifetime.',
 'frozen_before_execution':True,
}

# Symmetric quartic tensor: U=C+1/2 m_i² q_i²+1/24 K_ijkl q_i q_j q_k q_l.
M=np.array([PARAM['mphi2'],PARAM['mX2'],PARAM['mX2'],PARAM['mD2']])
SELF=np.array([PARAM['lambda_phi'],PARAM['lambda_X'],PARAM['lambda_X'],PARAM['lambda_D']])
PAIR=np.array([[0,PARAM['c'],PARAM['c'],PARAM['d']],[PARAM['c'],0,PARAM['lambda_X'],PARAM['alpha']],[PARAM['c'],PARAM['lambda_X'],0,PARAM['alpha']],[PARAM['d'],PARAM['alpha'],PARAM['alpha'],0.]])
K=np.zeros((4,4,4,4))
for i in range(4):K[i,i,i,i]=6*SELF[i]
for i in range(4):
 for j in range(i+1,4):
  for ix in set(itertools.permutations([i,i,j,j])):K[ix]=2*PAIR[i,j]


def serial(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
 raise TypeError(type(x).__name__)
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,default=serial)+'\n')
def log(obj):
 print(json.dumps(obj,default=serial),flush=True)
 with (OUT/'potential_execution.jsonl').open('a') as f:f.write(json.dumps(obj,default=serial)+'\n')


def tree(q):
 phi,x,y,D=np.asarray(q);r2=x*x+y*y
 return PARAM['C']+.5*PARAM['mphi2']*phi**2+.25*PARAM['lambda_phi']*phi**4+.5*PARAM['mX2']*r2+.25*PARAM['lambda_X']*r2**2+.5*PARAM['c']*phi**2*r2+.5*PARAM['mD2']*D**2+.25*PARAM['lambda_D']*D**4+.5*PARAM['d']*phi**2*D**2+.5*PARAM['alpha']*D**2*r2

def tree_grad(q):return M*q+np.einsum('iabc,a,b,c->i',K,q,q,q)/6

def tree_hessian(q):return np.diag(M)+.5*np.einsum('ijab,a,b->ij',K,q,q)


def cw(q,derivatives=False):
 q=np.asarray(q,dtype=float);H=tree_hessian(q);mass2,U=np.linalg.eigh(H)
 neg=mass2<0;nonzero=mass2!=0
 term=np.zeros(4,dtype=complex)
 term[nonzero]=mass2[nonzero]**2*(np.log(mass2[nonzero].astype(complex)/PARAM['mu']**2)-1.5)
 value=term.sum()/(64*np.pi**2)
 if not derivatives:return value,mass2
 if np.any(mass2<=0):raise ValueError('Analytic real derivatives require all four tree Hessian masses²>0')
 # g(t)=t(log(t/mu²)-1); dV1=Tr g(H)dH/(32pi²).
 logs=np.log(mass2/PARAM['mu']**2);g=mass2*(logs-1)
 T=np.einsum('ijab,b->ija',K,q)
 A=np.einsum('ki,klm,lj->ijm',U,T,U)
 B=np.einsum('ki,klmn,lj->ijmn',U,K,U)
 grad=np.einsum('i,iia->a',g,A)/(32*np.pi**2)
 divided=np.empty((4,4))
 for i in range(4):
  for j in range(4):
   if abs(mass2[i]-mass2[j])<=1e-10*max(mass2[i],mass2[j]):divided[i,j]=np.log((mass2[i]+mass2[j])/(2*PARAM['mu']**2))
   else:divided[i,j]=(g[i]-g[j])/(mass2[i]-mass2[j])
 hess=(np.einsum('i,iiab->ab',g,B)+np.einsum('ij,ija,jib->ab',divided,A,A))/(32*np.pi**2)
 return value.real,mass2,grad,hess

def effective(q,derivatives=False):
 if derivatives:
  v,m,g,h=cw(q,True);return tree(q)+v,m,tree_grad(q)+g,tree_hessian(q)+h
 v,m=cw(q);return tree(q)+v,m


def check_derivatives():
 out=[]
 for point in PROTOCOL['derivatives']['points']:
  q=np.array(point);value,m,grad,hess=cw(q,True)
  for step in PROTOCOL['derivatives']['finite_difference_steps']:
   fdg=np.zeros(4);fdh=np.zeros((4,4));fdhg=np.zeros((4,4));treefd=np.zeros((4,4))
   for a in range(4):
    ea=np.eye(4)[a]*step
    fp=cw(q+ea)[0].real;fm=cw(q-ea)[0].real
    fdg[a]=(fp-fm)/(2*step);fdh[a,a]=(fp-2*value+fm)/step**2
    fdhg[:,a]=(cw(q+ea,True)[2]-cw(q-ea,True)[2])/(2*step)
    treefd[a,a]=(tree(q+ea)-2*tree(q)+tree(q-ea))/step**2
    for b in range(a):
     eb=np.eye(4)[b]*step
     fdh[a,b]=fdh[b,a]=(cw(q+ea+eb)[0].real-cw(q+ea-eb)[0].real-cw(q-ea+eb)[0].real+cw(q-ea-eb)[0].real)/(4*step**2)
     treefd[a,b]=treefd[b,a]=(tree(q+ea+eb)-tree(q+ea-eb)-tree(q-ea+eb)+tree(q-ea-eb))/(4*step**2)
   out.append({'point':point,'step':step,'CW_gradient_max_abs_error':np.max(abs(fdg-grad)),'CW_hessian_from_values_max_abs_error':np.max(abs(fdh-hess)),'CW_hessian_from_gradients_max_abs_error':np.max(abs(fdhg-hess)),'tree_hessian_max_abs_error':np.max(abs(treefd-tree_hessian(q)))})
 return out


def main():
 start=time.perf_counter();(OUT/'potential_execution.jsonl').write_text('');dump(OUT/'potential_protocol.json',PROTOCOL)
 qtree=np.array([1.,0.,0.,0.]);V1tree,mvac,g1,h1=cw(qtree,True)
 delta=-np.linalg.solve(tree_hessian(qtree),g1)
 root=root_scalar(lambda phi:effective(np.array([phi,0.,0.,0.]),True)[2][0],bracket=(.98,1.02),xtol=1e-14)
 qvac=np.array([root.root,0.,0.,0.]);Vvac,_,gvac,Hvac=effective(qvac,True)
 strictH=tree_hessian(qtree)+h1+np.einsum('ijab,a,b->ij',K,qtree,delta)
 vacuum={
  'field_order':['phi','x','y','D'],'tree_point':qtree,'tree_curvature_diagonal':np.diag(tree_hessian(qtree)),'tree_energy':tree(qtree),'tree_masses_squared':mvac,'CW_at_tree_point':V1tree,'CW_tadpole_at_tree_point':g1,'CW_curvature_correction_at_tree_point':h1,
  'strict_one_loop_field_shift':delta,'strict_one_loop_vacuum_point':qtree+delta,'strict_one_loop_curvature_masses_squared':np.linalg.eigvalsh(strictH),'strict_one_loop_curvature_diagonal':np.diag(strictH),'strict_one_loop_vacuum_energy':tree(qtree)+V1tree,
  'minimum_of_truncated_V0_plus_V1':qvac,'truncated_vacuum_energy':Vvac,'truncated_vacuum_gradient':gvac,'truncated_curvature_masses_squared':np.linalg.eigvalsh(Hvac),'truncated_curvature_diagonal':np.diag(Hvac),'truncated_curvature_masses':np.sqrt(np.diag(Hvac)),
  'curvature_relative_shift_vs_tree':np.diag(Hvac)/np.diag(tree_hessian(qtree))-1,
  'truncated_minus_strict_vacuum_energy':Vvac-(tree(qtree)+V1tree),
  'field_shift_relative':qvac[0]-1,
  'note':'Curvatures of the constant-field effective potential, not momentum-dependent pole masses. Numerical minimum of V0+V1 contains induced higher-order pieces; it is not a two-loop result.',
 }
 log({'vacuum':vacuum})
 # Strict first-order bulk change and independent minimum of truncated potential.
 Ft=10**.25;qt=np.array([0.,Ft,0.,0.]);V1b,mbtree,gb,hb=cw(qt,True)
 h0=2*(tree(qt)-tree(qtree))/Ft**2;h1bulk=2*(V1b-V1tree)/Ft**2
 h1prime=2*gb[1]/Ft**2-4*(V1b-V1tree)/Ft**3;deltaF=-h1prime/(4*PARAM['lambda_X'])
 def bulk_value_grad(z):
  phi,F=z;q=np.array([phi,F,0.,0.]);V,m,g,h=effective(q,True);dv=V-Vvac
  value=2*dv/F**2
  grad=np.array([2*g[0]/F**2,2*g[1]/F**2-4*dv/F**3])
  hh=np.array([[2*h[0,0]/F**2,2*h[0,1]/F**2-4*g[0]/F**3],[2*h[0,1]/F**2-4*g[0]/F**3,2*h[1,1]/F**2-8*g[1]/F**3+12*dv/F**4]])
  return value,grad,hh
 opt=minimize_scalar(lambda F:bulk_value_grad([0,F])[0],bounds=(1.4,2.2),method='bounded',options={'xatol':1e-13,'maxiter':100})
 opts=[]
 for point in PROTOCOL['bulk']['two_dimensional_starts']:
  result=minimize(lambda z:bulk_value_grad(z)[0],point,jac=lambda z:bulk_value_grad(z)[1],bounds=[(-.15,.15),(1.4,2.2)],method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-12,'maxiter':100})
  opts.append({'start':point,'success':result.success,'message':str(result.message),'point':result.x,'value':result.fun,'gradient':bulk_value_grad(result.x)[1],'evaluations':result.nfev})
 min_tree_mass=np.inf;min_eff_curv=np.inf;min_ratio_curv=np.inf;max_imag=0.
 for phi in np.linspace(-.15,.15,31):
  for F in np.linspace(1.4,2.2,101):
   vv,mm,gg,hh=effective(np.array([phi,F,0.,0.]),True)
   min_tree_mass=min(min_tree_mass,mm.min());min_eff_curv=min(min_eff_curv,np.linalg.eigvalsh(hh).min());min_ratio_curv=min(min_ratio_curv,np.linalg.eigvalsh(bulk_value_grad([phi,F])[2]).min())
   max_imag=max(max_imag,abs(cw([phi,F,0,0])[0].imag))
 vb,mbmin,gbr,hbr=effective(np.array([0,opt.x,0.,0.]),True)
 bulk={
  'tree_F':Ft,'tree_ratio_mass_squared':h0,'tree_CW':V1b,'tree_background_masses_squared':mbtree,'truncated_minimum_background_tree_Hessian_masses_squared':mbmin,
  'strict_one_loop_ratio_correction':h1bulk,'strict_one_loop_ratio_relative_correction':h1bulk/h0,'strict_one_loop_ratio_mass_squared':h0+h1bulk,'strict_one_loop_F_shift':deltaF,'strict_one_loop_F':Ft+deltaF,
  'truncated_minimum_F':opt.x,'truncated_minimum_phi':0.,'truncated_minimum_ratio_mass_squared':opt.fun,'truncated_minimum_gradient':bulk_value_grad([0,opt.x])[1],'truncated_minimum_ratio_hessian':bulk_value_grad([0,opt.x])[2],
  'bulk_ratio_relative_correction':opt.fun/h0-1,'bulk_frequency_over_tree_X_mass':np.sqrt(opt.fun/.01),'bulk_frequency_over_truncated_vacuum_X_curvature_mass':np.sqrt(opt.fun/Hvac[1,1]),
  'one_dimensional_optimizer_success':opt.success,'two_dimensional_optimizer_checks':opts,
  'rectangle_min_tree_mass_squared_sampled':min_tree_mass,'rectangle_min_effective_curvature_squared_sampled':min_eff_curv,'rectangle_min_ratio_hessian_sampled':min_ratio_curv,'rectangle_max_CW_imaginary':max_imag,
  'rectangle_tree_convexity_certificate':{'A_lower':.0096,'B_lower':.00588,'offdiag_abs_upper':.0066,'determinant_lower':.0096*.00588-.0066**2,'phase_mass_lower':.00196,'D_mass_lower':.01},
  'note':'This is a homogeneous zero-derivative ratio diagnostic, not a finite-charge quantum energy or existence/stability proof.',
 }
 log({'bulk':bulk})
 checks=check_derivatives();log({'derivative_checks':checks})
 # Global-background map: complex values retained and flagged, never clipped.
 rows=[]
 for phi in np.linspace(0,1.2,61):
  for F in np.linspace(0,2.4,121):
   value,mm=cw([phi,F,0,0]);rows.append([phi,F,*mm,value.real,value.imag,int(np.sum(mm<0))])
 arr=np.array(rows);np.savetxt(OUT/'potential_background_scan.csv',arr,delimiter=',',header='phi,F,mass2_0,mass2_1,mass2_2,mass2_3,CW_real,CW_imaginary,negative_mode_count',comments='')
 neg=arr[:,-1]>0
 scan={'points':len(arr),'negative_Hessian_points':int(neg.sum()),'fraction_negative':neg.mean(),'minimum_mass_squared':arr[:,2:6].min(),'maximum_absolute_CW_imaginary':np.max(abs(arr[:,-2])),'negative_points_are_excluded_from_real_convex_minimizations':True,'principal_log_imaginary_sign_convention':'+i*pi; imaginary magnitude diagnoses unstable constant backgrounds, not a computed Q-ball lifetime'}
 # Flag, but do not integrate a local CW energy over the rotating profile.
 ref=OUT.parent/'grut_candidate_1_20261018/code/results_profile.npz';p=np.load(ref);spl=scipy.interpolate.CubicSpline(p['rho'],p['y'],axis=1);rows=[]
 for rho in np.linspace(0,20,1001):
  F,_,phi,_=spl(rho);value,mm=cw([phi,F,0,0]);rows.append([rho,phi,F,*mm,value.real,value.imag,int(np.sum(mm<0))])
 profile=np.array(rows);np.savetxt(OUT/'potential_profile_flags.csv',profile,delimiter=',',header='rho,phi,F,mass2_0,mass2_1,mass2_2,mass2_3,CW_real,CW_imaginary,negative_mode_count',comments='')
 pn=profile[:,-1]>0
 profile_flags={'points':len(profile),'negative_Hessian_points':int(pn.sum()),'rho_range_negative':[profile[pn,0].min(),profile[pn,0].max()] if pn.any() else None,'minimum_mass_squared':profile[:,3:7].min(),'maximum_absolute_CW_imaginary':np.max(abs(profile[:,-2])),'CW_profile_energy_integral_performed':False,'interpretation':'Constant-background Hessian evaluated on profile values only; derivative and nonlocal terms are absent.'}
 # Constant subtraction only defines relative matter energies. Gravity needs its own renormalized matching condition.
 rho_strict=vacuum['strict_one_loop_vacuum_energy'];rho_trunc=vacuum['truncated_vacuum_energy'];ratio=PARAM['v_over_reduced_Mpl']
 gravity={'fixed_C':PARAM['C'],'raw_vacuum_energy_strict_one_loop_over_v4':rho_strict,'raw_vacuum_energy_truncated_minimum_over_v4':rho_trunc,'required_constant_shift_for_zero_vacuum_strict_over_v4':-rho_strict,'required_constant_shift_for_zero_vacuum_truncated_over_v4':-rho_trunc,'strict_shift_relative_to_tree_C':-rho_strict/PARAM['C'],'strict_raw_vacuum_energy_over_reduced_Mpl4':rho_strict*ratio**4,'strict_induced_Lambda_over_reduced_Mpl2':rho_strict*ratio**4,'strict_induced_R_over_reduced_Mpl2':4*rho_strict*ratio**4,'strict_AdS_radius_times_v':np.sqrt(3/(abs(rho_strict)*ratio**2)),'gravity_calibration_applied':False,'note':'Subtracting Veff(vacuum) in the bulk ratio fixes a relative-energy reference only. It does not set the gravitating cosmological constant; no gravity matching value was fitted here.'}
 result={'protocol':PROTOCOL,'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'command':'python effective_potential.py','vacuum':vacuum,'bulk':bulk,'derivative_checks':checks,'background_scan':scan,'classical_profile_flags':profile_flags,'gravity_vacuum_shift':gravity,'elapsed_seconds':time.perf_counter()-start}
 dump(OUT/'potential_results.json',result);log({'background_scan':scan,'classical_profile_flags':profile_flags,'gravity_vacuum_shift':gravity,'complete':True,'elapsed_seconds':result['elapsed_seconds']})

if __name__=='__main__':main()
