#!/usr/bin/env python3
"""Bounded tree rotating Q-ball fluctuation spectrum. No quantum energy sum.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python modes.py
"""
from __future__ import annotations
import csv,hashlib,json,platform,sys,time
from pathlib import Path
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from scipy.linalg import eig,eigvalsh,eigh_tridiagonal
from scipy.sparse import diags,bmat,eye,csr_matrix
from scipy.sparse.linalg import spsolve,eigs,eigsh,ArpackNoConvergence,LinearOperator

OUT=Path(__file__).resolve().parent
REFERENCE=OUT.parent/'grut_candidate_1_20261018/code/results_profile.npz'
PROTOCOL={
 'action':'Same original four-scalar tree action and matched couplings. a=1,b=.1,g=.1,d/g²=.01,mD0²/(g²v²)=1,alpha0; D background zero.',
 'time_radius_units':'tau=g*v*t, rho=g*v*r; s generator exponent e^(s*tau), frequency nu=Im(s) for stable positive-frequency roots',
 'fields':'q=(rho*u,rho*w,rho*eta); qddot+G*qdot+K*q=0; G_uw=+2Omega,G_wu=-2Omega',
 'reference':str(REFERENCE),'Omega':'read unchanged from reference profile',
 'boundaries':'uniform interior rho=jL/(N+1); Dirichlet transformed fields at rho=0,L; S(L)=1',
 'background':'consistent fixed-Omega finite-difference stationary p=rho F,z=rho(S-1), Newton max15 iterations, max residual1e-10; no physical fit; retain failures',
 'dense_grids':[{'N':120,'L':30.},{'N':180,'L':30.}],
 'sparse_grids':[{'N':400,'L':30.},{'N':800,'L':30.},{'N':800,'L':40.}],
 'ell':[0,1,2,3,4],
 'sparse_shifts':[[0.,.08],[0.,.18],[.05,0.]],'sparse_k':18,'sparse_tol':1e-10,'sparse_maxiter':1500,
 'vacuum_check':{'N':80,'L':30.,'ell':[0,2],'expected':'sqrt(k_discrete²+1)±Omega, sqrt(k_discrete²+2), sqrt(k_discrete²+1.01)'},
 'controls':'full dense spectra at two spacings; sparse low-frequency/real-axis probes at two spacings and one larger box; phase and translation overlap/residuals; exact discrete charge-left-null check; amplitude constrained and phase energy eigenvalues',
 'internal_mode_window':[1e-5,.2],'all_positive_real_roots_retained':True,
 'higher_ell_scope':'K_ell=K_2+[ell(ell+1)-6]/rho² for ell>=2; numerical K2 positivity plus this algebraic comparison supports higher-ell energy test without infinite scan; continuum positivity is conditional on converged background/operator.',
 'budget_seconds':240,'time_budget_behavior':'do not start another grid once elapsed budget exceeded; record skipped grids; no paid compute or installs',
 'no_claims':['complete continuum spectrum','renormalized determinant or quantum energy','full soliton one-loop correction','quantum lifetime'],
 'frozen_before_execution':True,
}

def serial(o):
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 if isinstance(o,complex):return {'real':o.real,'imag':o.imag}
 raise TypeError(type(o).__name__)
def dump(path,data):path.write_text(json.dumps(data,indent=2,default=serial)+'\n')
def log(data):
 print(json.dumps(data,default=serial),flush=True)
 with (OUT/'mode_execution.jsonl').open('a') as f:f.write(json.dumps(data,default=serial)+'\n')

ROWS=[];SELECTED={};ATTEMPTS=[]

def background(N,L,O,spl,tol=1e-10):
 h=L/(N+1);r=h*np.arange(1,N+1);F,Fp,S,Sp=spl(r);p=r*F;z=r*(S-1)
 minusD=diags([-np.ones(N-1)/h**2,2*np.ones(N)/h**2,-np.ones(N-1)/h**2],[-1,0,1],format='csr')
 def residual(p,z):
  F=p/r;S=1+z/r
  return np.r_[minusD@p+(S*S+.1*F*F-O*O)*p,minusD@z+((S*S-1)+F*F)*S*r]
 start=time.perf_counter();history=[];success=False
 for it in range(15):
  R=residual(p,z);norm=np.max(abs(R));history.append(norm)
  if norm<tol:success=True;break
  F=p/r;S=1+z/r
  A=minusD+diags(S*S+.3*F*F-O*O);B=minusD+diags(3*S*S-1+F*F);C=diags(2*F*S)
  step=spsolve(bmat([[A,C],[C,B]],format='csr'),-R)
  # Fixed, bounded Newton line search; no profile architecture tuning.
  for weight in [1.,.5,.25,.125,.0625]:
   pn=p+weight*step[:N];zn=z+weight*step[N:]
   if np.max(abs(residual(pn,zn)))<norm:p,z=pn,zn;break
  else:break
 F=p/r;S=1+z/r
 rec={'operation':'discrete_background','N':N,'L':L,'h':h,'success':success,'residual_history':history,'max_residual':np.max(abs(residual(p,z))),'target_tolerance':tol,'minF':F.min(),'minS':S.min(),'Qbar_discrete':4*np.pi*h*O*np.dot(p,p),'seconds':time.perf_counter()-start}
 ATTEMPTS.append(rec);log(rec)
 if not success or F.min()<-1e-9 or S.min()<0:raise RuntimeError('Discrete background failed acceptance')
 return {'N':N,'L':L,'h':h,'r':r,'F':F,'S':S,'p':p,'z':z,'D':minusD,'record':rec}


def operators(bg,ell,O,vacuum=False):
 N=bg['N'];r=bg['r'];F=np.zeros(N) if vacuum else bg['F'];S=np.ones(N) if vacuum else bg['S']
 lap=bg['D']+diags(ell*(ell+1)/(r*r))
 Lu=lap+diags(S*S+.3*F*F-O*O);Lw=lap+diags(S*S+.1*F*F-O*O);Le=lap+diags(3*S*S-1+F*F);LD=lap+diags(1+.01*S*S)
 C=diags(2*F*S);zero=csr_matrix((N,N));I=eye(N,format='csr')
 K=bmat([[Lu,zero,C],[zero,Lw,zero],[C,zero,Le]],format='csr')
 G=bmat([[zero,2*O*I,zero],[-2*O*I,zero,zero],[zero,zero,zero]],format='csr')
 Z=csr_matrix((3*N,3*N));A=bmat([[Z,eye(3*N)],[-K,-G]],format='csr')
 return A,K,G,Lu,Lw,Le,LD,C


def energy_checks(bg,ell,O,ops):
 A,K,G,Lu,Lw,Le,LD,C=ops;N=bg['N'];h=bg['h'];p=bg['p'];H=bmat([[Lu,C],[C,Le]],format='csr')
 if ell==0:
  qq=np.r_[p,np.zeros(N)];fac=4*O*O/(h*np.dot(p,p));op=LinearOperator(H.shape,matvec=lambda x:H@x+fac*h*qq*np.dot(qq,x),dtype=float)
 else:op=H
 e=np.sort(eigsh(op,k=2,which='SA',tol=2e-10,return_eigenvectors=False))
 phase=eigh_tridiagonal(Lw.diagonal(),-np.ones(N-1)/h**2,select='i',select_range=(0,2),check_finite=False)[0]
 dm=eigh_tridiagonal(LD.diagonal(),-np.ones(N-1)/h**2,select='i',select_range=(0,2),check_finite=False)[0]
 rec={'real_amplitude_lowest_energy_eigenvalues':e,'phase_lowest_energy_eigenvalues':phase,'D_lowest_frequency':np.sqrt(dm),'rank_one_fixedQ':ell==0}
 if ell==0:
  # Exact discrete left-null charge and phase vector, plus its charge-changing Jordan partner.
  phasevec=np.r_[np.zeros(N),p,np.zeros(4*N)]
  charge=np.r_[2*O*p,np.zeros(3*N),-p,np.zeros(N)]
  rec['phase_generator_residual']=np.linalg.norm(A@phasevec)/np.linalg.norm(phasevec)
  rec['charge_left_null_residual']=np.linalg.norm(A.T@charge)/np.linalg.norm(charge)
  derivative=spsolve(H,np.r_[2*O*p,np.zeros(N)])
  qgen=np.r_[-derivative[:N],np.zeros(N),-derivative[N:]];vgen=phasevec[:3*N];gen=np.r_[qgen,vgen]
  rec['phase_Jordan_residual']=np.linalg.norm(A@gen-phasevec)/np.linalg.norm(phasevec)
  rec['dQbar_dOmega_discrete']=4*np.pi*h*(np.dot(p,p)+2*O*np.dot(p,derivative[:N]))
  rec['Jordan_partner_deltaQbar']=4*np.pi*h*np.dot(charge,gen)
 return rec


def classify(values,vectors,bg,ell,O,label,ops):
 N=bg['N'];h=bg['h'];r=bg['r'];p=bg['p'];A=ops[0]
 qphase=np.r_[np.zeros(N),p,np.zeros(N)];qphase/=np.linalg.norm(qphase)
 # Translation uses original smooth derivatives; consistent FD backgrounds converge to it.
 fp=SOURCE_SPLINE(r)[1];sp=SOURCE_SPLINE(r)[3]
 qt=np.r_[r*fp,np.zeros(N),r*sp];qt/=np.linalg.norm(qt)
 charge=np.r_[2*O*p,np.zeros(3*N),-p,np.zeros(N)];charge/=np.linalg.norm(charge)
 records=[];internal=[]
 for j,s in enumerate(values):
  v=vectors[:,j];q=v[:3*N];qn=np.linalg.norm(q)
  phaseover=abs(np.vdot(qphase,q))/qn if qn else 0.
  transover=abs(np.vdot(qt,q))/qn if qn else 0.
  symmetry='phase' if ell==0 and abs(s)<.03 and phaseover>.98 else ('translation' if ell==1 and abs(s)<.08 and transover>.90 else None)
  residual=np.linalg.norm(A@v-s*v)/((1+abs(s))*np.linalg.norm(v))
  dq=abs(np.vdot(charge,v))/np.linalg.norm(v) if ell==0 else 0.
  rec={'method':label,'N':N,'L':bg['L'],'ell':ell,'s_real':s.real,'s_imag':s.imag,'abs_s':abs(s),'symmetry':symmetry,'phase_overlap':phaseover,'translation_overlap':transover,'relative_generator_residual':residual,'deltaQ_normalized':dq,'deltaQ_zero_reason':'left-null charge test' if ell==0 else 'angular integral of Y_lm vanishes'}
  records.append(rec);ROWS.append(rec)
  if abs(s.real)<2e-6 and 1e-5<s.imag<.2-2e-6 and symmetry is None:
   internal.append(rec)
   key=f'{label}_N{N}_L{int(bg["L"])}_ell{ell}_nu{s.imag:.9f}'
   SELECTED[key+'_mode']=v;SELECTED[key+'_rho']=r
 positive=[x for x in records if x['s_real']>1e-7]
 nonzero=[x for x in records if x['abs_s']>1e-4 and x['symmetry']!='phase']
 summary={'N':N,'L':bg['L'],'ell':ell,'method':label,'eigenvalues_returned':len(values),'internal_positive_frequencies':[x['s_imag'] for x in sorted(internal,key=lambda x:x['s_imag'])],'positive_real_roots':sorted(positive,key=lambda x:-x['s_real'])[:12],'max_real_part':max(x['s_real'] for x in records),'max_real_part_excluding_symmetry':max([abs(x['s_real']) for x in records if x['symmetry'] is None],default=0.),'max_deltaQ_normalized_nonzero':max([x['deltaQ_normalized'] for x in nonzero],default=0.),'max_generator_residual':max(x['relative_generator_residual'] for x in records),'smallest_roots':sorted(records,key=lambda x:x['abs_s'])[:8]}
 return summary


def vacuum_check(O):
 N=80;L=30.;h=L/(N+1);r=h*np.arange(1,N+1);D=diags([-np.ones(N-1)/h**2,2*np.ones(N)/h**2,-np.ones(N-1)/h**2],[-1,0,1],format='csr');bg={'N':N,'L':L,'h':h,'r':r,'D':D}
 ans=[]
 for ell in [0,2]:
  ops=operators(bg,ell,O,vacuum=True);A=ops[0]
  vals=eig(A.toarray(),left=False,right=False,check_finite=False)
  actual=np.sort(vals.imag[vals.imag>1e-8]);k2=eigh_tridiagonal(2*np.ones(N)/h**2+ell*(ell+1)/(r*r),-np.ones(N-1)/h**2,check_finite=False)[0]
  expected=np.sort(np.r_[np.sqrt(k2+1)-O,np.sqrt(k2+1)+O,np.sqrt(k2+2)])
  dd=eigh_tridiagonal(ops[6].diagonal(),-np.ones(N-1)/h**2,check_finite=False)[0]
  rec={'ell':ell,'N':N,'L':L,'maximum_charged_eta_frequency_difference':np.max(abs(actual-expected)),'maximum_D_frequency_difference':np.max(abs(np.sqrt(dd)-np.sqrt(k2+1.01))),'maximum_absolute_generator_real_part':np.max(abs(vals.real)),'first_charged_lower_frequencies':np.sqrt(k2[:3]+1)-O,'first_charged_upper_frequencies':np.sqrt(k2[:3]+1)+O,'first_eta_frequencies':np.sqrt(k2[:3]+2),'first_D_frequencies':np.sqrt(k2[:3]+1.01)}
  ans.append(rec);log({'vacuum_check':rec})
 return ans


def save(result):
 result['attempts']=ATTEMPTS;dump(OUT/'mode_results.json',result)
 if ROWS:
  with (OUT/'mode_eigenvalues.csv').open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(ROWS[0]));writer.writeheader();writer.writerows(ROWS)
 np.savez(OUT/'mode_internal_vectors.npz',**SELECTED)


def main():
 global SOURCE_SPLINE
 start=time.perf_counter();(OUT/'mode_execution.jsonl').write_text('')
 raw=REFERENCE.read_bytes();p=np.load(REFERENCE);O=float(p['Omega']);SOURCE_SPLINE=CubicSpline(p['rho'],p['y'],axis=1)
 PROTOCOL['reference_sha256']=hashlib.sha256(raw).hexdigest();PROTOCOL['Omega']=O;dump(OUT/'mode_protocol.json',PROTOCOL)
 result={'protocol':PROTOCOL,'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'vacuum_checks':vacuum_check(O),'spectra':[],'skipped':[]}
 jobs=[('dense',g) for g in PROTOCOL['dense_grids']]+[('sparse',g) for g in PROTOCOL['sparse_grids']]
 for method,grid in jobs:
  if time.perf_counter()-start>240:result['skipped'].append({'method':method,**grid,'reason':'frozen time budget reached'});continue
  bg=background(grid['N'],grid['L'],O,SOURCE_SPLINE)
  for ell in PROTOCOL['ell']:
   if time.perf_counter()-start>240:result['skipped'].append({'method':method,**grid,'ell':ell,'reason':'frozen time budget reached'});continue
   ops=operators(bg,ell,O);A=ops[0];t0=time.perf_counter();energy=energy_checks(bg,ell,O,ops)
   if method=='dense':
    values,vectors=eig(A.toarray(),left=False,right=True,check_finite=False)
   else:
    allvals=[];allvec=[]
    for sr,si in PROTOCOL['sparse_shifts']:
     sigma=complex(sr,si)
     try:
      val,vec=eigs(A.astype(complex),k=18,sigma=sigma,which='LM',tol=1e-10,maxiter=1500)
      status='converged'
     except ArpackNoConvergence as exc:
      val=exc.eigenvalues;vec=exc.eigenvectors;status='partial_ARPACK_nonconvergence'
     except Exception as exc:
      ATTEMPTS.append({'operation':'sparse_shift','N':bg['N'],'L':bg['L'],'ell':ell,'sigma':sigma,'status':'exception','message':str(exc)});log(ATTEMPTS[-1]);continue
     ATTEMPTS.append({'operation':'sparse_shift','N':bg['N'],'L':bg['L'],'ell':ell,'sigma':sigma,'status':status,'count':len(val),'returned_eigenvalues':val});log(ATTEMPTS[-1])
     for i,s in enumerate(val):
      if all(abs(s-a)>2e-7 for a in allvals):allvals.append(s);allvec.append(vec[:,i])
    if not allvals:
     result['spectra'].append({'method':method,**grid,'ell':ell,'error':'no sparse eigenvalues returned'});save(result);continue
    values=np.array(allvals);vectors=np.column_stack(allvec)
   summary=classify(values,vectors,bg,ell,O,method,ops);summary['energy_checks']=energy;summary['seconds']=time.perf_counter()-t0
   result['spectra'].append(summary);log({'spectrum_summary':summary});save(result)
 result['elapsed_seconds']=time.perf_counter()-start;result['completed']=True;save(result);log({'complete':True,'elapsed_seconds':result['elapsed_seconds'],'skipped':result['skipped']})



def matched_box_check():
 global SOURCE_SPLINE,ATTEMPTS,ROWS,SELECTED
 import shutil
 start=time.perf_counter()
 amendment={'reason':'Separate outer-domain change from radial spacing; tighten a phase-zero diagnostic while preserving original roots.','matched_h_grids':[{'N':599,'L':30.},{'N':799,'L':40.}],'h':.05,'ell':[0,1,2,3,4],'sparse_shifts':[[0,.08],[0,.18],[.05,0]],'sparse_k':18,'maxiter':1500,'tol':1e-10,'phase_tightening':{'N':180,'L':30.,'background_tolerance':2e-13,'ell':0,'method':'full dense spectrum'},'extra_unconstrained_energy_diagnostic':{'N':799,'L':40.,'ell':0},'additional_budget_seconds':60,'frozen_before_execution':True,'no_physical_parameters_changed':True}
 dump(OUT/'mode_protocol_amendment_01.json',amendment)
 for name in ['mode_results.json','mode_execution.jsonl']:
  dest=OUT/name.replace('mode_','mode_initial_',1)
  if not dest.exists():shutil.copy2(OUT/name,dest)
 result=json.loads((OUT/'mode_results.json').read_text());ATTEMPTS=result['attempts']
 with (OUT/'mode_eigenvalues.csv').open() as f:ROWS=list(csv.DictReader(f))
 if (OUT/'mode_internal_vectors.npz').exists():
  pp=np.load(OUT/'mode_internal_vectors.npz');SELECTED={k:pp[k] for k in pp.files}
 pp=np.load(REFERENCE);O=float(pp['Omega']);SOURCE_SPLINE=CubicSpline(pp['rho'],pp['y'],axis=1)
 for grid in amendment['matched_h_grids']:
  bg=background(grid['N'],grid['L'],O,SOURCE_SPLINE)
  for ell in amendment['ell']:
   ops=operators(bg,ell,O);A=ops[0];t=time.perf_counter();energy=energy_checks(bg,ell,O,ops);values=[];vectors=[]
   for sr,si in amendment['sparse_shifts']:
    try:val,vec=eigs(A.astype(complex),k=18,sigma=complex(sr,si),which='LM',tol=1e-10,maxiter=1500);status='converged'
    except ArpackNoConvergence as exc:val=exc.eigenvalues;vec=exc.eigenvectors;status='partial_ARPACK_nonconvergence'
    ATTEMPTS.append({'operation':'matched_box_sparse','N':bg['N'],'L':bg['L'],'ell':ell,'sigma':complex(sr,si),'status':status,'count':len(val),'returned_eigenvalues':val})
    for j,s in enumerate(val):
     if all(abs(s-v)>2e-7 for v in values):values.append(s);vectors.append(vec[:,j])
   summary=classify(np.array(values),np.column_stack(vectors),bg,ell,O,'matched_box',ops);summary['energy_checks']=energy;summary['seconds']=time.perf_counter()-t
   result['spectra'].append(summary);log({'matched_box_summary':summary})
   if grid['N']==799 and ell==0:
    H=bmat([[ops[3],ops[7]],[ops[7],ops[5]]],format='csr')
    energy['unconstrained_amplitude_eigenvalues']=np.sort(eigsh(H,k=3,which='SA',tol=1e-10,return_eigenvectors=False))
   save(result)
 bg=background(180,30.,O,SOURCE_SPLINE,tol=2e-13);ops=operators(bg,0,O);t=time.perf_counter();val,vec=eig(ops[0].toarray(),left=False,right=True,check_finite=False)
 summary=classify(val,vec,bg,0,O,'dense_phase_tight',ops);summary['energy_checks']=energy_checks(bg,0,O,ops);summary['seconds']=time.perf_counter()-t
 result['spectra'].append(summary);log({'phase_tight_summary':summary})
 result['amendment']=amendment;result['amendment_elapsed_seconds']=time.perf_counter()-start
 result['summary']=summarize(result);save(result);log({'amendment_complete':True,'elapsed_seconds':result['amendment_elapsed_seconds'],'summary':result['summary']})


def summarize(result):
 compact=[]
 for s in result['spectra']:
  positive=[x for x in s['smallest_roots'] if x['s_imag']>1e-5 and x['symmetry'] is None]
  sym=[x for x in s['smallest_roots'] if x['symmetry'] is not None]
  compact.append({'method':s['method'],'N':s['N'],'L':s['L'],'ell':s['ell'],'lowest_positive_non_symmetry_frequency':min([x['s_imag'] for x in positive],default=None),'symmetry_roots':sym,'max_abs_real_non_symmetry':s['max_real_part_excluding_symmetry'],'fixedQ_deltaQ_nonzero_normalized_max':s['max_deltaQ_normalized_nonzero']})
 translation=[s for s in result['spectra'] if s['ell']==1 and s['L']==30. and s['method'] in ['sparse','matched_box']]
 xx=[];yy=[]
 for s in translation:
  modes=[x for x in s['smallest_roots'] if x['symmetry']=='translation' and x['s_imag']>0]
  if modes:xx.append((30/(s['N']+1))**2);yy.append(modes[0]['s_imag']**2)
 fit=np.polyfit(xx,yy,1)
 return {'compact_spectra':compact,'non_symmetry_internal_modes_found':sum(len(s['internal_positive_frequencies']) for s in result['spectra']),'max_abs_real_part_non_symmetry_all_returned_modes':max(s['max_real_part_excluding_symmetry'] for s in result['spectra']),'max_normalized_deltaQ_nonzero':max(s['max_deltaQ_normalized_nonzero'] for s in result['spectra']),'translation_nu_squared_fit_intercept':fit[1],'translation_nu_squared_fit_coefficient_h_squared':fit[0],'vacuum_continuum_thresholds':{'charged_lower':1-result['protocol']['Omega'],'charged_upper':1+result['protocol']['Omega'],'eta':np.sqrt(2),'D':np.sqrt(1.01)},'limitations':'Only dense finite grids and specified sparse probes, no full continuum spectral certification or quantum determinant.'}



def ell5_check():
 global SOURCE_SPLINE,ATTEMPTS,ROWS,SELECTED
 start=time.perf_counter()
 protocol={'reason':'Root numerical angular comparison B<31 conditionally excludes real subgap modes for ell>=6. Fill ell5 explicitly, without claiming a certified half-line bound.','ell':[5],'grids':[{'method':'dense_ell5','N':180,'L':30.},{'method':'sparse_ell5','N':800,'L':30.},{'method':'matched_ell5','N':599,'L':30.},{'method':'matched_ell5','N':799,'L':40.}],'shifts':[[0,.08],[0,.18],[.05,0]],'k':18,'tol':1e-10,'maxiter':1500,'eigenvalue_summary_dedup_distance':2e-7,'raw_shift_eigenvalues_retained':True,'budget_seconds':45,'frozen_before_execution':True,'no_physical_parameters_changed':True}
 dump(OUT/'mode_protocol_amendment_02.json',protocol)
 result=json.loads((OUT/'mode_results.json').read_text());ATTEMPTS=result['attempts']
 with (OUT/'mode_eigenvalues.csv').open() as f:ROWS=list(csv.DictReader(f))
 pp=np.load(OUT/'mode_internal_vectors.npz');SELECTED={k:pp[k] for k in pp.files}
 pp=np.load(REFERENCE);O=float(pp['Omega']);SOURCE_SPLINE=CubicSpline(pp['rho'],pp['y'],axis=1)
 for grid in protocol['grids']:
  bg=background(grid['N'],grid['L'],O,SOURCE_SPLINE);ops=operators(bg,5,O);A=ops[0];t=time.perf_counter();energy=energy_checks(bg,5,O,ops)
  if grid['method']=='dense_ell5':values,vectors=eig(A.toarray(),left=False,right=True,check_finite=False)
  else:
   values=[];vectors=[]
   for sr,si in protocol['shifts']:
    try:val,vec=eigs(A.astype(complex),k=18,sigma=complex(sr,si),which='LM',tol=1e-10,maxiter=1500);status='converged'
    except ArpackNoConvergence as exc:val=exc.eigenvalues;vec=exc.eigenvectors;status='partial_ARPACK_nonconvergence'
    rec={'operation':'ell5_sparse_shift','N':bg['N'],'L':bg['L'],'ell':5,'sigma':complex(sr,si),'status':status,'count':len(val),'returned_eigenvalues':val};ATTEMPTS.append(rec);log(rec)
    for j,s in enumerate(val):
     if all(abs(s-v)>2e-7 for v in values):values.append(s);vectors.append(vec[:,j])
   values=np.array(values);vectors=np.column_stack(vectors)
  summary=classify(values,vectors,bg,5,O,grid['method'],ops);summary['energy_checks']=energy;summary['seconds']=time.perf_counter()-t
  result['spectra'].append(summary);log({'ell5_summary':summary});save(result)
 result['ell5_amendment']=protocol;result['ell5_elapsed_seconds']=time.perf_counter()-start;result['summary']=summarize(result);save(result)
 log({'ell5_complete':True,'elapsed_seconds':result['ell5_elapsed_seconds'],'non_symmetry_internal_modes_found':result['summary']['non_symmetry_internal_modes_found'],'max_abs_real_non_symmetry':result['summary']['max_abs_real_part_non_symmetry_all_returned_modes']})

if __name__=='__main__':
 if '--box-check' in sys.argv:matched_box_check()
 elif '--ell5-check' in sys.argv:ell5_check()
 else:main()
