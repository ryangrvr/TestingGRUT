"""Exact integrity, raw accounting and optional same-author re-execution."""
from pathlib import Path
import argparse,hashlib,io,json,math,os,shutil,subprocess,sys,tarfile,tempfile
import numpy as np
P=Path(__file__).resolve().parents[1]
def need(ok,message):
 if not ok:raise AssertionError(message)
def leaves(x,p=()):
 if isinstance(x,dict):return {k:v for a,b in x.items() for k,v in leaves(b,p+(a,)).items()}
 if isinstance(x,list):return {k:v for a,b in enumerate(x) for k,v in leaves(b,p+(a,)).items()}
 return {p:x}
def projection(r):
 return dict(analytical=r['analytical'],full=[{k:v for k,v in q.items() if k not in ['seconds','nfev','max_Friedmann_relative','max_scalar_work_relative']} for q in r['full_action_runs']],late=[{k:v for k,v in q.items() if k not in ['numerical_nfev','max_scalar_work_relative']} for q in r['conditional_late_runs']+r['conditional_late_refinements']])
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--rerun',action='store_true');args=parser.parse_args()
 manifest=json.loads((P/'MANIFEST.json').read_text())['files']
 actual={str(f.relative_to(P)) for f in P.rglob('*') if f.is_file() and f.name!='MANIFEST.json'}
 need(actual==set(manifest),('strict packet leaf paths',sorted(actual^set(manifest))))
 for rel,digest in manifest.items():need(hashlib.sha256((P/rel).read_bytes()).hexdigest()==digest,('hash',rel))
 inputs=json.loads((P/'INPUTS.json').read_text());r=json.loads((P/'RESULTS.json').read_text())
 need(r['status']=='PASS' and len(r['analytic_controls'])==7 and len(r['numerical_gates'])==28,'status/counts')
 need(all(x['passed'] for x in r['analytic_controls']+r['numerical_gates']),'recorded failed gate')
 ref=json.loads((P/'RESTORE_REFERENCE.json').read_text())
 first=r['analytical'];late1=next(q for q in r['conditional_late_runs'] if q['b']==3 and q['c']==1);late2=next(q for q in r['conditional_late_runs'] if q['b']==3 and q['c']==2.5)
 restored=dict(rho_star_c1_GeV4=first['scales'][0]['rho_initial_GeV4'],portal_width_GeV=first['original_portal']['width_GeV'],minimum_efolds_b3_c1=first['derivative_completion_path_bounds'][1]['minimum_efolds_all_positive_energy'],original_final_h_over_v=r['full_action_runs'][0]['endpoint']['h_over_v'],completed_final_h_over_v=r['full_action_runs'][2]['endpoint']['h_over_v'],late_H_c1_b3=late1['H_today_kms_Mpc'],late_w_c1_b3=late1['w_today'],late_H_c2p5_b3=late2['H_today_kms_Mpc'],late_w_c2p5_b3=late2['w_today'])
 for k,v in restored.items():need(math.isclose(v,ref[k],rel_tol=1e-7,abs_tol=0),('restored headline mismatch',k,v,ref[k]))
 print('PASS restored headline continuity:',len(restored),'recorded values')
 raw=json.loads((P/'RAW_MANIFEST.json').read_text());trajectories={}
 expected={q['kind']+('_refined' if q['refined'] else '')+'.csv' for q in r['full_action_runs']}
 expected|={'late_b%g_c%g%s.csv'%(q['b'],q['c'],'_refined' if q['refined'] else '') for q in r['conditional_late_runs']+r['conditional_late_refinements']}
 need(len(expected)==16 and set(raw)==expected,'raw manifest members')
 with tarfile.open(P/'RAW_DATA.tar.gz','r:gz') as a:
  need(set(a.getnames())==expected and len(a.getmembers())==16,'archive members')
  for name in sorted(expected):
   b=a.extractfile(name).read()
   need(len(b)==raw[name]['bytes'] and hashlib.sha256(b).hexdigest()==raw[name]['sha256'],('raw hash',name))
   trajectories[name]=np.genfromtxt(io.BytesIO(b),delimiter=',',names=True)
 q=inputs['full'];M,m,v,lam,b=[q[k] for k in ['M','m','v','lam','b']];A=math.sqrt(2/3);V0=lam*v**4/4
 phi0=M/A*math.log1p(V0/(.75*m*m*M*M));residuals=[]
 for x in r['full_action_runs']:
  name=x['kind']+('_refined' if x['refined'] else '')+'.csv';d=trajectories[name]
  need(len(d)==inputs['numerics']['samples'],('rows',name))
  h=v*d['h_over_v'];phi=phi0*np.exp(d['log_phi_ratio']);up=d['u_per_N'];hp=d['h_per_N']
  extra=b*M/v*h*up if x['kind']=='original' else b*M/v*(h*up+(d['log_phi_ratio']+math.log(phi0/M))*hp)
  kinetic=(phi*up)**2+hp**2+extra**2
  U=.75*m*m*M*M*np.expm1(-A*phi/M)**2+np.exp(-2*q['c']*A*phi/M)*lam*(h*h-v*v)**2/4
  H2=np.exp(2*d['log_H']);rs=.5*H2*kinetic+U;rr=q['rad_to_Higgs_V0']*V0*np.exp(-4*d['N'])
  fc=float(np.max(abs(3*M*M*H2-rs-rr)/(rs+rr)));wc=float(np.max(abs(rs+d['expansion_work']-rs[0]))/rs[0])
  need(fc<1e-5 and wc<1e-6,('raw conservation',name,fc,wc))
  need(np.max(abs(rs-d['rho_scalar']))/rs[0]<1e-8,('raw scalar',name))
  need(np.max(abs(kinetic-d['field_speed_per_N']**2))<1e-8,('raw metric',name))
  need(abs(h[-1]/v-x['endpoint']['h_over_v'])<1e-12,('full endpoint',name))
  if x['kind']=='completed':need(abs(d['z'][-1]-d['z'][0])<=np.trapezoid(np.sqrt(kinetic),d['N'])+1e-8,('path bound',name))
  residuals.append(dict(name=name,Friedmann=fc,scalar_work=wc))
 p=inputs['physical'];H=p['H0_reference_kms_Mpc']/p['Mpc_km']*p['hbar_GeV_s'];C=3*p['M_GeV']**2*H**2
 for x in r['conditional_late_runs']+r['conditional_late_refinements']:
  name='late_b%g_c%g%s.csv'%(x['b'],x['c'],'_refined' if x['refined'] else '');d=trajectories[name]
  need(len(d)==inputs['late']['samples'],('late rows',name))
  r0=(x['c']*p['lam']*p['v_GeV']**4/4)**2/(p['alpha']**2*p['M_GeV']**4)
  V=r0/C*np.exp(2*d['canonical_displacement_over_M']/x['b']);E2=(d['H_kms_Mpc']/p['H0_reference_kms_Mpc'])**2;K=E2*d['velocity_per_N']**2/6
  wc=float(np.max(abs(K+V+d['expansion_work_over_reference_critical']-r0/C))/(r0/C))
  need(wc<1e-6,('late raw work',name,wc));need(abs((K[-1]-V[-1])/(K[-1]+V[-1])-x['w_today'])<1e-12,('late w',name))
  need(abs(d['H_kms_Mpc'][-1]-x['H_today_kms_Mpc'])<1e-10,('late H',name))
 print('PASS exact packet manifest:',len(manifest),'files; archive:16 raw CSV; raw energy/metric reconstruction')
 print(json.dumps(residuals,indent=2))
 if args.rerun:
  with tempfile.TemporaryDirectory(prefix='grut-freezing-') as tmp:
   t=Path(tmp);(t/'code').mkdir();shutil.copy2(P/'code/build.py',t/'code/build.py')
   env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
   proc=subprocess.run([sys.executable,str(t/'code/build.py')],capture_output=True,text=True,env=env)
   need(proc.returncode==0,('rerun failed',proc.stdout[-2000:],proc.stderr[-2000:]))
   new=json.loads((t/'RESULTS.json').read_text());need(set(leaves(new))==set(leaves(r)),'result leaf schema')
   need(new['status']=='PASS' and all(x['passed'] for x in new['analytic_controls']+new['numerical_gates']),'new failed gate')
   oldp=leaves(projection(r));newp=leaves(projection(new));need(set(oldp)==set(newp),'scientific leaf paths')
   for path,val in oldp.items():
    other=newp[path]
    if isinstance(val,(float,int)) and not isinstance(val,bool):need(math.isclose(val,other,rel_tol=2e-6,abs_tol=0),('scientific mismatch',path,val,other))
    else:need(val==other,('scientific mismatch',path,val,other))
   need(json.loads((t/'INPUTS.json').read_text())==inputs,'changed inputs')
   print('PASS fresh temporary rerun:',len(oldp),'scientific leaves; strict complete leaf schema; all 35 checks')
 print('PASS internal verification; no independent scientific replication')
if __name__=='__main__':main()
