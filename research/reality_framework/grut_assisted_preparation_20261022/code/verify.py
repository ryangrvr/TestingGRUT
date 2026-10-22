"""Read-only integrity, raw moments and optional fresh same-author calculation."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile
import numpy as np
P=Path(__file__).resolve().parents[1]
def need(ok,message):
 if not ok:raise AssertionError(message)
def leaves(x,path=()):
 if isinstance(x,dict):return {k:v for a,b in x.items() for k,v in leaves(b,path+(a,)).items()}
 if isinstance(x,list):return {k:v for a,b in enumerate(x) for k,v in leaves(b,path+(a,)).items()}
 return {path:x}
def projection(r):
 return dict(calibration=r['calibration'],snapshot=r['snapshot'],probe=r['probe'],probe_moments=r['probe_moments'],linear_finite_band=r['linear_finite_band'],fp=[dict(dx=q['dx'],L=q['L'],rows=[{k:v for k,v in x.items() if k not in ['norm','minimum','edge_mass']} for x in q['rows']]) for q in r['fp_runs']],modes=[dict(c=q['c'],xinitial=q['xinitial'],powers=[{k:v for k,v in x.items() if k.startswith('power_')} for x in q['rows']]) for q in r['mode_controls']])
def main():
 pa=argparse.ArgumentParser();pa.add_argument('--rerun',action='store_true');args=pa.parse_args()
 m=json.loads((P/'MANIFEST.json').read_text())['files']
 actual={str(f.relative_to(P)) for f in P.rglob('*') if f.is_file() and f.name!='MANIFEST.json'}
 need(actual==set(m),('payload leaf mismatch',sorted(actual^set(m))))
 for name,digest in m.items():need(hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,('hash',name))
 r=json.loads((P/'RESULTS.json').read_text());need(r['status']=='PASS' and len(r['checks'])==27 and all(q['passed'] for q in r['checks']),'recorded controls')
 for q in r['fp_runs']+[r['brownian_control']]:
  name='fp_dx%g_L%g.csv'%(q['dx'],q['L']) if q is not r['brownian_control'] else 'brownian_control.csv'
  d=np.genfromtxt(P/name,delimiter=',',names=True);x=d['canonical_over_Hdiv2pi'];pos=np.maximum(x,0);zero=np.argmin(abs(x))
  need(len(d)==q['points'],('grid length',name))
  for row in q['rows']:
   prob=d['prob_N%g'%row['N']];moments=dict(positive_probability=float(sum(prob[x>0])+.5*prob[zero]),positive_first_moment=float(np.dot(pos,prob)),positive_second_moment=float(np.dot(pos*pos,prob)))
   need(abs(sum(prob)-1)<1e-8 and min(prob)>=-1e-12,('raw probability',name))
   need(sum(prob[:5])+sum(prob[-5:])<1e-8,('raw boundaries',name))
   for k,v in moments.items():need(math.isclose(v,row[k],rel_tol=1e-12,abs_tol=1e-14),('raw moment',name,k,v,row[k]))
 print('PASS exact manifest:',len(m),'files; five raw probability grids, normalization and moments; 27 controls')
 if args.rerun:
  with tempfile.TemporaryDirectory(prefix='grut-assisted-') as tmp:
   t=Path(tmp);(t/'code').mkdir();shutil.copy2(P/'code/calculate.py',t/'code/calculate.py')
   env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
   p=subprocess.run([sys.executable,str(t/'code/calculate.py')],capture_output=True,text=True,env=env)
   need(p.returncode==0,('fresh calculation failed',p.stdout[-2000:],p.stderr[-2000:]))
   n=json.loads((t/'RESULTS.json').read_text());need(set(leaves(n))==set(leaves(r)),'complete leaf schema')
   need(n['status']=='PASS' and all(q['passed'] for q in n['checks']),'fresh controls')
   aa=leaves(projection(r));bb=leaves(projection(n));need(set(aa)==set(bb),'scientific leaf schema')
   for k,v in aa.items():
    if isinstance(v,(int,float)) and not isinstance(v,bool):need(math.isclose(v,bb[k],rel_tol=2e-6,abs_tol=0),('science comparison',k,v,bb[k]))
    else:need(v==bb[k],('value comparison',k,v,bb[k]))
   need(json.loads((P/'INPUTS.json').read_text())==json.loads((t/'INPUTS.json').read_text()),'changed inputs')
   print('PASS fresh calculation, all 27 checks and',len(aa),'scientific leaf comparisons; complete leaf schema')
 print('PASS internal verification only; no independent scientific replication')
if __name__=='__main__':main()
