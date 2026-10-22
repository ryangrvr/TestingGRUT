"""Read-only artifact and raw-physics checks; optional fresh re-execution."""
from pathlib import Path
import argparse, hashlib, json, math, os, shutil, subprocess, sys, tempfile
import numpy as np
from scipy.integrate import quad, simpson
from scipy.interpolate import CubicSpline

P=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def leaves(x,p=''):
    if isinstance(x,dict):
        for k,v in x.items():
            if k=='nfev':continue
            yield from leaves(v,p+'/'+k)
    elif isinstance(x,list):
        for j,v in enumerate(x):yield from leaves(v,p+'/'+str(j))
    else:yield p,x


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    errors=[];raw_checks=0
    m=json.loads((P/'MANIFEST.json').read_text())
    for f,h in m['files_sha256'].items():
        if not (P/f).is_file() or sha(P/f)!=h:errors.append('hash '+f)
    actual={str(f.relative_to(P)) for f in P.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name!='MANIFEST.json'}
    if actual!=set(m['files_sha256']):errors.append('manifest schema')
    d=json.loads((P/'RESULTS.json').read_text());I=json.loads((P/'INPUTS.json').read_text())
    if d['status']!='PASS' or len(d['checks'])!=42 or not all(x['passed'] for x in d['checks']):errors.append('internal dispositions')
    for x in d['calibration']:
        amp=x['pole_fraction_F']*x['inflaton_mass_GeV']**2/(2*math.sqrt(x['kappa'])*I['M_GeV'])
        rate=amp**2/(32*math.pi*x['inflaton_mass_GeV'])
        if abs(rate/x['partial_width_GeV']-1)>2e-14:errors.append('phase space width')
        raw_checks+=1
    for x in [*d['pulse'],d['pulse_refined']]:
        a=np.loadtxt(P/x['file'],delimiter=',',skiprows=1)
        if a.shape!=(x['npoints'],5):errors.append('pulse shape')
        if np.max(abs(a[:,1]-(a[:,3]**2+a[:,4]**2)))>1e-15:errors.append('mode occupancy')
        energy=float(simpson(a[:,0]**3*a[:,1],x=a[:,0])/(2*math.pi**2))
        if abs(energy/x['rho_over_m4']-1)>1e-12:errors.append('mode energy integral')
        raw_checks+=3
    for x in [*d['receivers'],d['receiver_refined']]:
        a=np.loadtxt(P/x['file'],delimiter=',',skiprows=1)
        t,N,h,u,rd,rv=a.T
        c=float(np.max(abs(h*h-u-rd-rv)/(h*h+abs(u)+abs(rd)+abs(rv))))
        if abs(c-x['constraint_relative'])>1e-12 or c>2e-7:errors.append('raw receiver constraint')
        # Independent retarded energy integral using rho_I a^3 exp(t)=constant.
        fit=CubicSpline(np.log(t),N)
        integral=quad(lambda tt:math.exp(float(fit(math.log(tt)))-tt),t[0],t[-1],epsabs=1e-7,epsrel=1e-10)[0]
        predicted=u[0]*math.exp(t[0]-4*N[-1])*integral
        err=abs(predicted/(rd[-1]+rv[-1])-1)
        if err>2e-6:errors.append('retarded receiver energy '+str(err))
        raw_checks+=2
    replay_leaves=0;replay_arrays=0
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='grut-inflaton-replay-') as td:
            q=Path(td);(q/'code').mkdir();(q/'data').mkdir()
            shutil.copy2(P/'INPUTS.json',q/'INPUTS.json');shutil.copy2(P/'code/calculate.py',q/'code/calculate.py')
            env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
            r=subprocess.run([sys.executable,str(q/'code/calculate.py')],cwd=q,env=env,text=True,capture_output=True,timeout=180)
            if r.returncode:raise RuntimeError(r.stdout+r.stderr)
            new=json.loads((q/'RESULTS.json').read_text());oldl=dict(leaves(d));newl=dict(leaves(new))
            if set(oldl)!=set(newl):errors.append('complete scientific leaf schema')
            for k,v in oldl.items():
                w=newl.get(k)
                if isinstance(v,(int,float)) and not isinstance(v,bool):
                    diagnostic=('/checks/' in k or 'error' in k or 'relative_quantum_Born' in k)
                    floor=2e-13 if diagnostic or v==0 else 0.0
                    if not isinstance(w,(int,float)) or not math.isclose(v,w,rel_tol=2e-8,abs_tol=floor):errors.append('replay leaf '+k)
                    replay_leaves+=1
                elif w!=v:errors.append('replay metadata '+k)
            oldfiles={f.name for f in (P/'data').glob('*.csv')};newfiles={f.name for f in (q/'data').glob('*.csv')}
            if oldfiles!=newfiles:errors.append('complete raw array schema')
            for f in sorted(oldfiles):
                aa=np.loadtxt(P/'data'/f,delimiter=',',skiprows=1);bb=np.loadtxt(q/'data'/f,delimiter=',',skiprows=1)
                if aa.shape!=bb.shape or not np.allclose(aa,bb,rtol=2e-8,atol=2e-13):errors.append('replay raw '+f)
                replay_arrays+=1
    print(json.dumps({'passed':not errors,'internal_checks':len(d['checks']),'raw_checks':raw_checks,
          'manifest_files':len(m['files_sha256']),'replay_scientific_leaves':replay_leaves,
          'replay_raw_arrays':replay_arrays,'errors':errors},indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
