"""Declared exploratory Gaussian controls; frozen integrator is included unchanged."""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_run',HERE/'run.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
EXE=R.BUILD/'gaussian'

def compile_solver():
    cmd=['g++','-std=c++17','-O3','-fopenmp','-ffp-contract=off',str(HERE/'gaussian.cpp'),'-o',str(EXE)]
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    return cmd

def noise(N,dx,seed):
    path=R.BUILD/f'gaussian_N{N}_s{seed}.noise.f64'
    rng=np.random.default_rng(seed)
    k=2*np.pi*np.fft.fftfreq(N,dx);kr=2*np.pi*np.fft.rfftfreq(N,dx)
    k2=k[:,None,None]**2+k[None,:,None]**2+kr[None,None,:]**2
    filt=np.exp(-.5*k2)
    means=[];rms=[]
    with path.open('wb') as out:
        for field in range(3):
            f=np.fft.rfftn(rng.standard_normal((N,N,N)))
            f*=filt;f.flat[0]=0
            q=np.fft.irfftn(f,s=(N,N,N),axes=(0,1,2))
            q-=q.mean();q/=q.std()
            means.append(float(q.mean()));rms.append(float(np.sqrt(np.mean(q*q))))
            q.astype('<f8',copy=False).tofile(out)
    q=np.memmap(path,dtype='<f8',mode='r',shape=(3,N,N,N))
    covariance=[[float(np.mean(q[i]*q[j])) for j in range(3)] for i in range(3)]
    rec={'N':N,'dx':dx,'L':N*dx,'seed':seed,'generator':'NumPy PCG64 standard_normal; rfftn, exp(-k^2/2), zero k=0, unit realization rms',
         'numpy_version':np.__version__,'means':means,'rms':rms,'covariance':covariance,'sha256':R.digest(path)}
    R.write(R.RAW/f'gaussian_N{N}_s{seed}.noise.json',rec)
    return path,rec

def run(name,N,dx,seed,dt=.06,T=48,mirror=1,state_in=None):
    path,init=noise(N,dx,seed)
    prefix=R.RAW/name
    state=R.BUILD/(name+'.state.f64')
    config={'N':N,'dx':dx,'dt':dt,'T':T,'noise-file':str(path),'mirror':mirror,'out':str(prefix),'threads':8,'state':str(state)}
    if state_in:config['state-in']=str(state_in)
    cmd=[str(EXE)]
    for k,v in config.items():cmd+=['--'+k,str(v)]
    start=time.monotonic();result=subprocess.run(cmd,check=True,capture_output=True,text=True)
    rec={'name':name,'arguments':config,'seed':seed,'command':cmd,'stdout':result.stdout.strip(),'wall_seconds':time.monotonic()-start,
         'frozen_source_sha256':R.digest(HERE/'expanding.cpp'),'wrapper_sha256':R.digest(HERE/'gaussian.cpp'),'noise':init}
    R.write(R.RAW/(name+'.config.json'),rec)
    print(result.stdout.strip(),flush=True)
    # Large full states are neither needed nor retained; scalar histories and slices are retained.
    if N>=128:state.unlink()
    return rec

def validate():
    compile_solver();checks=[]
    def check(name,ok,**v):checks.append({'name':name,'pass':bool(ok),**v})
    for sign in [1,-1]:run(f'g_v_mirror_{sign}',32,.5,101,T=36,mirror=sign)
    a=R.final_state('g_v_mirror_1',32);b=R.final_state('g_v_mirror_-1',32)
    err=float(max(np.max(np.abs(a[:,:2]-b[:,:2])),np.max(np.abs(a[:,2]+b[:,2]))))
    check('gaussian_exact_mirror',err==0,max_error=err)
    # Load exactly the initial state of the frozen implementation, then compare endpoints.
    R.run('g_v_frozen_initial',N=32,dx=.5,dt=.06,T=0,save_state=True,topology=0)
    R.run('g_v_frozen_evolved',N=32,dx=.5,dt=.06,T=36,save_state=True,topology=0)
    run('g_v_loaded',32,.5,101,T=36,state_in=R.BUILD/'g_v_frozen_initial.state.f64')
    f=R.final_state('g_v_frozen_evolved',32);g=R.final_state('g_v_loaded',32)
    err=float(np.max(np.abs(f-g)))
    check('identical_state_identical_evolution',err==0,max_error=err)
    meta=json.loads((R.RAW/'gaussian_N32_s101.noise.json').read_text())
    check('noise_zero_means',max(abs(x) for x in meta['means'])<1e-14,means=meta['means'])
    check('noise_unit_rms',max(abs(x-1) for x in meta['rms'])<1e-14,rms=meta['rms'])
    for name in ['g_v_mirror_1','g_v_loaded']:
        d=np.genfromtxt(R.RAW/(name+'.csv'),delimiter=',',names=True)
        drift=float(np.max(np.abs(d['Q_com']/d['Q_com'][0]-1)))
        check(name+'_charge',drift<1e-9,relative_drift=drift)
    x={'checks':checks,'passed':sum(c['pass'] for c in checks),'failed':sum(not c['pass'] for c in checks),
       'wrapper_sha256':R.digest(HERE/'gaussian.cpp'),'driver_sha256':R.digest(__file__),
       'frozen_source_sha256':R.digest(HERE/'expanding.cpp'),'numpy':np.__version__,'platform':platform.platform()}
    R.write(HERE/'GAUSSIAN_VALIDATION.json',x)
    if x['failed']:raise SystemExit('Gaussian validation failed')
    print(json.dumps(x),flush=True)

def production():
    x=json.loads((HERE/'GAUSSIAN_VALIDATION.json').read_text())
    if x['failed'] or x['wrapper_sha256']!=R.digest(HERE/'gaussian.cpp'):raise SystemExit('Validation required')
    compile_solver();start=time.monotonic();status=[]
    for name,N,seed in [('gaussian_s101',192,101),('gaussian_s202',192,202),('gaussian_large_s101',288,101)]:
        if time.monotonic()-start>360:status.append({'name':name,'status':'not_launched_cost_cap'});continue
        rec=run(name,N,.4,seed)
        status.append({'name':name,'status':'completed','wall_seconds':rec['wall_seconds']})
        R.write(HERE/'GAUSSIAN_STATUS.json',status)
    R.write(HERE/'GAUSSIAN_STATUS.json',status)

if __name__=='__main__':
    if sys.argv[1]=='validate':validate()
    elif sys.argv[1]=='production':production()
    else:raise SystemExit('Use validate or production')
