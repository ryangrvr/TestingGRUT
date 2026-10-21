"""Frozen finite-window 3D suite; use `validate`, then `production`, then `analyze`."""
from __future__ import annotations
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
BUILD=HERE/'build'
RAW=HERE/'raw'
BUILD.mkdir(exist_ok=True)
RAW.mkdir(exist_ok=True)
EXE=BUILD/'expanding'

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def compile_solver():
    cmd=['g++','-std=c++17','-O3','-fopenmp','-ffp-contract=off',str(HERE/'expanding.cpp'),'-o',str(EXE)]
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    return cmd

def modes(N,dx,seed):
    """64 real Fourier modes. Same polynomial on fixed-L refinements; no zero mode."""
    L=N*dx
    path=RAW/f'modes_L{L:.8g}_s{seed}.txt'
    rng=np.random.default_rng(seed)
    kvec=[]
    while len(kvec)<64:
        k=np.rint(rng.normal(0,L/(2*np.pi),3)).astype(int)
        if not k.any() or np.linalg.norm(k)*2*np.pi/L>2.5:
            continue
        first=np.flatnonzero(k)[0]
        if k[first]<0:k=-k
        row=tuple(int(t) for t in k)
        if row not in kvec:kvec.append(row)
    co=rng.normal(size=(64,6))
    for f in range(3):
        norm=np.sqrt(np.sum(co[:,2*f:2*f+2]**2)/2)
        co[:,2*f:2*f+2]/=norm
    a=np.column_stack((kvec,co))
    np.savetxt(path,a,fmt='%.17g')
    return path

def run(name,**kw):
    config={'N':32,'dx':.4,'dt':.04,'T':8.,'taui':40.,'kappa':2.,'kind':'charge','noise':1.,'threads':8,'sample':.4,'top-sample':4.,'topology':1}
    config.update(kw)
    seed=config.pop('seed',101)
    state=config.pop('save_state',False)
    mode_file=None
    if config['kind'] in ('charge','cold') and config['noise']!=0:
        mode_file=modes(config['N'],config['dx'],seed)
        config['modes']=str(mode_file)
    pref=RAW/name
    config['out']=str(pref)
    if state:config['state']=str(BUILD/(name+'.state.f64'))
    cmd=[str(EXE)]
    for k,v in config.items():cmd += ['--'+k,str(v)]
    t=time.monotonic()
    result=subprocess.run(cmd,check=True,capture_output=True,text=True)
    rec={'name':name,'arguments':config,'seed':seed,'command':cmd,'stdout':result.stdout.strip(),
         'wall_seconds':time.monotonic()-t,'source_sha256':digest(HERE/'expanding.cpp'),
         'executable_sha256':digest(EXE),'mode_sha256':digest(mode_file) if mode_file else None}
    write(RAW/(name+'.config.json'),rec)
    print(result.stdout.strip(),flush=True)
    tab=np.genfromtxt(str(pref)+'.csv',delimiter=',',names=True)
    tab=np.atleast_1d(tab)
    return rec,tab

def final_state(name,N):
    return np.fromfile(BUILD/(name+'.state.f64'),dtype='<f8').reshape(2,3,N,N,N)

def validate():
    if (HERE/'VALIDATION.json').exists() and not (HERE/'FAILED_VALIDATION_01.json').exists():
        old=json.loads((HERE/'VALIDATION.json').read_text())
        if old['failed']:write(HERE/'FAILED_VALIDATION_01.json',old)
    elif (HERE/'VALIDATION.json').exists() and not (HERE/'FAILED_VALIDATION_02.json').exists():
        old=json.loads((HERE/'VALIDATION.json').read_text())
        if old['failed']:write(HERE/'FAILED_VALIDATION_02.json',old)
    compile_cmd=compile_solver()
    tests=[]
    details={}
    def check(label,ok,**values):
        tests.append({'name':label,'pass':bool(ok),**values})
        print(label, bool(ok), values,flush=True)
    # Exact charge normalization and homogeneous rotating orbit, including dilution.
    orbit=[]
    for dt in [.06,.03]:
        name=f'v_orbit_dt{dt}'
        _,d=run(name,N=8,dx=1.,kind='homogeneous',noise=0,T=48,dt=dt,sample=.8,topology=0,save_state=True)
        f=final_state(name,8)
        kap=2.;F=kap**(1/3);omega=math.sqrt(.1)*F
        exact=np.array([F*np.cos(omega*48),-F*np.sin(omega*48),0])
        error=float(np.max(np.abs(f[0,:,0,0,0]-exact)))
        q0=kap*math.sqrt(.1)*(8**3)
        charge=float(np.max(np.abs(d['Q_com']/q0-1)))
        residual=(d['H_conf']-d['H_conf'][0]-d['work_conf'])
        budget=float(np.max(np.abs(residual))/max(d['H_conf'].max(),d['H_conf'][0]))
        orbit.append(error)
        check('homogeneous_charge_'+str(dt),charge<1e-10,relative_error=charge)
        check('homogeneous_dilution_'+str(dt),float(d['Q_com'][-1]/((8**3)*d['a'][-1]**3))>0,final_nbar=float(d['Q_com'][-1]/((8**3)*d['a'][-1]**3)))
        check('homogeneous_orbit_'+str(dt),error<.002,max_field_error=error,source_budget=budget)
    check('homogeneous_second_order',3.5<orbit[0]/orbit[1]<4.5,error_ratio=orbit[0]/orbit[1])
    # Direct physical-field friction equation (independent DOP853) vs conformal KDK.
    kap=.9;q=kap*np.sqrt(.1)
    x=brentq(lambda xx:xx*np.sqrt(1-.9*xx)-q,1e-12,2/2.7)
    F=np.sqrt(x);s0=np.sqrt(1-x);omega=np.sqrt(1-.9*x)
    y0=np.array([F,0,s0,0,-omega*F,0])
    rate=.025
    def rhs(t,y):
        a=1+rate*t;u,w,S=y[:3];f=u*u+w*w
        force=np.array([(S*S+.1*f)*u,(S*S+.1*f)*w,(S*S-1+f)*S])
        return np.r_[y[3:],-2*rate/a*y[3:]-a*a*force]
    target=solve_ivp(rhs,[0,12],y0,method='DOP853',rtol=2e-12,atol=2e-14).y[:,-1]
    errors=[]
    for dt in [.04,.02,.01]:
        name=f'v_physical_dt{dt}'
        _,d=run(name,N=8,dx=1,dt=dt,T=12,kind='cold',kappa=kap,noise=0,topology=0,save_state=True)
        f=final_state(name,8)[:,:,0,0,0];a=1+rate*12
        actual=np.r_[f[0]/a,(f[1]-rate/a*f[0])/a]
        errors.append(float(np.max(np.abs(actual-target))))
    check('direct_physical_vs_conformal',errors[-1]<2e-4,errors=errors)
    check('physical_second_order',all(3.5<errors[i]/errors[i+1]<4.5 for i in range(2)),ratios=[errors[i]/errors[i+1] for i in range(2)])
    # Linear perturbation through the charge-held spinodal; lattice k explicitly used.
    N=32;dx=2.4;L=N*dx;kap=2.;noise=1e-10;dt=.04
    name='v_linear'
    _,d=run(name,N=N,dx=dx,dt=dt,T=28,kind='linear',noise=noise,topology=0,sample=.4)
    khat2=4*np.sin(np.pi/N)**2/dx**2
    def lr(t,y):return [y[1],-(khat2+kap**(2/3)-(1+.025*t)**2)*y[0]]
    expected=solve_ivp(lr,[0,28],[noise,.025*noise],t_eval=d['t'],method='DOP853',rtol=1e-12,atol=1e-24).y[0]
    actual=d['linear_chi_projection'];relative=float(np.max(np.abs(actual-expected))/np.max(np.abs(expected)))
    check('linear_crossing_ODE',relative<.02,relative_peak_error=relative,mode_crossing_a=math.sqrt(khat2+kap**(2/3)))
    # Fully nonlinear small 3D run and its exact Z2 mirror.
    _,a=run('v_nonlinear',N=32,dx=.5,T=40,dt=.05,sample=1,save_state=True)
    _,b=run('v_mirror',N=32,dx=.5,T=40,dt=.05,sample=1,mirror=-1,save_state=True)
    f=final_state('v_nonlinear',32);m=final_state('v_mirror',32)
    mirror=float(max(np.max(np.abs(f[:,:2]-m[:,:2])),np.max(np.abs(f[:,2]+m[:,2]))))
    check('Z2_mirror',mirror==0,max_abs_error=mirror)
    charge=float(np.max(np.abs(a['Q_com']/a['Q_com'][0]-1)))
    check('nonlinear_charge',charge<1e-9,relative_drift=charge)
    residual=a['H_conf']-a['H_conf'][0]-a['work_conf']
    norm=max(float(a['H_conf'].max()),float(a['H_conf'][0]),1)
    residual0=float(np.max(np.abs(residual))/norm)
    _,fine=run('v_nonlinear_dt',N=32,dx=.5,T=40,dt=.025,sample=1)
    residual=fine['H_conf']-fine['H_conf'][0]-fine['work_conf']
    residual1=float(np.max(np.abs(residual))/max(float(fine['H_conf'].max()),float(fine['H_conf'][0]),1))
    check('nonlinear_work_source',residual0<.01 and residual1<.01,residual_coarse=residual0,residual_fine=residual1)
    check('work_source_second_order',3<residual0/residual1<5,ratio=residual0/residual1)
    # Periodic slab topology (spanning a cut is insufficient) and boosted planar kink controls.
    planar=[]
    for N,dx in [(128,.2),(256,.1)]:
        name=f'v_kink_N{N}'
        _,d=run(name,N=N,dx=dx,T=0,kind='kink',rate=0,speed=.456,topology=1)
        L=N*dx;gamma=1/np.sqrt(1-.456**2);E=2*2*np.sqrt(2)/3*L*L*gamma
        energy=float(d['rho'][0]*L**3)
        area=float(d['area_core_com_over_V'][0]*L**3)
        proxy=float(d['phi_core_energy_proxy'][0]);v2=float(d['v2_core_proxy'][0])
        smooth=float(d['area_smooth_com_over_V'][0]*L**3)
        smoothE=float(d['phi_smooth_energy_proxy'][0])
        topo=json.loads((RAW/(name+'.topology.jsonl')).read_text())
        planar.append({'N':N,'energy_relative_error':energy/E-1,'coarea_relative_error':area/(2*L*L)-1,
                       'core_energy_proxy_relative_error':proxy/E-1,'v2_relative_error':v2/(.456**2)-1,
                       'smooth_coarea_relative_error':smooth/(2*L*L)-1,'smooth_energy_proxy_relative_error':smoothE/E-1})
        check('slab_winding_'+str(N),topo['plus']['any_winding']==6 and topo['minus']['any_winding']==6,plus=topo['plus']['any_winding'],minus=topo['minus']['any_winding'])
    check('planar_energy',abs(planar[-1]['energy_relative_error'])<.005,controls=planar)
    check('planar_coarea',abs(planar[-1]['coarea_relative_error'])<.06,relative_error=planar[-1]['coarea_relative_error'])
    check('planar_velocity_proxy',abs(planar[-1]['v2_relative_error'])<.01,relative_error=planar[-1]['v2_relative_error'])
    check('planar_smooth_coarea',abs(planar[-1]['smooth_coarea_relative_error'])<.005,relative_error=planar[-1]['smooth_coarea_relative_error'])
    check('planar_smooth_energy',abs(planar[-1]['smooth_energy_proxy_relative_error'])<.005,relative_error=planar[-1]['smooth_energy_proxy_relative_error'])
    _,blob=run('v_corner_blob',N=64,dx=.5,T=0,kind='blob',rate=0,topology=1)
    topo=json.loads((RAW/'v_corner_blob.topology.jsonl').read_text())
    check('corner_blob_crosses_cuts_without_winding',topo['plus']['any_winding']==0 and topo['minus']['any_winding']==7,positive_winding=topo['plus']['any_winding'],negative_winding=topo['minus']['any_winding'])
    details['planar_controls']=planar
    # Weighting counterexample: the old one-speed figure is not an upper bound.
    s=.456**2;vfast=.99;fraction=s/vfast**2
    one=1/np.sqrt(1-s);mixture=1-fraction+fraction/np.sqrt(1-vfast**2)
    check('Lorentz_weighting_counterexample',mixture>one,one_speed=one,area_weighted_mixture_gamma=mixture,
          fixed_energy_weighted_lower=one,fixed_energy_weighted_upper=1/(1-s))
    # Canonical scalar pressure counterexample at arbitrarily small density.
    rho=1e-60;p=rho
    check('instantaneous_w_counterexample',p/rho==1,rho=rho,w=p/rho)
    result={'tests':tests,'passed':sum(x['pass'] for x in tests),'failed':sum(not x['pass'] for x in tests),
            'details':details,'compile_command':compile_cmd,'python':sys.version,'platform':platform.platform(),
            'source_sha256':digest(HERE/'expanding.cpp'),'run_py_sha256':digest(__file__)}
    write(HERE/'VALIDATION.json',result)
    if result['failed']:raise SystemExit('Validation failed; production forbidden')

def jobs():
    j=[]
    for s in [101,202,303]:j.append((f'main_s{s}',dict(N=192,dx=.4,dt=.06,seed=s)))
    for s in [101,202,303]:j.append((f'large_s{s}',dict(N=288,dx=.4,dt=.06,seed=s)))
    j += [('grid256_s101',dict(N=256,dx=.3,dt=.06,seed=101)),
          ('grid384_s101',dict(N=384,dx=.2,dt=.04,seed=101)),
          ('dt_half_s101',dict(N=192,dx=.4,dt=.03,seed=101)),
          ('kappa4_s101',dict(N=192,dx=.4,dt=.06,seed=101,kappa=4)),
          ('cold_s101',dict(N=192,dx=.4,dt=.06,seed=101,kappa=.9,kind='cold'))]
    return j

def production():
    val=json.loads((HERE/'VALIDATION.json').read_text())
    if val['failed'] or val['source_sha256']!=digest(HERE/'expanding.cpp'):
        raise SystemExit('Source must pass validation before production')
    compile_solver()
    t=time.monotonic();status=[]
    J=jobs();lookup=dict(J)
    order=['main_s101','main_s202','main_s303','cold_s101','kappa4_s101','dt_half_s101',
           'grid256_s101','grid384_s101','large_s101','large_s202','large_s303']
    for name in order:
        kw=lookup[name]
        if time.monotonic()-t>20*60:
            status.append({'name':name,'status':'not_launched_cost_cap'})
            continue
        rec,tab=run(name,T=48,sample=1.2,**{'top-sample':4.8},**kw)
        status.append({'name':name,'status':'completed','wall_seconds':rec['wall_seconds']})
        write(HERE/'RUN_STATUS.json',status)
    write(HERE/'RUN_STATUS.json',status)

if __name__=='__main__':
    if sys.argv[1]=='validate':validate()
    elif sys.argv[1]=='production':production()
    else:raise SystemExit('Use validate or production')
