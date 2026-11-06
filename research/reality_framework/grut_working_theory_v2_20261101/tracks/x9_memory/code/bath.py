#!/usr/bin/env python3
"""Signed lepton-Higgs vacuum spectral model and paired H55 preparation replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import ctypes
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, trapezoid
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import root
from scipy.special import roots_legendre

TRACKS=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(TRACKS/'x7_regulator/code'))
from overlap import initial_data, read_columns

LABELS=('norm','sterile_share','sterile_free_count','bath_positive','bath_negative',
        'bath_chiral_count','sterile_energy','bath_energy','interaction_energy','total_energy',
        'redshift_work','scalar_work')
SETTINGS={
    'reference': (70000.,24.,1.2e-5,'half'),
    'frequency': (70000.,12.,1.2e-5,'half'),
    'cutoff': (105000.,24.,1.2e-5,'half'),
    'timestep': (70000.,24.,6e-6,'half'),
    'background': (70000.,24.,1.2e-5,'radial'),
}
REFINED={
    'reference': (105000.,24.,1.2e-5,'half',2.),
    'frequency': (105000.,12.,1.2e-5,'half',1.),
    'timestep': (105000.,24.,6e-6,'half',2.),
    'cutoff': (140000.,24.,1.2e-5,'half',2.),
}
DENSE={**REFINED, 'timestep_fine': (105000.,24.,3e-6,'half',2.)}
BEGIN,END=0.16,0.42


def taper(x):
    u=np.clip(2*x-1,0,1)
    return 1-u*u*(3-2*u)


def principal(y):
    if np.any(np.real(y)<=0) or np.any(np.real(y)>=0.5):
        raise ValueError('Mass renormalization point outside the untapered band')
    d=((16*y-36)*y+24)*y-5
    return np.log(y/(1-y))+d*np.log((0.5-y)/(1-y))-(19/6-12*y+8*y*y)


def mass_matched(m,gamma,cutoff):
    return m*(1-gamma*principal((m/cutoff)**2))


@dataclass
class Star:
    k: float
    frequency: np.ndarray
    coupling: np.ndarray

    def matrix(self,m,delta=0.):
        n=len(self.frequency)
        h=np.diag(np.r_[self.k+delta,-self.k-delta,self.frequency+delta,self.frequency-delta])
        h[0,1]=h[1,0]=m
        for b in (0,1):
            sec=slice(2+b*n,2+(b+1)*n)
            h[b,sec]=h[sec,b]=self.coupling[b]
        return h

    def resolve(self,m,z,v,delta=0.):
        n=len(self.frequency)
        u=v[2:].reshape(2,n,-1)
        shifts=np.array([delta,-delta])
        d=1/(z-self.frequency[None,:]-shifts[:,None])
        g=self.coupling
        self_energy=np.sum(g*g*d,axis=1)
        a,b=z-self.k-delta-self_energy[0],z+self.k+delta-self_energy[1]
        target=v[:2]+np.sum((g*d)[:,:,None]*u,axis=1)
        s=np.array([b*target[0]+m*target[1],m*target[0]+a*target[1]])/(a*b-m*m)
        bath=d[:,:,None]*(u+g[:,:,None]*s[:,None,:])
        return np.vstack([s,bath.reshape(2*n,-1)])

    def positive(self,m,v,nodes=64):
        x,w=roots_legendre(nodes)
        angle=(x+1)*np.pi/4
        scale=np.hypot(self.k,m)
        values=scale*np.tan(angle)
        weight=scale*w/(4*np.cos(angle)**2)
        result=v/2
        for tau,d in zip(values,weight,strict=True):
            result-=d*self.resolve(m,1j*tau,v).real
        return result


def build(k,gamma,cutoff,spacing,inner=None):
    q=abs(k)
    edge=np.hypot(k,cutoff)
    panels=[(q,edge,spacing)] if inner is None else [(q,q+1024.,inner),(q+1024.,edge,spacing)]
    points=[]
    width=[]
    for lo,hi,goal in panels:
        n=int(np.ceil((hi-lo)/goal))
        d=(hi-lo)/n
        points.extend(lo+d*(np.arange(n)+0.5))
        width.extend(np.full(n,d))
    positive=np.array(points)
    widths=np.array(width)
    frequency=np.r_[positive,-positive]
    invariant=(frequency**2-k*k)/cutoff**2
    weight=gamma*taper(invariant)*np.tile(widths,2)
    coupling=np.sqrt(weight[None,:]*abs(frequency[None,:]-np.array([1,-1])[:,None]*k))
    return Star(k,frequency,coupling),float(max(widths))


def prepare(star,m,physical_mass,width=128.):
    E=np.hypot(star.k,physical_mass)
    n=len(star.frequency)
    source=np.zeros((2+2*n,3))
    source[:2,:2]=np.eye(2)
    packet=np.exp(-0.5*((star.frequency-E)/width)**2)*(star.frequency>0)
    source[2:,2]=(star.coupling*packet).ravel()
    projected=star.positive(m,source)
    u,s,_=np.linalg.svd(projected[:,:2],full_matrices=False)
    if s[-1]<1e-10:
        raise ValueError('Cannot isolate both projected sterile coordinates')
    spin=np.sqrt((1+np.array([1.,-1.])*star.k/E)/2)
    sterile=projected[:,:2] @ spin
    sterile/=np.linalg.norm(sterile)
    bath=projected[:,2]-u @ (u.T @ projected[:,2])
    bath/=np.linalg.norm(bath)
    check=star.positive(m,np.column_stack([sterile,bath]),nodes=96)
    error=float(np.max(abs(check-np.column_stack([sterile,bath]))))
    pair=np.array([sterile+1j*bath,sterile-1j*bath])/np.sqrt(2)
    return np.ascontiguousarray(pair),dict(positive_projection_error=error,
                                          orthogonality=float(abs(sterile @ bath)),
                                          bath_packet_sterile_leakage=float(np.linalg.norm(bath[:2])),
                                          sterile_covariance_difference=float(np.max(abs(
                                              np.outer(pair[0,:2],pair[0,:2].conj())-
                                              np.outer(pair[1,:2],pair[1,:2].conj())))))


def engine(work):
    binary=work/'star.so'
    subprocess.run(['gcc','-std=c11','-O3','-march=native','-fcx-limited-range','-Wall','-Wextra','-Werror',
                    '-shared','-fPIC',str(Path(__file__).with_name('star.c')),'-lm','-o',str(binary)],check=True)
    lib=ctypes.CDLL(str(binary))
    lib.propagate.argtypes=[ctypes.c_int]*3+[ctypes.c_double]*2+[ctypes.c_void_p]*5
    lib.propagate.restype=ctypes.c_int
    lib.propagate_pade.argtypes=lib.propagate.argtypes+[ctypes.c_void_p]
    lib.propagate_pade.restype=ctypes.c_int
    return lib


def evolve(lib,star,path,eta,pair,stride,pade=False):
    steps=len(eta)-1
    keep=np.r_[0,np.arange(stride,steps,stride),steps]
    data=np.empty((len(keep),2,len(LABELS)))
    state=pair.copy()
    args=(len(star.frequency),steps,stride,float(eta[1]-eta[0]),star.k,
          star.frequency.ctypes.data,star.coupling.ctypes.data,path.ctypes.data,
          state.ctypes.data,data.ctypes.data)
    checks=np.zeros(2)
    n=lib.propagate_pade(*args,checks.ctypes.data) if pade else lib.propagate(*args)
    if n!=len(keep) or not np.all(np.isfinite(data)):
        raise RuntimeError('Invalid propagated state or ledger')
    return keep,dict(zip(LABELS,data.transpose(2,0,1),strict=True)),state,checks


def relevant(history):
    t=history['t']
    keep=(t>=BEGIN-0.003)&(t<=END+0.002)
    names=('t','a','H','sigma','dsigma','chi','dchi')
    return {key:history[key][keep] for key in names}


def background(history,parents,scales,step):
    t=history['t']
    if len(t)<100 or not (t[0]<BEGIN and t[-1]>END):
        raise ValueError('Need the high-cadence production history through 0.42')
    a=history['a']
    scale=CubicHermiteSpline(t,a,a*history['H'])
    conformal=PchipInterpolator(t,cumulative_trapezoid(1/a,t,initial=0))
    low,high=conformal([BEGIN,END])
    steps=int(np.ceil((high-low)/step))
    eta=np.linspace(low,high,steps+1)
    time=PchipInterpolator(conformal(t),t)(eta)
    velocity=[]
    for field,vel in (('sigma','dsigma'),('chi','dchi')):
        spline=CubicHermiteSpline(t,history[field],history[vel])
        velocity.append(spline(time,1))
    a=scale(time)
    shift=np.array(velocity)[parents].T*a[:,None]/(2*scales)
    return time,eta-low,a,shift


def select(model,terminal,state,widths):
    q=state['q']
    a=terminal['a']
    out=[]
    for index in (0,1,3):
        parent,mass,scale,_=model['species'][index]
        k=np.array([-1,1])[:,None]*q/a
        E=np.hypot(k,mass)
        r=state['bloch'][index]
        occupation=(1+(mass*r[...,0]+k*r[...,2])/E)/2
        cum=np.cumsum(state['w']*np.sum(E*occupation,axis=0))
        i=int(np.searchsorted(cum,0.5*cum[-1]))
        for sign in (-1,1):
            out.append(dict(species=index+1,parent=parent,mass=mass,phase_scale=scale,
                            width=widths[index],q=float(q[i]),helicity=sign))
    return out


def window_average(time,values):
    low,high=END-0.02,END
    t=np.r_[low,time[(time>low)&(time<high)],high]
    return float(trapezoid(np.interp(t,time,values),t)/(high-low))


def summarize(time,data,prep,scale,dense):
    diff={k:v[:,0]-v[:,1] for k,v in data.items()}
    e=data['total_energy']
    budget=e-e[0]-data['redshift_work']-data['scalar_work']
    select=time>time[-1]-0.02
    interval=time[select][-1]-time[select][0]
    maxima={k:float(np.max(abs(y))) for k,y in diff.items()}
    return dict(
        initial_differences={k:float(y[0]) for k,y in diff.items()},
        terminal_differences={k:float(y[-1]) for k,y in diff.items()},
        final_window_differences={k:(window_average(time,y) if dense else
                                       float(trapezoid(y[select],time[select])/interval))
                                  for k,y in diff.items()},
        max_absolute_differences=maxima,
        mean_initial_energy=float(np.mean(e[0])),
        max_energy_work_error=float(np.max(abs(budget))/scale),
        max_norm_error=float(np.max(abs(data['norm']-1))),
        **prep,
    )


def stationary():
    rows=[]
    for M,G,q,cutoff in ((2.,0.025,0.7,16.),(2.,0.0025,0.7,16.),(250.,0.2256,3000.,70000.),
                         (400.,0.5776,15000.,70000.),(120.,0.052,8000.,105000.)):
        gamma=G/(2*np.pi*M)
        y=(M/cutoff)**2
        reference=quad(lambda x:-taper(x),0,1,weight='cauchy',wvar=y,epsabs=1e-11)[0]
        m=mass_matched(M,gamma,cutoff)
        E=np.hypot(q,M)
        born=G*M/E

        def residual(value):
            z=complex(*value)
            u=(z*z-q*q)/cutoff**2
            f=(z*z-q*q)*(1-gamma*(principal(u)-1j*np.pi))**2-m*m
            return f.real,f.imag

        pole=root(residual,[E,-born/2],tol=1e-11)
        if np.linalg.norm(residual(pole.x))>1e-6*M*M:
            raise RuntimeError('Unresolved continuum resonance')
        rows.append(dict(mass=M,rest_width=G,q=q,cutoff=cutoff,
                         independent_principal_error=float(abs(principal(y)-reference)),
                         pole_energy=float(pole.x[0]),on_shell_energy=float(E),
                         pole_width=float(-2*pole.x[1]),born_width=float(born),
                         born_width_relative_error=float(abs(-2*pole.x[1]/born-1))))
    return rows


def run(work,lib,name,mode_id,modes,paths,settings,dense):
    mode=modes[mode_id]
    cutoff,spacing,dt,profile,*extra=settings[name]
    inner=extra[0] if extra else None
    time,eta,a,delta=paths[name]
    m0=mode['mass']*a
    gamma=mode['width']/(2*np.pi*mode['mass'])
    mass=mass_matched(m0,gamma,cutoff)
    speed=delta[:,mode_id]
    path=np.column_stack([mass,speed,np.gradient(mass,eta,edge_order=2),m0])
    k=mode['q']*mode['helicity']
    star,actual=build(k,gamma,cutoff,spacing,inner)
    pair,prep=prepare(star,mass[0],m0[0])
    stride=max(1,int(2.4e-5/(eta[1]-eta[0]))) if dense else max(1,round(4e-4/(eta[1]-eta[0])))
    keep,obs,_,checks=evolve(lib,star,path,eta,pair,stride,pade=inner is not None)
    norm=np.hypot(k,m0[0])
    report=summarize(time[keep],obs,prep,norm,dense)
    if dense:
        tt=time[keep]
        idx=np.r_[np.arange(0,len(tt)-1,2),len(tt)-1]
        report['output_cadence_max_conformal_spacing']=float(np.max(np.diff(eta[keep])))
        report['decimated_window_absolute_differences']={
            k:abs(window_average(tt[idx],v[idx,0]-v[idx,1])-report['final_window_differences'][k])
            for k,v in obs.items()
        }
    e=np.hypot(k+speed,m0)
    branch=np.maximum(abs(e-speed),abs(e+speed))
    inv=np.maximum(branch**2-k*k,0)/cutoff**2
    report.update(actual_frequency_spacing=actual,conformal_duration=float(eta[-1]),
                  recurrence_time=2*np.pi/actual,dt=float(eta[1]-eta[0]),
                  points_per_chirality=len(star.frequency),cutoff=cutoff,
                  maximum_drive_branch_cutoff_ratio=float(np.max(inv)),inner_spacing=inner,
                  continuous_torque_difference=(float(np.max(abs(checks))/norm) if inner else None))
    np.savez_compressed(work/f'{name}_{mode_id}.npz',time=time[keep],**obs)
    return name,mode_id,report


def compare(base,variant):
    out={}
    for key in ('sterile_free_count','bath_positive','bath_energy','bath_chiral_count','scalar_work'):
        out[key]=max(abs(a['final_window_differences'][key]-b['final_window_differences'][key])/
                     max(abs(a['final_window_differences'][key]),a['max_absolute_differences'][key]*0.1,1e-12)
                     for a,b in zip(base,variant,strict=True))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('production',type=Path)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=6)
    phase=p.add_mutually_exclusive_group()
    phase.add_argument('--refined',action='store_true')
    phase.add_argument('--dense',action='store_true')
    args=p.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    model,terminal,_,state,sha,widths,_=initial_data()
    modes=select(model,terminal,state,widths)
    history={key:args.production/f'{key}_base.history' for key in ('half','radial')}
    histories={k:relevant(read_columns(path)) for k,path in history.items()}
    parents=np.array([row['parent'] for row in modes])
    scales=np.array([row['phase_scale'] for row in modes])
    settings=DENSE if args.dense else REFINED if args.refined else SETTINGS
    paths={name:background(histories[row[3]],parents,scales,row[2]) for name,row in settings.items()}
    lib=engine(work)
    tasks=[(name,i) for name in settings for i in range(len(modes))]
    results={name:[None]*len(modes) for name in settings}

    def call(task):
        return run(work,lib,*task,modes,paths,settings,args.dense)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for name,i,report in pool.map(call,tasks):
            results[name][i]=report
            print(name,i,'bath difference',report['final_window_differences']['bath_energy'],
                  'ledger',report['max_energy_work_error'],flush=True)
    summary=dict(
        scope='quadratic representation of vacuum Born self-energy on prescribed homogeneous H55 histories',
        preparation='one positive spectral excitation with opposite phase on a bath packet, common interacting vacuum',
        interval=[BEGIN,END],checkpoint_selection_sha256=sha,
        propagator='pade33' if args.refined or args.dense else 'split',dense_output=args.dense,
        source_trace_sha256={name:hashlib.sha256(np.column_stack(list(d.values())).astype('<f8').tobytes()).hexdigest()
                             for name,d in histories.items()},
        stationary_poles=stationary(),modes=modes,cases=results,
        final_window_refinement_differences={name:compare(results['reference'],v) for name,v in results.items()
                                             if name!='reference'},
        scientific_acceptance=False,
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
