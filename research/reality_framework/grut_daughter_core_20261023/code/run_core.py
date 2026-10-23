#!/usr/bin/env python3
"""Prepared spherical Gaussian patches; no cosmological Q distribution."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
from scipy.integrate import solve_bvp, simpson

ROOT=Path(__file__).resolve().parents[1]
CASES={
    'A':dict(amplitude=.18,width=60.,length=360.,h=.25,dt=.0625,tmax=12000.,probe=120.,absorber=270.,eta=.03),
    'A_fine':dict(amplitude=.18,width=60.,length=360.,h=.125,dt=.03125,tmax=12000.,probe=120.,absorber=270.,eta=.03),
    'A_time':dict(amplitude=.18,width=60.,length=360.,h=.25,dt=.03125,tmax=12000.,probe=120.,absorber=270.,eta=.03),
    'A_box':dict(amplitude=.18,width=60.,length=480.,h=.25,dt=.0625,tmax=12000.,probe=120.,absorber=360.,eta=.03),
    'A_far':dict(amplitude=.18,width=60.,length=480.,h=.25,dt=.0625,tmax=12000.,probe=240.,absorber=360.,eta=.03),
    'A_finer':dict(amplitude=.18,width=60.,length=360.,h=.0625,dt=.015625,tmax=12000.,probe=120.,absorber=270.,eta=.03),
    'A_remote':dict(amplitude=.18,width=60.,length=1200.,h=.25,dt=.0625,tmax=36000.,probe=600.,absorber=900.,eta=.03),
    'B':dict(amplitude=.09,width=120.,length=720.,h=.25,dt=.0625,tmax=48000.,probe=180.,absorber=540.,eta=.03),
}


def load_trace(run):
    if 'trace_parts' in run:
        return np.concatenate([np.load(ROOT/path)['trace'] for path in run['trace_parts']])
    return np.load(ROOT/run['trace_file'])['trace']


def store_trace(name,trace):
    paths=[]
    for part,start in enumerate(range(0,len(trace),8192)):
        path=f'data/{name}_trace_{part:04d}.npz'
        np.savez_compressed(ROOT/path,trace=trace[start:start+8192])
        paths.append(path)
    return dict(trace_parts=paths,trace_samples=len(trace),
                trace_float64_sha256=hashlib.sha256(trace.tobytes()).hexdigest())


def initial(c,output):
    amp,w,L=c['amplitude'],c['width'],c['length']
    def f(r):return amp*np.exp(-r*r/(2*w*w))
    def rhs(r,y):return np.array([y[1],(3*(y[0]**2-1)+f(r)**2)*y[0]])
    def bc(left,right):return [left[1],right[0]-1]
    x=np.linspace(0,L,601)
    guess=np.array([np.sqrt(1-f(x)**2/3),f(x)**2*x/(3*w*w)])
    s=solve_bvp(rhs,bc,x,guess,S=np.diag([0.,-2.]),tol=1e-10,bc_tol=1e-11,max_nodes=12000)
    assert s.success,s.message
    r=(np.arange(round(L/c['h']))+.5)*c['h'];F=f(r);S=s.sol(r)[0]
    # Exact full-field Q equals the NR norm at this map: pi_X=-i X.
    a=np.column_stack([r,F,np.zeros_like(F),S,np.zeros_like(F),-F,np.zeros_like(F)])
    np.savetxt(output,a,delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
    rr=np.linspace(0,L,20001);ff=f(rr);ss,sp=s.sol(rr)
    grad=3*(np.pi**1.5*amp**2*w**3)/(4*w*w)
    charge=np.pi**1.5*amp**2*w**3
    had=4*np.pi*simpson(rr*rr*(.5*(ff*rr/w**2)**2+.5*sp**2+
        .75*(ss*ss-1)**2+.5*(ss*ss-1)*ff*ff+.025*ff**4),x=rr)
    return dict(initial_NR_norm_Qbar=charge,initial_full_energy_continuum=charge+had,
                adiabatic_NR_Hamiltonian=had,NR_energy_including_rest=charge+had,
                s_elliptic_max_residual=float(max(s.rms_residuals)),
                initial_max_phi_depletion=float(1-min(S)),
                gaussian_focusing_parameter=(.01/.03-.1)*amp*amp*w*w,
                actual_cosmological_preparation=False)


def summarize(name,c,init,elapsed):
    p=ROOT/'data'/name
    hist=np.loadtxt(str(p)+'_history.csv',delimiter=',',skiprows=1)
    final=np.loadtxt(str(p)+'_final.csv',delimiter=',',skiprows=1)
    trace=np.loadtxt(str(p)+'_trace.csv',delimiter=',',skiprows=1)
    storage=store_trace(name,trace)
    Path(str(p)+'_trace.csv').unlink()
    e0,q0=hist[0,1:3]
    late=hist[:,0]>.8*c['tmax']
    row=dict(name=name,parameters=c,initial=init,elapsed_seconds=elapsed,
             Ebar_initial=float(e0),Qbar_initial=float(q0),Qphysical_initial=float(q0/.01),
             Ebar_final=float(hist[-1,1]),Qbar_final=float(hist[-1,2]),
             max_total_energy_relative_error=float(max(abs((hist[:,1]+hist[:,3])/e0-1))),
             max_total_charge_relative_error=float(max(abs((hist[:,2]+hist[:,4])/q0-1))),
             maximum_F=float(max(hist[:,7])),minimum_S=float(min(hist[:,8])),
             Qbar_absorbed=float(hist[-1,4]),Ebar_absorbed=float(hist[-1,3]),
             cumulative_Ebar_flux_probe=float(hist[-1,5]),cumulative_Qbar_flux_probe=float(hist[-1,6]),
             final_core_charge_fractions={str(r):float(hist[-1,10+2*i]/q0) for i,r in enumerate([40,60,80])},
             late_mean_core_charge_fractions={str(r):float(np.mean(hist[late,10+2*i])/q0) for i,r in enumerate([40,60,80])},
             late_core60_charge_fraction_range=[float(min(hist[late,12])/q0),float(max(hist[late,12])/q0)],
             final_core_energy_fractions={str(r):float(hist[-1,9+2*i]/e0) for i,r in enumerate([40,60,80])},
             **storage,history_file=f'data/{name}_history.csv',final_file=f'data/{name}_final.csv')
    return row


def main(names):
    dest=ROOT/'CORE_RESULTS.json'
    result=json.loads(dest.read_text()) if dest.exists() else dict(
        scope='Prepared spherical hierarchy-reduced full-field prototypes; no final cosmological dN/dQ, quantum decay, or 3D fragmentation claim.',runs=[])
    with tempfile.TemporaryDirectory(prefix='grut-radial-build-') as tmp:
        binary=Path(tmp)/'radial'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'code/radial.cpp'),'-o',str(binary)],check=True)
        for name in names:
            c=CASES[name];inputfile=Path(tmp)/'initial.csv';init=initial(c,inputfile)
            start=time.monotonic()
            subprocess.run([str(binary),str(inputfile),str(ROOT/'data'/name),
                 *[str(c[k]) for k in ['h','dt','length','tmax','probe','absorber','eta']]],check=True)
            row=summarize(name,c,init,time.monotonic()-start)
            result['runs']=[x for x in result['runs'] if x['name']!=name]+[row]
            result['source_sha256']=hashlib.sha256((ROOT/'code/radial.cpp').read_bytes()).hexdigest()
            dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            print(json.dumps(row,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('names',nargs='*',default=['A'])
    main(p.parse_args().names)
