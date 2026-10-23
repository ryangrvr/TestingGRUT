#!/usr/bin/env python3
"""Fourth-order time refinement at the unchanged finite-volume Hamiltonian."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from scipy.integrate import solve_bvp
import subprocess
import tempfile
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_daughter_core_20261023'
spec=importlib.util.spec_from_file_location('parent_core', PARENT/'code/run_core.py')
parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
BASE=dict(amplitude=.18,width=60.,length=360.,tmax=12000.,probe=120.,absorber=270.,eta=.03,order=4)
CASES={
 'H025':dict(h=.25,dt=.0625),
 'H025_time':dict(h=.25,dt=.03125),
 'H0125':dict(h=.125,dt=.03125),
 'H0125_time':dict(h=.125,dt=.015625),
 'H00625':dict(h=.0625,dt=.015625),
 'H003125':dict(h=.03125,dt=.0078125),
 'H0125_box':dict(h=.125,dt=.03125,length=480.,absorber=360.),
 'H0125_abs':dict(h=.125,dt=.03125,eta=.015),
 'H0125_late':dict(h=.125,dt=.03125,tmax=24000.),
 'N05':dict(h=.5,dt=.125,nodal=True),
 'N025':dict(h=.25,dt=.0625,nodal=True),
 'N0125':dict(h=.125,dt=.03125,nodal=True),
 'N00625':dict(h=.0625,dt=.015625,nodal=True),
 'N0125_time':dict(h=.125,dt=.015625,nodal=True),
 'N0125_box':dict(h=.125,dt=.03125,nodal=True,length=480.,absorber=360.),
 'N025_time':dict(h=.25,dt=.03125,nodal=True),
 'N0125_abs':dict(h=.125,dt=.03125,nodal=True,eta=.015),
 'N0125_far':dict(h=.125,dt=.03125,nodal=True,length=480.,absorber=360.,probe=240.),
}


def run(name):
    c={**BASE,**CASES[name]};prefix=ROOT/'data'/name
    with tempfile.TemporaryDirectory(prefix='grut-core-refine-') as tmp:
        tmp=Path(tmp);binary=tmp/'radial4';init=tmp/'initial.csv'
        source=ROOT/'code'/('radial_nodal.cpp' if c.get('nodal') else 'radial4.cpp')
        subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(binary)],check=True)
        initial=parent.initial(c,init)
        if c.get('nodal'):
            # Solve the same smooth neutral initial data, evaluated on nodes.
            amp,width,L=c['amplitude'],c['width'],c['length']
            f=lambda r:amp*np.exp(-r*r/(2*width*width))
            rhs=lambda r,y:np.array([y[1],(3*(y[0]**2-1)+f(r)**2)*y[0]])
            x=np.linspace(0,L,601)
            guess=np.array([np.sqrt(1-f(x)**2/3),f(x)**2*x/(3*width*width)])
            s=solve_bvp(rhs,lambda l,r:[l[1],r[0]-1],x,guess,S=np.diag([0.,-2.]),tol=1e-10,bc_tol=1e-11,max_nodes=12000)
            assert s.success
            r=np.arange(1,round(L/c['h']))*c['h'];F=f(r);zeros=np.zeros_like(r)
            np.savetxt(init,np.column_stack([r,F,zeros,s.sol(r)[0],zeros,-F,zeros]),delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
        start=time.monotonic()
        keys=['h','dt','length','tmax','probe','absorber','eta']+([] if c.get('nodal') else ['order'])
        proc=subprocess.run([str(binary),str(init),str(prefix),
             *[str(c[k]) for k in keys]],
             check=True,capture_output=True,text=True)
        elapsed=time.monotonic()-start
    history=np.loadtxt(str(prefix)+'_history.csv',delimiter=',',skiprows=1)
    trace=np.loadtxt(str(prefix)+'_trace.csv',delimiter=',',skiprows=1)
    # Keep a lossless trace in a single compressed file, retaining array hash.
    np.savez_compressed(str(prefix)+'_trace.npz',trace=trace)
    Path(str(prefix)+'_trace.csv').unlink()
    e0,q0=history[0,1:3];late=history[:,0]>.8*c['tmax']
    row=dict(name=name,parameters=c,initial=initial,seconds=elapsed,
        Ebar_initial=float(e0),Qbar_initial=float(q0),
        max_energy_error=float(np.max(abs((history[:,1]+history[:,3])/e0-1))),
        max_charge_error=float(np.max(abs((history[:,2]+history[:,4])/q0-1))),
        core60_late_mean_fraction=float(np.mean(history[late,12])/q0),
        core60_final_fraction=float(history[-1,12]/q0),
        Qbar_absorbed=float(history[-1,4]),Ebar_absorbed=float(history[-1,3]),
        trace_samples=len(trace),trace_sha256=hashlib.sha256(trace.tobytes()).hexdigest(),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        executable_stdout=proc.stdout.strip())
    (ROOT/f'{name}_RESULTS.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    print(json.dumps(row,indent=2),flush=True)
    return row


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('names',nargs='+',choices=list(CASES));args=p.parse_args()
    for name in args.names:run(name)
