#!/usr/bin/env python3
"""Stationary width controls and a locally constant rotating-kernel H55 sensitivity."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import roots_legendre

X7=Path(__file__).resolve().parents[2]/'x7_regulator'
sys.path.insert(0,str(X7/'code'))
from overlap import DURATION, GAPS, initial_data, measure, read_columns, summarize, write_input


def controls(engine):
    table=np.loadtxt(io.StringIO(subprocess.check_output([str(engine),'--controls'],text=True)))
    nodes,weights=roots_legendre(32)
    M,G=2.,0.04
    reference=[]
    for k,delta,omega,*_ in table:
        p=abs(k)
        u=np.array([M,0,k+delta])/omega
        r=np.array([.2,-.3,.6])
        b=float(r @ u)
        width=[]
        for sign in (1,-1):
            energy=omega-sign*delta
            rate=0.
            if energy>p:
                # Daughter lepton current after boosting its center-of-mass sphere.
                lepton_e=(energy+p*nodes)/2
                lepton_z=(p+energy*nodes)/2
                average=np.dot(weights,lepton_e-sign*np.sign(k)*lepton_z)
                rate=G/M*(1+sign*u[2])/2*average
            width.append(rate)
        gamma=sum(width)
        lost=-np.expm1(-2*gamma)*(1+b)/2
        n=np.exp(-gamma)*(r-b*u)+(b-2*lost)*u
        reference.append([*width,*n,lost])
    reference=np.array(reference)
    observed=table[:,[3,4,7,8,9,10]]
    zero=table[:,1]==0
    hflip=max(np.max(abs(row[[3,4]]-table[
        (table[:,0]==-row[0]) & (table[:,1]==-row[1]),:][0,[4,3]])) for row in table)
    anti=np.array([np.sign(row[0])*row[1]<0 and
                   abs(row[0])>M*M/(4*abs(row[1])) for row in table])
    return dict(probes=len(table),width_max_abs_error=float(np.max(abs(reference[:,:2]-observed[:,:2]))),
                bloch_max_abs_error=float(np.max(abs(reference[:,2:]-observed[:,2:]))),
                zero_drive_max_abs_error=float(np.max(abs(table[zero,3]+table[zero,4]-
                                                          G*M/np.hypot(table[zero,0],M)))),
                helicity_drive_reversal_error=float(hflip),
                anti_aligned_threshold_matches=bool(np.array_equal(
                    anti,np.count_nonzero(table[:,3:5],axis=1)==1)),
                max_total_width=float(np.max(table[:,3]+table[:,4])),
                max_bloch_norm=float(np.max(np.linalg.norm(observed[:,2:5],axis=1))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    model,terminal,_,state,sha,widths,_=initial_data()
    checkpoint=work/'checkpoint.bin'
    with checkpoint.open('wb') as f:
        for key in ('q','w','bloch'):
            np.asarray(state[key],dtype='=f8',order='C').tofile(f)
    engine=work/'rotating'
    source=Path(__file__).with_name('rotating.c')
    subprocess.run(['gcc','-std=c11','-O3','-march=native','-Wall','-Wextra','-Werror',
                    '-fopenmp',str(source),'-lm','-o',str(engine)],check=True)
    bench=controls(engine)
    print('Stationary controls',bench,flush=True)
    base=math.ceil(DURATION*model['nsteps']/(model['t1']-model['t0']))

    def run(multiplier):
        name='rotating' if multiplier==1 else 'rotating_half'
        species=[[*row[:3],width] for row,width in zip(model['species'],widths,strict=True)]
        cfg=dict(model,species=species,t0=terminal['t'],t1=terminal['t']+DURATION,
                 nsteps=multiplier*base,stride=4*multiplier,a0=terminal['a'],H0=terminal['H'],
                 fields=[terminal[k] for k in ('sigma','chi','s')],
                 velocities=[terminal[k] for k in ('dsigma','dchi','ds')],
                 rho_sm=terminal['rho_sm'],rho_dr=terminal['rho_dr'])
        paths=[work/f'{name}.{suffix}' for suffix in ('input','history','projection','audit')]
        write_input(paths[0],cfg)
        subprocess.run([str(engine),str(paths[0]),str(checkpoint),*(str(x) for x in paths[1:])],
                       check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
        d=read_columns(paths[1])
        projected=np.loadtxt(paths[2],ndmin=2)
        trace=measure(d,projected,[row[0] for row in species])
        summary=summarize(d,trace)
        summary['dt']=DURATION/cfg['nsteps']
        summary['timestep_restart_error']=max(abs(d[k][0]-terminal[k])/max(1.,abs(terminal[k]))
                                            for k in ('a','H','rho','pressure','dsigma','dchi','rho_1','rho_4'))
        audit=np.loadtxt(paths[3])
        summary['collision_substep']=dict(energy_error=float(audit[0]),charge_error=float(audit[1]),
                                          midpoint_bisections=int(audit[2]))
        volume=d['a']**4
        power=projected[:,6]
        total=float(trapezoid(power*volume,d['t']))
        indicators={}
        for key,idx in (('eigenbasis_and_secular',8),('bath_threshold',9),('all',10),('both_branches_open',11)):
            selected=projected[:,idx]
            indicators[key]=dict(initial=float(selected[0]/power[0]),
                                 integrated=float(trapezoid(selected*volume,d['t'])/total),
                                 min=float(np.min(selected/power)),max=float(np.max(selected/power)))
        radiation=(d['a'][-1]**4*d['rho_sm'][-1]-d['a'][0]**4*d['rho_sm'][0])
        summary['validity']=indicators
        summary['radiation_source_integral_error']=float(abs(total/radiation-1))
        summary['max_lab_width']=float(np.max(projected[:,13]))
        summary['net_branch_rate_initial_over_total_energy']=float(projected[0,7]/d['rho'][0])
        np.savez_compressed(work/f'{name}.npz',**trace)
        print(name,'gate',summary['gate_pass'],'radiation',summary['final']['radiation_sm_fraction'],
              'locally slow power fraction',indicators['all']['integrated'],flush=True)
        return name,(summary,trace)

    with ThreadPoolExecutor(max_workers=2) as pool:
        runs=dict(pool.map(run,(1,2)))
    a,b=(runs[name][1] for name in ('rotating','rotating_half'))
    old=json.loads((X7/'results/overlap.json').read_text())['cases']
    comparison={}
    for name in ('coherent_half','vacuum_decay_half'):
        comparison[name]={key:runs['rotating_half'][0]['final'][key]-old[name]['final'][key]
                          for key in ('fermion_fraction','radiation_sm_fraction',*GAPS)}
    output=dict(scope='diagonal locally constant vacuum Born sensitivity with prescribed mass-phase Yukawa frame',
                checkpoint_sha256=sha,rest_widths=widths,stationary_controls=bench,
                cases={k:s for k,(s,_) in runs.items()},validity_threshold=0.1,
                timestep_max_abs_differences={k:float(np.max(abs(a[k]-b[k]))) for k in (*GAPS,'fermion_fraction',
                                                                                       'radiation_sm_fraction')},
                final_differences_against_X7=comparison,scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
