#!/usr/bin/env python3
"""Compare earlier H55 phase preparations on the full homogeneous scalar ledger."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import trapezoid

TRACKS=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(TRACKS/'x5_coupled/code'))
from run import DT, digest, read_columns, write_input

X5=TRACKS/'x5_coupled'
SETTINGS={'coarse':(2880,DT,2),'radial':(5760,DT,2),'half':(2880,DT/2,4)}
OBSERVABLES=('fermion_fraction','fermion_energy','fermion_number','comoving_fermion_number',
             'n_1','n_2','n_3','n_4','a3_n_1','a3_n_2','a3_n_3','a3_n_4')


def diagnostics(d):
    rho=d['rho']
    t=d['t']
    energy=sum(d[f'rho_{j}'] for j in range(1,5))
    number=sum(d[f'n_{j}'] for j in range(1,5))
    power=sum(d[f'J_{j}']*d['dsigma' if j==1 else 'dchi'] for j in range(1,5))
    return dict(t=t,rho=rho,fermion_fraction=energy/rho,fermion_energy=energy,
                fermion_number=number,comoving_fermion_number=d['a']**3*number,
                normalized_work=power/rho,
                **{f'n_{j}':d[f'n_{j}'] for j in range(1,5)},
                **{f'a3_n_{j}':d['a']**3*d[f'n_{j}'] for j in range(1,5)})


def average(t,y,begin,end):
    interior=(t>begin)&(t<end)
    x=np.r_[begin,t[interior],end]
    return float(trapezoid(np.interp(x,t,y),x)/(end-begin))


def summarize(d,match):
    h=diagnostics(d)
    t=d['t']
    return dict(
        history_end=float(t[-1]),
        max_budget_error=float(np.max(abs(d['residual']))),
        max_bloch_norm=float(np.max(d['norm_max'])),
        pauli_min=float(np.min(d['f_min'])),pauli_max=float(np.max(d['f_max'])),
        change_at_copy={key:float(value[-1]-value[0]) for key,value in match.items()},
        normalized_work_at_copy_before=float(diagnostics(match)['normalized_work'][0]),
        normalized_work_at_copy_after=float(diagnostics(match)['normalized_work'][-1]),
        initial_gas_fraction_at_copy=float(diagnostics(match)['fermion_fraction'][0]),
        terminal={key:float(h[key][-1]) for key in OBSERVABLES},
        final_window={key:average(t,h[key],t[-1]-0.1,t[-1]) for key in OBSERVABLES},
        peak_fermion_fraction=float(np.max(h['fermion_fraction'])),
    ),h


def compare(baseline,reflect):
    a,x=baseline
    b,y=reflect
    final={key:b['terminal'][key]/a['terminal'][key]-1 for key in OBSERVABLES}
    window={key:b['final_window'][key]/a['final_window'][key]-1 for key in OBSERVABLES}
    take=(x['t']>=0.16)&(x['t']<=0.42)
    xt=x['t'][take]
    phase=abs(np.interp(xt,y['t'],y['normalized_work'])-x['normalized_work'][take])
    fraction=abs(np.interp(xt,y['t'],y['fermion_fraction'])-x['fermion_fraction'][take])
    return dict(terminal_relative_difference=final,final_window_relative_difference=window,
                early_max_work_difference_over_rho=float(np.max(phase)),
                early_max_fermion_fraction_difference=float(np.max(fraction)),
                early_fraction_difference_time=float(xt[np.argmax(fraction)]),
                final_window_within_one_percent=bool(max(map(abs,window.values()))<0.01))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--only',choices=SETTINGS)
    p.add_argument('--gather',action='store_true')
    args=p.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    cases=json.loads((X5/'results/h55_refined.json').read_text())['cases']
    base=next(row['model'] for row in cases if row['case']=='eight_bins')
    engine=work/'production'
    if not args.gather:
        subprocess.run(['gcc','-std=c11','-O3','-march=native','-fopenmp','-Wall','-Wextra','-Werror',
                        str(Path(__file__).with_name('production.c')),'-lm','-o',str(engine)],check=True)
    names=(args.only,) if args.only else SETTINGS
    tasks=[(name,reflect) for name in names for reflect in (False,True)]

    def run(task):
        name,reflect=task
        nq,target,fine_stride=SETTINGS[name]
        label=name+('_reflect' if reflect else '_base')
        model=dict(base,nq=nq,nsteps=math.ceil((base['t1']-base['t0'])/target),stride=64 if target==DT else 128)
        input_path,history,match=(work/f'{label}.{ext}' for ext in ('input','history','match'))
        if not args.gather:
            write_input(input_path,model)
            subprocess.run([str(engine),str(input_path),str(history),str(match),str(int(reflect)),str(fine_stride)],
                           check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
        d=read_columns(history)
        if abs(d['t'][-1]-model['t1'])>1e-12:
            raise ValueError(f'Truncated production: {label}')
        result,h=summarize(d,read_columns(match))
        result.update(radial_bins=nq,dt_target=target,history_sha256=digest(history))
        print(label,'final fermion energy fraction',result['terminal']['fermion_fraction'],flush=True)
        return task,(result,h)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=dict(pool.map(run,tasks))
    comparisons={name:compare(results[name,False],results[name,True]) for name in names}
    resolutions={}
    if 'coarse' in names:
        coarse=results['coarse',False][0]['final_window']
        resolutions={name:{key:results[name,False][0]['final_window'][key]/coarse[key]-1
                           for key in OBSERVABLES} for name in names if name!='coarse'}
    report=dict(scope='homogeneous H55 scalar and sterile feedback with collisionless sterile fields',
                copy_time=0.16,radial_reflection_endpoints=[3000.,6000.],final_window=0.10,
                cases={name+('_reflect' if reflect else '_base'):value[0]
                       for (name,reflect),value in results.items()},
                comparisons=comparisons,base_final_window_resolution_differences=resolutions,
                scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
