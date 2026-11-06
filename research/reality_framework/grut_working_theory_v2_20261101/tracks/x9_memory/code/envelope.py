#!/usr/bin/env python3
"""Bound sterile phase sensitivity for every sign pattern on one prescribed history."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import subprocess

from production import DT, X5, read_columns, write_input


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    models=json.loads((X5/'results/h55_refined.json').read_text())['cases']
    base=next(row['model'] for row in models if row['case']=='eight_bins')
    engine=work/'envelope'
    subprocess.run(['gcc','-std=c11','-O3','-march=native','-fopenmp','-Wall','-Wextra','-Werror',
                    str(Path(__file__).with_name('envelope.c')),'-lm','-o',str(engine)],check=True)
    settings={'coarse':(1440,DT),'fine':(2880,DT),'half':(1440,DT/2)}

    def run(name):
        nq,dt=settings[name]
        cfg=dict(base,nq=nq,nsteps=math.ceil((base['t1']-base['t0'])/dt))
        infile,outfile=(work/f'{name}.{ext}' for ext in ('input','output'))
        write_input(infile,cfg)
        subprocess.run([str(engine),str(infile),str(outfile)],check=True,
                       env=dict(os.environ,OMP_NUM_THREADS='1'))
        raw=read_columns(outfile)
        rows=[{key:float(value[i]) for key,value in raw.items()} for i in range(len(raw['species']))]
        bounds=rows[:-1]
        baseline=rows[-1]
        keys=('endpoint_number','endpoint_energy','mean_number','mean_energy')
        relative={key:sum(row[key] for row in bounds)/baseline[key] for key in keys}
        print(name,'final-window phase bounds',relative['mean_number'],relative['mean_energy'],flush=True)
        return name,dict(radial_bins=nq,dt_target=dt,per_species_absolute_bounds=bounds,
                         baseline=baseline,total_relative_bounds=relative)

    with ThreadPoolExecutor(max_workers=3) as pool:
        cases=dict(pool.map(run,settings))
    result=dict(scope='exact linear phase envelope with the common X5 midpoint scalar history held fixed',
                copy_time=0.16,final_window=0.10,cases=cases,scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
