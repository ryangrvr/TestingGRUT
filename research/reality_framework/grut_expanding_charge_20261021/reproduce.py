"""Replay a selected executed run into scratch, preserving archived reference files."""
from pathlib import Path
import argparse,importlib.util,json
import numpy as np
HERE=Path(__file__).resolve().parent

def module(name,file):
    s=importlib.util.spec_from_file_location(name,HERE/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    p=argparse.ArgumentParser();p.add_argument('--list',action='store_true');p.add_argument('--run')
    args=p.parse_args()
    ref=json.loads((HERE/'RESULTS.json').read_text())
    names=[x['name'] for x in ref['status'] if x['status']=='completed']
    if args.list or args.run is None:print('\n'.join(names));return
    selected=names if args.run=='all' else [args.run]
    if not set(selected)<=set(names):raise SystemExit('Run is not in the executed reference list')
    R=module('replay_frozen','run.py');G=module('replay_gaussian','gaussian_run.py')
    target=HERE/'scratch'/'reproduced';target.mkdir(parents=True,exist_ok=True)
    R.RAW=target;G.R.RAW=target
    J=dict(R.jobs());R.compile_solver();G.compile_solver()
    for name in selected:
        x=ref['runs'][name]
        if name.startswith('gaussian'):G.run(name,x['N'],x['dx'],x['seed'],dt=x['dt'])
        else:R.run(name,T=48,sample=1.2,**{'top-sample':4.8},**J[name])
        reference=HERE/'raw'/(name+'.csv')
        if reference.exists():
            a=np.atleast_1d(np.genfromtxt(reference,delimiter=',',names=True));b=np.atleast_1d(np.genfromtxt(target/(name+'.csv'),delimiter=',',names=True))
            if len(a)!=len(b):print(name,'time-row count differs');continue
            scaled={k:float(np.max(np.abs(a[k]-b[k])))/max(1.,float(np.max(np.abs(a[k])))) for k in a.dtype.names}
            print(name,'maximum column-scaled absolute difference',max(scaled.values()))

if __name__=='__main__':main()
