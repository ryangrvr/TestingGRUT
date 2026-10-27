"""Predeclared cold nonlinear controls, ordinary CPU workers only."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pilot import *

def task(spec):
    tag,N,L,seed,neutral,dt=spec
    try:return tag,run(tag,N=N,L=L,seed=seed,neutral=neutral,dt=dt,T=240.,snap=20.,sample=1.,threads=4,bits=64)
    except Exception as e:return tag,dict(error=repr(e),complete=False)

def main():
    checks=json.loads((ROOT/'ENGINE_CHECKS.json').read_text());assert all(v['passed'] for v in checks['groups'].values())
    build()
    specs=[('mixed_N96',96,48.,1,'mixed',.06),('onesign_N96',96,48.,1,'onesign',.06),
        ('mixed_N64',64,48.,1,'mixed',.06),('mixed_N128',128,48.,1,'mixed',.045),
        ('mixed_dt_half',96,48.,1,'mixed',.03),('mixed_seed2',96,48.,2,'mixed',.06),
        ('mixed_box64',128,64.,1,'mixed',.06)]
    results={}
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(task,s) for s in specs]
        for f in as_completed(futures):
            tag,r=f.result();results[tag]=r;save('RUNS.json',dict(cases=results,physical_formation_derived=False,physical_f_sym=None))
            print(tag,{k:v for k,v in r.items() if k not in ['parameters','snapshots']},flush=True)

if __name__=='__main__':main()
