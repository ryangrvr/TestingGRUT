"""Seven regulated scalar-limit histories, with unchanged matched seed functions."""
from en_common import *
from run import run
from concurrent.futures import ThreadPoolExecutor,as_completed

def main():
    jobs=[dict(tag='N32',N=32),dict(tag='N48',N=48),dict(tag='N64',N=64),
          dict(tag='halfstep48',N=48,hmax=.002,sample=.025),
          dict(tag='seed2_48',N=48,seed=2),dict(tag='cut10_48',N=48,cutoff=10.),dict(tag='box16_64',N=64,L=16.)]
    rows={}
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(run,**job):job['tag'] for job in jobs}
        for f in as_completed(futures):
            tag=futures[f]
            try:row=f.result()
            except Exception as exc:
                rows[tag]=dict(error=repr(exc));save('RUNS.json',dict(cases=rows));raise
            rows[tag]=row;save('RUNS.json',dict(cases=rows,physical_f_sym=None,full_model_evolved=False))
            print('RUN',tag,'seconds',row['seconds'],'energy',row['max_recorded_energy_relative_drift'],'Q',row['max_recorded_charge_relative_drift'],flush=True)

if __name__=='__main__':main()
