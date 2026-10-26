"""Debris-inclusive case-L outcomes and density-dependent numerical thresholds."""
import json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from run import ROOT,write
from relic import run


def outcome(job):
    label,kwargs=job;r=run(**kwargs);return label,r


def threshold(job):
    index,epsilon,ymult=job
    lo,hi=0.,400.;history=[]
    for _ in range(11):
        mid=(lo+hi)/2;r=run(epsilon,mid,ymult,steps=520)
        history.append({k:r[k] for k in ['R','f_extra','entropy_per_charge','max_cond_frac','max_thermo_constraint_error']})
        if r['f_extra']>.01:hi=mid
        else:lo=mid
    # Repeat around the first threshold at doubled step resolution.
    tightlo=max(0,lo*.9);tighthi=hi*1.1
    for _ in range(8):
        mid=(tightlo+tighthi)/2;r=run(epsilon,mid,ymult,steps=1040)
        history.append({k:r[k] for k in ['R','f_extra','entropy_per_charge','max_cond_frac','max_thermo_constraint_error']}|dict(steps=1040))
        if r['f_extra']>.01:tighthi=mid
        else:tightlo=mid
    return dict(candidate=index,epsilon=epsilon,ymult=ymult,R_L_520_bracket=[lo,hi],
                R_L_1040_bracket=[tightlo,tighthi],evaluations=history)


def main():
    source=json.loads((ROOT/'DECAY_RESULTS.json').read_text())['rows'];jobs=[]
    for s in source:
        for density in [1.,s['ymult_raw']]:
            for name in ['zero','R_Y','R_lower','R_central','R_upper']:
                R=0. if name=='zero' else s[name]
                jobs.append((f"{s['candidate']}_{density:.6g}_{name}",dict(epsilon=s['epsilon_net'],R=R,ymult=density,steps=1040)))
        for n in [260,520]:
            jobs.append((f"{s['candidate']}_step{n}",dict(epsilon=s['epsilon_net'],R=s['R_lower'],ymult=s['ymult_raw'],steps=n)))
        jobs.append((f"{s['candidate']}_quad128",dict(epsilon=s['epsilon_net'],R=s['R_lower'],ymult=s['ymult_raw'],steps=520,nodes=128)))
        jobs.append((f"{s['candidate']}_late",dict(epsilon=s['epsilon_net'],R=s['R_lower'],ymult=s['ymult_raw'],steps=1040,endT=1e-9)))
        jobs.append((f"{s['candidate']}_legacy",dict(epsilon=s['epsilon_net'],R=s['R_lower'],ymult=s['ymult_raw'],steps=260,endT=1e-6,scheme='review_endpoint',table='review')))
        for factor in [.1,10.]:
            jobs.append((f"{s['candidate']}_rate{factor}",dict(epsilon=s['epsilon_net'],R=s['R_lower'],ymult=s['ymult_raw'],steps=520,rate_factor=factor)))
    results={};t=time.perf_counter()
    with ProcessPoolExecutor(max_workers=8) as pool:
        for future in as_completed([pool.submit(outcome,j) for j in jobs]):
            key,r=future.result();results[key]=r
            write(ROOT/'BENCHMARK_RESULTS.json',results)
            print('BENCHMARK',len(results),'/',len(jobs),key,'f',r['f_extra'],flush=True)
        thjobs=[(s['candidate'],s['epsilon_net'],ym) for s in source for ym in [1.,s['ymult_raw']]]
        thresholds=[]
        for future in as_completed([pool.submit(threshold,j) for j in thjobs]):
            thresholds.append(future.result());write(ROOT/'THRESHOLD_RESULTS.json',thresholds)
            print('THRESHOLD',thresholds[-1]['candidate'],thresholds[-1]['ymult'],thresholds[-1]['R_L_1040_bracket'],flush=True)
    print('BENCHMARK_COMPLETE',time.perf_counter()-t,flush=True)


if __name__=='__main__':main()
