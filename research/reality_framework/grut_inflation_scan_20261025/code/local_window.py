"""Resolve the newly found supplied-phase counterexample to a no-window claim."""
import csv
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from run import ROOT,write
from scan import worker,compact


def main():
    center=11211167304108.865
    hs=center*np.linspace(.98,1.02,41)
    phases=np.pi/2+np.array([-.06,-.03,0.,.03,.06])
    jobs=[dict(H=float(h),cH=1.,delta=float(d),rtol=2e-12,atol=2e-14) for h in hs for d in phases]
    write(ROOT/'inputs/LOCAL_WINDOW_CONFIG.json',dict(H_center_GeV=center,
        H_factors=np.linspace(.98,1.02,41).tolist(),phases=phases.tolist(),n_runs=len(jobs),
        scope='Conditional supplied-phase map near the refined high-fraction point; excludes Y/vector debris and phase preparation.'))
    with ProcessPoolExecutor(max_workers=8) as pool:rows=[compact(r) for r in pool.map(worker,jobs)]
    with (ROOT/'LOCAL_WINDOW_RESULTS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    passing=[]
    for h in hs:
        rr=[r for r in rows if r['H_GeV']==h]
        if all(abs(r['epsX'])>=.19 and abs(r['YX_over_YQ'])>=1 for r in rr):passing.append(float(h))
    write(ROOT/'LOCAL_WINDOW_SUMMARY.json',dict(n_runs=len(rows),H_values_all_five_phases_above_019=passing,
        minimum_Y_pair_estimate_over_all_rows=min(r['approximate_Y_decay_pairs_per_abs_X_charge'] for r in rows),
        max_abs_eps=max(abs(r['epsX']) for r in rows),
        passing_count=int(sum(abs(r['epsX'])>=.19 for r in rows))))
    print('LOCAL_WINDOW_COMPLETE',len(rows),'H slices all phases above threshold',len(passing),flush=True)


if __name__=='__main__':main()
