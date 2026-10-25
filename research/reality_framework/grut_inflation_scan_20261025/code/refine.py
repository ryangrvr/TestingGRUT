"""Refine the regular/chaotic transition and every high-fraction release.

Adaptive rule written after the coarse scan: cH=1 grid has two candidates
near eps=0.15-0.19 at H=8.72e12. cH=10 loses phase robustness near 2e12.
Cover both transition regions broadly; do not optimize only one endpoint.
"""
import argparse,csv,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from run import ROOT,write
from scan import worker,compact
from build import build


def main():
    a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=8);args=a.parse_args()
    build()
    grids={1.:np.geomspace(5e11,1.4e13,31),10.:np.geomspace(1.2e11,3e12,31)}
    jobs=[dict(H=float(H),cH=c,delta=j*np.pi/24) for c,grid in grids.items() for H in grid for j in range(1,24)]
    config=dict(stage='refinement of coarse scan',H_grids={str(k):v.tolist() for k,v in grids.items()},
                phase_indices=list(range(1,24)),phase_unit='pi/24',n_unshifted=len(jobs),
                shift_rule='Add delta+1e-9 for every abs(epsX)>=0.14 and abs(YX/Y_Q)>=1.',
                additional_initial_condition='Same full aligned radial minimum and supplied source-phase mismatch as coarse scan.')
    write(ROOT/'inputs/REFINEMENT_CONFIG.json',config)
    rows=[];t0=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(worker,j) for j in jobs]
        for f in as_completed(futures):
            rows.append(compact(f.result()))
            if len(rows)%100==0:print('REFINE',len(rows),'/',len(jobs),'elapsed',round(time.perf_counter()-t0,1),flush=True)
        shifts=[dict(H=r['H_GeV'],cH=r['cH'],delta=r['delta']+1e-9) for r in rows
                if abs(r['epsX'])>=.14 and abs(r['YX_over_YQ'])>=1]
        shifted={}
        for f in as_completed([pool.submit(worker,j) for j in shifts]):
            r=compact(f.result());shifted[(r['cH'],r['H_GeV'],r['delta'])]=r
    candidates=[]
    for r in rows:
        s=shifted.get((r['cH'],r['H_GeV'],r['delta']+1e-9))
        r['phase_shift_abs_delta_epsX']=abs(r['epsX']-s['epsX']) if s else None
        r['shifted_epsX']=s['epsX'] if s else None
        r['shifted_YX_over_YQ']=s['YX_over_YQ'] if s else None
        r['provisional_regular']=bool(s and abs(r['epsX']-s['epsX'])<1e-4 and r['late_charge_drift_fraction']<.01)
        if r['provisional_regular']:candidates.append(r)
    rows.sort(key=lambda r:(r['cH'],r['H_GeV'],r['delta']))
    with (ROOT/'REFINEMENT_RESULTS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    with (ROOT/'REFINEMENT_SHIFT_RESULTS.csv').open('w',newline='') as f:
        sr=list(shifted.values());w=csv.DictWriter(f,fieldnames=list(sr[0]),lineterminator="\n");w.writeheader();w.writerows(sr)
    write(ROOT/'REFINED_CANDIDATES.json',dict(n_unshifted=len(rows),n_shifted=len(shifts),
        max_abs_epsX_over_all_rows=max(abs(r['epsX']) for r in rows),
        scope='Favorable supplied phase and numerical screening only. Seed ensembles and tighter tolerances still required.',
        rows=sorted(candidates,key=lambda r:abs(r['epsX']),reverse=True),runtime_seconds=time.perf_counter()-t0))
    print('REFINEMENT_COMPLETE',len(rows),len(shifts),'candidates',len(candidates),'seconds',time.perf_counter()-t0,flush=True)


if __name__=='__main__':main()
