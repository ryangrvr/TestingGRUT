"""Predeclared full-field H/phase grid, including the 1e-9 shift pairs.

Cached trajectories are temporary reproducible computation products.
The packet publishes every grid point's sufficient numerical diagnostics
and canonical endpoint, plus selected complete sampled histories.
"""
import argparse
import csv
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from build import build,WORK
from model import minimum, initial_z
from run import ROOT,evolve,write,NAMES

ENGINE_HASH=hashlib.sha256(b''.join((ROOT/'code'/p).read_bytes() for p in ['release.c','model.py','run.py'])).hexdigest()
CACHE=WORK/'cache'/ENGINE_HASH[:16]


def key(job):
    return hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:24]


def worker(job):
    path=CACHE/(key(job)+'.json')
    if path.exists():return json.loads(path.read_text())
    r=evolve(**job)
    write(path,r)
    return r


def compact(r):
    w=r['windows']['600-1000'];earlier=r['windows']['300-600']
    keep=['NX','NY','NS','Neff','NX_std','Neff_std','epsX','epsX_std','eps_from_mean_charge',
          'YX','Yeff','YX_over_YQ','Yeff_over_YQ','EX','EY','ES','EY_exc_quadratic','JY_approx',
          'approximate_Y_decay_pairs_per_abs_X_charge','NcircX','gauge_Higgs_norm_over_H','vector_mass_over_H']
    d={k:r[k] for k in ['H_GeV','cH','delta','a_end','rtol','atol','yield_factor','runtime_seconds']}
    d.update({k:w[k] for k in keep});d.update(r['controls'])
    d['NX_300_600']=earlier['NX'];d['epsX_300_600']=earlier['epsX']
    d['late_charge_drift_fraction']=abs(w['NX']-earlier['NX'])/max(abs(w['NX']),1e-30)
    d['EY_over_EX']=w['EY']/w['EX']
    d['absolute_Y_charge_mean_over_X']=abs(w['NY'])/max(abs(w['NX']),1e-30)
    d['Y_charge_rms_over_X']=w['Neff_std']/3/max(abs(w['NX']),1e-30)
    for j,z in enumerate(r['final_canonical_state']):d[f'final_state_{j}']=z
    return d


def main():
    a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=8);args=a.parse_args()
    build();CACHE.mkdir(parents=True,exist_ok=True)
    grids={str(c):np.geomspace(3e11/np.sqrt(c),4.7e13,13).tolist() for c in [1.,10.]}
    jobs=[]
    for c in [1.,10.]:
        for H in grids[str(c)]:
            for j in range(1,12):
                delta=j*np.pi/12
                for shift in [0.,1e-9]:jobs.append(dict(H=H,cH=c,delta=delta+shift))
    source=hashlib.sha256((ROOT/'code/release.c').read_bytes()).hexdigest()
    config=dict(grid=grids,phase_indices=list(range(1,12)),phase_unit='pi/12',shift=1e-9,
                n_runs=len(jobs),C_integrator_sha256=source,
                provisional_regular_filter='abs(delta epsX under 1e-9 phase shift)<1e-4; late NX drift<0.01; candidate still requires tighter tolerance and seed ensemble')
    write(ROOT/'inputs/SCAN_CONFIG.json',config)
    records=[];t0=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(worker,j):j for j in jobs}
        for f in as_completed(futures):
            r=f.result();records.append(r)
            if len(records)%22==0:
                print(f"SCAN {len(records)}/{len(jobs)} elapsed {time.perf_counter()-t0:.1f}s H={r['H_GeV']:.5g} cH={r['cH']:g}",flush=True)
    records.sort(key=lambda r:(r['cH'],r['H_GeV'],r['delta']))
    rows=[compact(r) for r in records]
    with (ROOT/'SCAN_RESULTS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    groups=[]
    for c in [1.,10.]:
        for H in grids[str(c)]:
            rr=[r for r in rows if r['cH']==c and r['H_GeV']==H]
            phase=[]
            for i in range(0,len(rr),2):
                x,y=rr[i:i+2];de=abs(x['epsX']-y['epsX'])
                phase.append(dict(delta=x['delta'],epsX=x['epsX'],YX_over_YQ=x['YX_over_YQ'],
                                  phase_shift_abs_delta_epsX=de,phase_shift_abs_delta_NX=abs(x['NX']-y['NX']),
                                  late_charge_drift_fraction=x['late_charge_drift_fraction'],
                                  provisional_regular=bool(de<1e-4 and x['late_charge_drift_fraction']<.01),
                                  numerical_tolerance_not_yet_certified=True))
            groups.append(dict(H_GeV=H,cH=c,initial_point=minimum(H,c),phases=phase,
                               max_abs_epsX=max(abs(p['epsX']) for p in phase)))
    write(ROOT/'GRID_RESULTS.json',dict(config=config,groups=groups,runtime_seconds=time.perf_counter()-t0))
    print('COMPLETE',len(records),'runs',time.perf_counter()-t0,flush=True)
    candidates=[p|{'H_GeV':g['H_GeV'],'cH':g['cH']} for g in groups for p in g['phases']
                if p['provisional_regular'] and abs(p['YX_over_YQ'])>=1 and abs(p['epsX'])>=.15]
    write(ROOT/'PROVISIONAL_CANDIDATES.json',dict(scope='Screening only; not a viability certificate.',rows=candidates))
    print('PROVISIONAL_CANDIDATES',len(candidates),flush=True)


if __name__=='__main__':main()
