#!/usr/bin/env python3
"""Compute the separate LOWER_H_TEST.md scalar-lattice control."""
import argparse
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
TRACKS=HERE.parents[1]
TAG='KFmid_H55_actual_nl'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('work',type=Path)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--budget',type=float,default=550.)
    a=ap.parse_args();work=a.work.resolve()
    verdict=json.loads((work/'analysis/blind_verdict.json').read_text())
    assert all(verdict['nodes'][n]['grid_points']==135 for n in ('hLF84_mid','hLF_mid'))
    assert verdict['preregistration']=='6d836ca5051b7bf3e9d677e26c3b7b2db0ae0f6c'
    sys.path.insert(0,str(TRACKS/'lf_fermion'/'code'))
    import lfcommon as L
    import export as E
    import numpy as np
    L.SP=work/'lower_H';L.RUNS=L.SP/'runs';L.BIN=L.SP/'build'/'lf3d'
    prior=json.loads((L.RES/'hist_HFa_LF.json').read_text())['alternates']['KFmid_H84_nl']
    cfg=prior['config'];hist=prior['lattice'];rd=L.RUNS/TAG
    if not (rd/'run.par').exists():
        L.make_run(TAG,H_GeV=5.5e10,N=cfg['N'],L=cfg['L'],dt=cfg['dt'],seed=cfg['seed'],
                   kcut_mchi=cfg['kcut_mchi'],gw=0,snaps=[],extra=['tracers 0']+cfg['kappa'])
    info=json.loads((rd/'info.json').read_text());assert info['H_GeV']==5.5e10
    for k in ('N','L','dt','seed','kcut_mchi'):assert info[k]==cfg[k]
    if not L.run(TAG,budget=a.budget,threads=a.threads):
        print('Checkpointed scalar lattice; rerun this command to finish.',flush=True)
        return 75
    data=L.load_diag(TAG);head=E.head(data);block=E.lattice_block(data)
    assert abs(data['t'][-1]-3)<cfg['dt']
    assert all(np.all(np.isfinite(x)) for x in data.values())
    vals={k:{'H55_fresh_scalar_lattice':head[k],'H84_reused_lattice':prior['head'][k]}
          for k in ('f_rad','f_grad_max','t_fcoh_lt_0p1','budget_absmax','m_s2_mean_3','R_cyc_3')}
    vals.update({k+'_t3':{'H55_fresh_scalar_lattice':block[k][-1],'H84_reused_lattice':hist[k][-1]}
                 for k in ('a','H','f_SM','f_chi','f_sig')})
    out=work/'analysis';p=out/'lower_H_diag.txt.gz'
    p.write_bytes(gzip.compress((rd/'diag.txt').read_bytes(),mtime=0))
    summary={'scope':'H55 scalar lattice on transferred H84 fermion drain; independent of strict verdict',
             'run':info,'held_kappa':cfg['kappa'],'comparison':vals,'lattice_block':block,
             'diag_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
             'full_3D_Dirac_evolved':False,'new_kinetic_splice_computed':False,'used_in_registered_verdict':False}
    (out/'lower_H_lattice.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(vals,indent=2),flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
