"""Partial-wave scattering strengths at selected stationary charges.

The same radial discretization and partial-wave cutoffs as the parent are used.
Files are new artifacts; the parent kernels are never edited. Resolving resonances
in momentum and interpolating between charges still require independent controls.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
from mass_analysis import stationary, D, Seed
import core

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_kinetics_20261028'

def use_q(q):
    sol=D.solve(q,.002,Seed(),length=160.,tol=2e-9)
    core.OMEGA=float(sol.p[0])
    def prof(r):
        r=np.asarray(r);out=sol.sol(np.minimum(r,160))
        return np.where(r>160,0,out[0]),np.where(r>160,1,out[2])
    core.profile=prof
    return sol

def main():
    p=argparse.ArgumentParser();p.add_argument('--q',type=float,required=True);p.add_argument('--scan',action='store_true')
    args=p.parse_args();started=time.monotonic();sol=use_q(args.q)
    if args.scan:
        ks=np.linspace(.27,.43,321);px=[]
        for k in ks:
            P,_,_=core.matrix(k,4)
            px.append(P[0,1])
        peak=int(np.argmax(px))
        out={'q':args.q,'omega':float(sol.p[0]),'peak_k':float(ks[peak]),'peak_P':float(px[peak]),
             'seconds':time.monotonic()-started}
        (ROOT/f'SCAN_Q{int(args.q)}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out),flush=True);return
    # The ell=4 resonance moves with Q, so the reference-only fine window
    # is insufficient. Use a shared dense interval spanning observed peaks.
    with np.load(PARENT/'CORE_KERNELS.npz') as d:ks=np.unique(np.r_[d['k'],np.linspace(.31,.405,381)])
    rows=[];checks=[]
    for i,k in enumerate(ks):
        w,u,rec,l=core.strengths(float(k))
        rows.append(w);checks.append([u,rec,l])
        if (i+1)%100==0:print(int(args.q),i+1,'/',len(ks),'s',round(time.monotonic()-started,1),flush=True)
    output=ROOT/f'core_Q{int(args.q)}.npz'
    np.savez_compressed(output,k=ks,W=np.array(rows),checks=np.array(checks),Q=args.q,Omega=float(sol.p[0]))
    print('done',int(args.q),'seconds',round(time.monotonic()-started,1),'max unitarity',np.max(np.array(checks)[:,0]),flush=True)

if __name__=='__main__':main()
