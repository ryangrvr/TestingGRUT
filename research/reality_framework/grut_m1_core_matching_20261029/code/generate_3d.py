"""Recreate a specified prefix of the parent periodic 3D object evolution."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_dynamics_20261027'
sys.path.insert(0,str(PARENT/'code'))
os.environ['GRUT_M1_DYNAMICS_WORK']=str(ROOT/'.work'/'parent_common')
import common

def main():
    a=argparse.ArgumentParser();a.add_argument('--N',type=int,required=True)
    a.add_argument('--perturbation',type=float,required=True);a.add_argument('--tag',required=True)
    a.add_argument('--T',type=float,default=80.);a.add_argument('--threads',type=int,default=2)
    a.add_argument('--bits',type=int,choices=[32,64],default=32)
    args=a.parse_args();d=ROOT/'.work'/args.tag;d.mkdir(parents=True,exist_ok=True)
    q,p,meta=common.object_state(args.N,64.,args.perturbation)
    common.write_initial(d/'initial.bin',q,p,64.)
    meta['initial_sha256']=hashlib.sha256((d/'initial.bin').read_bytes()).hexdigest()
    meta['initial_totals']=common.totals(q,p,64.)
    del q,p
    exe=ROOT/'.work'/'field3d'
    if not exe.exists():raise RuntimeError('Compile field3d.cpp as .work/field3d before starting concurrent runs')
    cmd=[str(exe),str(d/'initial.bin'),str(d/'state'),'0.04',str(args.T),'2','0.52','4',str(args.threads),str(args.bits)]
    with (d/'run.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    out={'configuration':vars(args),'initial':meta,'command':cmd,'snapshot_count':len(list(d.glob('state_s*.bin'))),
         'scope':'Periodic 64-wide box; 5% supplied perturbation is not the 10% radial cosmological pulse.'}
    (d/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
