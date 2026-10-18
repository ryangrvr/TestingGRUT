#!/usr/bin/env python3
"""Reproduce in a new directory, preserving all archived execution records."""
from pathlib import Path
import argparse, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True,help='New, nonexistent directory for the finite reproduction run')
    parser.add_argument('--prepare-only',action='store_true',help='Create the clean layout without executing physics')
    args=parser.parse_args();out=args.output.resolve()
    if out.exists():raise SystemExit('Output directory already exists; choose a fresh path to preserve existing records.')
    followup=out/ROOT.name
    followup.mkdir(parents=True)
    for path in ROOT.glob('*.py'):
        shutil.copy2(path,followup/path.name)
    for name in ['PROTOCOL.md','PROTOCOL_AMENDMENT_01.md','independent_protocol.json','independent_refinement_protocol.json']:
        shutil.copy2(ROOT/name,followup/name)
    background=ROOT.parent/'grut_candidate_1_20261018/code'
    newbackground=out/'grut_candidate_1_20261018/code'
    newbackground.mkdir(parents=True)
    for name in ['results_profile.npz','results.json','reporter_results.json']:
        shutil.copy2(background/name,newbackground/name)
    print('Prepared '+str(followup),flush=True)
    if args.prepare_only:return
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    commands=[['branch_scan.py'],['branch_scan.py','--refine'],['independent_minimize.py'],
              ['independent_refine.py'],['independent_analyze.py'],['evolve.py'],['analyze.py'],['compare_slices.py'],['verify.py'],['plot.py']]
    for command in commands:
        print('Running '+' '.join(command),flush=True)
        with (followup/(command[0].removesuffix('.py')+'_reproduction.log')).open('a') as log:
            subprocess.run([sys.executable,str(followup/command[0]),*command[1:]],cwd=followup,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Completed finite reproduction in '+str(followup),flush=True)


if __name__=='__main__':main()
