#!/usr/bin/env python3
"""Prepare and run one kg_radial case (spectrum track).  Usage:
  drive.py NAME --h H --dt DT --T T --L L [--rext ...] [--snaps ...]
Writes RUNDIR/NAME_* (binary outputs) and RUNDIR/NAME_meta.json.  Checks the no-reflection condition
2L - max(Rext) > T (a signal leaving r<0 at t=0 cannot return to any extraction radius inside the run)."""
import argparse, json, os, subprocess, sys, time, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('name'); ap.add_argument('--h', type=float, required=True); ap.add_argument('--dt', type=float, required=True)
ap.add_argument('--T', type=float, required=True); ap.add_argument('--L', type=float, required=True)
ap.add_argument('--rext', default='5,30,60,120,180,240,360,480,720,960')
ap.add_argument('--rin', default='30,45,60,90,120,180,240,360,480,720,960')
ap.add_argument('--snaps', default='3000,6000,12000,18000,24000')
ap.add_argument('--sdt', type=float, default=0.25); ap.add_argument('--hdt', type=float, default=1.0)
ap.add_argument('--rundir', required=True); ap.add_argument('--binary', required=True)
a = ap.parse_args()
rmax = max(float(x) for x in a.rext.split(','))
assert 2 * a.L - rmax > a.T, 'wall reflection could reach an extraction radius'
os.makedirs(a.rundir, exist_ok=True)
pre = os.path.join(a.rundir, a.name)
info = json.loads(subprocess.check_output([sys.executable, os.path.join(HERE, 'prep.py'), 'prepared', '--L', str(a.L), '--h', str(a.h), '--out', pre + '_init.bin']).decode().splitlines()[-1])
snaps = ','.join(s for s in a.snaps.split(',') if float(s) <= a.T) or 'none'
cmd = [a.binary, pre + '_init.bin', pre, str(a.dt), str(a.T), '1', str(a.sdt), a.rext, snaps, str(a.hdt), '900', a.rin]
t0 = time.time()
with open(pre + '.log', 'w') as lg:
    subprocess.run(cmd, stderr=lg, check=True)
meta = dict(name=a.name, h=a.h, dt=a.dt, T=a.T, L=a.L, rext=[float(x) for x in a.rext.split(',')], rin=[float(x) for x in a.rin.split(',')],
            snaps=snaps, sample_dt=a.sdt, hist_dt=a.hdt, prep=info, command=cmd, seconds=time.time() - t0,
            solver_sha256=hashlib.sha256(open(os.path.join(HERE, 'kg_radial.c'), 'rb').read()).hexdigest())
json.dump(meta, open(pre + '_meta.json', 'w'), indent=1)
print(json.dumps(meta))
