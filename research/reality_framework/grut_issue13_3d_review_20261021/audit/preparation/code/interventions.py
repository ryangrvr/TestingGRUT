"""Causal interventions on Codex's sparse seed s101 with Codex's unchanged expanding.cpp binary.

Usage: python3 -I interventions.py <codex_raw_dir> <scratch_runs_dir> <expanding_binary> [names...]
Writes modified mode tables (inputs [P]) and runs the frozen solver at Codex's main settings
(N=192, dx=0.4, dt=0.06, T=48, kappa=2) with OMP_NUM_THREADS=2. Codex's files are read only.
Columns of a mode table: nx ny nz | cu su | cw sw | cS sS  (u = Re X, w = Im X, S = phi/v).
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
import time
from pathlib import Path
import numpy as np

RAW, RUNS, EXE = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ONLY = sys.argv[4:]
BASE = np.loadtxt(RAW / 'modes_L76.8_s101.txt')


def variant(name):
    t = BASE.copy()
    if name == 'noX':          # remove the X perturbation: no S-X cross-correlation
        t[:, 3:7] = 0
    elif name == 'flipX':      # flip the X perturbation: cross-correlation sign reversed, S unchanged
        t[:, 3:7] *= -1
    elif name == 'noW':        # keep only the radial (u) X perturbation
        t[:, 5:7] = 0
    elif name == 'noU':        # keep only the phase (w) X perturbation
        t[:, 3:5] = 0
    elif name == 'gapfill':    # add one box-scale S mode (n=(1,0,0), k=0.0818), power 1/64, no X partner
        t = np.vstack([t, [1, 0, 0, 0, 0, 0, 0, 0.125, 0.125]])
    else:
        raise ValueError(name)
    path = RUNS / 'modes' / f'modes_L76.8_s101_{name}.txt'
    np.savetxt(path, t, fmt='%.17g')
    return path


def main():
    names = ONLY or ['flipX', 'noX', 'gapfill', 'noU', 'noW']
    env = dict(os.environ, OMP_NUM_THREADS='2')
    log = []
    for name in names:
        mp = variant(name)
        out = RUNS / f'int_{name}'
        cmd = [str(EXE), '--N', '192', '--dx', '0.4', '--dt', '0.06', '--T', '48', '--taui', '40.0', '--kappa', '2.0',
               '--kind', 'charge', '--noise', '1.0', '--threads', '2', '--sample', '1.2', '--top-sample', '4.8',
               '--topology', '1', '--modes', str(mp), '--out', str(out)]
        t0 = time.monotonic()
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        log.append({'name': name, 'cmd': cmd, 'rc': r.returncode, 'stdout': r.stdout.strip(), 'stderr': r.stderr.strip(),
                    'wall_s': time.monotonic() - t0})
        print(json.dumps(log[-1]), flush=True)
        (RUNS / 'interventions_log.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    main()
