"""Build s101 variant mode tables independently and run Codex's unchanged expanding binary (192^3, main settings).
Usage: python3 -I run_variants.py <codex_raw_dir> <runs_dir> <expanding_binary> name [name...]
names: flipX (negate X columns), Xscale<s> (scale X columns by s). S columns unchanged."""
import json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
raw, runs, exe = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
runs.mkdir(parents=True, exist_ok=True)
base = np.loadtxt(raw / 'modes_L76.8_s101.txt')
for name in sys.argv[4:]:
    t = base.copy()
    if name == 'flipX':
        t[:, 3:7] = -t[:, 3:7]
    elif name.startswith('Xscale'):
        t[:, 3:7] = t[:, 3:7] * float(name[6:])
    else:
        raise SystemExit(name)
    mp = runs / f'modes_v_{name}.txt'; np.savetxt(mp, t, fmt='%.17g')
    cmd = [str(exe), '--N', '192', '--dx', '0.4', '--dt', '0.06', '--T', '48', '--taui', '40.0', '--kappa', '2.0', '--kind', 'charge',
           '--noise', '1.0', '--threads', '2', '--sample', '1.2', '--top-sample', '4.8', '--topology', '1', '--modes', str(mp),
           '--out', str(runs / f'v_{name}')]
    t0 = time.monotonic(); r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='2'))
    print(json.dumps({'name': name, 'rc': r.returncode, 'stdout': r.stdout.strip(), 'stderr': r.stderr.strip(), 'wall_s': time.monotonic() - t0}), flush=True)
