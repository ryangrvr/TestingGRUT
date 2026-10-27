"""en_neutral item 3 driver: REPRODUCTION-type run of Codex's early3d.cpp (compiled from a git-archive COPY of the
7d9c93b packet, unmodified) on an initial file written by our own generator (en_lib; identical to Codex's
initial.py realization to 5e-16 in energy).  Only the snapshot cadence differs from Codex's production (denser,
32-bit), so the phi spectrum can be followed in time.  The dense trace (Delta tau = 0.05) is the engine's own.

usage: python3 -I run_codex_engine.py <published_dir> <early3d.cpp> <workdir> <N> <L> <snap_tau> <end_tau> [threads] [seed] [cutoff]
"""
import sys, json, struct, subprocess, hashlib, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import en_lib as E

pub, src, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
N, L, snap, end_tau = int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
threads = int(sys.argv[8]) if len(sys.argv) > 8 else 2
seed = int(sys.argv[9]) if len(sys.argv) > 9 else 1
cutoff = float(sys.argv[10]) if len(sys.argv) > 10 else 8.0
work.mkdir(parents=True, exist_ok=True)
M = np.load(pub / f'MATCH_L{L:g}.npz')
a0, C = float(M['a_match']), float(M['C_match'])
n, coef = E.realization_coefficients(M, L, seed=seed, cutoff=cutoff)
q, p = E.fields_from_coefficients(n, coef, M['light_mean_q'], M['light_mean_p'], N)
rho0 = E.energy_density(q, p, L)
init = work / 'initial.bin'
with open(init, 'wb') as f:
    f.write(b'GRUT3D1\n'); f.write(struct.pack('<id4d', N, L, 0.0, E.G2, E.LX, E.LP))
    np.asarray(q, '<f8').tofile(f); np.asarray(p, '<f8').tofile(f)
key = hashlib.sha256(src.read_bytes()).hexdigest()[:16]
exe = work.parent / ('early3d_' + key)
if not exe.exists():
    subprocess.run(['g++', '-std=c++17', '-O3', '-ffp-contract=off', '-fopenmp', str(src), '-o', str(exe)], check=True)
duration = 2 * (end_tau - 20) / (a0 + np.sqrt(a0 * a0 + 2 * C * (end_tau - 20)))
for old in work.glob('state_s*.bin'): old.unlink()
args = [str(exe), str(init), str(work / 'state'), '0.004', repr(float(duration)), '0.05', repr(float(snap)), '4', str(threads), '32', '20', repr(float(a0)), repr(float(C))]
t0 = time.time()
with open(work / 'run.log', 'w') as f: subprocess.run(args, check=True, stdout=f, stderr=subprocess.STDOUT)
tr = np.genfromtxt(work / 'state.tsv', names=True, delimiter='\t')
info = dict(N=N, L=L, seed=seed, cutoff=cutoff, snap_tau=snap, end_tau=end_tau, rho0_own=rho0, a0=a0, C=C, exe_source_sha=key,
            seconds=time.time() - t0, energy_drift=float(np.max(abs(tr['rho'] / tr['rho'][0] - 1))),
            charge_drift=float(np.max(abs(tr['nQ'] / tr['nQ'][0] - 1))), final_tau=float(tr['tau'][-1]),
            final_phi_plus_int=float((tr['Kphi'][-1] + tr['Gphi'][-1] + tr['Uphi'][-1] + tr['Umix'][-1]) / tr['rho'][-1]))
(work / 'info.json').write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
