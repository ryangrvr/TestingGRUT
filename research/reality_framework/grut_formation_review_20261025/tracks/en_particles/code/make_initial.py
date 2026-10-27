"""en_particles helper (reproduction-type): build a seeded initial state with Codex's UNMODIFIED initial.py and save it
in the engine's snapshot format so window_run.py can evolve it with the conformal-time clock.
Usage: python3 -B make_initial.py <packet_code_dir> <out.bin> N L seed [cutoff=8]
"""
import sys, struct
from pathlib import Path
import numpy as np
sys.path.insert(0, sys.argv[1])
from initial import initialize
N, L, seed = int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]); cut = float(sys.argv[6]) if len(sys.argv) > 6 else 8.
q, p, prm = initialize(N, L, seed, cut)
with open(sys.argv[2], 'wb') as f:
    f.write(struct.pack('<iddi', N, L, 0., 64)); np.asarray(q, '<f8').tofile(f); np.asarray(p, '<f8').tofile(f)
print(prm['nQ'], prm['rho'], prm['a_match'], prm['C_match'])
