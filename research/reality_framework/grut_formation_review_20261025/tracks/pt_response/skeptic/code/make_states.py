"""SKEPTIC helper (Codex engine from our copy; run from $PKC/code): advance the regenerated base48 states with the
production rule (one step() call per Delta s = 0.04 sample, hmax 0.002) to new probe phases s = 20.12 and 20.56 (default) or, with argument "late", s = 20.92 from 20.68."""
import sys
sys.path.insert(0, '.')
from pt_common import *  # noqa
set_box(12.); lib = library()
import sys as _s
JOBS = ((20., 3, 20.12), (20.44, 3, 20.56)) if len(_s.argv) < 2 else ((20.68, 6, 20.92),)
for src, n, dst in JOBS:
    z = np.ascontiguousarray(np.load(WORK / 'base48' / f'state_s{src:g}.npy'))
    for i in range(n): step(lib, z, .04, .002, threads=2)
    np.save(WORK / 'base48' / f'state_s{dst:g}.npy', z)
    m = z.mean(axis=(1, 2, 3)); print(dst, 'mean field', m[:2], 'r', float(np.hypot(m[0], m[1])))
