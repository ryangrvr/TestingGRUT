"""Own re-implementation of the 64-sparse-mode preparation recipe (stated in Codex's PROTOCOL / run.py
docstring): 64 distinct nonzero integer wavevectors k ~ rint(Normal(0, L/2pi)), |k|*2pi/L <= 2.5, sign
canonicalised by the first nonzero component; then a (64,6) standard-normal coefficient block, each
field's (cos,sin) pair normalised so sum(c^2)/2 = 1. Same seed => same file at fixed L for any N.
Usage: python3 -I make_modes.py L seed out.txt
"""
import hashlib
import sys

import numpy as np


def main():
    L, seed, out = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    rng = np.random.default_rng(seed)
    ks = []
    while len(ks) < 64:
        k = np.rint(rng.normal(0, L / (2 * np.pi), 3)).astype(int)
        if not k.any() or np.linalg.norm(k) * 2 * np.pi / L > 2.5:
            continue
        nz = np.flatnonzero(k)[0]
        if k[nz] < 0:
            k = -k
        t = tuple(int(x) for x in k)
        if t not in ks:
            ks.append(t)
    c = rng.normal(size=(64, 6))
    for f in range(3):
        c[:, 2 * f:2 * f + 2] /= np.sqrt(np.sum(c[:, 2 * f:2 * f + 2] ** 2) / 2)
    np.savetxt(out, np.column_stack((ks, c)), fmt='%.17g')
    print(out, hashlib.sha256(open(out, 'rb').read()).hexdigest())


if __name__ == '__main__':
    main()
