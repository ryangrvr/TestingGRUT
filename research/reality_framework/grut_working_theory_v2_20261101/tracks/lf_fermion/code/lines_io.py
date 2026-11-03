"""Reader of lf3d line dumps (lines.bin, trlines.bin). Claude, lf_fermion (issue #13 context)."""
import numpy as np
from pathlib import Path


def read_lines(rd):
    rd = Path(rd)
    with open(rd / 'lines.bin', 'rb') as f:
        N, nl = np.fromfile(f, dtype=np.int32, count=2)
        ax = np.fromfile(f, dtype=np.int32, count=nl)
        i1 = np.fromfile(f, dtype=np.int32, count=nl)
        i2 = np.fromfile(f, dtype=np.int32, count=nl)
        raw = np.fromfile(f, dtype=np.float64)
    rec = 3 + nl * N * 4
    nrec = raw.size // rec
    raw = raw[:nrec * rec].reshape(nrec, rec)
    tah = raw[:, :3]
    F = raw[:, 3:].reshape(nrec, nl, N, 4)      # sigma, sigma_dot, chi, chi_dot (v units, cosmic-time derivative)
    t = tah[:, 0]
    _, idx = np.unique(t, return_index=True)
    idx = np.sort(idx)
    return dict(N=int(N), nl=int(nl), ax=ax, i1=i1, i2=i2, t=tah[idx, 0], a=tah[idx, 1], H=tah[idx, 2], F=F[idx])


def read_trlines(rd, N, nl, NL=5, NCH=4):
    raw = np.fromfile(Path(rd) / 'trlines.bin', dtype=np.float64)
    rec = 3 + NL * NCH * nl * N
    nrec = raw.size // rec
    raw = raw[:nrec * rec].reshape(nrec, rec)
    t = raw[:, 0]
    _, idx = np.unique(t, return_index=True)
    idx = np.sort(idx)
    return dict(t=raw[idx, 0], a=raw[idx, 1], TE=raw[idx, 3:].reshape(len(idx), NL, NCH, nl, N))
