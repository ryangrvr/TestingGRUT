#!/usr/bin/env python3
"""Fine single-partial-wave scan to look for narrow high-l resonances (Claude review)."""
import json
import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import scatter as S

ell = int(sys.argv[1]); lo, hi, n = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
rows = []
for p in np.linspace(lo, hi, n):
    Sm, chk = S.smatrix(ell, np.sqrt(1 + p * p), rtol=1e-9)
    rows.append([float(p), float(abs(Sm[0, 1]) ** 2), chk['unitarity']])
(S.SCR / f'ellscan_{ell}.json').write_text(json.dumps(rows))
P = np.array(rows)
i = np.argmax(P[:, 1]); print(ell, 'max P', P[i, 1], 'at', P[i, 0], 'max unitarity', P[:, 2].max())
