"""Unit controls for topo.periodic_components (run: python3 -I test_topo.py)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
from topo import summarize
N = 24
x = np.arange(N)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
v = N ** 3
cases = {}
slab = (X >= 6) & (X < 12)                       # slab normal to x: winds y,z -> 6 ; complement also 6
cases['slab'] = (summarize(slab, v)['any_winding'], 6)
cases['slab_complement'] = (summarize(~slab, v)['any_winding'], 6)
d = np.minimum(X, N - X) ** 2 + np.minimum(Y, N - Y) ** 2 + np.minimum(Z, N - Z) ** 2
blob = d < 25                                     # corner blob crossing every cut, contractible -> 0
cases['corner_blob'] = (summarize(blob, v)['any_winding'], 0)
cases['corner_blob_complement'] = (summarize(~blob, v)['any_winding'], 7)
tube = ((Y - 5) ** 2 + (Z - 7) ** 2) < 9           # tube along x -> 1
cases['tube_x'] = (summarize(tube, v)['any_winding'], 1)
hel = np.zeros((N, N, N), bool)                   # diagonal staircase line along (1,1,0) -> bits x,y = 3
for i in range(N):
    hel[i, i, 3] = True; hel[i, (i + 1) % N, 3] = True
cases['diag_line'] = (summarize(hel, v)['any_winding'], 3)
ring = np.zeros((N, N, N), bool)                  # contractible ring in xy plane crossing x-boundary -> 0
for t in np.linspace(0, 2 * np.pi, 400):
    ring[int(round(5 * np.cos(t))) % N, int(round(10 + 5 * np.sin(t))) % N, 4] = True
cases['contractible_ring_across_boundary'] = (summarize(ring, v)['any_winding'], 0)
ok = all(a == b for a, b in cases.values())
print(cases, 'ALL_PASS' if ok else 'FAIL')
