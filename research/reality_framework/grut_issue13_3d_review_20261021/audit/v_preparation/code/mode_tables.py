"""V1/V3: mode-table identity (large_s101 vs s101) and exact translation antisymmetry of the S mode set.

Usage: python3 -I mode_tables.py <codex_raw_dir> <out_json>
A translation d (lattice cells) with n_m . d = N/2 (mod N) for every populated S mode m maps the S-only field
to its negative: S(x+d) = -S(x). With homogeneous X the (Z2 x translation)-symmetric evolution then keeps
f+ = f- exactly. This checks whether the noX 'exact balance' is enforced by symmetry.
"""
import json, sys
from pathlib import Path
import numpy as np

raw, out = Path(sys.argv[1]), Path(sys.argv[2])
res = {}
a = np.loadtxt(raw / 'modes_L76.8_s101.txt'); b = np.loadtxt(raw / 'modes_L115.2_s101.txt')
res['coef_max_abs_diff_L76.8_vs_L115.2'] = float(np.max(np.abs(a[:, 3:] - b[:, 3:])))
ka = a[:, :3] * 2 * np.pi / 76.8; kb = b[:, :3] * 2 * np.pi / 115.2
dk = np.linalg.norm(ka - kb, axis=1)
res['phys_k_diff_max'] = float(dk.max()); res['phys_k_diff_median'] = float(np.median(dk))
res['kmin_76.8'] = float(np.linalg.norm(ka, axis=1).min()); res['kmin_115.2'] = float(np.linalg.norm(kb, axis=1).min())
for s in (101, 202, 303):
    t = np.loadtxt(raw / f'modes_L76.8_s{s}.txt'); n = t[:, :3].astype(np.int64)
    N = 192
    # search all translations d in Z_N^3 for n.d == N/2 mod N for all modes (chunked brute force)
    found = []
    for dx_ in range(N):
        dy, dz = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        ph = (n[:, 0, None, None] * dx_ + n[:, 1, None, None] * dy[None] + n[:, 2, None, None] * dz[None]) % N
        ok = np.all(ph == N // 2, axis=0)
        if ok.any():
            for i, j in zip(*np.nonzero(ok)):
                found.append([dx_, int(dy[i, j]), int(dz[i, j])])
                if len(found) >= 4: break
        if len(found) >= 4: break
    res[f's{s}_antisymmetric_translations_found'] = found
    res[f's{s}_parity_counts'] = {f'{p}': int(np.sum((n.sum(1) % 2) == p)) for p in (0, 1)}
    res[f's{s}_nx_odd_all'] = bool(np.all(n[:, 0] % 2 == 1))
print(json.dumps(res, indent=1)); out.write_text(json.dumps(res, indent=1))
