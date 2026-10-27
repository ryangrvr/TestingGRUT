"""en_neutral_skeptic item 3b: shell-averaged classical phi occupation at late times vs the half quantum.
For each lattice-k shell (width 0.25) of the N48 reproduction snapshots: nbar+1/2 = <e_k V / w_k>_shell, with
w_k^2 = k_lat^2 + g2<|chi_X|^2> + 3 lamP <phi^2>. Above-vacuum energy = sum_shells max(nbar, 0) * sum_shell w_k / V.
OWN code. usage: python3 -I uv_occupation.py <rundir> <out_json>"""
import sys, json, struct
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_lib as K
rd, out = Path(sys.argv[1]), Path(sys.argv[2])
tr = K.read_trace(rd / 'state.tsv'); rho0 = float(tr['rho'][0])
a0 = float(tr['a'][0]); C = float((tr['a'][-1] - tr['a'][0]) / tr['eta'][-1])
N, L = 48, 12.0; V = L ** 3; K2 = K.lattice_k2_grid(N, L); kl = np.sqrt(K2)
edges = np.arange(0, kl.max() + 0.25, 0.25); ids = np.clip(np.digitize(kl.ravel(), edges) - 1, 0, len(edges) - 2)
res = {}
for path in sorted(rd.glob('state_s*.bin'))[-6:]:
    with open(path, 'rb') as f:
        n_, L_, s_, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f4').reshape(6, n_, n_, n_).astype(np.float64)
    tau = 20 + a0 * s_ + 0.5 * C * s_ * s_
    X2 = float(np.mean(a[0] ** 2 + a[1] ** 2)); P2 = float(np.mean(a[2] ** 2))
    ph = np.fft.fftn(a[2]) / N ** 3; pp = np.fft.fftn(a[5]) / N ** 3
    w = np.sqrt(K2 + K.G2 * X2 + 3 * K.LP * P2)
    e = 0.5 * (abs(pp) ** 2 + w ** 2 * abs(ph) ** 2); e[0, 0, 0] = 0
    occ = (e * V / w).ravel(); ww = w.ravel(); ee = e.ravel()
    nb = np.bincount(ids, weights=occ, minlength=len(edges) - 1) / np.maximum(np.bincount(ids, minlength=len(edges) - 1), 1)
    sw = np.bincount(ids, weights=ww, minlength=len(edges) - 1)
    se = np.bincount(ids, weights=ee, minlength=len(edges) - 1)
    above = np.maximum(nb - 0.5, 0) * sw / V
    sel8 = edges[:-1] >= 8; seeded = (edges[:-1] >= 1.75) & (edges[:-1] < 8)
    res[f'{tau:.0f}'] = dict(total_quadratic=float(e.sum() / rho0),
                             above_half_quantum_shell_clipped=float(above.sum() / rho0),
                             above_half_quantum_seeded_1p8_8=float(above[seeded].sum() / rho0),
                             above_half_quantum_k_gt_8=float(above[sel8].sum() / rho0),
                             shell_nbar_plus_half_max_k_gt_8=float(nb[sel8 & (sw > 0)].max()),
                             shells_k_gt_8_with_nbar_plus_half_above_half=int((nb[sel8 & (sw > 0)] > 0.5).sum()),
                             shells_k_gt_8=int((sel8 & (sw > 0)).sum()),
                             nbar_plus_half_by_shell={f'{edges[i]:.2f}': round(float(nb[i]), 3) for i in range(0, len(nb), 4) if sw[i] > 0})
    print(f'{tau:.0f}', {k: v for k, v in res[f"{tau:.0f}"].items() if k != 'nbar_plus_half_by_shell'})
print(res[f'{tau:.0f}']['nbar_plus_half_by_shell'])
json.dump(dict(label='OWN shell-averaged classical phi occupation, N48 reproduction snapshots', snapshots=res), open(out, 'w'), indent=1)
