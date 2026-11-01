"""SKEPTIC (pt_response_skeptic), INDEPENDENT re-analysis of pt_response's exact shell-1 projections (PROJ_exact1_s20.npz,
Codex tangent engine, pt_response probe patterns cos(k.x), sin(k.x)).  Our own reconstruction of the full within-shell
response matrix G_ab(k', k) from the cos/sin pairs (complex pattern e^{ik.x} = cos + i sin; e^{-ik.x} = cos - i sin),
then (i) the exact shell trace and the conjugate/normal RMS ratio over delays 0-6; (ii) a mode-resolved ratio
(RMS of |G_A(k,k)| vs |G_N(k,k)| per mode, then averaged), which is what a quasiparticle (diagonal) reading needs;
(iii) the forward vs scattered amplitude.
Usage: python3 -I exact_shell1.py <PROJ_exact1 npz> <out json>
"""
import sys, json
import numpy as np
d = np.load(sys.argv[1]); proj = d['proj']; srcs = d['srcs']; sh = d['shell']; allm = d['allmodes']; meta = [eval(m) for m in d['meta']]
N = 48
nt = proj.shape[0]
s1 = np.flatnonzero(sh == 1)           # positions (within allmodes) of shell-1 modes
idx = {int(allm[p]): p for p in range(len(allm))}
def neg(ix):
    i, j, k = np.unravel_index(ix, (N, N, N)); return int(np.ravel_multi_index(((-i) % N, (-j) % N, (-k) % N), (N, N, N)))
# patterns: meta[2m] = (ix,'cos'), meta[2m+1] = (ix,'sin'); tangent index = 2*pattern + b
pairs = [(meta[2 * m][0], 2 * m, 2 * m + 1) for m in range(len(meta) // 2)]
# G[t, a, b, kout(pos in allmodes), kin(pos in s1 list)]
cols = []; Gcols = []
for ix, pc, ps in pairs:
    for sign in (+1, -1):
        kin = ix if sign > 0 else neg(ix)
        # response to e^{sign i k.x}: cos + sign*i*sin ; tangent 2*p+b kicks pi_b
        resp = np.stack([proj[:, 2 * pc + b] + sign * 1j * proj[:, 2 * ps + b] for b in range(2)], axis=2)  # [t, a, b, mode]
        cols.append(idx[kin]); Gcols.append(resp)
G = np.stack(Gcols, axis=-1)  # [t, a, b, kout, kin]
cols = np.array(cols)
def chan(x):
    n = .5 * (x[:, 0, 0] + x[:, 1, 1] + 1j * (x[:, 1, 0] - x[:, 0, 1]))
    a = .5 * (x[:, 0, 0] - x[:, 1, 1] + 1j * (x[:, 1, 0] + x[:, 0, 1]))
    return n, a
diag = np.stack([G[:, :, :, cols[c], c] for c in range(len(cols))], axis=-1)  # [t,a,b,mode]
Nd, Ad = chan(diag)          # [t, mode]
trN, trA = Nd.mean(-1), Ad.mean(-1)
ratio_trace = float(np.linalg.norm(trA) / np.linalg.norm(trN))
per_mode = np.linalg.norm(Ad, axis=0) / np.linalg.norm(Nd, axis=0)
rms_mode = float(np.sqrt(np.sum(np.abs(Ad) ** 2) / np.sum(np.abs(Nd) ** 2)))
# initial slope check (normal ~ t)
times = np.arange(nt) * .04
slope = float((trN[1] / times[1]).real)
# scattered: all outputs (shells 1-4) except the kicked mode, both channels
Nall, Aall = chan(G)        # [t, kout, kin]
fwd = np.sqrt(np.sum(np.abs(Nd) ** 2 + np.abs(Ad) ** 2, axis=-1))
tot = np.sqrt(np.sum(np.abs(Nall) ** 2 + np.abs(Aall) ** 2, axis=(1, 2)))
sca = np.sqrt(np.maximum(tot ** 2 - fwd ** 2, 0))
bins = [(0, 1), (1, 2), (2, 4), (4, 6.01)]
scat = {f'{a}-{b}': float(np.sqrt(np.mean((sca ** 2)[(times >= a) & (times < b)]) / np.mean((fwd ** 2)[(times >= a) & (times < b)]))) for a, b in bins}
# the 'anomalous' channel partner: for a homogeneous condensate G_A couples k -> k (as defined); check if the anomalous
# response at the kicked mode k is larger at k or at -k
partner = np.stack([G[:, :, :, idx[neg(int(allm[cols[c]]))], c] for c in range(len(cols))], axis=-1)
Np, Ap = chan(partner)
res = dict(s0=float(d['s0']), ratio_exact_trace=ratio_trace, ratio_mode_rms=rms_mode,
           per_mode_ratio_min_max=[float(per_mode.min()), float(per_mode.max())], per_mode_ratio_median=float(np.median(per_mode)),
           initial_slope_normal=slope, scattered_over_forward_rms=scat,
           anomalous_at_minus_k_over_at_k=float(np.linalg.norm(Ap) / np.linalg.norm(Ad)),
           trace_by_window={f'{a}-{b}': float(np.linalg.norm(trA[(times >= a) & (times < b)]) / np.linalg.norm(trN[(times >= a) & (times < b)])) for a, b in bins})
print(json.dumps(res, indent=1))
json.dump(res, open(sys.argv[2], 'w'), indent=1)
