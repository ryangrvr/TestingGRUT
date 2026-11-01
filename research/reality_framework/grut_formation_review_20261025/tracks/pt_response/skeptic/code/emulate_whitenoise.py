"""SKEPTIC (pt_response_skeptic), INDEPENDENT: distribution of Codex's single white-noise shell-trace estimator of the
shell-1 conjugate/normal RMS ratio, emulated on the exact within-shell response matrix G(k', k) (pt_response exact1
projections, s0 = 20; outputs restricted to shell 1, so cross-shell leakage is NOT included -> a lower bound on the
scatter).  Real-pattern constraint src_{-k} = conj(src_k), Gaussian amplitudes as for a real white-noise field.
Usage: python3 -I emulate_whitenoise.py <PROJ_exact1 npz> <out json>
"""
import sys, json
import numpy as np
d = np.load(sys.argv[1]); proj = d['proj']; sh = d['shell']; allm = d['allmodes']; meta = [eval(m) for m in d['meta']]
N = 48
idx = {int(allm[p]): p for p in range(len(allm))}
def neg(ix):
    i, j, k = np.unravel_index(ix, (N, N, N)); return int(np.ravel_multi_index(((-i) % N, (-j) % N, (-k) % N), (N, N, N)))
kin = []; cols = []
for m in range(len(meta) // 2):
    ix = meta[2 * m][0]
    for sign in (1, -1):
        cols.append(np.stack([proj[:, 4 * m + b] + sign * 1j * proj[:, 4 * m + 2 + b] for b in range(2)], axis=2))  # [t,a,b,mode]
        kin.append(ix if sign > 0 else neg(ix))
G = np.stack(cols, -1)                       # [t, a, b, kout(all 388), kin(18)]
s1 = [idx[k] for k in kin]
G1 = G[:, :, :, s1, :]                       # outputs restricted to the 18 shell-1 modes, same order as kin
partner = [kin.index(neg(k)) for k in kin]
rng = np.random.default_rng(2026)
rat = []
for it in range(4000):
    z = (rng.normal(size=18) + 1j * rng.normal(size=18)) / np.sqrt(2)
    src = z.copy()
    for c, pc in enumerate(partner):
        if c < pc: src[pc] = np.conj(src[c])
    q = np.einsum('tabok,k->tabo', G1, src)
    x = np.einsum('tabo,o->tab', q, src.conj()) / np.sum(np.abs(src) ** 2)
    n = .5 * (x[:, 0, 0] + x[:, 1, 1] + 1j * (x[:, 1, 0] - x[:, 0, 1])); a = .5 * (x[:, 0, 0] - x[:, 1, 1] + 1j * (x[:, 1, 0] + x[:, 0, 1]))
    rat.append(np.linalg.norm(a) / np.linalg.norm(n))
rat = np.array(rat)
diag = np.stack([G1[:, :, :, c, c] for c in range(18)], -1)
n = .5 * (diag[:, 0, 0] + diag[:, 1, 1] + 1j * (diag[:, 1, 0] - diag[:, 0, 1])).mean(-1); a = .5 * (diag[:, 0, 0] - diag[:, 1, 1] + 1j * (diag[:, 1, 0] + diag[:, 0, 1])).mean(-1)
res = dict(exact_trace_ratio=float(np.linalg.norm(a) / np.linalg.norm(n)), emulated_percentiles_5_16_50_84_95=np.percentile(rat, [5, 16, 50, 84, 95]).tolist(),
           fraction_above_1p70=float(np.mean(rat > 1.7043)), fraction_below_1p028=float(np.mean(rat < 1.0278)), codex_single_probes=[1.0278, 1.7043])
print(json.dumps(res, indent=1)); json.dump(res, open(sys.argv[2], 'w'), indent=1)
