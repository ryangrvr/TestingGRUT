"""SKEPTIC (pt_response_skeptic), INDEPENDENT: build measured low-noise shell traces (normal, conjugate) at five orbital
phases from the probe projections (Codex tangent engine; probe designs: pt_response exact1/hutch, ours skp A/B).
s0 = 20     : shell 1 from the exact per-mode probes (36 tangents), shells 2-4 from pt_response Hutchinson M=8
s0 = 20.44  : pt_response Hutchinson M=8 (joint shells 1-4)
s0 = 20.12, 20.24, 20.68 : ours, family A (shell 1 only, M=8) and family B (shells 2-4, M=8)
Output npz: tr_normal[phase, shell(1..4), t], tr_anom[...], noise estimates (std of the probe mean) and the ratios.
Usage: python3 -I measured_traces.py <pt_response probes dir> <skeptic probes dir> <out npz> <out json>
"""
import sys, json
import numpy as np
A, B, OUTN, OUTJ = sys.argv[1:5]
N = 48

def chan(x):
    n = .5 * (x[..., 0, 0] + x[..., 1, 1] + 1j * (x[..., 1, 0] - x[..., 0, 1]))
    a = .5 * (x[..., 0, 0] - x[..., 1, 1] + 1j * (x[..., 1, 0] + x[..., 0, 1]))
    return n, a

def probe_traces(f, fam=None):
    d = np.load(f); proj = d['proj']; srcs = d['srcs']; sh = d['shell']
    fams = [str(x) for x in d['fams']] if 'fams' in d.files else ['J'] * srcs.shape[0]
    out = {}
    for n in (1, 2, 3, 4):
        sel = sh == n
        js = [j for j in range(srcs.shape[0]) if np.any(np.abs(srcs[j, sel]) > 1e-12)]
        x = np.zeros((len(js), proj.shape[0], 2, 2), complex)
        for c, j in enumerate(js):
            for b in range(2):
                q = proj[:, 2 * j + b][:, :, sel]
                x[c, :, :, b] = np.mean(q * srcs[j, sel].conj()[None, None, :], axis=-1)
        Nj, Aj = chan(x)
        out[n] = (Nj.mean(0), Aj.mean(0), Nj.std(0, ddof=1) / np.sqrt(len(js)), Aj.std(0, ddof=1) / np.sqrt(len(js)), len(js))
    return out

def exact1(f):
    d = np.load(f); proj = d['proj']; sh = d['shell']; allm = d['allmodes']; meta = [eval(m) for m in d['meta']]
    idx = {int(allm[p]): p for p in range(len(allm))}
    def neg(ix):
        i, j, k = np.unravel_index(ix, (N, N, N)); return int(np.ravel_multi_index(((-i) % N, (-j) % N, (-k) % N), (N, N, N)))
    diag = []
    for m in range(len(meta) // 2):
        ix = meta[2 * m][0]; pc, ps = 2 * m, 2 * m + 1
        for sign in (1, -1):
            kin = ix if sign > 0 else neg(ix)
            resp = np.stack([proj[:, 2 * pc + b] + sign * 1j * proj[:, 2 * ps + b] for b in range(2)], axis=2)
            diag.append(resp[:, :, :, idx[kin]])
    diag = np.stack(diag, -1)  # [t,a,b,mode]
    n, a = chan(np.moveaxis(diag, -1, 1))  # [t, mode]
    return n.mean(1), a.mean(1)

phases = [20., 20.12, 20.24, 20.44, 20.68]
nt = 151
TN = np.zeros((5, 4, nt), complex); TA = np.zeros_like(TN); EN = np.zeros((5, 4, nt)); EA = np.zeros_like(EN); src = {}
h20 = probe_traces(A + '/PROJ_hutch_s20_M8.npz'); e20 = exact1(A + '/PROJ_exact1_s20.npz')
h44 = probe_traces(A + '/PROJ_hutch_s20.44_M8.npz')
sk = {s: probe_traces(B + f'/PROJ_skp_s{s:g}_M8_seed424242.npz') for s in (20.12, 20.24, 20.68)}
for i, s in enumerate(phases):
    for n in (1, 2, 3, 4):
        if s == 20. and n == 1:
            TN[i, 0], TA[i, 0] = e20; src[f'{s}_{n}'] = 'exact'
        else:
            tab = h20 if s == 20. else h44 if s == 20.44 else sk[s]
            TN[i, n - 1], TA[i, n - 1], EN[i, n - 1], EA[i, n - 1], M = tab[n]; src[f'{s}_{n}'] = f'probe M={M}'
ratio = np.linalg.norm(TA, axis=-1) / np.linalg.norm(TN, axis=-1)
np.savez(OUTN, phases=phases, tr_normal=TN, tr_anom=TA, err_normal=EN, err_anom=EA, delays=np.arange(nt) * .04)
res = dict(phases=phases, source=src, ratio_rms_0_6=ratio.tolist(),
           ratio_range_shells_1_3=[float(ratio[:, :3].min()), float(ratio[:, :3].max())],
           ratio_range_shell_1=[float(ratio[:, 0].min()), float(ratio[:, 0].max())])
t = np.arange(nt) * .04
for a_, b_ in ((0, 2), (2, 6.01)):
    w = (t >= a_) & (t < b_)
    res[f'ratio_rms_{a_}_{b_:g}'] = (np.linalg.norm(TA[..., w], axis=-1) / np.linalg.norm(TN[..., w], axis=-1)).tolist()
print(json.dumps(res, indent=1))
json.dump(res, open(OUTJ, 'w'), indent=1)
