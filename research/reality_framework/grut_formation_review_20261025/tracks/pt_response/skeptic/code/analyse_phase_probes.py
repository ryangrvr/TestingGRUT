"""SKEPTIC (pt_response_skeptic), INDEPENDENT analysis of probe projections.

Reads our PROJ_skp_*.npz (probes_phase.py; families A = shell-1-only, B = shells 2-4) and, for a like-for-like
re-analysis, pt_response's PROJ_hutch_*.npz (joint shells 1-4) and PROJ_exact1_s20.npz.

For every probe j, shell n and kick component b we form the per-probe shell trace
    x_j[t, a, b] = mean_{k in shell} q^{(j,b)}_{a,k}(t) conj(src_{j,k})        (|src_k| = 1)
and the two channels  N = .5(x00+x11+i(x10-x01)),  A = .5(x00-x11+i(x10+x01)).
Estimators of the conjugate/normal RMS ratio over delays 0-6 (151 samples):
  'mean'     : ||mean_j A_j|| / ||mean_j N_j||          (what pt_response quotes; noise biases both norms upward)
  'unbiased' : sqrt( U(A) / U(N) ),  U(X) = sum_{i != j} Re<X_i, X_j> / (M(M-1))   (noise-free estimate of ||E X||^2)
  bootstrap over probes (2000 resamples) for both.
Usage: python3 -I analyse_phase_probes.py <out json> <npz> [<npz> ...]
"""
import sys, json, os
import numpy as np

OUT = sys.argv[1]; files = sys.argv[2:]

def chan(x):  # x[..., a, b]
    n = .5 * (x[..., 0, 0] + x[..., 1, 1] + 1j * (x[..., 1, 0] - x[..., 0, 1]))
    a = .5 * (x[..., 0, 0] - x[..., 1, 1] + 1j * (x[..., 1, 0] + x[..., 0, 1]))
    return n, a

def U(X):  # X[j, t] complex; unbiased ||E X||^2 from cross products
    M = X.shape[0]; S = X.sum(0)
    return float((np.vdot(S, S).real - np.sum(np.abs(X) ** 2)) / (M * (M - 1)))

def ratios(Nj, Aj, rng, nb=2000):
    M = Nj.shape[0]
    r_mean = np.linalg.norm(Aj.mean(0)) / np.linalg.norm(Nj.mean(0))
    un, ua = U(Nj), U(Aj)
    r_unb = np.sqrt(ua / un) if un > 0 and ua > 0 else float('nan')
    bm, bu = [], []
    for _ in range(nb):
        ix = rng.integers(0, M, M)
        bm.append(np.linalg.norm(Aj[ix].mean(0)) / np.linalg.norm(Nj[ix].mean(0)))
        # unbiased estimator on a resample: use distinct-index cross products only
        S = Nj[ix]; T = Aj[ix]
        num = np.vdot(T.sum(0), T.sum(0)).real - np.sum(np.abs(T) ** 2); den = np.vdot(S.sum(0), S.sum(0)).real - np.sum(np.abs(S) ** 2)
        bu.append(np.sqrt(num / den) if num > 0 and den > 0 else np.nan)
    bm = np.array(bm); bu = np.array(bu)
    single = np.linalg.norm(Aj, axis=1) / np.linalg.norm(Nj, axis=1)
    return dict(M=M, ratio_mean_trace=float(r_mean), ratio_unbiased=float(r_unb),
                boot_mean_trace_16_50_84=[float(x) for x in np.nanpercentile(bm, [16, 50, 84])],
                boot_unbiased_16_50_84=[float(x) for x in np.nanpercentile(bu, [16, 50, 84])],
                single_probe_ratio_min_max=[float(single.min()), float(single.max())],
                normal_noise_per_probe=float(np.mean(np.linalg.norm(Nj - Nj.mean(0), axis=1)) / np.linalg.norm(Nj.mean(0))))

res = {}
rng = np.random.default_rng(12345)
for f in files:
    d = np.load(f, allow_pickle=False); tag = os.path.basename(f)[5:-4]
    proj = d['proj']; srcs = d['srcs']; sh = d['shell']; s0 = float(d['s0'])
    if 'exact1' in tag:
        continue
    nprobe = srcs.shape[0]
    if 'fams' in d.files:
        fams = [str(x) for x in d['fams']]
    else:
        fams = ['J'] * nprobe  # pt_response joint shells 1-4
    rows = {}
    for fam in sorted(set(fams)):
        js = [j for j in range(nprobe) if fams[j] == fam]
        for n in (1, 2, 3, 4):
            sel = sh == n
            if np.all(np.abs(srcs[js][:, sel]) < 1e-12):
                continue
            x = np.zeros((len(js), proj.shape[0], 2, 2), complex)
            for c, j in enumerate(js):
                for b in range(2):
                    q = proj[:, 2 * j + b][:, :, sel]  # [t, a, mode]
                    x[c, :, :, b] = np.mean(q * srcs[j, sel].conj()[None, None, :], axis=-1)
            Nj, Aj = chan(x)
            rows[f'{fam}_shell{n}'] = ratios(Nj, Aj, rng)
    res[tag] = dict(s0=s0, file=f, rows=rows)
    print(tag, s0)
    for k, v in rows.items():
        print('   %-10s mean %.3f  unbiased %.3f  boot(unb) %s  single %s  noise %.2f' % (
            k, v['ratio_mean_trace'], v['ratio_unbiased'], np.round(v['boot_unbiased_16_50_84'], 3),
            np.round(v['single_probe_ratio_min_max'], 3), v['normal_noise_per_probe']))
json.dump(res, open(OUT, 'w'), indent=1)
