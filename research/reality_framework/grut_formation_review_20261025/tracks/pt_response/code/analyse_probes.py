"""INDEPENDENT analysis (Claude, pt_response) of our own probe projections (probe_study.py output).

exact1: builds G_ab(k', k) for all 18 shell-1 modes (k, k' in shell 1), the exact shell trace
        T_ab = (1/18) sum_k G_ab(k, k), the within-shell off-diagonal (scattering) weight, and the leakage of
        shell-1 kicks into shells 2-4.  Also simulates Codex's white-noise estimator restricted to shell-1 sources
        to show its probe-to-probe scatter.
hutch : Hutchinson trace per shell 1-4 from M unit-modulus random-phase probes; mean and standard error.
Writes MEAS_<tag>.npz (normal/anomalous [151, 4 shells], NaN where not measured) for hfb_driven.py, and a JSON.
Usage: python3 -I analyse_probes.py <probes dir> <out json>
"""
import sys, glob, json, os
import numpy as np

D = sys.argv[1].rstrip('/') + '/'; OUT = sys.argv[2]
N, L = 48, 12.
mm = np.fft.fftfreq(N) * N

def chan(T):  # T[t, a, b]
    nr = .5 * (T[:, 0, 0] + T[:, 1, 1] + 1j * (T[:, 1, 0] - T[:, 0, 1]))
    an = .5 * (T[:, 0, 0] - T[:, 1, 1] + 1j * (T[:, 1, 0] + T[:, 0, 1]))
    return nr, an

res = {}
for f in sorted(glob.glob(D + 'PROJ_exact1_*.npz')):
    d = np.load(f); proj = d['proj']; allm = d['allmodes']; sh = d['shell']; meta = [eval(m) for m in d['meta']]
    s0 = float(d['s0']); tag = os.path.basename(f)[5:-4]
    m1 = allm[sh == 1]; pos = {int(ix): j for j, ix in enumerate(allm)}
    def neg(ix):
        i, j, k = np.unravel_index(ix, (N, N, N)); return int(np.ravel_multi_index(((-i) % N, (-j) % N, (-k) % N), (N, N, N)))
    # complex-pattern responses: resp[kick mode kk (both signs)][t, a, b, mode]
    kicks = {}
    pats = {}
    for p, (ix, kind) in enumerate(meta):
        pats.setdefault(ix, {})[kind] = p
    for ix, pp in pats.items():
        rc = proj[:, [2 * pp['cos'], 2 * pp['cos'] + 1]]  # [t, b, a, mode]
        rs = proj[:, [2 * pp['sin'], 2 * pp['sin'] + 1]]
        kicks[ix] = np.transpose(rc + 1j * rs, (0, 2, 1, 3))        # response to e^{+ik.x}: [t, a, b, mode]
        kicks[neg(ix)] = np.transpose(rc - 1j * rs, (0, 2, 1, 3))   # response to e^{-ik.x}
    order = [int(ix) for ix in m1]
    G = np.zeros((151, 2, 2, len(order), len(order)), complex)  # G[t,a,b,k',k]
    for c, k in enumerate(order):
        G[:, :, :, :, c] = kicks[k][:, :, :, [pos[kk] for kk in order]]
    diag = np.einsum('tabkk->tabk', G)
    T = diag.mean(-1)
    nr, an = chan(T)
    offw = np.sqrt(np.sum(abs(G) ** 2) - np.sum(abs(diag) ** 2)) / np.sqrt(np.sum(abs(diag) ** 2))
    # per-mode spread of diagonal normal response
    nrm, anm = chan(np.moveaxis(diag, -1, 0).reshape(-1, 2, 2)); nrm = nrm.reshape(len(order), 151); anm = anm.reshape(len(order), 151)
    spread = float(np.mean([np.linalg.norm(nrm[c] - nr) for c in range(len(order))]) / np.linalg.norm(nr))
    # leakage into other shells
    leak = {}
    for s_ in (2, 3, 4):
        ms = [pos[int(ix)] for ix in allm[sh == s_]]
        tot = sum(np.sum(abs(kicks[k][:, :, :, ms]) ** 2) for k in order)
        leak[s_] = float(np.sqrt(tot / np.sum(abs(diag) ** 2)))
    # Codex-style white-noise estimator restricted to shell-1 sources: R = sum_{k',k} G(k',k) src_k conj(src_k') / sum|src|^2
    # time-resolved scattering: all response amplitude outside the kicked mode (shells 1-4 projected) vs diagonal
    tot_t = np.zeros(151); dia_t = np.zeros(151)
    for c, k in enumerate(order):
        allv = kicks[k]  # [t, a, b, mode]
        tot_t += np.sum(abs(allv) ** 2, axis=(1, 2, 3)); dia_t += np.sum(abs(allv[:, :, :, pos[k]]) ** 2, axis=(1, 2))
    blocks = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6.01)]
    scat_t = {f'{a}-{b:g}': float(np.sqrt(np.sum((tot_t - dia_t)[(np.arange(151) * .04 >= a) & (np.arange(151) * .04 < b)]) /
                                            np.sum(dia_t[(np.arange(151) * .04 >= a) & (np.arange(151) * .04 < b)]))) for a, b in blocks}
    rng = np.random.default_rng(5); ratios = []
    for _ in range(2000):
        a = rng.normal(size=len(order)) + 1j * rng.normal(size=len(order))
        src = np.zeros(len(order), complex)
        for c, k in enumerate(order):  # hermitian: src(-k) = conj src(k)
            src[c] = a[c]
        for c, k in enumerate(order):
            cneg = order.index(neg(k))
            if cneg > c: src[cneg] = np.conj(src[c])
        R = np.einsum('tabpk,k,p->tab', G, src, src.conj()) / np.sum(abs(src) ** 2)
        n_, a_ = chan(R); ratios.append(np.linalg.norm(a_) / np.linalg.norm(n_))
    ratios = np.array(ratios)
    res[tag] = dict(s0=s0, shell=1, exact_trace_ratio=float(np.linalg.norm(an) / np.linalg.norm(nr)),
                    offdiag_over_diag=float(offw), per_mode_normal_spread=spread, leakage_into_shells=leak,
                    whitenoise_estimator_ratio_percentiles={p: float(np.percentile(ratios, p)) for p in (5, 16, 50, 84, 95)},
                    late_over_early_normal=float(np.max(abs(nr[100:])) / np.max(abs(nr[:50]))),
                    scattered_over_diagonal_by_time=scat_t,
                    per_mode_late_over_early_median=float(np.median([np.max(abs(nrm[c][100:])) / np.max(abs(nrm[c][:50])) for c in range(len(order))])))
    out = np.full((151, 4), np.nan, complex); outA = out.copy(); out[:, 0] = nr; outA[:, 0] = an
    np.savez(D + f'MEAS_{tag}.npz', s0=s0, normal=out, anomalous=outA, per_mode_normal=nrm, per_mode_anomalous=anm)
    print(tag, json.dumps(res[tag], indent=1))

for f in sorted(glob.glob(D + 'PROJ_hutch_*.npz')):
    d = np.load(f); proj = d['proj']; srcs = d['srcs']; sh = d['shell']; s0 = float(d['s0']); tag = os.path.basename(f)[5:-4]
    M = proj.shape[1] // 2
    per = np.zeros((M, 151, 4, 2, 2), complex)
    for j in range(M):
        for b in range(2):
            q = proj[:, 2 * j + b]  # [t, a, mode]
            for s_ in range(1, 5):
                sel = sh == s_
                per[j, :, s_ - 1, :, b] = np.mean(q[:, :, sel] * srcs[j, sel].conj()[None, None, :], axis=-1)
    T = per.mean(0); t_ = np.arange(151) * .04
    rows = []
    out = np.zeros((151, 4), complex); outA = out.copy()
    for s_ in range(4):
        nr, an = chan(T[:, s_]); out[:, s_] = nr; outA[:, s_] = an
        single = [np.linalg.norm(chan(per[j, :, s_])[1]) / np.linalg.norm(chan(per[j, :, s_])[0]) for j in range(M)]
        # jackknife SE of the ratio
        jk = []
        for j in range(M):
            Tj = np.delete(per, j, axis=0).mean(0)[:, s_]; n_, a_ = chan(Tj); jk.append(np.linalg.norm(a_) / np.linalg.norm(n_))
        jk = np.array(jk); se = float(np.sqrt((M - 1) / M * np.sum((jk - jk.mean()) ** 2)))
        nr_single = [chan(per[j, :, s_])[0] for j in range(M)]
        noise = float(np.mean([np.linalg.norm(x - nr) for x in nr_single]) / np.linalg.norm(nr))
        # incoherent remainder per mode after removing the trace part (upper bound on scattering from shells 1-4 kicks)
        sel = sh == s_ + 1; scat = {}
        for a0, b0 in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6.01)]:
            tt = (t_ >= a0) & (t_ < b0); num = 0.; den = 0.
            for j in range(M):
                for b in range(2):
                    q = proj[tt][:, 2 * j + b][:, :, sel]  # [t, a, mode]
                    coh = T[tt][:, s_, :, b][:, :, None] * srcs[j, sel][None, None, :]
                    num += np.sum(abs(q - coh) ** 2); den += np.sum(abs(coh) ** 2)
            scat[f'{a0}-{b0:g}'] = float(np.sqrt(num / den))
        rows.append(dict(shell=s_ + 1, incoherent_over_trace_by_time=scat, ratio_mean_trace=float(np.linalg.norm(an) / np.linalg.norm(nr)), ratio_jackknife_se=se,
                         single_probe_ratio_range=[float(min(single)), float(max(single))], single_probe_normal_rel_noise=noise,
                         late_over_early_normal=float(np.max(abs(nr[100:])) / np.max(abs(nr[:50])))))
    res[tag] = dict(s0=s0, M=M, rows=rows)
    np.savez(D + f'MEAS_{tag}.npz', s0=s0, normal=out, anomalous=outA)
    print(tag, json.dumps(res[tag], indent=1))
json.dump(res, open(OUT, 'w'), indent=1)
