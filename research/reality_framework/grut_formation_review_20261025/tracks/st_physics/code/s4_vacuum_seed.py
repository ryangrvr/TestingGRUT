"""st_physics step 4: vacuum-seeded response estimate [E], using Codex's compiled tangent engine (response.c,
built by timing_common.build() in our own GRUT_TIMING_WORKDIR) with OUR seeds and OUR post-processing.

Label: engine = Codex code (reproduction-type use); seeds, covariance and diagnostics = own.

For each momentum k the engine propagates 12 unit seed vectors (the columns of the 12x12 transfer matrix) plus the
supplied preparation-phase seed (cross-check against the published amplitude arrays).  The initial state at a = 1
is the instantaneous adiabatic vacuum of the conformal canonical fields: in the eigenbasis of the initial conformal
mass matrix M0, <dq dq> = 1/(2 w), <dp dp> = w/2 with w^2 = k^2 + m_i^2 (all positive for k >= 1.8).  The
dimensionless power per ln k of an observable O = r . dv is k^3/(2 pi^2) r^T T Cov T^T r.

usage: python3 -B s4_vacuum_seed.py <codex timing code dir> <out.json>   (needs the GRUT_* env of BRIEF_TIMING.md)
"""
import sys, json, subprocess, hashlib
import numpy as np

code_dir, out = sys.argv[1:3]
sys.path.insert(0, code_dir)
import timing_common as T  # noqa: E402  (Codex module: background preparation, build of response.c)
M = T.M

KV = np.array([1.8, 2.0, 2.13, 2.4, 2.8])


def run(exe, tag, beta, decay, ks, seeds, times, e=1e-7, rtol=2e-10):
    v, res = M.prepared(beta)
    params = np.array([beta, decay, e, M.XI, M.ETA2, M.MU, M.KAP, times[-1], rtol], dtype='<f8')
    inp = T.WORK / (tag + '_stphys.in'); outp = T.WORK / (tag + '_stphys.bin')
    with inp.open('wb') as f:
        np.array([len(ks), len(times)], dtype='<i4').tofile(f); params.tofile(f); np.asarray(v, dtype='<f8').tofile(f)
        np.asarray(ks, dtype='<f8').tofile(f); np.asarray(times, dtype='<f8').tofile(f)
        np.asarray(seeds, dtype='<f8').tofile(f)
    proc = subprocess.run([str(exe), str(inp), str(outp)], capture_output=True, text=True, check=True)
    data = np.fromfile(outp, dtype='<f8').reshape(len(times), 20 + 12 * len(ks))
    assert np.allclose(data[:, 0], times, rtol=0, atol=1e-10)
    return data[:, 1:20], data[:, 20:].reshape(len(times), len(ks), 12), proc.stdout.strip()


def envelope(times, a, C, y, T_eta, lo=0.0):
    eta = a / C
    sel = times >= lo
    edges = np.arange(eta[sel][0], eta[sel][-1], T_eta)
    tc, ev = [], []
    for e0 in edges[:-1]:
        w = sel & (eta >= e0) & (eta < e0 + T_eta)
        if w.sum() >= 2:
            tc.append(times[w][-1]); ev.append(np.max(abs(y[w])))
    return np.array(tc), np.array(ev)


def main():
    exe = T.build()
    pkt = T.ROOT
    resj = {r['tag']: r for r in T.read('RESPONSE_RESULTS.json')['rows']}
    result = dict(label='engine: Codex response.c (reproduction-type use); seeds/covariance/analysis: own [E]', rows={})
    for tag, beta, decay, T_eta in (('quiet3', 3., .5, 0.8143), ('quiet5', 5., 1., 0.9208)):
        pub = np.load(pkt / f'{tag}.npz')
        times = pub['times']
        ks = np.repeat(KV, 13)
        seeds = np.zeros((len(ks), 12))
        sup = T.phase_seed(beta)
        for i in range(len(KV)):
            seeds[13 * i:13 * i + 12] = np.eye(12)
            seeds[13 * i + 12] = sup
        bg, vec, check = run(exe, tag, beta, decay, ks, seeds, times)
        bgdiff = float(np.max(abs(bg[:, 2:4] - pub['bg'][:, 2:4])) / np.max(abs(pub['bg'][:, 2:4])))
        print(tag, 'max X background difference (relative to max |x|) vs published:', bgdiff, flush=True)
        amp_pub, dq_pub, nx = T.responses(beta, bg, vec)  # amplitude/charge of every column (Codex's formula)
        a = np.exp(bg[:, 15]); C = resj[tag]['late_C']
        norm = abs(resj[tag]['late_NX_normalization_through80'])
        v0, _ = M.prepared(beta)
        M0, g0, gp0, chi0, chip0 = M.matrices(v0, beta)
        lam, U = np.linalg.eigh(M0)
        rows = []
        Pamp = np.zeros((len(times), len(KV))); Pchg = Pamp.copy(); PampX = Pamp.copy()
        for i, k in enumerate(KV):
            sl = slice(13 * i, 13 * i + 12)
            w = np.sqrt(k * k + lam)
            assert np.all(k * k + lam > 0)
            Cqq = (U / (2 * w)) @ U.T; Cpp = (U * (w / 2)) @ U.T
            Cov = np.block([[Cqq, np.zeros((6, 6))], [np.zeros((6, 6)), Cpp]])
            CovX = np.zeros((12, 12)); CovX[np.ix_([0, 1, 6, 7], [0, 1, 6, 7])] = Cov[np.ix_([0, 1, 6, 7], [0, 1, 6, 7])]
            # observable responses of the 12 unit columns: shape (t, 12)
            A = amp_pub[:, sl]; Q = dq_pub[:, sl] / norm
            Pamp[:, i] = k ** 3 / (2 * np.pi ** 2) * np.einsum('ti,ij,tj->t', A, Cov, A)
            PampX[:, i] = k ** 3 / (2 * np.pi ** 2) * np.einsum('ti,ij,tj->t', A, CovX, A)
            Pchg[:, i] = k ** 3 / (2 * np.pi ** 2) * np.einsum('ti,ij,tj->t', Q, Cov, Q)
            # supplied-seed column vs published arrays (only for k in Codex's grid)
            supamp = amp_pub[:, 13 * i + 12]
            j = np.flatnonzero(abs(pub['ks'] - k) < 1e-12)
            pubdiff = None
            if len(j):
                ref = pub['amplitude_per_radian'][:, j[0]]
                pubdiff = float(np.max(abs(supamp - ref)) / np.max(abs(ref)))
            # linearity check: T . seed from the unit columns vs the directly propagated supplied seed
            lin = float(np.max(abs(A @ sup - supamp)) / np.max(abs(supamp)))
            tc, e_vac = envelope(times, a, C, np.sqrt(Pamp[:, i]), T_eta, lo=300)
            _, e_sup = envelope(times, a, C, supamp, T_eta, lo=300)
            _, e_vacX = envelope(times, a, C, np.sqrt(PampX[:, i]), T_eta, lo=300)
            ratio = e_vac / e_sup
            rows.append(dict(k=float(k), omega_initial=w.tolist(),
                             naive_phase_seed_per_lnk=float(k / (2 * np.pi) / np.linalg.norm(chi0[:2])),
                             equivalent_supplied_radians_per_lnk_late=float(np.median(ratio)),
                             equivalent_supplied_radians_per_lnk_range=[float(ratio.min()), float(ratio.max())],
                             X_only_vacuum_fraction_of_power=float(np.median((e_vacX / e_vac) ** 2)),
                             published_amplitude_rel_diff=pubdiff, linearity_rel_diff=lin))
            print(tag, rows[-1], flush=True)
        # band-integrated RMS (trapezoid in ln k over 1.8-2.8; the band below 1.8 is omitted -> lower bound)
        lk = np.log(KV)
        s_amp = np.sqrt(np.trapezoid(Pamp, lk, axis=1)); s_chg = np.sqrt(np.trapezoid(Pchg, lk, axis=1))
        firsts = {}
        for name, s in (('amplitude_rms', s_amp), ('charge_rms', s_chg)):
            for lev in (0.1, 0.3, 1.0):
                wh = np.flatnonzero(s >= lev)
                firsts[f'{name}_{lev}'] = None if not len(wh) else dict(tau=float(times[wh[0]]), a=float(a[wh[0]]))
        result['rows'][tag] = dict(per_k=rows, background_X_rel_diff_vs_published=bgdiff, band_rms_first_sampled_crossings=firsts, engine_check=check,
                                   rms_at=dict((str(t), dict(amp=float(s_amp[np.argmin(abs(times - t))]),
                                                             charge=float(s_chg[np.argmin(abs(times - t))])))
                                               for t in (80, 134, 200, 300)),
                                   initial_mass_eigenvalues=lam.tolist())
        print(tag, firsts, flush=True)
    json.dump(result, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
