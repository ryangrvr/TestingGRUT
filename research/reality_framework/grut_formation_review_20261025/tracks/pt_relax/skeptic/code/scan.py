"""pt_relax skeptic (Claude): cutoff scan of the classical negative-frequency action at fixed box and seed.

Question: is the post-plateau decay of A_- that pt_relax reports on N32 (L = 12) and N48 a continuum classical
effect, or does it shrink toward zero as the lattice cutoff rises?  Same L = 12, seed 1, width 1 initial state
(Codex's prepare.initialize draws each seeded k <= 8 mode from its own (seed, L, i, j, k) SeedSequence, so the
infrared realization is identical on every N >= 32), N = 32, 40, 48, 64.

Labelling.
  * Dynamics: Codex's FD4/Yoshida engine (fieldstep.cpp via pt_common, PRIVATE git-archive copy of 5f5513b in
    $RF).  This is a reproduction tool, cross-checked by pt_relax against our own BM6 integrator to 2e-8.
  * Everything else here is own code, written fresh for this check (not copied from pt_relax/code):
    sampling, both window estimators, band/shell bookkeeping, energy fractions and occupations.

Estimators per 10-s window (M = 251 samples at ds = 0.04, inclusive ends), per initial net charge density nQ:
  A_-^{w} = sum_k sum_{omega_n < 0} |omega_n| |sum_j w_j psi_k(s_j) e^{-i omega_n s_j}|^2 / (M sum_j w_j^2),
  psi = chi1 + i chi2, forward-normalised spatial FFT, w = Hann (Codex definition) or 4-term Blackman-Harris.
  Bands: zero mode, 0 < |k_cont| <= 8, |k_cont| > 8.  A_+ likewise; A_+ - A_- = 1 checks charge.

Usage:  python3 -B scan.py N S_TARGET WALL_BUDGET_SECONDS     (restartable; state in $RUNDIR/L12_N<N>/)
"""
import os, sys, json, time
import numpy as np
from scipy.fft import fftn, fft

sys.path.insert(0, os.path.join(os.environ['PK'], 'code'))
import pt_common as P            # noqa: E402  Codex engine (reproduction tool)
from prepare import initialize   # noqa: E402  Codex initial state (same private copy)

L, SEED, WIDTH, HMAX = 12.0, 1, 1, 0.004
DS, M = 0.04, 251
LAMX, G2, LAMP = 1e-3, 1e-2, 3e-2


def window_bh(M):
    n = np.arange(M); x = 2 * np.pi * n / (M - 1)
    return 0.35875 - 0.48829 * np.cos(x) + 0.14128 * np.cos(2 * x) - 0.01168 * np.cos(3 * x)


class Lat:
    def __init__(self, N, nQ):
        self.N, self.V, self.nQ = N, L ** 3, nQ
        i = np.fft.fftfreq(N, d=1.0 / N)
        dk = 2 * np.pi / L
        k2 = dk ** 2 * (i[:, None, None] ** 2 + i[None, :, None] ** 2 + i[None, None, :] ** 2)
        self.kc = np.sqrt(k2)
        dx = L / N
        s2 = np.sin(np.pi * i / N) ** 2                       # sin^2(k dx / 2) per axis
        sym = 4 * s2 * (1 + s2 / 3) / dx ** 2                  # FD4 -Laplacian eigenvalue per axis
        self.klat2 = sym[:, None, None] + sym[None, :, None] + sym[None, None, :]
        self.zero = self.kc < 1e-12
        self.band = (~self.zero) & (self.kc <= 8 + 1e-9)
        self.uv = self.kc > 8 + 1e-9
        self.shell = np.rint(self.kc / dk).astype(int)
        self.nsh = int(self.shell.max()) + 1
        self.cnt = np.bincount(self.shell.ravel(), minlength=self.nsh)
        self.om = 2 * np.pi * np.fft.fftfreq(M, d=DS)
        self.wins = dict(hann=np.hanning(M), bh=window_bh(M))

    def sample(self, z):
        q, p = z[:3], z[3:]
        qk = fftn(q, axes=(1, 2, 3), norm='forward'); pk = fftn(p, axes=(1, 2, 3), norm='forward')
        r2 = q[0] ** 2 + q[1] ** 2; f2 = q[2] ** 2
        kin = 0.5 * (np.abs(pk) ** 2); grad = 0.5 * self.klat2 * np.abs(qk) ** 2
        quad = kin + grad                                       # (3, N, N, N)
        pot = np.array([LAMX / 4 * np.mean(r2 * r2), LAMP / 4 * np.mean(f2 * f2), G2 / 2 * np.mean(f2 * r2)])
        rho = float(quad.sum() + pot.sum())
        Q = float(np.mean(q[0] * p[1] - q[1] * p[0]))
        mx2 = 2 * LAMX * float(np.mean(r2)) + G2 * float(np.mean(f2))
        mp2 = G2 * float(np.mean(r2)) + 3 * LAMP * float(np.mean(f2))
        row = dict(rho=rho, Q=Q, mx2=mx2, mp2=mp2,
                   uv_quad_frac=float(quad[:, self.uv].sum() / rho),
                   phi_frac=float((quad[2].sum() + pot[1] + pot[2]) / rho))
        psi = qk[0] + 1j * qk[1]
        return row, psi, qk, pk

    def occupations(self, qk, pk, mx2, mp2):
        """mean classical occupation per real oscillator in each shell, Hartree masses, FD4 dispersion"""
        sh = self.shell.ravel(); nz = (~self.zero).ravel()
        out = {}
        for name, comps, m2 in [('X', [0, 1], mx2), ('phi', [2], mp2)]:
            w = np.sqrt(self.klat2 + m2)
            occ = sum(self.V * (np.abs(pk[c]) ** 2 + w ** 2 * np.abs(qk[c]) ** 2) / (2 * w) for c in comps) / len(comps)
            out[name] = (np.bincount(sh[nz], weights=occ.ravel()[nz], minlength=self.nsh) / np.maximum(self.cnt, 1)).tolist()
            out['E_' + name] = np.bincount(sh[nz], weights=sum(
                0.5 * (np.abs(pk[c]) ** 2 + self.klat2 * np.abs(qk[c]) ** 2) for c in comps).ravel()[nz],
                minlength=self.nsh).tolist()
        return out

    def window(self, buf):
        res = {}
        neg = self.om < 0; pos = self.om > 0
        for wn, w in self.wins.items():
            norm = M * np.sum(w * w)
            Am = np.zeros(buf.shape[1:]); Ap = np.zeros(buf.shape[1:])
            for i in range(buf.shape[1]):
                X = fft(w[:, None, None] * buf[:, i], axis=0)
                Pw = np.abs(X) ** 2 / norm
                Am[i] = np.tensordot(np.abs(self.om) * neg, Pw, axes=(0, 0))
                Ap[i] = np.tensordot(np.abs(self.om) * pos, Pw, axes=(0, 0))
            q = self.nQ
            r = dict(coh=float(Am[self.zero].sum() / q), band=float(Am[self.band].sum() / q),
                     uv=float(Am[self.uv].sum() / q))
            r['tot8'] = r['coh'] + r['band']; r['all'] = r['tot8'] + r['uv']
            r['charge'] = float((Ap.sum() - Am.sum()) / q)
            if wn == 'hann':
                r['shell_minus'] = (np.bincount(self.shell.ravel(), weights=Am.ravel(), minlength=self.nsh) / q).tolist()
            res[wn] = r
        return res


def main():
    N = int(sys.argv[1]); s_target = float(sys.argv[2]); budget = float(sys.argv[3])
    d = os.path.join(os.environ['RUNDIR'], f'L12_N{N}'); os.makedirs(d, exist_ok=True)
    rst = os.path.join(d, 'restart.npz')
    P.set_box(L); lib = P.library(); thr = int(os.environ.get('GRUT_PARTICLE_THREADS', '2'))
    if os.path.exists(rst):
        R = np.load(rst); z = np.ascontiguousarray(R['z']); i0 = int(R['i']); nQ = float(R['nQ']); rho0 = float(R['rho0'])
        part = int(R['part']) + 1
    else:
        z, par = initialize(N, L, SEED, WIDTH); z = np.ascontiguousarray(z); i0 = 0; part = 0
        nQ = float(np.mean(z[0] * z[4] - z[1] * z[3])); rho0 = None
    lat = Lat(N, nQ)
    row, psi, qk, pk = lat.sample(z)
    if rho0 is None: rho0 = row['rho']
    # window ends: 10, 20, 24, 34, 44, ...  Restarts only at ends >= 34, so a window never straddles a restart.
    def is_end(i):
        return i in (250, 500) or (i >= 600 and (i - 600) % 250 == 0)
    buf = np.zeros((M,) + psi.shape, complex); nb = 0
    hist = {} if i0 == 0 else None
    def push(ps):
        nonlocal nb
        if nb < M: buf[nb] = ps; nb += 1
        else: buf[:-1] = buf[1:]; buf[-1] = ps
    push(psi)
    rows = [dict(s=i0 * DS, **row)]; occs = []; wins = []
    if i0 == 0: occs.append(dict(s=0.0, **lat.occupations(qk, pk, row['mx2'], row['mp2'])))
    t0 = time.monotonic(); i = i0
    while True:
        P.step(lib, z, DS, HMAX, threads=thr); i += 1; s = i * DS
        row, psi, qk, pk = lat.sample(z)
        if i % 5 == 0: rows.append(dict(s=s, **row))
        if i % 250 == 0: occs.append(dict(s=s, **lat.occupations(qk, pk, row['mx2'], row['mp2'])))
        push(psi)
        if is_end(i) and nb == M:
            w = lat.window(buf); w.update(window=[round(s - 10, 6), round(s, 6)], dE=row['rho'] / rho0 - 1, dQ=row['Q'] / nQ - 1,
                                          uv_quad_frac=row['uv_quad_frac'], phi_frac=row['phi_frac'])
            wins.append(w)
            print(f"N{N} WIN {s-10:6.1f}-{s:6.1f} hann all {w['hann']['all']:.6f} tot8 {w['hann']['tot8']:.6f} "
                  f"coh {w['hann']['coh']:.5f} band {w['hann']['band']:.5f} uv {w['hann']['uv']:.5f} | bh all {w['bh']['all']:.6f} "
                  f"| chg {w['hann']['charge']:.8f} dE {w['dE']:.1e} Euv {row['uv_quad_frac']:.4f} wall {time.monotonic()-t0:.0f}", flush=True)
            if i >= 600:
                buf[0] = buf[nb - 1]; nb = 1                     # contiguous 10-s windows from s = 24 on
                if s >= s_target - 1e-9 or time.monotonic() - t0 > budget:
                    json.dump(dict(N=N, nQ=nQ, rho0=rho0, rows=rows, occs=occs, wins=wins),
                              open(os.path.join(d, f'part{part:03d}.json'), 'w'))
                    np.savez(rst, z=z, i=i, nQ=nQ, rho0=rho0, part=part)
                    print('SAVED', N, 'part', part, 's', s, flush=True)
                    return


if __name__ == '__main__':
    main()
