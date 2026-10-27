"""Track pt_relax (Claude): long continuation of one IR12 history and the windowed negative-frequency action.

Labelling.
  * Evolution: Codex's reviewed FD4/Yoshida engine (fieldstep.cpp, compiled by pt_common.library from a PRIVATE
    git-archive copy of 5f5513b at $RF). This is a reproduction tool, not independent dynamics. The initial field
    comes from Codex's prepare.initialize (same copy).  An independent BM6 cross-check segment is in xcheck.py.
  * Everything else is own code: sampling, the Hann-window frequency action, shell spectra, occupations, restart.

Estimator (Codex MATHEMATICS.md, re-implemented):  psi = chi1 + i chi2, forward-normalised spatial FFT psi_k(s_j),
  P_k(w_n) = |sum_j h_j psi_k(s_j) e^{-i w_n s_j}|^2 / (M sum h^2),   A_- = sum_{k, w<0} |w| P_k(w)   (per volume),
  reported per initial net charge density nQ.  Windows of 10 in s, sampled at ds = 0.04 (M = 251, Hann).
  Bands: zero mode ('coherent'), 0 < k_cont <= 8 (Codex's band, 'waves'), and all nonzero lattice modes.

Usage:  python3 -B long_run.py TAG S_TARGET [WALL_BUDGET_SECONDS]
"""
import os, sys, json, time
import numpy as np
from scipy.fft import fftn, fft

PK = os.environ['PK']
sys.path.insert(0, PK + '/code')
import pt_common as P  # noqa: E402  (Codex module from the private copy: engine + geometry only)
from prepare import initialize  # noqa: E402

LX, G2, LP = 1e-3, 1e-2, 3e-2
DS, M = 0.04, 251                     # sample spacing and samples per 10-s window (inclusive ends)
CASES = {'IR12_s1_w1': dict(N=32, L=12., seed=1, width=1, hmax=.004),
         'IR12_s2_w1': dict(N=32, L=12., seed=2, width=1, hmax=.004),
         'IR12_fine_s1_w1': dict(N=48, L=12., seed=1, width=1, hmax=.004)}
RUN = os.environ['RUNDIR']


def grids(N, L):
    m = np.fft.fftfreq(N) * N
    kc2 = (2 * np.pi / L) ** 2 * (m[:, None, None] ** 2 + m[None, :, None] ** 2 + m[None, None, :] ** 2)
    u = 2 * np.pi * m / N
    f = 4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3) / (L / N) ** 2
    kl = f[:, None, None] + f[None, :, None] + f[None, None, :]          # FD4 symbol (exact lattice -Laplacian)
    kc = np.sqrt(kc2)
    shell = np.floor(kc / (2 * np.pi / L) + .5).astype(int)            # Codex's shell convention
    return kc, kl, shell


class Obs:
    def __init__(self, N, L, nQ0):
        self.N, self.L, self.V, self.nQ0 = N, L, L ** 3, nQ0
        self.kc, self.kl, self.shell = grids(N, L)
        self.nsh = self.shell.max() + 1
        self.cnt = np.bincount(self.shell.ravel(), minlength=self.nsh)
        self.band = (self.kc <= 8 + 1e-12) & (self.kc > 0)
        self.zero = (self.kc == 0)
        self.uv = self.kc > 8 + 1e-12
        w = 2 * np.pi * np.fft.fftfreq(M, DS)
        self.w = w
        self.h = np.hanning(M)
        self.norm = M * np.sum(self.h ** 2)

    def equal_time(self, z, full=False):
        """Own energy/charge and equal-time diagnostics. Forward-normalised FFTs (Parseval: mean = sum)."""
        N3 = self.N ** 3
        q = z[:3]; p = z[3:]
        r2 = q[0] ** 2 + q[1] ** 2; ph2 = q[2] ** 2
        qh = fftn(q, axes=(-3, -2, -1), norm='forward'); ph = fftn(p, axes=(-3, -2, -1), norm='forward')
        K = 0.5 * np.sum(np.abs(ph) ** 2, axis=0)                      # per k, summed over fields
        Gx = 0.5 * self.kl * (np.abs(qh[0]) ** 2 + np.abs(qh[1]) ** 2)
        Gp = 0.5 * self.kl * np.abs(qh[2]) ** 2
        UX = LX / 4 * np.mean(r2 * r2); UP = LP / 4 * np.mean(ph2 * ph2); UM = G2 / 2 * np.mean(ph2 * r2)
        rho = float(K.sum() + Gx.sum() + Gp.sum() + UX + UP + UM)
        Q = float(np.mean(q[0] * p[1] - q[1] * p[0]))
        mr2 = float(np.mean(r2)); mph2 = float(np.mean(ph2))
        mH2 = 2 * LX * mr2 + G2 * mph2                                   # Codex's isotropic Hartree X mass^2
        mphi2 = G2 * mr2 + 3 * LP * mph2                                 # phi Hartree mass^2
        means = z.mean(axis=(1, 2, 3))
        E0 = .5 * np.sum(means[3:] ** 2) + LX / 4 * (means[0] ** 2 + means[1] ** 2) ** 2 + LP / 4 * means[2] ** 4 \
            + G2 / 2 * means[2] ** 2 * (means[0] ** 2 + means[1] ** 2)
        Ephi_quad = 0.5 * np.abs(ph[2]) ** 2 + Gp
        row = dict(rho=rho, Q=Q, r2=mr2, phi2=mph2, mH2=mH2, mphi2=mphi2,
                   outside_mean=1 - E0 / rho, phi_quad_frac=float(Ephi_quad[~self.zero].sum() / rho),
                   phi_total_frac=float((0.5 * np.sum(np.abs(ph[2]) ** 2) + Gp.sum() + UP + UM) / rho),
                   Q0_over_Q=float((means[0] * means[4] - means[1] * means[3]) / Q),
                   cond_r=float(np.hypot(means[0], means[1])))
        # quadratic energy outside k_cont > 8 (X and phi; kinetic + FD4 gradient), fraction of total energy
        EX = 0.5 * (np.abs(ph[0]) ** 2 + np.abs(ph[1]) ** 2) + Gx
        row['E_uv_frac'] = float((EX[self.uv].sum() + Ephi_quad[self.uv].sum()) / rho)
        # equal-time free-basis antiparticle reference at the Hartree mass (own form of Codex's ref_H)
        om = np.sqrt(self.kl + mH2)
        psi = qh[0] + 1j * qh[1]; pp = ph[0] + 1j * ph[1]
        B = 0.5 * (psi + 1j * pp / om); A = 0.5 * (psi - 1j * pp / om)
        nminus = om * np.abs(B) ** 2; nplus = om * np.abs(A) ** 2       # per volume, per mode
        row['ref_minus_all'] = float(nminus[~self.zero].sum() / self.nQ0)
        row['ref_minus_band'] = float(nminus[self.band].sum() / self.nQ0)
        row['ref_minus_uv'] = float(nminus[self.uv].sum() / self.nQ0)
        if not full:
            return row, psi
        # shell spectra: energy and mean classical occupation per real oscillator (V (|p|^2 + w^2|q|^2)/(2w))
        sh = self.shell.ravel()
        omp = np.sqrt(self.kl + mphi2)
        occX = self.V * (np.abs(ph[0]) ** 2 + np.abs(ph[1]) ** 2 + om ** 2 * (np.abs(qh[0]) ** 2 + np.abs(qh[1]) ** 2)) / (2 * om) / 2
        occP = self.V * (np.abs(ph[2]) ** 2 + omp ** 2 * np.abs(qh[2]) ** 2) / (2 * omp)
        nz = (~self.zero).ravel()
        spec = dict(
            EX=np.bincount(sh[nz], weights=EX.ravel()[nz], minlength=self.nsh),
            Ephi=np.bincount(sh[nz], weights=Ephi_quad.ravel()[nz], minlength=self.nsh),
            occX=np.bincount(sh[nz], weights=occX.ravel()[nz], minlength=self.nsh) / np.maximum(self.cnt, 1),
            occphi=np.bincount(sh[nz], weights=occP.ravel()[nz], minlength=self.nsh) / np.maximum(self.cnt, 1),
            nminus=np.bincount(sh[nz], weights=nminus.ravel()[nz], minlength=self.nsh) / self.nQ0,
            nplus=np.bincount(sh[nz], weights=nplus.ravel()[nz], minlength=self.nsh) / self.nQ0)
        return row, psi, spec

    def window(self, buf):
        """Hann-window frequency action of a (M, N, N, N) complex buffer; per initial net charge."""
        out = {}
        Pw = np.zeros((M,) + buf.shape[1:])
        # chunk over the first spatial axis to bound memory
        for i in range(buf.shape[1]):
            X = fft(self.h[:, None, None] * buf[:, i], axis=0)
            Pw[:, i] = np.abs(X) ** 2 / self.norm
        wneg = (self.w < 0); wpos = (self.w > 0)
        Am = np.tensordot(np.abs(self.w) * wneg, Pw, axes=(0, 0))       # per k
        Ap = np.tensordot(np.abs(self.w) * wpos, Pw, axes=(0, 0))
        q = self.nQ0
        out['coherent_minus'] = float(Am[self.zero].sum() / q)
        out['band_minus'] = float(Am[self.band].sum() / q)
        out['uv_minus'] = float(Am[self.uv].sum() / q)
        out['coherent_plus'] = float(Ap[self.zero].sum() / q)
        out['band_plus'] = float(Ap[self.band].sum() / q)
        out['uv_plus'] = float(Ap[self.uv].sum() / q)
        out['total_minus_band'] = out['coherent_minus'] + out['band_minus']        # Codex's headline definition
        out['total_minus_all'] = out['total_minus_band'] + out['uv_minus']
        out['charge_check_all'] = (out['coherent_plus'] + out['band_plus'] + out['uv_plus']
                                   - out['total_minus_all'])
        sh = self.shell.ravel()
        out['shell_minus'] = (np.bincount(sh, weights=Am.ravel(), minlength=self.nsh) / q).tolist()
        # power-weighted mean |w| of the negative-frequency part in the band (where the antiparticle action sits)
        pw_neg = np.tensordot(wneg.astype(float), Pw, axes=(0, 0))
        out['mean_abs_w_minus_band'] = float(Am[self.band].sum() / max(pw_neg[self.band].sum(), 1e-300))
        return out


def main():
    tag = sys.argv[1]; s_target = float(sys.argv[2]); budget = float(sys.argv[3]) if len(sys.argv) > 3 else 540.
    case = CASES[tag]; N, L = case['N'], case['L']
    d = os.path.join(RUN, tag); os.makedirs(d, exist_ok=True)
    rst = os.path.join(d, 'restart.npz')
    P.set_box(L); lib = P.library(); threads = int(os.environ.get('GRUT_PARTICLE_THREADS', '2'))
    if os.path.exists(rst):
        R = np.load(rst); z = np.ascontiguousarray(R['z']); s0 = float(R['s']); nQ0 = float(R['nQ0']); rho0 = float(R['rho0'])
        part = int(R['part']) + 1
    else:
        z, params = initialize(N, L, case['seed'], case['width']); z = np.ascontiguousarray(z)
        np.save(os.path.join(d, 'initial.npy'), z)
        s0 = 0.; part = 0
        nQ0 = float(np.mean(z[0] * z[4] - z[1] * z[3]))
        rho0 = None
    ob = Obs(N, L, nQ0)
    row, psi, spec = ob.equal_time(z, full=True)
    if rho0 is None: rho0 = row['rho']
    # the window grid: ends at 10, 20, 24, 34, 44, ...; each window needs the 250 previous samples.
    # restarts happen only at ends >= 34 that start a new contiguous window, so the buffer starts at s0.
    buf = np.zeros((M, N, N, N), complex)
    nbuf = 0
    rows, specs, wins = [], [], []
    def push(psi):
        nonlocal nbuf
        if nbuf < M:
            buf[nbuf] = psi; nbuf += 1
        else:
            buf[:-1] = buf[1:]; buf[-1] = psi
    i0 = int(round(s0 / DS)); push(psi)
    rows.append(dict(s=s0, **row)); specs.append(dict(s=s0, **{k: v.tolist() for k, v in spec.items()}))
    t0 = time.monotonic(); i = i0
    while True:
        i += 1; s = i * DS
        P.step(lib, z, DS, case['hmax'], threads=threads)
        full = (i % 10 == 0)
        res = ob.equal_time(z, full=full)
        row, psi = res[0], res[1]
        rows.append(dict(s=s, **row))
        if full: specs.append(dict(s=s, **{k: v.tolist() for k, v in res[2].items()}))
        push(psi)
        ends = (abs(s - 10) < 1e-9 or abs(s - 20) < 1e-9 or (s > 23.99 and abs(((s - 24) / 10) - round((s - 24) / 10)) < 1e-9))
        if ends and nbuf == M:
            w = ob.window(buf); w.update(window=[s - 10, s], rho_rel=row['rho'] / rho0 - 1, Q_rel=row['Q'] / nQ0 - 1)
            wins.append(w)
            print(f'WIN {s-10:7.2f}-{s:7.2f} coh {w["coherent_minus"]:.5f} band {w["band_minus"]:.5f} '
                  f'tot8 {w["total_minus_band"]:.5f} totall {w["total_minus_all"]:.5f} uv {w["uv_minus"]:.5f} '
                  f'Euv {row["E_uv_frac"]:.4e} dE {row["rho"]/rho0-1:.2e} wall {time.monotonic()-t0:.0f}', flush=True)
            stop = (s >= s_target - 1e-9) or (time.monotonic() - t0 > budget)
            if stop and s >= 34 - 1e-9:
                np.savez_compressed(os.path.join(d, f'part{part:03d}.npz'),
                                    rows=json.dumps(rows), specs=json.dumps(specs), wins=json.dumps(wins))
                np.savez(rst, z=z, s=s, nQ0=nQ0, rho0=rho0, part=part)
                print('SAVED', tag, 'part', part, 's', s, 'wall', round(time.monotonic() - t0, 1), flush=True)
                return
            if abs(s - 24) < 1e-9 or s > 24:
                pass
            # after a window ending at >= 24, keep only the last sample: windows from 24 on are contiguous
            if s > 23.99:
                last = buf[nbuf - 1].copy(); buf[0] = last; nbuf = 1
        if i % 2500 == 0:
            print('t', s, 'wall', round(time.monotonic() - t0, 1), flush=True)


if __name__ == '__main__':
    main()
