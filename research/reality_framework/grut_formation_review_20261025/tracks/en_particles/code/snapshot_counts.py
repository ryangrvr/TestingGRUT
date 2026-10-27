"""en_particles part 1b (own code): per-mode exact reference counts from a regenerated engine snapshot.

Formulation independent of Codex's analyze.py: full FFT of psi=chi1+i chi2 (forward norm, so sums are densities),
N_k(w) = (|pi_k|^2 + w^2 |psi_k|^2)/(2w), Q_k = Im(psi_k^* pi_k), n_pm = (N_k +- Q_k)/2 summed over k != 0.
Also rebuilds Codex's binned export from the snapshot with own code and compares it to SPECTRAL_EXPORT.npz.
Usage: python3 -I snapshot_counts.py <out.json> <export.npz or -> <state.bin> [<state.bin> ...]
"""
import sys, json, struct
import numpy as np

LX, G2, LP = 1e-3, 1e-2, 3e-2


def read_state(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    return N, L, s, np.fromfile(path, dtype='<f8', offset=24).reshape(6, N, N, N)


def fd4(N, L):
    a = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    f = lambda x: 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    return f(a)[:, None, None]+f(a)[None, :, None]+f(a)[None, None, :]


def counts(a, N, L):
    psi = np.fft.fftn(a[0]+1j*a[1], norm='forward'); pp = np.fft.fftn(a[3]+1j*a[4], norm='forward')
    phi = np.fft.fftn(a[2], norm='forward'); k2 = fd4(N, L); nz = k2 > 0
    Qk = np.imag(np.conj(psi)*pp); Q = float(Qk.sum()); Qnz = float(Qk[nz].sum())
    X2 = float((np.abs(psi)**2).sum()); f2 = float((np.abs(phi)**2).sum())
    mH2 = 2*LX*X2+G2*f2
    out = dict(Q=Q, Q_nonzero_over_Q=Qnz/Q, mH2=mH2, X2=X2, phi2=f2)
    kc = np.sqrt(k2)
    for name, m2 in [('massless', 0.), ('hartree', mH2)]:
        om = np.sqrt(np.where(nz, k2+m2, 1.))
        Nk = (np.abs(pp)**2/om+om*np.abs(psi)**2)/2
        nm = (.5*(Nk-Qk))[nz]; npl = (.5*(Nk+Qk))[nz]
        out[name] = dict(nminus_over_Q=float(nm.sum())/Q, nplus_over_Q=float(npl.sum())/Q)
        # per-mode occupation profile (numbers per mode, V=L^3): shells in FD4 |k|
        edges = np.arange(0, kc.max()+.5, .5); ids = np.digitize(kc[nz], edges)-1
        cnt = np.bincount(ids, minlength=len(edges)); occ = np.bincount(ids, weights=nm, minlength=len(edges))*L**3
        out[name]['occupation_minus_per_mode'] = [[float(edges[i]+.25), float(occ[i]/cnt[i])] for i in range(len(edges)) if cnt[i] > 0]
        out[name]['nminus_from_modes_with_occ_below_1_over_Q'] = float(
            sum(occ[i] for i in range(len(edges)) if cnt[i] > 0 and occ[i]/cnt[i] < 1)/L**3/Q)
    # per-mode Cauchy-Schwarz floor: N_k >= |pi_k||psi_k| for any positive w per mode
    Nmin = np.abs(pp)*np.abs(psi)
    out['per_mode_optimal_w_floor_nminus_over_Q'] = float((.5*(Nmin-Qk))[nz].sum())/Q
    return out


def rebuild_export(a, N, L):
    q = a[:3]; p = a[3:]
    qh = np.fft.rfftn(q, axes=(-3, -2, -1), norm='forward'); ph = np.fft.rfftn(p, axes=(-3, -2, -1), norm='forward')
    v = np.array([qh[0], qh[1], ph[0], ph[1], qh[2], ph[2]])
    f = np.fft.fftfreq(N)*2*np.pi*N/L; last = np.fft.rfftfreq(N)*2*np.pi*N/L
    kp = np.sqrt(f[:, None, None]**2+f[None, :, None]**2+last[None, None, :]**2)
    edges = np.linspace(0, kp.max()*1.00001, 51); ids = np.clip(np.digitize(kp.ravel(), edges)-1, 0, 49)
    w = np.full(N//2+1, 2.); w[[0, -1]] = 1.; w = np.broadcast_to(w[None, None, :], kp.shape).copy(); w[0, 0, 0] = 0
    cov = np.zeros((50, 6, 6))
    for i in range(6):
        for j in range(6):
            cov[:, i, j] = np.bincount(ids, weights=(w*np.real(np.conj(v[i])*v[j])).ravel(), minlength=50)
    return cov


def main(out, export, *paths):
    res = dict(source='own code on regenerated (reproduction-type) engine snapshots', rows=[])
    for p in paths:
        N, L, s, a = read_state(p)
        r = counts(a, N, L); r.update(path=p.split('/')[-1], s=s, N=N, L=L)
        if export != '-' and abs(s-23.915481855003456) < 1e-6:
            cov = rebuild_export(a, N, L); ref = np.load(export)['N48_covariance']
            r['export_rebuild_max_abs_diff_over_max'] = float(np.abs(cov-ref).max()/np.abs(ref).max())
        res['rows'].append(r)
        print(p.split('/')[-1], 's=%.4f' % s, 'Q=%.6f' % r['Q'], 'massless n-/Q=%.6f' % r['massless']['nminus_over_Q'],
              'hartree n-/Q=%.6f' % r['hartree']['nminus_over_Q'], 'floor=%.4f' % r['per_mode_optimal_w_floor_nminus_over_Q'],
              'mH2=%.3f' % r['mH2'], 'Qnz/Q=%.4f' % r['Q_nonzero_over_Q'], r.get('export_rebuild_max_abs_diff_over_max', ''), flush=True)
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
