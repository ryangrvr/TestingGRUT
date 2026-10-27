"""en_particles SKEPTIC (own code): re-estimate the measured-frequency counts with different estimators.

Input: a reproduction-type window.npz (uniform conformal-time samples psi_k, pi_k of psi = chi1 + i chi2, full FFT,
forward norm) written by en_particles' window_run.py driver around Codex's unmodified engine.

Estimators, all per net charge Q (Q measured from the data):
  A  'dpss'   : Slepian multitaper (NW=2.5, K=4) power, n_pm = sum_{w>0/<0} |w| P(w)   [different taper from Hann]
  B  'rect'   : rectangular window, same sum (leakage-prone; shows taper sensitivity)
  C  'virial' : per-mode time averages, n_- = 1/2 [ sqrt(<|pi|^2><|psi|^2>) - <Im psi* pi> ]  (no frequency split;
                exact for one line per branch at a common |w|, an UPPER bound on A-type counts by Cauchy-Schwarz)
  D  'centroid': per shell, branch centroid frequency w_c(k) = sum |w| P / sum P on each branch; reports the
                n_- weight carried by shells whose measured centroid deviates from the Hartree w_H(k) by > 10 %.
Also: n_- per shell group (box-scale k <= 1.6, 1.6 < k < 4, k >= 4), the zero-mode split, and the total.
Usage: python3 -I sk_windows.py <window.npz> <out.json> [SUB samples, default all] [STEP]
"""
import sys, json
import numpy as np
from scipy.signal.windows import dpss

LX, G2 = 1e-3, 1e-2


def k2_fd4(N, L):
    th = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    one = ((4/3)*(2-2*np.cos(th))-(2-2*np.cos(2*th))/12)/dx**2
    return (one[:, None, None]+one[None, :, None]+one[None, None, :]).ravel()


def nvec_shell(N):
    n = np.fft.fftfreq(N)*N
    return np.rint(np.sqrt(n[:, None, None]**2+n[None, :, None]**2+n[None, None, :]**2)).astype(int).ravel()


def spectra(x, ds, kind):
    """x: (M, modes) complex. Returns w (M,), P (M, modes) normalised so sum_w P = time-average |x|^2."""
    M = x.shape[0]; w = 2*np.pi*np.fft.fftfreq(M, ds)
    if kind == 'rect':
        X = np.fft.fft(x, axis=0); P = np.abs(X)**2/M**2
    elif kind == 'hann':
        h = np.hanning(M); X = np.fft.fft(x*h[:, None], axis=0); P = np.abs(X)**2/(M*(h**2).sum())
    else:
        tapers = dpss(M, 2.5, 4)  # unit-energy tapers
        P = 0.
        for t in tapers:
            X = np.fft.fft(x*t[:, None], axis=0); P = P+np.abs(X)**2/M
        P = P/len(tapers)
    return w, P


def analyse(s, psi, ppsi, phi2, N, L):
    M = len(s); ds = s[1]-s[0]
    k2 = k2_fd4(N, L); nz = k2 > 1e-12; shell = nvec_shell(N); kk = 2*np.pi/L*shell
    Q = float(np.mean(np.imag(np.conj(psi)*ppsi).sum(axis=1)))
    out = dict(s=[float(s[0]), float(s[-1])], samples=M, Q=Q)
    for kind in ['hann', 'dpss', 'rect']:
        w, P = spectra(psi, ds, kind); aw = np.abs(w)[:, None]
        nm = ((aw*P)*(w < 0)[:, None]).sum(axis=0); npl = ((aw*P)*(w > 0)[:, None]).sum(axis=0)
        ch = (w[:, None]*P).sum(axis=0)
        g = dict(nminus_nonzero=float(nm[nz].sum()/Q), nminus_zero=float(nm[~nz].sum()/Q),
                 nplus_nonzero=float(npl[nz].sum()/Q), nplus_zero=float(npl[~nz].sum()/Q),
                 charge_spectral=float(ch.sum()/Q))
        g['nminus_all'] = g['nminus_nonzero']+g['nminus_zero']
        # shell groups (nonzero modes)
        for name, sel in [('k_le_1p6', (kk <= 1.6) & nz), ('k_1p6_4', (kk > 1.6) & (kk < 4)), ('k_ge_4', kk >= 4)]:
            g['nminus_'+name] = float(nm[sel].sum()/Q)
        g['nminus_shell1'] = float(nm[shell == 1].sum()/Q); g['nminus_shell2'] = float(nm[shell == 2].sum()/Q)
        out[kind] = g
        if kind == 'dpss':
            # centroid frequencies per shell vs time-averaged Hartree
            X2 = (np.abs(psi)**2).sum(axis=1); mH2 = 2*LX*X2+G2*phi2; mH2bar = float(mH2.mean())
            rows = []; off = 0.; tot = float(nm[nz].sum())
            for n in range(1, shell.max()+1):
                sel = (shell == n) & nz
                if not sel.any(): continue
                Ps = P[:, sel].sum(axis=1)
                neg = w < 0; pos = w > 0
                wc_m = float((np.abs(w)*Ps)[neg].sum()/max(Ps[neg].sum(), 1e-300))
                wc_p = float((np.abs(w)*Ps)[pos].sum()/max(Ps[pos].sum(), 1e-300))
                ip = int(np.argmax(np.where(pos, Ps, -1))); im = int(np.argmax(np.where(neg, Ps, -1)))
                kf = float(np.sqrt(k2[sel].mean())); wH = float(np.sqrt(kf**2+mH2bar))
                nmn = float(nm[sel].sum())
                if abs(wc_m/wH-1) > .10: off += nmn
                rows.append(dict(shell=n, k=float(2*np.pi/L*n), k_fd4=kf, modes=int(sel.sum()), nminus=nmn/Q,
                                 centroid_neg=wc_m, centroid_pos=wc_p, peak_neg=float(-w[im]), peak_pos=float(w[ip]),
                                 omega_hartree=wH, m2_from_centroid_neg=wc_m**2-kf**2, m2_from_peak_neg=float(w[im]**2-kf**2)))
            out['dpss_shells'] = rows; out['mH2_mean'] = mH2bar
            out['nminus_fraction_in_shells_off_hartree_10pct'] = off/tot
    # C: per-mode virial (time averages), nonzero modes and zero mode
    A2 = (np.abs(psi)**2).mean(axis=0); B2 = (np.abs(ppsi)**2).mean(axis=0)
    q = np.imag(np.conj(psi)*ppsi).mean(axis=0)
    nv = .5*(np.sqrt(A2*B2)-q)
    out['virial'] = dict(nminus_nonzero=float(nv[nz].sum()/Q), nminus_zero=float(nv[~nz].sum()/Q),
                         nminus_all=float(nv.sum()/Q))
    return out


def main(path, out, sub=None, step=None):
    d = np.load(path); s = d['s']; N = int(d['N']); L = float(d['L'])
    psi = d['psi'].reshape(len(s), -1).astype(np.complex128); ppsi = d['ppsi'].reshape(len(s), -1).astype(np.complex128)
    a0 = np.arange(N//2+1); wr = np.where((a0 == 0) | (a0 == N//2), 1., 2.)
    phi2 = np.array([(wr[None, None, :]*np.abs(x)**2).sum() for x in d['phi']])
    M = len(s); sub = int(sub) if sub else M; step = int(step) if step else sub
    rows = []
    for i0 in range(0, M-sub+1, step):
        sl = slice(i0, i0+sub)
        r = analyse(s[sl], psi[sl], ppsi[sl], phi2[sl], N, L); rows.append(r)
        print('s %.2f-%.2f M=%d' % (r['s'][0], r['s'][1], r['samples']),
              ' '.join('%s: nz %.4f zero %.4f all %.4f' % (k, r[k]['nminus_nonzero'], r[k]['nminus_zero'], r[k]['nminus_all'])
                       for k in ['hann', 'dpss', 'rect', 'virial']), flush=True)
        g = r['dpss']
        print('   dpss groups: k<=1.6 %.4f  1.6-4 %.4f  k>=4 %.4f  shell1 %.4f shell2 %.4f  off-Hartree(>10%%) frac %.3f  mH2 %.2f' % (
            g['nminus_k_le_1p6'], g['nminus_k_1p6_4'], g['nminus_k_ge_4'], g['nminus_shell1'], g['nminus_shell2'],
            r['nminus_fraction_in_shells_off_hartree_10pct'], r['mH2_mean']), flush=True)
    json.dump(dict(source='en_particles_skeptic own code on reproduction-type window ' + path.split('/')[-2], rows=rows),
              open(out, 'w'), indent=1)
    return rows


if __name__ == '__main__':
    main(*sys.argv[1:])
