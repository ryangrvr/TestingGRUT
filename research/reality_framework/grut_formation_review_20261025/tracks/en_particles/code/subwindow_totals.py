"""en_particles: sliding sub-window totals of the measured-frequency counts (own code) [E].

For each sub-window (Hann taper) computes, per net charge, the lab-frame |w|-weighted spectral counts
n+-(zero mode) and n+-(nonzero modes), and the equal-time Hartree/massless references averaged with the same taper.
Usage: python3 -I subwindow_totals.py <out.json> <window.npz> [<window.npz> ...]  (env: SUB=100 STEP=20)
"""
import sys, os, json
import numpy as np

LX, G2 = 1e-3, 1e-2


def fd4(N, L):
    a = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    f = lambda x: 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    return f(a)[:, None, None]+f(a)[None, :, None]+f(a)[None, None, :]


def main(out, *paths):
    sub = int(os.environ.get('SUB', 100)); step = int(os.environ.get('STEP', 20)); rows = []
    for p in paths:
        d = np.load(p); s = d['s']; N = int(d['N']); L = float(d['L']); k2 = fd4(N, L).ravel(); nz = k2 > 0
        psi = d['psi'].reshape(len(s), -1); ppsi = d['ppsi'].reshape(len(s), -1); ds = s[1]-s[0]
        nQ = float(np.imag(np.conj(psi[0].astype(np.complex128))*ppsi[0].astype(np.complex128)).sum())  # conserved
        a0 = np.arange(N//2+1); wr = np.where((a0 == 0) | (a0 == N//2), 1., 2.)
        phi2 = np.array([(wr[None, None, :]*np.abs(x)**2).sum() for x in d['phi']])
        a0c, C = float(os.environ.get('A0', 8.869218738806)), float(os.environ.get('CC', 1.9857962750634526))
        for i0 in range(0, len(s)-sub+1, step):
            sl = slice(i0, i0+sub); M = sub; h = np.hanning(M); h2 = (h**2).sum(); wgt = h**2/h2
            X = np.fft.fft(psi[sl].astype(np.complex128)*h[:, None], axis=0)
            P = np.abs(X)**2/(M*h2); w = 2*np.pi*np.fft.fftfreq(M, ds); aw = np.abs(w)[:, None]
            neg = (w < 0)[:, None]; pos = (w > 0)[:, None]
            nm_nz = float((aw*P*neg)[:, nz].sum())/nQ; np_nz = float((aw*P*pos)[:, nz].sum())/nQ
            nm_0 = float((aw*P*neg)[:, ~nz].sum())/nQ; np_0 = float((aw*P*pos)[:, ~nz].sum())/nQ
            # equal-time references (same taper)
            refs = {'massless': 0., 'hartree': 0.}
            for j, wj in zip(range(i0, i0+sub), wgt):
                A = psi[j].astype(np.complex128); B = ppsi[j].astype(np.complex128)
                X2 = float((np.abs(A)**2).sum()); m2 = 2*LX*X2+G2*phi2[j]; Qk = np.imag(np.conj(A)*B)
                for name, mm in [('massless', 0.), ('hartree', m2)]:
                    om = np.sqrt(np.where(nz, k2+mm, 1.))
                    refs[name] += wj*float((.5*((np.abs(B)**2/om+om*np.abs(A)**2)/2-Qk/2))[nz].sum())/nQ
            rows.append(dict(window=p.split('/')[-2], N=N, L=L, nQ=nQ, s_start=float(s[i0]), s_end=float(s[i0+sub-1]),
                             tau_start=float(20+a0c*s[i0]+C*s[i0]**2/2),
                             tau_end=float(20+a0c*s[i0+sub-1]+C*s[i0+sub-1]**2/2),
                             n_minus_nonzero=nm_nz, n_plus_nonzero=np_nz, n_minus_zero=nm_0, n_plus_zero=np_0,
                             n_minus_all=nm_nz+nm_0, n_plus_all=np_nz+np_0, charge_all=np_nz+np_0-nm_nz-nm_0,
                             eq_massless=refs['massless'], eq_hartree=refs['hartree']))
            r = rows[-1]
            print(f"{r['window']:8s} s {r['s_start']:6.2f}-{r['s_end']:6.2f} tau {r['tau_start']:7.1f}-{r['tau_end']:7.1f} "
                  f"n-nz {nm_nz:.4f} n-0 {nm_0:.4f} n-all {r['n_minus_all']:.4f} n+all {r['n_plus_all']:.4f} Q {r['charge_all']:.5f} "
                  f"| eqH {refs['hartree']:.4f} eqM {refs['massless']:.4f}", flush=True)
    json.dump(dict(source='own code; Hann sub-windows of the regenerated N48 history', sub_samples=sub, step=step,
                   rows=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
