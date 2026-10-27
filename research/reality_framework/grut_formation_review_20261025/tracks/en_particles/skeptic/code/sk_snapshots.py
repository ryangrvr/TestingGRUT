"""en_particles SKEPTIC (own code): per-snapshot reference counts with an amplitude formulation.

Different from en_particles/snapshot_counts.py (full FFT of psi and N_k(w) +- Q_k) and from Codex analyze.py
(rfft, chosen-omega quadratic form): here each real field chi1, chi2, pi1, pi2 is rfft'd separately, the complex
mode amplitudes are rebuilt per half-space mode with BOTH +k and -k partners explicitly, and
    n_-(k) = |w psi_k + i pi_k|^2/(4w),   n_+(k) = |w psi_k - i pi_k|^2/(4w)
is summed. (psi_k ~ B e^{-iws} is the antiparticle branch; derivation in SKEPTIC.md.)
Also reports Q0/Q, the energy fraction outside the spatial mean and the Hartree mass for each snapshot.
Usage: python3 -I sk_snapshots.py <out.json> <state.bin> [...]
"""
import sys, json, struct
import numpy as np

LX, G2, LP = 1e-3, 1e-2, 3e-2


def load(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    return N, L, s, np.fromfile(path, dtype='<f8', offset=24).reshape(6, N, N, N)


def lattice_k2(N, L):
    # FD4 symbol for f''  ~ (4/3)(f+ + f- - 2f)/dx^2 - (f++ + f-- - 2f)/(12 dx^2), written directly from the stencil
    th = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    one = ((4/3)*(2-2*np.cos(th))-(2-2*np.cos(2*th))/12)/dx**2
    return one[:, None, None]+one[None, :, None]+one[None, None, :]


def amp_counts(q, p, N, L):
    # full complex psi = chi1 + i chi2 built from separate real FFTs (numpy fftn on real input, forward norm)
    F = lambda x: np.fft.fftn(x, norm='forward')
    c1, c2, d1, d2 = F(q[0]), F(q[1]), F(p[0]), F(p[1])
    psi = c1+1j*c2; pi = d1+1j*d2
    k2 = lattice_k2(N, L); nz = k2 > 1e-12
    X2 = float(np.mean(q[0]**2+q[1]**2)); f2 = float(np.mean(q[2]**2)); mH2 = 2*LX*X2+G2*f2
    Q = float(np.mean(q[0]*p[1]-q[1]*p[0]))
    out = dict(Q=Q, mH2=mH2)
    for name, m2 in [('massless', 0.), ('hartree', mH2)]:
        w = np.sqrt(np.where(nz, k2+m2, 1.))
        nm = np.abs(w*psi+1j*pi)**2/(4*w); npl = np.abs(w*psi-1j*pi)**2/(4*w)
        out[name] = dict(nminus_over_Q=float(nm[nz].sum()/Q), nplus_over_Q=float(npl[nz].sum()/Q),
                         charge_identity=float((npl[nz].sum()-nm[nz].sum())/Q))
    mq = q.reshape(3, -1).mean(1); mp = p.reshape(3, -1).mean(1)
    Q0 = mq[0]*mp[1]-mq[1]*mp[0]
    r0 = mq[0]**2+mq[1]**2
    E0 = .5*np.sum(mp**2)+LX/4*r0**2+G2/2*mq[2]**2*r0+LP/4*mq[2]**4
    # total energy with the same FD4 operator (spectral)
    lapE = 0.
    for c in range(3):
        fc = F(q[c]); lapE += .5*float((k2*np.abs(fc)**2).sum())
    r2 = q[0]**2+q[1]**2
    E = .5*float(np.mean(np.sum(p*p, axis=0)))+lapE+float(np.mean(LX/4*r2*r2+G2/2*q[2]**2*r2+LP/4*q[2]**4))
    out.update(Q0_over_Q=float(Q0/Q), E=E, outside_mean_fraction=float(1-E0/E))
    return out


def main(out, *paths):
    rows = []
    for pth in paths:
        N, L, s, a = load(pth)
        r = amp_counts(a[:3], a[3:], N, L); r.update(path=pth.split('/')[-1], s=s, N=N, L=L)
        rows.append(r)
        print(r['path'], 's=%.5f' % s, 'Q=%.6f' % r['Q'], 'massless n-/Q=%.6f' % r['massless']['nminus_over_Q'],
              'hartree n-/Q=%.6f' % r['hartree']['nminus_over_Q'], 'mH2=%.4f' % r['mH2'],
              'Q0/Q=%.6f' % r['Q0_over_Q'], 'outside=%.5f' % r['outside_mean_fraction'],
              'ident=%.1e' % (r['massless']['charge_identity']-1+r['Q0_over_Q']), flush=True)
    json.dump(dict(source='en_particles_skeptic own code (amplitude formulation) on reproduction-type snapshots', rows=rows),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
