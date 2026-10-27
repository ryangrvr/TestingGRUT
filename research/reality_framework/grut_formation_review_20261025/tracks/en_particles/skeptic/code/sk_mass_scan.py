"""en_particles SKEPTIC (own code): how strongly does the equal-time n_- reference depend on the assumed mass?

For a window of the regenerated history, computes the Hann^2-weighted time average of the equal-time count
    n_-(m^2) = sum_{k != 0} |w psi_k + i pi_k|^2/(4w),  w = sqrt(k_FD4^2 + m^2),  m^2 fixed,
for a grid of m^2, plus (i) the instantaneous isotropic Hartree mass m_H^2(s) and (ii) the per-mode time-averaged
virial count 1/2[sqrt(<|pi|^2><|psi|^2>) - <q>] with the same weights. If n_-(m^2) is flat over a broad m^2 range,
agreement of a measured count with the Hartree reference does not select the Hartree dispersion.
Usage: python3 -I sk_mass_scan.py <window.npz> <out.json> [i0] [i1]
"""
import sys, json
import numpy as np

LX, G2 = 1e-3, 1e-2


def k2_fd4(N, L):
    th = 2*np.pi*np.fft.fftfreq(N); dx = L/N
    one = ((4/3)*(2-2*np.cos(th))-(2-2*np.cos(2*th))/12)/dx**2
    return (one[:, None, None]+one[None, :, None]+one[None, None, :]).ravel()


def main(path, out, i0=0, i1=None):
    d = np.load(path); s = d['s']; i0 = int(i0); i1 = int(i1) if i1 else len(s); s = s[i0:i1]; M = len(s)
    N = int(d['N']); L = float(d['L'])
    psi = d['psi'][i0:i1].reshape(M, -1); ppsi = d['ppsi'][i0:i1].reshape(M, -1)
    a0 = np.arange(N//2+1); wr = np.where((a0 == 0) | (a0 == N//2), 1., 2.)
    phi2 = np.array([(wr[None, None, :]*np.abs(x)**2).sum() for x in d['phi'][i0:i1]])
    k2 = k2_fd4(N, L); nz = k2 > 1e-12; k2n = k2[nz]
    h = np.hanning(M); wt = h**2/(h**2).sum()
    grid = [0., 2., 5., 8., 10., 12., 14., 16., 18., 20., 22., 25., 30., 40., 60.]
    acc = np.zeros(len(grid)); accH = 0.; A2 = 0.; B2 = 0.; qa = 0.; Qs = []; mHs = []
    for j in range(M):
        a = psi[j].astype(np.complex128); b = ppsi[j].astype(np.complex128)
        Q = float(np.imag(np.conj(a)*b).sum()); Qs.append(Q)
        an = a[nz]; bn = b[nz]
        for g, m2 in enumerate(grid):
            w = np.sqrt(k2n+m2); acc[g] += wt[j]*float((np.abs(w*an+1j*bn)**2/(4*w)).sum())/Q
        mH2 = 2*LX*float((np.abs(a)**2).sum())+G2*phi2[j]; mHs.append(mH2)
        w = np.sqrt(k2n+mH2); accH += wt[j]*float((np.abs(w*an+1j*bn)**2/(4*w)).sum())/Q
        A2 = A2+wt[j]*np.abs(an)**2; B2 = B2+wt[j]*np.abs(bn)**2; qa = qa+wt[j]*np.imag(np.conj(an)*bn)
    Q = float(np.mean(Qs))
    vir = float((.5*(np.sqrt(A2*B2)-qa)).sum())/Q
    res = dict(source='en_particles_skeptic own code on reproduction-type window', window=[float(s[0]), float(s[-1])],
               samples=M, weights='Hann^2 (same as en_particles equal-time averages)',
               scan=[dict(m2=m2, nminus_over_Q=float(v)) for m2, v in zip(grid, acc)],
               instantaneous_hartree=dict(nminus_over_Q=accH, mH2_weighted_mean=float((wt*np.array(mHs)).sum()),
                                          mH2_min=float(min(mHs)), mH2_max=float(max(mHs))),
               virial_hann2=vir)
    json.dump(res, open(out, 'w'), indent=1)
    print('window s %.2f-%.2f' % (s[0], s[-1]), 'virial %.4f' % vir, 'instHartree %.4f' % accH,
          'mH2 mean %.2f [%.2f, %.2f]' % (res['instantaneous_hartree']['mH2_weighted_mean'], min(mHs), max(mHs)))
    for r in res['scan']: print('  m2 %5.1f  n-/Q %.4f' % (r['m2'], r['nminus_over_Q']))


if __name__ == '__main__':
    main(*sys.argv[1:])
