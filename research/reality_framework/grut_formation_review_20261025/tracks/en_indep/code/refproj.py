"""Own reference projections n_-(k) (massless / isotropic Hartree), MATHEMATICS.md formula, full complex FFT."""
import numpy as np
import indep3d as I

def refproj(q, p, L=12.0):
    N = q.shape[-1]; dx = L / N
    th = 2 * np.pi * np.fft.fftfreq(N); s2 = np.sin(th / 2) ** 2; f = 4 * s2 * (1 + s2 / 3) / dx ** 2
    K2 = f[:, None, None] + f[None, :, None] + f[None, None, :]
    qh = np.fft.fftn(q, axes=(1, 2, 3)) / N ** 3; ph = np.fft.fftn(p, axes=(1, 2, 3)) / N ** 3
    nz = K2 > 0
    Qtot = float(np.mean(q[0] * p[1] - q[1] * p[0]))
    Qk = (np.conj(qh[0]) * ph[1] - np.conj(qh[1]) * ph[0]).real
    r = q[0] ** 2 + q[1] ** 2
    o = {}
    for name, m2 in (('massless', 0.0), ('isotropic_Hartree', 2 * I.LX * float(r.mean()) + I.G2 * float(np.mean(q[2] ** 2)))):
        w = np.sqrt(K2[nz] + m2)
        nr = (abs(ph[0][nz]) ** 2 + abs(ph[1][nz]) ** 2 + w ** 2 * (abs(qh[0][nz]) ** 2 + abs(qh[1][nz]) ** 2)) / (2 * w)
        nplus = float(np.sum(nr + Qk[nz]) / 2); nminus = float(np.sum(nr - Qk[nz]) / 2)
        Q0 = float(q[0].mean() * p[1].mean() - q[1].mean() * p[0].mean())
        o[name] = dict(mass2=m2, nminus_over_net=nminus / Qtot, charge_identity=abs((nplus - nminus) - (Qtot - Q0)) / abs(Qtot))
    return o
