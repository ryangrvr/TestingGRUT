"""Free-mode (classical-to-quanta) spectral extraction for the spectrum track (Claude, 25 Oct 2026).

Units and normalisation (all dimensionless, Qbar units): U = r psi, Z = r(S-1),
  charge   Q = 4 pi int dr (-Im U* U_t),   charge current j_Q = 4 pi Im(U* U_r),
  energy current j_E = -4 pi [Re(U_t* U_r) + Z_t Z_r].
A free outgoing s-wave of the charged field has U(R,t) = int dw/2pi A(w) exp(-i w t),
with k(w) = sign(w) sqrt(w^2-1) (outgoing for both signs).  Then
  dN_X/dw    = 2 k |A(w)|^2           (w > 1;  particles, charge +1 each)
  dN_Xbar/dw = 2 |k| |A(-w)|^2        (w > 1;  antiparticles, charge -1 each)
  dE = |w| dN.   Physical numbers are N/g_X^2.
Real neutral field (mass^2 = 6): dN_s/dw = 4 k |Z(w)|^2 for w > sqrt 6.
Spatial projection at time T on the box sine modes s_n(r) (exact normal modes of the 6th-order
grid operator with odd continuation at r=0 and at the wall):
  alpha_n = (u_n + i v_n/w_n)/2, beta_n = (u_n - i v_n/w_n)/2,
  N_X(n) = 4 pi h w_n |alpha_n|^2,  N_Xbar(n) = 4 pi h w_n |beta_n|^2  (exact for the linear grid system).
"""
import numpy as np
from scipy.fft import dst

FOURPI = 4 * np.pi


def load_ext(prefix, ne):
    raw = np.fromfile(prefix + '_ext.bin')
    rec = 1 + 9 * ne
    m = len(raw) // rec
    raw = raw[:m * rec].reshape(m, rec)
    t = raw[:, 0]
    d = raw[:, 1:].reshape(m, ne, 9)
    return t, d


def load_hist(prefix, nin):
    raw = np.fromfile(prefix + '_hist.bin')
    rec = 3 + 2 * nin
    m = len(raw) // rec
    return raw[:m * rec].reshape(m, rec)


def load_snap(path):
    with open(path, 'rb') as fi:
        N = int(np.fromfile(fi, np.int64, 1)[0]); h, t = np.fromfile(fi, np.float64, 2)
        a = np.fromfile(fi, np.float64, 6 * (N + 1)).reshape(6, N + 1)
    return dict(N=N, h=h, t=t, x=a[0], y=a[1], z=a[2], vx=a[3], vy=a[4], vz=a[5])


def series(d, k):
    """Return U, Ut, Ur, Z, Zt, Zr at extraction radius index k."""
    s = d[:, k, :]
    return s[:, 0] + 1j * s[:, 1], s[:, 3] + 1j * s[:, 4], s[:, 6] + 1j * s[:, 7], s[:, 2], s[:, 5], s[:, 8]


def time_fluxes(t, U, Ut, Ur, Z, Zt, Zr):
    jQ = FOURPI * np.imag(np.conj(U) * Ur)
    jE = -FOURPI * (np.real(np.conj(Ut) * Ur) + Zt * Zr)
    jEs = -FOURPI * Zt * Zr
    cum = lambda j: np.concatenate([[0], np.cumsum(0.5 * (j[1:] + j[:-1]) * np.diff(t))])
    return dict(jQ=jQ, jE=jE, Qpass=cum(jQ), Epass=cum(jE), Es_pass=cum(jEs))


def taper(t, t0, t1, width):
    """1 on [t0+width, t1-width], cos^2 ramps at both ends, 0 outside."""
    w = np.ones_like(t)
    a = (t - t0) / width; b = (t1 - t) / width
    w = np.where(a < 1, np.sin(0.5 * np.pi * np.clip(a, 0, 1))**2, w)
    w = np.where(b < 1, w * np.sin(0.5 * np.pi * np.clip(b, 0, 1))**2, w)
    return w


def ft(t, x, pad=8):
    """A(w) = int dt x(t) e^{i w t} on a zero-padded grid; returns w (sorted) and A."""
    dt = t[1] - t[0]; M = len(t); Mp = int(2**np.ceil(np.log2(M * pad)))
    xp = np.zeros(Mp, complex); xp[:M] = x
    A = dt * Mp * np.fft.ifft(xp) * np.exp(1j * 2 * np.pi * np.fft.fftfreq(Mp, dt) * t[0])
    w = 2 * np.pi * np.fft.fftfreq(Mp, dt)
    o = np.argsort(w)
    return w[o], A[o]


def kofw(w, m2=1.0):
    """Outgoing branch: k = sign(w) sqrt(w^2-m2) for |w|>m; +i sqrt(m2-w^2) (decaying) below."""
    k = np.where(w * w > m2, np.sign(w) * np.sqrt(np.abs(w * w - m2)), 1j * np.sqrt(np.abs(m2 - w * w)))
    return k


def charged_spectrum(t, U, Ur, wt):
    """In/out-separated charged spectra from the tapered time series at one radius."""
    w, A = ft(t, wt * U); _, B = ft(t, wt * Ur)
    k = kofw(w)
    prop = np.abs(w) > 1
    aout = np.where(prop, (1j * k * A + B) / (2j * k + 1e-300), 0)
    ain = np.where(prop, (1j * k * A - B) / (2j * k + 1e-300), 0)
    kr = np.real(k)
    dNout = np.where(prop, 2 * np.abs(kr) * np.abs(aout)**2, 0)
    dNin = np.where(prop, 2 * np.abs(kr) * np.abs(ain)**2, 0)
    dw = w[1] - w[0]
    pos = w > 1; neg = w < -1; sub = ~prop
    # Parseval form of the tapered charge flux, split into its frequency bands
    fluxdens = 2 * np.imag(np.conj(A) * B)    # (4pi/2pi) Im(A* B)
    out = dict(w=w, A=A, aout=aout, ain=ain, dNout=dNout, dNin=dNin, dw=dw,
               NX=float(np.sum(dNout[pos]) * dw), NXbar=float(np.sum(dNout[neg]) * dw),
               NX_in=float(np.sum(dNin[pos]) * dw), NXbar_in=float(np.sum(dNin[neg]) * dw),
               EX=float(np.sum((np.abs(w) * dNout)[pos]) * dw), EXbar=float(np.sum((np.abs(w) * dNout)[neg]) * dw),
               flux_total=float(np.sum(fluxdens) * dw), flux_subthreshold=float(np.sum(fluxdens[sub]) * dw),
               power_subthreshold_fraction=float(np.sum(np.abs(A[sub])**2) / np.sum(np.abs(A)**2)))
    return out


def neutral_spectrum(t, Z, Zr, wt, m2=6.0):
    w, A = ft(t, wt * Z); _, B = ft(t, wt * Zr)
    k = kofw(w, m2); prop = (w > np.sqrt(m2))
    aout = np.where(prop, (1j * k * A + B) / (2j * k + 1e-300), 0)
    ain = np.where(prop, (1j * k * A - B) / (2j * k + 1e-300), 0)
    dw = w[1] - w[0]
    dN = np.where(prop, 4 * np.real(k) * np.abs(aout)**2, 0)
    dNin = np.where(prop, 4 * np.real(k) * np.abs(ain)**2, 0)
    fluxE = -2 * np.real(np.conj(-1j * w * A) * B)  # (4pi/2pi) * (-Re(Zt* Zr)) density, both signs of w
    return dict(w=w, dN=dN, Ns=float(np.sum(dN) * dw), Es=float(np.sum(w * dN) * dw), Ns_in=float(np.sum(dNin) * dw),
                Eflux_total=float(np.sum(fluxE) * dw), Eflux_subthreshold=float(np.sum(fluxE[np.abs(w) <= np.sqrt(m2)]) * dw))


def K2_6(k, h):
    return (49 / 18 - 3 * np.cos(k * h) + 0.3 * np.cos(2 * k * h) - (1 / 45) * np.cos(3 * k * h)) / h**2


def ramp(r, R0, width):
    s = np.clip((r - R0) / width, 0, 1)
    return np.sin(0.5 * np.pi * s)**2


def spatial_projection(snap, Rw, width, m2=1.0, field='charged'):
    N, h = snap['N'], snap['h']; L = N * h
    r = np.arange(1, N) * h; n = np.arange(1, N); k = n * np.pi / L
    w = np.sqrt(m2 + K2_6(k, h))
    win = ramp(r, Rw, width)
    T = lambda a: dst(a, type=1, norm='ortho')
    if field == 'charged':
        U = (snap['x'] + 1j * snap['y'])[1:N] * win; V = (snap['vx'] + 1j * snap['vy'])[1:N] * win
        u = T(U.real) + 1j * T(U.imag); v = T(V.real) + 1j * T(V.imag)
        a = 0.5 * (u + 1j * v / w); b = 0.5 * (u - 1j * v / w)
        NX = FOURPI * h * w * np.abs(a)**2; NXb = FOURPI * h * w * np.abs(b)**2
        Qdirect = FOURPI * h * np.sum(-np.imag(np.conj(U) * V))
        return dict(k=k, w=w, NXn=NX, NXbn=NXb, NX=float(NX.sum()), NXbar=float(NXb.sum()),
                    EX=float((w * NX).sum()), EXbar=float((w * NXb).sum()), Q_window_direct=float(Qdirect), dk=np.pi / L)
    Zf = snap['z'][1:N] * win; Zv = snap['vz'][1:N] * win
    u = T(Zf); v = T(Zv)
    En = FOURPI * h * 0.5 * (v**2 + w**2 * u**2)
    return dict(k=k, w=w, En=En, Nn=En / w, Es=float(En.sum()), Ns=float((En / w).sum()), dk=np.pi / L)


def propagate(t, aout_w, w, dR, M):
    """Free outgoing propagation of a measured outgoing spectrum by dR; returns U(t) on the first M samples."""
    k = kofw(w)
    A = aout_w * np.exp(1j * k * dR)
    A = np.where(np.isfinite(A), A, 0)
    # invert ft(): x(t_j) = (1/2pi) int dw A e^{-i w t}
    dt = t[1] - t[0]
    wu = np.fft.ifftshift(w); Au = np.fft.ifftshift(A)
    Mp = len(w)
    x = np.fft.fft(Au * np.exp(-1j * wu * t[0])) / (dt * Mp)
    return x[:M]
