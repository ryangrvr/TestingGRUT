"""e2e_lattice post-processing: percolation of the -v vacuum, field spectra, TT-projected GW spectra and their
redshift to today, summary of one run.  Labels: [C] on the lattice model; the redshift uses the run's own K
(from the fluid continuation) and standard entropy conservation (g_* = 106.75, g_s0 = 3.931).
"""
import math
import numpy as np
from scipy import ndimage
import latcommon as LC

GEV_TO_HZ = 1.0 / 6.582119569e-25          # 1 GeV = 1.519e24 s^-1
T0_GEV = 2.7255 * 8.617333e-14
GS0 = 3.931
GSTAR = 106.75
OMEGA_GAMMA_H2 = 2.473e-5                  # photons today (T0 = 2.7255 K)


def clusters(mask):
    """Clusters of a boolean field in the periodic box (6-connectivity).  Open-boundary labels are joined across the
    periodic faces by a weighted union-find that carries the box shift between label frames; a periodic cluster
    'wraps' (percolates on the torus) when two joins disagree on the shift.  Returns the set fraction, the largest
    periodic cluster's volume fraction, the number of periodic clusters, the wrapping volume fraction, and wraps."""
    if not mask.any():
        return dict(frac=0.0, largest=0.0, n=0, wrap_frac=0.0, wraps=False)
    lab, n = ndimage.label(mask)
    N = mask.shape[0]
    parent = list(range(n + 1))
    off = [np.zeros(3, dtype=int) for _ in range(n + 1)]
    wraps_root = set()

    def find(a):
        d = np.zeros(3, dtype=int)
        path = []
        while parent[a] != a:
            path.append(a)
            d = d + off[a]
            a = parent[a]
        # path compression with offsets
        acc = d.copy()
        for b in path:
            ob = off[b].copy()
            off[b] = acc.copy()
            parent[b] = a
            acc = acc - ob
        return a, d
    for ax in range(3):
        a0 = np.take(lab, 0, axis=ax).ravel(); a1 = np.take(lab, N - 1, axis=ax).ravel()
        sel = (a0 > 0) & (a1 > 0)
        pairs = np.unique(np.stack([a1[sel], a0[sel]], axis=1), axis=0) if sel.any() else []
        e = np.zeros(3, dtype=int); e[ax] = 1
        for p_, q_ in pairs:                        # p_ at face N-1, q_ at face 0: frame(q) = frame(p) + e
            rp, dp = find(int(p_)); rq, dq = find(int(q_))
            if rp == rq:
                if not np.array_equal(dq - dp, e):
                    wraps_root.add(rp)
            else:
                parent[rq] = rp
                off[rq] = dp + e - dq
                if rq in wraps_root:
                    wraps_root.add(rp)
    roots = np.array([find(i)[0] for i in range(n + 1)])
    wr = set(find(r)[0] for r in wraps_root)
    so = np.bincount(lab.ravel()); so[0] = 0
    per = np.bincount(roots[1:], weights=so[1:], minlength=n + 1)
    big = int(np.argmax(per))
    wfrac = float(sum(per[r] for r in wr) / mask.size)
    return dict(frac=float(mask.mean()), largest=float(per[big] / mask.size), n=int(len(set(roots[1:].tolist()))),
                wrap_frac=wfrac, wraps=bool(len(wr) > 0))


def power(field, L):
    """Shell-averaged power: d<f^2>/d ln k = k^3 P(k)/(2 pi^2), with <f^2> = sum over shells."""
    N = field.shape[0]
    F = np.fft.fftn(field - field.mean())
    kk = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing='ij')
    K = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    kf = 2 * np.pi / L
    nb = int(N * math.sqrt(3) / 2) + 1
    b = np.floor(K / kf + 0.5).astype(int)
    p = np.abs(F) ** 2 / N ** 6
    s = np.bincount(b.ravel(), weights=p.ravel(), minlength=nb)
    c = np.bincount(b.ravel(), minlength=nb)
    k = np.arange(nb) * kf
    # per unit ln k: shell sum / (dk/k) with dk = kf
    with np.errstate(divide='ignore', invalid='ignore'):
        dlnk = np.where(k > 0, kf / k, np.nan)
    return k[1:], (s / dlnk)[1:], c[1:]


def gw_spectrum(hdot, L, q):
    """TT projection with the central-difference lattice momentum k_eff,i = sin(k_i dx)/dx (cf. Figueroa et al.
    1110.0337), then d rho_GW / d ln k in rho~ = rho/v^2 units: (1/4q) <hdot_ij hdot_ij>^TT per ln k."""
    N = hdot[0].shape[0]
    dx = L / N
    F = [np.fft.fftn(h.astype(np.float64)) for h in hdot]        # 0 xx 1 yy 2 zz 3 xy 4 xz 5 yz
    kk = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    ke = np.sin(kk * dx) / dx
    KX, KY, KZ = np.meshgrid(ke, ke, ke, indexing='ij')
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        nx, ny, nz = [np.where(K2 > 0, X / np.sqrt(K2), 0.0) for X in (KX, KY, KZ)]
    n = [nx, ny, nz]
    idx = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3, (1, 0): 3, (0, 2): 4, (2, 0): 4, (1, 2): 5, (2, 1): 5}
    Pm = [[(1.0 if i == j else 0.0) - n[i] * n[j] for j in range(3)] for i in range(3)]
    tot = np.zeros(F[0].shape)
    for i in range(3):
        for j in range(3):
            acc = np.zeros(F[0].shape, dtype=complex)
            for l in range(3):
                for m in range(3):
                    lam = Pm[i][l] * Pm[j][m] - 0.5 * Pm[i][j] * Pm[l][m]
                    acc += lam * F[idx[(l, m)]]
            tot += np.abs(acc) ** 2
    tot[K2 == 0] = 0.0
    kc = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KXc, KYc, KZc = np.meshgrid(kc, kc, kc, indexing='ij')
    K = np.sqrt(KXc ** 2 + KYc ** 2 + KZc ** 2)
    kf = 2 * np.pi / L
    nb = int(N * math.sqrt(3) / 2) + 1
    b = np.floor(K / kf + 0.5).astype(int)
    s = np.bincount(b.ravel(), weights=tot.ravel(), minlength=nb) / N ** 6 / (4 * q)
    k = np.arange(nb) * kf
    with np.errstate(divide='ignore', invalid='ignore'):
        spec = np.where(k > 0, s * k / kf, 0.0)
    # non-projected trace for comparison
    raw = sum((1.0 if c < 3 else 2.0) * np.mean(h.astype(np.float64) ** 2) for c, h in enumerate(hdot)) / (4 * q)
    return k[1:], spec[1:], float(s.sum()), float(raw)


def today(k_com, Omega_rad_era, H_GeV, K, g=GSTAR):
    """Comoving k (H_d units, a = 1 at the slice) -> frequency today [Hz]; radiation-era fraction -> Omega_GW h^2."""
    MP = LC.MP_GEV
    Ta = (90.0 * K * K * MP ** 2 * H_GeV ** 2 / (math.pi ** 2 * g)) ** 0.25
    a0 = (g / GS0) ** (1.0 / 3.0) * Ta / T0_GEV
    f = np.asarray(k_com) * H_GeV * GEV_TO_HZ / (2 * math.pi * a0)
    fac = OMEGA_GAMMA_H2 * (g / 2.0) * (GS0 / g) ** (4.0 / 3.0)
    return f, fac * np.asarray(Omega_rad_era), dict(a0=a0, Ta=Ta, factor=fac)
