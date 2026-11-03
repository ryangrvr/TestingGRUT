"""lf1d: axially coupled Majorana singlets on a 1+1D line, exact free streaming plus local mass mixing, Fermi statistics
kept through the evolved Dirac sea (Claude, lf_fermion track; GitHub issue #13 context).

Model (mf_hybrid / e2e_sm couplings): -(M/2) Nbar exp(i gamma5 phi/f) N, studied in the UNROTATED basis: the chiral
components R, L of one helicity sector along the line are free massless (R moves at +1, L at -1 in comoving units and
conformal time eta) and are mixed locally by the complex mass a M exp(i theta(s, eta)), theta = phi/f.  Per comoving
momentum k the homogeneous problem is H(k) = [[k, c], [c*, -k]], c = a M e^{i theta}; rotating the phase away gives
mf_hybrid's two-level form [[p, M], [M, -p]] with p = k - mu5, mu5 = theta'/2 (so k > 0 and k < 0 on the line are the
two helicities of the 3D radial problem).
  * Free part: exact in k space (spectral derivative on a periodic grid, so no fermion doublers: the 'SLAC' dispersion
    E = k up to k_Nyq; the time-step and Nyquist conditions below keep every mode below the aliasing threshold).
  * Mass part: exact 2x2 exponential per site with theta at the step midpoint (Strang splitting).
  * Fermi statistics: the many-body state is the Slater determinant of the evolved negative-energy states.  We evolve
    n_r stochastic superpositions psi = sum_k xi_k v_k (|xi_k| = 1 random phases): E[psi psi^dagger] is the exact
    one-body density matrix C (low-cost-fermion idea of Borsanyi and Hindmarsh, arXiv:0809.4711).  For a homogeneous
    background C is diagonal in k, so the estimator is exact with n_r = 1; in an inhomogeneous background the noise of
    the occupation per mode is O(n_k), and averages down over modes and realisations.
  * Occupations: n(k) = u_k^dagger C(k, k) u_k with the positive-energy spinor u_k of H(k) for a reference phase.
  * 3D weighting (ray reduction): a homogeneous 3D problem decouples into these two-level systems, so its energy density
    per Majorana species is eps_3D = (1/2 pi^2) int_{-inf}^{inf} k^2 E_k n(k) dk with n(k) from the line (k > 0: one
    helicity, k < 0: the other).  For an inhomogeneous background the line is a ray sample [A, ray reduction].
Every choice (basis, splitting, estimator, reference phase, ray reduction) is a model choice.
"""
import math
import numpy as np


def spinors(k, c):
    """negative (v) and positive (u) energy eigenvectors of [[k, c], [c*, -k]] (c complex array or scalar)."""
    c = np.broadcast_to(c, k.shape)
    w = np.sqrt(k * k + np.abs(c) ** 2)
    nv = np.sqrt(np.abs(c) ** 2 + (w + k) ** 2)
    nu = np.sqrt(np.abs(c) ** 2 + (w - k) ** 2)
    # for c -> 0 the formulas degenerate on one side; regularise
    v = np.stack([c / np.maximum(nv, 1e-300), -(w + k) / np.maximum(nv, 1e-300)])
    u = np.stack([c / np.maximum(nu, 1e-300), (w - k) / np.maximum(nu, 1e-300)])
    small = np.abs(c) < 1e-14 * np.maximum(np.abs(k), 1)
    if np.any(small):
        pos = k > 0
        u[:, small & pos] = np.array([[1.0], [0.0]])
        u[:, small & ~pos] = np.array([[0.0], [1.0]])
        v[:, small & pos] = np.array([[0.0], [-1.0]])
        v[:, small & ~pos] = np.array([[1.0], [0.0]])
    return v, u, w


class Line:
    """n_s species (masses Ms, same theta) times n_r realisations on one periodic line of N sites, length L."""

    def __init__(self, N, L, Ms, n_r=1, seed=0, theta0=0.0, a0=1.0, dtype=np.complex128, phases='random'):
        self.N, self.L = N, L
        self.k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
        self.dk = 2 * np.pi / L
        self.Ms = np.asarray(Ms, dtype=float)
        self.ns, self.nr = len(Ms), n_r
        self.dtype = dtype
        rng = np.random.default_rng(seed)
        psi = np.zeros((self.ns, n_r, 2, N), dtype=complex)
        th0 = np.mean(theta0) if np.ndim(theta0) else theta0
        for j, M in enumerate(self.Ms):
            v, _, _ = spinors(self.k, a0 * M * np.exp(1j * th0))
            for r in range(n_r):
                xi = np.exp(2j * np.pi * rng.random(N)) if phases == 'random' else np.ones(N)
                psi[j, r] = v * xi
        # if theta0 is a profile, the sea is defined for the mean phase and then the profile is switched on
        # adiabatically by the caller (theta0 profile only used at t0 where gradients are tiny)
        self.psik = psi.astype(dtype)                     # state in k space

    def free(self, deta):
        ph = np.exp(-1j * self.k * deta).astype(self.dtype)
        self.psik[:, :, 0, :] *= ph
        self.psik[:, :, 1, :] *= np.conj(ph)

    def mass(self, deta, a, theta):
        """local kick exp(-i deta [[0, c], [c*, 0]]), c = a M e^{i theta(x)}; theta scalar or array(N)."""
        x = np.fft.ifft(self.psik, axis=-1)
        e = np.exp(1j * np.asarray(theta))
        for j, M in enumerate(self.Ms):
            g = a * M * deta
            cs, sn = math.cos(g), math.sin(g)
            R, Lc = x[j, :, 0, :], x[j, :, 1, :]
            Rn = cs * R - 1j * sn * e * Lc
            Ln = cs * Lc - 1j * sn * np.conj(e) * R
            x[j, :, 0, :], x[j, :, 1, :] = Rn, Ln
        self.psik = np.fft.fft(x, axis=-1).astype(self.dtype)

    def step(self, eta, deta, a_mid, theta_mid):
        self.free(0.5 * deta)
        self.mass(deta, a_mid, theta_mid)
        self.free(0.5 * deta)

    def occupation(self, a, theta_ref):
        """n(k) per species (mean over realisations), against the positive-energy spinor of H(k) at phase theta_ref.
        The k-space state is normalised so that a filled mode has |psi_k|^2 = 1."""
        out = []
        for j, M in enumerate(self.Ms):
            _, u, w = spinors(self.k, a * M * np.exp(1j * theta_ref))
            amp = np.conj(u[0]) * self.psik[j, :, 0, :] + np.conj(u[1]) * self.psik[j, :, 1, :]
            out.append(np.mean(np.abs(amp) ** 2, axis=0))
        return np.array(out)

    def eps3d(self, n, a, kmin_over_M=0.0):
        """3D-weighted energy density per Majorana species (comoving k, physical energy density) from n(k)."""
        out = []
        for j, M in enumerate(self.Ms):
            w = np.sqrt(self.k ** 2 + (a * M) ** 2)
            msk = np.abs(self.k) >= kmin_over_M * a * M
            out.append(np.sum((self.k ** 2 * w * n[j])[msk]) * self.dk / (2 * np.pi ** 2) / a ** 4)
        return np.array(out)

    def local_occ(self, a, theta_line, width, centres, kbins):
        """windowed (Gabor) occupation n(s, k) for |k| above the mass scale, chiral measure: R occupation for k > 0,
        L occupation for k < 0 (exact for |k| >> aM).  Returns array (ns, ncent, nbins) of mean n in each |k| bin, and
        the 3D-weighted energy density per window."""
        x = np.fft.ifft(self.psik, axis=-1)                     # site amplitudes (normalised: filled mode -> 1/N)
        s = np.arange(self.N) * self.L / self.N
        res_n, res_e = [], []
        for j, M in enumerate(self.Ms):
            nj, ej = [], []
            for c0 in centres:
                d = (s - c0 + 0.5 * self.L) % self.L - 0.5 * self.L
                win = np.exp(-0.5 * (d / width) ** 2)
                win /= math.sqrt(np.sum(win ** 2) / self.N)     # unit mean square
                Rk = np.fft.fft(x[j, :, 0, :] * win, axis=-1)
                Lk = np.fft.fft(x[j, :, 1, :] * win, axis=-1)
                nR = np.mean(np.abs(Rk) ** 2, axis=0)
                nL = np.mean(np.abs(Lk) ** 2, axis=0)
                n = np.where(self.k > 0, nR, nL)
                w = np.abs(self.k)
                e = np.sum((self.k ** 2 * w * n)[w > 3 * a * M]) * self.dk / (2 * np.pi ** 2) / a ** 4
                hist = [float(np.mean(n[(w >= lo) & (w < hi)])) if np.any((w >= lo) & (w < hi)) else 0.0
                        for lo, hi in zip(kbins[:-1], kbins[1:])]
                nj.append(hist)
                ej.append(e)
            res_n.append(nj)
            res_e.append(ej)
        return np.array(res_n), np.array(res_e)
