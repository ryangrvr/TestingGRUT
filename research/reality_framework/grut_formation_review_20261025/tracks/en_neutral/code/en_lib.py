"""en_neutral: own helper library (Claude review of Codex 7d9c93b, issue #13).

INDEPENDENT code: own Fourier-realization generator written from MATHEMATICS.md / the documented seed key, own
real-space FD4 energy (np.roll stencils, not Codex's Fourier symbol), own fixed-step RK4 for the homogeneous reduced
orbit and the delta-phi transfer matrices, own Floquet monodromy.  Codex's published MATCH_L12.npz (git 7d9c93b) is
read as data only.

Reduced conformal Hamiltonian (Codex MATHEMATICS.md):
  H = 1/2 pi^2 + 1/2 |grad chi|^2 + lamX/4 (chi1^2+chi2^2)^2 + g2/2 chi_phi^2 (chi1^2+chi2^2) + lamP/4 chi_phi^4
"""
import numpy as np

LX, G2, LP = 1e-3, 1e-2, 3e-2


def fd4_symbol_1d(kx, dx):
    """Fourier symbol of the symmetric 5-point fourth-order second difference (own derivation:
    -(4/3)(2cos x - 2) + (1/12)(2cos 2x - 2) = 4 s^2 + 4/3 s^4 with s = sin(x/2))."""
    s = np.sin(kx * dx / 2)
    return (4 * s * s + 4 / 3 * s ** 4) / dx ** 2


def lattice_k2(n_vec, N, L):
    """FD4 effective k^2 for integer wavevectors n_vec (..., 3)."""
    dx = L / N
    kk = 2 * np.pi / L * np.asarray(n_vec, float)
    return fd4_symbol_1d(kk[..., 0], dx) + fd4_symbol_1d(kk[..., 1], dx) + fd4_symbol_1d(kk[..., 2], dx)


def band_pairs(L, cutoff=8.0, lower=1.8):
    """One representative n of each +-n pair with lower <= 2 pi |n| / L <= cutoff.
    Representative: first nonzero component positive (the documented half-space)."""
    nmax = int(np.ceil(cutoff * L / (2 * np.pi)))
    a = np.arange(-nmax, nmax + 1)
    I, J, K = np.meshgrid(a, a, a, indexing='ij')
    n = np.stack([I.ravel(), J.ravel(), K.ravel()], axis=1)
    first = np.where(n[:, 0] != 0, n[:, 0], np.where(n[:, 1] != 0, n[:, 1], n[:, 2]))
    nn = (n * n).sum(1)
    k = 2 * np.pi / L * np.sqrt(nn)
    keep = (first > 0) & (k >= lower) & (k <= cutoff)
    return n[keep]


def realization_coefficients(match, L, seed=1, cutoff=8.0, fields='all'):
    """Complex Fourier coefficients (per representative n) of the documented finite-band realization.
    Returns n (P,3) and coeff (P,6) = (chi1, chi2, phi, pi1, pi2, pi_phi) / sqrt(V).
    The random key (seed, int(100 L), i+128, j+128, k+128) and the order (14 real normals then 14 imaginary
    normals, /sqrt2; 12 scalar columns then 2 phi columns) follow the documented recipe so that the same
    realization is reproduced; this key is the one non-independent ingredient."""
    keys = match['keys']; lookup = {int(v): i for i, v in enumerate(keys)}
    n = band_pairs(L, cutoff)
    out = np.zeros((len(n), 6), complex)
    for row, (i, j, k) in enumerate(n):
        ix = lookup[int(i * i + j * j + k * k)]
        rng = np.random.default_rng(np.random.SeedSequence([seed, int(L * 100), int(i) + 128, int(j) + 128, int(k) + 128]))
        re = rng.normal(size=14); im = rng.normal(size=14); z = (re + 1j * im) / np.sqrt(2)
        sc = match['evolved_scalar'][ix] @ z[:12]
        ph = match['evolved_phi'][ix] @ z[12:]
        out[row] = np.r_[sc[0], sc[1], ph[0], sc[6], sc[7], ph[1]] / np.sqrt(L ** 3)
    return n, out


def fields_from_coefficients(n, coeff, mean_q, mean_p, N):
    """Real fields q (3,N,N,N), p (3,N,N,N) on the grid x_j = j L / N from the half-space coefficients."""
    C = np.zeros((6, N, N, N), complex)
    for c in range(6):
        C[c][n[:, 0] % N, n[:, 1] % N, n[:, 2] % N] = coeff[:, c]
        C[c][(-n[:, 0]) % N, (-n[:, 1]) % N, (-n[:, 2]) % N] = np.conj(coeff[:, c])
    C[:3, 0, 0, 0] = mean_q; C[3:, 0, 0, 0] = mean_p
    f = np.fft.ifftn(C, axes=(1, 2, 3)) * N ** 3  # f(x) = sum_n c_n exp(i k.x)
    assert np.max(abs(f.imag)) < 1e-8 * max(1.0, np.max(abs(f.real)))
    return f.real[:3].copy(), f.real[3:].copy()


def fd4_gradient_energy_density(f, dx):
    """Mean of 1/2 sum_axes [ (4/3)(f(x+1)-f(x))^2 - (1/12)(f(x+2)-f(x))^2 ] / dx^2 (real space, np.roll)."""
    tot = 0.0
    for ax in range(f.ndim):
        d1 = np.roll(f, -1, axis=ax) - f
        d2 = np.roll(f, -2, axis=ax) - f
        tot += np.mean(4 / 3 * d1 * d1 - d2 * d2 / 12)
    return 0.5 * tot / dx ** 2


def energy_density(q, p, L):
    N = q.shape[-1]; dx = L / N
    r = q[0] ** 2 + q[1] ** 2
    E = 0.5 * np.mean(np.sum(p * p, axis=0))
    E += sum(fd4_gradient_energy_density(q[c], dx) for c in range(3))
    E += np.mean(LX / 4 * r * r + G2 / 2 * q[2] ** 2 * r + LP / 4 * q[2] ** 4)
    return float(E)


def orbit_rhs(x, v):
    r = x @ x
    return v, -LX * r * x


def rk4_orbit_and_transfer(x0, v0, k2, s_out, h=1e-3, mass_extra=None):
    """Fixed-step RK4 for the homogeneous reduced orbit and the 2x2 delta-phi transfer matrices
    T(s) for each k2 (phi'' = -(k2 + g2 |chi|^2 [+ mass_extra(s)]) phi).
    s_out must be increasing from 0.  Returns orbit (len(s_out),4) and T (len(s_out),nk,2,2)."""
    k2 = np.asarray(k2, float); nk = len(k2)
    y = np.zeros(4 + 4 * nk)
    y[:2] = x0; y[2:4] = v0
    Tm = np.zeros((nk, 2, 2)); Tm[:, 0, 0] = 1; Tm[:, 1, 1] = 1
    y[4:] = Tm.ravel()

    def f(s, y):
        x = y[:2]; v = y[2:4]; r = x @ x
        T = y[4:].reshape(nk, 2, 2)
        m2 = k2 + G2 * r + (0.0 if mass_extra is None else mass_extra(s))
        dT = np.empty_like(T)
        dT[:, 0, :] = T[:, 1, :]
        dT[:, 1, :] = -m2[:, None] * T[:, 0, :]
        out = np.empty_like(y)
        out[:2] = v; out[2:4] = -LX * r * x; out[4:] = dT.ravel()
        return out

    orbit = np.zeros((len(s_out), 4)); Ts = np.zeros((len(s_out), nk, 2, 2))
    s = 0.0; j = 0
    while j < len(s_out) and s_out[j] <= 1e-15:
        orbit[j] = y[:4]; Ts[j] = y[4:].reshape(nk, 2, 2); j += 1
    while j < len(s_out):
        target = s_out[j]
        nstep = max(1, int(np.ceil((target - s) / h - 1e-9)))
        hh = (target - s) / nstep
        for _ in range(nstep):
            k1 = f(s, y); k2_ = f(s + hh / 2, y + hh / 2 * k1); k3 = f(s + hh / 2, y + hh / 2 * k2_); k4 = f(s + hh, y + hh * k3)
            y = y + hh / 6 * (k1 + 2 * k2_ + 2 * k3 + k4); s += hh
        s = target
        orbit[j] = y[:4]; Ts[j] = y[4:].reshape(nk, 2, 2); j += 1
    return orbit, Ts


def orbit_turning_points(E, Q):
    """r_min, r_max of the planar quartic orbit with energy E and angular momentum (charge) Q:
    2E = r'^2 + Q^2/r^2 + lamX r^4 / 2  ->  turning points solve lamX/2 u^3 - 2E u + Q^2 = 0, u = r^2."""
    roots = np.roots([LX / 2, 0.0, -2 * E, Q * Q])
    u = np.sort(roots[np.isreal(roots)].real); u = u[u > 0]
    return np.sqrt(u[0]), np.sqrt(u[-1])


def circular_energy(Q):
    """Minimum E at fixed Q (circular orbit): r^6 = Q^2/lamX, E = 3/4 lamX^(1/3) |Q|^(4/3)."""
    return 0.75 * LX ** (1 / 3) * abs(Q) ** (4 / 3)


def radial_period(E, Q, extra_m2=0.0, n=4000):
    """Radial period by Gauss-Chebyshev-type quadrature in r between turning points, for the potential
    lamX/4 r^4 + extra_m2/2 r^2 (Hartree mass for the mean)."""
    roots = np.roots([LX / 2, extra_m2, -2 * E, Q * Q])
    u = np.sort(roots[np.isreal(roots)].real); u = u[u > 0]
    rmin, rmax = np.sqrt(u[0]), np.sqrt(u[-1])
    c, hw = (rmin + rmax) / 2, (rmax - rmin) / 2
    th = (np.arange(n) + 0.5) * np.pi / n - np.pi / 2
    r = c + hw * np.sin(th)
    pr2 = 2 * E - Q * Q / r ** 2 - LX * r ** 4 / 2 - extra_m2 * r * r
    return 2 * np.sum(hw * np.cos(th) / np.sqrt(np.maximum(pr2, 1e-300))) * np.pi / n, rmin, rmax


def floquet_on_rseries(r_of_s, s_grid, k2, extra=0.0):
    """Monodromy of phi'' = -(k2 + g2 r(s) + extra) phi over the given s-grid (one period), RK4 on the grid
    with r linearly interpolated at half steps.  Returns mu = ln(max|eig|)/T for each k2."""
    k2 = np.asarray(k2, float); nk = len(k2)
    T = np.zeros((nk, 2, 2)); T[:, 0, 0] = 1; T[:, 1, 1] = 1

    def f(r, T):
        d = np.empty_like(T); d[:, 0, :] = T[:, 1, :]; d[:, 1, :] = -(k2 + G2 * r + extra)[:, None] * T[:, 0, :]; return d
    for i in range(len(s_grid) - 1):
        h = s_grid[i + 1] - s_grid[i]; r0, r1 = r_of_s[i], r_of_s[i + 1]; rm = 0.5 * (r0 + r1)
        a = f(r0, T); b = f(rm, T + h / 2 * a); c = f(rm, T + h / 2 * b); d = f(r1, T + h * c)
        T = T + h / 6 * (a + 2 * b + 2 * c + d)
    tr = T[:, 0, 0] + T[:, 1, 1]
    lam = np.abs(tr) / 2 + np.sqrt(np.maximum(tr * tr / 4 - 1, 0))
    period = s_grid[-1] - s_grid[0]
    return np.log(np.maximum(lam, 1.0)) / period, T
