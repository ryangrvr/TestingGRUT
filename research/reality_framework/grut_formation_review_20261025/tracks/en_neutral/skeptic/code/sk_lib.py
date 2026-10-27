"""en_neutral_skeptic: own helper library (Claude adversarial check of track en_neutral, Codex 7d9c93b, issue #13).

INDEPENDENT of both Codex's and en_neutral's analysis code:
  * FD4 symbol written from the 5-point stencil (-1, 16, -30, 16, -1)/12 (en_neutral used the sin^2 form);
  * seeds taken by full complex FFT of the tau-20 initial field file (all lattice modes, both +-n);
  * homogeneous orbits from scipy DOP853 (rtol 1e-12) with dense output; measured pumps from cubic splines of the
    engine trace (en_neutral used linear interpolation);
  * delta-phi transfer matrices by RK4 on a grid refined x2 relative to the trace grid, with the pump pre-evaluated
    at every RK4 stage point.
Data read (never modified): the en_neutral reproduction run of Codex's unmodified engine (state.tsv, initial.bin),
which en_neutral and we both checked against the published N48 curve.

Reduced conformal Hamiltonian (Codex MATHEMATICS.md):
  H = 1/2 pi^2 + 1/2 |grad chi|^2 + lamX/4 |chi_X|^4 + g2/2 chi_phi^2 |chi_X|^2 + lamP/4 chi_phi^4
"""
import struct
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

LX, G2, LP = 1e-3, 1e-2, 3e-2


def read_initial(path):
    with open(path, 'rb') as f:
        magic = f.read(8)
        N, L, _s, g2, lx, lp = struct.unpack('<id4d', f.read(44))
        a = np.fromfile(f, dtype='<f8').reshape(6, N, N, N)
    assert magic == b'GRUT3D1\n' and abs(g2 - G2) < 1e-15 and abs(lx - LX) < 1e-15 and abs(lp - LP) < 1e-15
    return N, L, a[:3].copy(), a[3:].copy()


def fd4_k2_axis(m, N, L):
    """Symbol of -(f(x+2dx)(-1) + 16 f(x+dx) - 30 f + 16 f(x-dx) - f(x-2dx))/(12 dx^2) for integer index m."""
    dx = L / N
    th = 2 * np.pi * np.asarray(m, float) / N
    return (30 - 32 * np.cos(th) + 2 * np.cos(2 * th)) / (12 * dx * dx)


def lattice_k2_grid(N, L):
    m = np.fft.fftfreq(N) * N
    a = fd4_k2_axis(m, N, L)
    return a[:, None, None] + a[None, :, None] + a[None, None, :]


def phi_seed_by_k2(path, cont_band=None):
    """Full complex FFT (norm 1/N^3) of phi and pi_phi; group nonzero modes by unique lattice k^2.
    Returns dict with unique k2 values, summed |phi|^2, |pi|^2, Re(phi* pi) per k2, plus per-mode arrays."""
    N, L, q, p = read_initial(path)
    K2 = lattice_k2_grid(N, L)
    ph = np.fft.fftn(q[2]) / N ** 3
    pp = np.fft.fftn(p[2]) / N ** 3
    use = (abs(ph) + abs(pp)) > 1e-12
    k2 = K2[use]
    vals, ids = np.unique(np.round(k2, 10), return_inverse=True)
    nk = len(vals)
    qq = np.bincount(ids, weights=abs(ph[use]) ** 2, minlength=nk)
    PP = np.bincount(ids, weights=abs(pp[use]) ** 2, minlength=nk)
    qp = np.bincount(ids, weights=np.real(np.conj(ph[use]) * pp[use]), minlength=nk)
    cnt = np.bincount(ids, minlength=nk)
    means = dict(q=q.mean(axis=(1, 2, 3)), p=p.mean(axis=(1, 2, 3)))
    return dict(N=N, L=L, k2=vals, qq=qq, pp=PP, qp=qp, count=cnt, nmodes=int(use.sum()), means=means, q=q, p=p)


def energy_density_realspace(q, p, L):
    """Own real-space energy density with the 5-point FD4 gradient energy written as -1/2 f . Lap f."""
    N = q.shape[-1]; dx = L / N
    def lap(f):
        out = np.zeros_like(f)
        for ax in range(3):
            out = out - 30 * f + 16 * (np.roll(f, 1, ax) + np.roll(f, -1, ax)) - (np.roll(f, 2, ax) + np.roll(f, -2, ax))
        return out / (12 * dx * dx)
    r = q[0] ** 2 + q[1] ** 2
    E = 0.5 * np.mean(np.sum(p * p, axis=0))
    E += sum(-0.5 * np.mean(q[c] * lap(q[c])) for c in range(3))
    E += np.mean(LX / 4 * r * r + G2 / 2 * q[2] ** 2 * r + LP / 4 * q[2] ** 4)
    return float(E)


def free_orbit_dense(x0, v0, s1, rtol=1e-12, atol=1e-10):
    def rhs(s, y):
        r = y[0] ** 2 + y[1] ** 2
        return [y[2], y[3], -LX * r * y[0], -LX * r * y[1]]
    sol = solve_ivp(rhs, (0.0, s1), np.r_[x0, v0], method='DOP853', rtol=rtol, atol=atol, dense_output=True)
    assert sol.success
    return sol


def turning_points(E, Q):
    """r_min, r_max of the planar quartic orbit: 2E = r'^2 + Q^2/r^2 + lamX r^4/2 -> lamX/2 u^3 - 2E u + Q^2 = 0."""
    roots = np.roots([LX / 2, 0.0, -2 * E, Q * Q])
    u = np.sort(roots[np.isreal(roots)].real); u = u[u > 0]
    return np.sqrt(u[0]), np.sqrt(u[-1])


def one_free_period(E, Q):
    """Free orbit from pericentre (r_min, 0) with tangential velocity Q/r_min over one radial period
    (period found as the next pericentre by an event on the radial velocity)."""
    rmin, rmax = turning_points(E, Q)
    def rhs(s, y):
        r = y[0] ** 2 + y[1] ** 2
        return [y[2], y[3], -LX * r * y[0], -LX * r * y[1]]
    def peri(s, y):  # radial velocity x.v crosses zero from below at pericentre
        return y[0] * y[2] + y[1] * y[3]
    peri.direction = 1.0
    sol = solve_ivp(rhs, (0.0, 3.0), [rmin, 0.0, 0.0, Q / rmin], method='DOP853', rtol=1e-12, atol=1e-10,
                    dense_output=True, events=peri)
    ev = [t for t in sol.t_events[0] if t > 0.05]
    T = float(ev[0])
    return sol, T, rmin, rmax


def propagate(k2, m2_left, m2_mid, m2_right, s_grid, record_idx):
    """RK4 transfer matrices for phi'' = -(k2 + m2(s)) phi over s_grid, with m2 pre-evaluated at the
    left / mid / right stage points of each step. Returns T at record_idx (list of grid indices)."""
    nk = len(k2)
    T = np.zeros((nk, 2, 2)); T[:, 0, 0] = 1; T[:, 1, 1] = 1
    rec = {}
    rset = set(int(i) for i in record_idx)
    if 0 in rset: rec[0] = T.copy()
    def f(m, T):
        d = np.empty_like(T); d[:, 0, :] = T[:, 1, :]; d[:, 1, :] = -(k2 + m)[:, None] * T[:, 0, :]; return d
    for i in range(len(s_grid) - 1):
        h = s_grid[i + 1] - s_grid[i]
        a = f(m2_left[i], T); b = f(m2_mid[i], T + 0.5 * h * a); c = f(m2_mid[i], T + 0.5 * h * b); d = f(m2_right[i], T + h * c)
        T = T + h / 6 * (a + 2 * b + 2 * c + d)
        if i + 1 in rset: rec[i + 1] = T.copy()
    return rec


def energy_from_T(T, seed, m2_energy):
    """phi + interaction energy density sum_n 1/2(|pi_n|^2 + (k2 + m2)|phi_n|^2), per unique k2 and total."""
    t00, t01, t10, t11 = T[:, 0, 0], T[:, 0, 1], T[:, 1, 0], T[:, 1, 1]
    Qn = t00 ** 2 * seed['qq'] + 2 * t00 * t01 * seed['qp'] + t01 ** 2 * seed['pp']
    Pn = t10 ** 2 * seed['qq'] + 2 * t10 * t11 * seed['qp'] + t11 ** 2 * seed['pp']
    e = 0.5 * (Pn + (seed['k2'] + m2_energy) * Qn)
    return e


def read_trace(path):
    return np.genfromtxt(path, names=True, delimiter='\t')
