"""Independent 3D evolution of the reduced conformal X/phi Hamiltonian (Claude review track en_indep).

Own code: nothing is imported from Codex's packet (no early3d.cpp, no en_common lap/gradient/energy).
Only the published *input* data are read: the tau=20 initial field file (built by Codex's initial.py,
hash-checked) and the clock constants a_m, C from MATCH_L12.npz.

Model (MATHEMATICS.md of grut_early_nonlinear_20261027):
  H = sum_x [ 1/2 pi^2 + 1/2 chi(-Lap)chi + lX/4 r^2 + g2/2 phi^2 r + lP/4 phi^4 ],  r = chi1^2 + chi2^2,
  lX = 1e-3, g2 = 1e-2, lP = 3e-2; s = conformal time (H_i units); Hamiltonian autonomous in s.
Lattice operator: symmetric FD4 Laplacian, implemented spectrally with its exact symbol
  sum_axes 4 sin^2(k dx/2) (1 + sin^2(k dx/2)/3) / dx^2   (identical operator to the 5-point FD4 stencil).
Time integrator (different from Codex's Yoshida/Forest-Ruth triple-jump):
  'bm6'  : Blanes-Moan (2002) SRKN_6^b, 4th-order symmetric symplectic Runge-Kutta-Nystrom splitting (7 kicks, FSAL -> 6/step)
  'rk4'  : classical (non-symplectic) Runge-Kutta 4 on (q, p)       -- cross-check only
  'vv'   : velocity Verlet (2nd order)                               -- cross-check only
Clock: a(s) = a_m + C s (radiation, a''=0), d tau = a ds  =>  tau(s) = 20 + a_m s + C s^2/2,
  s(tau) = 2 (tau-20) / (a_m + sqrt(a_m^2 + 2 C (tau-20)))  (cancellation-free root).
Output grid: uniform in tau (spacing dtau_out); each output interval is split into n equal substeps with
  n = ceil(ds_interval / h0) (so the actual step never exceeds h0 and lands exactly on the output time).
"""
import sys, os, json, time, struct, hashlib, argparse
import numpy as np
from scipy import fft as sfft

LX, G2, LP = 1e-3, 1e-2, 3e-2
WORKERS = int(os.environ.get("INDEP_WORKERS", "2"))


def read_initial(path):
    with open(path, 'rb') as f:
        magic = f.read(8)
        assert magic == b'GRUT3D1\n', magic
        N, L, nu, g2, lx, lp = struct.unpack('<id4d', f.read(4 + 32 + 8))
        raw = np.fromfile(f, dtype='<f8')
    assert raw.size == 6 * N ** 3
    assert nu == 0.0 and abs(g2 - G2) < 1e-15 and abs(lx - LX) < 1e-15 and abs(lp - LP) < 1e-15
    raw = raw.reshape(6, N, N, N)
    return N, L, raw[:3].copy(), raw[3:].copy()


class Lattice:
    def __init__(self, N, L):
        self.N, self.L, self.dx = N, L, L / N
        th_full = 2 * np.pi * np.fft.fftfreq(N)          # k*dx on the full axes
        th_half = 2 * np.pi * np.fft.rfftfreq(N)         # k*dx on the last (real) axis
        def sym(th):
            s2 = np.sin(th / 2) ** 2
            return 4 * s2 * (1 + s2 / 3) / self.dx ** 2
        self.K2 = sym(th_full)[:, None, None] + sym(th_full)[None, :, None] + sym(th_half)[None, None, :]
        w = np.full(N // 2 + 1, 2.0); w[0] = 1.0
        if N % 2 == 0: w[-1] = 1.0
        self.wt = np.broadcast_to(w[None, None, :], self.K2.shape)   # rfft Parseval weights

    def minus_lap(self, q):
        qk = sfft.rfftn(q, axes=(-3, -2, -1), workers=WORKERS)
        qk *= self.K2
        return sfft.irfftn(qk, s=q.shape[-3:], axes=(-3, -2, -1), workers=WORKERS)

    def force(self, q):
        F = -self.minus_lap(q)
        r = q[0] * q[0] + q[1] * q[1]
        ph2 = q[2] * q[2]
        cX = LX * r + G2 * ph2
        F[0] -= cX * q[0]
        F[1] -= cX * q[1]
        F[2] -= (LP * ph2 + G2 * r) * q[2]
        return F

    def grad_energy_density(self, q):
        """(1/V) sum 1/2 chi(-Lap)chi by Parseval, per field."""
        qk = sfft.rfftn(q, axes=(-3, -2, -1), workers=WORKERS)
        n3 = float(self.N) ** 3
        # mean_x f(-Lap)f = (1/N^6) sum_k K2 |f_k|^2 (unnormalized FFT)
        return 0.5 * np.einsum('xyz,cxyz->c', self.wt * self.K2, (qk.real ** 2 + qk.imag ** 2)) / n3 ** 2


# Blanes & Moan, J. Comput. Appl. Math. 142 (2002) 313, SRKN_6^b (B = kick first)
BM_b = [0.0829844064174052, 0.396309801498368, -0.0390563049223486]
BM_b.append(1 - 2 * sum(BM_b))
BM_a = [0.245298957184271, 0.604872665711080]
BM_a.append(0.5 - sum(BM_a))
BM_KICK = [BM_b[0], BM_b[1], BM_b[2], BM_b[3], BM_b[2], BM_b[1], BM_b[0]]
BM_DRIFT = [BM_a[0], BM_a[1], BM_a[2], BM_a[2], BM_a[1], BM_a[0]]


def step_bm6(lat, q, p, h, nsub, F=None):
    """nsub BM6 steps of size h; FSAL: the closing kick of one step merges with the next opening kick."""
    if F is None: F = lat.force(q)
    for n in range(nsub):
        for i in range(6):
            p += (BM_KICK[i] * h) * F
            q += (BM_DRIFT[i] * h) * p
            F = lat.force(q)
        p += (BM_KICK[6] * h) * F
    return q, p, F


def step_vv(lat, q, p, h, nsub, F=None):
    if F is None: F = lat.force(q)
    for n in range(nsub):
        p += (0.5 * h) * F
        q += h * p
        F = lat.force(q)
        p += (0.5 * h) * F
    return q, p, F


def step_rk4(lat, q, p, h, nsub, F=None):
    for n in range(nsub):
        k1q, k1p = p, lat.force(q)
        k2q, k2p = p + 0.5 * h * k1p, lat.force(q + 0.5 * h * k1q)
        k3q, k3p = p + 0.5 * h * k2p, lat.force(q + 0.5 * h * k2q)
        k4q, k4p = p + h * k3p, lat.force(q + h * k3q)
        q = q + (h / 6) * (k1q + 2 * k2q + 2 * k3q + k4q)
        p = p + (h / 6) * (k1p + 2 * k2p + 2 * k3p + k4p)
    return q, p, None


STEPPERS = dict(bm6=step_bm6, vv=step_vv, rk4=step_rk4)
COLS = ['s', 'tau', 'a', 'E', 'KX', 'Kphi', 'GX', 'Gphi', 'UX', 'Uphi', 'Umix', 'nQ', 'jvar',
        'm1', 'm2', 'mphi', 'n1', 'n2', 'nphi', 'Xsq', 'phisq', 'maxXsq', 'work', 'wall']


def diagnostics(lat, q, p, s, tau, a, wall):
    r = q[0] * q[0] + q[1] * q[1]
    ph2 = q[2] * q[2]
    j = q[0] * p[1] - q[1] * p[0]
    KX = 0.5 * float(np.mean(p[0] ** 2 + p[1] ** 2)); Kp = 0.5 * float(np.mean(p[2] ** 2))
    G = lat.grad_energy_density(q)
    UX = LX / 4 * float(np.mean(r * r)); Up = LP / 4 * float(np.mean(ph2 * ph2)); Um = G2 / 2 * float(np.mean(ph2 * r))
    E = KX + Kp + G[0] + G[1] + G[2] + UX + Up + Um
    nQ = float(np.mean(j)); jvar = float(np.mean((j - nQ) ** 2))
    m = q.mean(axis=(1, 2, 3)); n = p.mean(axis=(1, 2, 3))
    work = G2 * float(np.mean(q[2] * p[2] * r))
    return [s, tau, a, E, KX, Kp, G[0] + G[1], G[2], UX, Up, Um, nQ, jvar, m[0], m[1], m[2], n[0], n[1], n[2],
            float(np.mean(r)), float(np.mean(ph2)), float(r.max()), work, wall]


def derived(rec):
    """Codex-defined diagnostics, implemented from MATHEMATICS.md definitions (own code)."""
    d = {c: rec[:, i] for i, c in enumerate(COLS)}
    E = d['E']
    mr = d['m1'] ** 2 + d['m2'] ** 2
    E0 = 0.5 * (d['n1'] ** 2 + d['n2'] ** 2 + d['nphi'] ** 2) + LX / 4 * mr ** 2 + G2 / 2 * d['mphi'] ** 2 * mr + LP / 4 * d['mphi'] ** 4
    Ephi = d['Kphi'] + d['Gphi'] + d['Uphi']
    Q0 = d['m1'] * d['n2'] - d['m2'] * d['n1']
    return dict(tau=d['tau'], s=d['s'], E=E, nQ=d['nQ'],
                energy_drift=E / E[0] - 1, charge_drift=d['nQ'] / d['nQ'][0] - 1,
                outside_mean=1 - E0 / E, mean_charge_fraction=Q0 / d['nQ'],
                Ephi_fraction=Ephi / E, phi_plus_mix_fraction=(Ephi + d['Umix']) / E,
                charge_contrast=np.sqrt(np.maximum(d['jvar'], 0)) / np.abs(d['nQ']),
                Xvar_fraction=(d['Xsq'] - mr) / d['Xsq'], maxXsq=d['maxXsq'], work=d['work'])


def tau_of_s(s, am, C, tau0=20.0):
    return tau0 + am * s + 0.5 * C * s * s


def s_of_tau(tau, am, C, tau0=20.0):
    x = tau - tau0
    return 2 * x / (am + np.sqrt(am * am + 2 * C * x))


def run(args):
    out = args.out
    os.makedirs(out, exist_ok=True)
    am, C = args.am, args.C
    step = STEPPERS[args.method]
    ck = os.path.join(out, 'checkpoint.npz')
    recf = os.path.join(out, 'records.npy')
    if os.path.exists(ck) and not args.fresh:
        z = np.load(ck)
        q, p, s, tau, N, L = z['q'], z['p'], float(z['s']), float(z['tau']), int(z['N']), float(z['L'])
        recs = list(np.load(recf))
        lat = Lattice(N, L)
        print('resume at tau', tau, flush=True)
    else:
        N, L, q, p = read_initial(args.init)
        lat = Lattice(N, L)
        s, tau = 0.0, 20.0
        recs = [diagnostics(lat, q, p, s, tau, am, 0.0)]
    t0 = time.monotonic(); F = None; nsteps = 0; maxh = 0.0
    last_ck = time.monotonic()
    while tau < args.until - 1e-9 and time.monotonic() - t0 < args.budget:
        dt_out = args.dtau if tau < args.fine_until - 1e-9 else args.dtau_late
        tnext = min(round((tau + dt_out) / 1e-9) * 1e-9, args.until)
        if tau < args.fine_until - 1e-9: tnext = min(tnext, args.fine_until)
        snext = s_of_tau(tnext, am, C)
        nsub = int(np.ceil((snext - s) / args.h0 - 1e-12))
        h = (snext - s) / nsub
        maxh = max(maxh, h)
        q, p, F = step(lat, q, p, h, nsub, F)
        nsteps += nsub
        s, tau = snext, tnext
        recs.append(diagnostics(lat, q, p, s, tau, am + C * s, time.monotonic() - t0))
        if not np.isfinite(recs[-1][3]): raise RuntimeError('nonfinite')
        if time.monotonic() - last_ck > 120:
            np.savez(ck, q=q, p=p, s=s, tau=tau, N=N, L=L); np.save(recf, np.array(recs)); last_ck = time.monotonic()
        if args.snap_taus and any(abs(tau - t) < 1e-7 for t in args.snap_taus):
            np.savez(os.path.join(out, f'snap_tau{tau:.2f}.npz'), q=q, p=p, s=s, tau=tau)
    np.savez(ck, q=q, p=p, s=s, tau=tau, N=N, L=L); np.save(recf, np.array(recs))
    print(f'reached tau {tau:.4f} s {s:.6f} steps {nsteps} maxh {maxh:.5f} seconds {time.monotonic()-t0:.1f}', flush=True)
    return tau


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--init', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--am', type=float, required=True)
    ap.add_argument('--C', type=float, required=True)
    ap.add_argument('--method', default='bm6')
    ap.add_argument('--h0', type=float, default=0.004)
    ap.add_argument('--dtau', type=float, default=0.05)
    ap.add_argument('--dtau_late', type=float, default=0.05)
    ap.add_argument('--fine_until', type=float, default=1e9)
    ap.add_argument('--until', type=float, default=800.0)
    ap.add_argument('--budget', type=float, default=540.0)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--snap_taus', type=float, nargs='*', default=[])
    run(ap.parse_args())
