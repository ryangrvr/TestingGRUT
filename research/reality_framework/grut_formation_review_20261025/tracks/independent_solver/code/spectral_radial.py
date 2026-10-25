#!/usr/bin/env python3
"""Independent radial solver for the prepared GRUT two-field collapse (Claude review track
`independent_solver`, 25 Oct 2026).

Written from the equations in Codex's MATHEMATICS.md (c4b99f3 section 4; bb211d1 section 2), not
from Codex's code.  Dimensionless equations (tau = g_X v t, rho = g_X v r, X = v psi/sqrt2, phi = v S):

    psi_tt = lap psi - (S^2 + a|psi|^2) psi,      S_tt = lap S - [b(S^2-1) + |psi|^2] S,
    a = lambda_X/g_X^2 = 0.1,  b = lambda_phi/g_X^2 = 3.

Radial variables U = r psi (complex, mass 1) and Z = r (S - 1) (real, mass sqrt 6).

Method (deliberately different from Codex's finite-difference + RK4 + CFS layer):
  * Space: Fourier-sine (DST-I) pseudo-spectral representation on [0, L] with U = Z = 0 at r = 0
    (regularity, odd continuation) and at r = L (hard wall).
  * Time: symplectic splitting.  The free Klein-Gordon part (k^2 + m^2) is propagated EXACTLY per sine
    mode (rotation by omega_k tau); the nonlinear remainder is a position-only kick.  Strang
    (2nd order) or Yoshida triple-jump (4th order) composition.  The kick conserves Noether charge
    exactly (real multiplier), the free rotation conserves it exactly, so charge is conserved to
    round-off; energy error is bounded (symplectic).
  * Boundary: (a) "big box": L chosen so that no signal reflected from the wall can return to the
    measured region within the run (cleanest control; no absorber at all);
    (b) optional smooth velocity sponge (gamma(r) p damping) for long runs, applied exactly
    inside the kick, with an energy/charge receiver.  It is tested separately with free packets.

Units of Q, E: the barred dimensionless quantities (Qbar = Q g_X^2), as in Codex's tables.
"""
import argparse, json, time, hashlib
from pathlib import Path
import numpy as np
from scipy.fft import dst

A_X = 0.1   # lambda_X / g_X^2
B_PHI = 3.0  # lambda_phi / g_X^2
M2 = np.array([1.0, 1.0, 6.0])  # vacuum mass^2 of rows x, y (charged) and Z (neutral)
W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0)); W0 = 1.0 - 2.0 * W1


def T(a):
    """Orthonormal DST-I along the last axis; it is its own inverse."""
    return dst(a, type=1, norm='ortho', axis=-1)


class Radial:
    def __init__(self, h, L, dt, scheme='strang', nonlinear=True, sponge=None, masses=None):
        self.h, self.L, self.dt, self.scheme, self.nonlinear = h, L, dt, scheme, nonlinear
        n = int(round(L / h)) - 1
        self.n = n
        self.r = h * np.arange(1, n + 1)
        self.k = np.pi * np.arange(1, n + 1) / L
        m2 = M2 if masses is None else np.asarray(masses, float)
        self.m2 = m2
        self.omega = np.sqrt(self.k[None, :] ** 2 + m2[:, None])  # (3, n)
        self.wmax = float(self.omega.max())
        self.Qh = np.zeros((3, n)); self.Ph = np.zeros((3, n))  # spectral positions / momenta
        self.tau = 0.0
        self.Erec = 0.0; self.Qrec = 0.0
        # sponge: dict(start=R_s, width=D, gmax=..., power=p) -> gamma(r) = gmax*((r-R_s)/D)^p (clipped)
        self.sponge = sponge
        if sponge:
            xi = np.clip((self.r - sponge['start']) / sponge['width'], 0.0, 1.0)
            self.gamma = sponge['gmax'] * xi ** sponge['power']
            self.gmask = self.gamma > 0
        else:
            self.gamma = None
        self._rot = {}
        if scheme == 'strang':
            self.drifts = [1.0]; self.kicks = [1.0]  # merged interior kick weights
            self.first_kick = 0.5
        elif scheme == 'yoshida4':
            self.drifts = [W1, W0, W1]
            self.kicks = [(W1 + W0) / 2, (W0 + W1) / 2, W1]
            self.first_kick = W1 / 2
        else:
            raise ValueError(scheme)
        theta = self.wmax * dt * max(abs(d) for d in self.drifts)
        assert theta < 0.9 * np.pi, f'step-size resonance risk: omega_max*dt={theta}'

    # ------------------------------------------------------------------ dynamics
    def rot(self, s):
        key = round(s, 15)
        if key not in self._rot:
            c = np.cos(self.omega * s); sn = np.sin(self.omega * s)
            self._rot[key] = (c, sn / self.omega, -sn * self.omega)
        return self._rot[key]

    def drift(self, s):
        c, so, wo = self.rot(s)
        q, p = self.Qh, self.Ph
        self.Qh, self.Ph = c * q + so * p, wo * q + c * p

    def force(self, q):
        """Nonlinear (non-free) force per unit 4 pi h; q rows = (x, y, Z) at nodes."""
        r = self.r
        x, y, z = q
        u2 = x * x + y * y
        S = 1.0 + z / r
        g = (S * S - 1.0) + A_X * u2 / (r * r)
        fz = -(9.0 * z * z / r + 3.0 * z ** 3 / (r * r) + u2 * S / r)
        return np.array([-g * x, -g * y, fz])

    def kick(self, s, want_sync=None):
        """Momentum update over time s.  With sponge: p -> (p + s/2 F) e^{-gamma s} + s/2 F.
        If want_sync is a float w, also return synchronized physical (q, p) at the point where the
        first half of this kick (weight w of s... see run loop) has been applied."""
        q = T(self.Qh)
        F = self.force(q) if self.nonlinear else np.zeros_like(q)
        if self.gamma is None:
            self.Ph = self.Ph + s * T(F)
            return q, F
        p = T(self.Ph)
        p1 = p + 0.5 * s * F
        damp = np.exp(-self.gamma * s)
        p2 = p1 * damp
        self._receiver(q, p1, p2)
        p3 = p2 + 0.5 * s * F
        self.Ph = T(p3)
        return q, F

    def _receiver(self, q, p_before, p_after):
        m = self.gmask
        c = 4 * np.pi * self.h
        self.Erec += c * 0.5 * np.sum(p_before[:, m] ** 2 - p_after[:, m] ** 2)
        x, y = q[0, m], q[1, m]
        qb = np.sum(y * p_before[0, m] - x * p_before[1, m])
        qa = np.sum(y * p_after[0, m] - x * p_after[1, m])
        self.Qrec += c * (qb - qa)

    def advance_block(self, nsteps):
        """Advance nsteps steps; state synchronized (positions and momenta at the same time) at
        entry and exit.  Interior kicks are merged (first-same-as-last).  With a sponge the kick is
        F(s/2) D(s) F(s/2) (D = exact exponential damping), only with the Strang scheme, so that no
        negative damping sub-step occurs."""
        dt = self.dt
        # With the Yoshida composition one kick weight is negative; the damping factor exp(-gamma s)
        # is then applied for negative s (exact backward damping flow).  The composition remains
        # consistent; the receiver ledger is signed and integrates exactly what is removed.
        self.kick(self.first_kick * dt)
        for i in range(nsteps):
            for j, (d, kw) in enumerate(zip(self.drifts, self.kicks)):
                self.drift(d * dt)
                last = (i == nsteps - 1) and (j == len(self.drifts) - 1)
                self.kick((self.first_kick if last else kw) * dt)
        self.tau += nsteps * dt

    # ------------------------------------------------------------------ diagnostics
    def physical(self):
        return T(self.Qh), T(self.Ph)

    def derivative(self, Qh):
        """d/dr of the sine series at the nodes, via the cosine series (DCT-I on the extended grid)."""
        # U(r) = sqrt(2/L) sum_k Qh_k sin(k r) with ortho DST-I normalisation sqrt(2/(n+1)) on nodes:
        # node values q_j = sqrt(2/(n+1)) sum_m Qh_m sin(pi m j/(n+1)).  Derivative at node j:
        # sqrt(2/(n+1)) sum_m Qh_m k_m cos(pi m j/(n+1)).  Compute with a DCT-I of length n+2
        # (nodes j = 0..n+1) on coefficients [0, k*Qh..., 0].
        from scipy.fft import dct
        n = self.n
        coef = np.zeros(Qh.shape[:-1] + (n + 2,))
        coef[..., 1:n + 1] = Qh * self.k
        # DCT-I (unnormalised): y_j = c_0 + (-1)^j c_{N-1} + 2 sum_{m=1}^{N-2} c_m cos(pi m j/(N-1)), N=n+2
        y = dct(coef, type=1, axis=-1)
        return np.sqrt(2.0 / (n + 1)) * 0.5 * y  # includes j = 0 and j = n+1

    def diagnostics(self, radii=(30.0, 60.0, 120.0, 240.0, 480.0)):
        q, p = self.physical()
        r, h = self.r, self.h
        c = 4 * np.pi * h
        x, y, z = q; px, py, pz = p
        u2 = x * x + y * y
        S = 1.0 + z / r
        dq = self.derivative(self.Qh)[:, 1:-1]
        # physical (covariant) local densities per unit r (times 4 pi): |U' - U/r|^2 etc.
        gU = (dq[0] - x / r) ** 2 + (dq[1] - y / r) ** 2
        gZ = (dq[2] - z / r) ** 2
        kinU = 0.5 * (px * px + py * py); kinZ = 0.5 * pz * pz
        potU = 0.5 * S * S * u2 + 0.25 * A_X * u2 * u2 / (r * r)
        potZ = 0.75 * z * z * (2.0 + z / r) ** 2
        e_ch = kinU + 0.5 * gU + 0.5 * u2        # charged free-like part (mass term m^2|U|^2/2)
        e_int = 0.5 * (S * S - 1.0) * u2 + 0.25 * A_X * u2 * u2 / (r * r)
        e_ne = kinZ + 0.5 * gZ + potZ
        qd = y * px - x * py
        dp = self.derivative(self.Ph)[:, 1:-1]
        # d/dr of the charge density (spectral derivatives) for the Euler-Maclaurin end correction
        dqd = dq[1] * px + y * dp[0] - dq[0] * py - x * dp[1]
        # total energy from spectral Parseval (exact discrete Hamiltonian)
        Ekin = 0.5 * np.sum(self.Ph ** 2)
        Egrad = 0.5 * np.sum((self.k[None, :] * self.Qh) ** 2)
        E = c * (Ekin + Egrad + np.sum(potU + potZ))
        Q = c * np.sum(qd)
        out = dict(tau=self.tau, E=E, Q=Q, Erec=self.Erec, Qrec=self.Qrec)
        for R in radii:
            j = int(round(R / h))  # node index of r = R is j-1 in zero-based arrays
            inside_lt = slice(0, j - 1)          # r < R
            node = j - 1
            trap = c * (np.sum(qd[inside_lt]) + 0.5 * qd[node])
            # Euler-Maclaurin: int_0^R q dr = h*sum' q - h^2/12 (q'(R) - q'(0)), q'(0) = 0 (q even)
            out[f'Q{int(R)}'] = trap - 4 * np.pi * h * h / 12.0 * dqd[node]
            out[f'Q{int(R)}_trap'] = trap
            out[f'Q{int(R)}_le'] = c * np.sum(qd[:node + 1])
            ein = e_ch + e_int + e_ne
            out[f'E{int(R)}'] = c * (np.sum(ein[inside_lt]) + 0.5 * ein[node])
            out[f'Ech_out{int(R)}'] = c * (np.sum(e_ch[node + 1:]) + 0.5 * e_ch[node])
            out[f'Ene_out{int(R)}'] = c * (np.sum(e_ne[node + 1:]) + 0.5 * e_ne[node])
            out[f'Eint_out{int(R)}'] = c * (np.sum(e_int[node + 1:]) + 0.5 * e_int[node])
        # local density check: sum of local energy densities vs Parseval total
        out['E_local_sum'] = c * np.sum(e_ch + e_int + e_ne)
        # nonlinear mass shift outside radius 240
        o = r > 240
        out['max_dm2_out240'] = float(np.max(np.abs(S[o] ** 2 - 1 + A_X * u2[o] / r[o] ** 2))) if o.any() else 0.0
        out['S0'] = float(S[0]); out['psi0'] = float(np.sqrt(u2[0]) / r[0])
        return out

    def mode_numbers(self, rmin, taper=20.0):
        """Free-field particle/antiparticle numbers of the exterior field r > rmin (smooth mask),
        projected on the exact sine eigenmodes of the box (free KG, mass 1 / sqrt 6)."""
        q, p = self.physical()
        r = self.r
        w = 0.5 * (1 + np.tanh((r - rmin) / taper))
        U = (q[0] + 1j * q[1]) * w; V = (p[0] + 1j * p[1]) * w
        Uh, Vh = T(U), T(V)
        om = self.omega[0]
        c = 4 * np.pi * self.h
        ap = (om * Uh + 1j * Vh) / (2 * om)
        am = (om * Uh - 1j * Vh) / (2 * om)
        Np = c * om * np.abs(ap) ** 2; Nm = c * om * np.abs(am) ** 2
        Zh, Wh = T(q[2] * w), T(p[2] * w)
        oz = self.omega[2]
        Ns = c * (np.abs(oz * Zh + 1j * Wh) ** 2) / (2 * oz)
        res = dict(rmin=rmin, taper=taper, N_plus=float(Np.sum()), N_minus=float(Nm.sum()),
                   E_plus=float((om * Np).sum()), E_minus=float((om * Nm).sum()),
                   N_s=float(Ns.sum()), E_s=float((oz * Ns).sum()),
                   Q_masked=float(c * np.sum((U.conj() * V).imag) * -1.0))
        # momentum spectra (binned) for record
        edges = np.array([0, .005, .01, .02, .05, .1, .2, .3, .4, .5, .7, 1, 1.5, 2, 3, 100])
        idx = np.digitize(self.k, edges) - 1
        res['k_edges'] = edges.tolist()
        res['N_plus_binned'] = np.bincount(idx, Np, len(edges) - 1).tolist()
        res['N_minus_binned'] = np.bincount(idx, Nm, len(edges) - 1).tolist()
        res['E_s_binned'] = np.bincount(idx, oz * Ns, len(edges) - 1).tolist()
        kmean_p = float((self.k * Np).sum() / max(Np.sum(), 1e-300))
        res['mean_k_plus'] = kmean_p
        return res

    # ------------------------------------------------------------------ initial data
    def prepared_gaussian(self, A=0.18, R=60.0, tol=1e-14):
        r = self.r
        f = A * np.exp(-r * r / (2 * R * R))
        U = r * f
        # static neutral equation, spectral fixed point:  (d^2 - 6) Z = N(Z) with
        # N(Z) = 9 Z^2/r + 3 Z^3/r^2 + f^2 r S  (exact rearrangement of Z'' = r[3(S^2-1)+f^2] S)
        lam = -(self.k ** 2 + 6.0)
        z = np.zeros_like(r)
        for it in range(200):
            S = 1 + z / r
            N = 9 * z * z / r + 3 * z ** 3 / (r * r) + f * f * r * S
            znew = T(T(N) / lam)
            dz = np.max(np.abs(znew - z)); z = znew
            if dz < tol:
                break
        self.static_iterations = it + 1
        self.static_residual = float(dz)
        q = np.array([U, np.zeros_like(r), z]); p = np.array([np.zeros_like(r), -U, np.zeros_like(r)])
        self.Qh, self.Ph = T(q), T(p)
        return f, z

    def free_packet(self, k0, width, center, row):
        """Outgoing Gaussian packet (positive frequency) in row 0 (complex charged: rows 0,1) or
        neutral row 2 (real field), exact free-field initial velocity."""
        r = self.r
        u = np.exp(-(r - center) ** 2 / (2 * width ** 2) + 1j * k0 * (r - center))
        if getattr(self, 'odd_packet', True):
            # subtract the mirror image so that the odd continuation through r = 0 is smooth
            u = u - np.exp(-(r + center) ** 2 / (2 * width ** 2) - 1j * k0 * (r + center))
        uh = T(u)
        if row in (0, 1):
            vh = -1j * self.omega[0] * uh
            self.Qh = np.array([uh.real, uh.imag, np.zeros(self.n)])
            self.Ph = np.array([vh.real, vh.imag, np.zeros(self.n)])
        else:
            vh = -1j * self.omega[2] * uh
            self.Qh = np.array([np.zeros(self.n), np.zeros(self.n), uh.real])
            self.Ph = np.array([np.zeros(self.n), np.zeros(self.n), vh.real])

    def state_hash(self):
        return hashlib.sha256(np.ascontiguousarray(np.r_[self.Qh, self.Ph]).tobytes()).hexdigest()


def run_formation(args):
    sponge = None
    if args.sponge_start is not None:
        sponge = dict(start=args.sponge_start, width=args.sponge_width, gmax=args.sponge_gmax,
                      power=args.sponge_power)
    sim = Radial(args.h, args.L, args.dt, args.scheme, True, sponge)
    f, z = sim.prepared_gaussian()
    d0 = sim.diagnostics()
    E0, Q0 = d0['E'], d0['Q']
    every = int(round(args.diag / args.dt))
    hist = [d0]
    modes = []
    out = Path(args.out)
    t0 = time.time()
    snap_times = sorted(set(float(s) for s in args.modes.split(',') if s)) if args.modes else []
    nsteps_total = int(round(args.tmax / args.dt))
    done = 0
    while done < nsteps_total:
        nb = min(every, nsteps_total - done)
        sim.advance_block(nb); done += nb
        d = sim.diagnostics(); hist.append(d)
        for st in snap_times:
            if abs(sim.tau - st) < 0.5 * args.dt:
                for rmin in (240.0, 480.0):
                    m = sim.mode_numbers(rmin); m['tau'] = sim.tau; modes.append(m)
                q, p = sim.physical()
                np.savez_compressed(str(out) + f'_snap{int(st)}.npz', r=sim.r, q=q, p=p, tau=sim.tau)
        if len(hist) % 200 == 0:
            print(f"tau={sim.tau:.1f} E/E0-1={(d['E']+d['Erec'])/E0-1:.3e} Q/Q0-1={(d['Q']+d['Qrec'])/Q0-1:.3e} "
                  f"Q60/Q0={d['Q60']/Q0:.6f} wall={time.time()-t0:.0f}s", flush=True)
    cols = list(hist[0].keys())
    H = np.array([[h_[c] for c in cols] for h_ in hist])
    np.savez_compressed(str(out) + '_history.npz', H=H, cols=np.array(cols))
    q, p = sim.physical()
    np.savez_compressed(str(out) + '_final.npz', r=sim.r, q=q, p=p, tau=sim.tau)
    tau = H[:, 0]
    i60 = cols.index('Q60')
    res = dict(args=vars(args), n=sim.n, omega_max_dt=sim.wmax * args.dt,
               static_iterations=sim.static_iterations, static_residual=sim.static_residual,
               E0=E0, Q0=Q0, seconds=time.time() - t0,
               max_rel_energy_error=float(np.max(np.abs((H[:, cols.index('E')] + H[:, cols.index('Erec')]) / E0 - 1))),
               max_rel_charge_error=float(np.max(np.abs((H[:, cols.index('Q')] + H[:, cols.index('Qrec')]) / Q0 - 1))),
               max_rel_local_energy_mismatch=float(np.max(np.abs(H[:, cols.index('E_local_sum')] / H[:, cols.index('E')] - 1))),
               final_hash=sim.state_hash(), modes=modes)
    windows = [(9600, 12000), (19200, 24000), (28800, 36000), (38400, 48000)]
    for a, b in windows:
        m = (tau > a) & (tau <= b + 1e-9)
        if m.sum() > 2 and tau[-1] >= b - 1e-6:
            y = H[m, i60] / Q0
            slope = np.polyfit(tau[m], y, 1)[0]
            res[f'Q60_mean_{a}_{b}'] = float(y.mean())
            res[f'Q60le_mean_{a}_{b}'] = float(np.mean(H[m, cols.index('Q60_le')]) / Q0)
            res[f'Q60_slope_per_1000_{a}_{b}'] = float(slope * 1000)
    res['final'] = {c: float(H[-1, k]) for k, c in enumerate(cols)}
    Path(str(out) + '_RESULT.json').write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ('modes', 'final', 'args')}, indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--h', type=float, default=0.25)
    ap.add_argument('--L', type=float, default=6400.0)
    ap.add_argument('--dt', type=float, default=0.04)
    ap.add_argument('--tmax', type=float, default=12000.0)
    ap.add_argument('--scheme', default='strang')
    ap.add_argument('--diag', type=float, default=4.0)
    ap.add_argument('--modes', default='')
    ap.add_argument('--sponge-start', type=float, default=None)
    ap.add_argument('--sponge-width', type=float, default=3000.0)
    ap.add_argument('--sponge-gmax', type=float, default=0.05)
    ap.add_argument('--sponge-power', type=float, default=3.0)
    ap.add_argument('--out', required=True)
    run_formation(ap.parse_args())
