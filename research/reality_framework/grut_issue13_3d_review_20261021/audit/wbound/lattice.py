#!/usr/bin/env python3
"""Classical lattice evolution for the w-bound audit.  numpy only, deterministic.  Units v = 1.

Fields: real phi, complex X = (x1 + i x2)/sqrt(2); D = 0 (invariant classical slice).
U = lam/4 (phi^2-1)^2 + g^2 phi^2 |X|^2 + lam_X |X|^4, benchmark lam = g^2 = 0.01, lam_X = 0.001.
dim = 1 means planar-symmetric 3D configurations (fields depend on z only); the isotropic pressure is still
p = K - Grad/3 - U with K = (phidot^2 + xdot^2)/2, Grad = (grad phi)^2/2 + (grad x)^2/2 (physical).

Background: flat (a = 1) or prescribed radiation a(t) = (t/t0)^(1/2) in cosmic time.  Canonical momenta per
comoving volume P_f = a^3 fdot; KDK leapfrog; the lattice gradient energy (forward differences) is paired with the
(2d+1)-point Laplacian, so sum f Lap f = -2 sum Grad exactly and the scalar virial identity
    d/dt [a^3 <delta . fdot>] = a^3 <fdot^2 - (grad f)^2/a^2 - delta . dU/df>
holds on the lattice in continuous time (O(dt^2) in KDK).
Recorded every `every` steps: rho, p, <Phi_0> = lam <phi^2 - 1>, g0 = <phi phidot + x . xdot>, phi extrema, the
fraction of sites with |phi| > 1, and (3D flat, optional) the dilatation moment Dm = sum (x - x0) . T^{0i} dV.
"""
import math

import numpy as np

LAM, G2, LX = 0.01, 0.01, 0.001


class Lattice:
    def __init__(self, dim, N, dx, frw=None, with_X=True):
        self.dim, self.N, self.dx = dim, N, dx
        shape = (N,) * dim
        self.phi = np.ones(shape)
        self.x1 = np.zeros(shape)
        self.x2 = np.zeros(shape)
        self.Pphi = np.zeros(shape)
        self.P1 = np.zeros(shape)
        self.P2 = np.zeros(shape)
        self.with_X = with_X
        self.frw = frw            # None or dict(t0=..): a = (t/t0)^(1/2), t starts at t0
        self.t = frw["t0"] if frw else 0.0

    # ---- background ----
    def a(self, t=None):
        t = self.t if t is None else t
        return math.sqrt(t / self.frw["t0"]) if self.frw else 1.0

    def H(self, t=None):
        t = self.t if t is None else t
        return 0.5 / t if self.frw else 0.0

    # ---- lattice operators ----
    def lap(self, f):
        out = -2.0 * self.dim * f
        for ax in range(self.dim):
            out = out + np.roll(f, 1, ax) + np.roll(f, -1, ax)
        return out / self.dx**2

    def grad_energy(self, f):
        s = np.zeros_like(f)
        for ax in range(self.dim):
            d = (np.roll(f, -1, ax) - f) / self.dx
            s += 0.5 * d * d
        return s

    # ---- potential ----
    def U(self):
        X2 = 0.5 * (self.x1**2 + self.x2**2)
        return LAM / 4 * (self.phi**2 - 1) ** 2 + G2 * self.phi**2 * X2 + LX * X2**2

    def dU(self):
        X2 = 0.5 * (self.x1**2 + self.x2**2)
        dphi = LAM * self.phi * (self.phi**2 - 1) + 2 * G2 * self.phi * X2
        dX2 = G2 * self.phi**2 + 2 * LX * X2
        return dphi, dX2 * self.x1, dX2 * self.x2

    # ---- evolution ----
    def kick(self, h, t):
        a = self.a(t)
        dphi, d1, d2 = self.dU()
        self.Pphi += h * (a * self.lap(self.phi) - a**3 * dphi)
        if self.with_X:
            self.P1 += h * (a * self.lap(self.x1) - a**3 * d1)
            self.P2 += h * (a * self.lap(self.x2) - a**3 * d2)

    def drift(self, h, t_mid):
        a3 = self.a(t_mid) ** 3
        self.phi += h * self.Pphi / a3
        if self.with_X:
            self.x1 += h * self.P1 / a3
            self.x2 += h * self.P2 / a3

    def step(self, dt):
        t = self.t
        self.kick(0.5 * dt, t)
        self.drift(dt, t + 0.5 * dt)
        self.kick(0.5 * dt, t + dt)
        self.t = t + dt

    # ---- diagnostics ----
    def diag(self, center=None):
        a = self.a()
        a3 = a**3
        fd = self.Pphi / a3
        x1d, x2d = self.P1 / a3, self.P2 / a3
        K = 0.5 * (fd**2 + x1d**2 + x2d**2)
        Gr = (self.grad_energy(self.phi) + self.grad_energy(self.x1) + self.grad_energy(self.x2)) / a**2
        Uu = self.U()
        rho = float(np.mean(K + Gr + Uu))
        p = float(np.mean(K - Gr / 3 - Uu))
        d = dict(t=self.t, a=a, H=self.H(), rho=rho, p=p,
                 Phi0=float(LAM * np.mean(self.phi**2 - 1)),
                 g0=float(np.mean(self.phi * fd + self.x1 * x1d + self.x2 * x2d)),
                 gv=float(np.mean((self.phi - 1) * fd + self.x1 * x1d + self.x2 * x2d)),
                 phi_min=float(self.phi.min()), phi_max=float(self.phi.max()),
                 frac_abs_gt1=float(np.mean(np.abs(self.phi) > 1)),
                 frac_neg=float(np.mean(self.phi < 0)),
                 rho_loc_max=float((K + Gr + Uu).max()),
                 X2mean=float(np.mean(0.5 * (self.x1**2 + self.x2**2))))
        if center is not None and self.dim == 3 and not self.frw:
            # dilatation moment Dm = sum_x (x - c) . T^{0i} dV, T^{0i} = -sum_f fdot d_i f (centered differences)
            N, dx = self.N, self.dx
            Dm = 0.0
            for ax in range(3):
                coord = (np.arange(N) - center[ax] + N / 2) % N - N / 2      # minimal image
                shp = [1, 1, 1]
                shp[ax] = N
                coord = coord.reshape(shp) * dx
                Ti = np.zeros_like(self.phi)
                for f, fdot in ((self.phi, fd), (self.x1, x1d), (self.x2, x2d)):
                    Ti -= fdot * (np.roll(f, -1, ax) - np.roll(f, 1, ax)) / (2 * dx)
                Dm += float(np.sum(coord * Ti)) * dx**3
            d["Dm"] = Dm
            d["E_tot"] = rho * (N * dx) ** 3
        return d


def gaussian_field(rng, dim, N, dx, kc, rms, kmin=0.0):
    """Real Gaussian random field with spectrum ~ exp(-k^2/(2 kc^2)) (k > kmin), zero mean, given rms."""
    shape = (N,) * dim
    noise = rng.normal(size=shape)
    fk = np.fft.rfftn(noise) if dim > 1 else np.fft.rfft(noise)
    ks = [np.fft.fftfreq(N, d=dx) * 2 * np.pi for _ in range(dim - 1)] + [np.fft.rfftfreq(N, d=dx) * 2 * np.pi]
    grids = np.meshgrid(*ks, indexing="ij")
    k2 = sum(gk**2 for gk in grids)
    filt = np.exp(-k2 / (4 * kc**2)) * (k2 > kmin**2)
    f = np.fft.irfftn(fk * filt, s=shape) if dim > 1 else np.fft.irfft(fk * filt, n=N)
    f -= f.mean()
    sd = f.std()
    return f * (rms / sd) if sd > 0 else f


def run(lat, dt, T, every=1, center=None, callback=None):
    rows = [lat.diag(center)]
    nsteps = int(round(T / dt))
    for i in range(1, nsteps + 1):
        lat.step(dt)
        if i % every == 0 or i == nsteps:
            rows.append(lat.diag(center))
            if callback:
                callback(lat, rows)
    return rows


def _trap(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, "trapezoid") else float(np.trapz(y, x))


def summarize_flat(rows, t_from=None, t_to=None):
    """Flat-space time averages over [t_from, t_to]; checks <3p-rho> = <Phi_0> + (g0(T)-g0(0))/T."""
    t = np.array([r["t"] for r in rows])
    sel = np.ones_like(t, bool)
    if t_from is not None:
        sel &= t >= t_from
    if t_to is not None:
        sel &= t <= t_to
    t = t[sel]
    rr = [r for r, s in zip(rows, sel) if s]
    rho = np.array([r["rho"] for r in rr])
    p = np.array([r["p"] for r in rr])
    Phi0 = np.array([r["Phi0"] for r in rr])
    g0 = np.array([r["g0"] for r in rr])
    T = t[-1] - t[0]
    av = lambda y: _trap(y, t) / T
    lhs = av(3 * p - rho)
    rhs = av(Phi0)
    bnd = (g0[-1] - g0[0]) / T
    return dict(T=float(T), wbar=av(p) / av(rho), rho_mean=av(rho), rho_drift=float(rho.std() / rho.mean()),
                avg_3p_minus_rho=lhs, avg_Phi0=rhs, boundary=bnd,
                identity_residual_rel=float(abs(lhs - rhs - bnd) / av(rho)),
                w_inst_first=float(p[0] / rho[0]), w_inst_max=float((p / rho).max()),
                phi_min=float(min(r["phi_min"] for r in rr)), phi_max=float(max(r["phi_max"] for r in rr)),
                one_sign=bool(min(r["phi_min"] for r in rr) > 0))


def running_wbar(rows):
    """Cumulative wbar_T = int_0^T p / int_0^T rho (flat space)."""
    t = np.array([r["t"] for r in rows])
    rho = np.array([r["rho"] for r in rows])
    p = np.array([r["p"] for r in rows])
    cp = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(t))])
    cr = np.concatenate([[0], np.cumsum(0.5 * (rho[1:] + rho[:-1]) * np.diff(t))])
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(cr > 0, cp / cr, p / rho)
    w[0] = p[0] / rho[0]
    return t, w
