#!/usr/bin/env python3
"""Independent lattice code (verifier v_network) for the real Z2 scalar
    U(phi) = (lam/4)(phi^2 - v^2)^2
in a radiation background, PRS form
    phi'' + alpha (a'/a) phi' - lap phi + a^beta U'(phi) = 0,   a = tau (a'/a = 1/tau).

Written from scratch by the verifier; it does not import or copy the module code.

Units (verifier's own choice, different from the module): v = 1, lam = 2, so the planar
kink is phi = tanh(z/w) with w = sqrt(2/lam)/v = 1, m^2 = U''(1) = 2 lam v^2 = 4,
sigma = (2 sqrt(2 lam)/3) v^3 = 4/3.  (Module units: length unit 1/m_phi = w/2, so
tau_module = 2 tau_here.)

Integrator: exact treatment of the friction term through P = tau^alpha phi'
    phi' = P / tau^alpha,   P' = tau^alpha (lap phi - a^beta U'(phi)),
staggered leapfrog (P at half steps).  This differs from the PRS (1-delta)/(1+delta)
discretisation used by the module but is consistent to O(dt^2).

Wall-size estimators (all comoving):
  raw   : number of sign-changing nearest-neighbour links x dx^(d-1)  ("Manhattan" measure)
  iso   : raw x pi/4 (2D) or raw x 2/3 (3D)  (orientation-averaged correction)
  grad  : per link weight |grad phi| / sum_i |d_i phi| (link-direction derivative = link
          difference, transverse derivatives = mean of the central differences at the link ends)
  ms    : 2D only, marching-squares contour length with linear interpolation of zero crossings
          (saddles resolved with the cell-centre average)
"""
import numpy as np

LAM = 2.0
SIGMA = 4.0 / 3.0


def dU(phi):
    return LAM * phi * (phi * phi - 1.0)


def U(phi):
    return 0.25 * LAM * (phi * phi - 1.0) ** 2


def laplacian(f, dx, out):
    d = f.ndim
    out.fill(0.0)
    for ax in range(d):
        # out += f(x+1) + f(x-1) along ax with periodic wrap, using slices
        sl = [slice(None)] * d
        sr = [slice(None)] * d
        # interior
        o = [slice(None)] * d
        o[ax] = slice(1, -1)
        sl[ax] = slice(2, None)
        sr[ax] = slice(None, -2)
        out[tuple(o)] += f[tuple(sl)] + f[tuple(sr)]
        # boundaries
        o0 = [slice(None)] * d; o0[ax] = 0
        a1 = [slice(None)] * d; a1[ax] = 1
        am = [slice(None)] * d; am[ax] = -1
        out[tuple(o0)] += f[tuple(a1)] + f[tuple(am)]
        oe = [slice(None)] * d; oe[ax] = -1
        a0 = [slice(None)] * d; a0[ax] = 0
        a2 = [slice(None)] * d; a2[ax] = -2
        out[tuple(oe)] += f[tuple(a0)] + f[tuple(a2)]
    out -= 2 * d * f
    out *= 1.0 / (dx * dx)
    return out


def gaussian_ic(N, d, dx, s_c, std, seed, mean_shift_frac=None):
    """Gaussian random field: white noise filtered with exp(-k^2 s_c^2 / 2) in Fourier space
    (equivalent to real-space convolution with a Gaussian kernel of std s_c), zero mean, given std.
    If mean_shift_frac = p is given, a constant is subtracted so that exactly a fraction p of the
    cells is positive."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((N,) * d)
    fk = np.fft.rfftn(noise)
    ks = [2 * np.pi * np.fft.fftfreq(N, d=dx)] * (d - 1) + [2 * np.pi * np.fft.rfftfreq(N, d=dx)]
    grids = np.meshgrid(*ks, indexing="ij")
    k2 = sum(g * g for g in grids)
    fk *= np.exp(-0.5 * k2 * s_c * s_c)
    f = np.fft.irfftn(fk, s=(N,) * d)
    f -= f.mean()
    f *= std / f.std()
    if mean_shift_frac is not None:
        q = np.quantile(f, 1.0 - mean_shift_frac)
        f -= q
    return f


def link_masks(phi):
    """Boolean masks of sign-changing links in each direction (link between x and x+e_i)."""
    s = phi > 0
    return [s != np.roll(s, -1, axis=ax) for ax in range(phi.ndim)]


def wall_measures(phi, dx):
    d = phi.ndim
    masks = link_masks(phi)
    raw_count = sum(int(m.sum()) for m in masks)
    raw = raw_count * dx ** (d - 1)
    iso = raw * (np.pi / 4 if d == 2 else 2.0 / 3.0)
    # gradient-weighted estimator
    cen = [(np.roll(phi, -1, axis=ax) - np.roll(phi, 1, axis=ax)) / (2 * dx) for ax in range(d)]
    gsum = 0.0
    for ax in range(d):
        m = masks[ax]
        if not m.any():
            continue
        comps = []
        for bx in range(d):
            if bx == ax:
                comps.append(((np.roll(phi, -1, axis=ax) - phi) / dx)[m])
            else:
                comps.append(0.5 * (cen[bx] + np.roll(cen[bx], -1, axis=ax))[m])
        comps = np.array(comps)
        g = np.sqrt((comps ** 2).sum(axis=0))
        l1 = np.abs(comps).sum(axis=0)
        gsum += float((g / np.maximum(l1, 1e-300)).sum())
    grad = gsum * dx ** (d - 1)
    out = {"raw": raw, "iso": iso, "grad": grad}
    if d == 2:
        out["ms"] = marching_squares_length(phi, dx)
    return out


def marching_squares_length(phi, dx):
    f00 = phi
    f10 = np.roll(phi, -1, axis=0)
    f01 = np.roll(phi, -1, axis=1)
    f11 = np.roll(f10, -1, axis=1)
    s00, s10, s01, s11 = f00 > 0, f10 > 0, f01 > 0, f11 > 0
    # edges: B (y=0, from 00 to 10), R (x=1, from 10 to 11), T (y=1, from 01 to 11), L (x=0, from 00 to 01)
    mB, mR, mT, mL = s00 != s10, s10 != s11, s01 != s11, s00 != s01
    with np.errstate(divide="ignore", invalid="ignore"):
        tB = f00 / (f00 - f10)
        tR = f10 / (f10 - f11)
        tT = f01 / (f01 - f11)
        tL = f00 / (f00 - f01)
    # local coordinates (x along axis 0, y along axis 1) of crossing points
    P = {"B": (tB, 0.0), "R": (1.0, tR), "T": (tT, 1.0), "L": (0.0, tL)}
    M = {"B": mB, "R": mR, "T": mT, "L": mL}
    ncross = mB.astype(np.int8) + mR + mT + mL
    two = ncross == 2
    four = ncross == 4

    def dist(e1, e2, mask):
        x1, y1 = P[e1]
        x2, y2 = P[e2]
        x1 = x1[mask] if np.ndim(x1) else x1
        y1 = y1[mask] if np.ndim(y1) else y1
        x2 = x2[mask] if np.ndim(x2) else x2
        y2 = y2[mask] if np.ndim(y2) else y2
        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    total = 0.0
    edges = ["B", "R", "T", "L"]
    for i in range(4):
        for j in range(i + 1, 4):
            m = two & M[edges[i]] & M[edges[j]]
            if m.any():
                total += float(dist(edges[i], edges[j], m).sum())
    if four.any():
        c = 0.25 * (f00 + f10 + f01 + f11)
        same = four & ((c > 0) == s00)
        other = four & ~((c > 0) == s00)
        if same.any():  # 00 and 11 connected through centre: cut corners 10 and 01
            total += float(dist("B", "R", same).sum() + dist("T", "L", same).sum())
        if other.any():  # cut corners 00 and 11
            total += float(dist("B", "L", other).sum() + dist("R", "T", other).sum())
    return total * dx


class Sim:
    def __init__(self, phi0, dx, alpha, beta, tau0, dphi0=None):
        self.phi = phi0.astype(np.float64).copy()
        self.d = phi0.ndim
        self.dx = dx
        self.alpha = alpha
        self.beta = beta
        self.tau = tau0
        self.lap = np.empty_like(self.phi)
        dphi0 = np.zeros_like(self.phi) if dphi0 is None else dphi0
        self.P = (tau0 ** alpha) * dphi0  # P = tau^alpha phi' (will be offset to half step at first step)
        self.half_started = False

    def force(self):
        laplacian(self.phi, self.dx, self.lap)
        a_b = self.tau ** self.beta if self.beta != 0 else 1.0
        # in-place: lap -= a_b * LAM * phi * (phi^2 - 1)
        if not hasattr(self, "_tmp"):
            self._tmp = np.empty_like(self.phi)
        t = self._tmp
        np.multiply(self.phi, self.phi, out=t)
        t -= 1.0
        t *= self.phi
        t *= a_b * LAM
        self.lap -= t
        return self.lap

    def step(self, dt):
        if not self.half_started:
            # first half kick: P(tau0 + dt/2) = P(tau0) + dt/2 tau0^alpha F
            F = self.force()
            self.P += 0.5 * dt * (self.tau ** self.alpha) * F
            self.half_started = True
        # drift with P at tau + dt/2
        th = self.tau + 0.5 * dt
        if not hasattr(self, "_tmp2"):
            self._tmp2 = np.empty_like(self.phi)
        np.multiply(self.P, dt / th ** self.alpha, out=self._tmp2)
        self.phi += self._tmp2
        self.tau += dt
        F = self.force()
        np.multiply(F, dt * (self.tau ** self.alpha), out=self._tmp2)
        self.P += self._tmp2  # now P at tau + dt/2 (staggered)

    def velocity(self, dt):
        # phi' at integer time approx from P at tau + dt/2 and tau - dt/2: use P(tau+dt/2) - dt/2 tau^a F
        F = self.lap  # force at current tau already computed at end of step
        Pn = self.P - 0.5 * dt * (self.tau ** self.alpha) * F
        return Pn / self.tau ** self.alpha


def wall_v2(phi, dphi, dx):
    """v^2 estimator near walls: sum phi'^2 / sum |grad phi|^2 over |phi| < 0.5 (rough)."""
    m = np.abs(phi) < 0.5
    if not m.any():
        return float("nan")
    g2 = np.zeros_like(phi)
    for ax in range(phi.ndim):
        g2 += ((np.roll(phi, -1, axis=ax) - np.roll(phi, 1, axis=ax)) / (2 * dx)) ** 2
    return float((dphi[m] ** 2).sum() / g2[m].sum())
