"""Independent lattice code for verifying modules/charged_exit_sim (verifier key v_charged_exit).

Written from the brief only; it does not import or copy the module's code.

Formulation (deliberately different from the module's conformal chi = a*f rescaling):
comoving-coordinate fields f in {S = phi/v, q1, q2} with X = v (q1 + i q2)/sqrt(2), and their CANONICAL
momenta p_f = a^2 f' (prime = conformal time).  Units: lengths/times in ell = 1/(sqrt(lambda_phi) v),
energies in lambda_phi v^4.  With Y = X/v, |Y|^2 = (q1^2+q2^2)/2:

  Ubar = (S^2-1)^2/4 + A S^2 |Y|^2 + B |Y|^4,   A = g_X^2/lambda_phi = 1,  B = lambda_X/lambda_phi = 0.1
  L    = a^2/2 (f'^2 - |grad f|^2) summed over f  -  a^4 Ubar
  H(t) = sum_f p_f^2/(2 a^2) + a^2 |grad f|^2/2 + a^4 Ubar
  p_S' = a^2 lap S - a^4 [S(S^2-1) + A S (q1^2+q2^2)]
  p_i' = a^2 lap q_i - a^4 q_i [A S^2 + B (q1^2+q2^2)]

Radiation era a = tau/tau_i.  Integrator: Strang splitting in extended phase space (kick at fixed tau,
exact drift q += p * int dtau/a^2 = p * tau_i^2 dt/(tau (tau+dt))).  Each sub-flow is U(1) invariant, so the
comoving charge Q = a^2 int (q2 q1' - q1 q2') d^2x = int (q2 p1 - q1 p2) d^2x is conserved to round-off.
Physical charge density n = (q2 p1 - q1 p2)/a^3 in units sqrt(lambda_phi) v^3.
"""
import numpy as np
from scipy import ndimage

A_G = 1.0
B_X = 0.1
NBAR_C = np.sqrt(8.0 * B_X) * 0.5 ** 1.5   # |Y_c|^2 = 1/(2A); n_c = sqrt(8B)|Y_c|^3


def lap(f, out, inv_dx2):
    """5-point periodic Laplacian (2D) or 7-point (3D), in place into out."""
    np.multiply(f, -2.0 * f.ndim, out=out)
    for ax in range(f.ndim):
        sl = [slice(None)] * f.ndim
        def S(s):
            sl2 = list(sl); sl2[ax] = s; return tuple(sl2)
        out[S(slice(1, None))] += f[S(slice(None, -1))]
        out[S(slice(0, 1))] += f[S(slice(-1, None))]
        out[S(slice(None, -1))] += f[S(slice(1, None))]
        out[S(slice(-1, None))] += f[S(slice(0, 1))]
    out *= inv_dx2


class Lattice:
    def __init__(self, N, dx, tau_i, ndim=2, with_X=True):
        self.N, self.dx, self.tau_i, self.ndim = N, dx, tau_i, ndim
        shp = (N,) * ndim
        self.S = np.zeros(shp); self.pS = np.zeros(shp)
        self.with_X = with_X
        self.q1 = np.zeros(shp); self.q2 = np.zeros(shp)
        self.p1 = np.zeros(shp); self.p2 = np.zeros(shp)
        self.tau = tau_i
        self._t1 = np.empty(shp); self._t2 = np.empty(shp); self._t3 = np.empty(shp)
        self.inv_dx2 = 1.0 / dx ** 2

    @property
    def a(self):
        return self.tau / self.tau_i

    def kick(self, a, h):
        """p += h * force evaluated at scale factor a (fields fixed)."""
        a2, a4 = a * a, a ** 4
        t1, t2, t3 = self._t1, self._t2, self._t3
        S = self.S
        if self.with_X:
            np.multiply(self.q1, self.q1, out=t3); t3 += self.q2 * self.q2      # F2 = q1^2+q2^2
        # phi
        lap(S, t1, self.inv_dx2)
        np.multiply(S, S, out=t2); t2 -= 1.0
        if self.with_X:
            t2 += A_G * t3
        t2 *= S
        t1 *= a2; t2 *= a4; t1 -= t2
        self.pS += h * t1
        if self.with_X:
            # mass^2 for q: A S^2 + B F2  (stored in t3)
            t3 *= B_X; t3 += A_G * S * S; t3 *= a4
            for q, p in ((self.q1, self.p1), (self.q2, self.p2)):
                lap(q, t1, self.inv_dx2)
                t1 *= a2
                np.multiply(q, t3, out=t2)
                t1 -= t2
                p += h * t1

    def drift(self, w):
        self.S += w * self.pS
        if self.with_X:
            self.q1 += w * self.p1
            self.q2 += w * self.p2

    def evolve(self, tau_end, dt):
        """KDK steps until tau_end (last step shortened to land exactly)."""
        nsteps = int(np.ceil((tau_end - self.tau) / dt - 1e-9))
        if nsteps <= 0:
            return
        h = (tau_end - self.tau) / nsteps
        ti = self.tau_i
        self.kick(self.tau / ti, 0.5 * h)
        for k in range(nsteps):
            t0 = self.tau
            t1 = t0 + h
            self.drift(ti * ti * h / (t0 * t1))
            self.tau = t1
            if k < nsteps - 1:
                self.kick(t1 / ti, h)
            else:
                self.kick(t1 / ti, 0.5 * h)

    # ---------- observables ----------
    def charge(self):
        """Comoving conserved charge per comoving cell volume, summed (brief sign)."""
        if not self.with_X:
            return 0.0
        return float(np.sum(self.q2 * self.p1 - self.q1 * self.p2)) * self.dx ** self.ndim

    def nphys(self):
        return (self.q2 * self.p1 - self.q1 * self.p2) / self.a ** 3

    def F2(self):
        return self.q1 ** 2 + self.q2 ** 2

    def save(self, path):
        np.savez(path, S=self.S, pS=self.pS, q1=self.q1, q2=self.q2, p1=self.p1, p2=self.p2,
                 tau=self.tau, tau_i=self.tau_i, dx=self.dx, N=self.N, ndim=self.ndim, with_X=self.with_X)

    @classmethod
    def load(cls, path):
        d = np.load(path)
        L = cls(int(d['N']), float(d['dx']), float(d['tau_i']), int(d['ndim']), bool(d['with_X']))
        for k in ('S', 'pS', 'q1', 'q2', 'p1', 'p2'):
            getattr(L, k)[...] = d[k]
        L.tau = float(d['tau'])
        return L


# ---------------- initial conditions ----------------

def smooth_noise(rng, N, dx, ndim, sc, std):
    """Gaussian random field: white noise times exp(-k^2 sc^2/2) in Fourier space, normalised to std."""
    w = rng.standard_normal((N,) * ndim)
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    if ndim == 2:
        k2 = k1[:, None] ** 2 + k1[None, :] ** 2
    else:
        k2 = k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + k1[None, None, :] ** 2
    f = np.real(np.fft.ifftn(np.fft.fftn(w) * np.exp(-0.5 * k2 * sc * sc)))
    f -= f.mean()
    return f * (std / f.std())


def init_charged_exit(lat, kappa, seed, phi_std=1e-2, x_std=1e-3, sc=1.0):
    """phi = 0 rotating condensate with n_i = kappa n_c at a = 1 plus noise.
    chi = a Y rotates as Y0 e^{-i Omega tau}, Omega^2 = 2 B Y0^2, Y0 = (n_i/sqrt(8B))^{1/3}."""
    rng = np.random.default_rng(seed)
    N, dx, nd = lat.N, lat.dx, lat.ndim
    n_i = kappa * NBAR_C
    Y0 = (n_i / np.sqrt(8 * B_X)) ** (1.0 / 3.0)
    Om = np.sqrt(2 * B_X) * Y0
    lat.S[...] = smooth_noise(rng, N, dx, nd, sc, phi_std)
    lat.pS[...] = 0.0
    e1 = smooth_noise(rng, N, dx, nd, sc, x_std)
    e2 = smooth_noise(rng, N, dx, nd, sc, x_std)
    c = np.sqrt(2.0) * Y0
    # chi = c (1 + e) at theta = 0 with complex e = e1 + i e2 (relative perturbation, per real component ~ x_std)
    chi1 = c * (1.0 + e1)
    chi2 = c * e2
    # co-rotating velocity chi' = -i Om chi  (chi = chi1 + i chi2)
    dchi1 = Om * chi2
    dchi2 = -Om * chi1
    a, ap = lat.a, 1.0 / lat.tau_i
    lat.q1[...] = chi1 / a; lat.q2[...] = chi2 / a
    lat.p1[...] = a * dchi1 - ap * chi1
    lat.p2[...] = a * dchi2 - ap * chi2
    return dict(n_i=n_i, Y0=Y0, Omega=Om)


def init_neutral(lat, seed, phi_std=1e-2, sc=1.0):
    rng = np.random.default_rng(seed)
    lat.S[...] = smooth_noise(rng, lat.N, lat.dx, lat.ndim, sc, phi_std)
    lat.pS[...] = 0.0
    for f in (lat.q1, lat.q2, lat.p1, lat.p2):
        f[...] = 0.0


def branch1_state(nbar):
    """Homogeneous phi != 0 rotating state: S^2 = 1 - F^2 (A=1), omega^2 = A S^2 + B F^2, n = omega F^2."""
    from scipy.optimize import brentq
    def g(F2):
        S2 = 1.0 - A_G * F2
        return np.sqrt(A_G * S2 + B_X * F2) * F2 - nbar
    F2 = brentq(g, 1e-12, 1.0 / A_G - 1e-12)   # lower root (n increasing on it)
    S0 = np.sqrt(1.0 - A_G * F2)
    om = np.sqrt(A_G * S0 ** 2 + B_X * F2)
    return F2, S0, om


def init_subcritical_plus(lat, frac_nc, seed, phi_std=1e-2, x_std=1e-3, sc=1.0):
    rng = np.random.default_rng(seed)
    N, dx, nd = lat.N, lat.dx, lat.ndim
    F2, S0, om = branch1_state(frac_nc * NBAR_C)
    lat.S[...] = S0 + smooth_noise(rng, N, dx, nd, sc, phi_std)
    a, ap = lat.a, 1.0 / lat.tau_i
    # physical S' = 0 -> p_S = a^2 S' = 0
    lat.pS[...] = 0.0
    e1 = smooth_noise(rng, N, dx, nd, sc, x_std)
    e2 = smooth_noise(rng, N, dx, nd, sc, x_std)
    c = np.sqrt(F2)
    q1 = c * (1 + e1); q2 = c * e2
    lat.q1[...] = q1; lat.q2[...] = q2
    # physical rotation dq/dt = -i om q  -> q' = a dq/dt ; p = a^2 q'
    lat.p1[...] = a ** 3 * om * q2
    lat.p2[...] = -a ** 3 * om * q1
    return dict(F2=F2, S0=S0, omega=om, nbar=frac_nc * NBAR_C)


# ---------------- diagnostics ----------------

def periodic_label(mask):
    """4-connected (2D) / 6-connected (3D) labelling on a torus."""
    lab, n = ndimage.label(mask)
    if n == 0:
        return lab, 0
    parent = np.arange(n + 1)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    pairs = []
    for ax in range(mask.ndim):
        a0 = np.take(lab, 0, axis=ax).ravel(); a1 = np.take(lab, -1, axis=ax).ravel()
        m = (a0 > 0) & (a1 > 0)
        if m.any():
            pairs.append(np.unique(np.stack([a0[m], a1[m]], 1), axis=0))
    for pr in pairs:
        for u, v in pr:
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[max(ru, rv)] = min(ru, rv)
    root = np.array([find(i) for i in range(n + 1)])
    uniq, newid = np.unique(root, return_inverse=True)
    lab = newid[lab]
    return lab, len(uniq) - 1


def edt_periodic(nonzero_mask, pad, return_indices=False):
    """Euclidean distance (in cells) to nearest False cell on a torus (approximate: wrap padding of 'pad' cells)."""
    m = np.pad(nonzero_mask, pad, mode='wrap')
    if return_indices:
        d, idx = ndimage.distance_transform_edt(m, return_indices=True)
        sl = tuple(slice(pad, pad + s) for s in nonzero_mask.shape)
        idx = tuple(ix[sl] for ix in idx)
        return d[sl], idx, m.shape
    d = ndimage.distance_transform_edt(m)
    sl = tuple(slice(pad, pad + s) for s in nonzero_mask.shape)
    return d[sl]


def link_count(lbl):
    """Number of nearest-neighbour links (periodic) across which the +/-1 label changes."""
    c = 0
    for ax in range(lbl.ndim):
        c += int(np.count_nonzero(lbl != np.roll(lbl, -1, axis=ax)))
    return c


def diagnostics(lat, qexcl_ell=1.0, pad=96, radii=(2.0, 4.0, 6.0)):
    a = lat.a
    S = lat.S
    nd = lat.ndim
    dx = lat.dx
    ncell = S.size
    out = dict(tau=lat.tau, a=a)
    plus = S > 0.5
    minus = S < -0.5
    out['f_plus'] = float(plus.mean()); out['f_minus'] = float(minus.mean())
    fp, fm = out['f_plus'], out['f_minus']
    out['p_plus'] = fp / (fp + fm) if fp + fm > 0 else float('nan')
    out['S_mean'] = float(S.mean())
    # raw sign-change length (contaminated by any sign oscillation), Crofton factor for isotropic curves
    sg = np.where(S >= 0, 1, -1).astype(np.int8)
    crof = (np.pi / 4.0) if nd == 2 else (2.0 / 3.0)
    Lraw = crof * link_count(sg) / (ncell * dx)
    out['L_raw'] = Lraw
    if lat.with_X:
        F2 = lat.F2()
        n = lat.nphys()
        Qtot = float(n.sum())
        out['Q_comoving'] = lat.charge()
        out['F2_mean'] = float(F2.mean())
        qm = F2 >= 1.0
        out['f_Qmatter'] = float(qm.mean())
    else:
        F2 = None; n = None; Qtot = 0.0; qm = np.zeros(S.shape, bool)
    # ---- vacuum-basin interface with Q-matter excluded from seeds ----
    if qm.any():
        dq = edt_periodic(~qm, pad)          # distance (cells) to nearest Q-matter cell
        nearQ = dq * dx * a <= qexcl_ell     # within qexcl_ell physical of Q-matter
    else:
        nearQ = np.zeros(S.shape, bool)
    seedp = plus & ~nearQ
    seedm = minus & ~nearQ
    seeds = seedp | seedm
    out['f_seed_plus'] = float(seedp.mean()); out['f_seed_minus'] = float(seedm.mean())
    res = {}
    for tag, sp, sm in (('basin', seedp, seedm), ('basin_noexcl', plus, minus)):
        sd = sp | sm
        if sp.any() and sm.any():
            d, idx, _ = edt_periodic(~sd, pad, return_indices=True)
            # map padded indices back to torus
            ii = [(ix - pad) % s for ix, s in zip(idx, S.shape)]
            basin = np.where(sp[tuple(ii)], 1, -1).astype(np.int8)
            nlinks = link_count(basin)
            res[tag] = (crof * nlinks / (ncell * dx), basin)
        else:
            res[tag] = (0.0, None)
    out['L_basin'] = res['basin'][0]
    out['L_basin_noexcl'] = res['basin_noexcl'][0]
    basin = res['basin'][1]
    # vacuum cluster statistics (periodic)
    for tag, m in (('plus', plus), ('minus', minus)):
        lab, nc = periodic_label(m)
        if nc > 0:
            sizes = np.bincount(lab.ravel())[1:]
            out['ncl_' + tag] = int(nc)
            out['maxcl_' + tag] = float(sizes.max() / ncell)
        else:
            out['ncl_' + tag] = 0; out['maxcl_' + tag] = 0.0
    # ---- charge near walls ----
    if lat.with_X and basin is not None and Qtot != 0:
        iface = np.zeros(S.shape, bool)
        for ax in range(nd):
            ch = basin != np.roll(basin, -1, axis=ax)
            iface |= ch | np.roll(ch, 1, axis=ax)
        dI = edt_periodic(~iface, pad) * dx * a   # physical distance to interface
        for r in radii:
            out['qfrac_within_%g' % r] = float(n[dI <= r].sum() / Qtot)
            out['area_within_%g' % r] = float((dI <= r).mean())
        # component definition: core = |S|<=0.5 or F2>=1 ; wall-type if adjacent to both seed signs
        core = (np.abs(S) <= 0.5) | qm
        lab, nc = periodic_label(core)
        if nc > 0:
            def adj(mm):
                dd = mm.copy()
                for ax in range(nd):
                    dd |= np.roll(mm, 1, axis=ax) | np.roll(mm, -1, axis=ax)
                return dd
            lp = np.unique(lab[adj(seedp) & core]); lm = np.unique(lab[adj(seedm) & core])
            wl = np.intersect1d(lp, lm); wl = wl[wl > 0]
            wallmask = np.isin(lab, wl)
            out['qfrac_wallcomp'] = float(n[wallmask].sum() / Qtot)
            out['qfrac_core'] = float(n[core].sum() / Qtot)
            out['qfrac_Qmatter'] = float(n[qm].sum() / Qtot)
            out['qfrac_Qmatter_on_wallcomp'] = float(n[qm & wallmask].sum() / Qtot)
            out['f_wallcomp'] = float(wallmask.mean())
        out['qfrac_vac'] = float(n[(np.abs(S) > 0.5) & ~qm].sum() / Qtot)
    return out
