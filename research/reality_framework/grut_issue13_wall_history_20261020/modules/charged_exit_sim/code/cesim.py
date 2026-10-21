"""Classical lattice simulation of the charged exit of the unchanged Candidate 1 scalar sector (phi real, X complex,
D = 0) in a radiation-dominated FRW background.  numpy/scipy only.

Units (charge_phase conventions): field in v, lengths and times in ell = 1/(sqrt(lambda_phi) v) = 1/(g_X v)
= 1/(100 GeV) (comoving lengths equal physical lengths at a = 1), so with S = phi/v, Y = X/v

    Ubar = (S^2 - 1)^2/4 + A_G S^2 |Y|^2 + B_X |Y|^4,   A_G = g_X^2/lambda_phi = 1,  B_X = lambda_X/lambda_phi = 0.1,

in units lambda_phi v^4.  Kinetic terms -1/2 (dS)^2 - |dY|^2.  Charge density unit sqrt(lambda_phi) v^3 = 0.1 v^3;
energy per area unit sqrt(lambda_phi) v^3 ell^... (see NOTES.md).

Background: ds^2 = a^2(tau) (-dtau^2 + dx^2), radiation era a = tau / tau_i (a = 1 at the start, a'' = 0).
The requested comoving equations  f'' + 2 (a'/a) f' - lap f + a^2 dU/df* = 0  are solved EXACTLY in the conformal
variables chi_S = a S, chi_Y = a Y (because a'' = 0):

    chi_S'' = lap chi_S - chi_S (chi_S^2 - a^2 + 2 A_G |chi_Y|^2)
    chi_Y'' = lap chi_Y - chi_Y (A_G chi_S^2 + 2 B_X |chi_Y|^2)

i.e. flat-space dynamics with Hamiltonian H(tau) = sum dx^d [ pS^2/2 + |grad chi_S|^2/2 + |pY|^2 + |grad chi_Y|^2
+ Ut(chi, tau) ],  Ut = (chi_S^2 - a^2)^2/4 + A_G chi_S^2 |chi_Y|^2 + B_X |chi_Y|^4,  pS = chi_S', pY = chi_Y'.
Energy budget: dH/dtau = dUt/dtau = - (a a') sum dx^d (chi_S^2 - a^2),  a' = 1/tau_i.
Comoving charge: Q = a^2 int 2 (y_r y_i' - y_i y_r') d^dx = int 2 (chi_r pi_i - chi_i pi_r) d^dx  (exactly conserved;
the brief's form a^2 int (q2 q1' - q1 q2') with X = (q1 + i q2)/sqrt2 equals -Q: sign convention only).
Physical charge density n = 2 (chi_r pi_i - chi_i pi_r) / a^3  [units sqrt(lambda_phi) v^3]; for a z-invariant
2D lattice the comoving charge per comoving z-length is conserved, so n dilutes as a^-3.

Integrator: kick-drift-kick (velocity Verlet) for the time-dependent Hamiltonian (Strang splitting; symplectic;
each sub-flow is U(1)-invariant so Q is conserved to round-off), periodic (2d+1)-point Laplacian.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

A_G = 1.0          # g_X^2 / lambda_phi
B_X = 0.1          # lambda_X / lambda_phi
NBAR_C = math.sqrt(8 * B_X) * (1.0 / (2 * A_G)) ** 1.5   # n_c / (sqrt(lam) v^3) = 0.316228 (charge_phase)
NBAR_S1 = 2 * A_G / (3 * (A_G ** 2 - B_X)) * math.sqrt(A_G - (A_G ** 2 - B_X) * 2 * A_G / (3 * (A_G ** 2 - B_X)))
NBAR_STAR = (1.0 / B_X) ** 0.25 * 1.0   # n_* = sqrt(lam) v^3 (lam_phi/lam_X)^(1/4) -> 1.778 in nbar units
SIGMA_WALL = 2 * math.sqrt(2) / 3       # tree wall tension in sqrt(lam) v^3 units (= 0.9428)
WALL_HALF_WIDTH = math.sqrt(2.0)        # phi = tanh(z / delta), delta = sqrt(2) ell


# ------------------------------------------------------------------------------------------------------------
# homogeneous initial states
# ------------------------------------------------------------------------------------------------------------
def branch2_state(nbar):
    """phi = 0 branch (equilibrium rotating condensate): |Y|^2 = (n^2/(8 B))^(1/3), omega = n/(2|Y|^2)."""
    y2 = (nbar ** 2 / (8 * B_X)) ** (1.0 / 3.0)
    om = nbar / (2 * y2)
    return 0.0, math.sqrt(y2), om


def branch1_state(nbar):
    """phi != 0 branch at fixed density (charge_phase part 1): S^2 = 1 - A x, Omega^2 = A - (A^2-B) x, x = 2|Y|^2,
    nbar = x Omega.  Small-x root (stable-radial part, x < x_s1)."""
    from scipy.optimize import brentq
    xs1 = 2 * A_G / (3 * (A_G ** 2 - B_X))
    f = lambda x: x * math.sqrt(A_G - (A_G ** 2 - B_X) * x) - nbar
    x = brentq(f, 1e-14, xs1)
    S0 = math.sqrt(max(0.0, 1 - A_G * x))
    om = math.sqrt(A_G - (A_G ** 2 - B_X) * x)
    return S0, math.sqrt(x / 2), om


def fixed_vacuum_state(nbar, S0=1.0):
    """Rotating X condensate in a fixed phi = S0 background: omega^2 = A S0^2 + 2 B |Y|^2, nbar = 2 omega |Y|^2."""
    from scipy.optimize import brentq
    f = lambda y2: 2 * math.sqrt(A_G * S0 ** 2 + 2 * B_X * y2) * y2 - nbar
    y2 = brentq(f, 0.0, 100.0)
    om = math.sqrt(A_G * S0 ** 2 + 2 * B_X * y2)
    return S0, math.sqrt(y2), om


# ------------------------------------------------------------------------------------------------------------
# Gaussian random fields (same continuum field on nested grids, for resolution checks)
# ------------------------------------------------------------------------------------------------------------
def gaussian_field(d, N, L, s_c, amp, seed, N_base=None):
    """Zero-mean Gaussian random field on an N^d periodic grid of side L: white noise drawn on a base grid N_base,
    smoothed by a Gaussian kernel of standard deviation s_c (spectrum P(k) ~ exp(-k^2 s_c^2), correlation
    C(r) ~ exp(-r^2/(4 s_c^2))), normalised to standard deviation amp ON THE BASE GRID, and band-limited
    zero-padded to N (identical continuum field for every N that is a multiple of N_base)."""
    if N_base is None:
        N_base = N
    assert N % N_base == 0
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((N_base,) * d)
    fk = np.fft.fftn(white)
    dxb = L / N_base
    k1 = 2 * np.pi * np.fft.fftfreq(N_base, d=dxb)
    k2 = np.zeros((N_base,) * d)
    for ax in range(d):
        shp = [1] * d
        shp[ax] = N_base
        k2 = k2 + k1.reshape(shp) ** 2
    fk *= np.exp(-0.5 * k2 * s_c ** 2)
    # remove Nyquist planes so that the zero-padded embedding is exactly real
    if N_base % 2 == 0:
        for ax in range(d):
            sl = [slice(None)] * d
            sl[ax] = N_base // 2
            fk[tuple(sl)] = 0.0
    fk.flat[0] = 0.0
    base = np.real(np.fft.ifftn(fk))
    norm = amp / base.std()
    if N == N_base:
        return base * norm
    # zero-pad in Fourier space
    big = np.zeros((N,) * d, dtype=complex)
    h = N_base // 2
    idx_small = np.r_[0:h, N_base - h:N_base] if True else None
    idx_big = np.r_[0:h, N - h:N]
    big[np.ix_(*([idx_big] * d))] = fk[np.ix_(*([idx_small] * d))]
    big *= (N / N_base) ** d
    out = np.real(np.fft.ifftn(big)) * norm
    return out


# ------------------------------------------------------------------------------------------------------------
# periodic connected components with wrapping (percolation on the torus)
# ------------------------------------------------------------------------------------------------------------
def plabel(mask, want_wrap=True):
    """Connected components (nearest-neighbour connectivity) of a boolean mask on a periodic lattice.
    Returns labels (0 = background, 1..n), n, sizes (length n), wraps (n x d bool: component winds around
    the torus along that axis, i.e. a non-contractible, box-spanning cluster), spans_proj (n bool: projection
    onto some axis covers all N coordinates, the network_sim convention)."""
    d = mask.ndim
    N = mask.shape[0]
    st = ndimage.generate_binary_structure(d, 1)
    lab, n0 = ndimage.label(mask, structure=st)
    if n0 == 0:
        return np.zeros(mask.shape, np.int32), 0, np.zeros(0, np.int64), np.zeros((0, d), bool), np.zeros(0, bool)
    ea, eb, eax = [], [], []
    for ax in range(d):
        a_ = np.take(lab, N - 1, axis=ax).ravel()
        b_ = np.take(lab, 0, axis=ax).ravel()
        sel = (a_ > 0) & (b_ > 0)
        if sel.any():
            pr = np.unique(np.stack([a_[sel], b_[sel]], 1), axis=0)
            ea.append(pr[:, 0]); eb.append(pr[:, 1]); eax.append(np.full(pr.shape[0], ax))
    if ea:
        ea = np.concatenate(ea); eb = np.concatenate(eb); eax = np.concatenate(eax)
    else:
        ea = eb = eax = np.zeros(0, np.int64)
    g = coo_matrix((np.ones(ea.size, np.int8), (ea, eb)), shape=(n0 + 1, n0 + 1))
    ncomp, comp = connected_components(g, directed=False)
    used = np.unique(comp[1:])
    remap = np.zeros(ncomp, np.int32)
    remap[used] = np.arange(1, used.size + 1, dtype=np.int32)
    lut = remap[comp]
    lut[0] = 0
    labels = lut[lab]
    n = int(used.size)
    sizes = np.bincount(labels.ravel(), minlength=n + 1)[1:]
    wraps = np.zeros((n, d), bool)
    if want_wrap and ea.size:
        # BFS on the boundary graph with integer shift potentials; a conflict = winding
        adj = {}
        for u, v_, ax in zip(ea.tolist(), eb.tolist(), eax.tolist()):
            adj.setdefault(u, []).append((v_, ax, +1))
            adj.setdefault(v_, []).append((u, ax, -1))
        pot = {}
        for start in adj:
            if start in pot:
                continue
            pot[start] = (0,) * d
            stack = [start]
            while stack:
                u = stack.pop()
                pu = pot[u]
                for v_, ax, sgn in adj[u]:
                    pv = list(pu)
                    pv[ax] += sgn
                    pv = tuple(pv)
                    if v_ not in pot:
                        pot[v_] = pv
                        stack.append(v_)
                    elif pot[v_] != pv:
                        diff = [pv[k] - pot[v_][k] for k in range(d)]
                        c = lut[u] - 1
                        for k in range(d):
                            if diff[k] != 0:
                                wraps[c, k] = True
    spans = np.zeros(n, bool)
    cand = np.nonzero(sizes >= N)[0] + 1
    if cand.size:
        lut2 = np.zeros(n + 1, np.int32)
        lut2[cand] = np.arange(1, cand.size + 1, dtype=np.int32)
        cl = lut2[labels]
        sel = cl > 0
        for ax in range(d):
            shp = [1] * d
            shp[ax] = N
            coord = np.broadcast_to(np.arange(N).reshape(shp), labels.shape)
            P = np.zeros((cand.size + 1, N), bool)
            P[cl[sel], coord[sel]] = True
            spans[cand - 1] |= P[1:].all(axis=1)
    return labels, n, sizes, wraps, spans


# ------------------------------------------------------------------------------------------------------------
# simulation state
# ------------------------------------------------------------------------------------------------------------
class Sim:
    def __init__(self, d, N, dx, tau_i, block=16):
        self.d, self.N, self.dx, self.tau_i = d, N, dx, tau_i
        self.B = block if d == 2 else max(1, block // 8)
        shp = (N + 2,) + (N,) * (d - 1)
        self.Sp = np.zeros(shp); self.Rp = np.zeros(shp); self.Ip = np.zeros(shp)
        self.PS = np.zeros((N,) * d); self.PR = np.zeros((N,) * d); self.PI = np.zeros((N,) * d)
        B = self.B
        bs = (B,) + (N,) * (d - 1)
        self._t = [np.empty(bs) for _ in range(6)]
        self.sum_s2 = 0.0

    @property
    def S(self):
        return self.Sp[1:-1]

    @property
    def R(self):
        return self.Rp[1:-1]

    @property
    def I(self):
        return self.Ip[1:-1]

    def ghosts(self):
        N = self.N
        for fp in (self.Sp, self.Rp, self.Ip):
            fp[0] = fp[N]
            fp[N + 1] = fp[1]

    def _lap(self, fp, r0, r1, out, tmp, c):
        """out = c * lapraw (without 1/dx^2 factor applied by caller via c)."""
        d = self.d
        cc = fp[r0 + 1:r1 + 1]
        np.add(fp[r0:r1], fp[r0 + 2:r1 + 2], out=out)
        for ax in range(1, d):
            def s(a, b):
                t = [slice(None)] * d
                t[ax] = slice(a, b)
                return tuple(t)
            out[s(1, None)] += cc[s(None, -1)]
            out[s(0, 1)] += cc[s(-1, None)]
            out[s(None, -1)] += cc[s(1, None)]
            out[s(-1, None)] += cc[s(0, 1)]
        np.multiply(cc, 2.0 * d, out=tmp)
        out -= tmp
        out *= c
        return out

    def kick(self, c, a, with_X=True):
        """P += c * F(fields, a).  Also accumulates sum(chi_S^2) for the energy budget."""
        N, B = self.N, self.B
        lS, lR, lI, s2, y2, t = self._t
        cdx = c / self.dx ** 2
        a2 = a * a
        tot = 0.0
        for r0 in range(0, N, B):
            r1 = min(N, r0 + B)
            nb = r1 - r0
            lS_, lR_, lI_, s2_, y2_, t_ = (x[:nb] for x in (lS, lR, lI, s2, y2, t))
            cS = self.Sp[r0 + 1:r1 + 1]
            self._lap(self.Sp, r0, r1, lS_, t_, cdx)
            np.multiply(cS, cS, out=s2_)
            tot += float(s2_.sum())
            if with_X:
                cR = self.Rp[r0 + 1:r1 + 1]
                cI = self.Ip[r0 + 1:r1 + 1]
                self._lap(self.Rp, r0, r1, lR_, t_, cdx)
                self._lap(self.Ip, r0, r1, lI_, t_, cdx)
                np.multiply(cR, cR, out=y2_)
                np.multiply(cI, cI, out=t_)
                y2_ += t_
                np.multiply(y2_, 2.0 * A_G, out=t_)
                t_ += s2_
            else:
                np.copyto(t_, s2_)
            t_ -= a2
            t_ *= cS
            t_ *= c
            lS_ -= t_
            self.PS[r0:r1] += lS_
            if with_X:
                # m2 = A chi_S^2 + 2 B |chi_Y|^2
                y2_ *= 2.0 * B_X
                s2_ *= A_G
                y2_ += s2_
                y2_ *= c
                np.multiply(cR, y2_, out=t_)
                lR_ -= t_
                self.PR[r0:r1] += lR_
                np.multiply(cI, y2_, out=t_)
                lI_ -= t_
                self.PI[r0:r1] += lI_
        self.sum_s2 = tot
        return tot

    def drift(self, dt, with_X=True):
        N, B = self.N, self.B
        t = self._t[0]
        pairs = ((self.Sp, self.PS), (self.Rp, self.PR), (self.Ip, self.PI)) if with_X else ((self.Sp, self.PS),)
        for fp, P in pairs:
            for r0 in range(0, N, B):
                r1 = min(N, r0 + B)
                nb = r1 - r0
                np.multiply(P[r0:r1], dt, out=t[:nb])
                fp[r0 + 1:r1 + 1] += t[:nb]
            fp[0] = fp[N]
            fp[N + 1] = fp[1]

    # --------------------------------------------------------------------------------------------------------
    def grad2_sum(self, f):
        tot = 0.0
        for ax in range(self.d):
            df = np.roll(f, -1, axis=ax) - f
            tot += float(np.einsum('...,...->', df, df) if False else (df * df).sum())
        return tot / self.dx ** 2

    def energy(self, a):
        S, R, I = self.S, self.R, self.I
        dv = self.dx ** self.d
        kin = 0.5 * float((self.PS ** 2).sum()) + float((self.PR ** 2).sum()) + float((self.PI ** 2).sum())
        grad = 0.5 * self.grad2_sum(S) + self.grad2_sum(R) + self.grad2_sum(I)
        s2 = S * S
        y2 = R * R + I * I
        pot = float((0.25 * (s2 - a * a) ** 2 + A_G * s2 * y2 + B_X * y2 * y2).sum())
        return dv * (kin + grad + pot), dv * kin, dv * grad, dv * pot

    def charge(self):
        return 2.0 * float((self.R * self.PI - self.I * self.PR).sum()) * self.dx ** self.d

    def dHdtau(self, a, sum_s2, adot=None):
        nc = self.N ** self.d
        if adot is None:
            adot = 1.0 / self.tau_i
        return -(a * adot) * (sum_s2 - nc * a * a) * self.dx ** self.d


# ------------------------------------------------------------------------------------------------------------
# diagnostics
# ------------------------------------------------------------------------------------------------------------
def link_weights(f, dx, ax, j):
    """PRS / Garagounis-Hindmarsh weight |grad f| / sum_i |f_,i| on the links (j, j + e_ax)."""
    d = f.ndim
    N = f.shape[0]
    jp = list(j)
    jp[ax] = (j[ax] + 1) % N
    jp = tuple(jp)
    g = (f[jp] - f[j]) / dx
    gsq = g * g
    gabs = np.abs(g)
    for b in range(d):
        if b == ax:
            continue
        def sh(base, k):
            t = list(base)
            t[b] = (base[b] + k) % N
            return tuple(t)
        c0 = (f[sh(j, 1)] - f[sh(j, -1)]) / (2 * dx)
        c1 = (f[sh(jp, 1)] - f[sh(jp, -1)]) / (2 * dx)
        gb = 0.5 * (c0 + c1)
        gsq = gsq + gb * gb
        gabs = gabs + np.abs(gb)
    return np.sqrt(gsq) / np.maximum(gabs, 1e-300)


def measure(sim, a, tau, qlist=False):
    """All diagnostics at a synchronous state.  Physical fields: S = chi_S/a, F^2 = 2|chi_Y|^2/a^2,
    physical charge density n = 2 (chi_r pi_i - chi_i pi_r)/a^3 (units sqrt(lam) v^3)."""
    d, N, dx = sim.d, sim.N, sim.dx
    V = (N * dx) ** d
    dv = dx ** d
    S = sim.S / a
    R, I = sim.R, sim.I
    F2 = 2.0 * (R * R + I * I) / (a * a)
    qc = 2.0 * (R * sim.PI - I * sim.PR)          # comoving charge density per comoving volume
    Qtot = float(qc.sum()) * dv
    Qabs = float(np.abs(qc).sum()) * dv
    plus = S > 0.5
    minus = S < -0.5
    core = ~(plus | minus)
    out = {"tau": tau, "a": a}
    nc = S.size
    out["f_plus"] = float(plus.sum()) / nc
    out["f_minus"] = float(minus.sum()) / nc
    out["f_core"] = float(core.sum()) / nc
    out["mean_S"] = float(S.mean())
    out["mean_abs_S"] = float(np.abs(S).mean())
    out["f_F2gt1"] = float((F2 > 1).sum()) / nc
    out["mean_n_phys"] = Qtot / V / a ** 3
    Qn = Qtot if abs(Qtot) > 0 else 1.0
    out["qfrac_core"] = float(qc[core].sum()) * dv / Qn
    out["qfrac_F2gt1"] = float(qc[F2 > 1].sum()) * dv / Qn
    out["qfrac_vac"] = float(qc[~core].sum()) * dv / Qn

    # --- core components and their vacuum contacts
    lc, ncomp, csz, cwrap, cspan = plabel(core, want_wrap=False)
    touchP = np.zeros(ncomp + 1, bool)
    touchM = np.zeros(ncomp + 1, bool)
    for ax in range(d):
        for sh in (1, -1):
            nbp = np.roll(plus, sh, axis=ax)
            nbm = np.roll(minus, sh, axis=ax)
            touchP |= np.bincount(lc[core & nbp], minlength=ncomp + 1) > 0
            touchM |= np.bincount(lc[core & nbm], minlength=ncomp + 1) > 0
    ctype = np.zeros(ncomp + 1, np.int8)   # 0 none (restored), 1 + only, 2 - only, 3 both (wall-type)
    ctype[touchP & ~touchM] = 1
    ctype[touchM & ~touchP] = 2
    ctype[touchP & touchM] = 3
    ctype[0] = -1
    cell_type = ctype[lc]                   # -1 for vacuum cells
    wallcomp = cell_type == 3
    qcomp = (cell_type == 1) | (cell_type == 2)
    restored = cell_type == 0
    out["f_core_walltype"] = float(wallcomp.sum()) / nc
    out["f_core_single"] = float(qcomp.sum()) / nc
    out["f_core_restored"] = float(restored.sum()) / nc
    out["qfrac_core_walltype"] = float(qc[wallcomp].sum()) * dv / Qn
    out["qfrac_core_walltype_F2lt1"] = float(qc[wallcomp & (F2 < 1)].sum()) * dv / Qn
    out["qfrac_core_single"] = float(qc[qcomp].sum()) * dv / Qn
    out["qfrac_core_restored"] = float(qc[restored].sum()) * dv / Qn
    out["n_core_comp"] = int(ncomp)
    out["n_core_walltype"] = int((ctype[1:] == 3).sum())
    out["n_core_single"] = int(((ctype[1:] == 1) | (ctype[1:] == 2)).sum())

    # --- sign-change links: walls vs internal (Q-ball) zero crossings
    pos = S > 0
    Lall = Lwall = Lbare = Lint = 0.0
    nlinks = nwl = 0
    for ax in range(d):
        sc = pos != np.roll(pos, -1, axis=ax)
        j = np.nonzero(sc)
        if j[0].size == 0:
            continue
        w = link_weights(S, dx, ax, j)
        jp = list(j)
        jp[ax] = (j[ax] + 1) % N
        jp = tuple(jp)
        t0, t1 = cell_type[j], cell_type[jp]
        vac_vac = (t0 == -1) & (t1 == -1)
        iswall = (t0 == 3) | (t1 == 3) | vac_vac
        bare = iswall & (F2[j] < 1) & (F2[jp] < 1)
        Lall += float(w.sum())
        Lwall += float(w[iswall].sum())
        Lbare += float(w[bare].sum())
        Lint += float(w[~iswall].sum())
        nlinks += j[0].size
        nwl += int(iswall.sum())
    sl = dx ** (d - 1)
    out["L_all"] = Lall * sl / V        # comoving length (area) per comoving area (volume)
    out["L_wall"] = Lwall * sl / V
    out["L_bare"] = Lbare * sl / V
    out["L_internal"] = Lint * sl / V
    out["A_wall"] = out["L_wall"] * tau / 2.0     # network_sim area parameter A = (A_c/V_c) tau / 2
    out["A_all"] = out["L_all"] * tau / 2.0
    out["xi_wall_comoving"] = (1.0 / out["L_wall"]) if out["L_wall"] > 0 else None
    out["xi_wall_physical"] = (a / out["L_wall"]) if out["L_wall"] > 0 else None

    # --- vacuum percolation
    for name, m in (("plus", plus), ("minus", minus)):
        lab, n, sz, wr, sp = plabel(m)
        out[f"n_comp_{name}"] = int(n)
        out[f"largest_frac_{name}"] = float(sz.max()) / nc if n else 0.0
        out[f"wraps_{name}"] = bool(wr.any()) if n else False
        out[f"wraps_both_axes_{name}"] = bool(wr.all(axis=1).any()) if n else False
        out[f"spans_{name}"] = bool(sp.any()) if n else False
        out[f"wrapping_volume_frac_{name}"] = float(sz[wr.any(axis=1)].sum()) / nc if n else 0.0

    # --- charge-restored lumps (F^2 > 1): Q-balls / Q-matter
    lq, nq, qsz, qwrap, qspan = plabel(F2 > 1)
    if nq:
        zf = a if d == 2 else 1.0   # 2D: charge per physical ell of z-length; 3D: total charge
        qlump = np.bincount(lq.ravel(), weights=qc.ravel(), minlength=nq + 1)[1:] * dv / zf
        att = np.bincount(lq[wallcomp], minlength=nq + 1)[1:] > 0
        corefrac = np.bincount(lq[core], minlength=nq + 1)[1:] / np.maximum(qsz, 1)
        big = np.abs(qlump) > 1.0
        out["n_lumps"] = int(big.sum())
        out["n_lumps_attached"] = int((big & att).sum())
        out["n_lumps_wrapping"] = int(qwrap.any(axis=1).sum())
        Qphys_tot = Qtot / zf
        out["qfrac_lumps"] = float(qlump[big].sum()) / (Qphys_tot if Qphys_tot else 1.0)
        out["qfrac_lumps_attached"] = float(qlump[big & att].sum()) / (Qphys_tot if Qphys_tot else 1.0)
        if big.any():
            qb = np.sort(qlump[big])
            out["lump_Q_quantiles"] = [float(np.quantile(qb, p)) for p in (0.1, 0.5, 0.9)] + [float(qb.max())]
            out["lump_mean_core_frac"] = float(np.average(corefrac[big], weights=np.abs(qlump[big])))
            # mean physical radius from area (2D: pi R^2) or volume (3D)
            area_phys = qsz[big] * dv * a ** d
            if d == 2:
                rr = np.sqrt(area_phys / np.pi)
            else:
                rr = (3 * area_phys / (4 * np.pi)) ** (1 / 3)
            out["lump_R_phys_median"] = float(np.median(rr))
            if qlist:
                out["lump_Q_list"] = [round(float(x), 3) for x in qlump[big]]
                out["lump_attached_list"] = [bool(x) for x in att[big]]
    else:
        out["n_lumps"] = 0
        out["n_lumps_attached"] = 0
        out["qfrac_lumps"] = 0.0
        out["qfrac_lumps_attached"] = 0.0
    out["Q_comoving"] = Qtot
    out["Q_abs_comoving"] = Qabs
    return out


# ------------------------------------------------------------------------------------------------------------
# run driver
# ------------------------------------------------------------------------------------------------------------
def setup(cfg):
    d, N, dx, tau_i = cfg["d"], cfg["N"], cfg["dx"], cfg["tau_i"]
    L = N * dx
    sim = Sim(d, N, dx, tau_i, block=cfg.get("block", 16))
    kind = cfg["kind"]
    seed = cfg["seed"]
    Nb = cfg.get("N_base", N)
    s_c = cfg.get("s_c", 1.0)
    noise_S = gaussian_field(d, N, L, s_c, cfg.get("eps_S", 1e-2), seed, Nb)
    info = {}
    if kind == "charged_exit":           # (A) homogeneous phi = 0 branch-2 condensate at n_i = kappa n_c
        nbar = cfg["kappa"] * NBAR_C
        S0, Y0, om = branch2_state(nbar)
        radial = 0.0                       # conformal: chi_Y circular orbit is exact
        sdot = 0.0
    elif kind == "neutral":              # (B)
        S0, Y0, om, radial, sdot, nbar = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    elif kind == "memory_sub":           # (C1) phi = +v side, branch 1, n below n_c
        nbar = cfg["kappa"] * NBAR_C
        S0, Y0, om = branch1_state(nbar)
        radial = -0.5 / tau_i
        sdot = S0 / tau_i
    elif kind == "memory_forced":        # (C2) phi = +v, X orbit in the fixed +v vacuum, n above n_s1
        nbar = cfg["kappa"] * NBAR_C
        S0, Y0, om = fixed_vacuum_state(nbar, 1.0)
        radial = -0.5 / tau_i
        sdot = S0 / tau_i
    else:
        raise ValueError(kind)
    info.update(S0=S0, Y0=Y0, omega=om, nbar_i=nbar, F2_i=2 * Y0 ** 2)
    sim.S[...] = S0 + noise_S
    sim.PS[...] = sdot
    if Y0 > 0:
        eX = cfg.get("eps_X", 1e-3)
        nr = gaussian_field(d, N, L, s_c, eX, seed + 1000, Nb)
        ni = gaussian_field(d, N, L, s_c, eX, seed + 2000, Nb)
        sim.R[...] = Y0 + nr
        sim.I[...] = ni
        # pY = (i om + radial) chi_Y
        sim.PR[...] = radial * sim.R - om * sim.I
        sim.PI[...] = radial * sim.I + om * sim.R
    sim.ghosts()
    return sim, info


def run(cfg, outpath, ckpt_dir=None, max_wall=None, verbose=True, snap_dir=None):
    t_start = time.time()
    d, N, dx, tau_i = cfg["d"], cfg["N"], cfg["dx"], cfg["tau_i"]
    dt = cfg.get("dt", dx / 5.0)
    a_f = cfg["a_final"]
    tau_f = a_f * tau_i
    nsteps = int(round((tau_f - tau_i) / dt))
    out_every = max(1, int(round(cfg.get("dtau_out", 1.0) / dt)))
    qlist_at = set(cfg.get("qlist_at_a", []))
    with_X = cfg["kind"] != "neutral" and not cfg.get("zero_X", False)
    a_fr = cfg.get("a_freeze")          # optional: stop the expansion at a = a_freeze (Minkowski continuation)

    def a_of(t):
        a_ = t / tau_i
        return min(a_, a_fr) if a_fr else a_

    def adot_of(t):
        return 0.0 if (a_fr and t / tau_i >= a_fr) else 1.0 / tau_i

    if a_fr:
        tau_f = cfg["tau_final"]
        nsteps = int(round((tau_f - tau_i) / dt))
    ck = os.path.join(ckpt_dir, cfg["name"] + ".npz") if ckpt_dir else None
    if ck and os.path.exists(ck):
        z = np.load(ck, allow_pickle=True)
        sim = Sim(d, N, dx, tau_i, block=cfg.get("block", 16))
        sim.Sp[...] = z["Sp"]; sim.Rp[...] = z["Rp"]; sim.Ip[...] = z["Ip"]
        sim.PS[...] = z["PS"]; sim.PR[...] = z["PR"]; sim.PI[...] = z["PI"]
        state = json.loads(str(z["state"]))
        step0 = state["step"]; info = state["info"]; series = state["series"]; budget = state["budget"]
        prev_dH = state["prev_dH"]; H0 = state["H0"]; Q0 = state["Q0"]
        wall_used = state.get("wall_used", 0.0)
        if verbose:
            print(f"resumed {cfg['name']} at step {step0}", flush=True)
        # PS etc. are synchronous at step0; re-open with a half kick below
    else:
        if cfg.get("init_from"):
            # continue from a saved synchronous state of another run (optionally with X and its momentum removed:
            # 'decharged' control with the identical phi configuration)
            z = np.load(cfg["init_from"], allow_pickle=True)
            sim = Sim(d, N, dx, tau_i, block=cfg.get("block", 16))
            sim.Sp[...] = z["Sp"]; sim.Rp[...] = z["Rp"]; sim.Ip[...] = z["Ip"]
            sim.PS[...] = z["PS"]; sim.PR[...] = z["PR"]; sim.PI[...] = z["PI"]
            info = json.loads(str(z["info"]))
            tau0 = float(z["tau"])
            step0 = int(round((tau0 - tau_i) / dt))
            info["init_from"] = os.path.basename(cfg["init_from"])
            info["tau_start"] = tau0
            if cfg.get("zero_X"):
                for arr_ in (sim.Rp, sim.Ip, sim.PR, sim.PI):
                    arr_[...] = 0.0
                info["zero_X"] = True
        else:
            sim, info = setup(cfg)
            step0 = 0
        series = []
        tau0 = tau_i + step0 * dt
        a0 = a_of(tau0)
        H0 = sim.energy(a0)[0]
        Q0 = sim.charge()
        budget = 0.0
        sum_s2 = float((sim.S * sim.S).sum())
        prev_dH = sim.dHdtau(a0, sum_s2, adot_of(tau0))
        wall_used = 0.0
        m = measure(sim, a0, tau0)
        E = sim.energy(a0)
        m.update(H=E[0], H_kin=E[1], H_grad=E[2], H_pot=E[3], budget_resid=0.0, dQ_rel=0.0, step=0,
                 wall_s=time.time() - t_start)
        series.append(m)
    # opening half kick at step0
    tau = tau_i + step0 * dt
    sim.kick(0.5 * dt, a_of(tau), with_X)
    step = step0
    last_budget_tau = tau
    while step < nsteps:
        sim.drift(dt, with_X)
        step += 1
        tau = tau_i + step * dt
        a = a_of(tau)
        is_out = (step % out_every == 0) or step == nsteps
        if is_out:
            ss = sim.kick(0.5 * dt, a, with_X)
        else:
            ss = sim.kick(dt, a, with_X)
        dH = sim.dHdtau(a, ss, adot_of(tau))
        budget += 0.5 * dt * (prev_dH + dH)
        prev_dH = dH
        if is_out:
            E = sim.energy(a)
            Q = sim.charge()
            want_q = any(abs(a - x) < 0.5 * dt * out_every / tau_i for x in qlist_at) or step == nsteps
            m = measure(sim, a, tau, qlist=want_q)
            if snap_dir and any(abs(a - x) < 0.5 * dt / tau_i * out_every for x in cfg.get("snap_at_a", [])):
                np.savez_compressed(os.path.join(snap_dir, f"{cfg['name']}_a{a:.3f}.npz"),
                                    S=(sim.S / a).astype(np.float32),
                                    F2=(2 * (sim.R ** 2 + sim.I ** 2) / a ** 2).astype(np.float32),
                                    n=(2 * (sim.R * sim.PI - sim.I * sim.PR) / a ** 3).astype(np.float32))
            m.update(H=E[0], H_kin=E[1], H_grad=E[2], H_pot=E[3], budget_resid=E[0] - H0 - budget,
                     budget_flow=budget, dQ_rel=(Q - Q0) / (abs(Q0) if Q0 else 1.0), Q_abs_drift=Q - Q0,
                     step=step, wall_s=wall_used + time.time() - t_start)
            series.append(m)
            if verbose:
                print(f"{cfg['name']} tau={tau:.2f} a={a:.3f} f+={m['f_plus']:.3f} f-={m['f_minus']:.3f} "
                      f"fc={m['f_core']:.3f} Lw={m['L_wall']:.4f} A={m['A_wall']:.3f} qcore={m['qfrac_core']:.3f} "
                      f"qwall={m['qfrac_core_walltype']:.3f} nl={m['n_lumps']} dQ={m['dQ_rel']:.2e} "
                      f"res={m['budget_resid'] / max(abs(E[1]) + abs(E[2]), 1e-30):.2e} t={m['wall_s']:.0f}s",
                      flush=True)
            if step < nsteps:
                sim.kick(0.5 * dt, a, with_X)
            # checkpoint if we are running out of wall time
            if ck and max_wall and (time.time() - t_start) > max_wall and step < nsteps:
                # state is synchronous after the measurement only if we undo the last half kick
                sim.kick(-0.5 * dt, a, with_X)
                state = dict(step=step, info=info, series=series, budget=budget, prev_dH=prev_dH, H0=H0, Q0=Q0,
                             wall_used=wall_used + time.time() - t_start)
                np.savez(ck, Sp=sim.Sp, Rp=sim.Rp, Ip=sim.Ip, PS=sim.PS, PR=sim.PR, PI=sim.PI,
                         state=json.dumps(state))
                if verbose:
                    print(f"checkpointed {cfg['name']} at step {step}", flush=True)
                return "incomplete"
    if cfg.get("save_state"):
        np.savez(cfg["save_state"], Sp=sim.Sp, Rp=sim.Rp, Ip=sim.Ip, PS=sim.PS, PR=sim.PR, PI=sim.PI,
                 tau=tau, info=json.dumps(info))
    res = dict(config=cfg, info=info, dt=dt, nsteps=nsteps, H0=H0, Q0=Q0, series=series,
               wall_s=wall_used + time.time() - t_start,
               derived=dict(n_c_bar=NBAR_C, n_s1_bar=NBAR_S1, n_star_bar=NBAR_STAR,
                            a_cross_nc=(cfg.get("kappa", 0) ** (1 / 3) if cfg["kind"] == "charged_exit" else None)))
    with open(outpath, "w") as f:
        json.dump(res, f, separators=(",", ":"))
    if ck and os.path.exists(ck):
        os.remove(ck)
    return "done"
