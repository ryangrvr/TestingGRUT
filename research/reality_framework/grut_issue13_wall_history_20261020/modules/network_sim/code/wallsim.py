"""Classical lattice evolution of Z2 domain walls of the unchanged Candidate 1 phi potential
(sector X = D = 0) in a radiation-dominated FRW background, plus wall-population diagnostics.

Units: v = 1 (field), m_phi = sqrt(2 lambda_phi) v = 1 (inverse length).  Then
    U(phi) = (lam/4) (phi^2 - 1)^2 with lam = 1/2,
    planar wall phi = tanh(z / delta), delta = 2 / m_phi = 2,  sigma = 2/3 (in v^2 m_phi units).
Physical benchmark conversion: m_phi = 141.42 GeV, v = 1000 GeV, sigma = 9.428e7 GeV^3 (not used in the dynamics,
which only depend on the ratios fixed above).

Equation (Press-Ryden-Spergel family, conformal time tau, comoving x, a = tau / tau_a in radiation era):
    phi'' + alpha (dln a/dln tau) (1/tau) phi' - lap(phi) + a^beta dU/dphi = 0.
alpha = 2, beta = 2 is the physical 3+1 equation (2D runs = fields independent of z).
beta = 0 (constant comoving wall width) with alpha = D is the PRS prescription.

Discretization: second-order (2d+1)-point Laplacian, leapfrog with PRS/standard friction discretization
    pi_{n+1/2} = [(1 - delta_n) pi_{n-1/2} + dt F_n] / (1 + delta_n),  delta_n = alpha dt / (2 tau_n)
    phi_{n+1}  = phi_n + dt pi_{n+1/2}.
"""
import json
import math
import os
import time

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

LAM = 0.5          # lambda_phi in units v = m_phi = 1
DELTA = 2.0        # tanh half-width of the planar wall in m_phi^-1
SIGMA = 2.0 / 3.0  # tension in units v^2 m_phi


def dU(phi, out=None, tmp=None):
    """dU/dphi = lam phi (phi^2 - 1)."""
    if out is None:
        return LAM * phi * (phi * phi - 1.0)
    np.multiply(phi, phi, out=out)
    out -= 1.0
    out *= phi
    out *= LAM
    return out


def U(phi):
    return 0.25 * LAM * (phi * phi - 1.0) ** 2


def lap_inplace(phi, out):
    """Periodic (2d+1)-point Laplacian times dx^2, written into out."""
    d = phi.ndim
    np.multiply(phi, -2.0 * d, out=out)
    for ax in range(d):
        def s(a, b):
            t = [slice(None)] * d
            t[ax] = slice(a, b)
            return tuple(t)
        out[s(1, None)] += phi[s(None, -1)]
        out[s(0, 1)] += phi[s(-1, None)]
        out[s(None, -1)] += phi[s(1, None)]
        out[s(-1, None)] += phi[s(0, 1)]
    return out


# ----------------------------------------------------------------------------------------------
# initial conditions
# ----------------------------------------------------------------------------------------------
def gaussian_ic(d, N, L, s_c, amp, seed, N_base=None, dtype=np.float32):
    """Zero-mean Gaussian random field: white noise on a base grid N_base (same box L) smoothed by a
    Gaussian kernel of standard deviation s_c (comoving), band-limited and zero-padded to N.
    Correlation function C(r) = amp^2 exp(-r^2 / (4 s_c^2)) (up to lattice effects).
    Using the same (seed, N_base, L, s_c) on different N gives the same continuum field (resolution checks)."""
    if N_base is None:
        N_base = N
    assert N % N_base == 0
    rng = np.random.default_rng(seed)
    white = rng.standard_normal((N_base,) * d)
    fk = np.fft.rfftn(white)
    dxb = L / N_base
    ks = [2 * np.pi * np.fft.fftfreq(N_base, d=dxb)] * (d - 1) + [2 * np.pi * np.fft.rfftfreq(N_base, d=dxb)]
    k2 = np.zeros(fk.shape)
    for ax, k in enumerate(ks):
        shp = [1] * d
        shp[ax] = k.size
        k2 = k2 + (k.reshape(shp)) ** 2
    fk *= np.exp(-0.5 * k2 * s_c ** 2)
    fk.flat[0] = 0.0
    # remove Nyquist planes (negligible after smoothing; keeps the embedding exactly real)
    h = N_base // 2
    for ax in range(d):
        sl = [slice(None)] * d
        sl[ax] = h
        fk[tuple(sl)] = 0.0
    if N == N_base:
        f = np.fft.irfftn(fk, s=(N,) * d)
    else:
        big = np.zeros((N,) * (d - 1) + (N // 2 + 1,), dtype=complex)
        # embed: for full axes take [0:h] and [-h:], last axis [0:h+1]
        idx_full = np.r_[0:h, N - h:N]
        idx_src = np.r_[0:h, N_base - h:N_base]
        if d == 2:
            big[np.ix_(idx_full, np.arange(h + 1))] = fk[np.ix_(idx_src, np.arange(h + 1))]
        else:
            big[np.ix_(idx_full, idx_full, np.arange(h + 1))] = fk[np.ix_(idx_src, idx_src, np.arange(h + 1))]
        f = np.fft.irfftn(big, s=(N,) * d) * (N / N_base) ** d
    # normalise with the base-grid standard deviation so that all resolutions share the same field
    fb = np.fft.irfftn(fk, s=(N_base,) * d)
    f *= amp / fb.std()
    return f.astype(dtype)


# ----------------------------------------------------------------------------------------------
# periodic connected-component labelling
# ----------------------------------------------------------------------------------------------
def plabel(mask, structure):
    """Connected components of a boolean mask on a periodic lattice.
    Returns (labels, n) with labels in 0..n (0 = background), compact."""
    d = mask.ndim
    N = mask.shape[0]
    m = np.pad(mask, 1, mode="wrap")
    lab, n = ndimage.label(m, structure=structure)
    if n == 0:
        return np.zeros(mask.shape, dtype=np.int32), 0
    pa, pb = [], []
    for ax in range(d):
        for p_idx, q_idx in ((0, N), (N + 1, 1)):
            a = np.take(lab, p_idx, axis=ax).ravel()
            b = np.take(lab, q_idx, axis=ax).ravel()
            sel = a > 0
            pa.append(a[sel])
            pb.append(b[sel])
    pa = np.concatenate(pa)
    pb = np.concatenate(pb)
    g = coo_matrix((np.ones(pa.size, dtype=np.int8), (pa, pb)), shape=(n + 1, n + 1))
    ncomp, comp = connected_components(g, directed=False)
    inner = tuple(slice(1, -1) for _ in range(d))
    li = lab[inner]
    # compact relabel of components that actually occur among labels 1..n
    used = np.unique(comp[1:])
    remap = np.zeros(ncomp, dtype=np.int32)
    remap[used] = np.arange(1, used.size + 1, dtype=np.int32)
    lut = remap[comp]
    lut[0] = 0
    return lut[li], int(used.size)


def spans(labels, ell):
    """True if component ell has cells at every coordinate along some axis (box-spanning)."""
    idx = np.nonzero(labels == ell)
    N = labels.shape[0]
    for ax in range(labels.ndim):
        if np.unique(idx[ax]).size == N:
            return True
    return False


def spanning_flags(labels, n, sizes):
    """Vectorised box-spanning test for all labels 1..n: a component spans if its projection onto some axis
    covers all N coordinates. Only components with >= N cells can span. Returns bool array of length n."""
    N = labels.shape[0]
    d = labels.ndim
    flags = np.zeros(n, dtype=bool)
    cand = np.nonzero(sizes >= N)[0] + 1
    if cand.size == 0:
        return flags
    lut = np.zeros(n + 1, dtype=np.int32)
    lut[cand] = np.arange(1, cand.size + 1, dtype=np.int32)
    cl = lut[labels]
    for ax in range(d):
        P = np.zeros((cand.size + 1, N), dtype=bool)
        shp = [1] * d
        shp[ax] = N
        coord = np.broadcast_to(np.arange(N).reshape(shp), labels.shape)
        sel = cl > 0
        P[cl[sel], coord[sel]] = True
        flags[cand - 1] |= P[1:].all(axis=1)
    return flags


# ----------------------------------------------------------------------------------------------
# wall diagnostics
# ----------------------------------------------------------------------------------------------
def wall_links(phi):
    """For each axis, the indices j with sign change between j and j+e_ax (periodic)."""
    pos = phi > 0
    out = []
    for ax in range(phi.ndim):
        s = pos != np.roll(pos, -1, axis=ax)
        out.append(np.nonzero(s))
    return pos, out


def prs_area(phi, dx, links):
    """Comoving wall area (length in 2D) by the PRS weighted link count:
    A = dx^(d-1) sum_links |grad phi| / sum_i |phi_,i| ; link-direction derivative is the link difference,
    transverse derivatives are the average of central differences at both ends (Garagounis-Hindmarsh style).
    Also returns the isotropic-corrected raw count (factor pi/4 in 2D, 2/3 in 3D) and per-link weights."""
    d = phi.ndim
    N = phi.shape[0]
    total_w = 0.0
    total_n = 0
    weights = []
    for ax in range(d):
        idx = links[ax]
        n = idx[0].size
        total_n += n
        if n == 0:
            weights.append(np.zeros(0))
            continue
        j = [i.astype(np.int64) for i in idx]
        jp = list(j)
        jp[ax] = (j[ax] + 1) % N
        g_ax = (phi[tuple(jp)].astype(np.float64) - phi[tuple(j)]) / dx
        gsq = g_ax ** 2
        gabs = np.abs(g_ax)
        for b in range(d):
            if b == ax:
                continue
            def sh(base, k):
                t = list(base)
                t[b] = (base[b] + k) % N
                return tuple(t)
            c0 = (phi[sh(j, 1)].astype(np.float64) - phi[sh(j, -1)]) / (2 * dx)
            c1 = (phi[sh(jp, 1)].astype(np.float64) - phi[sh(jp, -1)]) / (2 * dx)
            g_b = 0.5 * (c0 + c1)
            gsq += g_b ** 2
            gabs += np.abs(g_b)
        w = np.sqrt(gsq) / gabs
        weights.append(w)
        total_w += w.sum()
    iso = math.pi / 4 if d == 2 else 2.0 / 3.0
    return total_w * dx ** (d - 1), total_n * iso * dx ** (d - 1), weights


def measure(phi, pi, dx, tau, a, beta, clusters=True):
    d = phi.ndim
    N = phi.shape[0]
    V = (N * dx) ** d
    pos, links = wall_links(phi)
    area_prs, area_iso, weights = prs_area(phi, dx, links)
    nlinks = int(sum(l[0].size for l in links))
    res = {
        "tau": float(tau),
        "a": float(a),
        "nlinks": nlinks,
        "area_prs": float(area_prs),
        "area_iso": float(area_iso),
        "area_density": float(area_prs / V),
        # A = rho_w t / sigma = (A_c/V_c) / a * t, with t = a tau / 2 in radiation era
        "A": float(area_prs / V * tau / 2.0),
        "A_iso": float(area_iso / V * tau / 2.0),
        "xi_over_tau": float(V / area_prs / tau) if area_prs > 0 else None,
        "f_plus": float(pos.mean()),
        "phi_min": float(phi.min()),
        "phi_max": float(phi.max()),
        "phi_mean": float(phi.mean(dtype=np.float64)),
        "frac_absphi_lt_half": float((np.abs(phi) < 0.5).mean()),
    }
    # energies (comoving densities, per comoving volume) and wall velocity estimator
    kin = 0.5 * float(np.mean(pi.astype(np.float64) ** 2))
    grad2 = np.zeros(phi.shape, dtype=np.float64)
    gradc2 = np.zeros(phi.shape, dtype=np.float32)
    for ax in range(d):
        fw = (np.roll(phi, -1, axis=ax) - phi) / dx
        grad2 += fw.astype(np.float64) ** 2
        c = (np.roll(phi, -1, axis=ax) - np.roll(phi, 1, axis=ax)) / (2 * dx)
        gradc2 += c * c
    res["E_kin"] = kin
    res["E_grad"] = 0.5 * float(grad2.mean())
    res["E_pot"] = float(a ** beta * np.mean(U(phi.astype(np.float64))))
    core = np.abs(phi) < 0.5
    if core.any():
        num = float(np.sum(pi[core].astype(np.float64) ** 2))
        den = float(np.sum(gradc2[core].astype(np.float64)))
        res["v2_wall"] = num / den if den > 0 else None
    else:
        res["v2_wall"] = None
    if not clusters:
        return res
    face = ndimage.generate_binary_structure(d, 1)
    full = ndimage.generate_binary_structure(d, d)
    lp, n_p = plabel(pos, face)
    lm, n_m = plabel(~pos, face)
    res["n_domains_plus"] = n_p
    res["n_domains_minus"] = n_m
    res["n_domains"] = n_p + n_m
    for name, lab, n in (("plus", lp, n_p), ("minus", lm, n_m)):
        if n > 0:
            sizes = np.bincount(lab.ravel(), minlength=n + 1)[1:]
            fl = spanning_flags(lab, n, sizes)
            big = int(np.argmax(sizes))
            res[f"largest_{name}_frac"] = float(sizes.max() / phi.size)
            res[f"largest_{name}_spans"] = bool(fl[big])
            res[f"n_{name}_domains_spanning"] = int(fl.sum())
        else:
            res[f"largest_{name}_frac"] = 0.0
            res[f"largest_{name}_spans"] = False
            res[f"n_{name}_domains_spanning"] = 0
    del lp, lm
    # wall clusters: cells adjacent (face) to an opposite-sign cell, connected with full connectivity
    if nlinks > 0:
        wmask = np.zeros(phi.shape, dtype=bool)
        for ax in range(d):
            s = pos != np.roll(pos, -1, axis=ax)
            wmask |= s
            wmask |= np.roll(s, 1, axis=ax)
        lw, n_w = plabel(wmask, full)
        area_by = np.zeros(n_w + 1)
        for ax in range(d):
            lab_link = lw[links[ax]]
            area_by += np.bincount(lab_link, weights=weights[ax], minlength=n_w + 1)
        area_by = area_by[1:] * dx ** (d - 1)
        cells = np.bincount(lw.ravel(), minlength=n_w + 1)[1:]
        order = np.argsort(area_by)[::-1]
        fl = spanning_flags(lw, n_w, cells)
        span_area = float(area_by[fl].sum())
        n_span = int(fl.sum())
        tot = area_by.sum()
        res["n_wall_clusters"] = int(n_w)
        res["largest_cluster_area_frac"] = float(area_by[order[0]] / tot)
        res["spanning_cluster_area_frac"] = float(span_area / tot)
        res["n_spanning_clusters"] = int(n_span)
        res["closed_wall_area_frac"] = float(1.0 - span_area / tot)
    else:
        res["n_wall_clusters"] = 0
        res["largest_cluster_area_frac"] = None
        res["spanning_cluster_area_frac"] = None
        res["n_spanning_clusters"] = 0
        res["closed_wall_area_frac"] = None
    return res


# ----------------------------------------------------------------------------------------------
# evolution
# ----------------------------------------------------------------------------------------------
def evolve(cfg, outdir=None, verbose=True):
    """cfg keys: name, d, N, L, alpha, beta, tau_a (a = tau/tau_a), tau_i, tau_end, courant (dt = courant dx),
    s_c, amp, seed, N_base, p (positive fraction; None = unbiased zero mean), n_out, dtype, clusters,
    stop_when_empty, energy_check."""
    t_start = time.time()
    d, N, L = cfg["d"], cfg["N"], cfg["L"]
    dx = L / N
    alpha, beta = cfg["alpha"], cfg["beta"]
    tau_a = cfg.get("tau_a", 1.0)
    tau_i, tau_end = cfg["tau_i"], cfg["tau_end"]
    dt = cfg["courant"] * dx
    dtype = np.float64 if cfg.get("dtype", "f32") == "f64" else np.float32
    phi = gaussian_ic(d, N, L, cfg["s_c"], cfg["amp"], cfg["seed"], cfg.get("N_base"), dtype=dtype)
    p = cfg.get("p")
    shift = 0.0
    if p is not None and abs(p - 0.5) > 1e-12:
        shift = -float(np.quantile(phi.astype(np.float64), 1.0 - p))
        phi += dtype(shift)
    else:
        phi -= dtype(phi.mean(dtype=np.float64))
    pi = np.zeros_like(phi)
    buf = np.empty_like(phi)
    buf2 = np.empty_like(phi)
    nsteps = int(math.ceil((tau_end - tau_i) / dt - 1e-9))
    n_out = cfg.get("n_out", 32)
    out_taus = np.unique(np.round((np.geomspace(tau_i, tau_end, n_out) - tau_i) / dt).astype(int))
    out_steps = set(int(s) for s in out_taus) | {0, nsteps}
    energy_check = cfg.get("energy_check", False)
    series = []
    echeck = []
    W = 0.0
    prev_rate = None
    E0 = None
    empty_count = 0
    inv_dx2 = 1.0 / dx ** 2

    def a_of(t):
        return t / tau_a

    def energy(phi_, pi_mid, t):
        g2 = 0.0
        for ax in range(d):
            fw = (np.roll(phi_, -1, axis=ax) - phi_).astype(np.float64) / dx
            g2 += float(np.sum(fw * fw))
        kin = 0.5 * float(np.sum(pi_mid.astype(np.float64) ** 2))
        pot = a_of(t) ** beta * float(np.sum(U(phi_.astype(np.float64))))
        return (kin + 0.5 * g2 + pot) * dx ** d, pot * dx ** d

    pi_prev = pi.copy() if energy_check else None
    stopped_early = False
    for n in range(nsteps + 1):
        tau = tau_i + n * dt
        a = a_of(tau)
        if n in out_steps:
            m = measure(phi, pi, dx, tau, a, beta, clusters=cfg.get("clusters", True))
            m["step"] = n
            if beta == 2:
                m["wall_cells_2delta"] = float(2 * DELTA / (a * dx))
            else:
                m["wall_cells_2delta"] = float(2 * DELTA / dx)
            m["horizon_over_L"] = float(tau / L)
            series.append(m)
            if verbose:
                print(f"[{cfg['name']}] tau={tau:9.3f} A={m['A']:.4f} f+={m['f_plus']:.4f} "
                      f"nd={m.get('n_domains')} span={m.get('spanning_cluster_area_frac')} "
                      f"v2={m['v2_wall']} t={time.time()-t_start:.0f}s", flush=True)
            if cfg.get("stop_when_empty", False):
                if m["nlinks"] == 0:
                    empty_count += 1
                    if empty_count >= 2:
                        stopped_early = True
                        break
                else:
                    empty_count = 0
        if n == nsteps:
            break
        if time.time() - t_start > cfg.get("max_seconds", 24 * 60):
            stopped_early = "timeout"
            if verbose:
                print(f"[{cfg['name']}] TIMEOUT at tau={tau}", flush=True)
            break
        # force
        lap_inplace(phi, buf)
        buf *= inv_dx2
        dU(phi, out=buf2)
        if beta != 0:
            buf2 *= a ** beta
        buf -= buf2
        delta = 0.5 * alpha * dt / tau
        if tau < cfg.get("damp_until", -1.0):
            # optional numerical preparation phase (Garagounis-Hindmarsh style): extra friction damp_gamma * phi'
            delta += 0.5 * cfg["damp_gamma"] * dt
        if energy_check:
            pi_prev[...] = pi
        pi *= (1.0 - delta)
        buf *= dt
        pi += buf
        pi *= 1.0 / (1.0 + delta)
        if energy_check:
            # energy at integer step n uses pi_n = (pi_{n-1/2} + pi_{n+1/2}) / 2
            pim = 0.5 * (pi_prev.astype(np.float64) + pi)
            E, Epot = energy(phi, pim, tau)
            rate = -alpha / tau * float(np.sum(pim * pim)) * dx ** d + beta / tau * Epot
            if E0 is None:
                E0 = E
            else:
                W += 0.5 * dt * (rate + prev_rate)
            prev_rate = rate
            if n % max(1, nsteps // 200) == 0:
                echeck.append({"tau": tau, "E": E, "E0_plus_W": E0 + W,
                               "rel_resid": (E - E0 - W) / E})
        np.multiply(pi, dt, out=buf2)
        phi += buf2
    out = {"config": {k: (v if not isinstance(v, np.generic) else v.item()) for k, v in cfg.items()},
           "dx": dx, "dt": dt, "nsteps_planned": nsteps, "stopped_early": stopped_early,
           "mean_shift": shift, "wall_seconds": time.time() - t_start, "series": series}
    if energy_check:
        out["energy_check"] = echeck
        out["energy_max_abs_rel_resid"] = max(abs(e["rel_resid"]) for e in echeck) if echeck else None
        if echeck:
            emax = max(e["E"] for e in echeck)
            out["energy_max_abs_resid_over_maxE"] = max(abs(e["E"] - e["E0_plus_W"]) for e in echeck) / emax
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, cfg["name"] + ".json"), "w") as f:
            json.dump(out, f, indent=1)
    return out
