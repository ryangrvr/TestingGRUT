"""Isolated closed bubble control: a -v domain of comoving radius R0 inside the +v vacuum, released from rest at
conformal time tau_i in a radiation-dominated background (a = tau), unchanged phi potential, X = D = 0.

Three levels:
  (1) Nambu-Goto thin-wall ODE in comoving coordinates (derived in NOTES.md):
          R'' = -(1 - R'^2) [kappa / R + D_h R' / tau],
      kappa = number of curved directions (1: circle/cylinder, 2: sphere),
      D_h = alpha + beta/2 = effective Hubble damping of the wall (3 for the physical 3+1 equation; alpha for PRS).
      In radiation era the ODE is invariant under (tau, R) -> (s tau, s R): tau_c / tau_i depends only on x = R0/tau_i.
  (2) Field-theory radial solver (finite-volume Laplacian, leapfrog) for the same symmetric bubble.
  (3) Cartesian 2D/3D lattice checks of (2) with the network code's integrator.
Writes results/bubbles.json.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wallsim as ws  # noqa: E402

RES = os.path.join(HERE, "..", "results")
LAM, DELTA = ws.LAM, ws.DELTA


# ------------------------------------------------------------------------------------------------
# (1) Nambu-Goto
# ------------------------------------------------------------------------------------------------
def ng_collapse(x, kappa, D, tau_i=1.0, flat=False, samples=None):
    """Collapse conformal time of a thin wall released from rest at tau_i with comoving radius R0 = x tau_i.
    flat=True drops the Hubble term (Minkowski). Returns tau_c (and R at requested sample times)."""
    R0 = x * tau_i

    def rhs(t, y):
        R, u = y
        damp = 0.0 if flat else D * u / t
        return [u, -(1.0 - u * u) * (kappa / R + damp)]

    def ev(t, y):
        return y[0] - 1e-7 * R0
    ev.terminal = True
    ev.direction = -1

    def ev10(t, y):
        return y[0] - 0.1 * R0
    ev10.direction = -1
    sol = solve_ivp(rhs, (tau_i, tau_i + 50 * R0 + 50 * tau_i), [R0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13 * R0,
                    events=[ev, ev10], dense_output=True)
    tc = float(sol.t_events[0][0]) if sol.t_events[0].size else float("nan")
    t10 = float(sol.t_events[1][0]) if sol.t_events[1].size else float("nan")
    out = {"tau_c": tc, "tau_10": t10}
    if samples is not None:
        out["R_at"] = {str(s): (float(sol.sol(s)[0]) if s < tc else 0.0) for s in samples}
    return out


# ------------------------------------------------------------------------------------------------
# (2) radial field solver
# ------------------------------------------------------------------------------------------------
def radial_bubble(kappa, alpha, beta, tau_i, R0, tau_end, cells_per_delta=5.0, courant=0.25, tau_a=1.0,
                  n_track=3000):
    """Evolve phi(r) with phi'' + alpha/tau phi' - lap_kappa phi + a^beta dU = 0, a = tau/tau_a.
    Wall width (comoving) = DELTA / a^(beta/2); dr chosen to resolve it with cells_per_delta at tau_end."""
    a_end = tau_end / tau_a
    w_end = DELTA / a_end ** (beta / 2)
    dr = w_end / cells_per_delta
    r_max = R0 + (tau_end - tau_i) + 6 * DELTA / (tau_i / tau_a) ** (beta / 2) + 2.0
    M = int(math.ceil(r_max / dr))
    rf = (np.arange(M) + 1.0) * dr            # outer face of cell j
    rc = (np.arange(M) + 0.5) * dr            # cell centre
    Af = rf ** kappa
    Vc = (rf ** (kappa + 1) - (rf - dr) ** (kappa + 1)) / (kappa + 1)
    a_i = tau_i / tau_a
    phi = np.tanh((rc - R0) * a_i ** (beta / 2) / DELTA)
    pi = np.zeros(M)
    flux = np.zeros(M)
    lap = np.zeros(M)
    dt = courant * dr
    nsteps = int(math.ceil((tau_end - tau_i) / dt))
    every = max(1, nsteps // n_track)
    track = []
    tau_c = None
    reform_max_R = 0.0
    last_neg_tau = None
    failed = None
    t0 = time.time()
    for n in range(nsteps + 1):
        tau = tau_i + n * dt
        if n % every == 0 or n == nsteps:
            neg = np.nonzero(phi < 0)[0]
            if neg.size:
                j = neg[-1]
                last_neg_tau = tau
                # linear interpolation of the outermost zero crossing
                if j + 1 < M:
                    R = rc[j] + dr * phi[j] / (phi[j] - phi[j + 1])
                else:
                    R = rc[j]
            else:
                R = 0.0
            track.append((tau, R))
            w_now = DELTA / (tau / tau_a) ** (beta / 2)
            # first arrival of the wall at the centre: outermost zero within one comoving wall width of r = 0
            if tau_c is None and R < w_now:
                tau_c = tau
            elif tau_c is not None and R > reform_max_R:
                reform_max_R = R
        if n == nsteps:
            break
        a = tau / tau_a
        flux[:-1] = Af[:-1] * (phi[1:] - phi[:-1]) / dr
        flux[-1] = 0.0
        lap[0] = flux[0] / Vc[0]
        lap[1:] = (flux[1:] - flux[:-1]) / Vc[1:]
        F = lap - (a ** beta) * LAM * phi * (phi * phi - 1.0)
        dl = 0.5 * alpha * dt / tau
        pi = ((1.0 - dl) * pi + dt * F) / (1.0 + dl)
        phi = phi + dt * pi
        if n % 1000 == 0 and not np.all(np.isfinite(phi)):
            failed = tau
            break
    # refine tau_c by linear interpolation of the last two tracked points before R = 0
    tr = np.array(track)
    below = np.nonzero(tr[:, 1] <= 0.1 * R0)[0]
    tau_10 = float(tr[below[0], 0]) if below.size else None
    return {"kappa": kappa, "alpha": alpha, "beta": beta, "tau_i": tau_i, "R0": R0, "x": R0 / tau_i,
            "numerical_failure_tau": failed, "tau_10": tau_10,
            "tau_end": tau_end, "dr": dr, "M": M, "dt": dt, "tau_c": tau_c, "collapsed": tau_c is not None,
            "negative_cells_at_end": bool(np.any(phi < 0)),
            "tau_last_negative": last_neg_tau,
            "phi_centre_end": float(phi[0]), "phi_min_end": float(phi.min()),
            "wall_cells_per_delta_end": cells_per_delta, "seconds": time.time() - t0,
            "track_tau": tr[:, 0].tolist()[:: max(1, len(tr) // 60)],
            "track_R": tr[:, 1].tolist()[:: max(1, len(tr) // 60)],
            "_track": tr}


def rebound(res, tc_ng):
    """Largest negative-region radius after the thin-wall collapse time (+0.1 R0), i.e. collapse debris."""
    tr = res["_track"]
    sel = tr[:, 0] > tc_ng + 0.1 * res["R0"]
    mx = float(tr[sel, 1].max()) if sel.any() else 0.0
    out = {"debris_max_R_after_NG_collapse": mx, "debris_max_R_over_R0": mx / res["R0"]}
    if res["tau_last_negative"] is not None:
        out["debris_lifetime_over_R0"] = (res["tau_last_negative"] - tc_ng) / res["R0"]
    return out


def R_at(res, tau):
    tr = res["_track"]
    if tau <= tr[0, 0]:
        return float(tr[0, 1])
    return float(np.interp(tau, tr[:, 0], tr[:, 1]))


# ------------------------------------------------------------------------------------------------
# (3) Cartesian lattice check
# ------------------------------------------------------------------------------------------------
def cart_bubble(d, alpha, beta, tau_i, R0, tau_end, N, L, courant=1 / 3, tau_a=1.0, dtype=np.float32):
    dx = L / N
    x = (np.arange(N) - N / 2 + 0.5) * dx
    g = np.meshgrid(*([x] * d), indexing="ij")
    r = np.sqrt(sum(gg * gg for gg in g))
    del g
    a_i = tau_i / tau_a
    phi = np.tanh((r - R0) * a_i ** (beta / 2) / DELTA).astype(dtype)
    del r
    pi = np.zeros_like(phi)
    buf = np.empty_like(phi)
    buf2 = np.empty_like(phi)
    dt = courant * dx
    nsteps = int(math.ceil((tau_end - tau_i) / dt))
    track = []
    tau_c = None
    last_neg = None
    cell = dx ** d
    for n in range(nsteps + 1):
        tau = tau_i + n * dt
        if n % 2 == 0 or n == nsteps:
            vneg = float(np.count_nonzero(phi < 0)) * cell
            Reff = (vneg / math.pi) ** 0.5 if d == 2 else (3 * vneg / (4 * math.pi)) ** (1 / 3)
            track.append((tau, Reff))
            if vneg > 0:
                last_neg = tau
            if tau_c is None and Reff < DELTA / (tau / tau_a) ** (beta / 2):
                tau_c = tau
        if n == nsteps:
            break
        a = tau / tau_a
        ws.lap_inplace(phi, buf)
        buf *= 1.0 / dx ** 2
        ws.dU(phi, out=buf2)
        if beta:
            buf2 *= a ** beta
        buf -= buf2
        dl = 0.5 * alpha * dt / tau
        pi *= (1 - dl)
        buf *= dt
        pi += buf
        pi *= 1 / (1 + dl)
        np.multiply(pi, dt, out=buf2)
        phi += buf2
    tr = np.array(track)
    below = np.nonzero(tr[:, 1] <= 0.1 * R0)[0]
    return {"d": d, "alpha": alpha, "beta": beta, "tau_i": tau_i, "R0": R0, "N": N, "L": L, "dx": dx,
            "tau_10": float(tr[below[0], 0]) if below.size else None,
            "wall_cells_per_delta_end": DELTA / (tau_end / tau_a) ** (beta / 2) / dx, "tau_c": tau_c,
            "tau_last_negative": last_neg,
            "track_tau": tr[:, 0].tolist()[:: max(1, len(tr) // 40)],
            "track_R": tr[:, 1].tolist()[:: max(1, len(tr) // 40)]}


def main(quick=False):
    out = {"conventions": {
        "units": "v = 1, m_phi = 1, a = tau (radiation era), wall tanh half-width DELTA = 2/m_phi physical",
        "x": "R0/tau_i = comoving bubble radius over comoving Hubble radius (a/a')=tau at release",
        "D_h": "effective Hubble damping of the wall = alpha + beta/2",
        "kappa": "number of curved directions (1 circle or z-invariant cylinder, 2 sphere)"}}
    # ---- (1) NG
    ng = {}
    flat = {}
    for kappa in (1, 2):
        flat[str(kappa)] = (ng_collapse(1.0, kappa, 0, tau_i=1e6, flat=True)["tau_c"] - 1e6) / 1e6
    out["ng_flat_collapse_over_R0"] = {"kappa1": flat["1"], "kappa1_exact": math.pi / 2,
                                       "kappa2": flat["2"], "kappa2_exact": 1.3110287771460600}
    xs = np.geomspace(1e-3, 1e3, 25)
    for kappa, D in ((1, 2), (1, 3), (2, 3)):
        rows = []
        for x in xs:
            tc = ng_collapse(float(x), kappa, D)["tau_c"]
            rows.append({"x": float(x), "tau_c_over_tau_i": tc, "dtau_over_R0": (tc - 1.0) / x, "tau_c_over_R0": tc / x})
        ng[f"kappa{kappa}_D{D}"] = rows
    out["ng_rd_scan"] = ng
    # super-horizon asymptote and freezing diagnostics
    asym = {}
    for kappa, D in ((1, 2), (1, 3), (2, 3)):
        big = ng_collapse(1e4, kappa, D, samples=[0.5e4, 1e4])
        asym[f"kappa{kappa}_D{D}"] = {"tau_c_over_R0_x1e4": big["tau_c"] / 1e4,
                                      "R_over_R0_at_tau_eq_R0_over_2": big["R_at"]["5000.0"] / 1e4,
                                      "R_over_R0_at_tau_eq_R0": big["R_at"]["10000.0"] / 1e4}
    out["ng_superhorizon"] = asym
    # ---- (2) field radial, physical equation (alpha=2, beta=2 -> D_h = 3)
    rad = []
    taus_i = [4.0] if quick else [4.0, 8.0]
    xlist = [0.5, 1, 2, 4] if quick else [0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    for kappa in (1, 2):
        for tau_i in taus_i:
            for x in xlist:
                R0 = x * tau_i
                if R0 < 4 * DELTA / tau_i:   # require R0 >= 4 comoving wall widths at release
                    continue
                tcng = ng_collapse(x, kappa, 3, tau_i=tau_i)["tau_c"]
                if tcng > 110:
                    continue
                tau_end = 1.6 * tcng + 4.0
                r = radial_bubble(kappa, 2, 2, tau_i, R0, tau_end, cells_per_delta=5.0)
                r["tau_c_NG"] = tcng
                r["tau_10_NG"] = ng_collapse(x, kappa, 3, tau_i=tau_i)["tau_10"]
                r.update(rebound(r, tcng))
                r["R_over_R0_at_tau_R0"] = R_at(r, R0) / R0 if R0 > tau_i else None
                r["R_over_R0_at_tau_half_R0"] = R_at(r, 0.5 * R0) / R0 if 0.5 * R0 > tau_i else None
                r.pop("_track")
                print(f"radial phys kappa={kappa} tau_i={tau_i} x={x} tau_10={r['tau_10']} NG10={r['tau_10_NG']:.3f} "
                      f"M={r['M']} {r['seconds']:.1f}s", flush=True)
                rad.append(r)
    # resolution check of the radial solver: one case at 2x cells
    for kappa in (1, 2):
        tau_i, x = 4.0, 2.0
        tcng = ng_collapse(x, kappa, 3, tau_i=tau_i)["tau_c"]
        for cpd in (2.5, 10.0):
            r = radial_bubble(kappa, 2, 2, tau_i, x * tau_i, 1.6 * tcng + 4.0, cells_per_delta=cpd)
            r["tau_c_NG"] = tcng
            r["tau_10_NG"] = ng_collapse(x, kappa, 3, tau_i=tau_i)["tau_10"]
            r.update(rebound(r, tcng))
            r["resolution_check"] = True
            r.pop("_track")
            rad.append(r)
    out["radial_physical"] = rad
    # ---- (2') field radial, PRS (beta=0): sphere alpha=3, circle alpha=2 and 3
    prs = []
    xl = [0.4, 1, 4] if quick else [0.2, 0.4, 1, 2, 4, 8, 16]
    for kappa, alpha in ((2, 3), (1, 2), (1, 3)):
        tau_i = 50.0
        for x in xl:
            R0 = x * tau_i
            tcng = ng_collapse(x, kappa, alpha, tau_i=tau_i)["tau_c"]
            r = radial_bubble(kappa, alpha, 0, tau_i, R0, 1.6 * tcng + 10.0, cells_per_delta=4.0, courant=0.4)
            r["tau_c_NG"] = tcng
            r["tau_10_NG"] = ng_collapse(x, kappa, alpha, tau_i=tau_i)["tau_10"]
            r.update(rebound(r, tcng))
            r["R_over_R0_at_tau_R0"] = R_at(r, R0) / R0 if R0 > tau_i else None
            r.pop("_track")
            print(f"radial PRS kappa={kappa} alpha={alpha} x={x} tau_10={r['tau_10']} NG10={r['tau_10_NG']:.2f} {r['seconds']:.1f}s",
                  flush=True)
            prs.append(r)
    out["radial_prs"] = prs
    # ---- (3) Cartesian checks (physical equation)
    cart = []
    cases = [(2, 4.0, 4.0, 1024, 64.0), (3, 4.0, 4.0, 160, 16.0), (3, 4.0, 2.0, 128, 12.8)]
    if not quick:
        cases.append((2, 4.0, 16.0, 2048, 48.0))
    for d, tau_i, R0, N, L in cases:
        kappa = d - 1
        tcng = ng_collapse(R0 / tau_i, kappa, 3, tau_i=tau_i)["tau_c"]
        tau_end = 1.15 * tcng + 1.0
        c = cart_bubble(d, 2, 2, tau_i, R0, tau_end, N, L)
        rr = radial_bubble(kappa, 2, 2, tau_i, R0, tau_end, cells_per_delta=5.0)
        c["tau_c_radial"] = rr["tau_c"]
        c["tau_10_radial"] = rr["tau_10"]
        c["tau_c_NG"] = tcng
        c["tau_10_NG"] = ng_collapse(R0 / tau_i, kappa, 3, tau_i=tau_i)["tau_10"]
        print(f"cart d={d} R0={R0} tau_10={c['tau_10']} radial={rr['tau_10']} NG={c['tau_10_NG']:.3f} cells/delta_end="
              f"{c['wall_cells_per_delta_end']:.2f}", flush=True)
        cart.append(c)
    out["cartesian_checks"] = cart
    fn = "bubbles_quick.json" if quick else "bubbles.json"
    with open(os.path.join(RES, fn), "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    main(quick=(len(sys.argv) > 1 and sys.argv[1] == "quick"))
