#!/usr/bin/env python3
"""Independent check (verifier v_charge) of the T=0 homogeneous charged-condensate phase structure.

Written from scratch from the action in BRIEF.md; does not import or copy module code.
Physical units with v = 1 (v = 1 TeV benchmark): couplings g=gX=0.1, lp=lam_phi=0.01, lX=lam_X=0.001.
Fields: phi real; X = (q1 + i q2)/sqrt2 ; |X|^2 = y ; R^2 = q1^2+q2^2 = 2y.
U = lp/4 (phi^2-1)^2 + g^2 phi^2 y + lX y^2 .
Rotating condensate X = |X| e^{i w t}: n = 2 w y, rho = n^2/(4y) + U.
"""
import json, os
import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar

os.environ.setdefault("OMP_NUM_THREADS", "2")
g, lp, lX = 0.1, 0.01, 0.001
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_t0.json")


def U(phi, y):
    return lp / 4 * (phi**2 - 1) ** 2 + g**2 * phi**2 * y + lX * y**2


def rho_full(n, phi, y):
    return n**2 / (4 * y) + U(phi, y)


# ---- method 1: brute force 2D minimisation over (phi, y) at fixed n ----
def rho_bruteforce(n):
    best = None
    # dense grid start
    phis = np.linspace(0, 1.05, 211)
    ys = np.exp(np.linspace(np.log(1e-4), np.log(10), 600))
    P, Y = np.meshgrid(phis, ys, indexing="ij")
    R = rho_full(n, P, Y)
    idx = np.unravel_index(np.argmin(R), R.shape)
    starts = [(phis[idx[0]], ys[idx[1]]), (0.0, 1.0), (0.9, 0.05), (0.5, 0.3)]
    for p0, y0 in starts:
        res = minimize(lambda z: rho_full(n, z[0], np.exp(z[1])), [p0, np.log(y0)],
                       method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000))
        if best is None or res.fun < best.fun:
            best = res
    return best.fun, abs(best.x[0]), np.exp(best.x[1])


# ---- method 2: my own reduced analysis ----
c1 = g**4 / lp - lX  # effective attractive quartic after integrating phi out on branch 1


def W(y):
    """min over phi of U at fixed y."""
    y = np.asarray(y, float)
    yb = lp / (2 * g**2)
    return np.where(y < yb, g**2 * y - c1 * y**2, lp / 4 + lX * y**2)


def dW(y):
    yb = lp / (2 * g**2)
    return np.where(y < yb, g**2 - 2 * c1 * y, 2 * lX * y)


def d2W(y):
    yb = lp / (2 * g**2)
    return np.where(y < yb, -2 * c1, 2 * lX)


yb = lp / (2 * g**2)                   # |X|_c^2 where phi_min -> 0 at fixed |X|
# branch 2 (phi = 0) stationary: w^2 = 2 lX y ; n = 2 y w
n_c = 2 * yb * np.sqrt(2 * lX * yb)
# branch 1 radial stability: d2/dy2 [n^2/(4y) + W] = 2 W'/y + W'' > 0
y_s1 = brentq(lambda y: 2 * dW(y) / y + d2W(y), 1e-6, yb - 1e-9)
n_s1 = 2 * y_s1 * np.sqrt(dW(y_s1))


def branch1_rho(n):
    # stationary y on stable part (y < y_s1): n = 2 y sqrt(g^2 - 2 c1 y)
    y = brentq(lambda y: 2 * y * np.sqrt(g**2 - 2 * c1 * y) - n, 1e-14, y_s1)
    w = np.sqrt(g**2 - 2 * c1 * y)
    return n**2 / (4 * y) + W(y), y, w


def branch2_rho(n):
    y = (n / np.sqrt(8 * lX)) ** (2.0 / 3.0)
    w = np.sqrt(2 * lX * y)
    return n**2 / (4 * y) + lp / 4 + lX * y**2, y, w


n_t = brentq(lambda n: branch1_rho(n)[0] - branch2_rho(n)[0], n_c * 1.0000001, n_s1 * 0.9999999)
# Q-matter: minimise rho/n on branch 2 ; equivalently p(mu)=0 at mu = (lp lX)^(1/4)
res = minimize_scalar(lambda n: branch2_rho(n)[0] / n, bounds=(0.05, 1.0), method="bounded",
                      options=dict(xatol=1e-12))
n_star_num = res.x
mstar_num = res.fun
mstar_closed = (lp * lX) ** 0.25
y_star = mstar_closed**2 / (2 * lX)
n_star_closed = 2 * mstar_closed * y_star


# sound speed from mu(n) = d rho/dn by finite differences on each branch
def cs2_branch(fn, n, h=1e-7):
    mu = lambda nn: (fn(nn + h * nn)[0] - fn(nn - h * nn)[0]) / (2 * h * nn)
    m0 = mu(n)
    dm = (mu(n * (1 + 1e-4)) - mu(n * (1 - 1e-4))) / (2e-4 * n)
    return n * dm / m0, m0


cs2_b1 = {}
for n in [1e-4, 1e-3, 0.005, 0.01, 0.02, 0.0220085, 0.03, 0.035, 0.04, 0.042, 0.0425]:
    if n < n_s1:
        cs2_b1[f"{n:g}"] = cs2_branch(branch1_rho, n)[0]
cs2_b2 = {}
for n in [n_c * 1.001, 0.05, 0.1, n_star_closed, 1.0]:
    cs2_b2[f"{n:g}"] = cs2_branch(branch2_rho, n)[0]


# ---- linear dispersion: rotating-frame real-field perturbations (derived independently) ----
def growth_rates(y, phi0, kvals):
    R2 = 2 * y
    R = np.sqrt(R2)
    w = np.sqrt(g**2 * phi0**2 + lX * R2)
    out = []
    for k in kvals:
        # state (da, db, dphi, da', db', dphi')
        Mmat = np.zeros((6, 6))
        Mmat[0:3, 3:6] = np.eye(3)
        # da'' = 2 w db' - (k^2 + 2 lX R^2) da - 2 g^2 phi0 R dphi
        Mmat[3, 4] = 2 * w
        Mmat[3, 0] = -(k**2 + 2 * lX * R2)
        Mmat[3, 2] = -2 * g**2 * phi0 * R
        # db'' = -2 w da' - k^2 db
        Mmat[4, 3] = -2 * w
        Mmat[4, 1] = -k**2
        # dphi'' = -(k^2 + lp(3 phi0^2 - 1) + g^2 R^2) dphi - 2 g^2 phi0 R da
        Mmat[5, 2] = -(k**2 + lp * (3 * phi0**2 - 1) + g**2 * R2)
        Mmat[5, 0] = -2 * g**2 * phi0 * R
        ev = np.linalg.eigvals(Mmat)
        out.append(ev.real.max())
    return np.array(out)


def branch1_state_from_x(x):
    # Codex variable x = F^2 = 2 y (v=1)
    y = x / 2
    phi0 = np.sqrt(1 - 2 * g**2 * y / lp)
    return y, phi0


unit_k = np.sqrt(lp)  # Codex/charge_phase unit: sqrt(lp) v = 0.1
lin = {}
for x in [0.01, 0.1, 0.25, 0.5]:
    y, phi0 = branch1_state_from_x(x)
    kb = np.linspace(1e-4, 1.2, 4001) * unit_k
    gr = growth_rates(y, phi0, kb)
    i = np.argmax(gr)
    unstable = kb[gr > 1e-12]
    # small-k slope: gamma ~ sqrt(-cs2) k
    ks = np.array([1e-4, 2e-4]) * unit_k
    gs = growth_rates(y, phi0, ks)
    lin[f"x={x}"] = dict(gamma_max_units=float(gr[i] / unit_k), k_at_max_units=float(kb[i] / unit_k),
                         k_unstable_max_units=float(unstable.max() / unit_k) if unstable.size else 0.0,
                         cs2_from_small_k=float(-(gs[1] / ks[1]) ** 2))
# branch 2 check: stable?
br2 = {}
for y in [0.5, 0.75, 1.58114, 3.0]:
    kb = np.linspace(1e-4, 2.0, 2001) * unit_k
    gr = growth_rates(y, 0.0, kb)
    # sound speed: smallest |imag| eigenfrequency at small k
    k = 1e-4 * unit_k
    R2 = 2 * y
    w = np.sqrt(lX * R2)
    Mmat = np.zeros((6, 6)); Mmat[0:3, 3:6] = np.eye(3)
    Mmat[3, 4] = 2 * w; Mmat[3, 0] = -(k**2 + 2 * lX * R2)
    Mmat[4, 3] = -2 * w; Mmat[4, 1] = -k**2
    Mmat[5, 2] = -(k**2 - lp + g**2 * R2)
    ev = np.linalg.eigvals(Mmat)
    om = np.sort(np.abs(ev.imag))
    om_small = om[om > 0][0] if np.any(om > 0) else 0.0
    br2[f"y={y}"] = dict(max_growth=float(gr.max()), cs2_from_small_k=float((om_small / k) ** 2),
                         phi_mass2_over_lp=float((g**2 * R2 - lp) / lp))

# Codex production state
xP = 0.25
yP, phiP = branch1_state_from_x(xP)
wP = np.sqrt(g**2 - 2 * c1 * yP)
nP = 2 * yP * wP
cs2P = cs2_branch(branch1_rho, nP)[0]
kb = np.linspace(1e-4, 1.0, 20001) * unit_k
grP = growth_rates(yP, phiP, kb)
iP = np.argmax(grP)
unstP = kb[grP > 1e-12].max() / unit_k
L_box = 64.0
n_harm = int(np.floor(unstP * L_box / (2 * np.pi)))

# brute-force comparison
bf = []
for n in np.concatenate([np.linspace(0.002, 0.06, 15), [0.1, n_star_closed, 0.3, 1.0]]):
    rb, pb, ybf = rho_bruteforce(n)
    cands = [branch2_rho(n)[0] if n >= n_c else np.inf]
    if n < n_s1:
        cands.append(branch1_rho(n)[0])
    ra = min(cands)
    bf.append(dict(n=float(n), rho_bruteforce=float(rb), rho_reduced=float(ra), rel=float(abs(rb - ra) / ra),
                   phi=float(pb), y=float(ybf)))

# rho/n table
def rho_over_n(n, br):
    r = (branch1_rho if br == 1 else branch2_rho)(n)[0]
    return r / n / g  # in units of m_X = g v
ron = dict(prod_b1=rho_over_n(nP, 1), nc_b1=rho_over_n(n_c, 1), nc_b2=rho_over_n(n_c * (1 + 1e-12), 2),
           nt=rho_over_n(n_t, 2), n0p1_b2=rho_over_n(0.1, 2), nstar=rho_over_n(n_star_closed, 2),
           n1_b2=rho_over_n(1.0, 2))
# S jump at n_t
y1t = branch1_rho(n_t)[1]
y2t = branch2_rho(n_t)[1]

out = dict(
    note="verifier v_charge independent T=0 computation; v=1 units; densities in v^3",
    yb_absX2_c=yb, n_c=n_c, y_s1=y_s1, x_s1=2 * y_s1, n_s1=n_s1, n_t=n_t,
    S_jump_at_nt=dict(S_branch1=float(np.sqrt(1 - 2 * g**2 * y1t / lp)), x_branch1=2 * y1t, x_branch2=2 * y2t),
    n_star_numeric=n_star_num, n_star_closed=n_star_closed, mstar_numeric_v=mstar_num,
    mstar_closed_v=mstar_closed, mstar_GeV_at_v1TeV=mstar_closed * 1000, x_star=2 * y_star,
    absX_star=np.sqrt(y_star),
    cs2_branch1=cs2_b1, cs2_branch2=cs2_b2, linear_branch1=lin, linear_branch2=br2,
    codex_production=dict(x=xP, S=phiP, Omega_units=wP / unit_k, nbar_units=nP / unit_k**3 / 10 * 10,
                          n_v3=nP, n_over_nc=nP / n_c, cs2_thermo=cs2P,
                          gamma_max_units=float(grP[iP] / unit_k), k_max_units=float(kb[iP] / unit_k),
                          wavelength=float(2 * np.pi / (kb[iP] / unit_k)), k_unstable_edge=float(unstP),
                          harmonics_per_axis_L64=n_harm,
                          t_nonlin_ln200=float(np.log(1 / 0.005) / (grP[iP] / unit_k))),
    bruteforce=bf, max_rel_bruteforce=max(b["rel"] for b in bf), rho_over_n_mX=ron,
)
with open(OUT, "w") as f:
    json.dump(out, f, indent=1, default=float)
for k in ["yb_absX2_c", "n_c", "x_s1", "n_s1", "n_t", "S_jump_at_nt", "n_star_closed", "n_star_numeric",
          "mstar_closed_v", "x_star", "max_rel_bruteforce", "rho_over_n_mX"]:
    print(k, out[k])
print("cs2 b1", cs2_b1)
print("cs2 b2", cs2_b2)
print("lin b1", lin)
print("lin b2", br2)
print("prod", out["codex_production"])
