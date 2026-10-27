"""st_physics skeptic, check k4: physical-variable linear response on the full homogeneous background.

Independence label: the background equations and the potential are the model definition, imported from Codex's
Python reheating.py (rhs, potential_force, prepared).  Everything else is ours:
  * the perturbation equations are written in PHYSICAL variables (delta z, delta zdot) per complex field,
        delta zddot = -3 H delta zdot - (k/a)^2 delta z - J_F(z) delta z - xi R P_X delta z,
    with J_F the real 6x6 Jacobian of the model force F = dV/dzbar, obtained by central finite differences
    (no hand-written Hessian, no Codex canonical/conformal Hessian, no Ward projection, no response.c);
  * the integrator is scipy DOP853 on the combined background + tangent system;
  * the gauge (Coulomb, e = 1e-7) terms are dropped: e^2/k^2 ~ 2e-15 at k ~ 2.
Conversion to Codex's conformal canonical variables for comparison: chi = sqrt2 a z, chi' = sqrt2 a^2 (zdot + H z).

Columns per k: the 12 unit vectors of the conformal canonical (q, p) basis at a = 1, mapped to physical initial data
(delta z, delta zdot) = (q / sqrt2, (p - H q) / sqrt2); plus the supplied preparation-phase seed integrated
explicitly at k = 2.13 as a linearity check.

usage: python3 -B k4_phys_tangent.py <codex reality_framework dir> <tag quiet3|quiet5> <tau_end> <out.npz>
"""
import sys, time
import numpy as np
from scipy.integrate import solve_ivp

RF, TAG, TEND, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
sys.path.insert(0, RF + '/grut_reheating_preparation_20261026/code')
import reheating as Rh  # noqa: E402  (model definition only)

BETA, DECAY = {'quiet3': (3., .5), 'quiet5': (5., 1.)}[TAG]
XI, ETA2 = Rh.XI, Rh.ETA2
KV = np.array([1.8, 2.0, 2.13, 2.4, 2.8])
NCOL_K = 12
H_FD = 1e-2
# 4-point central stencil: exact (up to rounding) for polynomials of degree <= 4; the model force is cubic.
_E = np.eye(6) * H_FD
_SHIFTS = np.concatenate([_E, -_E, 2 * _E, -2 * _E])   # (24, 6)


def jac_force(zr):
    Z = zr[None, :] + _SHIFTS                          # (24, 6)
    z = (Z[:, 0::2] + 1j * Z[:, 1::2]).T               # (3, 24): potential_force is elementwise in its rows
    _, F = Rh.potential_force(z)
    Fr = np.empty((24, 6)); Fr[:, 0::2] = F.real.T; Fr[:, 1::2] = F.imag.T
    p1, m1, p2, m2 = Fr[:6], Fr[6:12], Fr[12:18], Fr[18:24]
    return (8 * (p1 - m1) - (p2 - m2)).T / (12 * H_FD)  # J[i, j] = dF_i / dz_j


def ricci(v):
    z, vel, pot, F, f, fd, R, hd = Rh.geometric(v, BETA, XI, ETA2)
    return R


PX = np.diag([1., 1., 0, 0, 0, 0])


def make_rhs(ks):
    nk = len(ks)
    k2 = np.repeat(ks ** 2, 1)

    def rhs(t, y):
        v = y[:19]
        out = np.empty_like(y)
        out[:19] = Rh.rhs(t, v, BETA, DECAY, XI, ETA2)
        a = np.exp(v[15]); H = v[16]
        R = ricci(v)
        J = jac_force(v[2:8]) + XI * R * PX
        T = y[19:].reshape(-1, 12)            # rows: columns of all k; [dz(6), dzdot(6)]
        dz, dzd = T[:, :6], T[:, 6:]
        acc = -3 * H * dzd - dz @ J.T - (KCOL[:, None] ** 2 / a ** 2) * dz
        out[19:] = np.hstack([dzd, acc]).ravel()
        return out
    return rhs


def main():
    v0, res = Rh.prepared(BETA)
    # supplied seed: our own central difference of the model preparation in its phase argument
    vp, _ = Rh.prepared(BETA, delta=1e-5); vm, _ = Rh.prepared(BETA, delta=-1e-5)
    dz_sup = (vp[2:8] - vm[2:8]) / 2e-5          # d z / d delta (real components), zdot derivative is 0
    a0, H0 = 1.0, v0[16]
    cols, colk = [], []
    for k in KV:
        for j in range(12):
            e = np.zeros(12); e[j] = 1.0
            q, p = e[:6], e[6:]
            cols.append(np.r_[q / np.sqrt(2) / a0, (p / a0 - H0 * q) / np.sqrt(2) / a0]); colk.append(k)
    cols.append(np.r_[dz_sup, np.zeros(6)]); colk.append(2.13)
    xrot = np.zeros(12); xrot[0] = -v0[3]; xrot[1] = v0[2]     # i x : delta z_X = i z_X, zdot unchanged (0)
    cols.append(xrot); colk.append(2.13)
    global KCOL
    KCOL = np.array(colk)
    y0 = np.r_[v0, np.array(cols).ravel()]
    times = np.r_[np.arange(0, 80, 0.05), np.arange(80, TEND + 1e-9, 0.5)]
    times = times[times <= TEND]
    t0 = time.time()
    sol = solve_ivp(make_rhs(KV), (0, TEND), y0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=times,
                    max_step=0.08 / BETA)
    print(TAG, 'status', sol.status, sol.message, 'nfev', sol.nfev, f'{time.time() - t0:.0f}s', flush=True)
    bg = sol.y[:19].T
    T = sol.y[19:].T.reshape(len(sol.t), len(cols), 12)
    np.savez_compressed(OUT, times=sol.t, bg=bg, T=T, colk=KCOL, KV=KV, dz_sup=dz_sup, xrot=xrot,
                        prep_residual=res, v0=v0)


if __name__ == '__main__':
    main()
