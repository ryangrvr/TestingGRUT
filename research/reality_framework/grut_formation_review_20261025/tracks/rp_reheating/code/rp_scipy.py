#!/usr/bin/env python3
"""rp_reheating track (Claude, own code): a second, separately written formulation of the same equations in Python,
with H EVOLVED by the Raychaudhuri equation (as Codex does) and the Friedmann constraint monitored, integrated with
SciPy DOP853 and Radau.  Used (i) to cross-check my C engine (which reconstructs H from the constraint) and (ii) to
compare with Codex's recorded 19-variable endpoints at tau = 10 (INDEPENDENT_REHEATING.json short_C / Radau_endpoint).
State: [ln a, i, j, X re, X im, Y re, Y im, S re, S im, vX re, vX im, vY re, vY im, vS re, vS im, r, H].
Usage: python3 -B rp_scipy.py out.json [end=10] [b=2] [gamma=0.1]
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rp

P, XI, ETA2 = rp.P, rp.XI, rp.ETA2


def rhs(t, u, b, gam, xi=XI, eta2=ETA2):
    lna, i, j = u[0], u[1], u[2]
    z = u[3:9:2] + 1j * u[4:9:2]; v = u[9:15:2] + 1j * u[10:15:2]
    r, H = u[15], u[16]
    V, F = rp.potential(z)
    K = float(np.sum(np.abs(v)**2))
    x2 = abs(z[0])**2
    f = 1 - 2 * xi * eta2 * x2
    fd = -4 * xi * eta2 * np.real(np.conj(z[0]) * v[0])
    rhoI = 0.5 * j * j + 0.5 * b * b * i * i; pI = 0.5 * j * j - 0.5 * b * b * i * i
    R = ((rhoI - 3 * pI) + eta2 * (4 * V - 2 * K + 12 * xi * abs(v[0])**2 - 12 * xi * np.real(np.conj(z[0]) * F[0]))) \
        / (f + 12 * xi**2 * eta2 * x2)
    acc = -3 * H * v - F
    acc[0] -= xi * R * z[0]
    fdd = -4 * xi * eta2 * (abs(v[0])**2 + np.real(np.conj(z[0]) * acc[0]))
    Hd = -(j * j + 2 * eta2 * K + 4 * r / 3 + fdd - H * fd) / (2 * f)
    d = np.empty(17)
    d[0] = H; d[1] = j; d[2] = -(3 * H + gam) * j - b * b * i
    d[3:9:2] = v.real; d[4:9:2] = v.imag; d[9:15:2] = acc.real; d[10:15:2] = acc.imag
    d[15] = gam * j * j - 4 * H * r; d[16] = Hd
    return d


def friedmann_residual(u, b, xi=XI, eta2=ETA2):
    i, j = u[1], u[2]
    z = u[3:9:2] + 1j * u[4:9:2]; v = u[9:15:2] + 1j * u[10:15:2]
    r, H = u[15], u[16]
    V, F = rp.potential(z); K = float(np.sum(np.abs(v)**2))
    f = 1 - 2 * xi * eta2 * abs(z[0])**2; fd = -4 * xi * eta2 * np.real(np.conj(z[0]) * v[0])
    rho = 0.5 * j * j + 0.5 * b * b * i * i + eta2 * (K + V) + r
    return (3 * f * H * H + 3 * H * fd - rho) / rho


def initial(b, delta=0.0):
    z, info = rp.prepare(b, delta)
    u = np.zeros(17)
    u[3:9:2] = z.real; u[4:9:2] = z.imag
    u[1] = np.sqrt(2 * info['rhoI0']) / b; u[16] = 1.0
    return u, info


def to_codex_order(u):
    """Codex's state order: [i, j, X, Y, S (re/im), velocities, r, ln a, h, W_total, W_hidden]; the two work integrals
    are not evolved here."""
    return [u[1], u[2]] + list(u[3:15]) + [u[15], u[0], u[16]]


if __name__ == '__main__':
    out = Path(sys.argv[1]); end = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    b = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0; gam = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
    u0, info = initial(b)
    res = dict(b=b, gamma=gam, end=end, prep=info, u0=u0.tolist())
    for method, kw in (('DOP853', dict(rtol=1e-12, atol=1e-14)), ('Radau', dict(rtol=1e-10, atol=1e-13))):
        t0 = time.time()
        sol = solve_ivp(rhs, (0, end), u0, method=method, args=(b, gam), dense_output=False, **kw)
        u = sol.y[:, -1]
        fr = max(abs(friedmann_residual(sol.y[:, k], b)) for k in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 400)))
        res[method] = dict(success=bool(sol.success), nfev=int(sol.nfev), seconds=time.time() - t0, final=u.tolist(),
                           final_codex_order=to_codex_order(u), max_Friedmann_residual_sampled=float(fr))
        print(method, 'nfev', sol.nfev, '%.1fs' % (time.time() - t0), 'Friedmann', fr, flush=True)
    out.write_text(json.dumps(res, indent=1))
