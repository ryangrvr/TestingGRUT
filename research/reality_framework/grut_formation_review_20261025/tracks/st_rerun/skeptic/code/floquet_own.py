"""st_rerun_skeptic: own Floquet check of the auxiliary quartic orbit (independent of Codex's orbit_diagnostic).

Input data only: Codex's matched initial chi (QUARTIC_SCAN.json 'initial_chi'; lambda = 1e-3).  Own method:
orbit invariants (E, L); radial period from quadrature between the turning points; monodromy integrated from
apoapsis with scipy Radau (Codex used DOP853 from an event-located start); precession removed by rotating with the
integrated end angle; mu = ln|max multiplier| / T.  Edges located by bisection on mu > 1e-6.
Usage: python3 -I floquet_own.py <QUARTIC_SCAN.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
pub = json.load(open(sys.argv[1])); LAM = 1e-3
x1, x2, v1, v2 = pub['initial_chi']
E = .5 * (v1 * v1 + v2 * v2) + LAM / 4 * (x1 * x1 + x2 * x2) ** 2; L = x1 * v2 - x2 * v1
pr2 = lambda r: 2 * E - L * L / r ** 2 - LAM * r ** 4 / 2
r0 = np.hypot(x1, x2)
rmin = brentq(pr2, 1e-6, r0 * (1 + 1e-12)) if pr2(r0) >= 0 else None
rmax = brentq(pr2, r0, 10 * r0)
c, h = (rmin + rmax) / 2, (rmax - rmin) / 2
T = 2 * quad(lambda s: h * np.cos(s) / np.sqrt(max(pr2(c + h * np.sin(s)), 1e-300)), -np.pi / 2, np.pi / 2, epsabs=0, epsrel=1e-13, limit=400)[0]


def rhs(t, y, k):
    q, M = y[:2], y[4:].reshape(4, 4); rr = q @ q
    hs = k * k * np.eye(2) + LAM * (rr * np.eye(2) + 2 * np.outer(q, q))
    return np.r_[y[2:4], -LAM * rr * q, np.vstack((M[2:], -hs @ M[:2])).ravel()]


def mu(k):
    s = solve_ivp(rhs, (0, T), np.r_[rmax, 0, 0, L / rmax, np.eye(4).ravel()], args=(k,), method='Radau', rtol=1e-12, atol=1e-13)
    y = s.y[:, -1]; ang = np.arctan2(y[1], y[0])
    R = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
    Mon = np.kron(np.eye(2), R) @ y[4:].reshape(4, 4)
    return max(0., float(np.log(abs(np.linalg.eigvals(Mon)).max()) / T))


out = dict(E=E, L=L, rmin=rmin, rmax=rmax, period=T, published_period=pub['radial_period'])
out['IR'] = {str(k): mu(k) / k for k in (.01, .02, .05)}
ks = np.round(np.arange(2.0, 2.26, .01), 3); m = [mu(k) for k in ks]
out['peak'] = dict(k=float(ks[int(np.argmax(m))]), mu=float(max(m)), mu_213=mu(2.13))
f = lambda k: mu(k) - 1e-6
out['edge1'] = brentq(f, 3.10, 3.16, xtol=1e-5)
out['band2_low'] = brentq(f, 5.50, 5.60, xtol=1e-5)
out['band2_high'] = brentq(f, 5.90, 6.00, xtol=1e-5)
k2 = np.round(np.arange(5.74, 5.80, .005), 3); m2 = [mu(k) for k in k2]
out['peak2'] = dict(k=float(k2[int(np.argmax(m2))]), mu=float(max(m2)))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print(json.dumps(out, indent=1))
