"""Independent Floquet check of Codex 2cbc7f1's isolated quartic X orbit (QUARTIC_ORBIT.json).

Own construction: the orbit is rebuilt from its two invariants (E, L) of the matched tau = 10 state. The radial
period and the precession angle come from quadrature between the turning points, not from an event search. The
monodromy is integrated with Radau from apoapsis, and the precession is removed by a rotation.
usage: python quartic_floquet.py <QUARTIC_ORBIT.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

src, out = sys.argv[1:3]
pub = json.load(open(src))
LAM = 1e-3
x1, x2, v1, v2 = pub['initial_chi']
E = 0.5 * (v1 * v1 + v2 * v2) + LAM / 4 * (x1 * x1 + x2 * x2) ** 2
Lz = x1 * v2 - x2 * v1


def pr2(r):  # twice the radial kinetic energy: 2E - L^2/r^2 - lam r^4/2
    return 2 * E - Lz * Lz / r ** 2 - LAM * r ** 4 / 2


r0 = np.hypot(x1, x2)
rmin = brentq(pr2, 1e-6, r0 if pr2(r0) > 0 else r0 * (1 + 1e-12))
rmax = brentq(pr2, r0, 10 * r0)
# substitution r = (rmin+rmax)/2 + (rmax-rmin)/2 sin(s) removes the endpoint singularities
c, h = (rmin + rmax) / 2, (rmax - rmin) / 2


def integrand(s, f):
    r = c + h * np.sin(s)
    return f(r) * h * np.cos(s) / np.sqrt(max(pr2(r), 1e-300))


T = 2 * quad(integrand, -np.pi / 2, np.pi / 2, args=(lambda r: 1.0,), epsabs=0, epsrel=1e-13, limit=400)[0]
dth = 2 * quad(integrand, -np.pi / 2, np.pi / 2, args=(lambda r: Lz / r ** 2,), epsabs=0, epsrel=1e-13, limit=400)[0]


def rhs(t, y, k):
    q, v, M = y[0:2], y[2:4], y[4:].reshape(4, 4)
    rr = q @ q
    hess = k * k * np.eye(2) + LAM * (rr * np.eye(2) + 2 * np.outer(q, q))
    dM = np.vstack((M[2:], -hess @ M[:2]))
    return np.r_[v, -LAM * rr * q, dM.ravel()]


def monodromy(k):
    y0 = np.r_[rmax, 0.0, 0.0, Lz / rmax, np.eye(4).ravel()]
    s = solve_ivp(rhs, (0, T), y0, args=(k,), method='Radau', rtol=1e-12, atol=1e-12)
    yT = s.y[:, -1]
    ang = np.arctan2(yT[1], yT[0])
    Rm = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
    close = np.kron(np.eye(2), Rm)
    Mon = close @ yT[4:].reshape(4, 4)
    return Mon, ang, np.hypot(yT[0], yT[1]) / rmax - 1, abs(yT[2] * np.cos(ang) + yT[3] * np.sin(ang))


res = dict(E=E, L=Lz, rmin=rmin, rmax=rmax, radial_period_quadrature=T, precession_quadrature=dth,
           published_period=pub['radial_period'], published_rotation_angle=pub['rotation_angle'])
_, ang, rerr, vr = monodromy(1.0)
res['orbit_closure'] = dict(angle_from_integration=ang, radius_mismatch=rerr, radial_velocity_at_T=vr)
J4 = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
rows, symp = [], 0.0
for r in pub['rows']:
    Mon, *_ = monodromy(r['k'])
    symp = max(symp, float(np.max(abs(Mon.T @ J4 @ Mon - J4))))
    m = float(abs(np.linalg.eigvals(Mon)).max())
    rows.append(dict(k=r['k'], multiplier=m, mu=np.log(m) / T, published_mu=r['mu_conformal'],
                     published_multiplier=r['multiplier_max']))
fine = []
for k in np.linspace(0.05, 6.0, 120):
    Mon, *_ = monodromy(k)
    fine.append(dict(k=float(k), mu=float(np.log(abs(np.linalg.eigvals(Mon)).max()) / T)))
res.update(rows=rows, fine_scan=fine, symplectic_error=symp,
           max_rel_mu_diff=max(abs(x['mu'] - x['published_mu']) / max(x['published_mu'], 1e-3) for x in rows),
           peak_on_published_nodes=max(rows, key=lambda x: x['mu']),
           peak_fine=max(fine, key=lambda x: x['mu']),
           band_fine=[x['k'] for x in fine if x['mu'] > 1e-6][:1] + [x['k'] for x in fine if x['mu'] > 1e-6][-1:])
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('rows', 'fine_scan')}, indent=1))
