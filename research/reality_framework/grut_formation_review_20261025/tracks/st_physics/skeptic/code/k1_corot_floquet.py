"""st_physics skeptic, check k1 (own code, independent method).

Floquet spectrum of the isolated complex quartic orbit written in the CO-ROTATING radial/tangential frame:
    chi = r e^{i theta},  delta chi = (u + i w) e^{i theta},  theta' = L / r^2,
    u'' = 2 theta' w' + theta'' w + (theta'^2 - k^2 - 3 lam r^2) u
    w'' = -2 theta' u' - theta'' u + (theta'^2 - k^2 - lam r^2) w
    r''  = L^2 / r^3 - lam r^3.
The coefficients depend only on r(eta), so they are exactly periodic with the radial period T; the monodromy is taken
from pericentre to pericentre (event r' = 0, ascending) with no precession bookkeeping.  This is a different
formulation from the reviewer's lab-frame monodromy (apoapsis start, rotation undone afterwards).

The invariants (E, L) are recomputed from the published backgrounds with our own conversion and, separately, taken
from Codex's QUARTIC_ORBIT.json initial_chi.  Late growth is a property of the matched orbit only.

usage: python3 -I k1_corot_floquet.py <codex timing packet dir> <QUARTIC_ORBIT.json> <out.json>
"""
import sys, json, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq

LAM = 1e-3


def invariants_from_bg(bg):
    a = np.exp(bg[:, 15]); H = bg[:, 16]
    x = bg[:, 2] + 1j * bg[:, 3]; xd = bg[:, 8] + 1j * bg[:, 9]
    chi = np.sqrt(2) * a * x                       # conformal canonical field
    chip = np.sqrt(2) * a * a * (xd + H * x)       # d chi / d eta = a d(chi)/d tau, a_tau = a H
    E = 0.5 * abs(chip) ** 2 + LAM / 4 * abs(chi) ** 4
    L = np.imag(np.conj(chi) * chip)
    return a, chi, chip, E, L


class CoRot:
    def __init__(self, E, L, lam=LAM):
        self.E, self.L, self.lam = float(E), float(L), lam
        f = lambda r: 2 * E - L * L / r ** 2 - lam * r ** 4 / 2
        rc = (L * L / lam) ** (1 / 6)
        self.rmin = brentq(f, 1e-6 * rc, rc, xtol=1e-14, rtol=1e-15)
        hi = rc
        while f(hi) > 0:
            hi *= 1.5
        self.rmax = brentq(f, rc, hi, xtol=1e-14, rtol=1e-15)
        # period by integration pericentre -> pericentre
        ev = lambda t, y: y[1]
        ev.direction = 1.0
        s = solve_ivp(self._bg, (0, 50), [self.rmin, 0.0, 0.0], method='DOP853', rtol=1e-13, atol=1e-12,
                      events=ev, dense_output=True)
        # first ascending zero of r' after leaving pericentre is the next pericentre? r' goes + then -, so the
        # ascending crossing (from - to +) is at the next pericentre.
        tz = [t for t in s.t_events[0] if t > 1e-6]
        self.T = float(tz[0])
        self.dtheta = float(s.sol(self.T)[2])

    def _bg(self, t, y):
        r, rp, th = y
        return [rp, self.L ** 2 / r ** 3 - self.lam * r ** 3, self.L / r ** 2]

    def _rhs(self, t, y, k):
        r, rp = y[0], y[1]
        tp = self.L / r ** 2
        tpp = -2 * self.L * rp / r ** 3
        Y = y[2:].reshape(4, 4)        # rows: u, w, u', w'
        u, w, up, wp = Y
        upp = 2 * tp * wp + tpp * w + (tp * tp - k * k - 3 * self.lam * r * r) * u
        wpp = -2 * tp * up - tpp * u + (tp * tp - k * k - self.lam * r * r) * w
        return np.r_[rp, self.L ** 2 / r ** 3 - self.lam * r ** 3, np.vstack([up, wp, upp, wpp]).ravel()]

    def monodromy(self, k, rtol=1e-12, atol=1e-13):
        y0 = np.r_[self.rmin, 0.0, np.eye(4).ravel()]
        s = solve_ivp(self._rhs, (0, self.T), y0, args=(k,), method='DOP853', rtol=rtol, atol=atol)
        return s.y[2:, -1].reshape(4, 4), s.y[:2, -1]

    def mu(self, k):
        Mo, end = self.monodromy(k)
        ev = np.linalg.eigvals(Mo)
        m = float(abs(ev).max())
        return np.log(m) / self.T, Mo, ev


def scan(orb, KS):
    return [dict(k=float(k), mu=orb.mu(k)[0]) for k in KS]


def main():
    pkt, qorb, out = sys.argv[1:4]
    t0 = time.time()
    res = dict(method='own co-rotating-frame Floquet (independent formulation); invariants recomputed by own code')
    orbits = {}
    for tag in ('quiet3', 'quiet5'):
        d = np.load(f'{pkt}/{tag}.npz')
        times, bg = d['times'], d['bg']
        a, chi, chip, E, L = invariants_from_bg(bg)
        i10 = int(np.argmin(abs(times - 10)))
        late = times >= 20
        res[tag + '_invariants'] = dict(E10=float(E[i10]), L10=float(L[i10]), E_late=float(E[late].mean()),
                                        L_late=float(L[late].mean()),
                                        E_late_ptp_rel=float(np.ptp(E[late]) / E[late].mean()),
                                        L_late_ptp_rel=float(np.ptp(L[late]) / L[late].mean()),
                                        L_64_80_mean=float(L[(times >= 64) & (times <= 80)].mean()),
                                        L_800=float(L[-1]),
                                        L_200_800_ptp_rel=float(np.ptp(L[times >= 200]) / L[times >= 200].mean()),
                                        E_200_800_ptp_rel=float(np.ptp(E[times >= 200]) / E[times >= 200].mean()))
        orbits[tag + '_tau10'] = (E[i10], L[i10])
        orbits[tag + '_late'] = (E[late].mean(), L[late].mean())
    q = json.load(open(qorb))
    c0 = np.array(q['initial_chi'])
    Eq = 0.5 * (c0[2] ** 2 + c0[3] ** 2) + LAM / 4 * (c0[0] ** 2 + c0[1] ** 2) ** 2
    Lq = c0[0] * c0[3] - c0[1] * c0[2]
    orbits['codex_aux'] = (Eq, Lq)
    KS = np.array([.02, .05, .1, .2, .4, .8, 1.2, 1.6, 2., 2.13, 2.4, 2.8, 3.2, 4., 5.6, 5.8, 6., 8.])
    for name, (E0, L0) in orbits.items():
        orb = CoRot(E0, L0)
        row = dict(E=float(E0), L=float(L0), rmin=orb.rmin, rmax=orb.rmax, T_eta=orb.T, precession=orb.dtheta)
        row['KS'] = scan(orb, KS)
        grid = np.arange(0.05, 7.0001, 0.05)
        mus = np.array([orb.mu(k)[0] for k in grid])
        j = int(np.argmax(mus))
        pk = minimize_scalar(lambda k: -orb.mu(k)[0], bounds=(grid[j] - 0.06, grid[j] + 0.06), method='bounded',
                             options=dict(xatol=1e-5))
        row['peak'] = dict(k=float(pk.x), mu=float(-pk.fun))
        thr = 1e-6
        on = mus > thr
        edges = []
        for i in range(len(grid) - 1):
            if on[i] != on[i + 1]:
                g = lambda k: orb.mu(k)[0] - thr
                edges.append(float(brentq(g, grid[i], grid[i + 1], xtol=1e-6)))
        row['band_edges_mu_gt_1e-6'] = edges
        if len(edges) >= 3:
            sel = (grid > edges[1]) & (grid < (edges[2] if len(edges) > 2 else 7))
            jj = int(np.argmax(np.where(sel, mus, -1)))
            pk2 = minimize_scalar(lambda k: -orb.mu(k)[0], bounds=(grid[jj] - 0.06, grid[jj] + 0.06),
                                  method='bounded', options=dict(xatol=1e-5))
            row['second_band_peak'] = dict(k=float(pk2.x), mu=float(-pk2.fun))
        row['ir_mu_over_k'] = [dict(k=k, mu_over_k=orb.mu(k)[0] / k) for k in (0.002, 0.005, 0.01, 0.02, 0.05)]
        mu0, M0, ev0 = orb.mu(0.0)
        # Jordan structure at k = 0: rank of (M - I) and the size of its nilpotent part
        D = M0 - np.eye(4)
        sv = np.linalg.svd(D, compute_uv=False)
        row['k0'] = dict(mu=mu0, eig_abs=np.abs(ev0).tolist(), singular_values_M_minus_I=sv.tolist(),
                         M_minus_I_squared_norm=float(np.linalg.norm(D @ D)), M_minus_I_norm=float(np.linalg.norm(D)))
        res[name] = row
        print(name, {k: row[k] for k in ('E', 'L', 'rmin', 'rmax', 'T_eta', 'peak', 'band_edges_mu_gt_1e-6')},
              row.get('second_band_peak'), row['ir_mu_over_k'][2], 'k0 sv', [f'{x:.3g}' for x in sv],
              f'{time.time() - t0:.0f}s', flush=True)
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
