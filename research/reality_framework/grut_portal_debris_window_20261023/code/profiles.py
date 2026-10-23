#!/usr/bin/env python3
"""Classical stationary test profiles; their charges are not formation outputs."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid, simpson

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent/'grut_unified_working_theory_20261022'
G, LX = .1, .001
LAMBDAS = [.005, .0075, .01, .015, .019, .02, .02005, .0202, .021, .025, .03, .04]


class Seed:
    def __init__(self, q=79988):
        a = np.loadtxt(PARENT/'data'/f'profile_Q{q}.csv', delimiter=',', skiprows=1)
        self.x = a[:, 0]
        self.y = a[:, 1:].T
        self.j = cumulative_trapezoid(self.x**2*self.y[0]**2, self.x, initial=0.)
        d = json.loads((PARENT/'RESULTS.json').read_text())
        self.p = np.array([next(r for r in d['profiles'] if r['Q'] == q)['Omega']])

    def sol(self, r):
        a = [np.interp(r, self.x, row, left=row[0], right=row[-1]) for row in self.y]
        a.append(np.interp(r, self.x, self.j, left=0., right=self.j[-1]))
        return np.array(a)


def solve(q, lam, seed, length=80., tolerance=2e-7):
    def rhs(r, y, p):
        f, fp, s, sp, j = y
        return np.array([fp, (s*s+LX/G**2*f*f-p[0]**2)*f,
                         sp, (lam/G**2*(s*s-1)+f*f)*s, r*r*f*f])

    def bc(l, r, p):
        return [l[1], l[3], r[0], r[2]-1, l[4], 4*np.pi*p[0]*r[4]-G**2*q]

    r = np.linspace(0, length, 901)
    guess = seed.sol(np.minimum(r, seed.x[-1]))
    outside = r > seed.x[-1]
    guess[:, outside] = np.array([0, 0, 1, 0, guess[4, ~outside][-1]])[:, None]
    sol = solve_bvp(rhs, bc, r, guess, p=seed.p, S=np.diag([0., -2., 0., -2., 0.]),
                    tol=tolerance, bc_tol=tolerance, max_nodes=30000)
    if not sol.success:
        raise RuntimeError(sol.message)
    test = sol.sol(np.linspace(0, length, 4001))
    if not (0 < sol.p[0] < 1 and np.min(test[0]) >= -1e-7
            and np.min(test[2]) >= -1e-7 and np.max(test[2]) <= 1.000001):
        raise RuntimeError('Non-admissible branch')
    return sol


def quantities(sol, q, lam):
    r = np.linspace(0, sol.x[-1], 16001)
    f, fp, s, sp, _ = sol.sol(r)
    omega = float(sol.p[0])
    grad = .5*(fp**2+sp**2)
    kinetic = .5*omega**2*f**2
    u = lam/G**2/4*(s*s-1)**2+.5*s*s*f*f+LX/G**2/4*f**4
    integral = lambda y: float(4*np.pi*simpson(r*r*y, x=r))
    eb = integral(grad+kinetic+u)
    depletion = integral(1-s*s)
    return dict(Q=q, lambda_phi=lam, Omega=omega, energy_per_v_GeV=eb/G,
                energy_per_free_mass=eb/(G**2*q), F0=float(f[0]), S0=float(s[0]),
                charge_relative_residual=omega*integral(f*f)/(G**2*q)-1,
                virial_relative_residual=(integral(grad)+3*(integral(u)-integral(kinetic)))/eb,
                mass_depletion_identity_relative_residual=(eb-lam/G**2*depletion)/eb,
                max_BVP_residual=float(max(sol.rms_residuals)),
                radial_domain=float(sol.x[-1]), tolerance=2e-7)


def save(sol, q, lam, refined=False):
    tag = f'Q{int(q)}_lambda{lam:g}'+('_refined' if refined else '')
    r = np.linspace(0, sol.x[-1], 1601)
    a = sol.sol(r)
    np.savetxt(ROOT/'data'/f'profile_{tag}.csv', np.column_stack([r, a[:4].T]),
               delimiter=',', header='rho,F,Fprime,S,Sprime', comments='')
    row = quantities(sol, q, lam)
    row['profile_file'] = f'data/profile_{tag}.csv'
    return row


def main():
    (ROOT/'data').mkdir(exist_ok=True, parents=True)
    rows, attempts = [], []
    seeds = {}
    seed = Seed()
    for q in np.geomspace(79988, 1e6, 19)[1:]:
        seed = solve(float(q), .01, seed)
    seeds[.01] = seed
    row = save(seed, 1000000, .01)
    rows.append(row)
    print('profile', row, flush=True)
    for direction in ([.0075, .005], [.015, .019, .02, .02005, .0202, .021, .025, .03, .04]):
        last_lam, seed = .01, seeds[.01]
        for target in direction:
            for lam in np.linspace(last_lam, target, max(2, int(abs(target-last_lam)/.001)+1))[1:]:
                seed = solve(1000000, float(lam), seed)
            last_lam = target
            seeds[target] = seed
            row = save(seed, 1000000, target)
            rows.append(row)
            print('profile', row, flush=True)
    # A continuation attempt at the old test charge; failure is not a
    # mathematical nonexistence result and never leads to a forced profile.
    original = Seed(39944)
    original = solve(39944, .01, original)
    rows.append(save(original, 39944, .01))
    last = .01
    for target in [.015, .019, .02, .021, .025, .03, .04]:
        try:
            for lam in np.linspace(last, target, max(2, int((target-last)/.0005)+1))[1:]:
                original = solve(39944, float(lam), original)
            last = target
            rows.append(save(original, 39944, target))
            attempts.append(dict(Q=39944, lambda_phi=target, converged=True))
        except RuntimeError as e:
            attempts.append(dict(Q=39944, lambda_phi=target, converged=False, message=str(e)))
            print('old-charge continuation stops', attempts[-1], flush=True)
            break
    refined = solve(1000000, .03, seeds[.03], length=100., tolerance=2e-9)
    row = save(refined, 1000000, .03, True)
    row['tolerance'] = 2e-9
    result = dict(scope='Stationary test profiles; not a formed charge distribution or a quantum lifetime.',
                  profiles=rows, old_charge_attempts=attempts, refinement=row)
    (ROOT/'PROFILE_RESULTS.json').write_text(json.dumps(result, indent=2)+'\n')
    print('completed', len(rows), 'profiles plus refinement', flush=True)


if __name__ == '__main__':
    main()
