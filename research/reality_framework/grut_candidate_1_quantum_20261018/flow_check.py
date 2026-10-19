#!/usr/bin/env python3
"""Bounded one-loop scalar RG checks; fixed original matching, v(mu0)=1.

The independent exact tensor/polynomial derivation is in rg_symbolic.py.
These numerical tests evaluate the actual scalar mass eigenvalues and logarithms.
They do not calculate a rotating soliton fluctuation determinant or a pole mass.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
NAMES = ['lambda_phi', 'lambda_X', 'c', 'd', 'lambda_D', 'alpha',
         'mphi2', 'mX2', 'mD2', 'C']
THETA0 = np.array([.01, .001, .01, .0001, .001, 0., -.01, 0., .01, .0025])
LOOP = 16 * np.pi**2
MU0 = .1


def beta(theta):
    lp, lx, c, d, ld, alpha, mp, mx, md, C = theta
    return np.array([
        18*lp**2+4*c**2+2*d**2,
        20*lx**2+2*c**2+2*alpha**2,
        c*(6*lp+8*lx+8*c)+2*d*alpha,
        d*(6*lp+6*ld+8*d)+4*c*alpha,
        18*ld**2+2*d**2+4*alpha**2,
        2*c*d+alpha*(8*lx+6*ld+8*alpha),
        6*lp*mp+4*c*mx+2*d*md,
        2*c*mp+8*lx*mx+2*alpha*md,
        2*d*mp+4*alpha*mx+6*ld*md,
        (mp**2+2*mx**2+md**2)/2,
    ])/LOOP


def parameter_gradient(q):
    phi, x, y, D = np.asarray(q)
    R = x*x+y*y
    return np.array([phi**4/4, R**2/4, phi**2*R/2,
                     phi**2*D**2/2, D**4/4, R*D**2/2,
                     phi**2/2, R/2, D**2/2, 1.])


def tree(q, theta):
    return parameter_gradient(q) @ theta


def mass_matrix(q, theta):
    """Direct differentiated polynomial; no effective_potential.py imports."""
    lp, lx, c, d, ld, alpha, mp, mx, md, C = theta
    phi, x, y, D = q
    R = x*x+y*y
    mat = np.diag([
        mp+3*lp*phi**2+c*R+d*D**2,
        mx+lx*(3*x*x+y*y)+c*phi**2+alpha*D**2,
        mx+lx*(x*x+3*y*y)+c*phi**2+alpha*D**2,
        md+3*ld*D**2+d*phi**2+alpha*R,
    ])
    mat[0, 1] = mat[1, 0] = 2*c*phi*x
    mat[0, 2] = mat[2, 0] = 2*c*phi*y
    mat[0, 3] = mat[3, 0] = 2*d*phi*D
    mat[1, 2] = mat[2, 1] = 2*lx*x*y
    mat[1, 3] = mat[3, 1] = 2*alpha*x*D
    mat[2, 3] = mat[3, 2] = 2*alpha*y*D
    return mat


def cw(q, theta, mu):
    mass2 = np.linalg.eigvalsh(mass_matrix(q, theta))
    if np.any(mass2 <= 0):
        raise ValueError(f'Nonpositive constant-background Hessian: {mass2}')
    return float(np.sum(mass2**2*(np.log(mass2/mu**2)-1.5))/(64*np.pi**2))


def integrate(ratio, rtol, atol):
    t = np.log(ratio)
    if t == 0:
        return THETA0.copy(), 0
    sol = solve_ivp(lambda t, z: beta(z), (0, t), THETA0,
                    method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1], sol.nfev


def main():
    start = time.perf_counter()
    protocol_bytes = (HERE/'flow_protocol.json').read_bytes()
    protocol = json.loads(protocol_bytes)
    assert protocol['parameters_in_order'] == NAMES
    assert np.array_equal(protocol['matching_values'], THETA0)
    ratios = protocol['scale_ratios']
    points = np.array(protocol['backgrounds_phi_x_y_D'])
    b0 = beta(THETA0)
    trajectories = []
    solver_differences = []
    records = []
    worst_identity = 0.
    min_mass = np.inf
    for ratio in ratios:
        th, nfev = integrate(ratio, 1e-11, 1e-15)
        thfine, _ = integrate(ratio, 1e-12, 1e-16)
        solver_differences.append(float(np.max(np.abs(th-thfine))))
        lp, lx, c, d, ld, alpha, mp, mx, md, C = th
        bounded = bool(min(lp, lx, ld) > 0 and min(c, d) >= 0
                       and alpha > -np.sqrt(lx*ld))
        trajectories.append({
            'mu_over_mu0': ratio, 'mu_over_v': MU0*ratio,
            'parameters': dict(zip(NAMES, map(float, th))),
            'leading_log_parameters': dict(zip(NAMES, map(float, THETA0+np.log(ratio)*b0))),
            'ode_function_evaluations': nfev,
            'quartic_sufficient_boundedness': bounded,
            'X_D_copositivity_margin': float(alpha+np.sqrt(lx*ld)),
            'tree_vacuum_point_phi': float(np.sqrt(-mp/lp)),
            'tree_vacuum_value': float(C-mp*mp/(4*lp)),
            'binding_coupling_ratio_c_squared_over_lambda_phi_lambda_X': float(c*c/(lp*lx)),
        })
        for i, q in enumerate(points):
            masses = np.linalg.eigvalsh(mass_matrix(q, th))
            min_mass = min(min_mass, masses.min())
            lhs = beta(th) @ parameter_gradient(q)
            rhs = np.sum(masses**2)/(32*np.pi**2)
            residual = abs(lhs-rhs)/max(abs(lhs), abs(rhs), 1e-30)
            worst_identity = max(worst_identity, residual)
            vv0 = tree(q, th)
            vv1 = cw(q, th, MU0*ratio)
            strict = tree(q, THETA0)+np.log(ratio)*(b0 @ parameter_gradient(q))+cw(q, THETA0, MU0*ratio)
            records.append({
                'background_index': i, 'mu_over_mu0': ratio,
                'V0_over_v4': float(vv0), 'CW_over_v4': vv1,
                'V0_plus_CW_over_v4': float(vv0+vv1),
                'strict_one_loop_over_v4': float(strict),
                'minimum_tree_mass_squared_over_v2': float(masses.min()),
                'RG_identity_relative_residual': float(residual),
            })
    scale_variations = []
    strict_worst = 0.
    for i, q in enumerate(points):
        rs = [r for r in records if r['background_index'] == i]
        v0spread = float(np.ptp([r['V0_over_v4'] for r in rs]))
        fullspread = float(np.ptp([r['V0_plus_CW_over_v4'] for r in rs]))
        strictspread = float(np.ptp([r['strict_one_loop_over_v4'] for r in rs]))
        strict_worst = max(strict_worst, strictspread)
        scale_variations.append({
            'background_index': i, 'point': q.tolist(),
            'tree_spread_over_v4': v0spread,
            'one_loop_running_spread_over_v4': fullspread,
            'strict_one_loop_spread_over_v4': strictspread,
            'reduction_tree_spread_divided_by_one_loop_spread': v0spread/fullspread,
        })
    bookkeeping = []
    worst_adjacent = 0.
    for i, q in enumerate(points):
        for ratio in [.25, .5, 2., 4.]:
            values = []
            for eta in [1., .5, .25, .125]:
                theta = THETA0+eta*np.log(ratio)*b0
                value = tree(q, theta)+eta*cw(q, theta, MU0*ratio)
                reference = tree(q, THETA0)+eta*cw(q, THETA0, MU0)
                residual = float(value-reference)
                values.append({'eta': eta, 'residual_over_v4': residual,
                               'residual_divided_by_eta_squared': residual/eta**2})
            normalized = np.array([v['residual_divided_by_eta_squared'] for v in values])
            adjacent = float(np.max(np.abs(normalized[:-1]/normalized[1:]-1)))
            worst_adjacent = max(worst_adjacent, adjacent)
            bookkeeping.append({'background_index': i, 'mu_over_mu0': ratio,
                                'values': values, 'worst_adjacent_relative_difference': adjacent})
    checks = {
        'positive_Hessians_on_prescribed_backgrounds': bool(min_mass > 0),
        'sufficient_quartic_boundedness_all_sampled_scales': all(t['quartic_sufficient_boundedness'] for t in trajectories),
        'numerical_RG_identity': bool(worst_identity < 1e-11),
        'strict_order_scale_cancellation': bool(strict_worst < 1e-16),
        'residual_scales_as_loop_order_squared': bool(worst_adjacent < .03),
        'ODE_tolerance_reproduction': bool(max(solver_differences) < 1e-12),
    }
    result = {
        'protocol_sha256': hashlib.sha256(protocol_bytes).hexdigest(),
        'command': 'python flow_check.py',
        'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform()},
        'beta_at_matching': dict(zip(NAMES, map(float, b0))),
        'trajectories': trajectories, 'scale_variations': scale_variations,
        'bookkeeping': bookkeeping,
        'maximum_RG_identity_relative_residual': float(worst_identity),
        'maximum_strict_order_scale_spread_over_v4': float(strict_worst),
        'maximum_bookkeeping_adjacent_relative_change': float(worst_adjacent),
        'maximum_ODE_tolerance_absolute_difference': max(solver_differences),
        'minimum_tree_mass_squared_on_flow_grid': float(min_mass),
        'checks': checks, 'all_checks_pass': all(checks.values()),
        'interpretation': 'Same matching data at every scale. Numerical RG resummation and V0+V1 leave higher-order residuals; these are neither a two-loop result nor a rigorous theory-error bound.',
        'elapsed_seconds': time.perf_counter()-start,
    }
    (HERE/'flow_results.json').write_text(json.dumps(result, indent=2)+'\n')
    with (HERE/'flow_scale_values.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (HERE/'flow_execution.jsonl').write_text(json.dumps({
        'checks': checks, 'elapsed_seconds': result['elapsed_seconds'],
        'protocol_sha256': result['protocol_sha256'],
        'maximum_bookkeeping_adjacent_relative_change': worst_adjacent,
    })+'\n')
    print(json.dumps({k: result[k] for k in ['checks', 'all_checks_pass', 'maximum_RG_identity_relative_residual', 'maximum_strict_order_scale_spread_over_v4', 'maximum_bookkeeping_adjacent_relative_change', 'maximum_ODE_tolerance_absolute_difference', 'elapsed_seconds']}, indent=2))
    if not result['all_checks_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
