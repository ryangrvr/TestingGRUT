#!/usr/bin/env python3
"""Tree-level boundedness of the full X, Y, Sigma potential (own code).

1. H=0 (units M_Y): numerical global search over the six canonical real components, log-uniform
   amplitudes 1e-3..1e4 M_Y with random phases; Young-inequality and Sigma>M_Y bound checks.
2. The exact 'integrated-out' potential V_eff(|X|) = min over Y, Sigma (aligned phases, full quartics),
   compared with the truncated quadratic-elimination EFT lam_X x^4 - (|beta6| + 2|alpha6|) x^6 and with the
   derived asymptotic quartic lam_eff = lam_X - (3/2) kap^(4/3) (2 lam_Y)^(-1/3).
3. With the supplied -c_H H^2|X|^2 at H=4.7e13 GeV: analytic floor -c_H^2 H^4/(2 lam_X) versus the found minimum.
Outputs results/boundedness.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import build, radial_V, GEV, COUP

OUT = Path(__file__).resolve().parents[1] / 'results'


def global_search_H0(n=400, seed=11):
    m = build(c_H=0.0, Hinf=GEV['M'])          # units of M_Y, no Hubble mass
    rng = np.random.default_rng(seed)
    worst = []
    for _ in range(n):
        amp = 10**rng.uniform(-3, 4, 3); ph = rng.uniform(-np.pi, np.pi, 3)
        z0 = np.ravel([[np.sqrt(2) * a * np.cos(t), np.sqrt(2) * a * np.sin(t)] for a, t in zip(amp, ph)])
        # minimise V/(1+|z|^4) so that runaway directions would show up as negative values at large |z|
        f = lambda z: m['V'](z) / (1 + np.dot(z, z)**2)
        r = minimize(f, z0, method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-10, fatol=1e-16))
        worst.append((r.fun, m['V'](r.x), np.linalg.norm(r.x) / np.sqrt(2)))
    worst.sort()
    # direct sampling of V/|z|^4 on large spheres (aligned and random phases)
    samp = []
    for Rr in [1e2, 1e3, 1e4, 1e5]:
        for _ in range(4000):
            v = rng.normal(size=6); v /= np.linalg.norm(v)
            samp.append(m['V'](Rr * v) / Rr**4)
    return dict(min_scaled=worst[0][0], V_at_min_scaled=worst[0][1], modulus_at_min=worst[0][2],
                min_V_over_r4_on_spheres=float(min(samp)), vacuum_value=m['V'](np.array([0, 0, 0, 0, m['h'], 0])))   # X=Y=0, |Sigma|=h/sqrt2 -> s1=h


def young_and_sigma():
    c = COUP
    kmax = (c['lam_X']**3 * c['lam_Y'] / 27)**0.25
    x = np.logspace(-4, 4, 2001)[:, None]; y = np.logspace(-6, 4, 2001)[None, :]
    lhs = 2 * c['kap'] * y * x**3; rhs = c['lam_X'] / 2 * x**4 + c['lam_Y'] / 2 * y**4
    s = np.linspace(1, 1e3, 200001)
    f = (s**2 - .005)**2 - .5 * (s - 1)**2 - .490025 * s**2
    return dict(kappa=c['kap'], kappa_max_sufficient=kmax, young_min_margin=float((rhs - lhs).min()),
                young_min_ratio=float((rhs / lhs).min()), sigma_bound_min=float(f.min()))


def eliminated_curve():
    m = build(c_H=0.0, Hinf=GEV['M'])
    c = COUP; M = 1.0; h = GEV['h'] / GEV['M']; mu = GEV['mu'] / GEV['M']
    D = 1 - 4 * mu**2 * h**2 / 2
    beta6 = -c['kap']**2 / D; alpha6 = mu * c['kap']**2 * (h / np.sqrt(2)) / D
    lam_eff = c['lam_X'] - 1.5 * c['kap']**(4 / 3) * (2 * c['lam_Y'])**(-1 / 3)
    rows = []
    for x in np.logspace(-2, 5, 57):
        best = None
        for y0, s0 in [(c['kap'] * x**3, h / np.sqrt(2)), (0.08 * x, h / np.sqrt(2)), (0.08 * x, (mu * (0.08 * x)**2 / 2)**(1 / 3))]:
            r = minimize(lambda q: radial_V([x, abs(q[0]), abs(q[1])], m, c_H=0.0), [y0, s0], method='Nelder-Mead',
                         options=dict(xatol=1e-13, fatol=1e-30, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
        trunc = c['lam_X'] * x**4 - (abs(beta6) + 2 * abs(alpha6)) * x**6
        rows.append(dict(x_over_MY=x, V_eff=best.fun, V_eff_over_x4=best.fun / x**4, truncated_EFT=trunc,
                         y=abs(best.x[0]), s=abs(best.x[1]), kap_x2=c['kap'] * x**2))
    xc = np.sqrt(c['lam_X'] / (abs(beta6) + 2 * abs(alpha6)))
    return dict(beta6_M2=beta6, alpha6_M2=alpha6, truncated_zero_crossing_x_over_MY=xc,
                truncated_zero_crossing_GeV=xc * GEV['M'], kap_x2_over_M2_at_crossing=c['kap'] * xc**2,
                lam_eff_derived=lam_eff, min_V_eff=min(r['V_eff'] for r in rows),
                V_eff_over_x4_at_largest=rows[-1]['V_eff_over_x4'], rows=rows)


def main():
    out = dict(young=young_and_sigma(), H0_search=global_search_H0(), eliminated=eliminated_curve())
    out['with_Hubble'] = dict(analytic_floor_over_H4=-1.0 / (2 * COUP['lam_X']))
    pt = json.loads((OUT / 'point_full.json').read_text())['baseline']
    out['with_Hubble']['found_min_over_H4'] = pt['V_min_over_H4']
    OUT.mkdir(exist_ok=True)
    (OUT / 'boundedness.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    e = out['eliminated']
    print('young', out['young']); print('H0', out['H0_search'])
    print('trunc crossing x/M', e['truncated_zero_crossing_x_over_MY'], 'kap x^2/M^2', e['kap_x2_over_M2_at_crossing'],
          'minVeff', e['min_V_eff'], 'Veff/x4 large', e['V_eff_over_x4_at_largest'], 'lam_eff', e['lam_eff_derived'])
    for r in e['rows'][::8]:
        print(f"x={r['x_over_MY']:.3e} Veff/x4={r['V_eff_over_x4']:.5e} trunc/x4={r['truncated_EFT'] / r['x_over_MY']**4:.4e} y/x={r['y'] / r['x_over_MY']:.4e}")


if __name__ == '__main__':
    main()
