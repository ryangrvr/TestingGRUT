#!/usr/bin/env python3
"""Robustness of the release angular instability to the supplied background (own code).

Linear-regime growth of the imaginary (angular) components from delta0=1e-14 along the light mode:
  RD, H=4.7e13 (baseline); MD (inflaton-oscillation era, no reheating before a=1000), H=4.7e13;
  RD at H=1e13 (its own inflationary point).
Conformal equation u'' - (a''/a) u + a^3 dV/dz*(u/a) = 0; RD: a=tau (tau0=1); MD: a=tau^2/4 (tau0=2), H_i=1.
Outputs results/variants.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import build, GEV, COUP
from point_full import analyse
from release import initial

OUT = Path(__file__).resolve().parents[1] / 'results'


def rhs_bg(p, bg):
    lX, lY, lS, k, M, h, mu = p['lam_X'], p['lam_Y'], p['lam_S'], p['kap'], p['M'], p['h'], p['mu']

    def rhs(t, s):
        if bg == 'RD':
            a, app = t, 0.0
        else:
            a, app = t * t / 4, 0.5
        uX, uY, uS = s[0] + 1j * s[1], s[2] + 1j * s[3], s[4] + 1j * s[5]
        fX = 2 * lX * abs(uX)**2 * uX + 3 * k * uY * np.conj(uX)**2
        fY = 2 * lY * abs(uY)**2 * uY + a * a * M * M * uY + k * uX**3 + 2 * mu * a * uS * np.conj(uY)
        fS = 2 * lS * (abs(uS)**2 - a * a * h * h / 2) * uS + mu * a * uY**2
        g = app / a
        fX, fY, fS = fX - g * uX, fY - g * uY, fS - g * uS
        return [s[6], s[7], s[8], s[9], s[10], s[11], -fX.real, -fX.imag, -fY.real, -fY.imag, -fS.real, -fS.imag]
    return rhs


def growth(Hinf, bg):
    pt = analyse(build(Hinf=Hinf), f'H={Hinf}')
    p = dict(M=GEV['M'] / Hinf, h=GEV['h'] / Hinf, mu=GEV['mu'] / Hinf, **COUP)
    d0 = 1e-14
    s0 = np.array(initial(pt, d0))
    if bg == 'MD':     # u = a z, u' = a' z + a^2 zdot ; at tau0=2: a=1, a'=1
        t0, t1 = 2.0, 2 * np.sqrt(1000.)
    else:
        t0, t1 = 1.0, 1000.
    a_of = (lambda t: t) if bg == 'RD' else (lambda t: t * t / 4)
    a_targets = [10, 20, 35, 50, 70, 100, 150, 200, 300, 500, 1000]
    t_eval = [t0 + 0] + [(a if bg == 'RD' else 2 * np.sqrt(a)) for a in a_targets]
    sol = solve_ivp(rhs_bg(p, bg), (t0, t1), s0, method='DOP853', rtol=1e-11, atol=1e-30, t_eval=t_eval)
    im = np.sqrt(sum(sol.y[i]**2 for i in [1, 3, 5, 7, 9, 11])) / d0
    return dict(H_GeV=Hinf, background=bg, light_mass_over_H=pt['light_mode']['mass_over_H'],
                Im_norm_over_delta0={str(a): float(v) for a, v in zip([1] + a_targets, im)},
                total_growth_to_a1000=float(im[-1] / im[0]), max_growth=float(im.max() / im[0]))


def main():
    rows = [growth(4.7e13, 'RD'), growth(4.7e13, 'MD'), growth(1e13, 'RD')]
    for r in rows:
        print(r['background'], r['H_GeV'], 'light m/H', round(r['light_mass_over_H'], 4), 'max growth', f"{r['max_growth']:.3e}",
              {k: f'{v:.2e}' for k, v in r['Im_norm_over_delta0'].items()})
    (OUT / 'variants.json').write_text(json.dumps(rows, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
