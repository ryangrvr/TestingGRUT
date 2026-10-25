#!/usr/bin/env python3
"""Charge budget of core-catalysed erasure (bb211d1 conversion Xbar + core(Q) -> X + core(Q-2)) (own code).

Each conversion destroys one X-Xbar pair and two bound charges and emits one fast carrier of the
object's sign.  Populations after fragmentation (NR quanta per net charge): N+- = (1/eps +- 1)/2;
a fraction beta of each sign is bound in same-sign lumps, the rest free.
Analytic end state when conversion is fast (Gamma/H >> 1 to T_vis ~ 0.6 keV, channels.json):
  minority fully erased iff  beta N+ >= (2-beta) N-  <=>  eps >= 1-beta;
  final bound fraction of net charge c_f = 1-(1-beta)/eps; otherwise remaining Xbar per net charge
  p_f >= ((1-beta)/eps - 1)/2  (all bound charge spent).
A mass-action integration (both signs, conversions only, daughters cooled instantly) checks the algebra.
Late matter per net charge = c_f E(Q)/Q + (1-c_f) m_X  (+ 2 p_f m_X), compared with the calibrated 0.98593 m_X.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from common import *

OUT = Path(__file__).resolve().parents[1] / 'results'


def analytic(eps, beta):
    Np, Nm = (1 / eps + 1) / 2, (1 / eps - 1) / 2
    if beta * Np >= (2 - beta) * Nm:
        return dict(erased=True, c_f=1 - (1 - beta) / eps, p_f=0.0)
    return dict(erased=False, c_f=0.0, p_f=((1 - beta) / eps - 1) / 2)


def mass_action(eps, beta):
    Np, Nm = (1 / eps + 1) / 2, (1 / eps - 1) / 2
    y0 = [(1 - beta) * Np, (1 - beta) * Nm, beta * Np, beta * Nm]  # F+, F-, C+, C-
    def rhs(t, y):
        Fp, Fm, Cp, Cm = np.maximum(y, 0)
        r1 = Fm * Cp * (Cp > 1e-12)  # Xbar on core
        r2 = Fp * Cm * (Cm > 1e-12)  # X on anticore
        return [r1 - r2, r2 - r1, -2 * r1, -2 * r2]
    sol = solve_ivp(rhs, (0, 1e6), y0, method='LSODA', rtol=1e-10, atol=1e-13)
    Fp, Fm, Cp, Cm = np.maximum(sol.y[:, -1], 0)
    return dict(c_f=Cp - Cm, p_f=Fm + Cm, net=Fp + Cp - Fm - Cm)


def main():
    c = MX * YQ / BDM
    rows = []
    for beta in (0.35145, 0.5, 0.8, 0.9, 0.97, 1.0):
        for eps in (0.99, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0.03):
            a = analytic(eps, beta)
            ma = mass_action(eps, beta) if beta < 1 else None
            late = a['c_f'] * E_PER_FREE_MASS + (1 - a['c_f']) * 1.0 + 2 * a['p_f']
            rows.append(dict(beta=beta, eps=eps, **a, mass_action=ma,
                             late_matter_over_calibrated=late * MX / B_PER_CHARGE,
                             f_extra=late * MX / B_PER_CHARGE - 1))
    # eps thresholds
    thr = []
    for beta in (0.35145, 0.5, 0.8, 0.9, 0.97):
        # 1% allowance: (c E/Q + 1 - c) m_X <= 1.01 B  -> c >= (1 - 1.01 B/m)/(1 - E/Q)
        cmin = max(0.0, (1 - 1.01 * B_PER_CHARGE / MX) / (1 - E_PER_FREE_MASS))
        thr.append(dict(beta=beta, eps_erasure=1 - beta, c_min_for_1pct=cmin,
                        eps_1pct=(1 - beta) / (1 - cmin)))
    out = dict(rows=rows, thresholds=thr, note='beta=0.35145: bb211d1 prepared-collapse retention (35.145%), used only as a benchmark')
    (OUT / 'catalysis.json').write_text(json.dumps(out, indent=1, default=float))
    for r in rows:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k != 'mass_action'},
              {k: round(v, 5) for k, v in (r['mass_action'] or {}).items()})
    print(json.dumps(thr, indent=1))


if __name__ == '__main__':
    main()
