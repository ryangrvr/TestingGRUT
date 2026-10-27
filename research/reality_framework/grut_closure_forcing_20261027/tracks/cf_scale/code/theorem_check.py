#!/usr/bin/env python3
"""cf_scale: numerical check of the one-loop monotonicity statement [D] used in NOTES.md (own code).

Statement: in the fermion-free GRUT hidden sector (scalars + one abelian gauge field), below M_P,
  d lambda_min / d ln mu = beta_V(n*) / 16pi^2 >= 0  whenever  lambda_min = min_{|n|=1} V4(n) <= 0,
because at the constrained minimiser n*, sum_i q_i^2 n_i dV4/dn_i = 4 lambda_min sum_i q_i^2 n_i^2, so
  16pi^2 beta_V(n*) = (1/2)Tr[(V4'')^2] + (3/2) M_A(n*)^4 - 12 e^2 lambda_min sum_i q_i^2 n_i*^2 >= 0.
Check: random coupling sets with lambda_min < 0; compare the finite-difference derivative of lambda_min along the
one-loop flow with zero. Usage: python3 -B theorem_check.py <out.json>
"""
import json, sys
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rg_common as R


def main(out):
    rng = np.random.default_rng(5)
    rows = []
    n_neg = 0
    while n_neg < 40:
        lam = rng.normal(size=R.NC) * 10**rng.uniform(-3, -0.5, size=R.NC)
        e = 10**rng.uniform(-3, -0.4)
        lmin, _ = R.V4_unit_min(lam, n=24, seed=int(rng.integers(1e6)))
        if lmin > 0:
            continue
        n_neg += 1
        y = np.append(lam, e)
        dt = 1e-4
        b = R.beta(y)
        lp, _ = R.V4_unit_min(y[:R.NC] + dt * b[:R.NC], n=24, seed=1)
        lm, _ = R.V4_unit_min(y[:R.NC] - dt * b[:R.NC], n=24, seed=1)
        d = (lp - lm) / (2 * dt)
        rows.append(dict(lambda_min=lmin, e=e, dlambda_min_dt=d))
    worst = min(r['dlambda_min_dt'] for r in rows)
    print('cases with lambda_min <= 0:', len(rows), ' min d(lambda_min)/dt =', worst,
          ' (all >= 0 up to minimiser noise?)', all(r['dlambda_min_dt'] > -1e-9 for r in rows))
    Path(out).write_text(json.dumps(dict(method=__doc__, n=len(rows), min_derivative=worst, rows=rows), indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
