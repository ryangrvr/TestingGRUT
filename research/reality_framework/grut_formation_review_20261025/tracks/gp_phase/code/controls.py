#!/usr/bin/env python3
"""gp_phase (c) controls for the chaotic full-model ensembles:
  1. late-time convergence: tau_end = 1e4 versus 1e3 (a = tau);
  2. integrator tolerance: rtol 1e-12 versus 1e-10;
  3. seed amplitude sigma = 1e-6, 1e-2, 0.1 versus 1e-3 (baseline model, no extra operator);
  4. DRT-type comparator: a Hubble-induced A-term present only during inflation (e.g. the
     curvature coupling xi R Y* X^3 / Lambda^2 with R = 12 H^2, zero in radiation) fixes the
     initial minimum, then is absent after release (Dine-Randall-Thomas mechanism, covariant form).
     Implemented with the op-4 shape (Y* X^3 |X|^2) during inflation only, its strength r set
     relative to the kappa term at the inflationary point; r = 12 xi H^2/(Lambda^2 kappa) for the
     curvature form (4.5e-6 at xi = 1, Lambda = M_red).
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import gauge_point as g
from phase_scan import one, Z_ALIGNED, EPS_REQ

OUT = Path(__file__).resolve().parents[1] / 'results'


def ensemble(z, op, C, N, sigma, tau_end, rtol, atol, seed, pool):
    rng = np.random.default_rng(seed)
    jobs = [(z * np.exp(1j * sigma * rng.standard_normal(3)), op, C, tau_end, rtol, atol) for _ in range(N)]
    res = [x for x in pool.map(one, jobs) if x['ok']]
    e = np.array([x['eps'] for x in res]); nX = np.array([x['nX'] for x in res]); nm = np.array([x['nmax'] for x in res])
    return dict(N=N, n_ok=len(res), sigma=sigma, tau_end=tau_end, rtol=rtol, eps_mean=float(e.mean()),
                eps_sem=float(e.std(ddof=1) / np.sqrt(len(e))), eps_std=float(e.std(ddof=1)), eps_abs_mean=float(abs(e).mean()),
                eps_abs_max=float(abs(e).max()), eps_net=float(nX.mean() / nm.mean()), frac_meeting_requirement=float(np.mean(e >= EPS_REQ)),
                max_gauss=float(max(x['gauss'] for x in res)), max_ledger=float(max(x['ledger'] for x in res)))


def main():
    t0 = time.time(); out = dict(eps_requirement=EPS_REQ)
    op2 = 2; r2 = g.ratio_for_coefficient(2, g.coefficient_from_wilson(2, 1.0), Z_ALIGNED)
    C2 = g.SIGN[2] * g.coefficient_for_ratio(2, r2, Z_ALIGNED) * np.exp(1j * np.pi / 2)
    z2 = g.inflation_minimum(2, C2)[0]
    with Pool(4) as pool:
        out['late_time'] = dict(
            baseline_1e3=ensemble(Z_ALIGNED, 0, 0j, 32, 1e-3, 1e3, 1e-10, 1e-12, 7, pool),
            baseline_1e4=ensemble(Z_ALIGNED, 0, 0j, 32, 1e-3, 1e4, 1e-10, 1e-12, 7, pool),
            op2_natural_1e3=ensemble(z2, op2, C2, 32, 1e-3, 1e3, 1e-10, 1e-12, 7, pool),
            op2_natural_1e4=ensemble(z2, op2, C2, 32, 1e-3, 1e4, 1e-10, 1e-12, 7, pool))
        print('late', {k: (round(v['eps_mean'], 4), round(v['eps_std'], 4)) for k, v in out['late_time'].items()}, time.time() - t0, flush=True)
        out['tolerance'] = {'op2_rtol1e-10': ensemble(z2, op2, C2, 64, 1e-3, 1e3, 1e-10, 1e-12, 11, pool),
                            'op2_rtol1e-12': ensemble(z2, op2, C2, 64, 1e-3, 1e3, 1e-12, 1e-14, 11, pool)}
        print('tol', {k: (round(v['eps_mean'], 4), round(v['eps_std'], 4)) for k, v in out['tolerance'].items()}, flush=True)
        out['seed_amplitude'] = {repr(s): ensemble(Z_ALIGNED, 0, 0j, 64, s, 1e3, 1e-10, 1e-12, 13, pool) for s in (1e-6, 1e-3, 1e-2, 0.1)}
        print('sigma', {k: (round(v['eps_mean'], 4), round(v['eps_std'], 4)) for k, v in out['seed_amplitude'].items()}, flush=True)
        drt = []
        r_xi1 = 12.0 * (g.HINF / g.MRED) ** 2 / g.P['kap']  # xi R Y*X^3/M_red^2 with xi=1, R=12H^2, relative to kappa
        for r in (r_xi1, 0.1, 0.3):
            for d in (np.pi / 2, -np.pi / 2):
                C = g.SIGN[4] * g.coefficient_for_ratio(4, r, Z_ALIGNED) * np.exp(1j * d)
                z = g.inflation_minimum(4, C)[0]
                rec = ensemble(z, 0, 0j, 64, 1e-3, 1e3, 1e-10, 1e-12, 17, pool)
                det = one((z, 0, 0j, 1e3, 1e-10, 1e-12))
                rec.update(r=r, delta=d, A_initial=float(np.angle(z[0] ** 3 * np.conj(z[1]))), deterministic_eps=det.get('eps'))
                drt.append(rec)
                print('DRT', '%.3g' % r, '%.3f' % d, round(rec['eps_mean'], 4), round(rec['eps_sem'], 4), round(rec['eps_std'], 4), rec['deterministic_eps'], flush=True)
        out['DRT_hubble_A_term'] = drt
        out['DRT_r_for_xi1_Mred'] = r_xi1
    out['runtime_s'] = time.time() - t0
    (OUT / 'controls.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
