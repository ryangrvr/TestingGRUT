#!/usr/bin/env python3
"""gp_phase (c): full homogeneous X, Y, Sigma release at the 7737487 gauge point with one extra
complex gauge-invariant operator and an irremovable CP phase delta. The field starts at rest at its
OWN inflationary minimum (Hubble mass -c_H H^2|X|^2, c_H=1, H_inf=4.7e13 GeV, extra operator
included), so there is no assumed initial phase mismatch. Release: Hubble mass off, radiation era
(the same supplied release as Codex's full_release.py and our affleck_dine W1 case).

The homogeneous dynamics are chaotic (see NOTES.md), so each configuration is run as an ensemble
with small random initial phase perturbations (sigma rad on each field phase: a stand-in for
inflationary phase fluctuations, which are larger). Reported: deterministic run, ensemble mean,
standard deviation, max |eps|, net eps = <n_X>/<n_max>.

eps_X = a^3 n_X / (2 sqrt(2 lambda_X) (a^4 rho_X / 3 lambda_X)^(3/4)), the charge fraction of the X
condensate with V = lambda_X|X|^4 (same definition as tracks/affleck_dine). Requirement quoted by
Codex: eps >= 0.99175 for a 1% leftover allowance.
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import gauge_point as g

OUT = Path(__file__).resolve().parents[1] / 'results'
Z_ALIGNED = np.array([23.824169767032206 + 0j, -1.891386373567112 + 0j, -0.26697595666823387 + 0j])
EPS_REQ = (1 + 0.01 / 0.9) ** (-0.75)


def one(args):
    z, op, C, tau_end, rtol, atol = args
    rows, end = g.run_release(z, op, C, tau_end=tau_end, rtol=rtol, atol=atol, nsamp=12)
    if not rows or 'runaway_at' in (end or {}) or 'failed_at' in (end or {}):
        return dict(ok=False, end=end)
    r = rows[-1]
    nmax = r['a3nX'] / r['epsX'] if r['epsX'] != 0 else float('nan')
    return dict(ok=True, eps=r['epsX'], nX=r['a3nX'], nY=r['a3nY'], nS=r['a3nS'], nmax=nmax,
                gauss=end['max_gauss'], ledger=end['max_ledger'], eps_tau316=[x['epsX'] for x in rows if abs(x['tau'] - 316.2) < 1][0] if tau_end == 1000 else None)


def configuration(op, r, delta, N=128, sigma=1e-3, tau_end=1000.0, rtol=1e-10, atol=1e-12, seed=1, pool=None):
    if op:
        Cabs = g.coefficient_for_ratio(op, r, Z_ALIGNED)
        C = g.SIGN[op] * Cabs * np.exp(1j * delta)
    else:
        Cabs, C = 0.0, 0j
    z, v, gn, nd = g.inflation_minimum(op, C) if op else (Z_ALIGNED.copy(), None, 0.0, 1)
    rec = dict(op=op, operator=g.NAME[op], r=r, delta=delta, sigma=sigma, N=N, tau_end=tau_end, rtol=rtol)
    if z is None:
        rec['status'] = 'no local inflationary minimum near the renormalizable one'
        return rec
    if op:
        rec.update(C_over_Hunits=[C.real, C.imag], wilson=g.wilson(op, Cabs), V_min=v, grad_norm=gn,
                   min_hessian=float(g.hessian_eigs(z, op, C).min()))
    rec['z_min'] = [[float(x.real), float(x.imag)] for x in z]
    rec['angles_min'] = dict(A=float(np.angle(z[0] ** 3 * np.conj(z[1]))), B=float(np.angle(z[1] ** 2 * np.conj(z[2]))),
                             X6Sstar_phase=float(np.angle(z[0] ** 6 * np.conj(z[2]))))
    rng = np.random.default_rng(seed)
    jobs = [(z, op, C, tau_end, rtol, atol)]
    for _ in range(N):
        jobs.append((z * np.exp(1j * sigma * rng.standard_normal(3)), op, C, tau_end, rtol, atol))
    res = pool.map(one, jobs) if pool else list(map(one, jobs))
    det, ens = res[0], [x for x in res[1:] if x['ok']]
    rec['deterministic'] = det
    rec['n_ok'] = len(ens); rec['n_runaway_or_fail'] = N - len(ens)
    if ens:
        e = np.array([x['eps'] for x in ens]); nX = np.array([x['nX'] for x in ens]); nm = np.array([x['nmax'] for x in ens])
        rec.update(eps_mean=float(e.mean()), eps_sem=float(e.std(ddof=1) / np.sqrt(len(e))), eps_std=float(e.std(ddof=1)),
                   eps_abs_max=float(abs(e).max()), eps_abs_mean=float(abs(e).mean()), eps_net=float(nX.mean() / nm.mean()),
                   frac_meeting_requirement=float(np.mean(e >= EPS_REQ)),
                   max_gauss=float(max(x['gauss'] for x in ens)), max_ledger=float(max(x['ledger'] for x in ens)))
        if ens[0].get('eps_tau316') is not None:
            e3 = np.array([x['eps_tau316'] for x in ens]); rec['eps_mean_tau316'] = float(e3.mean())
    return rec


def main():
    quick = 'quick' in sys.argv
    N = 16 if quick else 128
    t0 = time.time()
    deltas = [np.pi / 4, np.pi / 2, 3 * np.pi / 4, -np.pi / 2]
    nat = {op: g.ratio_for_coefficient(op, g.coefficient_from_wilson(op, 1.0), Z_ALIGNED) for op in (1, 2, 3, 4)}
    grid = {4: [nat[4], 0.01, 0.03, 0.1, 0.3, 1.0], 2: [0.1, 1.0, nat[2], 100.0], 3: [0.1, nat[3], 10.0], 1: [nat[1], 0.1, 1.0, 10.0]}
    out = dict(scope='Conditional homogeneous classical model; supplied release; not a cosmological prediction.',
               eps_requirement=EPS_REQ, natural_ratio_c1_Mred=nat, configs=[])
    with Pool(4) as pool:
        out['configs'].append(configuration(0, 0.0, 0.0, N=2 * N, pool=pool))
        for op in (4, 2, 3, 1):
            for r in grid[op]:
                for d in deltas:
                    rec = configuration(op, r, d, N=N, pool=pool)
                    out['configs'].append(rec)
                    print(op, '%.3g' % r, '%.3f' % d, rec.get('status', ''), 'mean=%.4f+-%.4f std=%.3f max|e|=%.3f net=%.4f fail=%s det=%.4f' % (
                        rec.get('eps_mean', np.nan), rec.get('eps_sem', np.nan), rec.get('eps_std', np.nan), rec.get('eps_abs_max', np.nan),
                        rec.get('eps_net', np.nan), rec.get('n_runaway_or_fail'), rec.get('deterministic', {}).get('eps', np.nan)), flush=True)
    out['runtime_s'] = time.time() - t0
    (OUT / ('phase_scan_quick.json' if quick else 'phase_scan.json')).write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
